// apps/sjs_run.cpp
//
// HighDims single experiment run:
//   - infer/resolve dim
//   - load/generate dataset
//   - create baseline (method + variant)
//   - run protocol with controlled RNG seeds
//   - write raw results CSV (+ optional sampled pairs)
//
// Dim handling:
//   - synthetic: uses --dim (default 2)
//   - binary: reads dim from binary header (R/S must match)
//   - csv: infers dim from first data row (R/S must match)
//
// Requires:
//   - sjs/dispatch/dim_dispatch.h : sjs::dispatch::DispatchDim<R>(dim, fn, out, err)
//   - sjs/baselines/baseline_factory.h : CreateBaseline<Dim>(method, variant, err)

#include "sjs/baselines/baseline_factory.h"

#include "sjs/baselines/runners/adaptive_runner.h"
#include "sjs/baselines/runners/enum_sampling_runner.h"
#include "sjs/baselines/runners/sampling_runner.h"

#include "sjs/core/logging.h"
#include "sjs/core/stats.h"
#include "sjs/core/timer.h"
#include "sjs/core/types.h"

#include "sjs/dispatch/dim_dispatch.h"

#include "sjs/data/synthetic/generator_factory.h"
#include "sjs/io/binary_io.h"
#include "sjs/io/csv_io.h"
#include "sjs/io/dataset.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>


// -----------------------------------------------------------------------------
// Minimal CLI parsing utilities.
// HighDims keeps sjs::Config in baseline_api.h; apps parse flags directly.
// -----------------------------------------------------------------------------

namespace sjs {

class ArgMap {
 public:
  static ArgMap FromArgv(int argc, char** argv) {
    ArgMap out;
    if (argc > 0 && argv && argv[0]) out.program_ = argv[0];

    bool end_of_opts = false;
    for (int i = 1; i < argc; ++i) {
      std::string_view tok = argv[i] ? std::string_view(argv[i]) : std::string_view{};
      if (!end_of_opts && tok == "--") {
        end_of_opts = true;
        continue;
      }

      if (!end_of_opts && IsOptionToken(tok)) {
        std::string_view key;
        std::string_view val;
        if (SplitOption(tok, &key, &val)) {
          const std::string k = NormalizeKey(key);
          if (!k.empty()) out.kv_[k] = std::string(val);
          continue;
        }

        // --key value OR --flag
        const std::string k = NormalizeKey(tok);
        if (k.empty()) continue;

        if (i + 1 < argc) {
          std::string_view nxt = argv[i + 1] ? std::string_view(argv[i + 1]) : std::string_view{};
          // Treat "-1" (or any numeric token) as a value, even though it starts with '-'.
          if (!IsOptionToken(nxt) || LooksLikeNumber(nxt)) {
            out.kv_[k] = std::string(nxt);
            ++i;
            continue;
          }
        }
        out.kv_[k] = "1";
        continue;
      }

      out.pos_.push_back(std::string(tok));
    }

    return out;
  }

  bool Has(std::string_view key) const {
    const std::string k = NormalizeKey(key);
    return kv_.find(k) != kv_.end();
  }

  std::optional<std::string_view> Get(std::string_view key) const {
    const std::string k = NormalizeKey(key);
    auto it = kv_.find(k);
    if (it == kv_.end()) return std::nullopt;
    return std::string_view(it->second);
  }

  const std::vector<std::string>& Positional() const { return pos_; }

  const std::unordered_map<std::string, std::string>& Items() const { return kv_; }

 private:
  static bool IsOptionToken(std::string_view s) noexcept {
    return s.size() >= 2 && s[0] == '-' && s != "-";
  }

  static bool SplitOption(std::string_view tok,
                          std::string_view* key_out,
                          std::string_view* val_out) noexcept {
    const size_t eq = tok.find('=');
    if (eq == std::string_view::npos) return false;
    if (key_out) *key_out = tok.substr(0, eq);
    if (val_out) *val_out = tok.substr(eq + 1);
    return true;
  }

  static std::string NormalizeKey(std::string_view key) {
    while (!key.empty() && key.front() == '-') key.remove_prefix(1);
    return std::string(key);
  }

  static bool LooksLikeNumber(std::string_view s) noexcept {
    if (s.empty()) return false;
    size_t i = 0;
    if (s[i] == '+' || s[i] == '-') {
      if (s.size() == 1) return false;
      ++i;
    }
    bool has_digit = false;
    for (; i < s.size(); ++i) {
      const char c = s[i];
      if (c >= '0' && c <= '9') {
        has_digit = true;
        continue;
      }
      if (c == '.' || c == 'e' || c == 'E' || c == '+' || c == '-') continue;
      return false;
    }
    return has_digit;
  }

  std::string program_;
  std::unordered_map<std::string, std::string> kv_;
  std::vector<std::string> pos_;
};

inline bool ParseDataSource(std::string_view s, DataSource* out) noexcept {
  if (!out) return false;
  if (detail::EqualsIgnoreCase(s, "synthetic")) { *out = DataSource::Synthetic; return true; }
  if (detail::EqualsIgnoreCase(s, "binary")) { *out = DataSource::Binary; return true; }
  if (detail::EqualsIgnoreCase(s, "csv")) { *out = DataSource::CSV; return true; }
  if (detail::EqualsIgnoreCase(s, "unknown")) { *out = DataSource::Unknown; return true; }
  return false;
}

namespace detail {
inline bool ParseLogLevel(std::string_view s, LogLevel* out) noexcept {
  if (!out) return false;
  if (EqualsIgnoreCase(s, "trace")) { *out = LogLevel::Trace; return true; }
  if (EqualsIgnoreCase(s, "debug")) { *out = LogLevel::Debug; return true; }
  if (EqualsIgnoreCase(s, "info")) { *out = LogLevel::Info; return true; }
  if (EqualsIgnoreCase(s, "warn") || EqualsIgnoreCase(s, "warning")) { *out = LogLevel::Warn; return true; }
  if (EqualsIgnoreCase(s, "error")) { *out = LogLevel::Error; return true; }
  if (EqualsIgnoreCase(s, "off")) { *out = LogLevel::Off; return true; }
  return false;
}
}  // namespace detail

namespace cli {

inline bool TryParseBool(std::string_view s, bool* out) noexcept {
  if (!out) return false;
  if (detail::EqualsIgnoreCase(s, "1") || detail::EqualsIgnoreCase(s, "true") ||
      detail::EqualsIgnoreCase(s, "yes") || detail::EqualsIgnoreCase(s, "y") ||
      detail::EqualsIgnoreCase(s, "on")) {
    *out = true;
    return true;
  }
  if (detail::EqualsIgnoreCase(s, "0") || detail::EqualsIgnoreCase(s, "false") ||
      detail::EqualsIgnoreCase(s, "no") || detail::EqualsIgnoreCase(s, "n") ||
      detail::EqualsIgnoreCase(s, "off")) {
    *out = false;
    return true;
  }
  return false;
}

inline bool TryParseU64(std::string_view s, u64* out) noexcept {
  if (!out) return false;
  try {
    std::size_t idx = 0;
    unsigned long long v = std::stoull(std::string(s), &idx, 10);
    if (idx != s.size()) return false;
    *out = static_cast<u64>(v);
    return true;
  } catch (...) {
    return false;
  }
}

inline bool TryParseI32(std::string_view s, i32* out) noexcept {
  if (!out) return false;
  try {
    std::size_t idx = 0;
    long long v = std::stoll(std::string(s), &idx, 10);
    if (idx != s.size()) return false;
    if (v < static_cast<long long>(std::numeric_limits<i32>::min()) ||
        v > static_cast<long long>(std::numeric_limits<i32>::max())) {
      return false;
    }
    *out = static_cast<i32>(v);
    return true;
  } catch (...) {
    return false;
  }
}

inline bool TryParseDouble(std::string_view s, double* out) noexcept {
  if (!out) return false;
  try {
    std::size_t idx = 0;
    double v = std::stod(std::string(s), &idx);
    if (idx != s.size()) return false;
    *out = v;
    return true;
  } catch (...) {
    return false;
  }
}

inline bool IsReservedKey(std::string_view k) noexcept {
  // Core config keys.
  if (k == "dataset_source" || k == "dataset" || k == "dim" || k == "path_r" || k == "path_s") return true;
  if (k == "gen" || k == "n_r" || k == "n_s" || k == "alpha" || k == "gen_seed") return true;
  if (k == "method" || k == "variant" || k == "t" || k == "seed" || k == "repeats" ||
      k == "j_star" || k == "budget" || k == "w_small" || k == "prefetch" ||
      k == "enum_cap" || k == "write_samples" || k == "verify") return true;
  if (k == "out_dir" || k == "run_tag" || k == "threads") return true;
  if (k == "log_level" || k == "with_timestamp" || k == "with_thread_id") return true;
  if (k == "csv_sep") return true;

  // Common app-only keys (do not forward to extra).
  if (k == "help" || k == "h") return true;
  if (k == "results_file" || k == "config" || k == "raw_file" || k == "summary_file") return true;
  if (k == "write_csv" || k == "out_r" || k == "out_s" || k == "csv_r" || k == "csv_s") return true;
  if (k == "oracle_max_checks" || k == "oracle_collect_limit" || k == "oracle_cap") return true;

  return false;
}

inline bool ParseConfigFromArgs(const ArgMap& args, Config* out, std::string* err) {
  if (!out) {
    if (err) *err = "ParseConfigFromArgs: out is null";
    return false;
  }

  Config cfg;  // defaults from baseline_api.h

  auto fail = [&](std::string_view msg) {
    if (err) *err = std::string(msg);
    return false;
  };

  // dataset
  if (auto v = args.Get("dataset_source")) {
    DataSource src = DataSource::Unknown;
    if (!ParseDataSource(*v, &src)) {
      return fail("Invalid --dataset_source (expected synthetic|binary|csv).");
    }
    cfg.dataset.source = src;
  }
  if (auto v = args.Get("dataset")) cfg.dataset.name = std::string(*v);
  if (auto v = args.Get("dim")) {
    i32 d = 0;
    if (!TryParseI32(*v, &d)) return fail("Invalid --dim (expected i32).");
    cfg.dataset.dim = d;
  }
  if (auto v = args.Get("path_r")) cfg.dataset.path_r = std::string(*v);
  if (auto v = args.Get("path_s")) cfg.dataset.path_s = std::string(*v);

  // synthetic
  if (auto v = args.Get("gen")) cfg.dataset.synthetic.generator = std::string(*v);
  if (auto v = args.Get("n_r")) {
    u64 x = 0;
    if (!TryParseU64(*v, &x)) return fail("Invalid --n_r (expected u64).");
    cfg.dataset.synthetic.n_r = x;
  }
  if (auto v = args.Get("n_s")) {
    u64 x = 0;
    if (!TryParseU64(*v, &x)) return fail("Invalid --n_s (expected u64).");
    cfg.dataset.synthetic.n_s = x;
  }
  if (auto v = args.Get("alpha")) {
    double x = 0.0;
    if (!TryParseDouble(*v, &x)) return fail("Invalid --alpha (expected number).");
    cfg.dataset.synthetic.alpha = x;
  }
  if (auto v = args.Get("gen_seed")) {
    u64 x = 0;
    if (!TryParseU64(*v, &x)) return fail("Invalid --gen_seed (expected u64).");
    cfg.dataset.synthetic.seed = x;
  }

  // run
  if (auto v = args.Get("method")) {
    Method m = Method::Unknown;
    if (!ParseMethod(*v, &m) || m == Method::Unknown) {
      return fail("Invalid --method (expected ours|range_tree|kd_tree).");
    }
    cfg.run.method = m;
  }
  if (auto v = args.Get("variant")) {
    Variant vv = Variant::Sampling;
    if (!ParseVariant(*v, &vv)) {
      return fail("Invalid --variant (expected sampling|enum_sampling|adaptive).");
    }
    cfg.run.variant = vv;
  }
  if (auto v = args.Get("t")) {
    u64 x = 0;
    if (!TryParseU64(*v, &x)) return fail("Invalid --t (expected u64).");
    cfg.run.t = x;
  }
  if (auto v = args.Get("seed")) {
    u64 x = 0;
    if (!TryParseU64(*v, &x)) return fail("Invalid --seed (expected u64).");
    cfg.run.seed = x;
  }
  if (auto v = args.Get("repeats")) {
    u64 x = 0;
    if (!TryParseU64(*v, &x)) return fail("Invalid --repeats (expected u64).");
    cfg.run.repeats = x;
  }
  // Framework III (Variant::Adaptive) knobs.
  //
  // IMPORTANT (HighDims alignment):
  //   - "budget B" controls how many join pairs we may cache/prefetch (performance-only).
  //   - Historical flag name "--j_star" is kept as a DEPRECATED alias of "--budget".
  //
  // Parse --budget / --j_star (alias) with conflict checking.
  bool has_j_star = false;
  u64 j_star_val = cfg.run.j_star;
  if (auto v = args.Get("j_star")) {
    u64 x = 0;
    if (!TryParseU64(*v, &x)) return fail("Invalid --j_star (expected u64).");
    has_j_star = true;
    j_star_val = x;
  }

  bool has_budget = false;
  u64 budget_val = 0;
  if (auto v = args.Get("budget")) {
    u64 x = 0;
    if (!TryParseU64(*v, &x)) return fail("Invalid --budget (expected u64).");
    has_budget = true;
    budget_val = x;
  }

  if (has_budget && has_j_star && budget_val != j_star_val) {
    return fail("Conflicting --budget and --j_star: they must match (prefer --budget). ");
  }

  if (has_budget) {
    // Canonical path: --budget wins.
    cfg.run.budget = budget_val;
    cfg.run.j_star = budget_val;  // legacy alias (kept consistent for JSON/logging)
    cfg.run.extra["budget"] = std::to_string(budget_val);
  } else if (has_j_star) {
    // Legacy alias: still forward as canonical budget key.
    cfg.run.budget = j_star_val;
    cfg.run.j_star = j_star_val;
    cfg.run.extra["budget"] = std::to_string(j_star_val);
  }


  // Small-block full-cache threshold (w_i <= w_small). 0 disables full-cache.
  if (auto v = args.Get("w_small")) {
    u64 x = 0;
    if (!TryParseU64(*v, &x)) return fail("Invalid --w_small (expected u64).");
    cfg.run.w_small = x;
    cfg.run.extra["w_small"] = std::to_string(x);
  }


  // Prefetch toggle (performance-only). Default: enabled.
  // NOTE: Prefetch affects only performance; disabling it should not change the output distribution.
  if (auto v = args.Get("prefetch")) {
    bool b = true;
    if (!TryParseBool(*v, &b)) return fail("Invalid --prefetch (expected 0/1/true/false).");
    cfg.run.prefetch = b;
    cfg.run.extra["prefetch"] = b ? "1" : "0";
  }

  if (auto v = args.Get("enum_cap")) {
    u64 x = 0;
    if (!TryParseU64(*v, &x)) return fail("Invalid --enum_cap (expected u64).");
    cfg.run.enum_cap = x;
  }
  if (auto v = args.Get("write_samples")) {
    bool b = false;
    if (!TryParseBool(*v, &b)) return fail("Invalid --write_samples (expected 0/1/true/false).");
    cfg.run.write_samples = b;
  }
  if (auto v = args.Get("verify")) {
    bool b = false;
    if (!TryParseBool(*v, &b)) return fail("Invalid --verify (expected 0/1/true/false).");
    cfg.run.verify = b;
  }

  // output
  if (auto v = args.Get("out_dir")) cfg.output.out_dir = std::string(*v);
  if (auto v = args.Get("run_tag")) cfg.output.run_tag = std::string(*v);

  // sys
  if (auto v = args.Get("threads")) {
    i32 x = 0;
    if (!TryParseI32(*v, &x) || x <= 0) return fail("Invalid --threads (expected i32>0).");
    cfg.sys.threads = x;
  }

  // logging
  if (auto v = args.Get("log_level")) {
    LogLevel lvl = LogLevel::Info;
    if (!detail::ParseLogLevel(*v, &lvl)) return fail("Invalid --log_level (trace|debug|info|warn|error|off).");
    cfg.logging.level = lvl;
  }
  if (auto v = args.Get("with_timestamp")) {
    bool b = true;
    if (!TryParseBool(*v, &b)) return fail("Invalid --with_timestamp (expected 0/1/true/false).");
    cfg.logging.with_timestamp = b;
  }
  if (auto v = args.Get("with_thread_id")) {
    bool b = false;
    if (!TryParseBool(*v, &b)) return fail("Invalid --with_thread_id (expected 0/1/true/false).");
    cfg.logging.with_thread_id = b;
  }

  // Special: csv_sep is used by CSV IO and is stored under run.extra["csv_sep"].
  if (auto v = args.Get("csv_sep")) cfg.run.extra["csv_sep"] = std::string(*v);

  // Forward unknown flags into:
  //   - cfg.run.extra (baseline knobs)
  //   - cfg.dataset.synthetic.extra (generator knobs) when source=synthetic
  for (const auto& kv : args.Items()) {
    const std::string& k = kv.first;
    const std::string& v = kv.second;
    if (IsReservedKey(k)) continue;
    cfg.run.extra[k] = v;
    if (cfg.dataset.source == DataSource::Synthetic) cfg.dataset.synthetic.extra[k] = v;
  }

  *out = std::move(cfg);
  return true;
}

}  // namespace cli

inline bool ParseConfigFromArgs(const ArgMap& args, Config* out, std::string* err) {
  return cli::ParseConfigFromArgs(args, out, err);
}

}  // namespace sjs

namespace fs = std::filesystem;

namespace sjs {
namespace apps {

namespace {

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

inline bool IsHelpRequested(const sjs::ArgMap& args) {
  return args.Has("help") || args.Has("h") || args.Has("-h") || args.Has("--help");
}

inline std::string SanitizeFilename(std::string_view s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    const bool ok =
        (std::isalnum(static_cast<unsigned char>(c)) != 0) || c == '_' || c == '-' || c == '.';
    out.push_back(ok ? c : '_');
  }
  if (out.empty()) out = "x";
  return out;
}

inline bool EnsureDir(const fs::path& p, std::string* err) {
  try {
    if (p.empty()) return true;
    fs::create_directories(p);
    return true;
  } catch (const std::exception& e) {
    SetErr(err, std::string("create_directories failed: ") + e.what());
    return false;
  }
}

inline fs::path ResolvePathByWalkingUp(const fs::path& p, int max_up = 6) {
  if (p.empty()) return p;
  if (p.is_absolute()) return p;
  if (fs::exists(p)) return p;

  fs::path base = fs::current_path();
  for (int i = 0; i < max_up; ++i) {
    const fs::path cand = base / p;
    if (fs::exists(cand)) return cand;

    if (!base.has_parent_path()) break;
    const fs::path parent = base.parent_path();
    if (parent == base) break;
    base = parent;
  }
  return p;
}

inline bool PeekBinaryDim(const fs::path& path, int* out_dim, std::string* err) {
  if (!out_dim) {
    SetErr(err, "PeekBinaryDim: out_dim is null");
    return false;
  }
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    SetErr(err, "PeekBinaryDim: cannot open file: " + path.string());
    return false;
  }
  sjs::binary::FileHeader h{};
  in.read(reinterpret_cast<char*>(&h), static_cast<std::streamsize>(sizeof(h)));
  if (!in) {
    SetErr(err, "PeekBinaryDim: failed to read header: " + path.string());
    return false;
  }
  if (std::memcmp(h.magic, sjs::binary::kMagic, sizeof(sjs::binary::kMagic)) != 0) {
    SetErr(err, "PeekBinaryDim: bad magic: " + path.string());
    return false;
  }
  if (h.endian != sjs::binary::kEndianMarker) {
    SetErr(err, "PeekBinaryDim: endian marker mismatch (file not little-endian?): " + path.string());
    return false;
  }
  if (h.dim == 0 || h.dim > static_cast<sjs::u32>(sjs::kMaxSupportedDim)) {
    SetErr(err, "PeekBinaryDim: invalid dim in header: " + std::to_string(h.dim));
    return false;
  }
  *out_dim = static_cast<int>(h.dim);
  return true;
}

inline std::string_view Trim(std::string_view s) {
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.remove_suffix(1);
  return s;
}

inline std::vector<std::string_view> SplitLine(std::string_view line, char sep) {
  std::vector<std::string_view> out;
  sjs::usize start = 0;
  while (start <= line.size()) {
    const sjs::usize pos = line.find(sep, start);
    if (pos == std::string_view::npos) {
      out.push_back(line.substr(start));
      break;
    }
    out.push_back(line.substr(start, pos - start));
    start = pos + 1;
  }
  return out;
}

inline bool InferCsvDim(const fs::path& path,
                        char sep,
                        bool has_header,
                        int* out_dim,
                        std::string* err) {
  if (!out_dim) {
    SetErr(err, "InferCsvDim: out_dim is null");
    return false;
  }
  std::ifstream in(path);
  if (!in) {
    SetErr(err, "InferCsvDim: cannot open file: " + path.string());
    return false;
  }
  std::string line;
  bool skipped_header = false;
  while (std::getline(in, line)) {
    std::string_view sv = Trim(std::string_view(line));
    if (sv.empty()) continue;
    if (!sv.empty() && sv.front() == '#') continue;

    if (has_header && !skipped_header) {
      skipped_header = true;
      continue;
    }

    const auto cols = SplitLine(sv, sep);
    const sjs::usize c = cols.size();
    if (c < 4) {
      SetErr(err, "InferCsvDim: too few columns in first data row: " + path.string());
      return false;
    }

    // Shapes supported by csv::ReadBoxesSimple:
    //   - without id: 2*Dim
    //   - with id: 1 + 2*Dim
    if (c % 2 == 0) {
      *out_dim = static_cast<int>(c / 2);
      return true;
    }
    if ((c - 1) % 2 == 0) {
      *out_dim = static_cast<int>((c - 1) / 2);
      return true;
    }
    SetErr(err, "InferCsvDim: column count not compatible with (2*Dim) or (1+2*Dim): " +
                    std::to_string(static_cast<unsigned long long>(c)));
    return false;
  }
  SetErr(err, "InferCsvDim: no data rows found: " + path.string());
  return false;
}

inline bool DetermineDatasetDim(sjs::Config* cfg, const sjs::ArgMap& args, std::string* err) {
  if (!cfg) {
    SetErr(err, "DetermineDatasetDim: cfg is null");
    return false;
  }

  // For synthetic: dim is taken from cfg.dataset.dim (default 2).
  if (cfg->dataset.source == sjs::DataSource::Synthetic) {
    if (cfg->dataset.dim <= 0) {
      SetErr(err, "Synthetic dataset requires a positive --dim");
      return false;
    }
    return true;
  }

  // Non-synthetic: infer dim from inputs (binary header or CSV first data row).
  if (cfg->dataset.path_r.empty() || cfg->dataset.path_s.empty()) {
    SetErr(err, "Non-synthetic dataset requires --path_r and --path_s");
    return false;
  }

  const fs::path path_r = ResolvePathByWalkingUp(fs::path(cfg->dataset.path_r));
  const fs::path path_s = ResolvePathByWalkingUp(fs::path(cfg->dataset.path_s));

  int dim_r = 0, dim_s = 0;
  std::string local_err;

  if (cfg->dataset.source == sjs::DataSource::Binary) {
    if (!PeekBinaryDim(path_r, &dim_r, &local_err)) {
      SetErr(err, "Failed to read dim from binary R: " + local_err);
      return false;
    }
    if (!PeekBinaryDim(path_s, &dim_s, &local_err)) {
      SetErr(err, "Failed to read dim from binary S: " + local_err);
      return false;
    }
  } else if (cfg->dataset.source == sjs::DataSource::CSV) {
    char sep = ',';
    if (auto it = cfg->run.extra.find("csv_sep"); it != cfg->run.extra.end()) {
      if (!it->second.empty()) {
        if (it->second == "tab" || it->second == "\\t") sep = '\t';
        else sep = it->second[0];
      }
    }
    if (!InferCsvDim(path_r, sep, /*has_header=*/true, &dim_r, &local_err)) {
      SetErr(err, "Failed to infer dim from CSV R: " + local_err);
      return false;
    }
    if (!InferCsvDim(path_s, sep, /*has_header=*/true, &dim_s, &local_err)) {
      SetErr(err, "Failed to infer dim from CSV S: " + local_err);
      return false;
    }
  } else {
    SetErr(err, "DetermineDatasetDim: unsupported dataset_source");
    return false;
  }

  if (dim_r != dim_s) {
    SetErr(err, "R/S dim mismatch: dim(R)=" + std::to_string(dim_r) +
                    " dim(S)=" + std::to_string(dim_s));
    return false;
  }

  // If user explicitly passed --dim, enforce it matches inferred.
  if (args.Has("dim")) {
    if (cfg->dataset.dim != dim_r) {
      SetErr(err, "User-specified --dim=" + std::to_string(cfg->dataset.dim) +
                      " does not match inferred dim=" + std::to_string(dim_r));
      return false;
    }
  } else {
    cfg->dataset.dim = dim_r;
  }

  // Rewrite paths to resolved ones for downstream load.
  cfg->dataset.path_r = path_r.string();
  cfg->dataset.path_s = path_s.string();

  return true;
}

template <int Dim>
bool LoadOrGenerateDataset(const sjs::Config& cfg,
                           sjs::Dataset<Dim, sjs::Scalar>* out,
                           sjs::synthetic::Report* gen_report,
                           std::string* err) {
  if (!out) {
    SetErr(err, "LoadOrGenerateDataset: out is null");
    return false;
  }
  if (cfg.dataset.dim != Dim) {
    SetErr(err, "Config dim mismatch: cfg.dataset.dim != compiled Dim");
    return false;
  }

  if (cfg.dataset.source == sjs::DataSource::Synthetic) {
    sjs::synthetic::DatasetSpec spec;
    spec.name = cfg.dataset.name;
    spec.n_r = cfg.dataset.synthetic.n_r;
    spec.n_s = cfg.dataset.synthetic.n_s;
    spec.alpha = cfg.dataset.synthetic.alpha;
    spec.seed = cfg.dataset.synthetic.seed;
    spec.params = cfg.dataset.synthetic.extra;

    sjs::synthetic::Report rep;
    std::string local_err;
    if (!sjs::synthetic::GenerateDataset<Dim, sjs::Scalar>(
            cfg.dataset.synthetic.generator, spec, out, &rep, &local_err)) {
      SetErr(err, "Synthetic generation failed: " + local_err);
      return false;
    }
    if (gen_report) *gen_report = rep;
    return true;
  }

  if (cfg.dataset.source == sjs::DataSource::Binary) {
    sjs::Relation<Dim, sjs::Scalar> R, S;
    std::string local_err;

    sjs::binary::RelationFileInfo infoR, infoS;
    sjs::binary::BinaryReadOptions opt;
    opt.generate_ids_if_missing = true;
    opt.drop_empty = false;

    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(
            cfg.dataset.path_r, &R, &infoR, opt, &local_err)) {
      SetErr(err, "Failed reading binary R from " + cfg.dataset.path_r + ": " + local_err);
      return false;
    }
    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(
            cfg.dataset.path_s, &S, &infoS, opt, &local_err)) {
      SetErr(err, "Failed reading binary S from " + cfg.dataset.path_s + ": " + local_err);
      return false;
    }

    sjs::Dataset<Dim, sjs::Scalar> ds;
    ds.name = cfg.dataset.name.empty() ? "binary" : cfg.dataset.name;
    ds.half_open = true;
    ds.R = std::move(R);
    ds.S = std::move(S);

    std::string v;
    if (!ds.Validate(/*require_proper=*/true, &v)) {
      SetErr(err, "Binary dataset validation failed: " + v);
      return false;
    }
    *out = std::move(ds);
    return true;
  }

  if (cfg.dataset.source == sjs::DataSource::CSV) {
    sjs::Relation<Dim, sjs::Scalar> R, S;
    std::string local_err;

    char sep = ',';
    if (auto it = cfg.run.extra.find("csv_sep"); it != cfg.run.extra.end()) {
      if (!it->second.empty()) {
        if (it->second == "tab" || it->second == "\\t") sep = '\t';
        else sep = it->second[0];
      }
    }

    if (!sjs::csv::ReadBoxesSimple<Dim, sjs::Scalar>(
            cfg.dataset.path_r, &R, sep, /*has_header=*/true, &local_err)) {
      SetErr(err, "Failed reading CSV R from " + cfg.dataset.path_r + ": " + local_err);
      return false;
    }
    if (!sjs::csv::ReadBoxesSimple<Dim, sjs::Scalar>(
            cfg.dataset.path_s, &S, sep, /*has_header=*/true, &local_err)) {
      SetErr(err, "Failed reading CSV S from " + cfg.dataset.path_s + ": " + local_err);
      return false;
    }

    sjs::Dataset<Dim, sjs::Scalar> ds;
    ds.name = cfg.dataset.name.empty() ? "csv" : cfg.dataset.name;
    ds.half_open = true;
    ds.R = std::move(R);
    ds.S = std::move(S);

    std::string v;
    if (!ds.Validate(/*require_proper=*/true, &v)) {
      SetErr(err, "CSV dataset validation failed: " + v);
      return false;
    }
    *out = std::move(ds);
    return true;
  }

  SetErr(err, "Unsupported dataset_source");
  return false;
}

inline std::vector<std::string> ResultHeader() {
  return {
      "dataset",
      "dim",
      "method",
      "variant",
      "rep",
      "seed",
      "n_r",
      "n_s",
      "t",
      "ok",
      "wall_ms",
      "count_value",
      "count_exact",
      "count_stderr",
      "count_ci_low",
      "count_ci_high",
      "used_enumeration",
      "enum_truncated",
      "enum_cap",
      "adaptive_branch",
      "adaptive_pilot_pairs",
      "phases_json",
      "note",
      "config_json",
      "gen_report_json",
  };
}

// -----------------------------------------------------------------------------
// Config JSON helpers
//
// Config::ToJsonLite() intentionally omits stringly-typed knobs (run.extra), but
// for HighDims experiments those keys often contain essential reproducibility
// knobs (e.g., Framework III budget/w_small/prefetch, and generator params).
//
// sjs_run extends the per-row config_json with a top-level "run_extra" object.
// -----------------------------------------------------------------------------

inline std::string JsonEscape(std::string_view s) {
  std::string out;
  out.reserve(s.size() + 8);
  static constexpr char kHex[] = "0123456789abcdef";

  for (char c : s) {
    switch (c) {
      case '\\': out += "\\\\"; break;
      case '"': out += "\\\""; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default: {
        const unsigned char uc = static_cast<unsigned char>(c);
        if (uc < 0x20U) {
          // Control character -> \u00XX
          out += "\\u00";
          out.push_back(kHex[(uc >> 4) & 0x0FU]);
          out.push_back(kHex[uc & 0x0FU]);
        } else {
          out.push_back(c);
        }
      } break;
    }
  }
  return out;
}

inline std::string MapToJsonObject(const std::unordered_map<std::string, std::string>& m) {
  if (m.empty()) return "{}";

  std::vector<std::pair<std::string, std::string>> items;
  items.reserve(m.size());
  for (const auto& kv : m) items.emplace_back(kv.first, kv.second);
  std::sort(items.begin(), items.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

  std::ostringstream oss;
  oss << "{";
  bool first = true;
  for (const auto& kv : items) {
    if (!first) oss << ",";
    first = false;
    oss << "\"" << JsonEscape(kv.first) << "\":\"" << JsonEscape(kv.second) << "\"";
  }
  oss << "}";
  return oss.str();
}

inline std::string ConfigToJsonWithRunExtra(const sjs::Config& cfg) {
  std::string base = cfg.ToJsonLite();
  if (cfg.run.extra.empty()) return base;
  if (!base.empty() && base.back() == '}') {
    base.pop_back();
    base += ",\"run_extra\":";
    base += MapToJsonObject(cfg.run.extra);
    base += "}";
  }
  return base;
}

template <int Dim>
inline bool WriteResultRow(sjs::csv::Writer& w,
                           const sjs::Config& cfg,
                           const sjs::Dataset<Dim, sjs::Scalar>& ds,
                           const sjs::baselines::RunReport& r,
                           int rep,
                           double wall_ms,
                           const sjs::synthetic::Report* gen_report) {
  const std::string gen_json = gen_report ? gen_report->ToJsonLite() : "{}";
  const std::string cfg_json = ConfigToJsonWithRunExtra(cfg);

  return w.WriteRowV(
      ds.name,
      Dim,
      sjs::ToString(r.method),
      sjs::ToString(r.variant),
      rep,
      static_cast<unsigned long long>(r.seed),
      static_cast<unsigned long long>(ds.R.Size()),
      static_cast<unsigned long long>(ds.S.Size()),
      static_cast<unsigned long long>(r.t),
      (r.ok ? 1 : 0),
      wall_ms,
      static_cast<double>(r.count.value),
      (r.count.exact ? 1 : 0),
      static_cast<double>(r.count.stderr),
      static_cast<double>(r.count.ci_low),
      static_cast<double>(r.count.ci_high),
      (r.used_enumeration ? 1 : 0),
      (r.enumeration_truncated ? 1 : 0),
      static_cast<unsigned long long>(r.enumeration_cap),
      r.adaptive_branch,
      static_cast<unsigned long long>(r.adaptive_pilot_pairs),
      r.phases.ToJsonMillis(),
      r.note,
      cfg_json,
      gen_json);
}

inline void PrintUsage() {
  std::cerr
      << "sjs_run: single experiment runner (HighDims)\n\n"
      << "Common flags:\n"
      << "  --method=<ours|range_tree|kd_tree>\n"
      << "  --variant=<sampling|enum_sampling|adaptive>\n"
      << "  --t=<num_samples>\n"
      << "  --seed=<seed>           (base seed; repeats add +rep)\n"
      << "  --repeats=<k>\n"
      << "  --out_dir=<dir>\n"
      << "  --results_file=<path>   (default: <out_dir>/run.csv)\n"
      << "  --write_samples=0|1\n"
      << "  --enum_cap=<cap>        (EnumSampling/Adaptive; 0 means no cap)\n"
      << "\nAdaptive (Framework III) knobs:\n"
      << "  --budget=<B>            (cache/prefetch budget in #pairs; preferred)\n"
      << "  --w_small=<w>           (full-cache threshold: cache blocks with w_i<=w; 0 disables)\n"
      << "  --prefetch=0|1          (enable prefetch sample cache; default 1)\n"
      << "  --j_star=<B>            (DEPRECATED alias of --budget; kept for compatibility)\n"
      << "\nDataset flags:\n"
      << "  --dataset_source=<synthetic|binary|csv>\n"
      << "  --dataset=<name>\n"
      << "  --dim=<d>               (synthetic required; binary/csv inferred unless provided)\n"
      << "  --path_r=<file> --path_s=<file>    (for binary/csv)\n"
      << "\nSynthetic flags:\n"
      << "  --gen=<stripe|uniform|clustered|hetero_sizes>\n"
      << "  --n_r=<N> --n_s=<N>\n"
      << "  --alpha=<float>\n"
      << "  --gen_seed=<seed>\n"
      << "  plus generator-specific params, e.g. --gap_factor=0.1 --core_lo=0.45 ...\n"
      << "\nNotes:\n"
      << "  - For binary datasets, dim is read from SJS binary header.\n"
      << "  - For CSV datasets, dim is inferred from the first data row.\n"
      << "  - If you compiled only a subset of dims, you must rebuild with those dims enabled.\n";
}

// Local runner helper: allow disabling Framework III prefetch without touching
// baseline internals. Prefetch is performance-only and should not affect the
// output distribution.
template <int Dim, class T = sjs::Scalar>
bool RunAdaptiveOnceMaybeNoPrefetch(sjs::baselines::IBaseline<Dim, T>* baseline,
                                   const sjs::Dataset<Dim, T>& dataset,
                                   const sjs::Config& cfg,
                                   sjs::u64 seed,
                                   bool enable_prefetch,
                                   sjs::baselines::RunReport* out,
                                   std::string* err = nullptr) {
  if (!baseline) {
    if (err) *err = "RunAdaptiveOnceMaybeNoPrefetch: baseline is null";
    return false;
  }
  if (!out) {
    if (err) *err = "RunAdaptiveOnceMaybeNoPrefetch: out is null";
    return false;
  }

  out->ok = false;
  out->error.clear();
  out->method = baseline->method();
  out->variant = sjs::Variant::Adaptive;
  out->baseline_name = std::string(baseline->Name());
  out->dataset_name = dataset.name;
  out->seed = seed;
  out->t = cfg.run.t;

  out->count = sjs::baselines::CountResult{};
  out->samples.Clear();

  out->used_enumeration = false;
  out->enumeration_truncated = false;
  out->enumeration_cap = cfg.run.enum_cap;

  out->enumeration_pairs_pass1 = 0;
  out->enumeration_pairs_pass2 = 0;
  out->enum_stats_pass1 = sjs::join::JoinStats{};
  out->enum_stats_pass2 = sjs::join::JoinStats{};

  out->adaptive_branch = "framework_III";
  out->adaptive_pilot_pairs = 0;
  out->note.clear();
  out->phases.Clear();

  std::string local_err;

  // SJS v3: separate RNG streams across phases.
  sjs::Rng rng_count(sjs::DeriveSeed(seed, 1));
  sjs::Rng rng_sample(sjs::DeriveSeed(seed, 2));

  // If prefetch is disabled, pass nullptr for Count() RNG.
  // Our Framework III baselines use rng only for prefetch sampling.
  sjs::Rng* rng_count_ptr = enable_prefetch ? &rng_count : nullptr;

  baseline->Reset();

  {
    auto _ = out->phases.Scoped("run_build");
    if (!baseline->Build(dataset, cfg, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  {
    auto _ = out->phases.Scoped("run_count");
    if (!baseline->Count(cfg, rng_count_ptr, &out->count, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  {
    auto _ = out->phases.Scoped("run_sample");
    if (!baseline->Sample(cfg, &rng_sample, &out->samples, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  // Defensive validation.
  {
    std::string v_err;
    if (!out->samples.Validate(&v_err)) {
      out->error = "SampleSet validation failed: " + v_err;
      if (err) *err = out->error;
      return false;
    }
  }

  // Record the prefetch choice (reproducibility).
  if (!enable_prefetch) {
    if (!out->note.empty()) out->note += "; ";
    out->note += "prefetch=0";
  }

  out->ok = true;
  return true;
}

template <int Dim>
int RunMainForDim(const sjs::Config& cfg_in, const sjs::ArgMap& args) {
  std::string err;

  // Copy cfg so we can tweak minor runtime fields without affecting caller.
  sjs::Config cfg = cfg_in;

  // Load/generate dataset.
  sjs::Dataset<Dim, sjs::Scalar> ds;
  sjs::synthetic::Report gen_report;
  sjs::synthetic::Report* gen_report_ptr = nullptr;

  {
    std::string local_err;
    if (!LoadOrGenerateDataset<Dim>(cfg, &ds, &gen_report, &local_err)) {
      SJS_LOG_ERROR(local_err);
      return 3;
    }
    if (cfg.dataset.source == sjs::DataSource::Synthetic) {
      gen_report_ptr = &gen_report;
      SJS_LOG_INFO("Generated dataset:", ds.name,
                   "Dim=", Dim,
                   "R=", static_cast<unsigned long long>(ds.R.Size()),
                   "S=", static_cast<unsigned long long>(ds.S.Size()),
                   "gen=", cfg.dataset.synthetic.generator,
                   "report=", gen_report.ToJsonLite());
    } else {
      SJS_LOG_INFO("Loaded dataset:", ds.name,
                   "Dim=", Dim,
                   "R=", static_cast<unsigned long long>(ds.R.Size()),
                   "S=", static_cast<unsigned long long>(ds.S.Size()));
    }
  }

  // Create baseline.
  auto baseline = sjs::baselines::CreateBaseline<Dim, sjs::Scalar>(cfg.run.method, cfg.run.variant, &err);
  if (!baseline) {
    SJS_LOG_ERROR("CreateBaseline failed:", err);
    PrintUsage();
    return 4;
  }
  SJS_LOG_INFO("Baseline:", baseline->Name(),
               "method=", sjs::ToString(baseline->method()),
               "variant=", sjs::ToString(baseline->variant()),
               "Dim=", Dim);

  // Output paths.
  const fs::path out_dir(cfg.output.out_dir);
  if (!EnsureDir(out_dir, &err)) {
    SJS_LOG_ERROR("Cannot create out_dir:", err);
    return 5;
  }

  std::string results_path = cfg.output.out_dir + "/run.csv";
  if (auto v = args.Get("results_file")) results_path = std::string(*v);

  if (!EnsureDir(fs::path(results_path).parent_path(), &err)) {
    SJS_LOG_ERROR("Cannot create results parent dir:", err);
    return 5;
  }

  sjs::csv::Writer writer(results_path, sjs::csv::Dialect{','}, &err);
  if (!writer.Ok()) {
    SJS_LOG_ERROR("Cannot open results file:", results_path, "err=", err);
    return 5;
  }

  if (!writer.WriteHeader(ResultHeader(), &err)) {
    SJS_LOG_ERROR("Failed writing CSV header:", err);
    return 5;
  }

  // Optional samples directory.
  fs::path samples_dir = out_dir / "samples";
  const bool write_samples = cfg.run.write_samples;
  if (write_samples) {
    if (!EnsureDir(samples_dir, &err)) {
      SJS_LOG_ERROR("Cannot create samples_dir:", err);
      return 5;
    }
  }

  // Deprecation notice: --j_star used to mean other things in older drafts.
  // In this repo, it is an alias of Framework III budget B.
  if (args.Has("j_star") && !args.Has("budget")) {
    SJS_LOG_WARN("--j_star is deprecated; use --budget. "
                 "In HighDims (Framework III / Variant::Adaptive), it is interpreted as budget B (#pairs). ");
  }

  // Framework III prefetch toggle (performance-only).
  // Default: enabled. If disabled, we pass nullptr for Count() RNG.
  bool enable_prefetch = cfg.run.prefetch;
  // Backward compatibility: some configs may only set run_extra["prefetch"].
  if (auto it = cfg.run.extra.find("prefetch"); it != cfg.run.extra.end()) {
    bool b = enable_prefetch;
    if (sjs::cli::TryParseBool(it->second, &b)) enable_prefetch = b;
  }


  // Run repeats.
  std::vector<double> wall_ms_all;
  std::vector<double> wall_ms_ok;
  sjs::u64 ok_reps = 0;
  wall_ms_all.reserve(static_cast<sjs::usize>(cfg.run.repeats));
  wall_ms_ok.reserve(static_cast<sjs::usize>(cfg.run.repeats));

  for (sjs::u64 rep = 0; rep < cfg.run.repeats; ++rep) {
    const sjs::u64 seed = cfg.run.seed + rep;

    sjs::baselines::RunReport report;
    std::string local_err;

    sjs::Stopwatch sw;
    bool ok = false;

    switch (cfg.run.variant) {
      case sjs::Variant::Sampling:
        ok = sjs::baselines::RunSamplingOnce<Dim, sjs::Scalar>(
            baseline.get(), ds, cfg, seed, &report, &local_err);
        break;
      case sjs::Variant::EnumSampling:
        ok = sjs::baselines::RunEnumSamplingOnce<Dim, sjs::Scalar>(
            baseline.get(), ds, cfg, seed, &report, &local_err);
        break;
      case sjs::Variant::Adaptive:
        ok = RunAdaptiveOnceMaybeNoPrefetch<Dim>(
            baseline.get(), ds, cfg, seed, enable_prefetch, &report, &local_err);
        break;
    }

    const double wall_ms = sw.ElapsedMillis();
    wall_ms_all.push_back(wall_ms);
    if (ok) {
      wall_ms_ok.push_back(wall_ms);
      ++ok_reps;
    }

    if (!ok) {
      SJS_LOG_ERROR("Run failed (rep=", rep, "):", local_err);
      report.ok = false;
      report.error = local_err;
    } else {
      SJS_LOG_INFO("Run ok (rep=", rep, ")",
                   "Dim=", Dim,
                   "wall_ms=", wall_ms,
                   "count=", static_cast<double>(report.count.value),
                   (report.count.exact ? "(exact)" : "(est)"),
                   "samples=", static_cast<unsigned long long>(report.samples.Size()));
    }

    if (!WriteResultRow<Dim>(writer, cfg, ds, report, static_cast<int>(rep), wall_ms, gen_report_ptr)) {
      SJS_LOG_ERROR("Failed writing result row to CSV:", results_path);
      return 6;
    }

    if (write_samples && report.ok && !report.samples.pairs.empty()) {
      // Write sampled pairs as TSV.
      const std::string base =
          SanitizeFilename(ds.name) + "__d" + std::to_string(Dim) + "__" +
          SanitizeFilename(sjs::ToString(report.method)) + "__" +
          SanitizeFilename(sjs::ToString(report.variant)) + "__t" +
          std::to_string(static_cast<unsigned long long>(report.t)) + "__seed" +
          std::to_string(static_cast<unsigned long long>(seed)) + "__rep" +
          std::to_string(static_cast<unsigned long long>(rep));

      const fs::path out_pairs = samples_dir / (base + ".tsv");

      std::string werr;
      if (!sjs::csv::WritePairs(out_pairs.string(),
                               sjs::Span<const sjs::PairId>(report.samples.pairs.data(),
                                                            report.samples.pairs.size()),
                               sjs::csv::Dialect{'\t'}, &werr)) {
        SJS_LOG_WARN("Failed writing samples:", out_pairs.string(), "err=", werr);
      }
    }
  }

  // Print a tiny summary.
  const sjs::Summary sum_all = sjs::Summarize(wall_ms_all);
  SJS_LOG_INFO("Wall-time summary (ms) [all]:", sum_all.ToJson());

  if (!wall_ms_ok.empty()) {
    const sjs::Summary sum_ok = sjs::Summarize(wall_ms_ok);
    const double ok_rate = (cfg.run.repeats > 0)
                               ? (static_cast<double>(ok_reps) / static_cast<double>(cfg.run.repeats))
                               : 0.0;
    SJS_LOG_INFO("Wall-time summary (ms) [ok-only]:", sum_ok.ToJson(), "ok_rate=", ok_rate);
  } else {
    SJS_LOG_WARN("No successful repetitions; ok_rate=0.0");
  }

  SJS_LOG_INFO("Wrote results to:", results_path);
  if (write_samples) {
    SJS_LOG_INFO("Wrote samples to:", (out_dir / "samples").string());
  }

  return 0;
}

}  // namespace

}  // namespace apps
}  // namespace sjs


// -----------------------------------------------------------------------------
// Dim-dispatch runner (namespace-scope: local classes cannot have member templates
// in C++17).
// -----------------------------------------------------------------------------
namespace {

struct RunDimRunner {
  const sjs::Config* cfg = nullptr;
  const sjs::ArgMap* args = nullptr;

  template <int Dim>
  int operator()() const {
    return sjs::apps::RunMainForDim<Dim>(*cfg, *args);
  }
};

}  // namespace

int main(int argc, char** argv) {
  sjs::ArgMap args = sjs::ArgMap::FromArgv(argc, argv);
  if (sjs::apps::IsHelpRequested(args)) {
    sjs::apps::PrintUsage();
    return 0;
  }

  sjs::Config cfg;
  std::string err;
  if (!sjs::ParseConfigFromArgs(args, &cfg, &err)) {
    SJS_LOG_ERROR("Config parse failed:", err);
    sjs::apps::PrintUsage();
    return 2;
  }

  if (!cfg.Validate(&err)) {
    SJS_LOG_ERROR("Config validation failed:", err);
    sjs::apps::PrintUsage();
    return 2;
  }

  // Determine/infer dim and rewrite paths if needed.
  if (!sjs::apps::DetermineDatasetDim(&cfg, args, &err)) {
    SJS_LOG_ERROR("Failed to determine dataset dim:", err);
    return 2;
  }

  // Now validate again with the possibly-updated dim.
  if (!cfg.Validate(&err)) {
    SJS_LOG_ERROR("Config validation failed after dim inference:", err);
    return 2;
  }

  sjs::Logger::Instance().SetConfig(cfg.logging);
  int rc = 0;
  std::string derr;
  if (!sjs::dispatch::DispatchDim<int>(cfg.dataset.dim, RunDimRunner{&cfg, &args}, &rc, &derr)) {
    SJS_LOG_ERROR("Unsupported dim=", cfg.dataset.dim, ". ", derr);
    return 2;
  }
  return rc;
}
