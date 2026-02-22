// apps/sjs_sweep.cpp
//
// HighDims parameter sweep harness.
//
// Extends the 2D sweep app with dim support:
//   - synthetic: can sweep dim via sweep.dim (list) or base.dataset.dim
//   - binary/csv: dim is inferred from files (R/S must match); sweep.dim is ignored
//
// Writes:
//   (1) sweep_raw.csv     : one row per run (repeat/seed)
//   (2) sweep_summary.csv : aggregated stats per (dim, dataset params, method, variant, t)
//
// Uses a minimal JSON subset parser (same spirit as 2D version) so no external deps.

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
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
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
      k == "j_star" || k == "enum_cap" || k == "write_samples" || k == "verify") return true;
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
  if (auto v = args.Get("j_star")) {
    u64 x = 0;
    if (!TryParseU64(*v, &x)) return fail("Invalid --j_star (expected u64).");
    cfg.run.j_star = x;
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

inline std::string FormatSci(double x, int prec = 3) {
  std::ostringstream oss;
  oss << std::scientific << std::setprecision(prec) << x;
  std::string s = oss.str();
  // Make it filename-friendly.
  for (char& c : s) {
    if (c == '+') c = 'p';
    if (c == '-') c = 'm';
    if (c == '.') c = '_';
  }
  return s;
}

inline std::string FirstLine(std::string_view s) {
  const size_t pos = s.find_first_of("\r\n");
  if (pos == std::string_view::npos) return std::string(s);
  return std::string(s.substr(0, pos));
}

inline bool ExistsNoThrow(const fs::path& p) noexcept {
  try {
    return fs::exists(p);
  } catch (...) {
    return false;
  }
}

inline std::string ResolvePathByProbingParents(std::string_view path_sv,
                                              const fs::path& sweep_dir,
                                              int max_parents = 8) {
  fs::path p{std::string(path_sv)};
  if (p.empty()) return std::string(path_sv);

  if (p.is_absolute()) return p.string();
  if (ExistsNoThrow(p)) return p.string();

  fs::path probe = sweep_dir;
  for (int i = 0; i < max_parents && !probe.empty(); ++i) {
    fs::path candidate = (probe / p).lexically_normal();
    if (ExistsNoThrow(candidate)) return candidate.string();
    fs::path parent = probe.parent_path();
    if (parent == probe) break;
    probe = parent;
  }
  return p.string();
}

inline void ResolveDatasetPathsFromSweep(const fs::path& sweep_json_path, sjs::Config* cfg) {
  if (!cfg) return;
  if (cfg->dataset.source == sjs::DataSource::Synthetic) return;

  fs::path sweep_abs;
  fs::path sweep_dir;
  try {
    sweep_abs = fs::absolute(sweep_json_path);
    sweep_dir = sweep_abs.has_parent_path() ? sweep_abs.parent_path() : fs::current_path();
  } catch (...) {
    sweep_dir = fs::current_path();
  }

  const std::string old_r = cfg->dataset.path_r;
  const std::string old_s = cfg->dataset.path_s;

  if (!cfg->dataset.path_r.empty()) {
    cfg->dataset.path_r = ResolvePathByProbingParents(cfg->dataset.path_r, sweep_dir);
  }
  if (!cfg->dataset.path_s.empty()) {
    cfg->dataset.path_s = ResolvePathByProbingParents(cfg->dataset.path_s, sweep_dir);
  }

  if (cfg->dataset.path_r != old_r) {
    SJS_LOG_INFO("Resolved path_r: '", old_r, "' -> '", cfg->dataset.path_r, "'");
  }
  if (cfg->dataset.path_s != old_s) {
    SJS_LOG_INFO("Resolved path_s: '", old_s, "' -> '", cfg->dataset.path_s, "'");
  }
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
    SetErr(err, "PeekBinaryDim: endian marker mismatch: " + path.string());
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
    if (c % 2 == 0) { *out_dim = static_cast<int>(c / 2); return true; }
    if ((c - 1) % 2 == 0) { *out_dim = static_cast<int>((c - 1) / 2); return true; }

    SetErr(err, "InferCsvDim: bad column count: " + std::to_string(static_cast<unsigned long long>(c)));
    return false;
  }
  SetErr(err, "InferCsvDim: no data rows found: " + path.string());
  return false;
}

inline bool DetermineFixedDimForNonSynthetic(sjs::Config* cfg, std::string* err) {
  if (!cfg) { SetErr(err, "DetermineFixedDimForNonSynthetic: cfg is null"); return false; }
  if (cfg->dataset.source == sjs::DataSource::Synthetic) return true;

  if (cfg->dataset.path_r.empty() || cfg->dataset.path_s.empty()) {
    SetErr(err, "Non-synthetic dataset requires path_r/path_s");
    return false;
  }
  const fs::path path_r = fs::path(cfg->dataset.path_r);
  const fs::path path_s = fs::path(cfg->dataset.path_s);

  int dr = 0, ds = 0;
  std::string local_err;

  if (cfg->dataset.source == sjs::DataSource::Binary) {
    if (!PeekBinaryDim(path_r, &dr, &local_err)) { SetErr(err, local_err); return false; }
    if (!PeekBinaryDim(path_s, &ds, &local_err)) { SetErr(err, local_err); return false; }
  } else if (cfg->dataset.source == sjs::DataSource::CSV) {
    char sep = ',';
    if (auto it = cfg->run.extra.find("csv_sep"); it != cfg->run.extra.end()) {
      if (!it->second.empty()) {
        if (it->second == "tab" || it->second == "\\t") sep = '\t';
        else sep = it->second[0];
      }
    }
    if (!InferCsvDim(path_r, sep, /*has_header=*/true, &dr, &local_err)) { SetErr(err, local_err); return false; }
    if (!InferCsvDim(path_s, sep, /*has_header=*/true, &ds, &local_err)) { SetErr(err, local_err); return false; }
  } else {
    SetErr(err, "Unsupported dataset_source");
    return false;
  }

  if (dr != ds) {
    SetErr(err, "R/S dim mismatch: " + std::to_string(dr) + " vs " + std::to_string(ds));
    return false;
  }
  cfg->dataset.dim = dr;
  return true;
}

template <int Dim>
bool LoadOrGenerateDataset(const sjs::Config& cfg,
                           sjs::Dataset<Dim, sjs::Scalar>* out,
                           sjs::synthetic::Report* gen_report,
                           std::string* err) {
  if (!out) { SetErr(err, "LoadOrGenerateDataset: out is null"); return false; }
  if (cfg.dataset.dim != Dim) { SetErr(err, "Config dim mismatch"); return false; }

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

    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(cfg.dataset.path_r, &R, &infoR, opt, &local_err)) {
      SetErr(err, "Failed reading binary R: " + local_err);
      return false;
    }
    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(cfg.dataset.path_s, &S, &infoS, opt, &local_err)) {
      SetErr(err, "Failed reading binary S: " + local_err);
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

    if (!sjs::csv::ReadBoxesSimple<Dim, sjs::Scalar>(cfg.dataset.path_r, &R, sep, /*has_header=*/true, &local_err)) {
      SetErr(err, "Failed reading CSV R: " + local_err);
      return false;
    }
    if (!sjs::csv::ReadBoxesSimple<Dim, sjs::Scalar>(cfg.dataset.path_s, &S, sep, /*has_header=*/true, &local_err)) {
      SetErr(err, "Failed reading CSV S: " + local_err);
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

inline std::vector<std::string> RawHeader() {
  return {
      "dataset",
      "dim",
      "generator",
      "alpha",
      "n_r",
      "n_s",
      "method",
      "variant",
      "t",
      "rep",
      "seed",
      "ok",
      "error",
      "wall_ms",
      "count_value",
      "count_exact",
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

inline std::vector<std::string> SummaryHeader() {
  return {
      "dataset",
      "dim",
      "generator",
      "alpha",
      "n_r",
      "n_s",
      "method",
      "variant",
      "t",
      "repeats",
      "ok_rate",
      "wall_mean_ms",
      "wall_stdev_ms",
      "wall_median_ms",
      "wall_p95_ms",
      "count_mean",
      "count_stdev",
      "exact_frac",
      "note",
  };
}

// --------------------------
// Minimal JSON parser (subset)
// --------------------------
namespace json {

enum class Type : sjs::u8 { Null, Bool, Number, String, Array, Object };

struct Value {
  Type type{Type::Null};
  bool b{false};
  double num{0.0};
  std::string str;
  std::vector<Value> arr;
  std::unordered_map<std::string, std::unique_ptr<Value>> obj;

  bool IsNull() const { return type == Type::Null; }
  bool IsBool() const { return type == Type::Bool; }
  bool IsNumber() const { return type == Type::Number; }
  bool IsString() const { return type == Type::String; }
  bool IsArray() const { return type == Type::Array; }
  bool IsObject() const { return type == Type::Object; }
};

class Parser {
 public:
  explicit Parser(std::string_view s) : s_(s) {}

  bool Parse(Value* out, std::string* err) {
    if (!out) { SetErr(err, "json::Parser: out is null"); return false; }
    SkipWs();
    if (!ParseValue(out, err)) return false;
    SkipWs();
    if (i_ != s_.size()) { SetErr(err, "json parse: trailing characters"); return false; }
    return true;
  }

 private:
  void SkipWs() {
    while (i_ < s_.size()) {
      char c = s_[i_];
      if (c == ' ' || c == '\n' || c == '\r' || c == '\t') ++i_;
      else break;
    }
  }

  bool Match(std::string_view kw) {
    if (s_.substr(i_, kw.size()) == kw) { i_ += kw.size(); return true; }
    return false;
  }

  bool ParseValue(Value* out, std::string* err) {
    SkipWs();
    if (i_ >= s_.size()) { SetErr(err, "json parse: unexpected end"); return false; }
    const char c = s_[i_];
    if (c == '{') return ParseObject(out, err);
    if (c == '[') return ParseArray(out, err);
    if (c == '"') return ParseString(out, err);
    if (c == 't') {
      if (!Match("true")) { SetErr(err, "json parse: expected true"); return false; }
      out->type = Type::Bool; out->b = true; return true;
    }
    if (c == 'f') {
      if (!Match("false")) { SetErr(err, "json parse: expected false"); return false; }
      out->type = Type::Bool; out->b = false; return true;
    }
    if (c == 'n') {
      if (!Match("null")) { SetErr(err, "json parse: expected null"); return false; }
      out->type = Type::Null; return true;
    }
    if (c == '-' || (c >= '0' && c <= '9')) return ParseNumber(out, err);

    SetErr(err, std::string("json parse: unexpected char '") + c + "'");
    return false;
  }

  bool ParseString(Value* out, std::string* err) {
    if (s_[i_] != '"') { SetErr(err, "json parse: expected '\"'"); return false; }
    ++i_;
    std::string buf;
    while (i_ < s_.size()) {
      const char c = s_[i_++];
      if (c == '"') {
        out->type = Type::String;
        out->str = std::move(buf);
        return true;
      }
      if (c == '\\') {
        if (i_ >= s_.size()) { SetErr(err, "json parse: bad escape"); return false; }
        const char e = s_[i_++];
        switch (e) {
          case '"': buf.push_back('"'); break;
          case '\\': buf.push_back('\\'); break;
          case '/': buf.push_back('/'); break;
          case 'b': buf.push_back('\b'); break;
          case 'f': buf.push_back('\f'); break;
          case 'n': buf.push_back('\n'); break;
          case 'r': buf.push_back('\r'); break;
          case 't': buf.push_back('\t'); break;
          default:
            SetErr(err, "json parse: unsupported escape"); return false;
        }
      } else {
        buf.push_back(c);
      }
    }
    SetErr(err, "json parse: unterminated string");
    return false;
  }

  bool ParseNumber(Value* out, std::string* err) {
    const size_t start = i_;
    if (s_[i_] == '-') ++i_;
    while (i_ < s_.size() && std::isdigit(static_cast<unsigned char>(s_[i_]))) ++i_;
    if (i_ < s_.size() && s_[i_] == '.') {
      ++i_;
      while (i_ < s_.size() && std::isdigit(static_cast<unsigned char>(s_[i_]))) ++i_;
    }
    if (i_ < s_.size() && (s_[i_] == 'e' || s_[i_] == 'E')) {
      ++i_;
      if (i_ < s_.size() && (s_[i_] == '+' || s_[i_] == '-')) ++i_;
      while (i_ < s_.size() && std::isdigit(static_cast<unsigned char>(s_[i_]))) ++i_;
    }

    const std::string num_str(s_.substr(start, i_ - start));
    try {
      out->type = Type::Number;
      out->num = std::stod(num_str);
      return true;
    } catch (...) {
      SetErr(err, "json parse: invalid number");
      return false;
    }
  }

  bool ParseArray(Value* out, std::string* err) {
    if (s_[i_] != '[') { SetErr(err, "json parse: expected '['"); return false; }
    ++i_;
    SkipWs();

    out->type = Type::Array;
    out->arr.clear();

    if (i_ < s_.size() && s_[i_] == ']') { ++i_; return true; }

    while (true) {
      Value v;
      if (!ParseValue(&v, err)) return false;
      out->arr.push_back(std::move(v));

      SkipWs();
      if (i_ >= s_.size()) { SetErr(err, "json parse: unterminated array"); return false; }
      const char c = s_[i_++];
      if (c == ']') break;
      if (c != ',') { SetErr(err, "json parse: expected ',' in array"); return false; }
      SkipWs();
    }
    return true;
  }

  bool ParseObject(Value* out, std::string* err) {
    if (s_[i_] != '{') { SetErr(err, "json parse: expected '{'"); return false; }
    ++i_;
    SkipWs();

    out->type = Type::Object;
    out->obj.clear();

    if (i_ < s_.size() && s_[i_] == '}') { ++i_; return true; }

    while (true) {
      Value key;
      if (!ParseString(&key, err)) return false;
      SkipWs();
      if (i_ >= s_.size() || s_[i_] != ':') { SetErr(err, "json parse: expected ':'"); return false; }
      ++i_;
      SkipWs();

      auto val = std::make_unique<Value>();
      if (!ParseValue(val.get(), err)) return false;
      out->obj.emplace(key.str, std::move(val));

      SkipWs();
      if (i_ >= s_.size()) { SetErr(err, "json parse: unterminated object"); return false; }
      const char c = s_[i_++];
      if (c == '}') break;
      if (c != ',') { SetErr(err, "json parse: expected ',' in object"); return false; }
      SkipWs();
    }
    return true;
  }

  std::string_view s_;
  size_t i_{0};
};

inline const Value* Get(const Value& obj, std::string_view key) {
  if (!obj.IsObject()) return nullptr;
  auto it = obj.obj.find(std::string(key));
  if (it == obj.obj.end()) return nullptr;
  return it->second.get();
}

inline bool GetString(const Value& v, std::string* out) {
  if (!out) return false;
  if (!v.IsString()) return false;
  *out = v.str;
  return true;
}

inline bool GetBool(const Value& v, bool* out) {
  if (!out) return false;
  if (!v.IsBool()) return false;
  *out = v.b;
  return true;
}

inline bool GetNumber(const Value& v, double* out) {
  if (!out) return false;
  if (!v.IsNumber()) return false;
  *out = v.num;
  return true;
}

inline bool GetU64(const Value& v, sjs::u64* out) {
  if (!out) return false;
  if (v.IsNumber()) {
    if (!(v.num >= 0.0)) return false;
    const long double x = static_cast<long double>(v.num);
    if (x > static_cast<long double>(std::numeric_limits<sjs::u64>::max())) return false;
    *out = static_cast<sjs::u64>(x + 0.5L);
    return true;
  }
  return false;
}

inline std::vector<sjs::u64> GetU64List(const Value* v, const std::vector<sjs::u64>& def) {
  if (!v || !v->IsArray()) return def;
  std::vector<sjs::u64> out;
  for (const auto& x : v->arr) {
    sjs::u64 u;
    if (GetU64(x, &u)) out.push_back(u);
  }
  return out.empty() ? def : out;
}

inline std::vector<double> GetDoubleList(const Value* v, const std::vector<double>& def) {
  if (!v || !v->IsArray()) return def;
  std::vector<double> out;
  for (const auto& x : v->arr) {
    double d;
    if (GetNumber(x, &d)) out.push_back(d);
  }
  return out.empty() ? def : out;
}

inline std::vector<std::string> GetStringList(const Value* v, const std::vector<std::string>& def) {
  if (!v || !v->IsArray()) return def;
  std::vector<std::string> out;
  for (const auto& x : v->arr) {
    if (x.IsString()) out.push_back(x.str);
  }
  return out.empty() ? def : out;
}

inline std::unordered_map<std::string, std::string> GetStringMap(const Value* v) {
  std::unordered_map<std::string, std::string> out;
  if (!v || !v->IsObject()) return out;
  for (const auto& kv : v->obj) {
    const Value& vv = *kv.second;
    if (vv.IsString()) out.emplace(kv.first, vv.str);
    else if (vv.IsNumber()) out.emplace(kv.first, std::to_string(vv.num));
    else if (vv.IsBool()) out.emplace(kv.first, vv.b ? "true" : "false");
  }
  return out;
}

}  // namespace json

// --------------------------
// Sweep spec
// --------------------------
struct SweepSpec {
  sjs::Config base;

  // sweep lists
  std::vector<int> dims;           // NEW: dim sweep (synthetic only)
  std::vector<double> alphas;
  // Equal-scale convenience sweep (N := |R| = |S|).
  // If sweep.n is provided, it overrides n_r/n_s lists.
  bool use_equal_n = false;
  std::vector<sjs::u64> n_list;
  std::vector<sjs::u64> n_r_list;
  std::vector<sjs::u64> n_s_list;
  std::vector<sjs::u64> gen_seeds;
  std::vector<sjs::u64> t_list;
  std::vector<sjs::Method> methods;
  std::vector<sjs::Variant> variants;
  std::vector<sjs::u64> seeds;     // optional explicit list (overrides repeats)

  // output
  std::string raw_file = "sweep_raw.csv";
  std::string summary_file = "sweep_summary.csv";
};

inline bool ApplyBaseConfigFromJson(const json::Value& root, sjs::Config* cfg, std::string* err) {
  if (!cfg) { SetErr(err, "ApplyBaseConfigFromJson: cfg is null"); return false; }

  const json::Value* base = json::Get(root, "base");
  if (!base || !base->IsObject()) return true;

  // dataset
  if (const json::Value* ds = json::Get(*base, "dataset"); ds && ds->IsObject()) {
    if (const json::Value* v = json::Get(*ds, "source")) {
      std::string s;
      if (json::GetString(*v, &s)) {
        sjs::DataSource src;
        if (sjs::ParseDataSource(s, &src)) cfg->dataset.source = src;
      }
    }
    if (const json::Value* v = json::Get(*ds, "name")) {
      std::string s;
      if (json::GetString(*v, &s)) cfg->dataset.name = s;
    }
    if (const json::Value* v = json::Get(*ds, "dim")) {
      sjs::u64 d;
      if (json::GetU64(*v, &d) && d <= static_cast<sjs::u64>(std::numeric_limits<sjs::i32>::max())) {
        cfg->dataset.dim = static_cast<sjs::i32>(d);
      }
    }
    if (const json::Value* v = json::Get(*ds, "path_r")) {
      std::string s; if (json::GetString(*v, &s)) cfg->dataset.path_r = s;
    }
    if (const json::Value* v = json::Get(*ds, "path_s")) {
      std::string s; if (json::GetString(*v, &s)) cfg->dataset.path_s = s;
    }

    if (const json::Value* syn = json::Get(*ds, "synthetic"); syn && syn->IsObject()) {
      if (const json::Value* v = json::Get(*syn, "generator")) {
        std::string s; if (json::GetString(*v, &s)) cfg->dataset.synthetic.generator = s;
      }
      if (const json::Value* v = json::Get(*syn, "n_r")) {
        sjs::u64 x; if (json::GetU64(*v, &x)) cfg->dataset.synthetic.n_r = x;
      }
      if (const json::Value* v = json::Get(*syn, "n_s")) {
        sjs::u64 x; if (json::GetU64(*v, &x)) cfg->dataset.synthetic.n_s = x;
      }
      if (const json::Value* v = json::Get(*syn, "alpha")) {
        double x; if (json::GetNumber(*v, &x)) cfg->dataset.synthetic.alpha = x;
      }
      if (const json::Value* v = json::Get(*syn, "seed")) {
        sjs::u64 x; if (json::GetU64(*v, &x)) cfg->dataset.synthetic.seed = x;
      }
      if (const json::Value* v = json::Get(*syn, "params")) {
        cfg->dataset.synthetic.extra = json::GetStringMap(v);
      }
    }
  }

  // run
  if (const json::Value* run = json::Get(*base, "run"); run && run->IsObject()) {
    if (const json::Value* v = json::Get(*run, "method")) {
      std::string s;
      if (json::GetString(*v, &s)) { sjs::Method m; if (sjs::ParseMethod(s, &m)) cfg->run.method = m; }
    }
    if (const json::Value* v = json::Get(*run, "variant")) {
      std::string s;
      if (json::GetString(*v, &s)) { sjs::Variant vv; if (sjs::ParseVariant(s, &vv)) cfg->run.variant = vv; }
    }
    if (const json::Value* v = json::Get(*run, "t")) { sjs::u64 x; if (json::GetU64(*v, &x)) cfg->run.t = x; }
    if (const json::Value* v = json::Get(*run, "seed")) { sjs::u64 x; if (json::GetU64(*v, &x)) cfg->run.seed = x; }
    if (const json::Value* v = json::Get(*run, "repeats")) { sjs::u64 x; if (json::GetU64(*v, &x)) cfg->run.repeats = x; }
    if (const json::Value* v = json::Get(*run, "enum_cap")) { sjs::u64 x; if (json::GetU64(*v, &x)) cfg->run.enum_cap = x; }
    if (const json::Value* v = json::Get(*run, "j_star")) { sjs::u64 x; if (json::GetU64(*v, &x)) cfg->run.j_star = x; }
    if (const json::Value* v = json::Get(*run, "write_samples")) { bool b; if (json::GetBool(*v, &b)) cfg->run.write_samples = b; }
    if (const json::Value* v = json::Get(*run, "extra")) { cfg->run.extra = json::GetStringMap(v); }
  }

  // output
  if (const json::Value* out = json::Get(*base, "output"); out && out->IsObject()) {
    if (const json::Value* v = json::Get(*out, "out_dir")) {
      std::string s; if (json::GetString(*v, &s)) cfg->output.out_dir = s;
    }
    if (const json::Value* v = json::Get(*out, "run_tag")) {
      std::string s; if (json::GetString(*v, &s)) cfg->output.run_tag = s;
    }
  }

  // logging
  if (const json::Value* log = json::Get(*base, "logging"); log && log->IsObject()) {
    if (const json::Value* v = json::Get(*log, "level")) {
      std::string s;
      if (json::GetString(*v, &s)) { sjs::LogLevel lvl; if (sjs::detail::ParseLogLevel(s, &lvl)) cfg->logging.level = lvl; }
    }
    if (const json::Value* v = json::Get(*log, "with_timestamp")) { bool b; if (json::GetBool(*v, &b)) cfg->logging.with_timestamp = b; }
    if (const json::Value* v = json::Get(*log, "with_thread_id")) { bool b; if (json::GetBool(*v, &b)) cfg->logging.with_thread_id = b; }
  }

  // sys
  if (const json::Value* sys = json::Get(*base, "sys"); sys && sys->IsObject()) {
    if (const json::Value* v = json::Get(*sys, "threads")) { sjs::u64 x; if (json::GetU64(*v, &x)) cfg->sys.threads = static_cast<sjs::i32>(x); }
  }

  return true;
}

inline bool LoadSweepSpecJson(const std::string& path, SweepSpec* spec, std::string* err) {
  if (!spec) { SetErr(err, "LoadSweepSpecJson: spec is null"); return false; }

  std::ifstream in(path);
  if (!in) { SetErr(err, "Cannot open sweep config: " + path); return false; }
  std::stringstream buffer;
  buffer << in.rdbuf();
  const std::string text = buffer.str();

  json::Value root;
  std::string jerr;
  if (!json::Parser(std::string_view(text)).Parse(&root, &jerr)) {
    SetErr(err, "JSON parse failed: " + jerr);
    return false;
  }
  if (!root.IsObject()) {
    SetErr(err, "Sweep JSON must be a top-level object");
    return false;
  }

  if (!ApplyBaseConfigFromJson(root, &spec->base, err)) return false;

  // sweep lists
  const json::Value* sweep = json::Get(root, "sweep");
  if (sweep && sweep->IsObject()) {
    // NEW: dim sweep
    {
      std::vector<sjs::u64> dims_u64 = json::GetU64List(json::Get(*sweep, "dim"), {});
      if (!dims_u64.empty()) {
        spec->dims.clear();
        for (sjs::u64 d : dims_u64) {
          if (d >= 2 && d <= static_cast<sjs::u64>(sjs::kMaxSupportedDim)) {
            spec->dims.push_back(static_cast<int>(d));
          }
        }
      }
    }

    spec->alphas    = json::GetDoubleList(json::Get(*sweep, "alpha"), spec->alphas);
    spec->n_r_list  = json::GetU64List(json::Get(*sweep, "n_r"), spec->n_r_list);
    spec->n_s_list  = json::GetU64List(json::Get(*sweep, "n_s"), spec->n_s_list);

    // Convenience: sweep "n" means equal sizes on both sides (N := |R| = |S|).
    // If present, it overrides n_r/n_s sweep lists.
    {
      std::vector<sjs::u64> ns = json::GetU64List(json::Get(*sweep, "n"), {});
      if (!ns.empty()) {
        spec->use_equal_n = true;
        spec->n_list = std::move(ns);
        spec->n_r_list = spec->n_list;
        spec->n_s_list = spec->n_list;
      }
    }
    spec->gen_seeds = json::GetU64List(json::Get(*sweep, "gen_seed"), spec->gen_seeds);
    spec->t_list    = json::GetU64List(json::Get(*sweep, "t"), spec->t_list);

    if (const json::Value* m = json::Get(*sweep, "method")) {
      std::vector<std::string> ms = json::GetStringList(m, {});
      if (!ms.empty()) {
        spec->methods.clear();
        for (const auto& s : ms) {
          sjs::Method mm;
          if (sjs::ParseMethod(s, &mm) && mm != sjs::Method::Unknown) spec->methods.push_back(mm);
        }
      }
    }
    if (const json::Value* v = json::Get(*sweep, "variant")) {
      std::vector<std::string> vs = json::GetStringList(v, {});
      if (!vs.empty()) {
        spec->variants.clear();
        for (const auto& s : vs) {
          sjs::Variant vv;
          if (sjs::ParseVariant(s, &vv)) spec->variants.push_back(vv);
        }
      }
    }

    spec->seeds = json::GetU64List(json::Get(*sweep, "seed"), spec->seeds);

    if (const json::Value* r = json::Get(*sweep, "repeats")) {
      sjs::u64 x;
      if (json::GetU64(*r, &x)) spec->base.run.repeats = x;
    }
  }

  // files
  if (const json::Value* files = json::Get(root, "files"); files && files->IsObject()) {
    if (const json::Value* v = json::Get(*files, "raw")) {
      std::string s; if (json::GetString(*v, &s)) spec->raw_file = s;
    }
    if (const json::Value* v = json::Get(*files, "summary")) {
      std::string s; if (json::GetString(*v, &s)) spec->summary_file = s;
    }
  }

  // Fill defaults if still empty
  if (spec->dims.empty()) spec->dims = {static_cast<int>(spec->base.dataset.dim)};
  if (spec->alphas.empty()) spec->alphas = {spec->base.dataset.synthetic.alpha};
  if (spec->n_r_list.empty()) spec->n_r_list = {spec->base.dataset.synthetic.n_r};
  if (spec->n_s_list.empty()) spec->n_s_list = {spec->base.dataset.synthetic.n_s};
  if (spec->gen_seeds.empty()) spec->gen_seeds = {spec->base.dataset.synthetic.seed};
  if (spec->t_list.empty()) spec->t_list = {spec->base.run.t};
  if (spec->methods.empty()) spec->methods = {spec->base.run.method};
  if (spec->variants.empty()) spec->variants = {spec->base.run.variant};

  return true;
}

inline void PrintUsage() {
  std::cerr
      << "sjs_sweep: run a parameter sweep from a JSON config (HighDims)\n\n"
      << "Required:\n"
      << "  --config=<path/to/sweep.json>    (or provide JSON path as a positional arg)\n\n"
      << "Common overrides (optional):\n"
      << "  --out_dir=<dir>                 (overrides base.output.out_dir)\n"
      << "  --raw_file=<path>               (overrides files.raw)\n"
      << "  --summary_file=<path>           (overrides files.summary)\n\n"
      << "Notes:\n"
      << "  - For synthetic sweeps, you can provide sweep.dim: [2,3,4,...].\n"
      << "  - For synthetic sweeps, sweep.n: [N1,N2,...] sets n_r=n_s=N (diagonal).\n"
      << "  - For binary/csv datasets, dim is inferred from the files.\n";
}

template <int Dim>
int RunSweepForDim(const sjs::apps::SweepSpec& spec,
                   sjs::csv::Writer& raw_writer,
                   sjs::csv::Writer& sum_writer) {
  using sjs::u64;
  using sjs::usize;

  std::string err;

  // For binary/csv, load dataset once (no dataset-param sweep).
  sjs::Dataset<Dim, sjs::Scalar> fixed_ds;
  sjs::synthetic::Report fixed_gen_report;
  bool has_fixed_ds = false;

  if (spec.base.dataset.source != sjs::DataSource::Synthetic) {
    std::string local_err;
    if (!LoadOrGenerateDataset<Dim>(spec.base, &fixed_ds, &fixed_gen_report, &local_err)) {
      SJS_LOG_ERROR(local_err);
      return 6;
    }
    has_fixed_ds = true;
    SJS_LOG_INFO("Loaded fixed dataset: ", fixed_ds.name,
                 " Dim=", Dim,
                 " R=", static_cast<unsigned long long>(fixed_ds.R.Size()),
                 " S=", static_cast<unsigned long long>(fixed_ds.S.Size()));
  }

  u64 total_runs = 0;
  u64 ok_runs = 0;

  const bool sweep_dataset = (spec.base.dataset.source == sjs::DataSource::Synthetic);

  const std::vector<u64> n_r_list = sweep_dataset ? spec.n_r_list : std::vector<u64>{0};
  const std::vector<u64> n_s_list = sweep_dataset ? spec.n_s_list : std::vector<u64>{0};
  const std::vector<double> alpha_list = sweep_dataset ? spec.alphas : std::vector<double>{0.0};
  const std::vector<u64> gen_seeds = sweep_dataset ? spec.gen_seeds : std::vector<u64>{0};

  // Dataset size combos:
  //   - default: cartesian product of n_r_list x n_s_list
  //   - if sweep.n was provided: diagonal pairs (n,n)
  std::vector<std::pair<u64, u64>> size_pairs;
  if (sweep_dataset && spec.use_equal_n && !spec.n_list.empty()) {
    size_pairs.reserve(spec.n_list.size());
    for (u64 n : spec.n_list) size_pairs.emplace_back(n, n);
  } else {
    size_pairs.reserve(n_r_list.size() * n_s_list.size());
    for (u64 n_r : n_r_list) {
      for (u64 n_s : n_s_list) {
        size_pairs.emplace_back(n_r, n_s);
      }
    }
  }

  for (const auto& ns_pair : size_pairs) {
    const u64 n_r = ns_pair.first;
    const u64 n_s = ns_pair.second;
    for (double alpha : alpha_list) {
      for (u64 gen_seed : gen_seeds) {
          // Dataset for this combo.
          sjs::Dataset<Dim, sjs::Scalar> ds;
          sjs::synthetic::Report gen_report;
          const sjs::synthetic::Report* gen_report_ptr = nullptr;

          sjs::Config cfg_ds = spec.base;
          cfg_ds.dataset.dim = Dim;

          if (sweep_dataset) {
            cfg_ds.dataset.synthetic.n_r = n_r;
            cfg_ds.dataset.synthetic.n_s = n_s;
            cfg_ds.dataset.synthetic.alpha = alpha;
            cfg_ds.dataset.synthetic.seed = gen_seed;

            // Auto-name dataset.
            {
              std::ostringstream nm;
              nm << (spec.base.dataset.name.empty() ? "synthetic" : spec.base.dataset.name)
                 << "__d" << Dim
                 << "__nr" << static_cast<unsigned long long>(n_r)
                 << "__ns" << static_cast<unsigned long long>(n_s)
                 << "__a" << sjs::apps::FormatSci(alpha, 3)
                 << "__g" << static_cast<unsigned long long>(gen_seed);
              cfg_ds.dataset.name = sjs::apps::SanitizeFilename(nm.str());
            }

            std::string local_err;
            if (!LoadOrGenerateDataset<Dim>(cfg_ds, &ds, &gen_report, &local_err)) {
              SJS_LOG_ERROR("Dataset generation failed: ", local_err);
              return 6;
            }
            gen_report_ptr = &gen_report;
          } else {
            (void)n_r; (void)n_s; (void)alpha; (void)gen_seed;
            if (!has_fixed_ds) {
              SJS_LOG_ERROR("Internal error: expected fixed dataset to be loaded.");
              return 6;
            }
            ds = fixed_ds;
            gen_report_ptr = nullptr;
          }

          for (u64 t : spec.t_list) {
            for (sjs::Method method : spec.methods) {
              for (sjs::Variant variant : spec.variants) {
                // cfg for this group.
                sjs::Config cfg = cfg_ds;
                cfg.run.t = t;
                cfg.run.method = method;
                cfg.run.variant = variant;

                // seeds
                std::vector<u64> seeds;
                if (!spec.seeds.empty()) {
                  seeds = spec.seeds;
                } else {
                  seeds.resize(static_cast<usize>(cfg.run.repeats));
                  for (u64 rep = 0; rep < cfg.run.repeats; ++rep) {
                    seeds[static_cast<usize>(rep)] = cfg.run.seed + rep;
                  }
                }

                // baseline
                std::string b_err;
                auto baseline = sjs::baselines::CreateBaseline<Dim, sjs::Scalar>(method, variant, &b_err);
                if (!baseline) {
                  SJS_LOG_ERROR("CreateBaseline failed (skipping this combo): ", b_err);

                  const double nan = std::numeric_limits<double>::quiet_NaN();
                  const std::string err_one_line = sjs::apps::FirstLine(b_err);
                  const std::string note = "unsupported baseline (skipped)";

                  for (usize i = 0; i < seeds.size(); ++i) {
                    const u64 seed = seeds[i];
                    total_runs++;

                    raw_writer.WriteRowV(
                        ds.name,
                        Dim,
                        (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.generator : ""),
                        (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.alpha : nan),
                        static_cast<unsigned long long>(ds.R.Size()),
                        static_cast<unsigned long long>(ds.S.Size()),
                        sjs::ToString(method),
                        sjs::ToString(variant),
                        static_cast<unsigned long long>(t),
                        static_cast<unsigned long long>(i),
                        static_cast<unsigned long long>(seed),
                        0,
                        err_one_line,
                        nan,
                        nan,
                        0,
                        0,
                        0,
                        static_cast<unsigned long long>(cfg.run.enum_cap),
                        "",
                        0ULL,
                        "{}",
                        note,
                        cfg.ToJsonLite(),
                        (gen_report_ptr ? gen_report_ptr->ToJsonLite() : "{}"));
                  }

                  sum_writer.WriteRowV(
                      ds.name,
                      Dim,
                      (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.generator : ""),
                      (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.alpha : nan),
                      static_cast<unsigned long long>(ds.R.Size()),
                      static_cast<unsigned long long>(ds.S.Size()),
                      sjs::ToString(method),
                      sjs::ToString(variant),
                      static_cast<unsigned long long>(t),
                      static_cast<unsigned long long>(seeds.size()),
                      0.0,   // ok_rate
                      -1.0,  // wall_mean_ms
                      -1.0,  // wall_stdev_ms
                      -1.0,  // wall_median_ms
                      -1.0,  // wall_p95_ms
                      -1.0,  // count_mean
                      -1.0,  // count_stdev
                      0.0,   // exact_frac
                      ("unsupported baseline: " + err_one_line));
                  continue;
                }

                // stats
                std::vector<double> wall_ms;
                wall_ms.reserve(seeds.size());
                std::vector<double> count_vals;
                count_vals.reserve(seeds.size());
                u64 exact_count = 0;
                u64 ok_in_group = 0;

                for (usize i = 0; i < seeds.size(); ++i) {
                  const u64 seed = seeds[i];
                  sjs::baselines::RunReport rep_out;
                  std::string local_err;

                  sjs::Stopwatch sw;
                  bool ok = false;
                  switch (variant) {
                    case sjs::Variant::Sampling:
                      ok = sjs::baselines::RunSamplingOnce<Dim, sjs::Scalar>(
                          baseline.get(), ds, cfg, seed, &rep_out, &local_err);
                      break;
                    case sjs::Variant::EnumSampling:
                      ok = sjs::baselines::RunEnumSamplingOnce<Dim, sjs::Scalar>(
                          baseline.get(), ds, cfg, seed, &rep_out, &local_err);
                      break;
                    case sjs::Variant::Adaptive:
                      ok = sjs::baselines::RunAdaptiveOnce<Dim, sjs::Scalar>(
                          baseline.get(), ds, cfg, seed, &rep_out, &local_err);
                      break;
                  }

                  const double wms = sw.ElapsedMillis();
                  wall_ms.push_back(wms);

                  const double count_value = ok ? static_cast<double>(rep_out.count.value)
                                                : std::numeric_limits<double>::quiet_NaN();
                  count_vals.push_back(count_value);

                  total_runs++;
                  if (!ok) {
                    rep_out.ok = false;
                    rep_out.error = local_err;
                  } else {
                    ok_runs++;
                    ok_in_group++;
                    if (rep_out.count.exact) exact_count++;
                  }

                  raw_writer.WriteRowV(
                      ds.name,
                      Dim,
                      (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.generator : ""),
                      (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.alpha
                                                                        : std::numeric_limits<double>::quiet_NaN()),
                      static_cast<unsigned long long>(ds.R.Size()),
                      static_cast<unsigned long long>(ds.S.Size()),
                      sjs::ToString(method),
                      sjs::ToString(variant),
                      static_cast<unsigned long long>(t),
                      static_cast<unsigned long long>(i),
                      static_cast<unsigned long long>(seed),
                      (rep_out.ok ? 1 : 0),
                      sjs::apps::FirstLine(rep_out.error),
                      wms,
                      count_value,
                      ((rep_out.ok && rep_out.count.exact) ? 1 : 0),
                      (rep_out.used_enumeration ? 1 : 0),
                      (rep_out.enumeration_truncated ? 1 : 0),
                      static_cast<unsigned long long>(rep_out.enumeration_cap),
                      rep_out.adaptive_branch,
                      static_cast<unsigned long long>(rep_out.adaptive_pilot_pairs),
                      rep_out.phases.ToJsonMillis(),
                      rep_out.note,
                      cfg.ToJsonLite(),
                      (gen_report_ptr ? gen_report_ptr->ToJsonLite() : "{}"));
                }

                const sjs::Summary wall_sum = sjs::Summarize(wall_ms);

                std::vector<double> count_vals_finite;
                count_vals_finite.reserve(count_vals.size());
                for (double v : count_vals) if (std::isfinite(v)) count_vals_finite.push_back(v);

                double cnt_mean = -1.0;
                double cnt_stdev = -1.0;
                if (!count_vals_finite.empty()) {
                  const sjs::Summary cnt_sum = sjs::Summarize(count_vals_finite);
                  cnt_mean = cnt_sum.mean;
                  cnt_stdev = cnt_sum.stdev;
                }

                const double ok_rate = (seeds.empty() ? 0.0 : static_cast<double>(ok_in_group) /
                                                          static_cast<double>(seeds.size()));
                const double exact_frac = (seeds.empty() ? 0.0 : static_cast<double>(exact_count) /
                                                          static_cast<double>(seeds.size()));

                sum_writer.WriteRowV(
                    ds.name,
                    Dim,
                    (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.generator : ""),
                    (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.alpha
                                                                      : std::numeric_limits<double>::quiet_NaN()),
                    static_cast<unsigned long long>(ds.R.Size()),
                    static_cast<unsigned long long>(ds.S.Size()),
                    sjs::ToString(method),
                    sjs::ToString(variant),
                    static_cast<unsigned long long>(t),
                    static_cast<unsigned long long>(seeds.size()),
                    ok_rate,
                    wall_sum.mean,
                    wall_sum.stdev,
                    wall_sum.median,
                    wall_sum.p95,
                    cnt_mean,
                    cnt_stdev,
                    exact_frac,
                    "");
              }
            }
      }
    }
  }

  }  // for (ns_pair : size_pairs)
  SJS_LOG_INFO("Sweep finished (Dim=", Dim, "). total_runs=", static_cast<unsigned long long>(total_runs),
               " ok_runs=", static_cast<unsigned long long>(ok_runs));
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

struct SweepDimRunner {
  const sjs::apps::SweepSpec* spec = nullptr;
  sjs::csv::Writer* raw_writer = nullptr;
  sjs::csv::Writer* sum_writer = nullptr;

  template <int Dim>
  int operator()() const {
    return sjs::apps::RunSweepForDim<Dim>(*spec, *raw_writer, *sum_writer);
  }
};

}  // namespace

int main(int argc, char** argv) {
  using sjs::u64;

  sjs::ArgMap args = sjs::ArgMap::FromArgv(argc, argv);
  if (sjs::apps::IsHelpRequested(args)) {
    sjs::apps::PrintUsage();
    return 0;
  }

  // Locate sweep JSON path.
  std::string sweep_path;
  if (auto v = args.Get("config")) {
    sweep_path = std::string(*v);
  } else if (!args.Positional().empty()) {
    sweep_path = args.Positional().front();
  }
  if (sweep_path.empty()) {
    SJS_LOG_ERROR("Missing --config=<sweep.json>.");
    sjs::apps::PrintUsage();
    return 2;
  }

  // Start from CLI config, then JSON overrides base.
  sjs::Config cfg0;

  std::string cfg_err;
  if (!sjs::ParseConfigFromArgs(args, &cfg0, &cfg_err)) {
    SJS_LOG_ERROR("Config parse failed: ", cfg_err);
    sjs::apps::PrintUsage();
    return 2;
  }

  sjs::apps::SweepSpec spec;
  spec.base = cfg0;

  std::string err;
  if (!sjs::apps::LoadSweepSpecJson(sweep_path, &spec, &err)) {
    SJS_LOG_ERROR(err);
    return 3;
  }

  // CLI overrides
  if (auto v = args.Get("out_dir")) spec.base.output.out_dir = std::string(*v);
  if (auto v = args.Get("raw_file")) spec.raw_file = std::string(*v);
  if (auto v = args.Get("summary_file")) spec.summary_file = std::string(*v);

  // Resolve dataset paths relative to sweep JSON for binary/csv portability.
  sjs::apps::ResolveDatasetPathsFromSweep(fs::path(sweep_path), &spec.base);

  // For non-synthetic datasets: infer fixed dim from files and override base/sweep dims.
  if (spec.base.dataset.source != sjs::DataSource::Synthetic) {
    if (!sjs::apps::DetermineFixedDimForNonSynthetic(&spec.base, &err)) {
      SJS_LOG_ERROR("Failed to infer dim for non-synthetic dataset: ", err);
      return 4;
    }
    spec.dims = {static_cast<int>(spec.base.dataset.dim)};
    SJS_LOG_INFO("Inferred dim from files: ", spec.base.dataset.dim);
  }

  if (!spec.base.Validate(&err)) {
    SJS_LOG_ERROR("Base config validation failed: ", err);
    return 4;
  }

  sjs::Logger::Instance().SetConfig(spec.base.logging);

  // Output files.
  const fs::path out_dir(spec.base.output.out_dir);
  if (!sjs::apps::EnsureDir(out_dir, &err)) {
    SJS_LOG_ERROR("Cannot create out_dir: ", err);
    return 5;
  }

  const std::string raw_path = (out_dir / spec.raw_file).string();
  const std::string summary_path = (out_dir / spec.summary_file).string();

  if (!sjs::apps::EnsureDir(fs::path(raw_path).parent_path(), &err) ||
      !sjs::apps::EnsureDir(fs::path(summary_path).parent_path(), &err)) {
    SJS_LOG_ERROR("Cannot create output directories: ", err);
    return 5;
  }

  sjs::csv::Writer raw_writer(raw_path, sjs::csv::Dialect{','}, &err);
  if (!raw_writer.Ok()) {
    SJS_LOG_ERROR("Cannot open raw CSV: ", raw_path, " err=", err);
    return 5;
  }
  sjs::csv::Writer sum_writer(summary_path, sjs::csv::Dialect{','}, &err);
  if (!sum_writer.Ok()) {
    SJS_LOG_ERROR("Cannot open summary CSV: ", summary_path, " err=", err);
    return 5;
  }

  raw_writer.WriteHeader(sjs::apps::RawHeader(), &err);
  sum_writer.WriteHeader(sjs::apps::SummaryHeader(), &err);

  // Run per dim (synthetic may sweep dims; non-synthetic dims is singleton).
  for (int dim : spec.dims) {
    sjs::Config base_dim_cfg = spec.base;
    base_dim_cfg.dataset.dim = dim;

    // Validate base for this dim (synthetic only matters; non-synthetic already fixed).
    if (!base_dim_cfg.Validate(&err)) {
      SJS_LOG_ERROR("Base config validation failed for dim=", dim, ": ", err);
      return 4;
    }

    sjs::apps::SweepSpec spec_dim = spec;
    spec_dim.base = base_dim_cfg;
    int rc = 0;
    std::string derr;
    if (!sjs::dispatch::DispatchDim<int>(dim, SweepDimRunner{&spec_dim, &raw_writer, &sum_writer}, &rc, &derr)) {
      SJS_LOG_ERROR("Unsupported dim=", dim, ". ", derr);
      return 4;
    }
    if (rc != 0) return rc;
  }

  SJS_LOG_INFO("Raw CSV: ", raw_path);
  SJS_LOG_INFO("Summary CSV: ", summary_path);

  // Copy sweep config next to outputs for reproducibility.
  try {
    const fs::path dst = out_dir / "sweep_config.json";
    std::error_code ec;
    fs::copy_file(fs::path(sweep_path), dst, fs::copy_options::overwrite_existing, ec);
    if (!ec) SJS_LOG_INFO("Saved sweep config to: ", dst.string());
  } catch (...) {
    // ignore
  }

  return 0;
}
