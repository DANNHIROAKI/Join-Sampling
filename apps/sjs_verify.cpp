// apps/sjs_verify.cpp
//
// HighDims correctness + sampling-quality verification.
//
// Steps:
//   1) infer/resolve dim
//   2) load/generate a small dataset
//   3) compute oracle join size |J| via brute-force
//   4) run one baseline to get count + samples
//   5) if universe small enough, evaluate sample uniformity on join pairs
//
// Dim handling:
//   - synthetic: uses --dim (default 2)
//   - binary: reads dim from binary header
//   - csv: infers dim from first data row
//
// Requires:
//   - sjs/dispatch/dim_dispatch.h : sjs::dispatch::DispatchDim<R>(dim, fn, out, err)
//   - sjs/baselines/baseline_factory.h : CreateBaseline<Dim>(method, variant, err)

#include "sjs/baselines/baseline_factory.h"

#include "sjs/baselines/runners/adaptive_runner.h"
#include "sjs/baselines/runners/enum_sampling_runner.h"
#include "sjs/baselines/runners/sampling_runner.h"

#include "sjs/core/logging.h"
#include "sjs/core/types.h"

#include "sjs/dispatch/dim_dispatch.h"

#include "sjs/data/synthetic/generator_factory.h"
#include "sjs/io/binary_io.h"
#include "sjs/io/csv_io.h"
#include "sjs/io/dataset.h"

#include "sjs/join/join_oracle.h"
#include "sjs/sampling/sample_quality.h"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <limits>
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

inline bool ParseU64(const std::optional<std::string_view>& v, sjs::u64* out) {
  if (!v || !out) return false;
  try {
    std::size_t idx = 0;
    const std::string sv(*v);
    unsigned long long x = std::stoull(sv, &idx, 10);
    if (idx != sv.size()) return false;
    *out = static_cast<sjs::u64>(x);
    return true;
  } catch (...) {
    return false;
  }
}

inline bool PeekBinaryDim(const fs::path& path, int* out_dim, std::string* err) {
  if (!out_dim) { SetErr(err, "PeekBinaryDim: out_dim is null"); return false; }
  std::ifstream in(path, std::ios::binary);
  if (!in) { SetErr(err, "PeekBinaryDim: cannot open file: " + path.string()); return false; }
  sjs::binary::FileHeader h{};
  in.read(reinterpret_cast<char*>(&h), static_cast<std::streamsize>(sizeof(h)));
  if (!in) { SetErr(err, "PeekBinaryDim: failed to read header"); return false; }
  if (std::memcmp(h.magic, sjs::binary::kMagic, sizeof(sjs::binary::kMagic)) != 0) {
    SetErr(err, "PeekBinaryDim: bad magic"); return false;
  }
  if (h.endian != sjs::binary::kEndianMarker) {
    SetErr(err, "PeekBinaryDim: endian mismatch"); return false;
  }
  if (h.dim == 0 || h.dim > static_cast<sjs::u32>(sjs::kMaxSupportedDim)) {
    SetErr(err, "PeekBinaryDim: invalid dim"); return false;
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
    if (pos == std::string_view::npos) { out.push_back(line.substr(start)); break; }
    out.push_back(line.substr(start, pos - start));
    start = pos + 1;
  }
  return out;
}

inline bool InferCsvDim(const fs::path& path, char sep, bool has_header, int* out_dim, std::string* err) {
  if (!out_dim) { SetErr(err, "InferCsvDim: out_dim is null"); return false; }
  std::ifstream in(path);
  if (!in) { SetErr(err, "InferCsvDim: cannot open file: " + path.string()); return false; }
  std::string line;
  bool skipped_header = false;
  while (std::getline(in, line)) {
    std::string_view sv = Trim(std::string_view(line));
    if (sv.empty()) continue;
    if (!sv.empty() && sv.front() == '#') continue;
    if (has_header && !skipped_header) { skipped_header = true; continue; }

    const auto cols = SplitLine(sv, sep);
    const sjs::usize c = cols.size();
    if (c % 2 == 0) { *out_dim = static_cast<int>(c / 2); return true; }
    if ((c - 1) % 2 == 0) { *out_dim = static_cast<int>((c - 1) / 2); return true; }

    SetErr(err, "InferCsvDim: bad column count");
    return false;
  }
  SetErr(err, "InferCsvDim: no data rows found");
  return false;
}

inline bool DetermineDatasetDim(sjs::Config* cfg, const sjs::ArgMap& args, std::string* err) {
  if (!cfg) { SetErr(err, "DetermineDatasetDim: cfg is null"); return false; }

  if (cfg->dataset.source == sjs::DataSource::Synthetic) {
    if (cfg->dataset.dim <= 0) { SetErr(err, "Synthetic requires --dim > 0"); return false; }
    return true;
  }

  if (cfg->dataset.path_r.empty() || cfg->dataset.path_s.empty()) {
    SetErr(err, "Non-synthetic requires --path_r and --path_s");
    return false;
  }

  const fs::path path_r = ResolvePathByWalkingUp(fs::path(cfg->dataset.path_r));
  const fs::path path_s = ResolvePathByWalkingUp(fs::path(cfg->dataset.path_s));

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

  if (dr != ds) { SetErr(err, "R/S dim mismatch"); return false; }

  if (args.Has("dim") && cfg->dataset.dim != dr) {
    SetErr(err, "User --dim mismatch with inferred dim");
    return false;
  }
  cfg->dataset.dim = dr;
  cfg->dataset.path_r = path_r.string();
  cfg->dataset.path_s = path_s.string();
  return true;
}

template <int Dim>
bool LoadOrGenerateDataset(const sjs::Config& cfg,
                           sjs::Dataset<Dim, sjs::Scalar>* out,
                           sjs::synthetic::Report* gen_report,
                           std::string* err) {
  if (!out) { SetErr(err, "LoadOrGenerateDataset: out is null"); return false; }

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
    sjs::binary::BinaryReadOptions opt;
    opt.generate_ids_if_missing = true;
    opt.drop_empty = false;

    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(cfg.dataset.path_r, &R, nullptr, opt, &local_err)) {
      SetErr(err, "Failed reading binary R: " + local_err);
      return false;
    }
    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(cfg.dataset.path_s, &S, nullptr, opt, &local_err)) {
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

inline void PrintUsage() {
  std::cerr
      << "sjs_verify: small-scale correctness & sampling-quality checks (HighDims)\n\n"
      << "This app understands the sjs_run flags, plus oracle controls:\n"
      << "  --oracle_max_checks=<u64>     (default 50000000)  limit on |R|*|S|\n"
      << "  --oracle_collect_limit=<u64>  (default 1000000)   max |J| to fully collect pairs\n"
      << "  --oracle_cap=<u64>            (default 0)         collect at most this many pairs\n"
      << "\nDim rules:\n"
      << "  - synthetic: use --dim\n"
      << "  - binary: dim from file header\n"
      << "  - csv: dim inferred from first data row\n";
}

inline double RelErr(double est, double truth) {
  if (truth == 0.0) return std::numeric_limits<double>::quiet_NaN();
  return (est - truth) / truth;
}

template <int Dim>
int RunVerifyForDim(const sjs::Config& cfg) {
  using sjs::u64;
  using sjs::usize;

  std::string err;

  u64 oracle_max_checks = 50'000'000ULL;
  u64 oracle_collect_limit = 1'000'000ULL;
  u64 oracle_cap = 0;

  // The oracle controls are parsed in main() and passed by cfg.run.extra? No:
  // here we keep them in main(); so this helper just assumes they are already set globally.
  // (We implement the parsing in main and pass as captured values instead.)
  (void)oracle_max_checks;
  (void)oracle_collect_limit;
  (void)oracle_cap;

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

struct VerifyDimRunner {
  const sjs::Config* cfg = nullptr;
  sjs::u64 oracle_max_checks = 0;
  sjs::u64 oracle_collect_limit = 0;
  sjs::u64 oracle_cap = 0;

  template <int Dim>
  int operator()() const {
    std::string err;

    // Load/generate dataset.
    sjs::Dataset<Dim, sjs::Scalar> ds;
    sjs::synthetic::Report gen_report;
    if (!sjs::apps::LoadOrGenerateDataset<Dim>(*cfg, &ds, &gen_report, &err)) {
      SJS_LOG_ERROR("Dataset load/generate failed:", err);
      return 3;
    }

    const sjs::u64 n_r = static_cast<sjs::u64>(ds.R.Size());
    const sjs::u64 n_s = static_cast<sjs::u64>(ds.S.Size());
    const __uint128_t checks128 = static_cast<__uint128_t>(n_r) * static_cast<__uint128_t>(n_s);
    if (checks128 > static_cast<__uint128_t>(oracle_max_checks)) {
      SJS_LOG_ERROR("Oracle would require |R|*|S|=",
                    static_cast<unsigned long long>(checks128),
                    " checks > oracle_max_checks=",
                    static_cast<unsigned long long>(oracle_max_checks),
                    ". Reduce n or raise --oracle_max_checks.");
      return 4;
    }

    SJS_LOG_INFO("Dataset:", ds.name, " Dim=", Dim,
                 " R=", static_cast<unsigned long long>(n_r),
                 " S=", static_cast<unsigned long long>(n_s));

    // Oracle count.
    sjs::join::JoinStats oracle_stats;
    const sjs::u64 oracle_count = sjs::join::CountNaive<Dim, sjs::Scalar>(ds.R, ds.S, &oracle_stats);
    SJS_LOG_INFO("Oracle |J| =", static_cast<unsigned long long>(oracle_count),
                 " candidate_checks=", static_cast<unsigned long long>(oracle_stats.candidate_checks));

    // Oracle universe (optional).
    std::vector<sjs::PairId> universe;
    bool have_full_universe = false;

    if (oracle_count <= oracle_collect_limit) {
      universe = sjs::join::CollectNaivePairs<Dim, sjs::Scalar>(ds.R, ds.S, /*cap=*/0, nullptr);
      have_full_universe = (universe.size() == static_cast<sjs::usize>(oracle_count));
      SJS_LOG_INFO("Collected full universe pairs:", static_cast<unsigned long long>(universe.size()));
    } else if (oracle_cap > 0) {
      universe = sjs::join::CollectNaivePairs<Dim, sjs::Scalar>(ds.R, ds.S, oracle_cap, nullptr);
      have_full_universe = false;
      SJS_LOG_WARN("Join too large to fully collect (|J|=", static_cast<unsigned long long>(oracle_count),
                   "); collected cap=", static_cast<unsigned long long>(universe.size()));
    }

    // Baseline factory.
    std::string factory_err;
    auto baseline = sjs::baselines::CreateBaseline<Dim, sjs::Scalar>(cfg->run.method, cfg->run.variant, &factory_err);
    if (!baseline) {
      SJS_LOG_ERROR("CreateBaseline failed:", factory_err);
      return 5;
    }

    // Run repeats.
    for (sjs::u64 rep = 0; rep < cfg->run.repeats; ++rep) {
      const sjs::u64 seed = cfg->run.seed + 9973ULL * rep;

      sjs::baselines::RunReport out;
      std::string local_err;

      bool ok = false;
      switch (cfg->run.variant) {
        case sjs::Variant::Sampling:
          ok = sjs::baselines::RunSamplingOnce<Dim, sjs::Scalar>(baseline.get(), ds, *cfg, seed, &out, &local_err);
          break;
        case sjs::Variant::EnumSampling:
          ok = sjs::baselines::RunEnumSamplingOnce<Dim, sjs::Scalar>(baseline.get(), ds, *cfg, seed, &out, &local_err);
          break;
        case sjs::Variant::Adaptive:
          ok = sjs::baselines::RunAdaptiveOnce<Dim, sjs::Scalar>(baseline.get(), ds, *cfg, seed, &out, &local_err);
          break;
      }

      std::cout << "---- run rep=" << rep << " seed=" << seed << " dim=" << Dim << " ----\n";
      if (!ok) {
        std::cout << "FAILED: " << local_err << "\n";
        continue;
      }

      const double est = static_cast<double>(out.count.value);
      const double truth = static_cast<double>(oracle_count);
      const double rel = sjs::apps::RelErr(est, truth);

      std::cout << "method=" << sjs::ToString(out.method)
                << " variant=" << sjs::ToString(out.variant)
                << " t=" << out.t << "\n";
      std::cout << "count=" << std::setprecision(17) << est
                << (out.count.exact ? " (exact)" : " (est)")
                << "  oracle=" << truth
                << "  rel_err=" << rel << "\n";
      std::cout << "samples=" << out.samples.Size() << "\n";

      if (have_full_universe && !universe.empty() && !out.samples.pairs.empty()) {
        const auto uni = sjs::sampling::quality::EvaluatePairUniformity(
            sjs::Span<const sjs::PairId>(universe.data(), universe.size()),
            sjs::Span<const sjs::PairId>(out.samples.pairs.data(), out.samples.pairs.size()));

        const double ac1 = sjs::sampling::quality::AutocorrelationHashedPairs(
            sjs::Span<const sjs::PairId>(out.samples.pairs.data(), out.samples.pairs.size()), /*lag=*/1);

        const auto ks = sjs::sampling::quality::KSPairsHashUniform01RankJitter(
            sjs::Span<const sjs::PairId>(universe.data(), universe.size()),
            sjs::Span<const sjs::PairId>(out.samples.pairs.data(), out.samples.pairs.size()));

        std::cout << "quality:\n";
        std::cout << "  missing_in_universe=" << uni.missing_in_universe << "\n";
        std::cout << "  unique_fraction=" << uni.unique_fraction
                  << " duplicate_fraction=" << uni.duplicate_fraction << "\n";
        std::cout << "  l1=" << uni.l1 << " l_inf=" << uni.l_inf
                  << " max_rel_error=" << uni.max_rel_error << "\n";
        std::cout << "  chi2_stat=" << uni.chi2.statistic
                  << " df=" << uni.chi2.df
                  << " p_value=" << uni.chi2.p_value << "\n";
        std::cout << "  autocorr_hash_lag1=" << ac1 << "\n";
        std::cout << "  ks_hash_uniform01 D=" << ks.D << " p=" << ks.p_value << "\n";
      } else {
        std::cout << "quality: skipped (universe not collected)\n";
      }
    }

    return 0;
  }
};

}  // namespace

int main(int argc, char** argv) {
  using sjs::u64;
  using sjs::usize;

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

  // Determine/infer dim & resolve paths.
  if (!sjs::apps::DetermineDatasetDim(&cfg, args, &err)) {
    SJS_LOG_ERROR("Failed to determine dataset dim:", err);
    return 2;
  }
  if (!cfg.Validate(&err)) {
    SJS_LOG_ERROR("Config validation failed after dim inference:", err);
    return 2;
  }

  sjs::Logger::Instance().SetConfig(cfg.logging);

  u64 oracle_max_checks = 50'000'000ULL;
  u64 oracle_collect_limit = 1'000'000ULL;
  u64 oracle_cap = 0;

  sjs::apps::ParseU64(args.Get("oracle_max_checks"), &oracle_max_checks);
  sjs::apps::ParseU64(args.Get("oracle_collect_limit"), &oracle_collect_limit);
  sjs::apps::ParseU64(args.Get("oracle_cap"), &oracle_cap);
  int rc = 0;
  std::string derr;
  if (!sjs::dispatch::DispatchDim<int>(
          cfg.dataset.dim,
          VerifyDimRunner{&cfg, oracle_max_checks, oracle_collect_limit, oracle_cap},
          &rc,
          &derr)) {
    SJS_LOG_ERROR("Unsupported dim=", cfg.dataset.dim, ". ", derr);
    return 2;
  }
  return rc;
}
