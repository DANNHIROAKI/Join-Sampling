// apps/sjs_gen_dataset.cpp
//
// HighDims synthetic dataset generator + exporter.
//
// Default mode (C++ generators):
//   Generates a Dataset<Dim> (R/S) using the configured synthetic generator, then
//   writes it out in:
//     - SJS binary format (fast, robust), and optionally
//     - CSV (for debugging/visualization)
//
// Special mode (Python generator: alacarte-rectgen):
//   If --gen=alacarte_rectgen, we delegate generation to:
//     tools/alacarte_rectgen_generate.py
//   and let Python write the binary files directly. This avoids loading huge
//   datasets into C++ memory and supports *any* dimension d.
//
// Dim handling:
//   - requires --dataset_source=synthetic
//   - uses --dim=<d>
//   - for C++ generators, Dim must be in SJS_DIMS (see include/sjs/dispatch/dim_dispatch.h)
//   - for alacarte_rectgen, any d>0 is accepted (Python generates and writes files)

#include "sjs/baselines/baseline_api.h"
#include "sjs/core/logging.h"
#include "sjs/core/types.h"

#include "sjs/dispatch/dim_dispatch.h"

#include "sjs/data/synthetic/generator_factory.h"
#include "sjs/io/binary_io.h"
#include "sjs/io/csv_io.h"
#include "sjs/io/dataset.h"

#include <cctype>
#include <cstdlib>
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
    const double v = std::stod(std::string(s), &idx);
    if (idx != s.size()) return false;
    *out = v;
    return true;
  } catch (...) {
    return false;
  }
}

inline bool IsReservedKey(std::string_view key) {
  const std::string k(key);
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
      return fail("Invalid --dataset_source (expected synthetic|binary|csv). ");
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

  // run (unused here but we keep parser consistent across apps)
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

inline bool IsHelpRequested(const sjs::ArgMap& args) {
  return args.Has("help") || args.Has("h") || args.Has("--help") || args.Has("-h");
}

inline std::string SanitizeFilename(std::string_view s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.') out.push_back(c);
    else out.push_back('_');
  }
  if (out.empty()) out = "dataset";
  return out;
}

inline bool EnsureDir(const fs::path& p, std::string* err) {
  std::error_code ec;
  if (fs::exists(p, ec)) return true;
  if (fs::create_directories(p, ec)) return true;
  if (err) *err = ec.message();
  return false;
}

inline bool ParseBool(const std::optional<std::string_view>& v, bool def) {
  if (!v) return def;
  if (*v == "1" || *v == "true" || *v == "yes") return true;
  if (*v == "0" || *v == "false" || *v == "no") return false;
  return def;
}

inline std::string ShellQuote(std::string_view s) {
#if defined(_WIN32)
  // Minimal Windows quoting.
  std::string out;
  out.reserve(s.size() + 2);
  out.push_back('"');
  for (char c : s) {
    if (c == '"') out += "\\\"";
    else out.push_back(c);
  }
  out.push_back('"');
  return out;
#else
  // POSIX sh-safe single-quote escaping: ' -> '\''
  std::string out;
  out.reserve(s.size() + 2);
  out.push_back('\'');
  for (char c : s) {
    if (c == '\'') out += "'\\''";
    else out.push_back(c);
  }
  out.push_back('\'');
  return out;
#endif
}

inline std::string GetExtraOrEnv(const std::unordered_map<std::string, std::string>& extra,
                                 std::string_view key,
                                 const char* env_key,
                                 std::string def) {
  if (auto it = extra.find(std::string(key)); it != extra.end()) return it->second;
  if (env_key) {
    if (const char* v = std::getenv(env_key); v && *v) return std::string(v);
  }
  return def;
}

inline u64 GetExtraU64(const std::unordered_map<std::string, std::string>& extra,
                       std::string_view key,
                       u64 def) {
  auto it = extra.find(std::string(key));
  if (it == extra.end()) return def;
  try {
    std::size_t idx = 0;
    const unsigned long long v = std::stoull(it->second, &idx, 10);
    if (idx != it->second.size()) return def;
    return static_cast<u64>(v);
  } catch (...) {
    return def;
  }
}

inline bool IsAlacarteRectGenName(std::string_view s) {
  using sjs::detail::EqualsIgnoreCase;
  return EqualsIgnoreCase(s, "alacarte_rectgen") ||
         EqualsIgnoreCase(s, "alacarte-rectgen") ||
         EqualsIgnoreCase(s, "rectgen") ||
         EqualsIgnoreCase(s, "alacarte");
}

inline void PrintUsage() {
  std::cerr
      << "sjs_gen_dataset: synthetic dataset generator (HighDims)\n\n"
      << "Required:\n"
      << "  --dataset_source=synthetic\n"
      << "  --dim=<d>\n"
      << "  --gen=<stripe|uniform|clustered|hetero_sizes|alacarte_rectgen>\n"
      << "  --dataset=<name>\n"
      << "  --n_r=<N> --n_s=<N>\n"
      << "  --alpha=<float>\n"
      << "  --gen_seed=<seed>\n"
      << "\nOutput:\n"
      << "  --out_dir=<dir>\n"
      << "  --out_r=<path> --out_s=<path>            (override binary outputs)\n"
      << "  --write_csv=0|1                          (default 0; debug)\n"
      << "  --csv_r=<path> --csv_s=<path>            (override CSV outputs)\n"
      << "  --csv_sep=,|tab|\\t                       (default ',')\n"
      << "\nAlacarte RectGen (Python) extra knobs:\n"
      << "  --python=<python_exe>                    (default: env SJS_PYTHON or python3)\n"
      << "  --rectgen_script=<path>                  (default: tools/alacarte_rectgen_generate.py)\n"
      << "  --audit_pairs=<u64>                      (default: 2000000; 0 disables audit)\n"
      << "  --audit_seed=<u64>                       (default: 0)\n"
      << "\n";
}

// ----------------------------
// Special: Python generator (alacarte-rectgen)
// ----------------------------

int RunAlacarteRectGenPython(const sjs::Config& cfg, const sjs::ArgMap& args) {
  std::string err;

  if (cfg.dataset.source != sjs::DataSource::Synthetic) {
    SJS_LOG_ERROR("sjs_gen_dataset requires --dataset_source=synthetic");
    return 2;
  }
  if (cfg.dataset.dim <= 0) {
    SJS_LOG_ERROR("Invalid --dim. Must be > 0.");
    return 2;
  }

  // Output paths.
  fs::path out_dir(cfg.output.out_dir);
  if (!args.Has("out_dir")) {
    out_dir = fs::path("data/synthetic");
  }
  if (!EnsureDir(out_dir, &err)) {
    SJS_LOG_ERROR("Cannot create out_dir:", err);
    return 4;
  }

  const int dim = cfg.dataset.dim;
  const std::string base = SanitizeFilename(cfg.dataset.name) + "__d" + std::to_string(dim);

  std::string out_r = (out_dir / (base + "_R.bin")).string();
  std::string out_s = (out_dir / (base + "_S.bin")).string();
  if (auto v = args.Get("out_r")) out_r = std::string(*v);
  if (auto v = args.Get("out_s")) out_s = std::string(*v);

  const std::string rep_path = (out_dir / (base + "_gen_report.json")).string();

  // Optional CSV.
  const bool write_csv = ParseBool(args.Get("write_csv"), /*def=*/false);
  char sep = ',';
  if (auto v = args.Get("csv_sep")) {
    if (*v == "tab" || *v == "\\t") sep = '\t';
    else if (!v->empty()) sep = (*v)[0];
  }
  std::string csv_r = (out_dir / (base + "_R.csv")).string();
  std::string csv_s = (out_dir / (base + "_S.csv")).string();
  if (auto v = args.Get("csv_r")) csv_r = std::string(*v);
  if (auto v = args.Get("csv_s")) csv_s = std::string(*v);

  // Resolve Python + script.
  const auto& extra = cfg.dataset.synthetic.extra;
  const std::string python = GetExtraOrEnv(extra, "python", "SJS_PYTHON", "python3");
  std::string script = GetExtraOrEnv(extra, "rectgen_script", "SJS_ALACARTE_RECTGEN_SCRIPT",
                                     "tools/alacarte_rectgen_generate.py");
  if (!fs::exists(script)) {
    // Try a couple of common fallbacks.
    const fs::path p0 = fs::path("tools") / "alacarte_rectgen_generate.py";
    const fs::path p1 = fs::path("../tools") / "alacarte_rectgen_generate.py";
    const fs::path p2 = fs::path("../../tools") / "alacarte_rectgen_generate.py";
    const fs::path p3 = fs::path("../../../tools") / "alacarte_rectgen_generate.py";
    if (fs::exists(p0)) script = p0.string();
    else if (fs::exists(p1)) script = p1.string();
    else if (fs::exists(p2)) script = p2.string();
    else if (fs::exists(p3)) script = p3.string();
  }
  if (!fs::exists(script)) {
    SJS_LOG_ERROR("alacarte_rectgen script not found:", script,
                  "\n  Provide --rectgen_script=... or set env SJS_ALACARTE_RECTGEN_SCRIPT.");
    return 3;
  }

  const u64 audit_pairs = GetExtraU64(extra, "audit_pairs", 2'000'000ULL);
  const u64 audit_seed = GetExtraU64(extra, "audit_seed", 0ULL);

  // Build command.
  std::ostringstream cmd;
  cmd << ShellQuote(python)
      << " " << ShellQuote(script)
      << " --nR=" << cfg.dataset.synthetic.n_r
      << " --nS=" << cfg.dataset.synthetic.n_s
      << " --d=" << dim
      << " --alpha_out=" << cfg.dataset.synthetic.alpha
      << " --seed=" << cfg.dataset.synthetic.seed
      << " --out_r=" << ShellQuote(out_r)
      << " --out_s=" << ShellQuote(out_s)
      << " --dataset_name=" << ShellQuote(cfg.dataset.name)
      << " --report_path=" << ShellQuote(rep_path)
      << " --audit_pairs=" << audit_pairs
      << " --audit_seed=" << audit_seed;

  if (write_csv) {
    cmd << " --write_csv"
        << " --csv_r=" << ShellQuote(csv_r)
        << " --csv_s=" << ShellQuote(csv_s)
        << " --csv_sep=" << ShellQuote(std::string(1, sep));
  }

  SJS_LOG_INFO("Running alacarte_rectgen Python generator...\nCmd: ", cmd.str());
  const int rc = std::system(cmd.str().c_str());
  if (rc != 0) {
    SJS_LOG_ERROR("alacarte_rectgen generation failed (rc=", rc,
                  "). Make sure you installed the dependency: pip install alacarte-rectgen");
    return 3;
  }

  if (!fs::exists(out_r) || !fs::exists(out_s)) {
    SJS_LOG_ERROR("Python generator returned success, but output files are missing.");
    SJS_LOG_ERROR("Expected:", out_r, "and", out_s);
    return 3;
  }

  SJS_LOG_INFO("Wrote binary:", out_r, "and", out_s);
  if (write_csv) {
    SJS_LOG_INFO("Wrote CSV:", csv_r, "and", csv_s);
  }
  if (fs::exists(rep_path)) {
    SJS_LOG_INFO("Wrote generator report:", rep_path);
  } else {
    SJS_LOG_WARN("Generator report missing:", rep_path);
  }
  return 0;
}

// ----------------------------
// Default: C++ synthetic generators
// ----------------------------

template <int Dim>
int RunForDim(const sjs::Config& cfg, const sjs::ArgMap& args) {
  std::string err;

  if (cfg.dataset.source != sjs::DataSource::Synthetic) {
    SJS_LOG_ERROR("sjs_gen_dataset requires --dataset_source=synthetic");
    return 2;
  }
  if (cfg.dataset.dim != Dim) {
    SJS_LOG_ERROR("Dim dispatch mismatch (cfg.dataset.dim != Dim).");
    return 2;
  }

  // Generate dataset.
  sjs::Dataset<Dim, sjs::Scalar> ds;
  sjs::synthetic::DatasetSpec spec;
  spec.name = cfg.dataset.name;
  spec.n_r = cfg.dataset.synthetic.n_r;
  spec.n_s = cfg.dataset.synthetic.n_s;
  spec.alpha = cfg.dataset.synthetic.alpha;
  spec.seed = cfg.dataset.synthetic.seed;
  spec.params = cfg.dataset.synthetic.extra;

  sjs::synthetic::Report rep;
  if (!sjs::synthetic::GenerateDataset<Dim, sjs::Scalar>(
          cfg.dataset.synthetic.generator, spec, &ds, &rep, &err)) {
    SJS_LOG_ERROR("Generation failed:", err);
    return 3;
  }

  SJS_LOG_INFO("Generated dataset:", ds.name,
               "Dim=", Dim,
               "R=", static_cast<unsigned long long>(ds.R.Size()),
               "S=", static_cast<unsigned long long>(ds.S.Size()),
               "report=", rep.ToJsonLite());

  // Output paths.
  fs::path out_dir(cfg.output.out_dir);
  // Default output directory for this generator app is data/synthetic unless overridden.
  if (!args.Has("out_dir")) {
    out_dir = fs::path("data/synthetic");
  }

  if (!EnsureDir(out_dir, &err)) {
    SJS_LOG_ERROR("Cannot create out_dir:", err);
    return 4;
  }

  const std::string base = SanitizeFilename(ds.name) + "__d" + std::to_string(Dim);

  std::string out_r = (out_dir / (base + "_R.bin")).string();
  std::string out_s = (out_dir / (base + "_S.bin")).string();
  if (auto v = args.Get("out_r")) out_r = std::string(*v);
  if (auto v = args.Get("out_s")) out_s = std::string(*v);

  // Write binary.
  {
    sjs::binary::BinaryWriteOptions opt;
    opt.half_open = true;
    opt.write_ids = true;
    opt.scalar = sjs::binary::ScalarEncoding::Float64;
    opt.write_name = true;

    if (!sjs::binary::WriteDatasetBinaryPair<Dim, sjs::Scalar>(out_r, out_s, ds, opt, &err)) {
      SJS_LOG_ERROR("Binary write failed:", err);
      return 5;
    }
    SJS_LOG_INFO("Wrote binary:", out_r, "and", out_s);
  }

  // Optional CSV.
  const bool write_csv = ParseBool(args.Get("write_csv"), /*def=*/false);
  if (write_csv) {
    char sep = ',';
    if (auto v = args.Get("csv_sep")) {
      if (*v == "tab" || *v == "\\t") sep = '\t';
      else if (!v->empty()) sep = (*v)[0];
    }

    std::string csv_r = (out_dir / (base + "_R.csv")).string();
    std::string csv_s = (out_dir / (base + "_S.csv")).string();
    if (auto v = args.Get("csv_r")) csv_r = std::string(*v);
    if (auto v = args.Get("csv_s")) csv_s = std::string(*v);

    sjs::csv::Dialect d;
    d.sep = sep;
    d.write_header = true;

    if (!sjs::csv::WriteBoxes<Dim, sjs::Scalar>(csv_r, ds.R, d, &err)) {
      SJS_LOG_ERROR("CSV write R failed:", err);
      return 6;
    }
    if (!sjs::csv::WriteBoxes<Dim, sjs::Scalar>(csv_s, ds.S, d, &err)) {
      SJS_LOG_ERROR("CSV write S failed:", err);
      return 6;
    }
    SJS_LOG_INFO("Wrote CSV:", csv_r, "and", csv_s);
  }

  // Write generator report JSON-lite.
  {
    const std::string rep_path = (out_dir / (base + "_gen_report.json")).string();
    std::ofstream out(rep_path);
    if (out) {
      out << rep.ToJsonLite() << "\n";
      SJS_LOG_INFO("Wrote generator report:", rep_path);
    }
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

struct GenDatasetDimRunner {
  const sjs::Config* cfg = nullptr;
  const sjs::ArgMap* args = nullptr;

  template <int Dim>
  int operator()() const {
    return sjs::apps::RunForDim<Dim>(*cfg, *args);
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
  sjs::Logger::Instance().SetConfig(cfg.logging);

  if (cfg.dataset.source != sjs::DataSource::Synthetic) {
    SJS_LOG_ERROR("sjs_gen_dataset requires --dataset_source=synthetic");
    return 2;
  }

  // Special-case: alacarte-rectgen (Python) supports ANY dim.
  if (sjs::apps::IsAlacarteRectGenName(cfg.dataset.synthetic.generator)) {
    return sjs::apps::RunAlacarteRectGenPython(cfg, args);
  }

  // Default: C++ synthetic generators use compile-time dim dispatch.
  if (cfg.dataset.dim <= 0) {
    SJS_LOG_ERROR("Invalid --dim. Must be > 0.");
    return 2;
  }
  int rc = 0;
  std::string derr;
  if (!sjs::dispatch::DispatchDim<int>(cfg.dataset.dim, GenDatasetDimRunner{&cfg, &args}, &rc, &derr)) {
    SJS_LOG_ERROR("Unsupported dim=", cfg.dataset.dim, ". ", derr);
    return 2;
  }
  return rc;
}
