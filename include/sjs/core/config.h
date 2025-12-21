#pragma once
// sjs/core/config.h
//
// Experiment configuration (CLI / JSON-friendly).
//
// Design goals:
//  - Keep core free of heavy dependencies (JSON support is optional).
//  - Provide a single Config struct used by all apps (run/sweep/gen/verify).
//  - Be future-proof: supports Dim>2, real datasets, and baseline-specific knobs.
//
// Convention:
//  - CLI uses --key=value or --key value (e.g., --method=ours --variant sampling --t 10000).
//  - Unknown keys are stored into `extra` maps so you can add new knobs without breaking old runs.

#include "sjs/core/types.h"
#include "sjs/core/logging.h"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sjs {

enum class DataSource : u8 {
  Synthetic = 0,
  Binary = 1,
  CSV = 2,
  Unknown = 255,
};

inline constexpr std::string_view ToString(DataSource s) noexcept {
  switch (s) {
    case DataSource::Synthetic: return "synthetic";
    case DataSource::Binary: return "binary";
    case DataSource::CSV: return "csv";
    case DataSource::Unknown: return "unknown";
  }
  return "unknown";
}

inline bool ParseDataSource(std::string_view s, DataSource* out) noexcept {
  if (!out) return false;
  auto eq = detail::EqualsIgnoreCase;
  if (eq(s, "synthetic") || eq(s, "syn")) { *out = DataSource::Synthetic; return true; }
  if (eq(s, "binary") || eq(s, "bin")) { *out = DataSource::Binary; return true; }
  if (eq(s, "csv")) { *out = DataSource::CSV; return true; }
  *out = DataSource::Unknown;
  return false;
}

// --------------------------
// Small key/value argument map
// --------------------------
class ArgMap {
 public:
  ArgMap() = default;

  static ArgMap FromArgv(int argc, char** argv) {
    ArgMap m;
    for (int i = 1; i < argc; ++i) {
      std::string_view token(argv[i]);
      if (token.rfind("--", 0) == 0) {
        token.remove_prefix(2);
        const auto eq_pos = token.find('=');
        if (eq_pos != std::string_view::npos) {
          const std::string key(token.substr(0, eq_pos));
          const std::string val(token.substr(eq_pos + 1));
          m.kv_.emplace(key, val);
        } else {
          const std::string key(token);
          // If next arg exists and isn't another flag, treat as value; else as boolean flag = true.
          if (i + 1 < argc) {
            std::string_view nxt(argv[i + 1]);
            if (!(nxt.rfind("-", 0) == 0)) {
              m.kv_.emplace(key, std::string(nxt));
              ++i;
            } else {
              m.kv_.emplace(key, "true");
            }
          } else {
            m.kv_.emplace(key, "true");
          }
        }
      } else if (token.rfind("-", 0) == 0) {
        // Single-dash options: treat "-k value" as key="k"
        token.remove_prefix(1);
        const std::string key(token);
        if (i + 1 < argc) {
          std::string_view nxt(argv[i + 1]);
          if (!(nxt.rfind("-", 0) == 0)) {
            m.kv_.emplace(key, std::string(nxt));
            ++i;
          } else {
            m.kv_.emplace(key, "true");
          }
        } else {
          m.kv_.emplace(key, "true");
        }
      } else {
        m.positional_.push_back(std::string(token));
      }
    }
    return m;
  }

  bool Has(std::string_view key) const {
    return kv_.find(std::string(key)) != kv_.end();
  }

  std::optional<std::string_view> Get(std::string_view key) const {
    auto it = kv_.find(std::string(key));
    if (it == kv_.end()) return std::nullopt;
    return std::string_view(it->second);
  }

  const std::unordered_map<std::string, std::string>& KV() const { return kv_; }
  const std::vector<std::string>& Positional() const { return positional_; }

 private:
  std::unordered_map<std::string, std::string> kv_;
  std::vector<std::string> positional_;
};

namespace detail {

inline bool ParseBool(std::string_view s, bool* out) noexcept {
  if (!out) return false;
  if (s.empty()) return false;
  if (EqualsIgnoreCase(s, "1") || EqualsIgnoreCase(s, "true") || EqualsIgnoreCase(s, "yes") ||
      EqualsIgnoreCase(s, "y") || EqualsIgnoreCase(s, "on")) {
    *out = true;
    return true;
  }
  if (EqualsIgnoreCase(s, "0") || EqualsIgnoreCase(s, "false") || EqualsIgnoreCase(s, "no") ||
      EqualsIgnoreCase(s, "n") || EqualsIgnoreCase(s, "off")) {
    *out = false;
    return true;
  }
  return false;
}

inline bool ParseU64(std::string_view s, u64* out) {
  if (!out) return false;
  if (s.empty()) return false;
  // stod/stoull are fine for config; speed isn't critical here.
  try {
    std::size_t idx = 0;
    const unsigned long long v = std::stoull(std::string(s), &idx, 10);
    if (idx != s.size()) return false;
    *out = static_cast<u64>(v);
    return true;
  } catch (...) {
    return false;
  }
}

inline bool ParseI32(std::string_view s, i32* out) {
  if (!out) return false;
  if (s.empty()) return false;
  try {
    std::size_t idx = 0;
    const long v = std::stol(std::string(s), &idx, 10);
    if (idx != s.size()) return false;
    *out = static_cast<i32>(v);
    return true;
  } catch (...) {
    return false;
  }
}

inline bool ParseDouble(std::string_view s, double* out) {
  if (!out) return false;
  if (s.empty()) return false;
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

// Insert unknown keys into a map<string,string>.
inline void StoreExtras(const std::unordered_map<std::string, std::string>& all_kv,
                        const std::vector<std::string>& known_keys,
                        std::unordered_map<std::string, std::string>* out_extra) {
  if (!out_extra) return;
  out_extra->clear();
  out_extra->reserve(all_kv.size());

  auto is_known = [&](const std::string& k) -> bool {
    for (const auto& kk : known_keys) {
      if (k == kk) return true;
    }
    return false;
  };

  for (const auto& kv : all_kv) {
    if (!is_known(kv.first)) {
      (*out_extra)[kv.first] = kv.second;
    }
  }
}

}  // namespace detail

// --------------------------
// Dataset config
// --------------------------
struct SyntheticConfig {
  // Generator name (e.g., "stripe_ctrl_alpha", "uniform", "clustered"...)
  std::string generator = "stripe_ctrl_alpha";

  // Sizes for R and S
  u64 n_r = 100000;
  u64 n_s = 100000;

  // A high-level knob for density; exact semantics depend on generator.
  // For your stripe generator this maps naturally to alpha_out.
  double alpha = 1e-6;

  // Generator seed (separate from sampling seed)
  u64 seed = 1;

  // Extra generator-specific knobs, passed through.
  std::unordered_map<std::string, std::string> extra;
};

struct DatasetConfig {
  DataSource source = DataSource::Synthetic;

  // Semantic name for logs/results.
  std::string name = "synthetic";

  // Dimensionality (keep as runtime value; most algorithms are templated by Dim but
  // experiments choose Dim at runtime).
  i32 dim = kDefaultDim;

  // For binary/csv datasets
  std::string path_r;
  std::string path_s;

  // For synthetic datasets
  SyntheticConfig synthetic;
};

// --------------------------
// Run config
// --------------------------
struct RunConfig {
  Method method = Method::Ours;
  Variant variant = Variant::Sampling;

  // Sample size t
  u64 t = 10000;

  // Sampling seed
  u64 seed = 1;

  // Number of repeats (for median/CI)
  u64 repeats = 5;

  // Adaptive threshold (implementation-defined; e.g., J* or time budget).
  // Keep as generic knobs; the adaptive runner can interpret them.
  u64 j_star = 1000000;

  // Stop enumerate early after this many discovered pairs (0 = disabled).
  u64 enum_cap = 0;

  // Emit sampled pairs to disk (can be large!)
  bool write_samples = false;

  // Verify correctness on small datasets (oracle compare, chi2, etc.)
  bool verify = false;

  // Extra baseline-specific knobs (e.g., KD leaf size, PBSM grid, etc.)
  std::unordered_map<std::string, std::string> extra;
};

struct OutputConfig {
  std::string out_dir = "results/raw";
  std::string run_tag;  // optional, helps grouping sweeps
};

struct SystemConfig {
  i32 threads = 1;
};

struct Config {
  DatasetConfig dataset;
  RunConfig run;
  OutputConfig output;
  SystemConfig sys;
  LoggingConfig logging;

  // Validate basic constraints; returns false and sets err on failure.
  bool Validate(std::string* err = nullptr) const {
    auto fail = [&](std::string_view msg) {
      if (err) *err = std::string(msg);
      return false;
    };

    if (dataset.dim <= 0) return fail("dataset.dim must be > 0");
    if (dataset.dim > kMaxSupportedDim) return fail("dataset.dim too large");
    if (run.t == 0) return fail("run.t must be > 0");
    if (run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      return fail("run.t too large (must fit in u32); reduce t or update baselines to use u64");
    }
    if (run.repeats == 0) return fail("run.repeats must be > 0");
    if (sys.threads <= 0) return fail("sys.threads must be > 0");

    if (dataset.source == DataSource::Binary || dataset.source == DataSource::CSV) {
      if (dataset.path_r.empty() || dataset.path_s.empty()) {
        return fail("dataset.path_r and dataset.path_s must be set for non-synthetic datasets");
      }
    }
    if (dataset.source == DataSource::Synthetic) {
      if (dataset.synthetic.n_r == 0 || dataset.synthetic.n_s == 0) {
        return fail("synthetic.n_r and synthetic.n_s must be > 0");
      }
      if (!(dataset.synthetic.alpha >= 0.0)) {
        return fail("synthetic.alpha must be >= 0");
      }
    }
    return true;
  }

  std::string ToJsonLite() const {
    std::ostringstream oss;
    oss << "{"
        << "\"dataset\":{\"source\":\"" << ToString(dataset.source) << "\","
        << "\"name\":\"" << dataset.name << "\","
        << "\"dim\":" << dataset.dim << ","
        << "\"path_r\":\"" << dataset.path_r << "\","
        << "\"path_s\":\"" << dataset.path_s << "\","
        << "\"synthetic\":{\"generator\":\"" << dataset.synthetic.generator << "\","
        << "\"n_r\":" << dataset.synthetic.n_r << ","
        << "\"n_s\":" << dataset.synthetic.n_s << ","
        << "\"alpha\":" << dataset.synthetic.alpha << ","
        << "\"seed\":" << dataset.synthetic.seed << "}"
        << "},"
        << "\"run\":{\"method\":\"" << ToString(run.method) << "\","
        << "\"variant\":\"" << ToString(run.variant) << "\","
        << "\"t\":" << run.t << ","
        << "\"seed\":" << run.seed << ","
        << "\"repeats\":" << run.repeats << ","
        << "\"j_star\":" << run.j_star << ","
        << "\"enum_cap\":" << run.enum_cap << ","
        << "\"write_samples\":" << (run.write_samples ? "true" : "false") << ","
        << "\"verify\":" << (run.verify ? "true" : "false")
        << "},"
        << "\"output\":{\"out_dir\":\"" << output.out_dir << "\","
        << "\"run_tag\":\"" << output.run_tag << "\"},"
        << "\"sys\":{\"threads\":" << sys.threads << "},"
        << "\"logging\":{\"level\":\"" << ToString(logging.level) << "\","
        << "\"with_timestamp\":" << (logging.with_timestamp ? "true" : "false") << ","
        << "\"with_thread_id\":" << (logging.with_thread_id ? "true" : "false") << "}"
        << "}";
    return oss.str();
  }

  // Parse from CLI arguments. Unknown args go into `run.extra` / `synthetic.extra`.
  static Config FromArgs(int argc, char** argv) {
    const ArgMap args = ArgMap::FromArgv(argc, argv);
    Config cfg;

    // Known keys list (anything else will be stored into extra).
    std::vector<std::string> known = {
        "method","variant","t","seed","repeats","j_star","enum_cap","write_samples","verify",
        "dataset_source","dataset","dim","path_r","path_s",
        "gen","n_r","n_s","alpha","gen_seed",
        "out_dir","run_tag",
        "threads",
        "log_level","log_timestamp","log_thread",
    };

    // ------------- run -------------
    if (auto v = args.Get("method")) {
      Method m;
      if (ParseMethod(*v, &m) && m != Method::Unknown) cfg.run.method = m;
    }
    if (auto v = args.Get("variant")) {
      Variant var;
      if (ParseVariant(*v, &var)) cfg.run.variant = var;
    }
    if (auto v = args.Get("t")) { detail::ParseU64(*v, &cfg.run.t); }
    if (auto v = args.Get("seed")) { detail::ParseU64(*v, &cfg.run.seed); }
    if (auto v = args.Get("repeats")) { detail::ParseU64(*v, &cfg.run.repeats); }
    if (auto v = args.Get("j_star")) { detail::ParseU64(*v, &cfg.run.j_star); }
    if (auto v = args.Get("enum_cap")) { detail::ParseU64(*v, &cfg.run.enum_cap); }
    if (auto v = args.Get("write_samples")) { detail::ParseBool(*v, &cfg.run.write_samples); }
    if (auto v = args.Get("verify")) { detail::ParseBool(*v, &cfg.run.verify); }

    // ------------- dataset -------------
    if (auto v = args.Get("dataset_source")) {
      DataSource ds;
      if (ParseDataSource(*v, &ds) && ds != DataSource::Unknown) cfg.dataset.source = ds;
    }
    if (auto v = args.Get("dataset")) cfg.dataset.name = std::string(*v);
    if (auto v = args.Get("dim")) { detail::ParseI32(*v, &cfg.dataset.dim); }
    if (auto v = args.Get("path_r")) cfg.dataset.path_r = std::string(*v);
    if (auto v = args.Get("path_s")) cfg.dataset.path_s = std::string(*v);

    // synthetic knobs
    if (auto v = args.Get("gen")) cfg.dataset.synthetic.generator = std::string(*v);
    if (auto v = args.Get("n_r")) detail::ParseU64(*v, &cfg.dataset.synthetic.n_r);
    if (auto v = args.Get("n_s")) detail::ParseU64(*v, &cfg.dataset.synthetic.n_s);
    if (auto v = args.Get("alpha")) detail::ParseDouble(*v, &cfg.dataset.synthetic.alpha);
    if (auto v = args.Get("gen_seed")) detail::ParseU64(*v, &cfg.dataset.synthetic.seed);

    // ------------- output -------------
    if (auto v = args.Get("out_dir")) cfg.output.out_dir = std::string(*v);
    if (auto v = args.Get("run_tag")) cfg.output.run_tag = std::string(*v);

    // ------------- system -------------
    if (auto v = args.Get("threads")) detail::ParseI32(*v, &cfg.sys.threads);

    // ------------- logging -------------
    if (auto v = args.Get("log_level")) {
      LogLevel lvl;
      if (detail::ParseLogLevel(*v, &lvl)) cfg.logging.level = lvl;
    }
    if (auto v = args.Get("log_timestamp")) detail::ParseBool(*v, &cfg.logging.with_timestamp);
    if (auto v = args.Get("log_thread")) detail::ParseBool(*v, &cfg.logging.with_thread_id);

    // Store unknown keys into extras for forward compatibility.
    detail::StoreExtras(args.KV(), known, &cfg.run.extra);
    // Generator extras can share the same map; you can split by prefix later if you want.
    cfg.dataset.synthetic.extra = cfg.run.extra;

    return cfg;
  }
};

}  // namespace sjs
