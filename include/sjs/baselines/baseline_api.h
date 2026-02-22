#pragma once
// sjs/baselines/baseline_api.h
//
// Unified baseline interface for Join-Sampling experiments (HighDims).
//
// This header provides:
//   - A lightweight experiment Config (since HighDims repo removes core/config.h)
//   - CountResult, SampleSet
//   - IJoinEnumerator (deterministic join stream + stats)
//   - IBaseline<Dim,T> interface used by runners
//   - RunReport container for CSV/JSON writing
//
// Design goals:
//   - Keep dependencies light (no JSON library here).
//   - Keep interface stable across Dim >= 2.
//   - Baselines may interpret cfg.run.extra in their own way.

#include "sjs/core/types.h"
#include "sjs/core/rng.h"
#include "sjs/core/timer.h"
#include "sjs/core/logging.h"

#include "sjs/io/dataset.h"
#include "sjs/join/join_types.h"

#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sjs {

// --------------------------
// Data / run configuration
// --------------------------

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

struct SyntheticConfig {
  // Generator name (see sjs/data/synthetic/*)
  std::string generator = "stripe_ctrl_alpha";

  // Sizes
  u64 n_r = 100000;
  u64 n_s = 100000;

  // A generic density knob (generator-specific semantics).
  double alpha = 1e-6;

  // Generator seed (separate from sampling seed).
  u64 seed = 1;

  // Generator-specific passthrough knobs.
  std::unordered_map<std::string, std::string> extra;
};

struct DatasetConfig {
  DataSource source = DataSource::Synthetic;

  // Semantic dataset name for logs/results.
  std::string name = "synthetic";

  // Runtime Dim (apps dispatch to template Dim).
  i32 dim = kDefaultDim;

  // For binary/csv datasets:
  std::string path_r;
  std::string path_s;

  // For synthetic datasets:
  SyntheticConfig synthetic;
};

struct RunConfig {
  // Algorithm family and variant selector (see core/types.h).
  Method method = Method::Ours;
  Variant variant = Variant::Sampling;

  // Sample size t (number of i.i.d. outputs).
  // NOTE: t==0 is allowed (count-only run; sampling should return an empty SampleSet).
  u64 t = 10000;

  // Sampling seed.
  u64 seed = 1;

  // Number of repeats (for median/CI sweeps).
  u64 repeats = 5;

  // Framework III (Variant::Adaptive) knobs.
  //
  // These knobs affect performance only (caching/prefetch), not the target
  // i.i.d. uniform distribution over join pairs.
  //
  // - budget: total cache/prefetch budget B (in #partner ids). 0 disables caching/prefetch.
  // - w_small: full-cache threshold; if w_i <= w_small and budget allows, cache full partner list.
  // - prefetch: enable sample prefetch cache (may reduce/eliminate the second pass).
  u64 budget = 1000000;
  u64 w_small = 0;
  bool prefetch = true;

  // DEPRECATED legacy alias for budget (kept for compatibility with older configs).
  u64 j_star = 1000000;

  // Enumeration cap for EnumSampling/Adaptive (0 = no cap).
  u64 enum_cap = 0;

  // Emit sampled pairs to disk (can be large).
  bool write_samples = false;

  // Enable correctness verification (small datasets).
  bool verify = false;

  // Extra baseline-specific knobs (stringly typed).
  std::unordered_map<std::string, std::string> extra;
};

struct OutputConfig {
  std::string out_dir = "results/raw";
  std::string run_tag;  // optional
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

  bool Validate(std::string* err = nullptr) const {
    auto fail = [&](std::string_view msg) {
      if (err) *err = std::string(msg);
      return false;
    };

    if (dataset.dim <= 0) return fail("dataset.dim must be > 0");
    if (dataset.dim > kMaxSupportedDim) return fail("dataset.dim too large");
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
        << "\"budget\":" << run.budget << ","
        << "\"w_small\":" << run.w_small << ","
        << "\"prefetch\":" << (run.prefetch ? "true" : "false") << ","
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
};

}  // namespace sjs

namespace sjs {
namespace baselines {

// --------------------------
// CountResult
// --------------------------
struct CountResult {
  long double value = 0.0L;
  bool exact = false;

  long double stderr = std::numeric_limits<long double>::quiet_NaN();
  long double ci_low = std::numeric_limits<long double>::quiet_NaN();
  long double ci_high = std::numeric_limits<long double>::quiet_NaN();

  u64 aux_draws = 0;

  bool HasStdErr() const noexcept { return std::isfinite(static_cast<double>(stderr)); }
  bool HasCI() const noexcept {
    return std::isfinite(static_cast<double>(ci_low)) && std::isfinite(static_cast<double>(ci_high));
  }

  u64 RoundedU64() const noexcept {
    if (!(value > 0.0L)) return 0ULL;
    long double x = value;
    if (x > static_cast<long double>(std::numeric_limits<u64>::max())) {
      return std::numeric_limits<u64>::max();
    }
    return static_cast<u64>(x + 0.5L);
  }

  std::string ToJsonLite() const {
    std::ostringstream oss;
    oss << "{";
    oss << "\"value\":" << static_cast<double>(value) << ",";
    oss << "\"exact\":" << (exact ? "true" : "false") << ",";
    if (HasStdErr()) {
      oss << "\"stderr\":" << static_cast<double>(stderr) << ",";
    }
    if (HasCI()) {
      oss << "\"ci_low\":" << static_cast<double>(ci_low) << ",";
      oss << "\"ci_high\":" << static_cast<double>(ci_high) << ",";
    }
    oss << "\"aux_draws\":" << aux_draws;
    oss << "}";
    return oss.str();
  }
};

inline CountResult MakeExactCount(u64 v) {
  CountResult r;
  r.value = static_cast<long double>(v);
  r.exact = true;
  r.stderr = std::numeric_limits<long double>::quiet_NaN();
  r.ci_low = std::numeric_limits<long double>::quiet_NaN();
  r.ci_high = std::numeric_limits<long double>::quiet_NaN();
  r.aux_draws = 0;
  return r;
}

inline CountResult MakeEstimateCount(long double v,
                                    long double stderr = std::numeric_limits<long double>::quiet_NaN(),
                                    long double ci_low = std::numeric_limits<long double>::quiet_NaN(),
                                    long double ci_high = std::numeric_limits<long double>::quiet_NaN(),
                                    u64 aux_draws = 0) {
  CountResult r;
  r.value = v;
  r.exact = false;
  r.stderr = stderr;
  r.ci_low = ci_low;
  r.ci_high = ci_high;
  r.aux_draws = aux_draws;
  return r;
}

// --------------------------
// SampleSet
// --------------------------
struct SampleSet {
  std::vector<PairId> pairs;
  std::vector<double> weights;

  bool weighted = false;
  bool with_replacement = true;

  void Clear() {
    pairs.clear();
    weights.clear();
    weighted = false;
    with_replacement = true;
  }

  usize Size() const noexcept { return pairs.size(); }
  bool Empty() const noexcept { return pairs.empty(); }

  bool Validate(std::string* err = nullptr) const {
    if (!weighted) {
      if (!weights.empty()) {
        if (err) *err = "SampleSet::Validate: weights must be empty when weighted=false";
        return false;
      }
      return true;
    }
    if (weights.size() != pairs.size()) {
      if (err) *err = "SampleSet::Validate: weights.size() != pairs.size()";
      return false;
    }
    for (double w : weights) {
      if (!(w >= 0.0) || !std::isfinite(w)) {
        if (err) *err = "SampleSet::Validate: weight must be finite and >= 0";
        return false;
      }
    }
    return true;
  }
};

// --------------------------
// IJoinEnumerator
// --------------------------
class IJoinEnumerator : public join::IJoinStream {
 public:
  ~IJoinEnumerator() override = default;
  virtual const join::JoinStats& Stats() const noexcept = 0;
};

// --------------------------
// IBaseline
// --------------------------
template <int Dim, class T = Scalar>
class IBaseline {
 public:
  using BoxT = Box<Dim, T>;
  using DatasetT = Dataset<Dim, T>;

  virtual ~IBaseline() = default;

  virtual Method method() const noexcept = 0;
  virtual Variant variant() const noexcept = 0;
  virtual std::string_view Name() const noexcept = 0;

  virtual void Reset() = 0;

  virtual bool Build(const DatasetT& ds,
                     const sjs::Config& cfg,
                     PhaseRecorder* phases,
                     std::string* err) = 0;

  virtual bool Count(const sjs::Config& cfg,
                     Rng* rng,
                     CountResult* out,
                     PhaseRecorder* phases,
                     std::string* err) = 0;

  virtual bool Sample(const sjs::Config& cfg,
                      Rng* rng,
                      SampleSet* out,
                      PhaseRecorder* phases,
                      std::string* err) = 0;

  virtual std::unique_ptr<IJoinEnumerator> Enumerate(const sjs::Config& cfg,
                                                     PhaseRecorder* phases,
                                                     std::string* err) = 0;
};

// --------------------------
// Runner-facing report
// --------------------------
struct RunReport {
  bool ok = false;
  std::string error;

  Method method = Method::Unknown;
  Variant variant = Variant::Sampling;
  std::string baseline_name;
  std::string dataset_name;

  // Compile-time dimension (matches the templated Dim used by runners).
  // This is stored in the report so that writers (CSV/JSON/filenames) can
  // always include `dim` without relying on fragile `if constexpr` tricks.
  i32 dim = 0;

  u64 seed = 0;
  u64 t = 0;

  CountResult count;
  SampleSet samples;

  bool used_enumeration = false;
  bool enumeration_truncated = false;
  u64 enumeration_cap = 0;
  u64 enumeration_pairs_pass1 = 0;
  u64 enumeration_pairs_pass2 = 0;
  join::JoinStats enum_stats_pass1{};
  join::JoinStats enum_stats_pass2{};

  std::string adaptive_branch;
  u64 adaptive_pilot_pairs = 0;

  PhaseRecorder phases;

  std::string note;
};

}  // namespace baselines
}  // namespace sjs