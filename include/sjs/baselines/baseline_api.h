#pragma once
// sjs/baselines/baseline_api.h
//
// Unified baseline interface for Spatial Join Sampling experiments.
//
// This project evaluates multiple spatial-join baselines under three variants:
//   (1) Sampling
//   (2) Enumerate + Sampling
//   (3) Adaptive (threshold-based switching)
//
// To ensure experimental consistency, all baselines expose the same API:
//   - Reset()      : reset internal state (does not necessarily free memory)
//   - Build(...)   : preprocess / build indexes
//   - Count(...)   : compute join cardinality (exact or estimated)
//   - Sample(...)  : draw join samples (optionally weighted)
//   - Enumerate(...) : return a deterministic enumerator of all join pairs
//
// Runners (see baselines/runners/*) orchestrate the calls above and record
// timing breakdown via sjs::PhaseRecorder.

#include "sjs/core/types.h"
#include "sjs/core/config.h"
#include "sjs/core/rng.h"
#include "sjs/core/timer.h"
#include "sjs/io/dataset.h"
#include "sjs/join/join_types.h"

#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {

// --------------------------
// CountResult
// --------------------------
// Represents an exact or estimated join cardinality.
// - value: the best estimate of |J|
// - exact: true if value is exact
// - stderr/ci_*: optional uncertainty information (algorithms may leave as NaN)
struct CountResult {
  long double value = 0.0L;
  bool exact = false;

  long double stderr = std::numeric_limits<long double>::quiet_NaN();
  long double ci_low = std::numeric_limits<long double>::quiet_NaN();
  long double ci_high = std::numeric_limits<long double>::quiet_NaN();

  // Optional: number of random draws / probes used for this estimate.
  u64 aux_draws = 0;

  bool HasStdErr() const noexcept { return std::isfinite(static_cast<double>(stderr)); }
  bool HasCI() const noexcept {
    return std::isfinite(static_cast<double>(ci_low)) && std::isfinite(static_cast<double>(ci_high));
  }

  // Convert to u64 (rounded). Negative values clamp to 0.
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
// A set of sampled join pairs.
// - pairs: sampled (r_id, s_id)
// - weights: optional per-sample weights (for importance sampling)
//
// Conventions:
// - If weighted == false, weights should be empty.
// - If weighted == true, weights.size() must equal pairs.size().
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
// Deterministic join enumerator with optional stats.
//
// Baselines implementing Enumerate() should return an object whose Reset() puts
// the stream back to the start and whose Next() enumerates intersecting pairs
// in a deterministic order across runs.
class IJoinEnumerator : public join::IJoinStream {
 public:
  ~IJoinEnumerator() override = default;

  // May return a zeroed stats struct if the enumerator does not collect stats.
  virtual const join::JoinStats& Stats() const noexcept = 0;
};

// --------------------------
// IBaseline
// --------------------------
// Templated by Dim for compile-time dimensionality.
// Runners are also templated by Dim and will instantiate for Dim=2 now.
//
// Semantics:
//  - Reset() can be called between repeats (different seeds) without rebuilding
//    the dataset, but must leave the instance in a state where Build() and
//    subsequent Count/Sample/Enumerate calls are valid.
//  - Build() may cache pointers to the dataset relations; dataset must outlive
//    the baseline object.
//  - Count() and Sample() should use cfg.run.* and cfg.run.extra for knobs.
//  - Enumerate() should return a fresh enumerator (or a reset one) ready to
//    produce pairs from the start.
//
// Error handling:
//  - On failure, methods return false and set *err.
//  - Unsupported operations should also return false with a clear message.

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
                     const Config& cfg,
                     PhaseRecorder* phases,
                     std::string* err) = 0;

  virtual bool Count(const Config& cfg,
                     Rng* rng,
                     CountResult* out,
                     PhaseRecorder* phases,
                     std::string* err) = 0;

  virtual bool Sample(const Config& cfg,
                      Rng* rng,
                      SampleSet* out,
                      PhaseRecorder* phases,
                      std::string* err) = 0;

  virtual std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                                     PhaseRecorder* phases,
                                                     std::string* err) = 0;
};

// --------------------------
// Runner-facing report (optional convenience container)
// --------------------------
// Runners in baselines/runners fill this struct for downstream CSV export.
struct RunReport {
  bool ok = false;
  std::string error;

  Method method = Method::Unknown;
  Variant variant = Variant::Sampling;
  std::string baseline_name;
  std::string dataset_name;

  u64 seed = 0;
  u64 t = 0;

  CountResult count;
  SampleSet samples;

  // Enumeration-related metadata (used by EnumSampling / Adaptive)
  bool used_enumeration = false;
  bool enumeration_truncated = false;
  u64 enumeration_cap = 0;
  u64 enumeration_pairs_pass1 = 0;
  u64 enumeration_pairs_pass2 = 0;
  join::JoinStats enum_stats_pass1{};
  join::JoinStats enum_stats_pass2{};

  // Adaptive decision metadata
  std::string adaptive_branch;  // e.g., "enumerate_all" or "fallback_sampling"
  u64 adaptive_pilot_pairs = 0;

  // Timing breakdown
  PhaseRecorder phases;

  std::string note;
};

}  // namespace baselines
}  // namespace sjs
