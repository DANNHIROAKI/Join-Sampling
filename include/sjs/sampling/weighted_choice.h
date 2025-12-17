#pragma once
// sjs/sampling/weighted_choice.h
//
// Utilities for sampling from discrete weighted sets.
//
// This file complements alias_table.h:
//  - WeightedChoiceLinear: O(n) one-off sampling (no preprocessing).
//  - PrefixDistribution: O(n) build + O(log n) sampling via prefix sums.
//
// These utilities are frequently useful inside baselines where weights
// change between iterations (alias tables would be too expensive to rebuild).

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/core/rng.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace sjs {
namespace sampling {

// ---------- One-off (O(n)) weighted choice ----------
//
// Returns an index in [0,weights.size()) proportional to weights[i].
//
// Behavior:
//  - Negative weights are rejected (returns kInvalidUsize and sets err).
//  - If all weights are zero, falls back to uniform over indices.
//  - For floating weights, NaN/Inf are rejected.
//
// This function is convenient when you only need 1-2 samples, or weights
// change every time.
inline constexpr usize kInvalidUsize = static_cast<usize>(-1);

namespace detail {

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

inline bool IsFiniteNonNeg(double w) {
  return (w >= 0.0) && std::isfinite(w);
}

}  // namespace detail

// Double weights version.
inline usize WeightedChoiceLinear(Span<const double> weights, Rng* rng, std::string* err = nullptr) {
  SJS_ASSERT(rng != nullptr);
  const usize n = weights.size();
  if (n == 0) {
    detail::SetErr(err, "WeightedChoiceLinear: empty weights");
    return kInvalidUsize;
  }

  long double sum = 0.0L;
  for (usize i = 0; i < n; ++i) {
    const double w = weights[i];
    if (!detail::IsFiniteNonNeg(w)) {
      detail::SetErr(err, "WeightedChoiceLinear: weight must be finite and >= 0");
      return kInvalidUsize;
    }
    sum += static_cast<long double>(w);
  }

  if (!(sum > 0.0L)) {
    // Uniform fallback.
    return static_cast<usize>(rng->UniformU64(static_cast<u64>(n)));
  }

  // r in [0,sum)
  const long double r = static_cast<long double>(rng->NextDouble()) * sum;
  long double acc = 0.0L;
  for (usize i = 0; i < n; ++i) {
    acc += static_cast<long double>(weights[i]);
    if (r < acc) return i;
  }
  // Due to rounding, might fall through; return last.
  return n - 1;
}

// Integer weights version.
inline usize WeightedChoiceLinear(Span<const u64> weights, Rng* rng, std::string* err = nullptr) {
  SJS_ASSERT(rng != nullptr);
  const usize n = weights.size();
  if (n == 0) {
    detail::SetErr(err, "WeightedChoiceLinear(u64): empty weights");
    return kInvalidUsize;
  }

  __uint128_t sum128 = 0;
  for (usize i = 0; i < n; ++i) sum128 += static_cast<__uint128_t>(weights[i]);

  if (sum128 == 0) {
    return static_cast<usize>(rng->UniformU64(static_cast<u64>(n)));
  }

  // If sum fits in u64, use integer sampling for exactness.
  if (sum128 <= static_cast<__uint128_t>(std::numeric_limits<u64>::max())) {
    const u64 sum = static_cast<u64>(sum128);
    const u64 r = rng->UniformU64(sum);  // [0,sum)
    __uint128_t acc = 0;
    for (usize i = 0; i < n; ++i) {
      acc += static_cast<__uint128_t>(weights[i]);
      if (static_cast<__uint128_t>(r) < acc) return i;
    }
    return n - 1;
  }

  // Otherwise fall back to long double.
  const long double sum = static_cast<long double>(sum128);
  const long double r = static_cast<long double>(rng->NextDouble()) * sum;
  long double acc = 0.0L;
  for (usize i = 0; i < n; ++i) {
    acc += static_cast<long double>(weights[i]);
    if (r < acc) return i;
  }
  return n - 1;
}

// ---------- Prefix-sum distribution (O(log n) samples) ----------
class PrefixDistribution {
 public:
  PrefixDistribution() = default;

  void Clear() {
    prefix_.clear();
    total_ = 0.0L;
    uniform_fallback_ = false;
  }

  usize Size() const noexcept { return prefix_.size(); }
  bool Empty() const noexcept { return prefix_.empty(); }

  long double TotalWeight() const noexcept { return total_; }

  // Build from double weights.
  bool Build(Span<const double> weights, std::string* err = nullptr) {
    Clear();
    const usize n = weights.size();
    prefix_.resize(n);
    long double sum = 0.0L;

    for (usize i = 0; i < n; ++i) {
      const double w = weights[i];
      if (!detail::IsFiniteNonNeg(w)) {
        detail::SetErr(err, "PrefixDistribution::Build: weight must be finite and >=0");
        Clear();
        return false;
      }
      sum += static_cast<long double>(w);
      prefix_[i] = sum;
    }

    total_ = sum;
    uniform_fallback_ = !(sum > 0.0L);
    return true;
  }

  // Build from u64 weights.
  bool BuildFromU64(Span<const u64> weights, std::string* err = nullptr) {
    Clear();
    const usize n = weights.size();
    prefix_.resize(n);

    __uint128_t sum128 = 0;
    for (usize i = 0; i < n; ++i) {
      sum128 += static_cast<__uint128_t>(weights[i]);
      prefix_[i] = static_cast<long double>(sum128);
    }
    total_ = static_cast<long double>(sum128);
    uniform_fallback_ = (sum128 == 0);
    return true;
  }

  // Sample an index. If all weights are zero, falls back to uniform.
  usize Sample(Rng* rng) const noexcept {
    SJS_DASSERT(rng != nullptr);
    const usize n = prefix_.size();
    SJS_DASSERT(n > 0);

    if (uniform_fallback_) {
      return static_cast<usize>(rng->UniformU64(static_cast<u64>(n)));
    }

    const long double r = static_cast<long double>(rng->NextDouble()) * total_;
    // Find first prefix > r.
    auto it = std::upper_bound(prefix_.begin(), prefix_.end(), r);
    if (it == prefix_.end()) return n - 1;
    return static_cast<usize>(it - prefix_.begin());
  }

 private:
  std::vector<long double> prefix_;
  long double total_{0.0L};
  bool uniform_fallback_{false};
};

}  // namespace sampling
}  // namespace sjs
