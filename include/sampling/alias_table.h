#pragma once
// sampling/alias_table.h
//
// Alias table (Vose's alias method) for O(1) sampling from a discrete distribution.
//
// This is useful when you need to draw many samples from a fixed weight vector.
//
// API:
//  - Build(weights): preprocess O(n)
//  - Sample(rng): draw one index in O(1)
//
// Notes:
//  - Accepts non-negative weights. If all weights sum to 0, it falls back to uniform.
//  - Stores per-bucket probability threshold in [0,1] and alias index.
//  - Deterministic given input weights and RNG seed.
//
// References:
//  - A. J. Walker (1974), "New fast method for generating discrete random numbers..."
//  - M. D. Vose (1991), "A linear algorithm for generating random numbers with a given distribution"

#include "core/types.h"
#include "core/assert.h"
#include "core/rng.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace sampling {

class AliasTable {
 public:
  AliasTable() = default;

  // Constants for integer thresholds.
  static constexpr u32 kProbOne = std::numeric_limits<u32>::max();
  static constexpr u64 kProbScaleU64 = 4294967296ULL;         // 2^32
  static constexpr long double kProbScaleLd = 4294967296.0L;  // 2^32

  void Clear() {
    prob_.clear();
    alias_.clear();
    n_ = 0;
    total_weight_ = 0.0;
    uniform_fallback_ = false;
    bound_ = 0;
    lemire_threshold_ = 0;
  }

  usize Size() const noexcept { return n_; }
  bool Empty() const noexcept { return n_ == 0; }
  double TotalWeight() const noexcept { return total_weight_; }

  // For fixed-size distributions, sampling requires drawing:
  //  (1) i ~ Uniform{0..n-1}
  //  (2) u ~ Uniform[0,1)
  // The project-wide RNG already provides fast u64 draws, so we implement
  // Lemire's unbiased bounded sampling *inside* the alias table and precompute
  // the rejection threshold once per Build() to avoid an expensive modulo per
  // sample.
  //
  // References:
  //   - Lemire (2019): Fast Random Integer Generation in an Interval
  //     https://arxiv.org/abs/1805.10941
  u32 SampleUniformIndex(Rng* rng) const noexcept {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(n_ > 0);

#if defined(__SIZEOF_INT128__)
    // Unbiased in [0, bound_).
    while (true) {
      const u64 x = rng->NextU64();
      const __uint128_t m = static_cast<__uint128_t>(x) * static_cast<__uint128_t>(bound_);
      const u64 l = static_cast<u64>(m);
      if (l >= lemire_threshold_) {
        return static_cast<u32>(m >> 64);
      }
    }
#else
    return static_cast<u32>(rng->UniformU64(bound_));
#endif
  }

  // Build from double weights.
  // Returns false on invalid (negative/NaN/Inf) weights.
  bool Build(Span<const double> weights, std::string* err = nullptr) {
    Clear();
    n_ = weights.size();
    if (n_ == 0) return true;
    if (n_ > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "AliasTable::Build: n too large for u32 alias indices";
      Clear();
      return false;
    }

    // Validate + sum in long double for robustness.
    long double sum = 0.0L;
    for (usize i = 0; i < n_; ++i) {
      const double w = weights[i];
      if (!(w >= 0.0) || !std::isfinite(w)) {
        if (err) *err = "AliasTable::Build: weight must be finite and >= 0";
        Clear();
        return false;
      }
      sum += static_cast<long double>(w);
    }
    total_weight_ = static_cast<double>(sum);

    prob_.assign(n_, kProbOne);
    alias_.assign(n_, 0);

    if (!(sum > 0.0L)) {
      // All weights are 0 -> uniform.
      uniform_fallback_ = true;
      for (usize i = 0; i < n_; ++i) {
        prob_[i] = kProbOne;
        alias_[i] = static_cast<u32>(i);
      }

      bound_ = static_cast<u64>(n_);
      lemire_threshold_ = (bound_ == 0) ? 0 : (static_cast<u64>(-bound_) % bound_);
      return true;
    }

    uniform_fallback_ = false;

    // Scaled probabilities: p_i = w_i * n / sum, so avg is 1.
    std::vector<long double> scaled(n_);
    for (usize i = 0; i < n_; ++i) {
      scaled[i] = static_cast<long double>(weights[i]) * static_cast<long double>(n_) / sum;
    }

    std::vector<u32> small;
    std::vector<u32> large;
    small.reserve(n_);
    large.reserve(n_);

    for (usize i = 0; i < n_; ++i) {
      if (scaled[i] < 1.0L) small.push_back(static_cast<u32>(i));
      else large.push_back(static_cast<u32>(i));
    }

    // Main construction.
    while (!small.empty() && !large.empty()) {
      const u32 s = small.back(); small.pop_back();
      const u32 l = large.back(); large.pop_back();

      // Convert p in (0,1) to an integer threshold in [0, 2^32).
      // We store thresholds as u32, with a special sentinel for p==1.
      {
        const long double p = scaled[s];
        if (!(p > 0.0L)) {
          prob_[s] = 0;
        } else if (p >= 1.0L) {
          prob_[s] = kProbOne;
        } else {
          const long double x = p * kProbScaleLd;
          u64 thr = static_cast<u64>(x);
          if (thr >= kProbScaleU64) thr = kProbScaleU64 - 1;
          prob_[s] = static_cast<u32>(thr);
        }
      }
      alias_[s] = l;

      // Reduce l by the deficit of s.
      scaled[l] = (scaled[l] + scaled[s]) - 1.0L;

      if (scaled[l] < 1.0L) small.push_back(l);
      else large.push_back(l);
    }

    // Remaining buckets get probability 1.
    for (u32 idx : large) {
      prob_[idx] = kProbOne;
      alias_[idx] = idx;
    }
    for (u32 idx : small) {
      prob_[idx] = kProbOne;
      alias_[idx] = idx;
    }

    // Defensive clamps.
    for (usize i = 0; i < n_; ++i) {
      if (alias_[i] >= n_) alias_[i] = static_cast<u32>(i);
    }

    bound_ = static_cast<u64>(n_);
    lemire_threshold_ = (bound_ == 0) ? 0 : (static_cast<u64>(-bound_) % bound_);

    return true;
  }

  // Build from integer weights (u64).
  bool BuildFromU64(Span<const u64> weights, std::string* err = nullptr) {
    Clear();
    n_ = weights.size();
    if (n_ == 0) return true;
    if (n_ > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "AliasTable::BuildFromU64: n too large for u32 alias indices";
      Clear();
      return false;
    }

    __uint128_t sum128 = 0;
    for (usize i = 0; i < n_; ++i) sum128 += static_cast<__uint128_t>(weights[i]);

    prob_.assign(n_, kProbOne);
    alias_.resize(n_);

    if (sum128 == 0) {
      // uniform
      for (usize i = 0; i < n_; ++i) alias_[i] = static_cast<u32>(i);
      uniform_fallback_ = true;
      total_weight_ = 0.0;

      bound_ = static_cast<u64>(n_);
      lemire_threshold_ = (bound_ == 0) ? 0 : (static_cast<u64>(-bound_) % bound_);
      return true;
    }

    // Convert to long double for scaled computation.
    const long double sum = static_cast<long double>(sum128);
    total_weight_ = static_cast<double>(sum);
    uniform_fallback_ = false;

    std::vector<long double> scaled(n_);
    std::vector<u32> small;
    std::vector<u32> large;
    small.reserve(n_);
    large.reserve(n_);

    const long double n_ld = static_cast<long double>(n_);
    for (usize i = 0; i < n_; ++i) {
      scaled[i] = static_cast<long double>(weights[i]) * n_ld / sum;
      if (scaled[i] < 1.0L) small.push_back(static_cast<u32>(i));
      else large.push_back(static_cast<u32>(i));
    }

    while (!small.empty() && !large.empty()) {
      const u32 s = small.back(); small.pop_back();
      const u32 l = large.back(); large.pop_back();

      {
        const long double p = scaled[s];
        if (!(p > 0.0L)) {
          prob_[s] = 0;
        } else if (p >= 1.0L) {
          prob_[s] = kProbOne;
        } else {
          const long double x = p * kProbScaleLd;
          u64 thr = static_cast<u64>(x);
          if (thr >= kProbScaleU64) thr = kProbScaleU64 - 1;
          prob_[s] = static_cast<u32>(thr);
        }
      }
      alias_[s] = l;

      scaled[l] = (scaled[l] + scaled[s]) - 1.0L;

      if (scaled[l] < 1.0L) small.push_back(l);
      else large.push_back(l);
    }

    for (u32 idx : large) { prob_[idx] = kProbOne; alias_[idx] = idx; }
    for (u32 idx : small) { prob_[idx] = kProbOne; alias_[idx] = idx; }

    for (usize i = 0; i < n_; ++i) {
      if (alias_[i] >= n_) alias_[i] = static_cast<u32>(i);
    }

    bound_ = static_cast<u64>(n_);
    lemire_threshold_ = (bound_ == 0) ? 0 : (static_cast<u64>(-bound_) % bound_);
    return true;
  }

  // Sample an index according to the built distribution.
  // Precondition: Size() > 0 and Build() has been called.
  usize Sample(Rng* rng) const noexcept {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(n_ > 0);

    const u32 i = SampleUniformIndex(rng);
    if (uniform_fallback_) {
      return static_cast<usize>(i);
    }

    const u32 thr = prob_[static_cast<usize>(i)];
    if (thr == kProbOne) return static_cast<usize>(i);

    // u32 uniform in [0, 2^32).
    const u32 u = rng->NextU32();
    return (u < thr) ? static_cast<usize>(i) : static_cast<usize>(alias_[static_cast<usize>(i)]);
  }

  // NOTE: Alias table does not store the original normalized probabilities; it stores
  // per-bucket thresholds used for sampling. Keep your original weights if you need
  // exact probability mass values.
  double BucketProbThreshold(usize i) const noexcept {
    SJS_DASSERT(i < n_);
    const u32 thr = prob_[i];
    if (thr == kProbOne) return 1.0;
    return static_cast<double>(thr) * (1.0 / 4294967296.0);  // 2^32
  }

 private:
  // Integer thresholds for the alias coin flip.
  // For bucket i:
  //   choose i if u < prob_[i], else alias_[i]
  // where u is uniform on [0, 2^32).
  //
  // prob_[i] == kProbOne means probability 1 (always choose i).
  std::vector<u32> prob_;
  std::vector<u32> alias_;
  usize n_{0};
  double total_weight_{0.0};
  bool uniform_fallback_{false};

  // Precomputed for Lemire bounded sampling of the uniform bucket index.
  u64 bound_{0};
  u64 lemire_threshold_{0};

};

// -----------------------------------------------------------------------------
// AliasSmall<N>: stack-based alias table for small N (no heap allocation)
//
// Intended for hot inner loops where (1) the distribution has very small support
// (e.g., O(log m) segment tree cover nodes) and (2) we need to draw many samples
// from the same weights.
//
// This mirrors AliasTable but stores all arrays on the stack (std::array) and
// only supports BuildFromU64.
// -----------------------------------------------------------------------------

template <usize MaxN>
class AliasSmall {
 public:

  static constexpr u32 kProbOne = std::numeric_limits<u32>::max();
  static constexpr u64 kProbScaleU64 = 4294967296ULL;         // 2^32
  static constexpr long double kProbScaleLd = 4294967296.0L;  // 2^32
  AliasSmall() = default;

  void Clear() {
    n_ = 0;
    uniform_fallback_ = false;
  }

  u32 Size() const noexcept { return n_; }
  bool Empty() const noexcept { return n_ == 0; }

  bool BuildFromU64(Span<const u64> weights, std::string* err = nullptr) {
    Clear();
    const usize n = weights.size();
    if (n == 0) return true;
    if (n > MaxN) {
      if (err) *err = "AliasSmall::BuildFromU64: n exceeds MaxN";
      return false;
    }
    if (n > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "AliasSmall::BuildFromU64: n too large for u32 alias indices";
      return false;
    }

    n_ = static_cast<u32>(n);

    __uint128_t sum128 = 0;
    for (usize i = 0; i < n; ++i) sum128 += static_cast<__uint128_t>(weights[i]);

    if (sum128 == 0) {
      // Uniform fallback (all weights 0).
      uniform_fallback_ = true;
      for (u32 i = 0; i < n_; ++i) {
        prob_[i] = kProbOne;
        alias_[i] = i;
      }

      bound_ = static_cast<u64>(n_);
      lemire_threshold_ = (bound_ == 0) ? 0 : (static_cast<u64>(-bound_) % bound_);
      return true;
    }

    uniform_fallback_ = false;

    const long double sum = static_cast<long double>(sum128);
    const long double n_ld = static_cast<long double>(n_);

    std::array<long double, MaxN> scaled;

    std::array<u32, MaxN> small;
    std::array<u32, MaxN> large;
    u32 small_sz = 0;
    u32 large_sz = 0;

    for (u32 i = 0; i < n_; ++i) {
      scaled[i] = static_cast<long double>(weights[static_cast<usize>(i)]) * n_ld / sum;
      if (scaled[i] < 1.0L) {
        small[small_sz++] = i;
      } else {
        large[large_sz++] = i;
      }
    }

    while (small_sz > 0 && large_sz > 0) {
      const u32 s = small[--small_sz];
      const u32 l = large[--large_sz];

      {
        const long double p = scaled[s];
        if (!(p > 0.0L)) {
          prob_[s] = 0;
        } else if (p >= 1.0L) {
          prob_[s] = kProbOne;
        } else {
          const long double x = p * kProbScaleLd;
          u64 thr = static_cast<u64>(x);
          if (thr >= kProbScaleU64) thr = kProbScaleU64 - 1;
          prob_[s] = static_cast<u32>(thr);
        }
      }
      alias_[s] = l;

      scaled[l] = (scaled[l] + scaled[s]) - 1.0L;

      if (scaled[l] < 1.0L) small[small_sz++] = l;
      else large[large_sz++] = l;
    }

    while (large_sz > 0) {
      const u32 idx = large[--large_sz];
      prob_[idx] = kProbOne;
      alias_[idx] = idx;
    }
    while (small_sz > 0) {
      const u32 idx = small[--small_sz];
      prob_[idx] = kProbOne;
      alias_[idx] = idx;
    }

    // Defensive clamps.
    for (u32 i = 0; i < n_; ++i) {
      if (alias_[i] >= n_) alias_[i] = i;
    }

    bound_ = static_cast<u64>(n_);
    lemire_threshold_ = (bound_ == 0) ? 0 : (static_cast<u64>(-bound_) % bound_);

    return true;
  }

  // Sample an index according to the built distribution.
  // Precondition: Size() > 0 and BuildFromU64() has been called.
  u32 Sample(Rng* rng) const noexcept {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(n_ > 0);

    // Sample bucket i ~ Uniform{0..n-1} using precomputed Lemire threshold.
#if defined(__SIZEOF_INT128__)
    u32 i;
    while (true) {
      const u64 x = rng->NextU64();
      const __uint128_t m = static_cast<__uint128_t>(x) * static_cast<__uint128_t>(bound_);
      const u64 l = static_cast<u64>(m);
      if (l >= lemire_threshold_) {
        i = static_cast<u32>(m >> 64);
        break;
      }
    }
#else
    const u32 i = rng->UniformU32(n_);
#endif

    if (uniform_fallback_) return i;
    const u32 thr = prob_[i];
    if (thr == kProbOne) return i;
    const u32 u = rng->NextU32();
    return (u < thr) ? i : alias_[i];
  }

 private:
  std::array<u32, MaxN> prob_{};
  std::array<u32, MaxN> alias_{};
  u32 n_{0};
  bool uniform_fallback_{false};

  // Precomputed for Lemire bounded sampling.
  u64 bound_{0};
  u64 lemire_threshold_{0};
};


}  // namespace sampling
}  // namespace sjs
