#pragma once
// sjs/sampling/alias_table.h
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

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/core/rng.h"

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

  void Clear() {
    prob_.clear();
    alias_.clear();
    n_ = 0;
    total_weight_ = 0.0;
    uniform_fallback_ = false;
  }

  usize Size() const noexcept { return n_; }
  bool Empty() const noexcept { return n_ == 0; }
  double TotalWeight() const noexcept { return total_weight_; }

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

    prob_.assign(n_, 1.0);
    alias_.assign(n_, 0);

    if (!(sum > 0.0L)) {
      // All weights are 0 -> uniform.
      uniform_fallback_ = true;
      for (usize i = 0; i < n_; ++i) {
        prob_[i] = 1.0;
        alias_[i] = static_cast<u32>(i);
      }
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

      prob_[s] = static_cast<double>(scaled[s]);   // in (0,1)
      alias_[s] = l;

      // Reduce l by the deficit of s.
      scaled[l] = (scaled[l] + scaled[s]) - 1.0L;

      if (scaled[l] < 1.0L) small.push_back(l);
      else large.push_back(l);
    }

    // Remaining buckets get probability 1.
    for (u32 idx : large) {
      prob_[idx] = 1.0;
      alias_[idx] = idx;
    }
    for (u32 idx : small) {
      prob_[idx] = 1.0;
      alias_[idx] = idx;
    }

    // Defensive clamps.
    for (usize i = 0; i < n_; ++i) {
      if (prob_[i] < 0.0) prob_[i] = 0.0;
      if (prob_[i] > 1.0) prob_[i] = 1.0;
      if (alias_[i] >= n_) alias_[i] = static_cast<u32>(i);
    }

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

    prob_.assign(n_, 1.0);
    alias_.resize(n_);

    if (sum128 == 0) {
      // uniform
      for (usize i = 0; i < n_; ++i) alias_[i] = static_cast<u32>(i);
      uniform_fallback_ = true;
      total_weight_ = 0.0;
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

      prob_[s] = static_cast<double>(scaled[s]);
      alias_[s] = l;

      scaled[l] = (scaled[l] + scaled[s]) - 1.0L;

      if (scaled[l] < 1.0L) small.push_back(l);
      else large.push_back(l);
    }

    for (u32 idx : large) { prob_[idx] = 1.0; alias_[idx] = idx; }
    for (u32 idx : small) { prob_[idx] = 1.0; alias_[idx] = idx; }

    for (usize i = 0; i < n_; ++i) {
      if (prob_[i] < 0.0) prob_[i] = 0.0;
      if (prob_[i] > 1.0) prob_[i] = 1.0;
      if (alias_[i] >= n_) alias_[i] = static_cast<u32>(i);
    }
    return true;
  }

  // Sample an index according to the built distribution.
  // Precondition: Size() > 0 and Build() has been called.
  usize Sample(Rng* rng) const noexcept {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(n_ > 0);

    if (uniform_fallback_) {
      return static_cast<usize>(rng->UniformU64(static_cast<u64>(n_)));
    }

    const usize i = static_cast<usize>(rng->UniformU64(static_cast<u64>(n_)));
    const double u = rng->NextDouble();
    // u in [0,1). Choose i with prob_[i], else alias.
    return (u < prob_[i]) ? i : static_cast<usize>(alias_[i]);
  }

  // NOTE: Alias table does not store the original normalized probabilities; it stores
  // per-bucket thresholds used for sampling. Keep your original weights if you need
  // exact probability mass values.
  double BucketProbThreshold(usize i) const noexcept {
    SJS_DASSERT(i < n_);
    return prob_[i];
  }

 private:
  std::vector<double> prob_;
  std::vector<u32> alias_;
  usize n_{0};
  double total_weight_{0.0};
  bool uniform_fallback_{false};
};

}  // namespace sampling
}  // namespace sjs
