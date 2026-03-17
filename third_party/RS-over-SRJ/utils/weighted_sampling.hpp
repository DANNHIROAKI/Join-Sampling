#pragma once
// utils/weighted_sampling.hpp (modified)
//
// RS-over-SRJ uses Walker's alias method for weighted sampling.
//
// In the upstream repository, the alias table is implemented locally and uses
// its own RNG stream (pcg32). In the main repository, we already have a
// project-wide alias table (sjs::sampling::AliasTable) and a project-wide RNG
// (sjs::Rng).
//
// For reproducibility across baselines (same seed => same sampling stream) and
// to avoid mixing RNG sources, we **default** to wrapping the project alias
// implementation.
//
// To address fairness / auditability concerns, we also provide an optional
// header-only "upstream-like" Walker alias implementation that is independent
// from the rest of the project alias code.
//
// Choose implementation:
//   - default (recommended):   RS_OVER_SRJ_USE_UPSTREAM_ALIAS = 0
//   - upstream-like (Walker):  RS_OVER_SRJ_USE_UPSTREAM_ALIAS = 1
//
// Both implementations:
//   - accept u64 weights (zeros allowed, but at least one weight must be > 0)
//   - sample using the provided sjs::Rng* (no hidden RNG state)

#include "core/rng.h"
#include "core/types.h"

#include <string>
#include <vector>

#ifndef RS_OVER_SRJ_USE_UPSTREAM_ALIAS
#define RS_OVER_SRJ_USE_UPSTREAM_ALIAS 0
#endif

#if RS_OVER_SRJ_USE_UPSTREAM_ALIAS

// -----------------------------------------------------------------------------
// Upstream-like Walker alias method (self-contained)
// -----------------------------------------------------------------------------

class alias {
 public:
  alias() = default;

  explicit alias(const std::vector<sjs::u64>& weights) { (void)Build(weights); }

  void Clear() {
    prob_.clear();
    alias_.clear();
  }

  sjs::usize Size() const noexcept { return prob_.size(); }
  bool Empty() const noexcept { return prob_.empty(); }

  bool Build(const std::vector<sjs::u64>& weights, std::string* err = nullptr) {
    Clear();
    if (weights.empty()) {
      if (err) *err = "alias::Build: empty weight array";
      return false;
    }

    // Sum weights; require at least one positive weight.
    long double sum = 0.0L;
    for (sjs::u64 w : weights) sum += static_cast<long double>(w);
    if (!(sum > 0.0L)) {
      if (err) *err = "alias::Build: all weights are zero";
      return false;
    }

    const sjs::usize n = weights.size();
    prob_.assign(n, 0.0);
    alias_.assign(n, 0);

    // Scale to average 1.0.
    std::vector<double> scaled(n);
    for (sjs::usize i = 0; i < n; ++i) {
      scaled[i] = static_cast<double>(static_cast<long double>(weights[i]) * static_cast<long double>(n) / sum);
    }

    std::vector<sjs::u32> small;
    std::vector<sjs::u32> large;
    small.reserve(n);
    large.reserve(n);

    for (sjs::usize i = 0; i < n; ++i) {
      if (scaled[i] < 1.0) {
        small.push_back(static_cast<sjs::u32>(i));
      } else {
        large.push_back(static_cast<sjs::u32>(i));
      }
    }

    while (!small.empty() && !large.empty()) {
      const sjs::u32 s = small.back();
      small.pop_back();
      const sjs::u32 l = large.back();
      large.pop_back();

      prob_[s] = scaled[s];
      alias_[s] = l;

      scaled[l] = (scaled[l] + scaled[s]) - 1.0;
      if (scaled[l] < 1.0) {
        small.push_back(l);
      } else {
        large.push_back(l);
      }
    }

    // Remaining entries have probability 1.
    for (sjs::u32 i : large) {
      prob_[i] = 1.0;
      alias_[i] = i;
    }
    for (sjs::u32 i : small) {
      prob_[i] = 1.0;
      alias_[i] = i;
    }

    return true;
  }

  // Sample one index using the provided RNG.
  sjs::usize Sample(sjs::Rng* rng) const noexcept {
    const sjs::u32 n = static_cast<sjs::u32>(prob_.size());
    if (n == 0) return 0;
    const sjs::u32 i = rng->UniformU32(n);
    const double u = rng->NextDouble();
    return (u < prob_[i]) ? static_cast<sjs::usize>(i) : static_cast<sjs::usize>(alias_[i]);
  }

 private:
  std::vector<double> prob_;
  std::vector<sjs::u32> alias_;
};

#else

// -----------------------------------------------------------------------------
// Project alias wrapper (default)
// -----------------------------------------------------------------------------

#include "sampling/alias_table.h"

// Keep the original class name `alias` to minimize downstream edits.
class alias {
 public:
  alias() = default;

  explicit alias(const std::vector<sjs::u64>& weights) { (void)Build(weights); }

  void Clear() { table_.Clear(); }

  sjs::usize Size() const noexcept { return table_.Size(); }
  bool Empty() const noexcept { return table_.Empty(); }

  // Build from u64 weights.
  bool Build(const std::vector<sjs::u64>& weights, std::string* err = nullptr) {
    return table_.BuildFromU64(sjs::Span<const sjs::u64>(weights), err);
  }

  // Sample one index using the provided RNG.
  sjs::usize Sample(sjs::Rng* rng) const noexcept { return table_.Sample(rng); }

 private:
  sjs::sampling::AliasTable table_;
};

#endif
