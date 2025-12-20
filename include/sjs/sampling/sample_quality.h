#pragma once
// sjs/sampling/sample_quality.h
//
// Statistical diagnostics for sampling quality on small datasets.
//
// Intended usage:
//  - For small-scale validation, compute the exact join universe U (oracle).
//  - Draw samples from an algorithm.
//  - Evaluate whether samples look uniform over U and independent across time.
//
// Provided metrics (lightweight, dependency-free):
//  - Chi-square goodness-of-fit to uniform discrete distribution.
//  - Kolmogorov–Smirnov (KS) test for continuous samples in [0,1) or two-sample KS.
//  - Autocorrelation (lag-k) on a numeric representation of samples.
//  - Convenience evaluator for PairId samples vs oracle universe.
//
// Notes:
//  - These tests have assumptions; for rigorous p-values you may need larger sample sizes.
//  - We implement common approximations (Numerical Recipes style) without external libs.

#include "sjs/core/types.h"
#include "sjs/core/assert.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sjs {
namespace sampling {
namespace quality {

namespace detail {

inline void Clamp01(double* p) {
  if (!p) return;
  if (*p < 0.0) *p = 0.0;
  if (*p > 1.0) *p = 1.0;
}

// --------------------------
// Regularized upper incomplete gamma Q(a,x)
// Q(a,x) = Γ(a,x) / Γ(a)
// Used for chi-square p-values: p = Q(df/2, chi2/2)
//
// Implementation adapted from Numerical Recipes (series + continued fraction).
// --------------------------
inline double GammaQ(double a, double x) {
  // Domain checks.
  if (!(a > 0.0) || x < 0.0 || !std::isfinite(a) || !std::isfinite(x)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (x == 0.0) return 1.0;

  constexpr int ITMAX = 200;
  constexpr double EPS = 3.0e-10;
  constexpr double FPMIN = 1.0e-300;

  const double gln = std::lgamma(a);

  if (x < a + 1.0) {
    // Use series for P(a,x), then Q=1-P
    double ap = a;
    double del = 1.0 / a;
    double sum = del;
    for (int n = 1; n <= ITMAX; ++n) {
      ap += 1.0;
      del *= x / ap;
      sum += del;
      if (std::fabs(del) < std::fabs(sum) * EPS) break;
    }
    const double gamser = sum * std::exp(-x + a * std::log(x) - gln);
    double q = 1.0 - gamser;
    if (q < 0.0) q = 0.0;
    if (q > 1.0) q = 1.0;
    return q;
  } else {
    // Use continued fraction for Q(a,x)
    double b = x + 1.0 - a;
    double c = 1.0 / FPMIN;
    double d = 1.0 / b;
    double h = d;

    for (int i = 1; i <= ITMAX; ++i) {
      const double an = -static_cast<double>(i) * (static_cast<double>(i) - a);
      b += 2.0;
      d = an * d + b;
      if (std::fabs(d) < FPMIN) d = FPMIN;
      c = b + an / c;
      if (std::fabs(c) < FPMIN) c = FPMIN;
      d = 1.0 / d;
      const double del = d * c;
      h *= del;
      if (std::fabs(del - 1.0) < EPS) break;
    }

    const double gammcf = std::exp(-x + a * std::log(x) - gln) * h;
    double q = gammcf;
    if (q < 0.0) q = 0.0;
    if (q > 1.0) q = 1.0;
    return q;
  }
}

// --------------------------
// Kolmogorov distribution approximation for KS p-value
// p ≈ 2 * Σ_{j=1..∞} (-1)^{j-1} exp(-2 j^2 λ^2)
// --------------------------
inline double KolmogorovProb(double lambda) {
  if (!(lambda >= 0.0) || !std::isfinite(lambda)) return std::numeric_limits<double>::quiet_NaN();
  if (lambda == 0.0) return 1.0;
  const double a2 = -2.0 * lambda * lambda;

  double sum = 0.0;
  double fac = 2.0;
  for (int j = 1; j <= 200; ++j) {
    const double term = fac * std::exp(a2 * static_cast<double>(j) * static_cast<double>(j));
    sum += term;
    if (std::fabs(term) < 1e-12) break;
    fac = -fac;
  }
  // Clamp.
  if (sum < 0.0) sum = 0.0;
  if (sum > 1.0) sum = 1.0;
  return sum;
}

// Deterministic 64-bit hash -> [0,1).
inline u64 SplitMix64(u64 x) noexcept {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  x = x ^ (x >> 31);
  return x;
}

inline u64 HashPair64(const PairId& p) noexcept {
  // Pack (r,s) into 64-bit (assumes Id is u32; see core/types.h).
  const u64 key = (static_cast<u64>(p.r) << 32) ^ static_cast<u64>(p.s);
  return SplitMix64(key);
}

inline double HashPairToUnit01(const PairId& p) noexcept {
  // Use 53 bits for a double mantissa.
  const u64 h = HashPair64(p);
  const u64 mant = (h >> 11);  // top 53 bits
  return static_cast<double>(mant) * (1.0 / 9007199254740992.0);  // 2^53
}

}  // namespace detail

// --------------------------
// Chi-square test (uniform discrete)
// --------------------------
struct ChiSquareResult {
  double statistic = 0.0;
  u64 df = 0;
  double p_value = std::numeric_limits<double>::quiet_NaN();  // upper tail
};

// Compute chi-square goodness-of-fit to uniform over k bins.
// obs_counts size = k. Expected count in each bin is n/k (n=sum obs).
inline ChiSquareResult ChiSquareUniform(const std::vector<u64>& obs_counts) {
  ChiSquareResult r;
  const usize k = obs_counts.size();
  if (k == 0) return r;

  long double n = 0.0L;
  for (u64 c : obs_counts) n += static_cast<long double>(c);

  if (!(n > 0.0L)) {
    r.statistic = 0.0;
    r.df = (k > 0) ? static_cast<u64>(k - 1) : 0;
    r.p_value = 1.0;
    return r;
  }

  const long double expected = n / static_cast<long double>(k);
  // If expected is extremely small, statistic becomes unstable; still compute.
  long double chi2 = 0.0L;
  for (u64 c : obs_counts) {
    const long double diff = static_cast<long double>(c) - expected;
    chi2 += (diff * diff) / expected;
  }

  r.statistic = static_cast<double>(chi2);
  r.df = (k > 0) ? static_cast<u64>(k - 1) : 0;

  // p-value = Q(df/2, chi2/2)
  if (r.df == 0) {
    r.p_value = 1.0;
  } else {
    r.p_value = detail::GammaQ(0.5 * static_cast<double>(r.df), 0.5 * r.statistic);
  }
  detail::Clamp01(&r.p_value);
  return r;
}

// --------------------------
// Kolmogorov–Smirnov tests
// --------------------------
struct KSTestResult {
  double D = 0.0;
  double p_value = std::numeric_limits<double>::quiet_NaN();
};

// One-sample KS test vs Uniform[0,1).
// Input samples should be in [0,1). Values outside are clamped for robustness.
inline KSTestResult KSOneSampleUniform01(std::vector<double> samples) {
  KSTestResult out;
  const usize n = samples.size();
  if (n == 0) {
    out.D = 0.0;
    out.p_value = 1.0;
    return out;
  }

  for (double& x : samples) {
    if (!std::isfinite(x)) x = 0.0;
    if (x < 0.0) x = 0.0;
    if (x >= 1.0) x = std::nextafter(1.0, 0.0);
  }
  std::sort(samples.begin(), samples.end());

  double D = 0.0;
  for (usize i = 0; i < n; ++i) {
    const double xi = samples[i];
    const double cdf = xi;
    const double emp_hi = static_cast<double>(i + 1) / static_cast<double>(n);
    const double emp_lo = static_cast<double>(i) / static_cast<double>(n);
    const double d1 = std::fabs(emp_hi - cdf);
    const double d2 = std::fabs(cdf - emp_lo);
    D = std::max(D, std::max(d1, d2));
  }

  out.D = D;

  const double en = std::sqrt(static_cast<double>(n));
  const double lambda = (en + 0.12 + 0.11 / en) * D;
  out.p_value = detail::KolmogorovProb(lambda);
  detail::Clamp01(&out.p_value);
  return out;
}

// Two-sample KS test (continuous).
inline KSTestResult KSTwoSample(std::vector<double> a, std::vector<double> b) {
  KSTestResult out;
  const usize n = a.size();
  const usize m = b.size();
  if (n == 0 || m == 0) {
    out.D = 0.0;
    out.p_value = 1.0;
    return out;
  }

  for (double& x : a) if (!std::isfinite(x)) x = 0.0;
  for (double& x : b) if (!std::isfinite(x)) x = 0.0;

  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());

  usize i = 0, j = 0;
  double cdf_a = 0.0, cdf_b = 0.0;
  double D = 0.0;

  while (i < n && j < m) {
    const double xa = a[i];
    const double xb = b[j];
    if (xa < xb) {
      const double x = xa;
      while (i < n && a[i] == x) ++i;
      cdf_a = static_cast<double>(i) / static_cast<double>(n);
    } else if (xb < xa) {
      const double x = xb;
      while (j < m && b[j] == x) ++j;
      cdf_b = static_cast<double>(j) / static_cast<double>(m);
    } else {
      const double x = xa;
      while (i < n && a[i] == x) ++i;
      while (j < m && b[j] == x) ++j;
      cdf_a = static_cast<double>(i) / static_cast<double>(n);
      cdf_b = static_cast<double>(j) / static_cast<double>(m);
    }
    D = std::max(D, std::fabs(cdf_a - cdf_b));
  }

  // Remaining tail(s)
  // If a has remaining values, then cdf_b == 1 and max difference is |cdf_a - 1|.
  if (i < n) D = std::max(D, std::fabs(1.0 - cdf_a));
  // If b has remaining values, then cdf_a == 1 and max difference is |cdf_b - 1|.
  if (j < m) D = std::max(D, std::fabs(1.0 - cdf_b));

  out.D = D;

  const double en = std::sqrt(static_cast<double>(n) * static_cast<double>(m) /
                              static_cast<double>(n + m));
  const double lambda = (en + 0.12 + 0.11 / en) * D;
  out.p_value = detail::KolmogorovProb(lambda);
  detail::Clamp01(&out.p_value);
  return out;
}

// --------------------------
// Autocorrelation
// --------------------------
inline double Autocorrelation(Span<const double> x, int lag = 1) {
  if (lag <= 0) return std::numeric_limits<double>::quiet_NaN();
  const usize n = x.size();
  if (n <= static_cast<usize>(lag)) return std::numeric_limits<double>::quiet_NaN();

  long double mean = 0.0L;
  for (double v : x) mean += static_cast<long double>(v);
  mean /= static_cast<long double>(n);

  long double var = 0.0L;
  for (double v : x) {
    const long double d = static_cast<long double>(v) - mean;
    var += d * d;
  }
  if (!(var > 0.0L)) return 0.0;

  long double cov = 0.0L;
  const usize m = n - static_cast<usize>(lag);
  for (usize i = 0; i < m; ++i) {
    const long double a = static_cast<long double>(x[i]) - mean;
    const long double b = static_cast<long double>(x[i + static_cast<usize>(lag)]) - mean;
    cov += a * b;
  }

  cov /= static_cast<long double>(m);
  var /= static_cast<long double>(n);

  return static_cast<double>(cov / var);
}

// Autocorrelation on PairId samples by hashing to [0,1).
inline double AutocorrelationHashedPairs(Span<const PairId> samples, int lag = 1) {
  if (samples.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x;
  x.reserve(samples.size());
  for (const auto& p : samples) {
    x.push_back(detail::HashPairToUnit01(p));
  }
  return Autocorrelation(Span<const double>(x.data(), x.size()), lag);
}

// KS sanity check for PairId samples against a known (small) oracle universe.
//
// IMPORTANT: One-sample KS vs continuous Uniform[0,1) assumes *continuous* support.
// If we simply hash PairId -> [0,1), the resulting distribution is discrete with
// at most |U| distinct values, which can lead to systematic rejection (tiny p-values)
// when the sample size t is large.
//
// To obtain a meaningful KS sanity check in EXP-1 settings, we:
//   1) Order the oracle universe deterministically by a 64-bit hash
//      (tie-break by (r,s) to make the order total).
//   2) Map each sampled pair to its rank in that order.
//   3) Add a small deterministic jitter in [0,1) based on sample index,
//      then scale by |U| to get values in [0,1).
//
// Under a truly uniform sampler over U, these values behave like Uniform[0,1).
inline KSTestResult KSPairsHashUniform01RankJitter(Span<const PairId> universe,
                                                   Span<const PairId> samples) {
  KSTestResult out;
  if (universe.empty() || samples.empty()) {
    out.D = 0.0;
    out.p_value = 1.0;
    return out;
  }

  struct Key {
    u64 h;
    PairId p;
  };

  std::vector<Key> keys;
  keys.reserve(universe.size());
  for (const auto& p : universe) {
    keys.push_back(Key{detail::HashPair64(p), p});
  }

  std::sort(keys.begin(), keys.end(), [](const Key& a, const Key& b) {
    if (a.h != b.h) return a.h < b.h;
    if (a.p.r != b.p.r) return a.p.r < b.p.r;
    return a.p.s < b.p.s;
  });

  std::unordered_map<PairId, u32, PairIdHash> rank;
  rank.reserve(universe.size() * 2);
  for (u32 i = 0; i < static_cast<u32>(keys.size()); ++i) {
    rank.emplace(keys[i].p, i);
  }

  const double U = static_cast<double>(universe.size());
  std::vector<double> u;
  u.reserve(samples.size());

  for (usize i = 0; i < samples.size(); ++i) {
    const auto it = rank.find(samples[i]);
    const double r = (it == rank.end()) ? 0.0 : static_cast<double>(it->second);

    // Deterministic jitter in [0,1), varying per sample index.
    const u64 j = detail::SplitMix64(0xD1B54A32D192ED03ULL + static_cast<u64>(i));
    const u64 mant = (j >> 11);  // top 53 bits
    const double jitter = static_cast<double>(mant) * (1.0 / 9007199254740992.0);  // 2^53

    u.push_back((r + jitter) / U);
  }

  return KSOneSampleUniform01(std::move(u));
}

// --------------------------
// Convenience: evaluate PairId samples vs oracle universe
// --------------------------
struct PairUniformityResult {
  u64 universe_size = 0;  // |U|
  u64 sample_size = 0;    // t

  u64 missing_in_universe = 0;  // sample pairs not found in U (should be 0)

  u64 unique_in_sample = 0;
  double unique_fraction = 0.0;
  double duplicate_fraction = 0.0;

  double expected_per_bin = 0.0;  // t / |U|
  double l1 = 0.0;                // L1 distance between empirical and uniform
  double l_inf = 0.0;             // L-infinity distance
  double max_rel_error = 0.0;     // max |obs-exp| / exp

  ChiSquareResult chi2{};
};

// Evaluate uniformity of sample pairs over a known universe U (oracle pairs).
// Universe must represent the support of the target distribution (typically all join pairs).
inline PairUniformityResult EvaluatePairUniformity(Span<const PairId> universe,
                                                   Span<const PairId> samples) {
  PairUniformityResult out;
  out.universe_size = static_cast<u64>(universe.size());
  out.sample_size = static_cast<u64>(samples.size());
  if (universe.empty() || samples.empty()) {
    out.expected_per_bin = (universe.empty() ? 0.0 : static_cast<double>(samples.size()) /
                                             static_cast<double>(universe.size()));
    out.unique_in_sample = 0;
    out.unique_fraction = 0.0;
    out.duplicate_fraction = 0.0;
    out.l1 = 0.0;
    out.l_inf = 0.0;
    out.max_rel_error = 0.0;
    out.chi2 = ChiSquareResult{};
    return out;
  }

  // Build mapping from pair -> index.
  std::unordered_map<PairId, u32, PairIdHash> idx;
  idx.reserve(universe.size() * 2);

  for (usize i = 0; i < universe.size(); ++i) {
    idx.emplace(universe[i], static_cast<u32>(i));
  }

  // Count occurrences.
  std::vector<u64> counts(universe.size(), 0);
  for (const auto& p : samples) {
    auto it = idx.find(p);
    if (it == idx.end()) {
      out.missing_in_universe++;
      continue;
    }
    counts[it->second] += 1;
  }

  // Unique / duplicates.
  u64 uniq = 0;
  for (u64 c : counts) if (c > 0) ++uniq;
  out.unique_in_sample = uniq;
  out.unique_fraction = static_cast<double>(uniq) / static_cast<double>(samples.size());
  out.duplicate_fraction = 1.0 - out.unique_fraction;

  // Distances to uniform.
  const double U = static_cast<double>(universe.size());
  const double t = static_cast<double>(samples.size());
  const double p0 = 1.0 / U;
  out.expected_per_bin = t / U;

  double l1 = 0.0;
  double l_inf = 0.0;
  double max_rel = 0.0;

  for (u64 c : counts) {
    const double phat = static_cast<double>(c) / t;
    const double diff = std::fabs(phat - p0);
    l1 += diff;
    l_inf = std::max(l_inf, diff);

    if (out.expected_per_bin > 0.0) {
      const double rel = std::fabs(static_cast<double>(c) - out.expected_per_bin) / out.expected_per_bin;
      if (rel > max_rel) max_rel = rel;
    }
  }
  out.l1 = l1;
  out.l_inf = l_inf;
  out.max_rel_error = max_rel;

  out.chi2 = ChiSquareUniform(counts);
  return out;
}

}  // namespace quality
}  // namespace sampling
}  // namespace sjs
