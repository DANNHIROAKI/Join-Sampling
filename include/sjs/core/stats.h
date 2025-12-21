#pragma once
// sjs/core/stats.h
//
// Simple statistics helpers for experiments.
// Provides:
//  - OnlineStats: Welford's algorithm for mean/variance, plus min/max.
//  - Summary: quantiles over a sample vector (requires sorting a copy).

#include "sjs/core/types.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace sjs {

class OnlineStats {
 public:
  OnlineStats() { Clear(); }

  void Clear() {
    n_ = 0;
    mean_ = 0.0L;
    m2_ = 0.0L;
    min_ = std::numeric_limits<long double>::infinity();
    max_ = -std::numeric_limits<long double>::infinity();
  }

  void Push(double x) {
    const long double v = static_cast<long double>(x);
    ++n_;
    const long double delta = v - mean_;
    mean_ += delta / static_cast<long double>(n_);
    const long double delta2 = v - mean_;
    m2_ += delta * delta2;
    if (v < min_) min_ = v;
    if (v > max_) max_ = v;
  }

  usize Count() const { return n_; }

  double Mean() const { return (n_ == 0) ? 0.0 : static_cast<double>(mean_); }

  // Unbiased sample variance by default (dividing by n-1).
  double Variance(bool unbiased = true) const {
    if (n_ < (unbiased ? 2 : 1)) return 0.0;
    const long double denom = unbiased ? static_cast<long double>(n_ - 1) : static_cast<long double>(n_);
    return static_cast<double>(m2_ / denom);
  }

  double Stddev(bool unbiased = true) const { return std::sqrt(Variance(unbiased)); }

  double Min() const { return (n_ == 0) ? 0.0 : static_cast<double>(min_); }
  double Max() const { return (n_ == 0) ? 0.0 : static_cast<double>(max_); }

 private:
  usize n_{0};
  long double mean_{0.0L};
  long double m2_{0.0L};
  long double min_{0.0L};
  long double max_{0.0L};
};

// Quantile on a *sorted* vector (ascending). q in [0,1].
// Uses linear interpolation between adjacent points (like numpy default).
inline double QuantileSorted(const std::vector<double>& sorted, double q) {
  if (sorted.empty()) return 0.0;
  if (q <= 0.0) return sorted.front();
  if (q >= 1.0) return sorted.back();

  const double pos = q * static_cast<double>(sorted.size() - 1);
  const usize idx = static_cast<usize>(pos);
  const double frac = pos - static_cast<double>(idx);

  const double a = sorted[idx];
  const double b = sorted[std::min(idx + 1, sorted.size() - 1)];
  return a + (b - a) * frac;
}

struct Summary {
  usize n{0};

  double mean{0.0};
  double stdev{0.0};
  double sem{0.0};  // standard error of the mean

  double min{0.0};
  double max{0.0};

  double median{0.0};
  double p90{0.0};
  double p95{0.0};
  double p99{0.0};

  std::string ToJson() const {
    std::ostringstream oss;
    oss << "{"
        << "\"n\":" << n << ","
        << "\"mean\":" << mean << ","
        << "\"stdev\":" << stdev << ","
        << "\"sem\":" << sem << ","
        << "\"min\":" << min << ","
        << "\"max\":" << max << ","
        << "\"median\":" << median << ","
        << "\"p90\":" << p90 << ","
        << "\"p95\":" << p95 << ","
        << "\"p99\":" << p99
        << "}";
    return oss.str();
  }
};

// Compute Summary from raw samples; makes a sorted copy for quantiles.
inline Summary Summarize(const std::vector<double>& samples) {
  Summary s;
  if (samples.empty()) return s;

  // Be robust to failed/placeholder measurements: ignore NaN/Inf values.
  std::vector<double> finite;
  finite.reserve(samples.size());
  for (double x : samples) {
    if (std::isfinite(x)) finite.push_back(x);
  }
  s.n = finite.size();
  if (finite.empty()) return s;

  OnlineStats os;
  for (double x : finite) os.Push(x);

  s.mean = os.Mean();
  s.stdev = os.Stddev(/*unbiased=*/true);
  s.sem = (s.n > 0) ? (s.stdev / std::sqrt(static_cast<double>(s.n))) : 0.0;
  s.min = os.Min();
  s.max = os.Max();

  std::vector<double> sorted = finite;
  std::sort(sorted.begin(), sorted.end());
  s.median = QuantileSorted(sorted, 0.5);
  s.p90 = QuantileSorted(sorted, 0.90);
  s.p95 = QuantileSorted(sorted, 0.95);
  s.p99 = QuantileSorted(sorted, 0.99);
  return s;
}
}  // namespace sjs
