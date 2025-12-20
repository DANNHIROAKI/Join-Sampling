// src/sampling/quality_checks.cpp
//
// Small-scale sampling quality diagnostics.
//
// The heavy lifting (chi-square, KS, autocorrelation, pair-uniformity) lives in
// the header-only implementation:
//   - include/sjs/sampling/sample_quality.h
//
// This file provides convenient wrappers and JSON-lite serialization helpers
// so apps can record diagnostics without duplicating formatting logic.

#include "sjs/sampling/sample_quality.h"

#include "sjs/core/types.h"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace sjs {
namespace sampling {

namespace {

inline void JsonNumOrNull(std::ostringstream& oss, double x) {
  if (std::isfinite(x)) {
    oss << x;
  } else {
    oss << "null";
  }
}

inline void JsonU64(std::ostringstream& oss, u64 x) {
  oss << static_cast<unsigned long long>(x);
}

}  // namespace

std::string ChiSquareToJsonLite(const quality::ChiSquareResult& r) {
  std::ostringstream oss;
  oss << "{";
  oss << "\"statistic\":"; JsonNumOrNull(oss, r.statistic);
  oss << ",\"df\":"; JsonU64(oss, r.df);
  oss << ",\"p_value\":"; JsonNumOrNull(oss, r.p_value);
  oss << "}";
  return oss.str();
}

std::string KSTestToJsonLite(const quality::KSTestResult& r) {
  std::ostringstream oss;
  oss << "{";
  oss << "\"D\":"; JsonNumOrNull(oss, r.D);
  oss << ",\"p_value\":"; JsonNumOrNull(oss, r.p_value);
  oss << "}";
  return oss.str();
}

std::string PairUniformityToJsonLite(const quality::PairUniformityResult& r) {
  std::ostringstream oss;
  oss << "{";
  oss << "\"universe_size\":"; JsonU64(oss, r.universe_size);
  oss << ",\"sample_size\":"; JsonU64(oss, r.sample_size);
  oss << ",\"missing_in_universe\":"; JsonU64(oss, r.missing_in_universe);
  oss << ",\"unique_in_sample\":"; JsonU64(oss, r.unique_in_sample);
  oss << ",\"unique_fraction\":"; JsonNumOrNull(oss, r.unique_fraction);
  oss << ",\"duplicate_fraction\":"; JsonNumOrNull(oss, r.duplicate_fraction);
  oss << ",\"expected_per_bin\":"; JsonNumOrNull(oss, r.expected_per_bin);
  oss << ",\"l1\":"; JsonNumOrNull(oss, r.l1);
  oss << ",\"l_inf\":"; JsonNumOrNull(oss, r.l_inf);
  oss << ",\"max_rel_error\":"; JsonNumOrNull(oss, r.max_rel_error);
  oss << ",\"chi2\":" << ChiSquareToJsonLite(r.chi2);
  oss << "}";
  return oss.str();
}

// Evaluate a bundle of common diagnostics for PairId samples.
//
// Returns a JSON-lite object:
// {
//   "pair_uniformity": {...},
//   "autocorr_hash_lag": <number|null>,
//   "ks_hash_uniform01": {...}
// }
std::string EvaluatePairSampleQualityJsonLite(Span<const PairId> universe,
                                              Span<const PairId> samples,
                                              int autocorr_lag) {
  const quality::PairUniformityResult uni = quality::EvaluatePairUniformity(universe, samples);

  const double ac = quality::AutocorrelationHashedPairs(samples, autocorr_lag);

  // KS sanity check: map samples into [0,1) using a universe-rank + jitter scheme
  // to avoid systematic rejection on small discrete universes (see sample_quality.h).
  const quality::KSTestResult ks = quality::KSPairsHashUniform01RankJitter(universe, samples);

  std::ostringstream oss;
  oss << "{";
  oss << "\"pair_uniformity\":" << PairUniformityToJsonLite(uni);
  oss << ",\"autocorr_hash_lag\":"; JsonNumOrNull(oss, ac);
  oss << ",\"ks_hash_uniform01\":" << KSTestToJsonLite(ks);
  oss << "}";
  return oss.str();
}

}  // namespace sampling
}  // namespace sjs
