// tests/test_sampling_quality.cpp
//
// Sampling/statistics sanity tests:
//  - Alias table behavior
//  - Chi-square uniform test helper
//  - KS one-sample uniform test helper
//  - Pair-uniformity evaluator vs a known universe

#include "sjs/core/rng.h"
#include "sjs/core/types.h"
#include "sjs/sampling/alias_table.h"
#include "sjs/sampling/sample_quality.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct TestContext {
  int fails = 0;

  void Check(bool ok, const char* expr, const char* file, int line) {
    if (ok) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line << "  CHECK(" << expr << ")\n";
  }

  template <class A, class B>
  void CheckEq(const A& a, const B& b, const char* ea, const char* eb,
               const char* file, int line) {
    if (a == b) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line << "  CHECK_EQ(" << ea << ", " << eb
              << ")  got " << a << " vs " << b << "\n";
  }

  void CheckNear(double a, double b, double rel_eps, const char* ea, const char* eb,
                 const char* file, int line) {
    const double diff = std::fabs(a - b);
    const double scale = std::max({1.0, std::fabs(a), std::fabs(b)});
    if (diff <= rel_eps * scale) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line
              << "  CHECK_NEAR(" << ea << ", " << eb << ")  got " << a << " vs " << b
              << "  diff=" << diff << "\n";
  }
};

#define CHECK(ctx, expr) (ctx).Check((expr), #expr, __FILE__, __LINE__)
#define CHECK_EQ(ctx, a, b) (ctx).CheckEq((a), (b), #a, #b, __FILE__, __LINE__)
#define CHECK_NEAR(ctx, a, b, eps) (ctx).CheckNear((a), (b), (eps), #a, #b, __FILE__, __LINE__)

using sjs::PairId;
using sjs::PairIdHash;
using sjs::Rng;

static void TestAliasTableBasic(TestContext& t) {
  // Two outcomes with weights 1:3.
  const std::vector<double> w = {1.0, 3.0};
  sjs::sampling::AliasTable tab;
  CHECK(t, tab.Build(w));

  // Empirical frequency should be close to 1/4 and 3/4.
  Rng rng(123);
  const int N = 8000;
  int c0 = 0, c1 = 0;
  for (int i = 0; i < N; ++i) {
    const sjs::u32 idx = static_cast<sjs::u32>(tab.Sample(&rng));
    CHECK(t, idx < 2);
    if (idx == 0) ++c0;
    else ++c1;
  }

  const double p0 = static_cast<double>(c0) / static_cast<double>(N);
  const double p1 = static_cast<double>(c1) / static_cast<double>(N);

  CHECK(t, p0 > 0.15 && p0 < 0.35);
  CHECK(t, p1 > 0.65 && p1 < 0.85);
}

static void TestChiSquareUniform(TestContext& t) {
  // Perfectly uniform counts => chi2=0, p=1.
  const std::vector<sjs::u64> counts = {10, 10, 10, 10};
  const auto r = sjs::sampling::quality::ChiSquareUniform(counts);

  CHECK_NEAR(t, r.statistic, 0.0, 1e-12);
  CHECK_NEAR(t, r.p_value, 1.0, 1e-12);
  CHECK_EQ(t, r.df, 3u);
}

static void TestKSUniform01(TestContext& t) {
  // Draw some pseudo-uniform samples; KS statistic should be small-ish.
  Rng rng(7);
  std::vector<double> x;
  x.reserve(200);
  for (int i = 0; i < 200; ++i) x.push_back(rng.NextDouble());

  const auto r = sjs::sampling::quality::KSOneSampleUniform01(x);
  CHECK(t, std::isfinite(r.D));
  CHECK(t, r.D >= 0.0);
  CHECK(t, r.D < 0.25);  // loose bound
  CHECK(t, r.p_value >= 0.0 && r.p_value <= 1.0);
}

static void TestPairUniformityEvaluator(TestContext& t) {
  // Universe of 5 join pairs.
  std::vector<PairId> U;
  for (int i = 0; i < 5; ++i) U.emplace_back(0, i);

  // Sample: repeat each element 10 times => exactly uniform.
  std::vector<PairId> S;
  for (int rep = 0; rep < 10; ++rep) {
    for (int i = 0; i < 5; ++i) S.emplace_back(0, i);
  }

  const auto res = sjs::sampling::quality::EvaluatePairUniformity(
      sjs::Span<const PairId>(U.data(), U.size()),
      sjs::Span<const PairId>(S.data(), S.size()));

  CHECK_EQ(t, res.universe_size, 5u);
  CHECK_EQ(t, res.sample_size, 50u);
  CHECK_EQ(t, res.missing_in_universe, 0u);
  CHECK_EQ(t, res.unique_in_sample, 5u);
  CHECK_NEAR(t, res.l1, 0.0, 1e-12);
  CHECK_NEAR(t, res.chi2.statistic, 0.0, 1e-12);
  CHECK_NEAR(t, res.chi2.p_value, 1.0, 1e-12);
}

static void TestAutocorrelation(TestContext& t) {
  // A monotonically increasing series should have a high lag-1 autocorrelation.
  std::vector<double> x;
  x.reserve(100);
  for (int i = 0; i < 100; ++i) x.push_back(static_cast<double>(i));

  const double r1 = sjs::sampling::quality::Autocorrelation(
      sjs::Span<const double>(x.data(), x.size()), 1);

  CHECK(t, std::isfinite(r1));
  CHECK(t, r1 > 0.9);
}

}  // namespace

int main() {
  TestContext t;

  TestAliasTableBasic(t);
  TestChiSquareUniform(t);
  TestKSUniform01(t);
  TestPairUniformityEvaluator(t);
  TestAutocorrelation(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_sampling_quality\n";
    return 0;
  }
  std::cerr << "[FAILED] test_sampling_quality: " << t.fails << " failure(s)\n";
  return 1;
}
