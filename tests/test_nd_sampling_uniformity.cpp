// tests/test_nd_sampling_uniformity.cpp
//
// Small deterministic regression tests for ND sampling quality. The goal is
// not to prove statistical correctness, but to catch obvious regressions in
// uniformity / i.i.d. behavior on a tiny oracle-known dataset.

#include "baselines/baseline_factory.h"

#include "core/config.h"
#include "core/rng.h"
#include "geometry/box.h"
#include "io/dataset.h"
#include "join/join_oracle.h"
#include "sampling/sample_quality.h"

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
};

#define CHECK(ctx, expr) (ctx).Check((expr), #expr, __FILE__, __LINE__)
#define CHECK_EQ(ctx, a, b) (ctx).CheckEq((a), (b), #a, #b, __FILE__, __LINE__)

using sjs::Box;
using sjs::Config;
using sjs::Dataset;
using sjs::Method;
using sjs::PairId;
using sjs::Point;
using sjs::Rng;
using sjs::Scalar;
using sjs::Variant;
using sjs::u64;

static bool PairLess(const PairId& a, const PairId& b) {
  if (a.r != b.r) return a.r < b.r;
  return a.s < b.s;
}

template <int Dim>
Dataset<Dim, Scalar> MakeUniformityDataset() {
  using P = Point<Dim, Scalar>;
  using B = Box<Dim, Scalar>;

  auto fill = [](double lo, double hi) {
    P p_lo{};
    P p_hi{};
    for (int i = 0; i < Dim; ++i) {
      p_lo[i] = lo;
      p_hi[i] = hi;
    }
    return std::pair<P, P>{p_lo, p_hi};
  };

  Dataset<Dim, Scalar> ds;
  ds.name = std::string("uniformity_nd_") + std::to_string(Dim);
  ds.R.name = "R";
  ds.S.name = "S";

  {
    auto [lo, hi] = fill(0.00, 0.55); ds.R.Add(B(lo, hi));
  }
  {
    auto [lo, hi] = fill(0.40, 0.95); ds.R.Add(B(lo, hi));
  }
  {
    auto [lo, hi] = fill(0.80, 1.35); ds.R.Add(B(lo, hi));
  }
  {
    auto [lo, hi] = fill(2.00, 2.40); ds.R.Add(B(lo, hi));
  }

  {
    auto [lo, hi] = fill(0.10, 0.65); ds.S.Add(B(lo, hi));
  }
  {
    auto [lo, hi] = fill(0.50, 1.05); ds.S.Add(B(lo, hi));
  }
  {
    auto [lo, hi] = fill(0.90, 1.45); ds.S.Add(B(lo, hi));
  }
  {
    auto [lo, hi] = fill(3.00, 3.50); ds.S.Add(B(lo, hi));
  }

  ds.EnsureIds();
  return ds;
}

template <int Dim>
void RunSamplingQualityCheck(TestContext& t, Method method, u64 seed) {
  const auto ds = MakeUniformityDataset<Dim>();
  std::string err;
  CHECK(t, ds.Validate(true, &err));

  auto universe = sjs::join::CollectNaivePairs<Dim, Scalar>(ds.R, ds.S, /*cap=*/0, nullptr);
  std::sort(universe.begin(), universe.end(), PairLess);
  CHECK_EQ(t, universe.size(), static_cast<size_t>(7));

  Config cfg;
  cfg.dataset.dim = Dim;
  cfg.run.method = method;
  cfg.run.variant = Variant::Sampling;
  cfg.run.repeats = 1;
  cfg.run.seed = seed;
  cfg.run.t = 14000;
  cfg.run.enum_cap = 0;

  std::string build_err;
  auto baseline = sjs::baselines::CreateBaseline<Dim>(method, Variant::Sampling, &build_err);
  CHECK(t, baseline != nullptr);
  if (!baseline) return;

  sjs::PhaseRecorder phases;
  CHECK(t, baseline->Build(ds, cfg, &phases, &build_err));

  sjs::baselines::CountResult cnt;
  Rng rng_cnt(seed);
  CHECK(t, baseline->Count(cfg, &rng_cnt, &cnt, &phases, &build_err));
  CHECK_EQ(t, cnt.RoundedU64(), static_cast<sjs::u64>(universe.size()));

  sjs::baselines::SampleSet samples;
  Rng rng(seed + 1);
  CHECK(t, baseline->Sample(cfg, &rng, &samples, &phases, &build_err));
  CHECK_EQ(t, samples.Size(), static_cast<sjs::usize>(cfg.run.t));

  const auto uni = sjs::sampling::quality::EvaluatePairUniformity(
      sjs::Span<const PairId>(universe.data(), universe.size()),
      sjs::Span<const PairId>(samples.pairs.data(), samples.pairs.size()));
  const auto ks = sjs::sampling::quality::KSPairsHashUniform01RankJitter(
      sjs::Span<const PairId>(universe.data(), universe.size()),
      sjs::Span<const PairId>(samples.pairs.data(), samples.pairs.size()));
  const double ac1 = sjs::sampling::quality::AutocorrelationHashedPairs(
      sjs::Span<const PairId>(samples.pairs.data(), samples.pairs.size()), 1);

  CHECK_EQ(t, uni.missing_in_universe, 0u);
  CHECK(t, uni.chi2.statistic < 30.0);
  CHECK(t, uni.l_inf < 0.03);
  CHECK(t, ks.D < 0.05);
  CHECK(t, std::fabs(ac1) < 0.08);
}

}  // namespace

int main() {
  TestContext t;

  RunSamplingQualityCheck<3>(t, Method::Ours, 17);
  RunSamplingQualityCheck<3>(t, Method::RSOverSRJ, 29);
  RunSamplingQualityCheck<5>(t, Method::Ours, 41);
  RunSamplingQualityCheck<5>(t, Method::RSOverSRJ, 53);

  if (t.fails == 0) {
    std::cout << "[OK] test_nd_sampling_uniformity\n";
    return 0;
  }
  std::cerr << "[FAILED] test_nd_sampling_uniformity: " << t.fails << " failure(s)\n";
  return 1;
}
