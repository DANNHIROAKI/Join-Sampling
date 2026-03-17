// tests/test_baselines_nd_smoke.cpp
//
// Small smoke tests for Dim=3/4/5 baselines.

#include "baselines/baseline_factory.h"

#include "core/config.h"
#include "core/rng.h"
#include "geometry/predicates.h"
#include "io/dataset.h"
#include "join/join_oracle.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <unordered_set>
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
using sjs::PairId;
using sjs::PairIdHash;
using sjs::Point;
using sjs::Rng;
using sjs::Scalar;
using sjs::Variant;
using sjs::Method;

static bool PairLess(const PairId& a, const PairId& b) {
  if (a.r != b.r) return a.r < b.r;
  return a.s < b.s;
}

template <int Dim>
Dataset<Dim, Scalar> MakeTinyDatasetND() {
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
  ds.name = std::string("tiny_nd_") + std::to_string(Dim);
  ds.R.name = "R";
  ds.S.name = "S";

  {
    auto [lo, hi] = fill(0.0, 0.60); ds.R.Add(B(lo, hi));
  }
  {
    auto [lo, hi] = fill(0.50, 1.10); ds.R.Add(B(lo, hi));
  }
  {
    auto [lo, hi] = fill(2.0, 2.5); ds.R.Add(B(lo, hi));
  }

  {
    auto [lo, hi] = fill(0.20, 0.80); ds.S.Add(B(lo, hi));
  }
  {
    auto [lo, hi] = fill(0.90, 1.40); ds.S.Add(B(lo, hi));
  }
  {
    auto [lo, hi] = fill(3.0, 3.5); ds.S.Add(B(lo, hi));
  }

  ds.EnsureIds();
  return ds;
}

template <int Dim>
void CheckAllSamplesAreJoinPairs(TestContext& t,
                                 const Dataset<Dim, Scalar>& ds,
                                 const std::unordered_set<PairId, PairIdHash>& oracle_set,
                                 const sjs::baselines::SampleSet& samples) {
  std::string err;
  CHECK(t, samples.Validate(&err));
  for (const auto& p : samples.pairs) {
    CHECK(t, oracle_set.find(p) != oracle_set.end());
    const auto& r = ds.R.boxes[static_cast<size_t>(p.r)];
    const auto& s = ds.S.boxes[static_cast<size_t>(p.s)];
    CHECK(t, sjs::IntersectsHalfOpen(r, s));
  }
}

template <int Dim>
void CheckEnumeratorMatchesOracle(TestContext& t,
                                  std::unique_ptr<sjs::baselines::IJoinEnumerator> it,
                                  const std::vector<PairId>& oracle_sorted) {
  CHECK(t, it != nullptr);
  std::vector<PairId> out;
  PairId p;
  while (it->Next(&p)) out.push_back(p);
  std::sort(out.begin(), out.end(), PairLess);
  CHECK_EQ(t, out.size(), oracle_sorted.size());
  for (size_t i = 0; i < out.size(); ++i) {
    CHECK(t, out[i] == oracle_sorted[i]);
  }
}

template <int Dim>
void TestDim(TestContext& t) {
  const auto ds = MakeTinyDatasetND<Dim>();
  std::string err;
  CHECK(t, ds.Validate(true, &err));

  const auto oracle_pairs = sjs::join::CollectNaivePairs<Dim, Scalar>(ds.R, ds.S, /*cap=*/0);
  std::vector<PairId> oracle_sorted = oracle_pairs;
  std::sort(oracle_sorted.begin(), oracle_sorted.end(), PairLess);

  std::unordered_set<PairId, PairIdHash> oracle_set;
  for (const auto& p : oracle_pairs) oracle_set.insert(p);
  CHECK(t, !oracle_pairs.empty());

  struct Combo { Method m; Variant v; } combos[] = {
      {Method::Ours, Variant::Sampling},
      {Method::Ours, Variant::EnumSampling},
      {Method::RSOverSRJ, Variant::Sampling},
  };

  Config cfg;
  cfg.dataset.dim = Dim;
  cfg.run.repeats = 1;
  cfg.run.seed = 17;
  cfg.run.t = 16;
  cfg.run.enum_cap = 0;

  for (const auto& c : combos) {
    cfg.run.method = c.m;
    cfg.run.variant = c.v;
    std::string err_local;
    auto baseline = sjs::baselines::CreateBaseline<Dim>(c.m, c.v, &err_local);
    CHECK(t, baseline != nullptr);
    if (!baseline) continue;

    sjs::PhaseRecorder phases;
    CHECK(t, baseline->Build(ds, cfg, &phases, &err_local));

    sjs::baselines::CountResult cnt;
    Rng rng_cnt(cfg.run.seed);
    CHECK(t, baseline->Count(cfg, &rng_cnt, &cnt, &phases, &err_local));
    CHECK(t, cnt.RoundedU64() == static_cast<sjs::u64>(oracle_pairs.size()));

    sjs::baselines::SampleSet samples;
    Rng rng_s(cfg.run.seed + 1);
    CHECK(t, baseline->Sample(cfg, &rng_s, &samples, &phases, &err_local));
    CHECK_EQ(t, samples.Size(), static_cast<sjs::usize>(cfg.run.t));
    CheckAllSamplesAreJoinPairs<Dim>(t, ds, oracle_set, samples);

    auto it = baseline->Enumerate(cfg, &phases, &err_local);
    CheckEnumeratorMatchesOracle<Dim>(t, std::move(it), oracle_sorted);
  }
}

}  // namespace

int main() {
  TestContext t;
  TestDim<3>(t);
  TestDim<4>(t);
  TestDim<5>(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_baselines_nd_smoke\n";
    return 0;
  }
  std::cerr << "[FAILED] test_baselines_nd_smoke: " << t.fails << " failure(s)\n";
  return 1;
}
