// tests/test_baselines_smoke.cpp
//
// Smoke tests for all registered baselines (Dim=2) and their three variants:
//   - Sampling
//   - Enumerate+Sampling
//   - Adaptive
//
// The goal is NOT to benchmark. It is to ensure:
//   - factory/registry wiring works
//   - Build/Count/Sample/Enumerate do not crash
//   - Samples are valid join pairs (belong to the join output)
//   - Enumerators output exactly the oracle join pairs on a tiny dataset

#include "baselines/baseline_factory_2d.h"  // lives in src/

#include "sjs/core/config.h"
#include "sjs/core/rng.h"
#include "sjs/core/timer.h"
#include "sjs/geometry/predicates.h"
#include "sjs/io/dataset.h"
#include "sjs/join/join_oracle.h"

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
using sjs::baselines::BaselineSpec2D;
using sjs::baselines::CountResult;
using sjs::baselines::CreateBaseline2D;
using sjs::baselines::IBaseline;
using sjs::baselines::SampleSet;
using sjs::baselines::BaselineRegistry2D;

static bool PairLess(const PairId& a, const PairId& b) {
  if (a.r != b.r) return a.r < b.r;
  return a.s < b.s;
}

static Dataset<2, Scalar> MakeTinyJoinDataset() {
  using P2 = Point<2, Scalar>;
  using B2 = Box<2, Scalar>;

  Dataset<2, Scalar> ds;
  ds.name = "tiny_smoke";

  // R: three rectangles
  ds.R.name = "R";
  ds.R.Add(B2(P2{0.0, 0.0}, P2{1.0, 1.0}));    // r0
  ds.R.Add(B2(P2{0.5, 0.5}, P2{1.5, 1.5}));    // r1
  ds.R.Add(B2(P2{2.0, 2.0}, P2{3.0, 3.0}));    // r2

  // S: three rectangles
  ds.S.name = "S";
  ds.S.Add(B2(P2{0.2, 0.2}, P2{0.8, 0.8}));     // s0 intersects r0,r1
  ds.S.Add(B2(P2{1.0, 1.0}, P2{2.5, 2.5}));     // s1 intersects r1,r2 (NOT r0)
  ds.S.Add(B2(P2{3.0, 0.0}, P2{4.0, 1.0}));     // s2 disjoint

  ds.EnsureIds();
  return ds;
}

static void CheckAllSamplesAreJoinPairs(TestContext& t,
                                       const Dataset<2, Scalar>& ds,
                                       const std::unordered_set<PairId, PairIdHash>& oracle_set,
                                       const SampleSet& samples) {
  std::string err;
  CHECK(t, samples.Validate(&err));

  for (const auto& p : samples.pairs) {
    // 1) must be in oracle set
    CHECK(t, oracle_set.find(p) != oracle_set.end());

    // 2) and must actually intersect geometrically
    const auto& r = ds.R.boxes[static_cast<size_t>(p.r)];
    const auto& s = ds.S.boxes[static_cast<size_t>(p.s)];
    CHECK(t, sjs::IntersectsHalfOpen(r, s));
  }
}

static void CheckEnumeratorMatchesOracle(TestContext& t,
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

static void TestAllBaselines(TestContext& t) {
  const auto ds = MakeTinyJoinDataset();
  std::string err;
  CHECK(t, ds.Validate(true, &err));

  // Oracle join pairs.
  const auto oracle_pairs = sjs::join::CollectNaivePairs<2, Scalar>(ds.R, ds.S, /*cap=*/0);
  std::vector<PairId> oracle_sorted = oracle_pairs;
  std::sort(oracle_sorted.begin(), oracle_sorted.end(), PairLess);

  std::unordered_set<PairId, PairIdHash> oracle_set;
  oracle_set.reserve(oracle_pairs.size() * 2);
  for (const auto& p : oracle_pairs) oracle_set.insert(p);

  CHECK(t, !oracle_pairs.empty());

  // Base config used for all baselines.
  Config cfg;
  cfg.dataset.dim = 2;
  cfg.run.repeats = 1;
  cfg.run.seed = 42;
  cfg.run.t = 32;
  cfg.run.enum_cap = 0;      // no truncation
  cfg.run.j_star = 1000000;  // treat this dataset as "small" for adaptive

  const auto reg = BaselineRegistry2D();
  CHECK(t, reg.size() > 0);

  for (sjs::usize i = 0; i < reg.size(); ++i) {
    const BaselineSpec2D spec = reg[i];

    cfg.run.method = spec.method;
    cfg.run.variant = spec.variant;

    std::string err_local;
    auto baseline = CreateBaseline2D(spec.method, spec.variant, &err_local);
    if (!baseline) {
      std::cerr << "[FAIL] Could not create baseline for key=" << spec.key
                << "  err=" << err_local << "\n";
      ++t.fails;
      continue;
    }

    CHECK_EQ(t, baseline->method(), spec.method);
    CHECK_EQ(t, baseline->variant(), spec.variant);

    sjs::PhaseRecorder phases;
    phases.Clear();

    // Build
    err_local.clear();
    CHECK(t, baseline->Build(ds, cfg, &phases, &err_local));
    if (!err_local.empty()) {
      // Baselines should not set err on success.
      std::cerr << "[WARN] Build returned success but err was non-empty: " << err_local << "\n";
    }

    // Count
    CountResult cnt;
    err_local.clear();
    Rng rng_cnt(cfg.run.seed);
    CHECK(t, baseline->Count(cfg, &rng_cnt, &cnt, &phases, &err_local));

    // For this tiny dataset, all baselines are expected to produce a sensible count.
    CHECK(t, cnt.value >= 0.0L);
    CHECK(t, std::isfinite(static_cast<double>(cnt.value)));
    CHECK(t, cnt.RoundedU64() > 0u);
    if (cnt.exact) {
      CHECK_EQ(t, cnt.RoundedU64(), static_cast<sjs::u64>(oracle_pairs.size()));
    }

    // Sample
    SampleSet samples;
    err_local.clear();
    Rng rng_samp(cfg.run.seed + 1);
    CHECK(t, baseline->Sample(cfg, &rng_samp, &samples, &phases, &err_local));
    CHECK_EQ(t, samples.Size(), static_cast<sjs::usize>(cfg.run.t));

    CheckAllSamplesAreJoinPairs(t, ds, oracle_set, samples);

    // Enumerate
    err_local.clear();
    auto it = baseline->Enumerate(cfg, &phases, &err_local);
    CheckEnumeratorMatchesOracle(t, std::move(it), oracle_sorted);

    if (t.fails == 0) {
      std::cout << "[OK] baseline " << spec.key << "\n";
    } else {
      std::cout << "[INFO] baseline " << spec.key << " done (failures so far=" << t.fails << ")\n";
    }
  }
}

}  // namespace

int main() {
  TestContext t;
  TestAllBaselines(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_baselines_smoke\n";
    return 0;
  }
  std::cerr << "[FAILED] test_baselines_smoke: " << t.fails << " failure(s)\n";
  return 1;
}
