// tests/test_baselines_smoke_nd.cpp
//
// Smoke tests for baseline factory + runners in ND.
//
// Goal: ensure each baseline can be created and can run Build/Count/Sample (and
// EnumSampling/Adaptive where supported) on a small synthetic dataset without crashing.
// This test is intentionally lightweight and should finish fast.
//
// Assumptions (consistent with repository layout):
//  - sjs::Method includes: Ours, RangeTree, KDTree
//  - sjs::Variant includes: Sampling, EnumSampling, Adaptive
//  - Factory: sjs::baselines::CreateBaseline<Dim,T>(Method,Variant,err)
//  - Runners: RunSamplingOnce / RunEnumSamplingOnce / RunAdaptiveOnce

#include "sjs/baselines/baseline_api.h"
#include "sjs/baselines/baseline_factory.h"
#include "sjs/baselines/runners/adaptive_runner.h"
#include "sjs/baselines/runners/enum_sampling_runner.h"
#include "sjs/baselines/runners/sampling_runner.h"

#include "sjs/core/rng.h"
#include "sjs/io/dataset.h"
#include "sjs/join/join_oracle.h"

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

using sjs::Method;
using sjs::Rng;
using sjs::Scalar;
using sjs::Variant;
using sjs::u32;
using sjs::u64;
using sjs::usize;

template <int Dim>
sjs::Dataset<Dim, Scalar> MakeSmallDataset(u32 seed) {
  using BoxT = sjs::Box<Dim, Scalar>;

  Rng rng(seed);
  sjs::Dataset<Dim, Scalar> ds;
  ds.name = "smoke_nd_dim" + std::to_string(Dim);
  ds.R.name = "R";
  ds.S.name = "S";

  const int nr = 20;
  const int ns = 25;

  ds.R.boxes.reserve(nr);
  ds.S.boxes.reserve(ns);

  // Make at least one guaranteed intersecting pair to avoid an empty join.
  {
    BoxT b = BoxT::Empty();
    for (int i = 0; i < Dim; ++i) {
      const usize k = static_cast<usize>(i);
      b.lo.v[k] = 0.1;
      b.hi.v[k] = 0.9;
    }
    ds.R.boxes.push_back(b);
    ds.S.boxes.push_back(b);
  }

  auto rand_box = [&]() -> BoxT {
    BoxT b = BoxT::Empty();
    for (int i = 0; i < Dim; ++i) {
      const usize k = static_cast<usize>(i);
      const double lo = rng.UniformDouble(0.0, 0.6);
      const double w = rng.UniformDouble(0.25, 0.7);  // wide => likely intersect
      const double hi = std::min(1.0, lo + w);
      b.lo.v[k] = lo;
      b.hi.v[k] = std::max(hi, lo + 1e-6);
    }
    return b;
  };

  for (int i = 1; i < nr; ++i) ds.R.boxes.push_back(rand_box());
  for (int i = 1; i < ns; ++i) ds.S.boxes.push_back(rand_box());

  ds.R.ForceSequentialIds();
  ds.S.ForceSequentialIds();
  return ds;
}

template <int Dim>
bool RunOneBaseline(TestContext& t, Method method, Variant variant) {
  // Build dataset.
  auto ds = MakeSmallDataset<Dim>(1337 + static_cast<u32>(method) * 17u + static_cast<u32>(variant));

  // Exact oracle count for basic sanity.
  const u64 oracle = sjs::join::CountNaive<Dim, Scalar>(ds.R, ds.S, nullptr);
  CHECK(t, oracle > 0);

  // Prepare a config; only set the most essential fields.
  sjs::Config cfg{};
  cfg.dataset.dim = Dim;
  cfg.run.method = method;
  cfg.run.variant = variant;
  cfg.run.t = 200;
  cfg.run.repeats = 1;
  cfg.run.enum_cap = 500000;  // plenty for small join
  cfg.run.j_star = 500000;    // adaptive threshold (safe large)
  cfg.sys.threads = 1;

  std::string err;
  auto baseline = sjs::baselines::CreateBaseline<Dim, Scalar>(method, variant, &err);

  if (!baseline) {
    // KDTree may not support EnumSampling/Adaptive.
    const bool expected_missing =
        (method == Method::KDTree && (variant != Variant::Sampling));
    CHECK(t, expected_missing);
    return expected_missing;
  }

  sjs::baselines::RunReport report;
  bool ok = false;

  if (variant == Variant::Sampling) {
    ok = sjs::baselines::RunSamplingOnce<Dim, Scalar>(baseline.get(), ds, cfg, /*seed=*/7, &report, &err);
  } else if (variant == Variant::EnumSampling) {
    ok = sjs::baselines::RunEnumSamplingOnce<Dim, Scalar>(baseline.get(), ds, cfg, /*seed=*/7, &report, &err);
  } else if (variant == Variant::Adaptive) {
    ok = sjs::baselines::RunAdaptiveOnce<Dim, Scalar>(baseline.get(), ds, cfg, /*seed=*/7, &report, &err);
  } else {
    CHECK(t, false);
    return false;
  }

  CHECK(t, ok);
  CHECK(t, report.ok);
  if (!ok || !report.ok) {
    std::cerr << "baseline failed: method=" << static_cast<int>(method)
              << " variant=" << static_cast<int>(variant)
              << " err=" << err << "\n";
    return false;
  }

  // Basic report sanity.
  CHECK(t, report.count.value >= 0.0L);
  CHECK(t, std::isfinite(static_cast<double>(report.count.value)));

  // If join is non-empty, sampling should typically produce t samples.
  CHECK(t, report.samples.Validate(&err));
  CHECK(t, report.samples.Size() == static_cast<usize>(cfg.run.t));

  return true;
}

}  // namespace

int main() {
  TestContext t;

  // Smoke for a couple dims (keep runtime low).
  for (Method m : {Method::Ours, Method::RangeTree, Method::KDTree}) {
    for (Variant v : {Variant::Sampling, Variant::EnumSampling, Variant::Adaptive}) {
      (void)RunOneBaseline<2>(t, m, v);
      (void)RunOneBaseline<3>(t, m, v);
    }
  }

  if (t.fails == 0) {
    std::cout << "[OK] test_baselines_smoke_nd\n";
    return 0;
  }
  std::cerr << "[FAILED] test_baselines_smoke_nd: " << t.fails << " failure(s)\n";
  return 1;
}
