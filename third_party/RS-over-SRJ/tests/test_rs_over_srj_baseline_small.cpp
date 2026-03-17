// third_party/RS-over-SRJ/tests/test_rs_over_srj_baseline_small.cpp
//
// Small correctness smoke test for the RS-over-SRJ KDTree (KDS) baseline.
//
// Validates:
//  - Count() equals the brute-force join size for a hand-made dataset
//  - Sample() returns only intersecting pairs

#include "isrjs_kds/isrjs_kds.hpp"

#include "core/config.h"
#include "core/rng.h"
#include "io/dataset.h"

#include <cstdint>
#include <iostream>
#include <string>

namespace {

template <int Dim, class T>
sjs::u64 BruteJoinSize(const sjs::Dataset<Dim, T>& ds) {
  sjs::u64 total = 0;
  for (sjs::usize i = 0; i < ds.R.boxes.size(); ++i) {
    for (sjs::usize j = 0; j < ds.S.boxes.size(); ++j) {
      if (ds.R.boxes[i].Intersects(ds.S.boxes[j])) ++total;
    }
  }
  return total;
}

}  // namespace

int main() {
  using T = sjs::Scalar;
  static constexpr int Dim = 2;

  using Baseline = sjs::baselines::rs_over_srj::RSOverSRJKDTreeSamplingBaseline<Dim, T>;

  // Hand-made dataset with boundary-touch cases.
  sjs::Dataset<Dim, T> ds;
  ds.name = "rs_over_srj_small";
  ds.half_open = true;

  auto P = [](T x, T y) { return sjs::Point<Dim, T>({x, y}); };
  auto B = [&](T x1, T y1, T x2, T y2) { return sjs::Box<Dim, T>(P(x1, y1), P(x2, y2)); };

  // R relation
  ds.R.Add(B(0.0, 0.0, 2.0, 2.0));  // r0
  ds.R.Add(B(5.0, 5.0, 6.0, 6.0));  // r1

  // S relation
  ds.S.Add(B(1.0, 1.0, 3.0, 3.0));  // intersects r0
  ds.S.Add(B(2.0, 2.0, 4.0, 4.0));  // touches r0 at corner only => NOT intersect
  ds.S.Add(B(5.5, 5.5, 7.0, 7.0));  // intersects r1

  // Use sequential IDs so PairId maps directly to box indices.
  ds.R.ForceSequentialIds();
  ds.S.ForceSequentialIds();

  std::string v_err;
  if (!ds.Validate(true, &v_err)) {
    std::cerr << "FAIL: dataset invalid: " << v_err << "\n";
    return 1;
  }

  const sjs::u64 brute = BruteJoinSize(ds);
  if (brute != 2) {
    std::cerr << "FAIL: brute join size expected 2, got " << brute << "\n";
    return 1;
  }

  sjs::Config cfg;
  cfg.dataset.dim = Dim;
  cfg.run.seed = 42;
  cfg.run.t = 1000;

  Baseline bl;
  std::string err;
  if (!bl.Build(ds, cfg, /*phases=*/nullptr, &err)) {
    std::cerr << "FAIL: Build failed: " << err << "\n";
    return 1;
  }

  sjs::baselines::CountResult cr;
  if (!bl.Count(cfg, /*rng=*/nullptr, &cr, /*phases=*/nullptr, &err)) {
    std::cerr << "FAIL: Count failed: " << err << "\n";
    return 1;
  }

  const sjs::u64 cnt = cr.RoundedU64();
  if (!cr.exact) {
    std::cerr << "FAIL: expected exact count, but cr.exact=false\n";
    return 1;
  }
  if (cnt != brute) {
    std::cerr << "FAIL: Count mismatch: got " << cnt << ", expected " << brute << "\n";
    return 1;
  }

  // Sampling correctness.
  sjs::Rng rng(cfg.run.seed);
  sjs::baselines::SampleSet ss;
  if (!bl.Sample(cfg, &rng, &ss, /*phases=*/nullptr, &err)) {
    std::cerr << "FAIL: Sample failed: " << err << "\n";
    return 1;
  }

  if (ss.pairs.size() != static_cast<sjs::usize>(cfg.run.t)) {
    std::cerr << "FAIL: expected " << cfg.run.t << " samples, got " << ss.pairs.size() << "\n";
    return 1;
  }

  for (const auto& p : ss.pairs) {
    const sjs::usize ri = static_cast<sjs::usize>(p.r);
    const sjs::usize si = static_cast<sjs::usize>(p.s);
    if (ri >= ds.R.boxes.size() || si >= ds.S.boxes.size()) {
      std::cerr << "FAIL: sampled id out of range: (" << p.r << "," << p.s << ")\n";
      return 1;
    }
    if (!ds.R.boxes[ri].Intersects(ds.S.boxes[si])) {
      std::cerr << "FAIL: sampled a non-intersecting pair: (" << p.r << "," << p.s << ")\n";
      return 1;
    }
  }

  std::cout << "PASS\n";
  return 0;
}
