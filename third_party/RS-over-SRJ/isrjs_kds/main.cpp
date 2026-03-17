// isrjs_kds/main.cpp (modified)
//
// NOTE:
// This file is kept ONLY as a tiny sanity-check demo.
// When vendoring RS-over-SRJ into the main repository as third_party code, you
// should generally NOT compile this translation unit into the main benchmark
// binaries (the project already has apps/sjs_run).
//
// To build this demo explicitly, compile with:
//   -DRS_OVER_SRJ_DEMO_MAIN

#ifdef RS_OVER_SRJ_DEMO_MAIN

#include "isrjs_kds.hpp"

#include "core/config.h"
#include "core/rng.h"
#include "io/dataset.h"

#include <iostream>
#include <string>

int main() {
  using Baseline = sjs::baselines::rs_over_srj::RSOverSRJKDTreeSamplingBaseline<2, sjs::Scalar>;

  // Tiny hand-made dataset.
  sjs::Dataset<2, sjs::Scalar> ds;
  ds.name = "demo";
  ds.half_open = true;

  // R: two rectangles
  ds.R.boxes.push_back(sjs::Box<2, sjs::Scalar>{sjs::Point<2, sjs::Scalar>({0.0, 0.0}),
                                               sjs::Point<2, sjs::Scalar>({2.0, 2.0})});
  ds.R.boxes.push_back(sjs::Box<2, sjs::Scalar>{sjs::Point<2, sjs::Scalar>({5.0, 5.0}),
                                               sjs::Point<2, sjs::Scalar>({6.0, 6.0})});

  // S: three rectangles
  ds.S.boxes.push_back(sjs::Box<2, sjs::Scalar>{sjs::Point<2, sjs::Scalar>({1.0, 1.0}),
                                               sjs::Point<2, sjs::Scalar>({3.0, 3.0})});  // intersects R0
  ds.S.boxes.push_back(sjs::Box<2, sjs::Scalar>{sjs::Point<2, sjs::Scalar>({2.0, 2.0}),
                                               sjs::Point<2, sjs::Scalar>({4.0, 4.0})});  // touches corner => no strict overlap
  ds.S.boxes.push_back(sjs::Box<2, sjs::Scalar>{sjs::Point<2, sjs::Scalar>({5.5, 5.5}),
                                               sjs::Point<2, sjs::Scalar>({7.0, 7.0})});  // intersects R1

  std::string v_err;
  if (!ds.Validate(true, &v_err)) {
    std::cerr << "Dataset invalid: " << v_err << "\n";
    return 1;
  }

  sjs::Config cfg;
  cfg.dataset.dim = 2;
  cfg.run.t = 10;
  cfg.run.seed = 1;

  Baseline bl;

  std::string err;
  if (!bl.Build(ds, cfg, /*phases=*/nullptr, &err)) {
    std::cerr << "Build failed: " << err << "\n";
    return 1;
  }

  sjs::baselines::CountResult cr;
  if (!bl.Count(cfg, /*rng=*/nullptr, &cr, /*phases=*/nullptr, &err)) {
    std::cerr << "Count failed: " << err << "\n";
    return 1;
  }
  std::cout << "|J| = " << cr.RoundedU64() << " (exact=" << (cr.exact ? "true" : "false") << ")\n";

  sjs::Rng rng(cfg.run.seed);
  sjs::baselines::SampleSet ss;
  if (!bl.Sample(cfg, &rng, &ss, /*phases=*/nullptr, &err)) {
    std::cerr << "Sample failed: " << err << "\n";
    return 1;
  }

  std::cout << "Samples (r_id, s_id):\n";
  for (const auto& p : ss.pairs) {
    std::cout << "(" << p.r << "," << p.s << ")\n";
  }

  return 0;
}

#endif  // RS_OVER_SRJ_DEMO_MAIN