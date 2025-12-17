// src/baselines/baseline_names.cpp
//
// Baseline registry + help strings.
//
// Keeping this in a .cpp (not a header) reduces compile-time and keeps the
// registry centralized for CLI/help/error messaging.

#include "baselines/baseline_factory_2d.h"

#include <algorithm>
// #include <ostringstream>
#include <sstream>
#include <string>
#include <vector>

namespace sjs {
namespace baselines {

Span<const BaselineSpec2D> BaselineRegistry2D() noexcept {
  // IMPORTANT: Keep this list consistent with CreateBaseline2D(...) in
  // baseline_factory_2d.cpp.
  //
  // Canonical key convention: "<method>/<variant>".
  static const BaselineSpec2D kReg[] = {
      // Ours
      {Method::Ours, Variant::Sampling, "ours/sampling",
       "Our method (Sampling): 2-pass sweep, exact i.i.d. uniform join samples"},
      {Method::Ours, Variant::EnumSampling, "ours/enum_sampling",
       "Our method (Enum+Sampling): enumerate join stream then rank-sample"},
      {Method::Ours, Variant::Adaptive, "ours/adaptive",
       "Our method (Adaptive): enumerate if small else fallback to sampling"},

      // AABB-tree baseline
      {Method::AABB, Variant::Sampling, "aabb/sampling", "Plane sweep + dynamic AABB-tree"},
      {Method::AABB, Variant::EnumSampling, "aabb/enum_sampling", "AABB baseline + enumerate then sample"},
      {Method::AABB, Variant::Adaptive, "aabb/adaptive", "AABB baseline + adaptive small/large join"},

      // Interval tree baseline
      {Method::IntervalTree, Variant::Sampling, "interval_tree/sampling", "Plane sweep + interval tree (y-dim)"},
      {Method::IntervalTree, Variant::EnumSampling, "interval_tree/enum_sampling",
       "Interval-tree baseline + enumerate then sample"},
      {Method::IntervalTree, Variant::Adaptive, "interval_tree/adaptive",
       "Interval-tree baseline + adaptive small/large join"},

      // KD-tree baseline
      {Method::KDTree, Variant::Sampling, "kd_tree/sampling", "Plane sweep + KD-tree (embedding)"},
      {Method::KDTree, Variant::EnumSampling, "kd_tree/enum_sampling", "KD-tree baseline + enumerate then sample"},
      {Method::KDTree, Variant::Adaptive, "kd_tree/adaptive", "KD-tree baseline + adaptive small/large join"},

      // R-tree baseline
      {Method::RTree, Variant::Sampling, "r_tree/sampling", "Plane sweep + dynamic R-tree"},
      {Method::RTree, Variant::EnumSampling, "r_tree/enum_sampling", "R-tree baseline + enumerate then sample"},
      {Method::RTree, Variant::Adaptive, "r_tree/adaptive", "R-tree baseline + adaptive small/large join"},

      // Range tree baseline
      {Method::RangeTree, Variant::Sampling, "range_tree/sampling", "Plane sweep + range tree"},
      {Method::RangeTree, Variant::EnumSampling, "range_tree/enum_sampling", "Range-tree baseline + enumerate then sample"},
      {Method::RangeTree, Variant::Adaptive, "range_tree/adaptive", "Range-tree baseline + adaptive small/large join"},

      // PBSM baseline
      {Method::PBSM, Variant::Sampling, "pbsm/sampling", "PBSM-style partition-based spatial merge"},
      {Method::PBSM, Variant::EnumSampling, "pbsm/enum_sampling", "PBSM baseline + enumerate then sample"},
      {Method::PBSM, Variant::Adaptive, "pbsm/adaptive", "PBSM baseline + adaptive small/large join"},

      // TLSOP baseline
      {Method::TLSOP, Variant::Sampling, "tlsop/sampling", "Two-layer SOP sampling baseline"},
      {Method::TLSOP, Variant::EnumSampling, "tlsop/enum_sampling", "TLSOP baseline + enumerate then sample"},
      {Method::TLSOP, Variant::Adaptive, "tlsop/adaptive", "TLSOP baseline + adaptive small/large join"},

      // SIRS baseline
      {Method::SIRS, Variant::Sampling, "sirs/sampling", "SIRS sampling baseline"},
      {Method::SIRS, Variant::EnumSampling, "sirs/enum_sampling", "SIRS baseline + enumerate then sample"},
      {Method::SIRS, Variant::Adaptive, "sirs/adaptive", "SIRS baseline + adaptive small/large join"},

      // Rejection baseline (KDS-style)
      {Method::Rejection, Variant::Sampling, "rejection/sampling", "Rejection sampling baseline"},
      {Method::Rejection, Variant::EnumSampling, "rejection/enum_sampling", "Rejection baseline + enumerate then sample"},
      {Method::Rejection, Variant::Adaptive, "rejection/adaptive", "Rejection baseline + adaptive small/large join"},

      // Tsunami baseline (NOTE: requires Method::Tsunami to exist in core/types.h)
      {Method::Tsunami, Variant::Sampling, "tsunami/sampling", "Tsunami-style sampling baseline"},
      {Method::Tsunami, Variant::EnumSampling, "tsunami/enum_sampling", "Tsunami baseline + enumerate then sample"},
      {Method::Tsunami, Variant::Adaptive, "tsunami/adaptive", "Tsunami baseline + adaptive small/large join"},
  };

  return Span<const BaselineSpec2D>(kReg, sizeof(kReg) / sizeof(kReg[0]));
}

bool IsBaselineSupported2D(Method method, Variant variant) noexcept {
  const auto reg = BaselineRegistry2D();
  for (usize i = 0; i < reg.size(); ++i) {
    if (reg[i].method == method && reg[i].variant == variant) return true;
  }
  return false;
}

std::string BaselineHelp2D() {
  std::ostringstream oss;

  oss << "  Methods (canonical):\n";
  // Gather unique method names in registry order.
  {
    std::vector<Method> seen;
    const auto reg = BaselineRegistry2D();
    seen.reserve(reg.size());
    for (usize i = 0; i < reg.size(); ++i) {
      const Method m = reg[i].method;
      if (std::find(seen.begin(), seen.end(), m) == seen.end()) seen.push_back(m);
    }
    for (Method m : seen) {
      oss << "    - " << ToString(m) << "\n";
    }
  }

  oss << "\n  Variants (canonical):\n";
  oss << "    - sampling\n";
  oss << "    - enum_sampling\n";
  oss << "    - adaptive\n";

  oss << "\n  Supported (Dim=2) method/variant combinations:\n";
  const auto reg = BaselineRegistry2D();
  for (usize i = 0; i < reg.size(); ++i) {
    const auto& r = reg[i];
    oss << "    - " << r.key << "  (" << ToString(r.method) << ", " << ToString(r.variant) << ")";
    if (!r.desc.empty()) {
      oss << ": " << r.desc;
    }
    oss << "\n";
  }

  oss << "\n  Example:\n";
  oss << "    --method=ours --variant=sampling --t=100000 --seed=1\n";

  return oss.str();
}

}  // namespace baselines
}  // namespace sjs
