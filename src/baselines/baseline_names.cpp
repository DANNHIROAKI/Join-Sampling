// src/baselines/baseline_names.cpp
//
// Baseline registry + help strings (HighDims).
//
// Centralized here so:
//  - CLI / error messages can print a stable list
//  - Factory can reference BaselineHelp() for diagnostics

#include "sjs/baselines/baseline_names.h"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace sjs {
namespace baselines {

Span<const BaselineSpec> BaselineRegistry() noexcept {
  // IMPORTANT:
  // Keep this list consistent with CreateBaseline<Dim>() in
  // src/baselines/baseline_factory_dispatch.cpp.

  static const BaselineSpec kReg[] = {
      // Ours (HighDims)
      {Method::Ours, Variant::Sampling, "ours/sampling",
       "Our method (Framework II): 2-pass sweep + HighDims event-block indices, exact i.i.d. uniform"},
      {Method::Ours, Variant::EnumSampling, "ours/enum_sampling",
       "Our method (Framework I): materialize full join then uniform index sampling"},
      {Method::Ours, Variant::Adaptive, "ours/adaptive",
       "Our method (Framework III): budgeted full-cache + (optional) adaptive sample-cache"},

      // Range-tree baseline
      {Method::RangeTree, Variant::Sampling, "range_tree/sampling",
       "Baseline (Framework II): sweep + event-block primitives via range-tree-like structure (correctness-first build may scan)"},
      {Method::RangeTree, Variant::EnumSampling, "range_tree/enum_sampling",
       "Range-tree baseline (Framework I): materialize full join then uniform index sampling"},
      {Method::RangeTree, Variant::Adaptive, "range_tree/adaptive",
       "Range-tree baseline (Framework III): budgeted full-cache + adaptive sampling"},

      // KD-tree baseline (sampling-only)
      {Method::KDTree, Variant::Sampling, "kd_tree/sampling",
       "KD-tree baseline: exact COUNT + exact SAMPLE per event-block (sampling only)"},
  };

  return Span<const BaselineSpec>(kReg, sizeof(kReg) / sizeof(kReg[0]));
}

bool IsBaselineSupported(Method method, Variant variant) noexcept {
  const auto reg = BaselineRegistry();
  for (usize i = 0; i < reg.size(); ++i) {
    if (reg[i].method == method && reg[i].variant == variant) return true;
  }
  return false;
}

std::string BaselineHelp() {
  std::ostringstream oss;

  oss << "  Methods (canonical):\n";
  {
    std::vector<Method> seen;
    const auto reg = BaselineRegistry();
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

  oss << "\n  Supported method/variant combinations:\n";
  const auto reg = BaselineRegistry();
  for (usize i = 0; i < reg.size(); ++i) {
    const auto& r = reg[i];
    oss << "    - " << r.key << "  (" << ToString(r.method) << ", " << ToString(r.variant) << ")";
    if (!r.desc.empty()) oss << ": " << r.desc;
    oss << "\n";
  }

  oss << "\n  Example:\n";
  oss << "    --method=ours --variant=sampling --t=100000 --seed=1\n";

  return oss.str();
}

}  // namespace baselines
}  // namespace sjs
