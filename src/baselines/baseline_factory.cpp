#include "baselines/baseline_factory.h"

#include <sstream>

namespace sjs {
namespace baselines {

std::string BaselineHelpSupportedDims() {
  std::ostringstream oss;
  oss << "  Supported dimensions: 2 / 3 / 4 / 5\n"
      << "  Methods:\n"
      << "    - ours\n"
      << "    - kd_tree\n"
      << "\n  Variants:\n"
      << "    - sampling\n"
      << "    - enum_sampling\n"
      << "\n  Supported combinations:\n"
      << "    - ours/sampling      : SJS Framework II (2D keeps the optimized fast path; 3D/4D/5D use fixed-mode local recursive structures aligned with the outlines)\n"
      << "    - ours/enum_sampling : SJS Framework I\n"
      << "    - kd_tree/sampling   : lifted KD-tree baseline on 2d dimensions\n";
  return oss.str();
}

}  // namespace baselines
}  // namespace sjs
