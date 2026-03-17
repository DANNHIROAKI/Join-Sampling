#pragma once
// core/dim_dispatch.h
//
// Small runtime -> compile-time dimensionality dispatcher used by apps and
// factories. Supported experiment dimensions are intentionally limited to the
// research scope of this repository: 2D/3D/4D/5D.

#include "core/types.h"

#include <string>

namespace sjs {

inline constexpr bool IsSupportedExperimentDim(int dim) noexcept {
  return dim >= 2 && dim <= 5;
}

inline std::string SupportedDimError(int dim) {
  return "Supported dimensions are 2/3/4/5; got --dim=" + std::to_string(dim);
}

template <class Fn>
decltype(auto) DispatchSupportedDim(int dim, Fn&& fn) {
  switch (dim) {
    case 2: return fn(DimTag<2>{});
    case 3: return fn(DimTag<3>{});
    case 4: return fn(DimTag<4>{});
    case 5: return fn(DimTag<5>{});
    default:
      // Caller should guard with IsSupportedExperimentDim / SupportedDimError.
      return fn(DimTag<2>{});
  }
}

}  // namespace sjs
