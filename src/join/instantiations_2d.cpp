// src/join/instantiations_2d.cpp
//
// Optional explicit instantiations for Dim=2 join components.
//
// In this project, many algorithms are header-only templates.
// Explicit instantiations can reduce compile times and keep object code
// size under control when many translation units use the same template
// instantiations.
//
// NOTE:
//  - We do not use 'extern template' declarations in headers (to keep
//    headers simple), so these instantiations are primarily a convenience
//    for libraries/binaries that link this TU.

#include "sjs/join/join_enumerator.h"
#include "sjs/join/join_oracle.h"

#include "sjs/core/types.h"

namespace sjs {
namespace join {

// Deterministic plane-sweep join stream for (Dim=2, Scalar).
template class PlaneSweepJoinStream<2, Scalar>;

// Naive streaming join oracle for (Dim=2, Scalar).
template class NaiveJoinStream<2, Scalar>;

}  // namespace join
}  // namespace sjs
