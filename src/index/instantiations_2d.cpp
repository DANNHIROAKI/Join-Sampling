// src/index/instantiations_2d.cpp
//
// Optional explicit instantiations for common Dim=2 index templates.
//
// Most index structures in this repository are header-only templates.
// This TU centralizes explicit instantiations for Dim=2 to reduce compile time
// and potential object-code bloat when many translation units use the same
// instantiations.

#include "sjs/index/aabb_tree.h"
#include "sjs/index/interval_tree.h"
#include "sjs/index/kd_tree.h"
#include "sjs/index/r_tree.h"
#include "sjs/index/range_tree.h"

#include "sjs/core/types.h"

namespace sjs {
namespace index {

template class AABBTree<2, Scalar>;
template class RTree<2, Scalar>;
template class KDTree<2, Scalar>;
template class RangeTree<2, Scalar>;
template class IntervalTreeOnAxis<2, Scalar>;

}  // namespace index
}  // namespace sjs
