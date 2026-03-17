// src/geometry/instantiations_2d.cpp
//
// Optional explicit instantiations for common Dim=2 geometry templates.
//
// Most geometry code is header-only templates. Explicit instantiation here is
// optional, but helps reduce compilation overhead when many translation units
// use the same Dim=2 specializations.

#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/predicates.h"
#include "geometry/embedding.h"

#include "core/types.h"

namespace sjs {

// Common primitives
template struct Point<2, Scalar>;
template struct Box<2, Scalar>;

// Embedding bounds helper used by KD/SIRS-style baselines.
template struct DomainBounds<2, Scalar>;

// Common free-function instantiation (optional).
template Scalar IntersectionVolume<2, Scalar>(const Box<2, Scalar>&, const Box<2, Scalar>&) noexcept;

}  // namespace sjs
