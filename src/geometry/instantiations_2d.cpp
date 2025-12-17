// src/geometry/instantiations_2d.cpp
//
// Optional explicit instantiations for common Dim=2 geometry templates.
//
// Most geometry code is header-only templates. Explicit instantiation here is
// optional, but helps reduce compilation overhead when many translation units
// use the same Dim=2 specializations.

#include "sjs/geometry/point.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/predicates.h"
#include "sjs/geometry/embedding.h"

#include "sjs/core/types.h"

namespace sjs {

// Common primitives
template struct Point<2, Scalar>;
template struct Box<2, Scalar>;

// Embedding bounds helper used by KD/RangeTree/SIRS-style baselines.
template struct DomainBounds<2, Scalar>;

// Common free-function instantiation (optional).
template Scalar IntersectionVolume<2, Scalar>(const Box<2, Scalar>&, const Box<2, Scalar>&) noexcept;

}  // namespace sjs
