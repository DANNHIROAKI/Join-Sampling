// src/data/synthetic/instantiations_2d.cpp
//
// Optional explicit instantiations for Dim=2 synthetic generator classes.
//
// This TU centralizes code generation for the 2D specializations of the
// synthetic generators. It is not strictly required (templates would be
// instantiated wherever used), but helps reduce compile time and object bloat
// in large builds.

#include "sjs/core/types.h"

// IMPORTANT: include generator.h (not individual generators) to avoid a
// circular-include edge case where including a generator header directly would
// pull in generator.h, which in turn defines the built-in factory.
//
// Including generator.h first ensures all built-in generators are visible
// before the factory is instantiated.
#include "sjs/data/synthetic/generator.h"

namespace sjs {
namespace synthetic {

template class StripeCtrlAlphaGenerator<2, Scalar>;
template class UniformGenerator<2, Scalar>;
template class ClusteredGenerator<2, Scalar>;
template class HeteroSizesGenerator<2, Scalar>;

}  // namespace synthetic
}  // namespace sjs
