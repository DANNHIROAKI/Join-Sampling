// src/baselines/range_tree/enum_sampling.cpp
//
/* Explicit template instantiation for range_tree::RangeTreeEnumSamplingBaseline in 2D.
 *
 * This project keeps baseline implementations as header-only templates for:
 *   - easy Dim extension (2D now, higher-D later)
 *   - convenient inlining for performance experiments
 *
 * By explicitly instantiating the Dim=2 specialization here, we centralize
 * code generation for the common case.
 */
#include "sjs/baselines/range_tree/enum_sampling.h"

namespace sjs::baselines::range_tree {

template class RangeTreeEnumSamplingBaseline<2, sjs::Scalar>;

}  // namespace sjs::baselines::range_tree
