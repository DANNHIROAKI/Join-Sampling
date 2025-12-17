// src/baselines/interval_tree/sampling.cpp
//
/* Explicit template instantiation for interval_tree::IntervalTreeSamplingBaseline in 2D.
 *
 * This project keeps baseline implementations as header-only templates for:
 *   - easy Dim extension (2D now, higher-D later)
 *   - convenient inlining for performance experiments
 *
 * By explicitly instantiating the Dim=2 specialization here, we centralize
 * code generation for the common case.
 */
#include "sjs/baselines/interval_tree/sampling.h"

namespace sjs::baselines::interval_tree {

template class IntervalTreeSamplingBaseline<2, sjs::Scalar>;

}  // namespace sjs::baselines::interval_tree
