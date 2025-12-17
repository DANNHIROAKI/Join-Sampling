// src/baselines/rejection/adaptive.cpp
//
/* Explicit template instantiation for rejection::RejectionAdaptiveBaseline in 2D.
 *
 * This project keeps baseline implementations as header-only templates for:
 *   - easy Dim extension (2D now, higher-D later)
 *   - convenient inlining for performance experiments
 *
 * By explicitly instantiating the Dim=2 specialization here, we centralize
 * code generation for the common case.
 */
#include "sjs/baselines/rejection/adaptive.h"

namespace sjs::baselines::rejection {

template class RejectionAdaptiveBaseline<2, sjs::Scalar>;

}  // namespace sjs::baselines::rejection
