// src/baselines/pbsm/adaptive.cpp
//
/* Explicit template instantiation for pbsm::PBSMAdaptiveBaseline in 2D.
 *
 * This project keeps baseline implementations as header-only templates for:
 *   - easy Dim extension (2D now, higher-D later)
 *   - convenient inlining for performance experiments
 *
 * By explicitly instantiating the Dim=2 specialization here, we centralize
 * code generation for the common case.
 */
#include "sjs/baselines/pbsm/adaptive.h"

namespace sjs::baselines::pbsm {

template class PBSMAdaptiveBaseline<2, sjs::Scalar>;

}  // namespace sjs::baselines::pbsm
