// src/baselines/ours/sampling.cpp
//
// Explicit template instantiation for our method (Sampling variant) in 2D.
//
// Why this file exists
// --------------------
// The baseline implementations are header-only templates (for easy extension
// to higher dimensions). In a larger build, that can lead to repeated template
// instantiations across translation units. By explicitly instantiating the
// 2D specialization here, we:
//   - centralize codegen for Dim=2
//   - (optionally) allow other TUs to use `extern template` declarations later
//
// NOTE:
//   The rest of the system (e.g., src/baselines/baseline_factory_2d.cpp)
//   can still include the headers; this TU simply ensures a concrete
//   instantiation is available in the final binary.

#include "sjs/baselines/ours/sampling.h"

namespace sjs::baselines::ours {

// Instantiate the baseline on the project-wide scalar type (double).
template class OursSamplingBaseline<2, sjs::Scalar>;

}  // namespace sjs::baselines::ours
