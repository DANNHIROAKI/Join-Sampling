// src/baselines/ours/enum_sampling.cpp
//
// Explicit template instantiation for our method (Enumerate+Sampling variant)
// in 2D.

#include "baselines/ours/enum_sampling.h"

namespace sjs::baselines::ours {

template class OursEnumSamplingBaseline<2, sjs::Scalar>;

}  // namespace sjs::baselines::ours
