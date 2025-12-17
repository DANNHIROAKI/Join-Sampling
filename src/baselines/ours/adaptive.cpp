// src/baselines/ours/adaptive.cpp
//
// Explicit template instantiation for our method (Adaptive variant) in 2D.

#include "sjs/baselines/ours/adaptive.h"

namespace sjs::baselines::ours {

template class OursAdaptiveBaseline<2, sjs::Scalar>;

}  // namespace sjs::baselines::ours
