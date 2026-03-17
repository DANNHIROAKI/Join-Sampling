#include "baselines/ours/enum_sampling.h"

namespace sjs::baselines::ours {

template class OursEnumSamplingBaseline<3, sjs::Scalar>;
template class OursEnumSamplingBaseline<4, sjs::Scalar>;
template class OursEnumSamplingBaseline<5, sjs::Scalar>;

}  // namespace sjs::baselines::ours
