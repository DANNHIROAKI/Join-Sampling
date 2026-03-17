#include "baselines/ours/sampling.h"

namespace sjs::baselines::ours {

template class OursSamplingBaseline<3, sjs::Scalar>;
template class OursSamplingBaseline<4, sjs::Scalar>;
template class OursSamplingBaseline<5, sjs::Scalar>;

}  // namespace sjs::baselines::ours
