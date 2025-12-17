// src/sampling/instantiations_2d.cpp
//
// Optional explicit instantiations for sampling helpers when used with
// (Dim=2) join streams.
//
// The codebase is largely header-only. This TU exists as a convenient place to
// centralize explicit instantiations to reduce compile times in large builds.

#include "sjs/sampling/rank_sampling.h"

#include "sjs/join/join_types.h"

#include "sjs/core/types.h"

#include <string>
#include <vector>

namespace sjs {
namespace sampling {

// Explicit instantiation of two-pass rank sampling when the stream is accessed
// through the polymorphic join::IJoinStream interface.
//
// This is the most common usage in Enum+Sampling runners:
//   std::unique_ptr<join::IJoinStream> stream = ...;
//   RankSampleWithReplacement(stream.get(), t, &rng, &samples, &info, &err);

template bool RankSampleWithReplacement<join::IJoinStream, PairId>(
    join::IJoinStream* stream,
    u64 t,
    Rng* rng,
    std::vector<PairId>* out_samples,
    RankSamplingInfo* info,
    std::string* err);

template bool CountStreamItems<join::IJoinStream, PairId>(
    join::IJoinStream* stream,
    u64* out_n,
    std::string* err);

}  // namespace sampling
}  // namespace sjs
