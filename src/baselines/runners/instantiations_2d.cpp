// src/baselines/runners/instantiations_2d.cpp
//
// Optional explicit template instantiations for runner functions at Dim=2.
//
// Why this file exists:
//   - Runner functions are templates defined in headers.
//   - If your build uses many apps/translation units, explicit instantiation
//     can reduce compile time *if* you also add corresponding 'extern template'
//     declarations in the runner headers.
//
// This file is safe to keep even without 'extern template' declarations; it
// will just instantiate these functions once here (compilers usually de-dup
// template instantiations via COMDAT/weak ODR).

#include "sjs/baselines/runners/adaptive_runner.h"
#include "sjs/baselines/runners/enum_sampling_runner.h"
#include "sjs/baselines/runners/sampling_runner.h"

namespace sjs {
namespace baselines {

// Dim=2 explicit instantiations

template bool RunSamplingOnce<2, Scalar>(IBaseline<2, Scalar>* baseline,
                                        const Dataset<2, Scalar>& dataset,
                                        const Config& cfg,
                                        u64 seed,
                                        RunReport* out,
                                        std::string* err);

template bool RunEnumSamplingOnce<2, Scalar>(IBaseline<2, Scalar>* baseline,
                                            const Dataset<2, Scalar>& dataset,
                                            const Config& cfg,
                                            u64 seed,
                                            RunReport* out,
                                            std::string* err);

template bool RunAdaptiveOnce<2, Scalar>(IBaseline<2, Scalar>* baseline,
                                        const Dataset<2, Scalar>& dataset,
                                        const Config& cfg,
                                        u64 seed,
                                        RunReport* out,
                                        std::string* err);

}  // namespace baselines
}  // namespace sjs
