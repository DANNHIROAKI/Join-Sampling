// src/data/instantiations_2d.cpp
//
// Optional explicit instantiations for Dim=2 data-generation entry points.
//
// Most data generation is implemented as header-only templates (see
// include/sjs/data/synthetic/*). This TU centralizes explicit instantiations for
// the most common Dim=2 usage to reduce compile time in large builds.

#include "sjs/core/types.h"
#include "sjs/io/dataset.h"
#include "sjs/data/synthetic/generator.h"

#include <memory>
#include <string>
#include <string_view>

namespace sjs {
namespace synthetic {

// Explicit instantiation of the high-level generation function.
template bool GenerateDataset<2, Scalar>(std::string_view generator_name,
                                        const DatasetSpec& spec,
                                        Dataset<2, Scalar>* out_ds,
                                        Report* report,
                                        std::string* err);

// Explicit instantiation of the generator factory.
template std::unique_ptr<ISyntheticGenerator<2, Scalar>> Create<2, Scalar>(std::string_view name);

}  // namespace synthetic
}  // namespace sjs
