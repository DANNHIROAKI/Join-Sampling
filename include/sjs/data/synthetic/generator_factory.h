#pragma once
// sjs/data/synthetic/generator_factory.h
//
// Factory + convenience wrappers for synthetic generators.
//
// This file centralizes the built-in generators:
//  - stripe_ctrl_alpha
//  - uniform
//  - clustered
//  - hetero_sizes
//
// Also provides a runtime-dim wrapper (uses sjs/dispatch/dim_dispatch.h).

#include "sjs/data/synthetic/generator.h"

// Built-in generators
#include "sjs/data/synthetic/stripe_ctrl_alpha.h"
#include "sjs/data/synthetic/uniform.h"
#include "sjs/data/synthetic/clustered.h"
#include "sjs/data/synthetic/hetero_sizes.h"

// Python-backed generator (alacarte-rectgen)
#include "sjs/data/synthetic/alacarte_rectgen.h"

// Optional (but recommended): runtime dim dispatch helper
#include "sjs/dispatch/dim_dispatch.h"

#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace sjs {
namespace synthetic {

// Create generator by name for a fixed compile-time Dim.
template <int Dim, class T = Scalar>
inline std::unique_ptr<ISyntheticGenerator<Dim, T>> CreateSyntheticGenerator(std::string_view generator_name) {
  using detail::EqualsIgnoreCase;

  // Accept some short aliases for convenience.
  if (EqualsIgnoreCase(generator_name, "stripe_ctrl_alpha") ||
      EqualsIgnoreCase(generator_name, "stripe") ||
      EqualsIgnoreCase(generator_name, "scirg")) {
    return std::make_unique<StripeCtrlAlphaGenerator<Dim, T>>();
  }
  if (EqualsIgnoreCase(generator_name, "uniform") ||
      EqualsIgnoreCase(generator_name, "uni")) {
    return std::make_unique<UniformGenerator<Dim, T>>();
  }
  if (EqualsIgnoreCase(generator_name, "clustered") ||
      EqualsIgnoreCase(generator_name, "clusters") ||
      EqualsIgnoreCase(generator_name, "hotspot")) {
    return std::make_unique<ClusteredGenerator<Dim, T>>();
  }
  if (EqualsIgnoreCase(generator_name, "hetero_sizes") ||
      EqualsIgnoreCase(generator_name, "hetero") ||
      EqualsIgnoreCase(generator_name, "mix_sizes")) {
    return std::make_unique<HeteroSizesGenerator<Dim, T>>();
  }

  if (EqualsIgnoreCase(generator_name, "alacarte_rectgen") ||
      EqualsIgnoreCase(generator_name, "alacarte-rectgen") ||
      EqualsIgnoreCase(generator_name, "rectgen") ||
      EqualsIgnoreCase(generator_name, "alacarte")) {
    return std::make_unique<AlacarteRectGenGenerator<Dim, T>>();
  }

  return nullptr;
}

// Convenience one-shot API (compile-time Dim).
template <int Dim, class T = Scalar>
inline bool GenerateDataset(std::string_view generator_name,
                            const DatasetSpec& spec,
                            Dataset<Dim, T>* out_ds,
                            Report* report = nullptr,
                            std::string* err = nullptr) {
  auto gen = CreateSyntheticGenerator<Dim, T>(generator_name);
  if (!gen) {
    if (err) *err = "Unknown generator: " + std::string(generator_name);
    return false;
  }
  return gen->Generate(spec, out_ds, report, err);
}

// Runtime Dim wrapper.
// Usage pattern (in apps):
//   Report rep;
//   std::string err;
//   DatasetSpec spec = ...
//   bool ok = synthetic::DispatchGenerateByDim<Scalar>(
//       dim, gen_name, spec,
//       [&](auto& ds){ /* ds is Dataset<Dim,Scalar>& */ write(ds); return true; },
//       &rep, &err);
//
// Fn signature: bool Fn(Dataset<Dim,T>&)

// Helper functor for dim_dispatch (namespace-scope: local classes cannot have
// member templates in C++17).
template <class T, class Fn>
struct DispatchGenerateByDimRunner {
  std::string_view generator_name;
  const DatasetSpec* spec = nullptr;
  Fn* fn = nullptr;
  Report* report = nullptr;
  std::string* err = nullptr;

  template <int Dim>
  bool operator()() const {
    Dataset<Dim, T> ds;
    if (!GenerateDataset<Dim, T>(generator_name, *spec, &ds, report, err)) return false;
    return (*fn)(ds);
  }
};

template <class T = Scalar, class Fn>
inline bool DispatchGenerateByDim(int dim,
                                  std::string_view generator_name,
                                  const DatasetSpec& spec,
                                  Fn&& fn,
                                  Report* report = nullptr,
                                  std::string* err = nullptr) {
  bool ok = false;
  std::string derr;
  if (!sjs::dispatch::DispatchDim<bool>(
          dim,
          DispatchGenerateByDimRunner<T, std::remove_reference_t<Fn>>{generator_name, &spec, &fn, report, err},
          &ok,
          &derr)) {
    if (err) *err = derr;
    return false;
  }
  return ok;
}

}  // namespace synthetic
}  // namespace sjs
