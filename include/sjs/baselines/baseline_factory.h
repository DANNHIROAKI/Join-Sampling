#pragma once
// sjs/baselines/baseline_factory.h
//
// Dimension-generic baseline factory declarations.
//
// This header only declares CreateBaseline<Dim>().
// The definition is provided in a single TU:
//   src/baselines/baseline_factory_dispatch.cpp
// which includes concrete baseline implementations and explicitly instantiates
// CreateBaseline<Dim>() for the supported Dim list (see sjs/dispatch/dim_dispatch.h).
//
// This design avoids including all baseline headers in every app TU.

#include "sjs/baselines/baseline_api.h"
#include "sjs/baselines/baseline_names.h"
#include "sjs/core/types.h"

#include <memory>
#include <sstream>
#include <string>
#include <string_view>

namespace sjs {
namespace baselines {

// Primary factory (DECLARATION ONLY).
// Defined + explicitly instantiated in src/baselines/baseline_factory_dispatch.cpp.
template <int Dim, class T = Scalar>
std::unique_ptr<IBaseline<Dim, T>> CreateBaseline(Method method,
                                                  Variant variant,
                                                  std::string* err = nullptr);

// Convenience overload: create baseline from a config.
template <int Dim, class T = Scalar>
inline std::unique_ptr<IBaseline<Dim, T>> CreateBaseline(const sjs::Config& cfg,
                                                         std::string* err = nullptr) {
  return CreateBaseline<Dim, T>(cfg.run.method, cfg.run.variant, err);
}

// Convenience overload: parse from strings.
template <int Dim, class T = Scalar>
inline std::unique_ptr<IBaseline<Dim, T>> CreateBaseline(std::string_view method,
                                                         std::string_view variant,
                                                         std::string* err = nullptr) {
  Method m = Method::Unknown;
  Variant v = Variant::Sampling;

  if (!ParseMethod(method, &m) || m == Method::Unknown) {
    if (err) {
      std::ostringstream oss;
      oss << "Unknown method: '" << method << "'.\n";
      oss << "Available baselines:\n" << BaselineHelp();
      *err = oss.str();
    }
    return nullptr;
  }
  if (!ParseVariant(variant, &v)) {
    if (err) {
      std::ostringstream oss;
      oss << "Unknown variant: '" << variant << "'.\n";
      oss << "Allowed variants: sampling | enum_sampling | adaptive.\n";
      oss << "Available baselines:\n" << BaselineHelp();
      *err = oss.str();
    }
    return nullptr;
  }

  return CreateBaseline<Dim, T>(m, v, err);
}

}  // namespace baselines
}  // namespace sjs
