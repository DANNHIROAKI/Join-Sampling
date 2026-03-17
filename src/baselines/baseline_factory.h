#pragma once
// src/baselines/baseline_factory.h
//
// Generic 2D/3D/4D/5D baseline factory. The original Dim=2 wrapper factory is
// kept intact for backward compatibility and tests; new apps should prefer this
// header.

#include "baselines/baseline_api.h"
#include "baselines/ours/enum_sampling.h"
#include "baselines/ours/sampling.h"
#include "core/types.h"
#include "third_party/RS-over-SRJ/isrjs_kds/isrjs_kds.hpp"

#include <memory>
#include <sstream>
#include <string>
#include <string_view>

namespace sjs {
namespace baselines {

std::string BaselineHelpSupportedDims();

namespace detail {

inline std::unique_ptr<IBaseline<2, Scalar>> MakeUnsupported2D(std::string_view method,
                                                               std::string_view variant,
                                                               std::string* err) {
  if (err) {
    std::ostringstream oss;
    oss << "Unsupported baseline: method='" << method << "' variant='" << variant << "'.\n"
        << BaselineHelpSupportedDims();
    *err = oss.str();
  }
  return nullptr;
}

template <int Dim, class T = Scalar>
inline std::unique_ptr<IBaseline<Dim, T>> MakeUnsupported(std::string_view method,
                                                          std::string_view variant,
                                                          std::string* err) {
  if (err) {
    std::ostringstream oss;
    oss << "Unsupported baseline: method='" << method << "' variant='" << variant
        << "' for Dim=" << Dim << ".\n"
        << BaselineHelpSupportedDims();
    *err = oss.str();
  }
  return nullptr;
}

}  // namespace detail

template <int Dim, class T = Scalar>
inline std::unique_ptr<IBaseline<Dim, T>> CreateBaseline(Method method,
                                                         Variant variant,
                                                         std::string* err = nullptr) {
  switch (method) {
    case Method::Ours:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<ours::OursSamplingBaseline<Dim, T>>();
        case Variant::EnumSampling:
          return std::make_unique<ours::OursEnumSamplingBaseline<Dim, T>>();
        default:
          break;
      }
      break;
    case Method::RSOverSRJ:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<rs_over_srj::RSOverSRJKDTreeSamplingBaseline<Dim, T>>();
        default:
          break;
      }
      break;
    case Method::Unknown:
    default:
      break;
  }
  return detail::MakeUnsupported<Dim, T>(ToString(method), ToString(variant), err);
}

template <int Dim, class T = Scalar>
inline std::unique_ptr<IBaseline<Dim, T>> CreateBaseline(const Config& cfg,
                                                         std::string* err = nullptr) {
  return CreateBaseline<Dim, T>(cfg.run.method, cfg.run.variant, err);
}

template <int Dim, class T = Scalar>
inline std::unique_ptr<IBaseline<Dim, T>> CreateBaseline(std::string_view method,
                                                         std::string_view variant,
                                                         std::string* err = nullptr) {
  Method m = Method::Unknown;
  Variant v = Variant::Sampling;
  if (!ParseMethod(method, &m) || m == Method::Unknown) {
    if (err) {
      std::ostringstream oss;
      oss << "Unknown method: '" << method << "'.\n" << BaselineHelpSupportedDims();
      *err = oss.str();
    }
    return nullptr;
  }
  if (!ParseVariant(variant, &v)) {
    if (err) {
      std::ostringstream oss;
      oss << "Unknown variant: '" << variant << "'.\n" << BaselineHelpSupportedDims();
      *err = oss.str();
    }
    return nullptr;
  }
  return CreateBaseline<Dim, T>(m, v, err);
}

}  // namespace baselines
}  // namespace sjs
