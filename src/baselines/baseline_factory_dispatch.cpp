// src/baselines/baseline_factory_dispatch.cpp
//
// Dim-generic baseline factory implementation + explicit instantiations.
// This TU intentionally includes baseline implementation headers.

#include "sjs/baselines/baseline_factory.h"

#include "sjs/baselines/baseline_names.h"
#include "sjs/core/types.h"
#include "sjs/dispatch/dim_dispatch.h"

// --------------------------
// Baseline implementations
// --------------------------
#include "sjs/baselines/ours/highdims/sampling.h"
#include "sjs/baselines/ours/highdims/enum_sampling.h"
#include "sjs/baselines/ours/highdims/adaptive.h"

#include "sjs/baselines/range_tree/sampling.h"
#include "sjs/baselines/range_tree/enum_sampling.h"
#include "sjs/baselines/range_tree/adaptive.h"

#include "sjs/baselines/kd_tree/sampling.h"

#include <sstream>
#include <string_view>
#include <utility>

namespace sjs {
namespace baselines {

namespace {

inline std::string MakeUnsupportedMsg(int dim, std::string_view method, std::string_view variant) {
  std::ostringstream oss;
  oss << "Unsupported baseline: dim=" << dim
      << " method='" << method << "' variant='" << variant << "'.\n\n";
  oss << "Available baselines:\n" << BaselineHelp();
  return oss.str();
}

}  // namespace

template <int Dim, class T>
std::unique_ptr<IBaseline<Dim, T>> CreateBaseline(Method method, Variant variant, std::string* err) {
  // IMPORTANT: Keep consistent with BaselineRegistry() in baseline_names.cpp.

  switch (method) {
    case Method::Ours:
      switch (variant) {
        case Variant::Sampling:
          // Our main HighDims implementation (Framework II)
          return std::make_unique<ours::highdims::OursHighDimsSamplingBaseline<Dim, T>>();
        case Variant::EnumSampling:
          // Framework I
          return std::make_unique<ours::highdims::OursHighDimsEnumSamplingBaseline<Dim, T>>();
        case Variant::Adaptive:
          // Framework III
          return std::make_unique<ours::highdims::OursHighDimsAdaptiveBaseline<Dim, T>>();
        default:
          break;
      }
      break;

    case Method::RangeTree:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<range_tree::RangeTreeSamplingBaseline<Dim, T>>();
        case Variant::EnumSampling:
          return std::make_unique<range_tree::RangeTreeEnumSamplingBaseline<Dim, T>>();
        case Variant::Adaptive:
          return std::make_unique<range_tree::RangeTreeAdaptiveBaseline<Dim, T>>();
        default:
          break;
      }
      break;

    case Method::KDTree:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<kd_tree::KdTreeSamplingBaseline<Dim, T>>();
        default:
          break;
      }
      break;

    case Method::Unknown:
    default:
      break;
  }

  if (err) *err = MakeUnsupportedMsg(Dim, ToString(method), ToString(variant));
  return nullptr;
}

// --------------------------
// Explicit instantiations for Scalar over supported dims
// --------------------------

// SJS_DIMS is the single source of truth for which compile-time dimensions are
// supported. It is defined (with a default list) in:
//   include/sjs/dispatch/dim_dispatch.h
//
// IMPORTANT: Do not re-define a separate default list here, otherwise it's easy
// for the explicit instantiation set to diverge from runtime dispatch.
#ifndef SJS_DIMS
#error "SJS_DIMS is not defined. Include sjs/dispatch/dim_dispatch.h (or define SJS_DIMS) before compiling this TU."
#endif

#define SJS_INSTANTIATE_BASELINE_FACTORY(D)                                                     \
  template std::unique_ptr<IBaseline<D, Scalar>> CreateBaseline<D, Scalar>(Method, Variant, std::string*);

SJS_DIMS(SJS_INSTANTIATE_BASELINE_FACTORY)

#undef SJS_INSTANTIATE_BASELINE_FACTORY

}  // namespace baselines
}  // namespace sjs