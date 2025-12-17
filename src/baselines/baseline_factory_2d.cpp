// src/baselines/baseline_factory_2d.cpp
//
// Implementation of the Dim=2 baseline factory.
//
// This translation unit intentionally includes *all* baseline headers.
// Doing so keeps the app-layer include surface small (apps only include
// baseline_factory_2d.h), at the cost of compiling this TU a bit heavier.

#include "baselines/baseline_factory_2d.h"

#include "sjs/core/types.h"  // ParseMethod/ParseVariant

// --------------------------
// Baseline implementations (each provides 3 variants)
// --------------------------
#include "sjs/baselines/ours/adaptive.h"
#include "sjs/baselines/ours/enum_sampling.h"
#include "sjs/baselines/ours/sampling.h"

#include "sjs/baselines/aabb/adaptive.h"
#include "sjs/baselines/aabb/enum_sampling.h"
#include "sjs/baselines/aabb/sampling.h"

#include "sjs/baselines/interval_tree/adaptive.h"
#include "sjs/baselines/interval_tree/enum_sampling.h"
#include "sjs/baselines/interval_tree/sampling.h"

#include "sjs/baselines/kd_tree/adaptive.h"
#include "sjs/baselines/kd_tree/enum_sampling.h"
#include "sjs/baselines/kd_tree/sampling.h"

#include "sjs/baselines/r_tree/adaptive.h"
#include "sjs/baselines/r_tree/enum_sampling.h"
#include "sjs/baselines/r_tree/sampling.h"

#include "sjs/baselines/range_tree/adaptive.h"
#include "sjs/baselines/range_tree/enum_sampling.h"
#include "sjs/baselines/range_tree/sampling.h"

#include "sjs/baselines/pbsm/adaptive.h"
#include "sjs/baselines/pbsm/enum_sampling.h"
#include "sjs/baselines/pbsm/sampling.h"

#include "sjs/baselines/tlsop/adaptive.h"
#include "sjs/baselines/tlsop/enum_sampling.h"
#include "sjs/baselines/tlsop/sampling.h"

#include "sjs/baselines/sirs/adaptive.h"
#include "sjs/baselines/sirs/enum_sampling.h"
#include "sjs/baselines/sirs/sampling.h"

#include "sjs/baselines/rejection/adaptive.h"
#include "sjs/baselines/rejection/enum_sampling.h"
#include "sjs/baselines/rejection/sampling.h"

#include "sjs/baselines/tsunami/adaptive.h"
#include "sjs/baselines/tsunami/enum_sampling.h"
#include "sjs/baselines/tsunami/sampling.h"

#include <sstream>

namespace sjs {
namespace baselines {

namespace {

inline std::unique_ptr<IBaseline<2, Scalar>> MakeUnsupported(std::string_view method,
                                                            std::string_view variant,
                                                            std::string* err) {
  if (err) {
    std::ostringstream oss;
    oss << "Unsupported baseline for Dim=2: method='" << method << "' variant='" << variant << "'.\n";
    oss << "Available baselines:\n" << BaselineHelp2D();
    *err = oss.str();
  }
  return nullptr;
}

}  // namespace

std::unique_ptr<IBaseline<2, Scalar>> CreateBaseline2D(Method method,
                                                       Variant variant,
                                                       std::string* err) {
  // We intentionally keep this as a switch (fast, explicit, compiler-friendly).
  // It must stay consistent with BaselineRegistry2D() in baseline_names.cpp.

  switch (method) {
    case Method::Ours:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<ours::OursSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<ours::OursEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<ours::OursAdaptiveBaseline<2, Scalar>>();
      }
      break;

    case Method::AABB:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<aabb::AABBSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<aabb::AABBEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<aabb::AABBAdaptiveBaseline<2, Scalar>>();
      }
      break;

    case Method::IntervalTree:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<interval_tree::IntervalTreeSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<interval_tree::IntervalTreeEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<interval_tree::IntervalTreeAdaptiveBaseline<2, Scalar>>();
      }
      break;

    case Method::KDTree:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<kd_tree::KDTreeSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<kd_tree::KDTreeEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<kd_tree::KDTreeAdaptiveBaseline<2, Scalar>>();
      }
      break;

    case Method::RTree:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<r_tree::RTreeSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<r_tree::RTreeEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<r_tree::RTreeAdaptiveBaseline<2, Scalar>>();
      }
      break;

    case Method::RangeTree:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<range_tree::RangeTreeSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<range_tree::RangeTreeEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<range_tree::RangeTreeAdaptiveBaseline<2, Scalar>>();
      }
      break;

    case Method::PBSM:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<pbsm::PBSMSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<pbsm::PBSMEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<pbsm::PBSMAdaptiveBaseline<2, Scalar>>();
      }
      break;

    case Method::TLSOP:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<tlsop::TLSOPSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<tlsop::TLSOPEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<tlsop::TLSOPAdaptiveBaseline<2, Scalar>>();
      }
      break;

    case Method::SIRS:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<sirs::SIRSSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<sirs::SIRSEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<sirs::SIRSAdaptiveBaseline<2, Scalar>>();
      }
      break;

    case Method::Rejection:
      switch (variant) {
        case Variant::Sampling:
          return std::make_unique<rejection::RejectionSamplingBaseline<2, Scalar>>();
        case Variant::EnumSampling:
          return std::make_unique<rejection::RejectionEnumSamplingBaseline<2, Scalar>>();
        case Variant::Adaptive:
          return std::make_unique<rejection::RejectionAdaptiveBaseline<2, Scalar>>();
      }
      break;

    case Method::Unknown:
      // fallthrough
    default:
      break;
  }

  return MakeUnsupported(ToString(method), ToString(variant), err);
}

std::unique_ptr<IBaseline<2, Scalar>> CreateBaseline2D(std::string_view method,
                                                       std::string_view variant,
                                                       std::string* err) {
  Method m = Method::Unknown;
  Variant v = Variant::Sampling;

  if (!ParseMethod(method, &m) || m == Method::Unknown) {
    if (err) {
      std::ostringstream oss;
      oss << "Unknown method: '" << method << "'.\n";
      oss << "Available baselines:\n" << BaselineHelp2D();
      *err = oss.str();
    }
    return nullptr;
  }
  if (!ParseVariant(variant, &v)) {
    if (err) {
      std::ostringstream oss;
      oss << "Unknown variant: '" << variant << "'.\n";
      oss << "Allowed variants: sampling | enum_sampling | adaptive.\n";
      oss << "Available baselines:\n" << BaselineHelp2D();
      *err = oss.str();
    }
    return nullptr;
  }

  return CreateBaseline2D(m, v, err);
}

}  // namespace baselines
}  // namespace sjs
