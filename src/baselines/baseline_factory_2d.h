#pragma once
// src/baselines/baseline_factory_2d.h
//
// Factory + registry helpers for Dim=2 baselines.
//
// Design goals:
//  - Apps should only include this header (and baseline_api.h transitively),
//    not individual baseline implementations.
//  - The corresponding .cpp includes all baseline headers and constructs the
//    requested implementation.
//  - Centralize the "what baselines exist" list for CLI/help.
//
// Note: This header lives under src/ (not include/). In CMake, add ${PROJECT_SOURCE_DIR}/src
// to your include directories so apps can:  #include "baselines/baseline_factory_2d.h".

#include "sjs/baselines/baseline_api.h"

#include <memory>
#include <string>
#include <string_view>

namespace sjs {
namespace baselines {

// A lightweight registry row.
struct BaselineSpec2D {
  Method method = Method::Unknown;
  Variant variant = Variant::Sampling;

  // Canonical key used in logs/CSV (stable across versions).
  // Convention: "<method>/<variant>" e.g., "ours/sampling".
  std::string_view key;

  // A short human-readable description.
  std::string_view desc;
};

// Return the list of baselines supported by this build for Dim=2.
//
// The returned Span points to static storage (no allocations).
Span<const BaselineSpec2D> BaselineRegistry2D() noexcept;

// Convenience predicates.
bool IsBaselineSupported2D(Method method, Variant variant) noexcept;

// Build a help string suitable for --help output.
std::string BaselineHelp2D();

// --------------------------
// Factory
// --------------------------
// Create a baseline implementation for Dim=2.
//
// Returns nullptr on error and sets *err.
std::unique_ptr<IBaseline<2, Scalar>> CreateBaseline2D(Method method,
                                                       Variant variant,
                                                       std::string* err = nullptr);

// Convenience overloads.
inline std::unique_ptr<IBaseline<2, Scalar>> CreateBaseline2D(const Config& cfg,
                                                              std::string* err = nullptr) {
  return CreateBaseline2D(cfg.run.method, cfg.run.variant, err);
}

std::unique_ptr<IBaseline<2, Scalar>> CreateBaseline2D(std::string_view method,
                                                       std::string_view variant,
                                                       std::string* err = nullptr);

}  // namespace baselines
}  // namespace sjs
