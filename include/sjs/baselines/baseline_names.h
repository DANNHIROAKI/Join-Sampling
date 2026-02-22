#pragma once
// sjs/baselines/baseline_names.h
//
// Baseline registry declarations (names + help strings).
//
// The actual registry list is defined in:
//   src/baselines/baseline_names.cpp
//
// Keeping the list centralized avoids inconsistencies between:
//   - CLI help
//   - factory supported combinations
//   - tests that iterate over all baselines

#include "sjs/core/types.h"

#include <string>
#include <string_view>

namespace sjs {
namespace baselines {

// A lightweight registry row.
struct BaselineSpec {
  Method method = Method::Unknown;
  Variant variant = Variant::Sampling;

  // Canonical key used in logs/CSV (stable across versions).
  // Convention: "<method>/<variant>" e.g., "ours/sampling".
  std::string_view key;

  // Short human-readable description.
  std::string_view desc;
};

// Return the list of baselines supported by this build.
//
// The returned Span points to static storage (no allocations).
Span<const BaselineSpec> BaselineRegistry() noexcept;

// Convenience predicate.
bool IsBaselineSupported(Method method, Variant variant) noexcept;

// Build a help string suitable for --help output / error messages.
std::string BaselineHelp();

}  // namespace baselines
}  // namespace sjs
