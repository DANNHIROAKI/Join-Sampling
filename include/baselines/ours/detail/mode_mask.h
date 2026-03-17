#pragma once
// baselines/ours/detail/mode_mask.h
//
// Fixed-mode helpers for the high-dimensional SJS extension.
// After sweeping axis 0, the remaining Dim-1 projection dimensions induce
// 2^(Dim-1) mutually exclusive modes using the A/B lower-endpoint rule.

#include "core/types.h"
#include "geometry/box.h"

#include <array>
#include <limits>

namespace sjs {
namespace baselines {
namespace ours {
namespace detail {

template <int Dim>
inline constexpr u32 kModeCountV = static_cast<u32>(1u << (Dim - 1));

template <int Dim>
using ModeWeights = std::array<u64, kModeCountV<Dim>>;

template <int Dim, class T>
inline bool ComputeModeMaskSkipAxis0(const Box<Dim, T>& q,
                                     const Box<Dim, T>& other,
                                     u32* out_mode_mask) {
  static_assert(Dim >= 2, "ComputeModeMaskSkipAxis0 requires Dim >= 2");
  if (!out_mode_mask) return false;

  u32 mask = 0;
  for (int axis = 1; axis < Dim; ++axis) {
    const T q_lo = q.lo.v[static_cast<usize>(axis)];
    const T q_hi = q.hi.v[static_cast<usize>(axis)];
    const T o_lo = other.lo.v[static_cast<usize>(axis)];
    const T o_hi = other.hi.v[static_cast<usize>(axis)];

    // Half-open overlap on this projection axis.
    if (!(o_lo < q_hi && q_lo < o_hi)) return false;

    // Pattern A (bit=0): other.lo <= q.lo < other.hi
    // Pattern B (bit=1): q.lo < other.lo < q.hi
    if (o_lo <= q_lo) {
      // Ties go to A. Overlap check above guarantees q_lo < other.hi.
      continue;
    }

    // Because intervals overlap, q.lo < other.lo implies other.lo < q.hi.
    mask |= static_cast<u32>(1u << static_cast<u32>(axis - 1));
  }

  *out_mode_mask = mask;
  return true;
}

}  // namespace detail
}  // namespace ours
}  // namespace baselines
}  // namespace sjs
