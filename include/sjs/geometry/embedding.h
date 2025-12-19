#pragma once
// sjs/geometry/embedding.h
//
// Embeddings that reduce box intersection conditions to orthogonal range
// conditions over points.
//
// Motivation:
//  - Several baselines (KD-tree / RangeTree / SIRS-style range sampling)
//    work over points with orthogonal constraints.
//  - A common trick is to embed each box r into a point p(r) = [lo(r), hi(r)]
//    in 2*Dim dimensions, then express intersection with query box q as
//    constraints on p(r):
//
//    For half-open boxes:
//      r intersects q  <=>  for all i:
//        lo_i(r) < hi_i(q)    AND    lo_i(q) < hi_i(r)
//
//    The first part constrains lo(r) by an upper bound; the second constrains
//    hi(r) by a lower bound. This is an orthogonal "2*Dim-dimensional range"
//    with some open bounds.
//  - In practice, range structures often implement half-open boxes [L,U).
//    We convert strict bounds to half-open by using nextafter().
//
// This header provides:
//  - DomainBounds<Dim>: finite coordinate bounds for turning half-infinite
//    queries into finite ones.
//  - EmbedLowerUpper(): box -> point in 2*Dim dims.
//  - MakeIntersectQueryRange(): query box -> half-open query range in embedded space.
//  - SkipDim0 variants (common when Dim0 is handled by sweep events).

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/geometry/box.h"

#include <cmath>
#include <limits>
#include <type_traits>

namespace sjs {

template <int Dim, class T = Scalar>
struct DomainBounds {
  static_assert(Dim >= 1, "DomainBounds<Dim>: Dim must be >= 1");
  static_assert(Dim <= kMaxSupportedDim, "DomainBounds<Dim>: Dim too large");

  Point<Dim, T> min_lo = Point<Dim, T>::Filled(std::numeric_limits<T>::infinity());
  Point<Dim, T> max_hi = Point<Dim, T>::Filled(-std::numeric_limits<T>::infinity());

  void ExpandToInclude(const Box<Dim, T>& b) noexcept {
    for (int i = 0; i < Dim; ++i) {
      const usize idx = static_cast<usize>(i);
      if (b.lo.v[idx] < min_lo.v[idx]) min_lo.v[idx] = b.lo.v[idx];
      if (b.hi.v[idx] > max_hi.v[idx]) max_hi.v[idx] = b.hi.v[idx];
    }
  }

  bool IsInitialized() const noexcept {
    for (int i = 0; i < Dim; ++i) {
      const usize idx = static_cast<usize>(i);
      if (!(min_lo.v[idx] <= max_hi.v[idx])) return false;
    }
    return true;
  }

  static DomainBounds FromBoxes(Span<const Box<Dim, T>> boxes) noexcept {
    DomainBounds d;
    if (boxes.empty()) {
      // Reasonable default for empty input.
      d.min_lo = Point<Dim, T>::Zero();
      d.max_hi = Point<Dim, T>::Zero();
      return d;
    }
    for (const auto& b : boxes) d.ExpandToInclude(b);
    // If somehow still uninitialized, fall back.
    if (!d.IsInitialized()) {
      d.min_lo = Point<Dim, T>::Zero();
      d.max_hi = Point<Dim, T>::Zero();
    }
    return d;
  }

  // Merge two bounds.
  static DomainBounds Merge(const DomainBounds& a, const DomainBounds& b) noexcept {
    DomainBounds d;
    for (int i = 0; i < Dim; ++i) {
      const usize idx = static_cast<usize>(i);
      d.min_lo.v[idx] = (a.min_lo.v[idx] < b.min_lo.v[idx]) ? a.min_lo.v[idx] : b.min_lo.v[idx];
      d.max_hi.v[idx] = (a.max_hi.v[idx] > b.max_hi.v[idx]) ? a.max_hi.v[idx] : b.max_hi.v[idx];
    }
    return d;
  }
};

// Next representable value greater than x.
template <class T>
inline T NextUp(T x) noexcept {
  if constexpr (std::numeric_limits<T>::has_infinity) {
    return std::nextafter(x, std::numeric_limits<T>::infinity());
  } else {
    return std::nextafter(x, std::numeric_limits<T>::max());
  }
}

// Next representable value smaller than x.
template <class T>
inline T NextDown(T x) noexcept {
  if constexpr (std::numeric_limits<T>::has_infinity) {
    return std::nextafter(x, -std::numeric_limits<T>::infinity());
  } else {
    return std::nextafter(x, std::numeric_limits<T>::lowest());
  }
}

template <int Dim, class T = Scalar>
using EmbeddedPoint = Point<2 * Dim, T>;

template <int Dim, class T = Scalar>
using EmbeddedBox = Box<2 * Dim, T>;

// Embed a box into a 2*Dim point: [lo0..lo(Dim-1), hi0..hi(Dim-1)].
template <int Dim, class T>
inline EmbeddedPoint<Dim, T> EmbedLowerUpper(const Box<Dim, T>& b) noexcept {
  EmbeddedPoint<Dim, T> p;
  for (int i = 0; i < Dim; ++i) {
    const usize idx = static_cast<usize>(i);
    const usize dim_usize = static_cast<usize>(Dim);
    p.v[idx] = b.lo.v[idx];
    p.v[idx + dim_usize] = b.hi.v[idx];
  }
  return p;
}

// Build a finite half-open query range in embedded space that captures the
// intersection predicate under half-open semantics.
//
// Given domain bounds (min_lo, max_hi), we create:
//
//   For lo-coordinates:   lo(r)_i in [min_lo_i, hi(q)_i)
//   For hi-coordinates:   hi(r)_i in (lo(q)_i, max_hi_i]    (strict lower bound)
//
// We convert the strict lower bound "hi(r)_i > lo(q)_i" into
// "hi(r)_i >= nextafter(lo(q)_i, +inf)".
// We also convert the inclusive upper bound "hi(r)_i <= max_hi_i" into a
// half-open upper bound "< nextafter(max_hi_i, +inf)" to ensure inclusion.
template <int Dim, class T>
inline EmbeddedBox<Dim, T> MakeIntersectQueryRange(const Box<Dim, T>& q,
                                                   const DomainBounds<Dim, T>& domain) noexcept {
  EmbeddedBox<Dim, T> r;
  for (int i = 0; i < Dim; ++i) {
    const usize idx = static_cast<usize>(i);
    const usize dim_usize = static_cast<usize>(Dim);
    // lo(r)_i < hi(q)_i
    r.lo.v[idx] = domain.min_lo.v[idx];
    r.hi.v[idx] = q.hi.v[idx];

    // hi(r)_i > lo(q)_i  ==> hi(r)_i >= nextafter(lo(q)_i, +inf)
    r.lo.v[idx + dim_usize] = NextUp(q.lo.v[idx]);
    // hi(r)_i <= max_hi_i  ==> hi(r)_i < nextafter(max_hi_i, +inf)
    r.hi.v[idx + dim_usize] = NextUp(domain.max_hi.v[idx]);
  }
  return r;
}

// ---------- SkipDim0 variants ----------
// Common when Dim0 (e.g., x) is handled by sweep and always overlaps.
// These embed only dimensions 1..Dim-1, producing 2*(Dim-1) space.

template <int Dim, class T = Scalar>
using EmbeddedPointSkipDim0 = Point<2 * (Dim - 1), T>;

template <int Dim, class T = Scalar>
using EmbeddedBoxSkipDim0 = Box<2 * (Dim - 1), T>;

template <int Dim, class T>
inline EmbeddedPointSkipDim0<Dim, T> EmbedLowerUpperSkipDim0(const Box<Dim, T>& b) noexcept {
  static_assert(Dim >= 2, "EmbedLowerUpperSkipDim0 requires Dim >= 2");
  EmbeddedPointSkipDim0<Dim, T> p;
  const usize k = static_cast<usize>(Dim - 1);
  for (int j = 1; j < Dim; ++j) {
    const usize idx = static_cast<usize>(j - 1);
    const usize j_idx = static_cast<usize>(j);
    p.v[idx] = b.lo.v[j_idx];
    p.v[idx + k] = b.hi.v[j_idx];
  }
  return p;
}

template <int Dim, class T>
inline EmbeddedBoxSkipDim0<Dim, T> MakeIntersectQueryRangeSkipDim0(
    const Box<Dim, T>& q, const DomainBounds<Dim, T>& domain) noexcept {
  static_assert(Dim >= 2, "MakeIntersectQueryRangeSkipDim0 requires Dim >= 2");
  EmbeddedBoxSkipDim0<Dim, T> r;
  const usize k = static_cast<usize>(Dim - 1);
  for (int j = 1; j < Dim; ++j) {
    const usize idx = static_cast<usize>(j - 1);
    const usize j_idx = static_cast<usize>(j);
    // lo(r)_j < hi(q)_j
    r.lo.v[idx] = domain.min_lo.v[j_idx];
    r.hi.v[idx] = q.hi.v[j_idx];

    // hi(r)_j > lo(q)_j
    r.lo.v[idx + k] = NextUp(q.lo.v[j_idx]);
    r.hi.v[idx + k] = NextUp(domain.max_hi.v[j_idx]);
  }
  return r;
}

}  // namespace sjs
