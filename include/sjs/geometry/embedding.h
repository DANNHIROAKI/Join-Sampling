#pragma once
// sjs/geometry/embedding.h
//
// Embeddings for reducing box intersection conditions to orthogonal range
// conditions over points.
//
// Key identity for half-open boxes (SJS-HighDims):
//   r intersects q  <=>  for all i:
//     lo(r)_i < hi(q)_i   AND   lo(q)_i < hi(r)_i
//
// Embed each box r into a point p(r) = [lo(r), hi(r)] in 2*Dim dimensions.
// Then intersect(q) becomes an orthogonal range query in embedded space with
// open/closed bounds converted into half-open bounds.
//
// Also provides SkipDim0 variants (for sweep on axis 0): embed dims 1..Dim-1 only.

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/geometry/box.h"

#include <cmath>
#include <limits>
#include <type_traits>

namespace sjs {

// Domain bounds used to turn half-infinite query constraints into finite ones.
template <int Dim, class T = Scalar>
struct DomainBounds {
  static_assert(Dim >= 1, "DomainBounds<Dim>: Dim must be >= 1");
  static_assert(Dim <= kMaxSupportedDim, "DomainBounds<Dim>: Dim too large");

  // Initialize to extreme finite values (works for both floats and integrals).
  Point<Dim, T> min_lo = Point<Dim, T>::Filled(std::numeric_limits<T>::max());
  Point<Dim, T> max_hi = Point<Dim, T>::Filled(std::numeric_limits<T>::lowest());

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
      d.min_lo = Point<Dim, T>::Zero();
      d.max_hi = Point<Dim, T>::Zero();
      return d;
    }
    for (const auto& b : boxes) d.ExpandToInclude(b);
    if (!d.IsInitialized()) {
      d.min_lo = Point<Dim, T>::Zero();
      d.max_hi = Point<Dim, T>::Zero();
    }
    return d;
  }

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

// Next representable > x.
template <class T>
inline T NextUp(T x) noexcept {
  if constexpr (std::is_floating_point_v<T>) {
    if constexpr (std::numeric_limits<T>::has_infinity) {
      return std::nextafter(x, std::numeric_limits<T>::infinity());
    } else {
      return std::nextafter(x, std::numeric_limits<T>::max());
    }
  } else if constexpr (std::is_integral_v<T>) {
    return (x == std::numeric_limits<T>::max()) ? x : static_cast<T>(x + 1);
  } else {
    return x;
  }
}

// Next representable < x.
template <class T>
inline T NextDown(T x) noexcept {
  if constexpr (std::is_floating_point_v<T>) {
    if constexpr (std::numeric_limits<T>::has_infinity) {
      return std::nextafter(x, -std::numeric_limits<T>::infinity());
    } else {
      return std::nextafter(x, std::numeric_limits<T>::lowest());
    }
  } else if constexpr (std::is_integral_v<T>) {
    return (x == std::numeric_limits<T>::lowest()) ? x : static_cast<T>(x - 1);
  } else {
    return x;
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
  const usize k = static_cast<usize>(Dim);
  for (int i = 0; i < Dim; ++i) {
    const usize idx = static_cast<usize>(i);
    p.v[idx] = b.lo.v[idx];
    p.v[idx + k] = b.hi.v[idx];
  }
  return p;
}

// Build a finite half-open query range in embedded space that captures the
// intersection predicate under half-open semantics.
//
// For each i:
//   lo(r)_i < hi(q)_i      ==> lo(r)_i in [domain.min_lo_i, hi(q)_i)
//   hi(r)_i > lo(q)_i      ==> hi(r)_i in [NextUp(lo(q)_i), ...)
//   hi(r)_i <= domain.max_hi_i  ==> hi(r)_i < NextUp(domain.max_hi_i)
template <int Dim, class T>
inline EmbeddedBox<Dim, T> MakeIntersectQueryRange(const Box<Dim, T>& q,
                                                   const DomainBounds<Dim, T>& domain) noexcept {
  EmbeddedBox<Dim, T> r;
  const usize k = static_cast<usize>(Dim);
  for (int i = 0; i < Dim; ++i) {
    const usize idx = static_cast<usize>(i);

    // lo(r)_i in [min_lo_i, hi(q)_i)
    r.lo.v[idx] = domain.min_lo.v[idx];
    r.hi.v[idx] = q.hi.v[idx];

    // hi(r)_i in [NextUp(lo(q)_i), NextUp(max_hi_i))
    r.lo.v[idx + k] = NextUp(q.lo.v[idx]);
    r.hi.v[idx + k] = NextUp(domain.max_hi.v[idx]);
  }
  return r;
}

// ---------- SkipDim0 variants ----------
// Used when axis 0 is already handled by sweep (always overlaps in that axis),
// so we only need dimensions 1..Dim-1 in the embedded space.

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

    // lo(r)_j in [min_lo_j, hi(q)_j)
    r.lo.v[idx] = domain.min_lo.v[j_idx];
    r.hi.v[idx] = q.hi.v[j_idx];

    // hi(r)_j in [NextUp(lo(q)_j), NextUp(max_hi_j))
    r.lo.v[idx + k] = NextUp(q.lo.v[j_idx]);
    r.hi.v[idx + k] = NextUp(domain.max_hi.v[j_idx]);
  }
  return r;
}

}  // namespace sjs
