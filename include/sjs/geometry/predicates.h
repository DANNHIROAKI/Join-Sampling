#pragma once
// sjs/geometry/predicates.h
//
// Common predicates for half-open boxes and points.
// We keep them as free functions so algorithms can use them without
// depending on Box's member methods (although Box also provides them).

#include "sjs/geometry/box.h"

#include <algorithm>

namespace sjs {

// 1D half-open interval intersection: [a0,a1) intersects [b0,b1) iff a0 < b1 && b0 < a1.
template <class T>
constexpr bool Intersects1DHalfOpen(T a0, T a1, T b0, T b1) noexcept {
  return (a0 < b1) && (b0 < a1);
}

// 1D half-open containment: x in [a0,a1) iff a0 <= x < a1.
template <class T>
constexpr bool Contains1DHalfOpen(T a0, T a1, T x) noexcept {
  return (a0 <= x) && (x < a1);
}

// Box/box intersection under half-open semantics.
template <int Dim, class T>
constexpr bool IntersectsHalfOpen(const Box<Dim, T>& a, const Box<Dim, T>& b) noexcept {
  return a.Intersects(b);
}

// Box/point containment under half-open semantics.
template <int Dim, class T>
constexpr bool ContainsHalfOpen(const Box<Dim, T>& b, const Point<Dim, T>& p) noexcept {
  return b.Contains(p);
}

// Box contains box (half-open).
template <int Dim, class T>
constexpr bool ContainsBoxHalfOpen(const Box<Dim, T>& outer, const Box<Dim, T>& inner) noexcept {
  return outer.ContainsBox(inner);
}

// Intersection volume (0 if disjoint).
template <int Dim, class T>
T IntersectionVolume(const Box<Dim, T>& a, const Box<Dim, T>& b) noexcept {
  const Box<Dim, T> inter = a.Intersection(b);
  return inter.Volume();
}

// Strict overlap test (equivalent to IntersectsHalfOpen).
// Provided as semantic sugar: "overlaps with positive measure".
// For half-open boxes, Intersects already implies positive measure in all dims.
template <int Dim, class T>
constexpr bool Overlaps(const Box<Dim, T>& a, const Box<Dim, T>& b) noexcept {
  return a.Intersects(b);
}

}  // namespace sjs
