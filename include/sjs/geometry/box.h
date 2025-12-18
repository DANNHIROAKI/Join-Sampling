#pragma once
// sjs/geometry/box.h
//
// Axis-aligned half-open boxes: [lo, hi) in each dimension.
//
// This half-open convention is important for deterministic spatial joins
// without boundary duplicates (e.g., when using plane sweep tie-break rules).
//
// Key semantics:
//  - A point x is inside box b iff for all i: b.lo[i] <= x[i] < b.hi[i].
//  - Two boxes a and b intersect iff for all i: a.lo[i] < b.hi[i] AND b.lo[i] < a.hi[i].
//    (strict < due to half-open upper bound)
//
// Box stores lo and hi as Points. No automatic normalization is performed; use
// IsValid() / IsEmpty() checks and construct carefully.

#include "sjs/geometry/point.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <sstream>

namespace sjs {

template <int Dim, class T = Scalar>
struct Box {
  static_assert(Dim >= 1, "Box<Dim>: Dim must be >= 1");
  static_assert(Dim <= kMaxSupportedDim, "Box<Dim>: Dim too large");

  using value_type = T;
  static constexpr int kDim = Dim;

  Point<Dim, T> lo{};
  Point<Dim, T> hi{};

  // -------- constructors --------
  constexpr Box() noexcept = default;

  constexpr Box(const Point<Dim, T>& lo_, const Point<Dim, T>& hi_) noexcept : lo(lo_), hi(hi_) {}

  static Box Empty() noexcept {
    // A canonical "empty" box: lo=+inf, hi=-inf.
    // This makes Union/Expand patterns simple.
    const T pinf = std::numeric_limits<T>::infinity();
    const T ninf = -std::numeric_limits<T>::infinity();
    return Box(Point<Dim, T>::Filled(pinf), Point<Dim, T>::Filled(ninf));
  }

  static Box FromCenterSize(const Point<Dim, T>& center, const Point<Dim, T>& size) noexcept {
    Point<Dim, T> lo_, hi_;
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      const T half = size.v[i] / T(2);
      lo_.v[i] = center.v[i] - half;
      hi_.v[i] = center.v[i] + half;
    }
    return Box(lo_, hi_);
  }

  // -------- validity / emptiness --------
  // Valid: lo <= hi for all dims (allows zero width).
  constexpr bool IsValid() const noexcept {
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      if (lo.v[i] > hi.v[i]) return false;
    }
    return true;
  }

  // Empty in half-open sense: any dimension has lo >= hi.
  constexpr bool IsEmpty() const noexcept {
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      if (!(lo.v[i] < hi.v[i])) return true;
    }
    return false;
  }

  // Proper: strictly positive width in all dims (non-empty).
  constexpr bool IsProper() const noexcept { return !IsEmpty(); }

  // -------- geometry helpers --------
  constexpr T Width(int axis) const noexcept {
    SJS_DASSERT(axis >= 0 && axis < Dim);
    return hi.v[static_cast<usize>(axis)] - lo.v[static_cast<usize>(axis)];
  }

  // Volume in Dim dimensions (area in 2D, length in 1D).
  T Volume() const noexcept {
    if (IsEmpty()) return T(0);
    T vol = T(1);
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      const T w = hi.v[i] - lo.v[i];
      if (!(w > T(0))) return T(0);
      vol *= w;
    }
    return vol;
  }

  Point<Dim, T> Center() const noexcept {
    Point<Dim, T> c;
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) c.v[i] = (lo.v[i] + hi.v[i]) / T(2);
    return c;
  }

  // Expand this box to include a point (treating the point as an infinitesimal).
  // For half-open semantics, we include p by ensuring lo<=p and hi>p.
  // We use nextafter(p, +inf) to keep half-open consistent.
  void ExpandToIncludePoint(const Point<Dim, T>& p) noexcept {
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      if (p.v[i] < lo.v[i]) lo.v[i] = p.v[i];
      if (p.v[i] >= hi.v[i]) {
        if constexpr (std::numeric_limits<T>::has_infinity) {
          hi.v[i] = std::nextafter(p.v[i], std::numeric_limits<T>::infinity());
        } else {
          hi.v[i] = std::nextafter(p.v[i], std::numeric_limits<T>::max());
        }
      }
    }
  }

  // Expand this box to include another box (union / bounding box).
  void ExpandToIncludeBox(const Box& b) noexcept {
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      if (b.lo.v[i] < lo.v[i]) lo.v[i] = b.lo.v[i];
      if (b.hi.v[i] > hi.v[i]) hi.v[i] = b.hi.v[i];
    }
  }

  // -------- membership tests (half-open) --------
  constexpr bool Contains(const Point<Dim, T>& p) const noexcept {
    // Empty box contains nothing.
    if (IsEmpty()) return false;
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      const T x = p.v[i];
      if (x < lo.v[i]) return false;
      if (!(x < hi.v[i])) return false;  // half-open upper bound
    }
    return true;
  }

  // Outer contains inner iff every point in inner is also in outer.
  // For half-open, this is lo_out <= lo_in and hi_in <= hi_out.
  constexpr bool ContainsBox(const Box& inner) const noexcept {
    if (inner.IsEmpty()) return true;  // empty is contained by any box
    if (IsEmpty()) return false;
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      if (inner.lo.v[i] < lo.v[i]) return false;
      if (inner.hi.v[i] > hi.v[i]) return false;
    }
    return true;
  }

  // Intersection test (half-open). See predicates.h for free function.
  constexpr bool Intersects(const Box& b) const noexcept {
    // Quick reject empties.
    if (IsEmpty() || b.IsEmpty()) return false;
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      if (!(lo.v[i] < b.hi.v[i] && b.lo.v[i] < hi.v[i])) return false;
    }
    return true;
  }

  // Compute intersection box (half-open). Returns Empty() if no intersection.
  Box Intersection(const Box& b) const noexcept {
    if (!Intersects(b)) return Empty();
    Point<Dim, T> lo2, hi2;
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      lo2.v[i] = (lo.v[i] > b.lo.v[i]) ? lo.v[i] : b.lo.v[i];
      hi2.v[i] = (hi.v[i] < b.hi.v[i]) ? hi.v[i] : b.hi.v[i];
    }
    Box out(lo2, hi2);
    if (out.IsEmpty()) return Empty();
    return out;
  }

  std::string ToString() const {
    std::ostringstream oss;
    oss << "[" << lo << " , " << hi << ")";
    return oss.str();
  }
};

template <int Dim, class T>
inline std::ostream& operator<<(std::ostream& os, const Box<Dim, T>& b) {
  os << b.ToString();
  return os;
}

}  // namespace sjs
