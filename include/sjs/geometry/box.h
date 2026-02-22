#pragma once
// sjs/geometry/box.h
//
// Axis-aligned half-open boxes: [lo, hi) in each dimension.
//
// Semantics (SJS-HighDims):
//  - Point containment: lo[i] <= x[i] < hi[i] for all i.
//  - Box intersection:  loA[i] < hiB[i] AND loB[i] < hiA[i] for all i.
//    (strict < due to half-open upper bound)
//
// Box does not auto-normalize. Use IsValid()/IsEmpty() as needed.

#include "sjs/geometry/point.h"

#include <cmath>
#include <limits>
#include <ostream>
#include <sstream>
#include <type_traits>

namespace sjs {
namespace detail {

// Next representable value greater than x (used to maintain half-open inclusions).
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

}  // namespace detail

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
    // Canonical empty box: lo=+INF (or max), hi=-INF (or lowest).
    // Works for both floating and integral T.
    const T big = std::numeric_limits<T>::max();
    const T low = std::numeric_limits<T>::lowest();
    return Box(Point<Dim, T>::Filled(big), Point<Dim, T>::Filled(low));
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

  // Proper: strictly positive width in all dims.
  constexpr bool IsProper() const noexcept { return !IsEmpty(); }

  // -------- geometry helpers --------
  constexpr T Width(int axis) const noexcept {
    SJS_DASSERT(axis >= 0 && axis < Dim);
    return hi.v[static_cast<usize>(axis)] - lo.v[static_cast<usize>(axis)];
  }

  // Volume (area in 2D).
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

  // Expand this box to include a point under half-open semantics.
  // Ensure lo<=p and hi>p in every dimension where p is extreme.
  void ExpandToIncludePoint(const Point<Dim, T>& p) noexcept {
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      if (p.v[i] < lo.v[i]) lo.v[i] = p.v[i];
      if (p.v[i] >= hi.v[i]) hi.v[i] = detail::NextUp(p.v[i]);
    }
  }

  // Expand this box to include another box (bounding box).
  void ExpandToIncludeBox(const Box& b) noexcept {
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      if (b.lo.v[i] < lo.v[i]) lo.v[i] = b.lo.v[i];
      if (b.hi.v[i] > hi.v[i]) hi.v[i] = b.hi.v[i];
    }
  }

  // -------- membership tests (half-open) --------
  constexpr bool Contains(const Point<Dim, T>& p) const noexcept {
    if (IsEmpty()) return false;
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      const T x = p.v[i];
      if (x < lo.v[i]) return false;
      if (!(x < hi.v[i])) return false;
    }
    return true;
  }

  // Half-open box containment: inner ⊆ outer iff lo_out<=lo_in and hi_in<=hi_out.
  constexpr bool ContainsBox(const Box& inner) const noexcept {
    if (inner.IsEmpty()) return true;
    if (IsEmpty()) return false;
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      if (inner.lo.v[i] < lo.v[i]) return false;
      if (inner.hi.v[i] > hi.v[i]) return false;
    }
    return true;
  }

  // Half-open intersection.
  constexpr bool Intersects(const Box& b) const noexcept {
    if (IsEmpty() || b.IsEmpty()) return false;
    for (usize i = 0; i < static_cast<usize>(Dim); ++i) {
      if (!(lo.v[i] < b.hi.v[i] && b.lo.v[i] < hi.v[i])) return false;
    }
    return true;
  }

  // Intersection box; returns Empty() if disjoint.
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
