#pragma once
// sjs/geometry/point.h
//
// A small fixed-size point/vector type used across geometry and indices.
// Designed for compile-time dimensionality (Dim) and memory efficiency.
//
// Notes:
//  - Uses std::array<T, Dim> for contiguous storage.
//  - Provides basic arithmetic and component-wise operations.
//  - Intended to be trivially copyable when T is trivially copyable.
//
// Half-open box semantics are defined in box.h/predicates.h; Point itself is
// purely numeric.

#include "sjs/core/types.h"
#include "sjs/core/assert.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <ostream>
#include <sstream>
#include <type_traits>

namespace sjs {

template <int Dim, class T = Scalar>
struct Point {
  static_assert(Dim >= 1, "Point<Dim>: Dim must be >= 1");
  static_assert(Dim <= kMaxSupportedDim, "Point<Dim>: Dim too large");

  using value_type = T;
  static constexpr int kDim = Dim;
  static constexpr usize kDimUsize = static_cast<usize>(Dim);

  // Suppress sign-conversion warning: Dim is compile-time constant and always positive
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsign-conversion"
  std::array<T, Dim> v{};
  #pragma GCC diagnostic pop

  // -------- constructors --------
  constexpr Point() noexcept = default;

  // Suppress sign-conversion warning: Dim is compile-time constant and always positive
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsign-conversion"
  explicit constexpr Point(const std::array<T, Dim>& a) noexcept : v(a) {}
  #pragma GCC diagnostic pop

  // Construct from initializer list (size must match Dim in debug builds).
  explicit Point(std::initializer_list<T> init) {
    SJS_DASSERT_MSG(static_cast<int>(init.size()) == Dim,
                    "Point initializer_list size must equal Dim");
    usize i = 0;
    for (T x : init) {
      if (i < kDimUsize) v[i++] = x;
    }
    // If init shorter in release, remaining are zero-initialized.
  }

  // Fill all coordinates with a constant.
  static constexpr Point Filled(T x) noexcept {
    Point p;
    for (usize i = 0; i < kDimUsize; ++i) p.v[i] = x;
    return p;
  }

  static constexpr Point Zero() noexcept { return Filled(T(0)); }

  // -------- element access --------
  constexpr T& operator[](int i) noexcept {
    SJS_DASSERT(i >= 0 && i < Dim);
    return v[static_cast<usize>(i)];
  }
  constexpr const T& operator[](int i) const noexcept {
    SJS_DASSERT(i >= 0 && i < Dim);
    return v[static_cast<usize>(i)];
  }

  constexpr T* data() noexcept { return v.data(); }
  constexpr const T* data() const noexcept { return v.data(); }

  constexpr auto begin() noexcept { return v.begin(); }
  constexpr auto end() noexcept { return v.end(); }
  constexpr auto begin() const noexcept { return v.begin(); }
  constexpr auto end() const noexcept { return v.end(); }

  // -------- comparisons --------
  friend constexpr bool operator==(const Point& a, const Point& b) noexcept {
    for (usize i = 0; i < kDimUsize; ++i) {
      if (!(a.v[i] == b.v[i])) return false;
    }
    return true;
  }
  friend constexpr bool operator!=(const Point& a, const Point& b) noexcept { return !(a == b); }

  // Lexicographic order (useful for sorting).
  friend constexpr bool operator<(const Point& a, const Point& b) noexcept {
    for (usize i = 0; i < kDimUsize; ++i) {
      if (a.v[i] < b.v[i]) return true;
      if (b.v[i] < a.v[i]) return false;
    }
    return false;
  }

  // -------- arithmetic --------
  friend constexpr Point operator+(const Point& a, const Point& b) noexcept {
    Point out;
    for (usize i = 0; i < kDimUsize; ++i) out.v[i] = a.v[i] + b.v[i];
    return out;
  }
  friend constexpr Point operator-(const Point& a, const Point& b) noexcept {
    Point out;
    for (usize i = 0; i < kDimUsize; ++i) out.v[i] = a.v[i] - b.v[i];
    return out;
  }
  friend constexpr Point operator*(const Point& a, T s) noexcept {
    Point out;
    for (usize i = 0; i < kDimUsize; ++i) out.v[i] = a.v[i] * s;
    return out;
  }
  friend constexpr Point operator*(T s, const Point& a) noexcept { return a * s; }
  friend constexpr Point operator/(const Point& a, T s) noexcept {
    Point out;
    for (usize i = 0; i < kDimUsize; ++i) out.v[i] = a.v[i] / s;
    return out;
  }

  Point& operator+=(const Point& b) noexcept {
    for (usize i = 0; i < kDimUsize; ++i) v[i] += b.v[i];
    return *this;
  }
  Point& operator-=(const Point& b) noexcept {
    for (usize i = 0; i < kDimUsize; ++i) v[i] -= b.v[i];
    return *this;
  }
  Point& operator*=(T s) noexcept {
    for (usize i = 0; i < kDimUsize; ++i) v[i] *= s;
    return *this;
  }
  Point& operator/=(T s) noexcept {
    for (usize i = 0; i < kDimUsize; ++i) v[i] /= s;
    return *this;
  }

  // -------- component-wise ops --------
  static constexpr Point Min(const Point& a, const Point& b) noexcept {
    Point out;
    for (usize i = 0; i < kDimUsize; ++i) out.v[i] = (a.v[i] < b.v[i]) ? a.v[i] : b.v[i];
    return out;
  }
  static constexpr Point Max(const Point& a, const Point& b) noexcept {
    Point out;
    for (usize i = 0; i < kDimUsize; ++i) out.v[i] = (a.v[i] < b.v[i]) ? b.v[i] : a.v[i];
    return out;
  }

  // L-infinity norm (max abs coordinate), helpful for quick magnitude checks.
  T NormInf() const noexcept {
    T m = T(0);
    for (usize i = 0; i < kDimUsize; ++i) {
      const T a = (v[i] < T(0)) ? -v[i] : v[i];
      if (a > m) m = a;
    }
    return m;
  }

  std::string ToString() const {
    std::ostringstream oss;
    oss << "(";
    for (usize i = 0; i < kDimUsize; ++i) {
      oss << v[i];
      if (i + 1 < kDimUsize) oss << ",";
    }
    oss << ")";
    return oss.str();
  }
};

template <int Dim, class T>
inline std::ostream& operator<<(std::ostream& os, const Point<Dim, T>& p) {
  os << p.ToString();
  return os;
}

}  // namespace sjs
