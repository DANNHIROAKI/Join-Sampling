#pragma once
// sjs/core/types.h
//
// Core, dependency-light types shared across the project.
// The project is currently focused on 2D, but these types are designed
// to scale to higher dimensions and very large datasets.

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ostream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace sjs {

// --------------------------
// Fixed-width integer aliases
// --------------------------
using i8  = std::int8_t;
using i16 = std::int16_t;
using i32 = std::int32_t;
using i64 = std::int64_t;

using u8  = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

using usize = std::size_t;

// --------------------------
// Scalar type (geometry uses this)
// --------------------------
// Use double by default for robustness; can be changed globally if needed.
using Scalar = double;

// --------------------------
// Object identifiers
// --------------------------
// 32-bit saves memory and is enough up to 4,294,967,295 objects per relation.
// If you need more, change to u64 project-wide.
using Id = u32;
inline constexpr Id kInvalidId = std::numeric_limits<Id>::max();

// --------------------------
// Dimensionality helpers
// --------------------------
inline constexpr int kDefaultDim = 2;
inline constexpr int kMaxSupportedDim = 32;  // soft cap; used for sanity checks

template <int Dim>
struct DimTag {
  static_assert(Dim >= 1, "Dim must be >= 1");
  static_assert(Dim <= kMaxSupportedDim, "Dim too large");
  static constexpr int value = Dim;
};

// --------------------------
// Pair identifiers (join output refers to (r_id, s_id))
// --------------------------
struct PairId {
  Id r{};
  Id s{};

  constexpr PairId() = default;
  constexpr PairId(Id r_, Id s_) : r(r_), s(s_) {}

  friend constexpr bool operator==(const PairId& a, const PairId& b) noexcept {
    return a.r == b.r && a.s == b.s;
  }
  friend constexpr bool operator!=(const PairId& a, const PairId& b) noexcept {
    return !(a == b);
  }
};

struct PairIdHash {
  usize operator()(const PairId& p) const noexcept {
    // 32-bit ids packed into 64-bit, then hashed.
    // This is stable and fast; sufficient for experiments.
    const u64 key = (static_cast<u64>(p.r) << 32) ^ static_cast<u64>(p.s);
    // Splitmix64 finalizer for better bit diffusion.
    u64 x = key + 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x = x ^ (x >> 31);
    return static_cast<usize>(x);
  }
};

// --------------------------
// Experiment selectors
// --------------------------
enum class Variant : u8 {
  Sampling = 0,
  EnumSampling = 1,
  Adaptive = 2,
};

enum class Method : u8 {
  Ours = 0,
  AABB = 1,
  IntervalTree = 2,
  KDTree = 3,
  RTree = 4,
  RangeTree = 5,
  PBSM = 6,
  TLSOP = 7,
  SIRS = 8,
  Rejection = 9,
  Tsunami = 10,   // <-- add
  Unknown = 255,
};


inline constexpr std::string_view ToString(Variant v) noexcept {
  switch (v) {
    case Variant::Sampling: return "sampling";
    case Variant::EnumSampling: return "enum_sampling";
    case Variant::Adaptive: return "adaptive";
  }
  return "unknown";
}

inline constexpr std::string_view ToString(Method m) noexcept {
  switch (m) {
    case Method::Ours: return "ours";
    case Method::AABB: return "aabb";
    case Method::IntervalTree: return "interval_tree";
    case Method::KDTree: return "kd_tree";
    case Method::RTree: return "r_tree";
    case Method::RangeTree: return "range_tree";
    case Method::PBSM: return "pbsm";
    case Method::TLSOP: return "tlsop";
    case Method::SIRS: return "sirs";
    case Method::Rejection: return "rejection";
    case Method::Unknown: return "unknown";
    case Method::Tsunami: return "tsunami";
  }
  return "unknown";
}

namespace detail {
inline constexpr char LowerAscii(char c) noexcept {
  return (c >= 'A' && c <= 'Z') ? static_cast<char>(c - 'A' + 'a') : c;
}
inline bool EqualsIgnoreCase(std::string_view a, std::string_view b) noexcept {
  if (a.size() != b.size()) return false;
  for (usize i = 0; i < a.size(); ++i) {
    if (LowerAscii(a[i]) != LowerAscii(b[i])) return false;
  }
  return true;
}
}  // namespace detail

inline bool ParseVariant(std::string_view s, Variant* out) noexcept {
  if (!out) return false;
  if (detail::EqualsIgnoreCase(s, "sampling")) {
    *out = Variant::Sampling;
    return true;
  }
  if (detail::EqualsIgnoreCase(s, "enum_sampling") || detail::EqualsIgnoreCase(s, "enumerate_sampling") ||
      detail::EqualsIgnoreCase(s, "enumerate+sampling") || detail::EqualsIgnoreCase(s, "enum+sampling")) {
    *out = Variant::EnumSampling;
    return true;
  }
  if (detail::EqualsIgnoreCase(s, "adaptive")) {
    *out = Variant::Adaptive;
    return true;
  }
  return false;
}

inline bool ParseMethod(std::string_view s, Method* out) noexcept {
  if (!out) return false;

  // Accept a few aliases to make CLI/config friendlier.
  if (detail::EqualsIgnoreCase(s, "ours")) { *out = Method::Ours; return true; }
  if (detail::EqualsIgnoreCase(s, "aabb") || detail::EqualsIgnoreCase(s, "aabb_tree")) {
    *out = Method::AABB; return true;
  }
  if (detail::EqualsIgnoreCase(s, "interval_tree") || detail::EqualsIgnoreCase(s, "it")) {
    *out = Method::IntervalTree; return true;
  }
  if (detail::EqualsIgnoreCase(s, "kd_tree") || detail::EqualsIgnoreCase(s, "kdtree") || detail::EqualsIgnoreCase(s, "kd")) {
    *out = Method::KDTree; return true;
  }
  if (detail::EqualsIgnoreCase(s, "r_tree") || detail::EqualsIgnoreCase(s, "rtree") || detail::EqualsIgnoreCase(s, "r")) {
    *out = Method::RTree; return true;
  }
  if (detail::EqualsIgnoreCase(s, "range_tree") || detail::EqualsIgnoreCase(s, "rangetree")) {
    *out = Method::RangeTree; return true;
  }
  if (detail::EqualsIgnoreCase(s, "pbsm")) { *out = Method::PBSM; return true; }
  if (detail::EqualsIgnoreCase(s, "tlsop") || detail::EqualsIgnoreCase(s, "two_layer_sop") ||
      detail::EqualsIgnoreCase(s, "two-layer-sop")) {
    *out = Method::TLSOP; return true;
  }
  if (detail::EqualsIgnoreCase(s, "sirs")) { *out = Method::SIRS; return true; }
  if (detail::EqualsIgnoreCase(s, "rejection") || detail::EqualsIgnoreCase(s, "kds_reject") ||
      detail::EqualsIgnoreCase(s, "kds-reject")) {
    *out = Method::Rejection; return true;
  }
  if (detail::EqualsIgnoreCase(s, "tsunami")) { *out = Method::Tsunami; return true; }
  *out = Method::Unknown;
  return false;
}

inline std::ostream& operator<<(std::ostream& os, Variant v) {
  os << ToString(v);
  return os;
}
inline std::ostream& operator<<(std::ostream& os, Method m) {
  os << ToString(m);
  return os;
}

// --------------------------
// Lightweight Span (C++17) to avoid copying vectors.
// Use sjs::Span<T> similarly to std::span<T> (C++20).
// --------------------------
template <class T>
class Span {
 public:
  using element_type = T;
  using value_type = std::remove_cv_t<T>;
  using pointer = T*;
  using reference = T&;
  using iterator = pointer;
  using const_iterator = const T*;

  constexpr Span() noexcept : data_(nullptr), size_(0) {}
  constexpr Span(pointer ptr, usize n) noexcept : data_(ptr), size_(n) {}

  template <class Alloc>
  /*implicit*/ Span(std::vector<value_type, Alloc>& v) noexcept : data_(v.data()), size_(v.size()) {}

  template <class Alloc>
  /*implicit*/ Span(const std::vector<value_type, Alloc>& v) noexcept : data_(v.data()), size_(v.size()) {}

  template <usize N>
  /*implicit*/ Span(std::array<value_type, N>& a) noexcept : data_(a.data()), size_(N) {}

  template <usize N>
  /*implicit*/ Span(const std::array<value_type, N>& a) noexcept : data_(a.data()), size_(N) {}

  constexpr pointer data() const noexcept { return data_; }
  constexpr usize size() const noexcept { return size_; }
  constexpr bool empty() const noexcept { return size_ == 0; }

  constexpr reference operator[](usize i) const noexcept { return data_[i]; }
  constexpr iterator begin() const noexcept { return data_; }
  constexpr iterator end() const noexcept { return data_ + size_; }

  constexpr Span subspan(usize offset, usize count) const noexcept {
    return (offset <= size_) ? Span(data_ + offset, (count <= (size_ - offset) ? count : (size_ - offset))) : Span();
  }

 private:
  pointer data_;
  usize size_;
};

}  // namespace sjs
