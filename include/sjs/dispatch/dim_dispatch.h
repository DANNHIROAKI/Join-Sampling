#pragma once
// sjs/dispatch/dim_dispatch.h
//
// Runtime-dim -> compile-time Dim dispatch.
//
// Usage pattern:
//   struct Fn {
//     template<int Dim> R operator()() const { ... }
//   };
//   R out;
//   std::string err;
//   bool ok = sjs::dispatch::DispatchDim<R>(dim, Fn{...}, &out, &err);
//
// Notes:
//  - Supported dims are controlled by the X-macro SJS_DIMS(X).
//  - You can override SJS_DIMS(X) at compile time, e.g. via a generated header
//    or compiler definitions, but the simplest is to edit the default list below.
//
// Default supported dims: 2..12
// (You can widen this list if you need higher dimensions.)

#include "sjs/core/types.h"

#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#ifndef SJS_MAX_DIM
// Default maximum dimension compiled into the library.
//
// You may reduce this value (e.g. -DSJS_MAX_DIM=6) to speed up compilation,
// or increase it if you need higher dimensions.
#define SJS_MAX_DIM 12
#endif

#if SJS_MAX_DIM < 2
#error "SJS_MAX_DIM must be >= 2"
#endif

#ifndef SJS_DIMS
// X-macro list of supported dimensions (2..SJS_MAX_DIM).
// You can override SJS_DIMS(X) directly if you want a non-contiguous set.
//
// Implementation detail: we use a small macro chain to conditionally include
// X(d) for each d based on SJS_MAX_DIM.

#if SJS_MAX_DIM >= 3
#define SJS_DIMS_3(X) X(3) SJS_DIMS_4(X)
#else
#define SJS_DIMS_3(X)
#endif

#if SJS_MAX_DIM >= 4
#define SJS_DIMS_4(X) X(4) SJS_DIMS_5(X)
#else
#define SJS_DIMS_4(X)
#endif

#if SJS_MAX_DIM >= 5
#define SJS_DIMS_5(X) X(5) SJS_DIMS_6(X)
#else
#define SJS_DIMS_5(X)
#endif

#if SJS_MAX_DIM >= 6
#define SJS_DIMS_6(X) X(6) SJS_DIMS_7(X)
#else
#define SJS_DIMS_6(X)
#endif

#if SJS_MAX_DIM >= 7
#define SJS_DIMS_7(X) X(7) SJS_DIMS_8(X)
#else
#define SJS_DIMS_7(X)
#endif

#if SJS_MAX_DIM >= 8
#define SJS_DIMS_8(X) X(8) SJS_DIMS_9(X)
#else
#define SJS_DIMS_8(X)
#endif

#if SJS_MAX_DIM >= 9
#define SJS_DIMS_9(X) X(9) SJS_DIMS_10(X)
#else
#define SJS_DIMS_9(X)
#endif

#if SJS_MAX_DIM >= 10
#define SJS_DIMS_10(X) X(10) SJS_DIMS_11(X)
#else
#define SJS_DIMS_10(X)
#endif

#if SJS_MAX_DIM >= 11
#define SJS_DIMS_11(X) X(11) SJS_DIMS_12(X)
#else
#define SJS_DIMS_11(X)
#endif

#if SJS_MAX_DIM >= 12
#define SJS_DIMS_12(X) X(12)
#else
#define SJS_DIMS_12(X)
#endif

#define SJS_DIMS(X) X(2) SJS_DIMS_3(X)

#endif  // SJS_DIMS

namespace sjs::dispatch {

inline std::vector<int> SupportedDims() {
  std::vector<int> v;
  v.reserve(16);
#define SJS_PUSH_DIM(D) v.push_back(D);
  SJS_DIMS(SJS_PUSH_DIM)
#undef SJS_PUSH_DIM
  return v;
}

inline std::string SupportedDimsString() {
  std::ostringstream oss;
  bool first = true;
#define SJS_APPEND_DIM(D)       \
  do {                          \
    if (!first) oss << ",";     \
    first = false;              \
    oss << (D);                 \
  } while (0);
  SJS_DIMS(SJS_APPEND_DIM)
#undef SJS_APPEND_DIM
  return oss.str();
}

inline bool IsSupportedDim(int dim) {
  switch (dim) {
#define SJS_CASE(D) case (D): return true;
    SJS_DIMS(SJS_CASE)
#undef SJS_CASE
    default: return false;
  }
}

// Dispatch for non-void return type R.
template <class R, class Fn>
inline std::enable_if_t<!std::is_void_v<R>, bool>
DispatchDim(int dim, Fn&& fn, R* out, std::string* err = nullptr) {
  switch (dim) {
#define SJS_CASE(D)                                                                 \
  case (D): {                                                                       \
    if (out) {                                                                      \
      *out = std::forward<Fn>(fn).template operator()<D>();                         \
    } else {                                                                        \
      (void)std::forward<Fn>(fn).template operator()<D>();                          \
    }                                                                               \
    return true;                                                                    \
  }
    SJS_DIMS(SJS_CASE)
#undef SJS_CASE
    default: {
      if (err) {
        *err = "Unsupported dim=" + std::to_string(dim) +
               ". Supported dims: [" + SupportedDimsString() + "].";
      }
      return false;
    }
  }
}

// Dispatch for void-return functor.
template <class Fn>
inline bool DispatchDimVoid(int dim, Fn&& fn, std::string* err = nullptr) {
  switch (dim) {
#define SJS_CASE(D)                                                                 \
  case (D): {                                                                       \
    std::forward<Fn>(fn).template operator()<D>();                                  \
    return true;                                                                    \
  }
    SJS_DIMS(SJS_CASE)
#undef SJS_CASE
    default: {
      if (err) {
        *err = "Unsupported dim=" + std::to_string(dim) +
               ". Supported dims: [" + SupportedDimsString() + "].";
      }
      return false;
    }
  }
}

}  // namespace sjs::dispatch