#pragma once
// sjs/core/assert.h
//
// Minimal, dependency-free assertions and checks.
// Use SJS_ASSERT for always-on checks, SJS_DASSERT for debug-only checks.
// Use SJS_CHECK_* macros for rich diagnostics with values.

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string_view>
#include <type_traits>

namespace sjs {
namespace detail {

[[noreturn]] inline void Fail(const char* kind,
                             std::string_view expr,
                             std::string_view msg,
                             const char* file,
                             int line,
                             const char* func) {
  std::ostringstream oss;
  oss << "[SJS][" << kind << "] " << file << ":" << line << " in " << func << "\n";
  if (!expr.empty()) oss << "  expr: " << expr << "\n";
  if (!msg.empty()) oss << "  msg : " << msg << "\n";
  std::cerr << oss.str() << std::flush;
  std::abort();
}

template <class T>
inline void PrintValue(std::ostream& os, const T& v) {
  if constexpr (std::is_same_v<T, bool>) {
    os << (v ? "true" : "false");
  } else {
    os << v;
  }
}

template <class A, class B>
[[noreturn]] inline void CheckFail(const char* a_expr,
                                  const char* b_expr,
                                  const char* op,
                                  const A& a_val,
                                  const B& b_val,
                                  const char* file,
                                  int line,
                                  const char* func) {
  std::ostringstream oss;
  oss << "[SJS][CHECK] " << file << ":" << line << " in " << func << "\n"
      << "  check: (" << a_expr << ") " << op << " (" << b_expr << ")\n"
      << "  lhs  : ";
  PrintValue(oss, a_val);
  oss << "\n  rhs  : ";
  PrintValue(oss, b_val);
  oss << "\n";
  std::cerr << oss.str() << std::flush;
  std::abort();
}

}  // namespace detail
}  // namespace sjs

// Branch prediction hints (portable enough for GCC/Clang; safe fallback elsewhere)
#if defined(__GNUC__) || defined(__clang__)
  #define SJS_LIKELY(x)   (__builtin_expect(!!(x), 1))
  #define SJS_UNLIKELY(x) (__builtin_expect(!!(x), 0))
#else
  #define SJS_LIKELY(x)   (x)
  #define SJS_UNLIKELY(x) (x)
#endif

#define SJS_ASSERT(cond)                                                                 \
  do {                                                                                   \
    if (SJS_UNLIKELY(!(cond))) {                                                         \
      ::sjs::detail::Fail("ASSERT", #cond, "", __FILE__, __LINE__, __func__);            \
    }                                                                                    \
  } while (0)

#define SJS_ASSERT_MSG(cond, msg)                                                        \
  do {                                                                                   \
    if (SJS_UNLIKELY(!(cond))) {                                                         \
      ::sjs::detail::Fail("ASSERT", #cond, (msg), __FILE__, __LINE__, __func__);         \
    }                                                                                    \
  } while (0)

#ifndef NDEBUG
  #define SJS_DASSERT(cond) SJS_ASSERT(cond)
  #define SJS_DASSERT_MSG(cond, msg) SJS_ASSERT_MSG(cond, msg)
#else
  #define SJS_DASSERT(cond)       do { (void)sizeof(cond); } while (0)
  #define SJS_DASSERT_MSG(cond, msg) do { (void)sizeof(cond); (void)sizeof(msg); } while (0)
#endif

#define SJS_UNREACHABLE()                                                                \
  do {                                                                                   \
    ::sjs::detail::Fail("UNREACHABLE", "", "", __FILE__, __LINE__, __func__);            \
  } while (0)

// Rich checks with value printing.
#define SJS_CHECK_OP(a, b, op)                                                           \
  do {                                                                                   \
    const auto& _sjs_a = (a);                                                            \
    const auto& _sjs_b = (b);                                                            \
    if (SJS_UNLIKELY(!(_sjs_a op _sjs_b))) {                                              \
      ::sjs::detail::CheckFail(#a, #b, #op, _sjs_a, _sjs_b, __FILE__, __LINE__, __func__);\
    }                                                                                    \
  } while (0)

#define SJS_CHECK_EQ(a, b) SJS_CHECK_OP(a, b, ==)
#define SJS_CHECK_NE(a, b) SJS_CHECK_OP(a, b, !=)
#define SJS_CHECK_LT(a, b) SJS_CHECK_OP(a, b, <)
#define SJS_CHECK_LE(a, b) SJS_CHECK_OP(a, b, <=)
#define SJS_CHECK_GT(a, b) SJS_CHECK_OP(a, b, >)
#define SJS_CHECK_GE(a, b) SJS_CHECK_OP(a, b, >=)
