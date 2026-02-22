#pragma once
// sjs/data/synthetic/generator.h
//
// Synthetic dataset generation interfaces + shared utilities.
//
// A generator produces Dataset<Dim> with two relations (R,S) of half-open boxes:
//   r = Π_i [L_i(r), R_i(r)) and always L_i < R_i.
//
// Design goals:
//  - Header-only generators, C++17, dependency-light.
//  - Dimension-agnostic: all generators are templated by Dim.
//  - Generator-specific params passed as string->string map (easy to wire from JSON/CLI).
//
// NOTE: The "factory/registry" lives in generator_factory.h (not here).

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/io/dataset.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sjs {
namespace synthetic {

// --------------------------
// Spec + Report
// --------------------------

using ParamMap = std::unordered_map<std::string, std::string>;

struct DatasetSpec {
  // Dataset tag used in output/logging.
  std::string name = "synthetic";

  // Sizes.
  u64 n_r = 100000;
  u64 n_s = 100000;

  // High-level knob; semantics depend on generator.
  // For stripe_ctrl_alpha: k = round(alpha * (n_r + n_s)) unless overridden by params["k_target"].
  double alpha = 1.0;

  // RNG seed for dataset generation (separate from algorithm sampling seed).
  u64 seed = 1;

  // Shared domain per axis: [domain_lo, domain_hi)^Dim.
  double domain_lo = 0.0;
  double domain_hi = 1.0;

  // Generator-specific parameters.
  ParamMap params;
};

struct Report {
  std::string generator;
  std::string dataset_name;

  u64 n_r = 0;
  u64 n_s = 0;

  // If generator can provide exact |J| by construction (stripe_ctrl_alpha),
  // these are populated.
  bool has_exact_k = false;
  u64 k_target = 0;
  u64 k_achieved = 0;

  double alpha_target = std::numeric_limits<double>::quiet_NaN();
  double alpha_achieved = std::numeric_limits<double>::quiet_NaN();

  std::string notes;

  std::string ToJsonLite() const {
    std::ostringstream oss;
    oss << "{"
        << "\"generator\":\"" << generator << "\","
        << "\"dataset\":\"" << dataset_name << "\","
        << "\"n_r\":" << n_r << ","
        << "\"n_s\":" << n_s << ","
        << "\"has_exact_k\":" << (has_exact_k ? "true" : "false") << ","
        << "\"k_target\":" << k_target << ","
        << "\"k_achieved\":" << k_achieved << ","
        << "\"alpha_target\":" << alpha_target << ","
        << "\"alpha_achieved\":" << alpha_achieved << ","
        << "\"notes\":\"" << notes << "\""
        << "}";
    return oss.str();
  }
};

// --------------------------
// Parameter parsing helpers
// --------------------------
namespace detail {

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

inline std::optional<std::string_view> FindParam(const ParamMap& m, std::string_view key) {
  auto it = m.find(std::string(key));
  if (it == m.end()) return std::nullopt;
  return std::string_view(it->second);
}

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

inline bool TryParseBool(std::string_view s, bool* out) {
  if (!out) return false;
  if (EqualsIgnoreCase(s, "1") || EqualsIgnoreCase(s, "true") || EqualsIgnoreCase(s, "yes") ||
      EqualsIgnoreCase(s, "y") || EqualsIgnoreCase(s, "on")) {
    *out = true;
    return true;
  }
  if (EqualsIgnoreCase(s, "0") || EqualsIgnoreCase(s, "false") || EqualsIgnoreCase(s, "no") ||
      EqualsIgnoreCase(s, "n") || EqualsIgnoreCase(s, "off")) {
    *out = false;
    return true;
  }
  return false;
}

inline bool TryParseI32(std::string_view s, i32* out) {
  if (!out) return false;
  try {
    std::size_t idx = 0;
    long v = std::stol(std::string(s), &idx, 10);
    if (idx != s.size()) return false;
    *out = static_cast<i32>(v);
    return true;
  } catch (...) {
    return false;
  }
}

inline bool TryParseU64(std::string_view s, u64* out) {
  if (!out) return false;
  try {
    std::size_t idx = 0;
    unsigned long long v = std::stoull(std::string(s), &idx, 10);
    if (idx != s.size()) return false;
    *out = static_cast<u64>(v);
    return true;
  } catch (...) {
    return false;
  }
}

inline bool TryParseDouble(std::string_view s, double* out) {
  if (!out) return false;
  try {
    std::size_t idx = 0;
    double v = std::stod(std::string(s), &idx);
    if (idx != s.size()) return false;
    *out = v;
    return true;
  } catch (...) {
    return false;
  }
}

inline double GetDouble(const ParamMap& m, std::string_view key, double def) {
  if (auto v = FindParam(m, key)) {
    double x;
    if (TryParseDouble(*v, &x)) return x;
  }
  return def;
}

inline u64 GetU64(const ParamMap& m, std::string_view key, u64 def) {
  if (auto v = FindParam(m, key)) {
    u64 x;
    if (TryParseU64(*v, &x)) return x;
  }
  return def;
}

inline i32 GetI32(const ParamMap& m, std::string_view key, i32 def) {
  if (auto v = FindParam(m, key)) {
    i32 x;
    if (TryParseI32(*v, &x)) return x;
  }
  return def;
}

inline bool GetBool(const ParamMap& m, std::string_view key, bool def) {
  if (auto v = FindParam(m, key)) {
    bool x;
    if (TryParseBool(*v, &x)) return x;
  }
  return def;
}

// Side-specific: first try "<prefix><key>", then "<key>".
inline double GetDoubleSide(const ParamMap& m, std::string_view prefix, std::string_view key, double def) {
  std::string pk;
  pk.reserve(prefix.size() + key.size());
  pk.append(prefix);
  pk.append(key);
  if (auto v = FindParam(m, pk)) {
    double x;
    if (TryParseDouble(*v, &x)) return x;
  }
  return GetDouble(m, key, def);
}

inline i32 GetI32Side(const ParamMap& m, std::string_view prefix, std::string_view key, i32 def) {
  std::string pk;
  pk.reserve(prefix.size() + key.size());
  pk.append(prefix);
  pk.append(key);
  if (auto v = FindParam(m, pk)) {
    i32 x;
    if (TryParseI32(*v, &x)) return x;
  }
  return GetI32(m, key, def);
}

inline bool GetBoolSide(const ParamMap& m, std::string_view prefix, std::string_view key, bool def) {
  std::string pk;
  pk.reserve(prefix.size() + key.size());
  pk.append(prefix);
  pk.append(key);
  if (auto v = FindParam(m, pk)) {
    bool x;
    if (TryParseBool(*v, &x)) return x;
  }
  return GetBool(m, key, def);
}

// In-place Fisher-Yates shuffle for relations (boxes + ids together).
template <int Dim, class T>
inline void ShuffleRelation(Relation<Dim, T>* rel, Rng* rng) {
  if (!rel) return;
  SJS_ASSERT(rng != nullptr);

  const usize n = rel->boxes.size();
  if (n <= 1) return;

  const bool has_ids = !rel->ids.empty();
  if (has_ids) {
    SJS_ASSERT(rel->ids.size() == n);
  }

  for (usize i = n - 1; i > 0; --i) {
    const usize j = static_cast<usize>(rng->UniformU64(static_cast<u64>(i + 1)));
    std::swap(rel->boxes[i], rel->boxes[j]);
    if (has_ids) std::swap(rel->ids[i], rel->ids[j]);
  }
}

inline bool CheckIdFits(u64 n, std::string* err, std::string_view who) {
  if (n > static_cast<u64>(std::numeric_limits<Id>::max())) {
    if (err) {
      std::ostringstream oss;
      oss << std::string(who) << ": n=" << n << " exceeds Id max="
          << static_cast<u64>(std::numeric_limits<Id>::max());
      *err = oss.str();
    }
    return false;
  }
  return true;
}

}  // namespace detail

// --------------------------
// Generator interface
// --------------------------

template <int Dim, class T = Scalar>
class ISyntheticGenerator {
 public:
  using DatasetT = Dataset<Dim, T>;

  virtual ~ISyntheticGenerator() = default;
  virtual std::string_view Name() const noexcept = 0;

  // Generate dataset into out_ds. Returns true on success.
  // If report != nullptr, fills generation metadata.
  virtual bool Generate(const DatasetSpec& spec, DatasetT* out_ds, Report* report, std::string* err) = 0;
};

}  // namespace synthetic
}  // namespace sjs
