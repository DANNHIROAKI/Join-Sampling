#pragma once
// sjs/data/synthetic/generator.h
//
// Synthetic dataset generation interfaces.
//
// Goals:
//  - C++17, header-only, dependency-light.
//  - Current focus: 2D experiments, but APIs are templated by Dim for easy extension.
//  - Pluggable generators via a small factory.
//
// A generator produces a Dataset<Dim> with two relations (R,S) of axis-aligned boxes.
// All boxes follow the project's half-open convention: [lo, hi) on every dimension.
//
// This module intentionally does NOT compute the true join size |J| for arbitrary
// generators (that could be expensive). For the stripe-controlled generator we
// can provide exact k=|J| by construction; other generators may leave it unset.

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/io/dataset.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
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
// Spec + report
// --------------------------

struct DatasetSpec {
  // Human-friendly tag used in logs / output file names.
  std::string name = "synthetic";

  // Sizes for R and S.
  u64 n_r = 100000;
  u64 n_s = 100000;

  // A high-level knob; exact semantics depend on the generator.
  // For stripe_ctrl_alpha, this is alpha where k = round(alpha * (n_r + n_s)).
  double alpha = 1e-6;

  // RNG seed for generation (separate from sampling seed).
  u64 seed = 1;

  // Domain bounds per axis. Most generators assume a hyper-rectangle domain:
  //   [domain_lo, domain_hi) ^ Dim
  // Default matches the paper note: [0,1)^2.
  double domain_lo = 0.0;
  double domain_hi = 1.0;

  // Generator-specific parameters (string -> string).
  // Example: {"core_lo":"0.45", "core_hi":"0.55", "gap_factor":"0.1"}.
  std::unordered_map<std::string, std::string> params;
};

// Optional metadata returned by generators.
struct Report {
  std::string generator;
  std::string dataset_name;

  u64 n_r = 0;
  u64 n_s = 0;

  // If known exactly (e.g., stripe generator), these are populated.
  // Otherwise they remain 0 / NaN.
  bool has_exact_k = false;
  u64 k_target = 0;
  u64 k_achieved = 0;

  double alpha_target = std::numeric_limits<double>::quiet_NaN();
  double alpha_achieved = std::numeric_limits<double>::quiet_NaN();

  // Free-form notes for debugging/diagnostics.
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

inline std::optional<std::string_view> FindParam(const std::unordered_map<std::string, std::string>& m,
                                                 std::string_view key) {
  auto it = m.find(std::string(key));
  if (it == m.end()) return std::nullopt;
  return std::string_view(it->second);
}

inline bool EqualsIgnoreCase(std::string_view a, std::string_view b) noexcept {
  if (a.size() != b.size()) return false;
  for (usize i = 0; i < a.size(); ++i) {
    const char ca = (a[i] >= 'A' && a[i] <= 'Z') ? static_cast<char>(a[i] - 'A' + 'a') : a[i];
    const char cb = (b[i] >= 'A' && b[i] <= 'Z') ? static_cast<char>(b[i] - 'A' + 'a') : b[i];
    if (ca != cb) return false;
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

// Get a typed param or default.
inline double GetDouble(const std::unordered_map<std::string, std::string>& m,
                        std::string_view key,
                        double def) {
  if (auto v = FindParam(m, key)) {
    double x;
    if (TryParseDouble(*v, &x)) return x;
  }
  return def;
}

inline u64 GetU64(const std::unordered_map<std::string, std::string>& m,
                  std::string_view key,
                  u64 def) {
  if (auto v = FindParam(m, key)) {
    u64 x;
    if (TryParseU64(*v, &x)) return x;
  }
  return def;
}

inline i32 GetI32(const std::unordered_map<std::string, std::string>& m,
                  std::string_view key,
                  i32 def) {
  if (auto v = FindParam(m, key)) {
    i32 x;
    if (TryParseI32(*v, &x)) return x;
  }
  return def;
}

inline bool GetBool(const std::unordered_map<std::string, std::string>& m,
                    std::string_view key,
                    bool def) {
  if (auto v = FindParam(m, key)) {
    bool x;
    if (TryParseBool(*v, &x)) return x;
  }
  return def;
}

// Side-specific parameters: first try "<prefix><key>", then "<key>".
inline double GetDoubleSide(const std::unordered_map<std::string, std::string>& m,
                            std::string_view prefix,
                            std::string_view key,
                            double def) {
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

inline bool GetBoolSide(const std::unordered_map<std::string, std::string>& m,
                        std::string_view prefix,
                        std::string_view key,
                        bool def) {
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

inline i32 GetI32Side(const std::unordered_map<std::string, std::string>& m,
                      std::string_view prefix,
                      std::string_view key,
                      i32 def) {
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

// Fisher-Yates shuffle using sjs::Rng.
template <class T>
inline void ShuffleInPlace(std::vector<T>* v, Rng* rng) {
  if (!v || v->size() <= 1) return;
  SJS_ASSERT(rng != nullptr);
  for (usize i = v->size() - 1; i > 0; --i) {
    const usize j = static_cast<usize>(rng->UniformU64(static_cast<u64>(i + 1)));
    std::swap((*v)[i], (*v)[j]);
  }
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

  // Generate into out_ds. Returns true on success.
  // If report != nullptr, fills metadata.
  virtual bool Generate(const DatasetSpec& spec,
                        DatasetT* out_ds,
                        Report* report,
                        std::string* err) = 0;
};

// Factory (defined in generator.h but implemented by including the known generators).
template <int Dim, class T = Scalar>
std::unique_ptr<ISyntheticGenerator<Dim, T>> Create(std::string_view generator_name);

// Convenience one-shot API.
template <int Dim, class T = Scalar>
inline bool GenerateDataset(std::string_view generator_name,
                            const DatasetSpec& spec,
                            Dataset<Dim, T>* out_ds,
                            Report* report = nullptr,
                            std::string* err = nullptr) {
  auto gen = Create<Dim, T>(generator_name);
  if (!gen) {
    if (err) *err = "Unknown generator: " + std::string(generator_name);
    return false;
  }
  return gen->Generate(spec, out_ds, report, err);
}

}  // namespace synthetic
}  // namespace sjs


// --------------------------
// Built-in generator registry
// --------------------------
//
// To keep the experiment code simple, we provide a header-only factory that
// knows about the built-in generators listed in this folder.
//
// If you want a leaner build, move this Create() implementation into a .cc file
// and only include the generators you need.

#include "sjs/data/synthetic/stripe_ctrl_alpha.h"
#include "sjs/data/synthetic/uniform.h"
#include "sjs/data/synthetic/clustered.h"
#include "sjs/data/synthetic/hetero_sizes.h"

namespace sjs {
namespace synthetic {

template <int Dim, class T>
inline std::unique_ptr<ISyntheticGenerator<Dim, T>> Create(std::string_view generator_name) {
  // Accept some short aliases for convenience.
  if (detail::EqualsIgnoreCase(generator_name, "stripe_ctrl_alpha") ||
      detail::EqualsIgnoreCase(generator_name, "stripe") ||
      detail::EqualsIgnoreCase(generator_name, "scirg")) {
    return std::make_unique<StripeCtrlAlphaGenerator<Dim, T>>();
  }
  if (detail::EqualsIgnoreCase(generator_name, "uniform") ||
      detail::EqualsIgnoreCase(generator_name, "uni")) {
    return std::make_unique<UniformGenerator<Dim, T>>();
  }
  if (detail::EqualsIgnoreCase(generator_name, "clustered") ||
      detail::EqualsIgnoreCase(generator_name, "hotspot") ||
      detail::EqualsIgnoreCase(generator_name, "clusters")) {
    return std::make_unique<ClusteredGenerator<Dim, T>>();
  }
  if (detail::EqualsIgnoreCase(generator_name, "hetero_sizes") ||
      detail::EqualsIgnoreCase(generator_name, "hetero") ||
      detail::EqualsIgnoreCase(generator_name, "mix_sizes")) {
    return std::make_unique<HeteroSizesGenerator<Dim, T>>();
  }
  return nullptr;
}

}  // namespace synthetic
}  // namespace sjs
