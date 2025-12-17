// src/data/synthetic/generator_factory_2d.cpp
//
// Dim=2 synthetic generator registry and factory.
//
// The generator implementations are templates defined in headers. For apps that
// operate only on 2D datasets, it is convenient to select a generator by string
// name at runtime without threading template parameters everywhere.
//
// This TU provides:
//  - a registry of available generators
//  - CreateSyntheticGenerator2D(name)
//  - GenerateSyntheticDataset2D(spec)
//
// NOTE
// ----
// This file is intentionally self-contained and does not depend on baseline
// code. It is primarily used by io/load_dataset_2d.cpp and by small tooling.

#include "sjs/core/logging.h"
#include "sjs/core/types.h"

#include "sjs/data/synthetic/generator.h"

#include <algorithm>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

namespace sjs {
namespace synthetic {

namespace {

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

inline bool EqualsIgnoreCase(std::string_view a, std::string_view b) {
  if (a.size() != b.size()) return false;
  for (usize i = 0; i < a.size(); ++i) {
    char ca = a[i];
    char cb = b[i];
    if (ca >= 'A' && ca <= 'Z') ca = static_cast<char>(ca - 'A' + 'a');
    if (cb >= 'A' && cb <= 'Z') cb = static_cast<char>(cb - 'A' + 'a');
    if (ca != cb) return false;
  }
  return true;
}

}  // namespace

// A lightweight registry row.
struct GeneratorSpec2D {
  std::string_view name;  // canonical name, e.g., "stripe_ctrl_alpha"
  std::string_view desc;  // short description
};

Span<const GeneratorSpec2D> GeneratorRegistry2D() noexcept {
  static const GeneratorSpec2D kReg[] = {
      {"stripe_ctrl_alpha", "Stripe-controlled alpha generator (exact k via degree sequence)"},
      {"uniform", "Uniform random rectangles"},
      {"clustered", "Clustered/hotspot rectangles"},
      {"hetero_sizes", "Heterogeneous sizes (big+small mixture)"},
  };
  return Span<const GeneratorSpec2D>(kReg, sizeof(kReg) / sizeof(kReg[0]));
}

bool IsGeneratorSupported2D(std::string_view name) noexcept {
  const auto reg = GeneratorRegistry2D();
  for (usize i = 0; i < reg.size(); ++i) {
    if (EqualsIgnoreCase(name, reg[i].name)) return true;
  }
  return false;
}

std::string GeneratorHelp2D() {
  std::ostringstream oss;
  oss << "  Synthetic generators (Dim=2):\n";
  const auto reg = GeneratorRegistry2D();
  for (usize i = 0; i < reg.size(); ++i) {
    oss << "    - " << reg[i].name;
    if (!reg[i].desc.empty()) oss << ": " << reg[i].desc;
    oss << "\n";
  }
  return oss.str();
}

std::unique_ptr<ISyntheticGenerator<2, Scalar>> CreateSyntheticGenerator2D(std::string_view name,
                                                                          std::string* err) {
  // Delegate to the header-only factory (Create<Dim,T>) to keep the mapping
  // consistent across the codebase.
  auto gen = Create<2, Scalar>(name);
  if (!gen) {
    std::ostringstream oss;
    oss << "Unknown synthetic generator: '" << name << "'.\n";
    oss << GeneratorHelp2D();
    SetErr(err, oss.str());
    return nullptr;
  }
  return gen;
}

bool GenerateSyntheticDataset2D(std::string_view generator_name,
                               const DatasetSpec& spec,
                               Dataset<2, Scalar>* out_ds,
                               Report* report,
                               std::string* err) {
  if (!out_ds) {
    SetErr(err, "GenerateSyntheticDataset2D: out_ds is null");
    return false;
  }

  auto gen = CreateSyntheticGenerator2D(generator_name, err);
  if (!gen) return false;

  const bool ok = gen->Generate(spec, out_ds, report, err);
  if (!ok) return false;

  // Defensive: ensure the dataset is labeled half-open (the project-wide convention).
  out_ds->half_open = true;

  return true;
}

}  // namespace synthetic
}  // namespace sjs
