// src/io/load_dataset_2d.cpp
//
// Dim=2 unified dataset loading / generation entry point.
//
// This translation unit provides a small runtime wrapper that reads the project
// Config and produces a concrete Dataset<2> instance.
//
// Supported sources (cfg.dataset.source):
//  - synthetic : calls a configured synthetic generator
//  - binary    : loads relation files in the project binary format
//  - csv       : loads simple CSV/TSV rectangle files (debug/small datasets)
//
// Notes
// -----
// - The core Dataset/Relation types live in include/sjs/io/dataset.h.
// - The actual loaders are header-only (binary_io.h/csv_io.h/generator.h).
// - This file is intentionally minimal; it mainly wires Config -> loader.

#include "sjs/core/config.h"
#include "sjs/core/logging.h"
#include "sjs/core/types.h"

#include "sjs/io/binary_io.h"
#include "sjs/io/csv_io.h"
#include "sjs/io/dataset.h"
#include "sjs/io/realdata_stub.h"

#include "sjs/data/synthetic/generator.h"

// Prefer the non-templated Dim=2 wrapper implemented in
// src/data/synthetic/generator_factory_2d.cpp. We forward-declare it here
// to avoid introducing an extra header under include/.
// (Apps can also call synthetic::GenerateDataset<2> directly if desired.)
namespace sjs {
namespace synthetic {
bool GenerateSyntheticDataset2D(std::string_view generator_name,
                               const DatasetSpec& spec,
                               Dataset<2, Scalar>* out_ds,
                               Report* report,
                               std::string* err);
}  // namespace synthetic
}  // namespace sjs

#include <cctype>
#include <cerrno>
#include <cmath>
#include <filesystem>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

namespace sjs {
namespace io {

namespace {

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

inline std::string_view Trim(std::string_view s) {
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.remove_suffix(1);
  return s;
}

inline bool TryParseDouble(std::string_view s, double* out) {
  if (!out) return false;
  s = Trim(s);
  if (s.empty()) return false;
  try {
    std::size_t idx = 0;
    const double v = std::stod(std::string(s), &idx);
    if (idx != s.size()) return false;
    if (!std::isfinite(v)) return false;
    *out = v;
    return true;
  } catch (...) {
    return false;
  }
}

inline char GuessCsvSepFromPath(const std::string& path) {
  // If user provides .tsv, default to tab; else comma.
  try {
    const std::string ext = std::filesystem::path(path).extension().string();
    if (sjs::detail::EqualsIgnoreCase(ext, ".tsv") || sjs::detail::EqualsIgnoreCase(ext, ".tab")) {
      return '\t';
    }
  } catch (...) {
    // ignore
  }
  return ',';
}

inline bool FillSyntheticSpecFromConfig(const Config& cfg, synthetic::DatasetSpec* spec, std::string* err) {
  if (!spec) {
    SetErr(err, "FillSyntheticSpecFromConfig: spec is null");
    return false;
  }

  synthetic::DatasetSpec s;
  s.name = cfg.dataset.name.empty() ? std::string("synthetic") : cfg.dataset.name;
  s.n_r = cfg.dataset.synthetic.n_r;
  s.n_s = cfg.dataset.synthetic.n_s;
  s.alpha = cfg.dataset.synthetic.alpha;
  s.seed = cfg.dataset.synthetic.seed;

  // Copy-through generator parameters.
  s.params = cfg.dataset.synthetic.extra;

  // Optional: allow domain bounds via params.
  {
    auto it = s.params.find("domain_lo");
    if (it != s.params.end()) {
      double v = 0.0;
      if (!TryParseDouble(it->second, &v)) {
        SetErr(err, "synthetic.domain_lo parse error: " + it->second);
        return false;
      }
      s.domain_lo = v;
    }
  }
  {
    auto it = s.params.find("domain_hi");
    if (it != s.params.end()) {
      double v = 0.0;
      if (!TryParseDouble(it->second, &v)) {
        SetErr(err, "synthetic.domain_hi parse error: " + it->second);
        return false;
      }
      s.domain_hi = v;
    }
  }

  if (!(s.domain_hi > s.domain_lo)) {
    SetErr(err, "synthetic domain must satisfy domain_hi > domain_lo");
    return false;
  }

  *spec = std::move(s);
  return true;
}

}  // namespace

// --------------------------
// Public API (2D)
// --------------------------

// Load/generate a 2D dataset according to cfg.
//
// - out_ds must be non-null.
// - If out_syn_report is non-null and source is synthetic, it is filled.
bool LoadDataset2D(const Config& cfg,
                   Dataset<2, Scalar>* out_ds,
                   synthetic::Report* out_syn_report,
                   std::string* err) {
  if (!out_ds) {
    SetErr(err, "LoadDataset2D: out_ds is null");
    return false;
  }
  if (cfg.dataset.dim != 2) {
    std::ostringstream oss;
    oss << "LoadDataset2D: cfg.dataset.dim=" << cfg.dataset.dim << " (expected 2)";
    SetErr(err, oss.str());
    return false;
  }

  std::string cfg_err;
  if (!cfg.Validate(&cfg_err)) {
    SetErr(err, "Config validation failed: " + cfg_err);
    return false;
  }

  Dataset<2, Scalar> ds;
  ds.half_open = true;

  switch (cfg.dataset.source) {
    case DataSource::Synthetic: {
      synthetic::DatasetSpec spec;
      if (!FillSyntheticSpecFromConfig(cfg, &spec, err)) return false;

      synthetic::Report rep;
      const std::string gen_name = cfg.dataset.synthetic.generator;

      SJS_LOG_INFO("Generating synthetic dataset (2D): gen=", gen_name,
                   " n_r=", spec.n_r, " n_s=", spec.n_s,
                   " alpha=", spec.alpha, " seed=", spec.seed);

      std::string gen_err;
      if (!synthetic::GenerateSyntheticDataset2D(gen_name, spec, &ds,
                                                 out_syn_report ? &rep : nullptr,
                                                 &gen_err)) {
        SetErr(err, "Synthetic generation failed: " + gen_err);
        return false;
      }

      // Fill output report if requested.
      if (out_syn_report) {
        *out_syn_report = std::move(rep);
      }

      // If generator didn't set a dataset name, use cfg.dataset.name.
      if (ds.name.empty()) ds.name = spec.name;
      break;
    }

    case DataSource::Binary: {
      binary::BinaryReadOptions opt;
      opt.generate_ids_if_missing = true;
      opt.drop_empty = false;

      std::string io_err;
      if (!binary::ReadDatasetBinaryPair<2, Scalar>(cfg.dataset.path_r,
                                                    cfg.dataset.path_s,
                                                    &ds,
                                                    opt,
                                                    &io_err)) {
        SetErr(err, "Binary dataset load failed: " + io_err);
        return false;
      }
      ds.name = cfg.dataset.name.empty() ? std::string("binary") : cfg.dataset.name;
      break;
    }

    case DataSource::CSV: {
      const char sep_r = GuessCsvSepFromPath(cfg.dataset.path_r);
      const char sep_s = GuessCsvSepFromPath(cfg.dataset.path_s);

      Relation<2, Scalar> R;
      Relation<2, Scalar> S;
      R.name = "R";
      S.name = "S";

      std::string io_err;
      if (!csv::ReadBoxesSimple<2, Scalar>(cfg.dataset.path_r, &R, sep_r,
                                           /*has_header=*/true,
                                           &io_err)) {
        SetErr(err, "CSV load failed for R: " + io_err);
        return false;
      }
      io_err.clear();
      if (!csv::ReadBoxesSimple<2, Scalar>(cfg.dataset.path_s, &S, sep_s,
                                           /*has_header=*/true,
                                           &io_err)) {
        SetErr(err, "CSV load failed for S: " + io_err);
        return false;
      }

      ds.name = cfg.dataset.name.empty() ? std::string("csv") : cfg.dataset.name;
      ds.R = std::move(R);
      ds.S = std::move(S);
      break;
    }

    case DataSource::Unknown:
    default: {
      SetErr(err,
             "LoadDataset2D: cfg.dataset.source is unknown. "
             "Set --dataset_source=synthetic|binary|csv.");
      return false;
    }
  }

  // Post-processing common to all sources.
  ds.half_open = true;
  ds.EnsureIds();

  // Optionally drop empties if requested via cfg.dataset.synthetic.extra etc.
  // (We keep this conservative here; most generators already ensure proper boxes.)
  // Users can always call Dataset::RemoveEmptyBoxes() in their app if desired.

  {
    std::string v_err;
    if (!ds.Validate(/*require_proper=*/true, &v_err)) {
      SetErr(err, "Loaded dataset is invalid: " + v_err);
      return false;
    }
  }

  SJS_LOG_INFO("Dataset ready: name=", ds.name,
               " |R|=", ds.R.Size(), " |S|=", ds.S.Size());

  *out_ds = std::move(ds);
  return true;
}

}  // namespace io
}  // namespace sjs
