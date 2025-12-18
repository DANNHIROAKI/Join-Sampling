// apps/sjs_gen_dataset.cpp
//
// Synthetic dataset generator + exporter.
//
// Generates a Dataset<2> (R/S) using the configured synthetic generator, then
// writes it out in:
//   - SJS binary format (fast, robust), and optionally
//   - CSV (for debugging/visualization)
//
// Examples:
//   ./sjs_gen_dataset --dataset_source=synthetic --gen=stripe --dataset=demo
//             --n_r=100000 --n_s=100000 --alpha=1e-6 --gen_seed=1
//             --out_dir=data/demo --write_csv=1
//
// Output default files:
//   <out_dir>/<dataset>_R.bin
//   <out_dir>/<dataset>_S.bin
//   (optional)
//   <out_dir>/<dataset>_R.csv
//   <out_dir>/<dataset>_S.csv

#include "sjs/core/config.h"
#include "sjs/core/logging.h"
#include "sjs/core/types.h"

#include "sjs/data/synthetic/generator.h"
#include "sjs/io/binary_io.h"
#include "sjs/io/csv_io.h"
#include "sjs/io/dataset.h"

#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>

namespace fs = std::filesystem;

namespace sjs {
namespace apps {

namespace {

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

inline bool IsHelpRequested(const sjs::ArgMap& args) {
  return args.Has("help") || args.Has("h") || args.Has("-h") || args.Has("--help");
}

inline std::string SanitizeFilename(std::string_view s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    const bool ok = (std::isalnum(static_cast<unsigned char>(c)) != 0) || c == '_' || c == '-' || c == '.';
    out.push_back(ok ? c : '_');
  }
  if (out.empty()) out = "x";
  return out;
}

inline bool EnsureDir(const fs::path& p, std::string* err) {
  try {
    if (p.empty()) return true;
    fs::create_directories(p);
    return true;
  } catch (const std::exception& e) {
    SetErr(err, std::string("create_directories failed: ") + e.what());
    return false;
  }
}

inline bool ParseBool(const std::optional<std::string_view>& v, bool def) {
  if (!v) return def;
  if (*v == "1" || *v == "true" || *v == "yes") return true;
  if (*v == "0" || *v == "false" || *v == "no") return false;
  return def;
}

inline void PrintUsage() {
  std::cerr
      << "sjs_gen_dataset: synthetic dataset generator\n\n"
      << "Required:\n"
      << "  --dataset_source=synthetic\n"
      << "  --gen=<stripe|uniform|clustered|hetero_sizes>\n"
      << "  --dataset=<name>\n"
      << "  --n_r=<N> --n_s=<N>\n"
      << "  --alpha=<float>\n"
      << "  --gen_seed=<seed>\n"
      << "\nOutput:\n"
      << "  --out_dir=<dir>\n"
      << "  --out_r=<path> --out_s=<path>            (override binary outputs)\n"
      << "  --write_csv=0|1                          (default 0)\n"
      << "  --csv_r=<path> --csv_s=<path>            (override CSV outputs)\n"
      << "  --csv_sep=,|tab|\\t                       (default ',')\n"
      << "\n";
}

}  // namespace

}  // namespace apps
}  // namespace sjs

int main(int argc, char** argv) {
  sjs::ArgMap args = sjs::ArgMap::FromArgv(argc, argv);
  if (sjs::apps::IsHelpRequested(args)) {
    sjs::apps::PrintUsage();
    return 0;
  }

  sjs::Config cfg = sjs::Config::FromArgs(argc, argv);
  std::string err;
  if (!cfg.Validate(&err)) {
    SJS_LOG_ERROR("Config validation failed:", err);
    sjs::apps::PrintUsage();
    return 2;
  }
  sjs::Logger::Instance().SetConfig(cfg.logging);

  if (cfg.dataset.dim != 2) {
    SJS_LOG_ERROR("This build supports only Dim=2 for now. Got --dim=", cfg.dataset.dim);
    return 2;
  }
  if (cfg.dataset.source != sjs::DataSource::Synthetic) {
    SJS_LOG_ERROR("sjs_gen_dataset requires --dataset_source=synthetic");
    return 2;
  }

  // Generate dataset.
  sjs::Dataset<2, sjs::Scalar> ds;
  sjs::synthetic::DatasetSpec spec;
  spec.name = cfg.dataset.name;
  spec.n_r = cfg.dataset.synthetic.n_r;
  spec.n_s = cfg.dataset.synthetic.n_s;
  spec.alpha = cfg.dataset.synthetic.alpha;
  spec.seed = cfg.dataset.synthetic.seed;
  spec.params = cfg.dataset.synthetic.extra;

  sjs::synthetic::Report rep;
  if (!sjs::synthetic::GenerateDataset<2, sjs::Scalar>(cfg.dataset.synthetic.generator, spec, &ds, &rep, &err)) {
    SJS_LOG_ERROR("Generation failed:", err);
    return 3;
  }

  SJS_LOG_INFO("Generated dataset:", ds.name,
               "R=", static_cast<unsigned long long>(ds.R.Size()),
               "S=", static_cast<unsigned long long>(ds.S.Size()),
               "report=", rep.ToJsonLite());

  // Output paths
  fs::path out_dir(cfg.output.out_dir);
  // Default output directory for this generator app is data/synthetic unless overridden.
  if (!args.Has("out_dir")) {
    out_dir = fs::path("data/synthetic");
  }

  if (!sjs::apps::EnsureDir(out_dir, &err)) {
    SJS_LOG_ERROR("Cannot create out_dir:", err);
    return 4;
  }

  const std::string base = sjs::apps::SanitizeFilename(ds.name);

  std::string out_r = (out_dir / (base + "_R.bin")).string();
  std::string out_s = (out_dir / (base + "_S.bin")).string();
  if (auto v = args.Get("out_r")) out_r = std::string(*v);
  if (auto v = args.Get("out_s")) out_s = std::string(*v);

  // Write binary.
  {
    sjs::binary::BinaryWriteOptions opt;
    opt.half_open = true;
    opt.write_ids = true;
    opt.scalar = sjs::binary::ScalarEncoding::Float64;
    opt.write_name = true;

    if (!sjs::binary::WriteDatasetBinaryPair<2, sjs::Scalar>(out_r, out_s, ds, opt, &err)) {
      SJS_LOG_ERROR("Binary write failed:", err);
      return 5;
    }
    SJS_LOG_INFO("Wrote binary:", out_r, "and", out_s);
  }

  // Optional CSV.
  const bool write_csv = sjs::apps::ParseBool(args.Get("write_csv"), /*def=*/false);
  if (write_csv) {
    char sep = ',';
    if (auto v = args.Get("csv_sep")) {
      if (*v == "tab" || *v == "\\t") sep = '\t';
      else if (!v->empty()) sep = (*v)[0];
    }

    std::string csv_r = (out_dir / (base + "_R.csv")).string();
    std::string csv_s = (out_dir / (base + "_S.csv")).string();
    if (auto v = args.Get("csv_r")) csv_r = std::string(*v);
    if (auto v = args.Get("csv_s")) csv_s = std::string(*v);

    sjs::csv::Dialect d;
    d.sep = sep;
    d.write_header = true;

    if (!sjs::csv::WriteBoxes<2, sjs::Scalar>(csv_r, ds.R, d, &err)) {
      SJS_LOG_ERROR("CSV write R failed:", err);
      return 6;
    }
    if (!sjs::csv::WriteBoxes<2, sjs::Scalar>(csv_s, ds.S, d, &err)) {
      SJS_LOG_ERROR("CSV write S failed:", err);
      return 6;
    }
    SJS_LOG_INFO("Wrote CSV:", csv_r, "and", csv_s);
  }

  // Write generator report JSON-lite.
  {
    const std::string rep_path = (out_dir / (base + "_gen_report.json")).string();
    std::ofstream out(rep_path);
    if (out) {
      out << rep.ToJsonLite() << "\n";
      SJS_LOG_INFO("Wrote generator report:", rep_path);
    }
  }

  return 0;
}
