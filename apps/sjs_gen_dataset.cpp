// apps/sjs_gen_dataset.cpp
//
// Synthetic dataset generator + exporter.
//
// Generates a Dataset<Dim> (R/S) using the configured synthetic generator, then
// writes it out in:
//   - SJS binary format (fast, robust), and optionally
//   - CSV (for debugging/visualization)

#include "core/config.h"
#include "core/dim_dispatch.h"
#include "core/logging.h"
#include "core/types.h"

#include "data/synthetic/generator.h"
#include "io/binary_io.h"
#include "io/csv_io.h"
#include "io/dataset.h"

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
      << "  --gen=<alacarte_rectgen|alacarte-rectgen|rectgen|alacarte>\n"
      << "  --dataset=<name>\n"
      << "  --dim=<2|3|4|5>\n"
      << "  --n_r=<N> --n_s=<N>\n"
      << "  --alpha=<float>\n"
      << "  --gen_seed=<seed>\n"
      << "  [optional generator params: --rectgen_script=... --audit_pairs=... --audit_seed=...]\n"
      << "\nOutput:\n"
      << "  --out_dir=<dir>\n"
      << "  --out_r=<path> --out_s=<path>            (override binary outputs)\n"
      << "  --write_csv=0|1                          (default 0)\n"
      << "  --csv_r=<path> --csv_s=<path>            (override CSV outputs)\n"
      << "  --csv_sep=,|tab|\\t                       (default ',')\n"
      << "\n";
}

template <int Dim>
int RunGenerateDataset(const sjs::Config& cfg, const sjs::ArgMap& args) {
  static_assert(Dim >= 2 && Dim <= 5, "RunGenerateDataset only supports Dim=2..5");

  std::string err;
  sjs::Dataset<Dim, sjs::Scalar> ds;
  sjs::synthetic::DatasetSpec spec;
  spec.name = cfg.dataset.name;
  spec.n_r = cfg.dataset.synthetic.n_r;
  spec.n_s = cfg.dataset.synthetic.n_s;
  spec.alpha = cfg.dataset.synthetic.alpha;
  spec.seed = cfg.dataset.synthetic.seed;
  spec.params = cfg.dataset.synthetic.extra;

  sjs::synthetic::Report rep;
  if (!sjs::synthetic::GenerateDataset<Dim, sjs::Scalar>(cfg.dataset.synthetic.generator, spec, &ds, &rep, &err)) {
    SJS_LOG_ERROR("Generation failed:", err);
    return 3;
  }

  if (ds.name.empty()) ds.name = spec.name.empty() ? cfg.dataset.name : spec.name;
  if (ds.name.empty()) ds.name = "synthetic";

  SJS_LOG_INFO("Generated dataset:", ds.name,
               "dim=", Dim,
               "R=", static_cast<unsigned long long>(ds.R.Size()),
               "S=", static_cast<unsigned long long>(ds.S.Size()),
               "report=", rep.ToJsonLite());

  fs::path out_dir(cfg.output.out_dir);
  if (!args.Has("out_dir")) {
    out_dir = fs::path("data/synthetic");
  }

  if (!EnsureDir(out_dir, &err)) {
    SJS_LOG_ERROR("Cannot create out_dir:", err);
    return 4;
  }

  const std::string base = SanitizeFilename(ds.name);

  std::string out_r = (out_dir / (base + "_R.bin")).string();
  std::string out_s = (out_dir / (base + "_S.bin")).string();
  if (auto v = args.Get("out_r")) out_r = std::string(*v);
  if (auto v = args.Get("out_s")) out_s = std::string(*v);

  {
    sjs::binary::BinaryWriteOptions opt;
    opt.half_open = true;
    opt.write_ids = true;
    opt.scalar = sjs::binary::ScalarEncoding::Float64;
    opt.write_name = true;

    if (!sjs::binary::WriteDatasetBinaryPair<Dim, sjs::Scalar>(out_r, out_s, ds, opt, &err)) {
      SJS_LOG_ERROR("Binary write failed:", err);
      return 5;
    }
    SJS_LOG_INFO("Wrote binary:", out_r, "and", out_s);
  }

  const bool write_csv = ParseBool(args.Get("write_csv"), /*def=*/false);
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

    if (!sjs::csv::WriteBoxes<Dim, sjs::Scalar>(csv_r, ds.R, d, &err)) {
      SJS_LOG_ERROR("CSV write R failed:", err);
      return 6;
    }
    if (!sjs::csv::WriteBoxes<Dim, sjs::Scalar>(csv_s, ds.S, d, &err)) {
      SJS_LOG_ERROR("CSV write S failed:", err);
      return 6;
    }
    SJS_LOG_INFO("Wrote CSV:", csv_r, "and", csv_s);
  }

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

  if (!sjs::IsSupportedExperimentDim(cfg.dataset.dim)) {
    SJS_LOG_ERROR(sjs::SupportedDimError(cfg.dataset.dim));
    return 2;
  }
  if (cfg.dataset.source != sjs::DataSource::Synthetic) {
    SJS_LOG_ERROR("sjs_gen_dataset requires --dataset_source=synthetic");
    return 2;
  }

  return sjs::DispatchSupportedDim(cfg.dataset.dim, [&](auto tag) {
    constexpr int Dim = decltype(tag)::value;
    return sjs::apps::RunGenerateDataset<Dim>(cfg, args);
  });
}
