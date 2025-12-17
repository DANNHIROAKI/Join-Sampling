// src/apps/sjs_run.cpp
//
// Single experiment run:
//   - load/generate dataset
//   - create baseline (method + variant)
//   - run protocol once (or multiple repeats) with controlled RNG seeds
//   - write raw results CSV (+ optional sampled pairs)
//
// This app is intended to be the basic building block for larger sweeps.
//
// Examples:
//   ./sjs_run --dataset_source=synthetic --gen=stripe --dataset=demo \
//             --n_r=100000 --n_s=100000 --alpha=1e-6 --gen_seed=1 \
//             --method=ours --variant=sampling --t=100000 --seed=1 --repeats=3 \
//             --out_dir=out/demo
//
//   ./sjs_run --dataset_source=binary --path_r=data/R.bin --path_s=data/S.bin \
//             --method=aabb --variant=enum_sampling --t=10000 --enum_cap=500000 \
//             --out_dir=out/bin_case

#include "baselines/baseline_factory_2d.h"

#include "sjs/baselines/runners/adaptive_runner.h"
#include "sjs/baselines/runners/enum_sampling_runner.h"
#include "sjs/baselines/runners/sampling_runner.h"

#include "sjs/core/config.h"
#include "sjs/core/logging.h"
#include "sjs/core/stats.h"
#include "sjs/core/timer.h"
#include "sjs/core/types.h"

#include "sjs/data/synthetic/generator.h"
#include "sjs/io/binary_io.h"
#include "sjs/io/csv_io.h"
#include "sjs/io/dataset.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

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

template <int Dim>
bool LoadOrGenerateDataset(const sjs::Config& cfg,
                           sjs::Dataset<Dim, sjs::Scalar>* out,
                           sjs::synthetic::Report* gen_report,
                           std::string* err) {
  if (!out) {
    SetErr(err, "LoadOrGenerateDataset: out is null");
    return false;
  }

  if (cfg.dataset.dim != Dim) {
    SetErr(err, "Config dim mismatch: cfg.dataset.dim != compiled Dim");
    return false;
  }

  if (cfg.dataset.source == sjs::DataSource::Synthetic) {
    sjs::synthetic::DatasetSpec spec;
    spec.name = cfg.dataset.name;
    spec.n_r = cfg.dataset.synthetic.n_r;
    spec.n_s = cfg.dataset.synthetic.n_s;
    spec.alpha = cfg.dataset.synthetic.alpha;
    spec.seed = cfg.dataset.synthetic.seed;
    spec.params = cfg.dataset.synthetic.extra;

    sjs::synthetic::Report rep;
    std::string local_err;
    if (!sjs::synthetic::GenerateDataset<Dim, sjs::Scalar>(
            cfg.dataset.synthetic.generator, spec, out, &rep, &local_err)) {
      SetErr(err, "Synthetic generation failed: " + local_err);
      return false;
    }
    if (gen_report) *gen_report = rep;
    return true;
  }

  if (cfg.dataset.source == sjs::DataSource::Binary) {
    if (cfg.dataset.path_r.empty() || cfg.dataset.path_s.empty()) {
      SetErr(err, "Binary dataset requires --path_r and --path_s");
      return false;
    }
    sjs::Relation<Dim, sjs::Scalar> R, S;
    std::string local_err;
    sjs::binary::RelationFileInfo infoR, infoS;
    sjs::binary::BinaryReadOptions opt;
    opt.generate_ids_if_missing = true;
    opt.drop_empty = false;
    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(cfg.dataset.path_r, &R, &infoR, opt, &local_err)) {
      SetErr(err, "Failed reading binary R: " + local_err);
      return false;
    }
    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(cfg.dataset.path_s, &S, &infoS, opt, &local_err)) {
      SetErr(err, "Failed reading binary S: " + local_err);
      return false;
    }
    sjs::Dataset<Dim, sjs::Scalar> ds;
    ds.name = cfg.dataset.name.empty() ? "binary" : cfg.dataset.name;
    ds.half_open = true;
    ds.R = std::move(R);
    ds.S = std::move(S);

    std::string v;
    if (!ds.Validate(/*require_proper=*/true, &v)) {
      SetErr(err, "Binary dataset validation failed: " + v);
      return false;
    }
    *out = std::move(ds);
    return true;
  }

  if (cfg.dataset.source == sjs::DataSource::CSV) {
    if (cfg.dataset.path_r.empty() || cfg.dataset.path_s.empty()) {
      SetErr(err, "CSV dataset requires --path_r and --path_s");
      return false;
    }
    sjs::Relation<Dim, sjs::Scalar> R, S;
    std::string local_err;

    // Simple reader uses sep=',' by default. You can change to TSV by passing --csv_sep=\t.
    char sep = ',';
    if (auto v = cfg.run.extra.find("csv_sep"); v != cfg.run.extra.end()) {
      if (!v->second.empty()) {
        if (v->second == "tab" || v->second == "\\t") sep = '\t';
        else sep = v->second[0];
      }
    }

    if (!sjs::csv::ReadBoxesSimple<Dim, sjs::Scalar>(cfg.dataset.path_r, &R, sep, /*has_header=*/true, &local_err)) {
      SetErr(err, "Failed reading CSV R: " + local_err);
      return false;
    }
    if (!sjs::csv::ReadBoxesSimple<Dim, sjs::Scalar>(cfg.dataset.path_s, &S, sep, /*has_header=*/true, &local_err)) {
      SetErr(err, "Failed reading CSV S: " + local_err);
      return false;
    }

    sjs::Dataset<Dim, sjs::Scalar> ds;
    ds.name = cfg.dataset.name.empty() ? "csv" : cfg.dataset.name;
    ds.half_open = true;
    ds.R = std::move(R);
    ds.S = std::move(S);

    std::string v;
    if (!ds.Validate(/*require_proper=*/true, &v)) {
      SetErr(err, "CSV dataset validation failed: " + v);
      return false;
    }
    *out = std::move(ds);
    return true;
  }

  SetErr(err, "Unsupported dataset_source");
  return false;
}

inline std::vector<std::string> ResultHeader() {
  return {
      "dataset",
      "method",
      "variant",
      "rep",
      "seed",
      "n_r",
      "n_s",
      "t",
      "ok",
      "wall_ms",
      "count_value",
      "count_exact",
      "count_stderr",
      "count_ci_low",
      "count_ci_high",
      "used_enumeration",
      "enum_truncated",
      "enum_cap",
      "adaptive_branch",
      "adaptive_pilot_pairs",
      "phases_json",
      "note",
      "config_json",
      "gen_report_json",
  };
}

inline std::string MaybeNum(double x) {
  if (!std::isfinite(x)) return "";
  std::ostringstream oss;
  oss << x;
  return oss.str();
}

inline bool WriteResultRow(sjs::csv::Writer& w,
                           const sjs::Config& cfg,
                           const sjs::Dataset<2, sjs::Scalar>& ds,
                           const sjs::baselines::RunReport& r,
                           int rep,
                           double wall_ms,
                           const sjs::synthetic::Report* gen_report) {
  const auto method = sjs::ToString(r.method);
  const auto variant = sjs::ToString(r.variant);

  const std::string gen_json = gen_report ? gen_report->ToJsonLite() : "{}";

  return w.WriteRowV(
      ds.name,
      method,
      variant,
      rep,
      static_cast<unsigned long long>(r.seed),
      static_cast<unsigned long long>(ds.R.Size()),
      static_cast<unsigned long long>(ds.S.Size()),
      static_cast<unsigned long long>(r.t),
      (r.ok ? 1 : 0),
      wall_ms,
      static_cast<double>(r.count.value),
      (r.count.exact ? 1 : 0),
      static_cast<double>(r.count.stderr),
      static_cast<double>(r.count.ci_low),
      static_cast<double>(r.count.ci_high),
      (r.used_enumeration ? 1 : 0),
      (r.enumeration_truncated ? 1 : 0),
      static_cast<unsigned long long>(r.enumeration_cap),
      r.adaptive_branch,
      static_cast<unsigned long long>(r.adaptive_pilot_pairs),
      r.phases.ToJsonMillis(),
      r.note,
      cfg.ToJsonLite(),
      gen_json);
}

inline void PrintUsage() {
  std::cerr
      << "sjs_run: single experiment runner\n\n"
      << "Common flags:\n"
      << "  --method=<ours|aabb|interval_tree|kd_tree|r_tree|range_tree|pbsm|tlsop|sirs|rejection|tsunami>\n"
      << "  --variant=<sampling|enum_sampling|adaptive>\n"
      << "  --t=<num_samples>\n"
      << "  --seed=<seed>           (base seed; repeats add +rep)\n"
      << "  --repeats=<k>\n"
      << "  --out_dir=<dir>\n"
      << "  --results_file=<path>   (default: <out_dir>/run.csv)\n"
      << "  --write_samples=0|1\n"
      << "  --enum_cap=<cap>        (EnumSampling/Adaptive; 0 means no cap)\n"
      << "  --j_star=<threshold>    (Adaptive; 0 disables pilot enum)\n"
      << "\nDataset flags:\n"
      << "  --dataset_source=<synthetic|binary|csv>\n"
      << "  --dataset=<name>\n"
      << "  --dim=2\n"
      << "  --path_r=<file> --path_s=<file>    (for binary/csv)\n"
      << "\nSynthetic flags:\n"
      << "  --gen=<stripe|uniform|clustered|hetero_sizes>\n"
      << "  --n_r=<N> --n_s=<N>\n"
      << "  --alpha=<float>\n"
      << "  --gen_seed=<seed>\n"
      << "  plus generator-specific params, e.g. --gap_factor=0.1 --core_lo=0.45 ...\n"
      << "\nBaselines supported in this build (Dim=2):\n"
      << sjs::baselines::BaselineHelp2D()
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
    SJS_LOG_ERROR("This build currently supports only Dim=2. Got --dim=", cfg.dataset.dim);
    return 2;
  }

  // Load/generate dataset.
  sjs::Dataset<2, sjs::Scalar> ds;
  sjs::synthetic::Report gen_report;
  sjs::synthetic::Report* gen_report_ptr = nullptr;

  {
    std::string local_err;
    if (!sjs::apps::LoadOrGenerateDataset<2>(cfg, &ds, &gen_report, &local_err)) {
      SJS_LOG_ERROR(local_err);
      return 3;
    }
    if (cfg.dataset.source == sjs::DataSource::Synthetic) {
      gen_report_ptr = &gen_report;
      SJS_LOG_INFO("Generated dataset:", ds.name,
                   "R=", static_cast<unsigned long long>(ds.R.Size()),
                   "S=", static_cast<unsigned long long>(ds.S.Size()),
                   "gen=", cfg.dataset.synthetic.generator,
                   "report=", gen_report.ToJsonLite());
    } else {
      SJS_LOG_INFO("Loaded dataset:", ds.name,
                   "R=", static_cast<unsigned long long>(ds.R.Size()),
                   "S=", static_cast<unsigned long long>(ds.S.Size()));
    }
  }

  // Create baseline.
  auto baseline = sjs::baselines::CreateBaseline2D(cfg, &err);
  if (!baseline) {
    SJS_LOG_ERROR("CreateBaseline2D failed:", err);
    sjs::apps::PrintUsage();
    return 4;
  }
  SJS_LOG_INFO("Baseline:", baseline->Name(),
               "method=", sjs::ToString(baseline->method()),
               "variant=", sjs::ToString(baseline->variant()));

  // Output paths.
  const fs::path out_dir(cfg.output.out_dir);
  if (!sjs::apps::EnsureDir(out_dir, &err)) {
    SJS_LOG_ERROR("Cannot create out_dir:", err);
    return 5;
  }

  std::string results_path = cfg.output.out_dir + "/run.csv";
  if (auto v = args.Get("results_file")) results_path = std::string(*v);

  // Ensure parent dir exists.
  if (!sjs::apps::EnsureDir(fs::path(results_path).parent_path(), &err)) {
    SJS_LOG_ERROR("Cannot create results parent dir:", err);
    return 5;
  }

  sjs::csv::Writer writer(results_path, sjs::csv::Dialect{','}, &err);
  if (!writer.Ok()) {
    SJS_LOG_ERROR("Cannot open results file:", results_path, "err=", err);
    return 5;
  }

  if (!writer.WriteHeader(sjs::apps::ResultHeader(), &err)) {
    SJS_LOG_ERROR("Failed writing CSV header:", err);
    return 5;
  }

  // Optional samples directory.
  fs::path samples_dir = out_dir / "samples";
  const bool write_samples = cfg.run.write_samples;
  if (write_samples) {
    if (!sjs::apps::EnsureDir(samples_dir, &err)) {
      SJS_LOG_ERROR("Cannot create samples_dir:", err);
      return 5;
    }
  }

  // Run repeats.
  std::vector<double> wall_ms_all;
  wall_ms_all.reserve(static_cast<sjs::usize>(cfg.run.repeats));

  for (sjs::u64 rep = 0; rep < cfg.run.repeats; ++rep) {
    const sjs::u64 seed = cfg.run.seed + rep;

    sjs::baselines::RunReport report;
    std::string local_err;

    sjs::Stopwatch sw;
    bool ok = false;

    switch (cfg.run.variant) {
      case sjs::Variant::Sampling:
        ok = sjs::baselines::RunSamplingOnce<2, sjs::Scalar>(
            baseline.get(), ds, cfg, seed, &report, &local_err);
        break;
      case sjs::Variant::EnumSampling:
        ok = sjs::baselines::RunEnumSamplingOnce<2, sjs::Scalar>(
            baseline.get(), ds, cfg, seed, &report, &local_err);
        break;
      case sjs::Variant::Adaptive:
        ok = sjs::baselines::RunAdaptiveOnce<2, sjs::Scalar>(
            baseline.get(), ds, cfg, seed, &report, &local_err);
        break;
    }

    const double wall_ms = sw.ElapsedMillis();
    wall_ms_all.push_back(wall_ms);

    if (!ok) {
      SJS_LOG_ERROR("Run failed (rep=", rep, "):", local_err);
      report.ok = false;
      report.error = local_err;
    } else {
      SJS_LOG_INFO("Run ok (rep=", rep, ")",
                   "wall_ms=", wall_ms,
                   "count=", static_cast<double>(report.count.value),
                   (report.count.exact ? "(exact)" : "(est)"),
                   "samples=", static_cast<unsigned long long>(report.samples.Size()));
    }

    if (!sjs::apps::WriteResultRow(writer, cfg, ds, report, static_cast<int>(rep),
                                   wall_ms, gen_report_ptr)) {
      SJS_LOG_ERROR("Failed writing result row to CSV:", results_path);
      return 6;
    }

    if (write_samples && report.ok && !report.samples.pairs.empty()) {
      // Write sampled pairs as TSV.
      const std::string base =
          sjs::apps::SanitizeFilename(ds.name) + "__" +
          sjs::apps::SanitizeFilename(sjs::ToString(report.method)) + "__" +
          sjs::apps::SanitizeFilename(sjs::ToString(report.variant)) + "__t" +
          std::to_string(static_cast<unsigned long long>(report.t)) + "__seed" +
          std::to_string(static_cast<unsigned long long>(seed)) + "__rep" +
          std::to_string(static_cast<unsigned long long>(rep));
      const fs::path out_pairs = samples_dir / (base + ".tsv");

      std::string werr;
      if (!sjs::csv::WritePairs(out_pairs.string(),
                               sjs::Span<const sjs::PairId>(report.samples.pairs.data(),
                                                            report.samples.pairs.size()),
                               sjs::csv::Dialect{'\t'}, &werr)) {
        SJS_LOG_WARN("Failed writing samples:", out_pairs.string(), "err=", werr);
      }
    }
  }

  // Print a tiny summary to stdout/stderr.
  const sjs::Summary sum = sjs::Summarize(wall_ms_all);
  SJS_LOG_INFO("Wall-time summary (ms):", sum.ToJson());

  SJS_LOG_INFO("Wrote results to:", results_path);
  if (write_samples) {
    SJS_LOG_INFO("Wrote samples to:", samples_dir.string());
  }

  return 0;
}
