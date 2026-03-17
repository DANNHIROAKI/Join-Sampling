// apps/sjs_run.cpp
//
// Single experiment run:
//   - load/generate dataset
//   - create baseline (method + variant)
//   - run protocol once (or multiple repeats) with controlled RNG seeds
//   - write raw results CSV (+ optional sampled pairs)

#include "baselines/baseline_factory.h"

#include "baselines/runners/enum_sampling_runner.h"
#include "baselines/runners/sampling_runner.h"

#include "core/config.h"
#include "core/dim_dispatch.h"
#include "core/logging.h"
#include "core/stats.h"
#include "core/timer.h"
#include "core/types.h"

#include "data/synthetic/generator.h"
#include "io/binary_io.h"
#include "io/csv_io.h"
#include "io/dataset.h"

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

inline fs::path ResolvePathByWalkingUp(const fs::path& p, int max_up = 6) {
  if (p.empty()) return p;
  if (p.is_absolute()) return p;
  if (fs::exists(p)) return p;

  fs::path base = fs::current_path();
  for (int i = 0; i < max_up; ++i) {
    const fs::path cand = base / p;
    if (fs::exists(cand)) return cand;

    if (!base.has_parent_path()) break;
    const fs::path parent = base.parent_path();
    if (parent == base) break;
    base = parent;
  }
  return p;
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

    const fs::path path_r = ResolvePathByWalkingUp(fs::path(cfg.dataset.path_r));
    const fs::path path_s = ResolvePathByWalkingUp(fs::path(cfg.dataset.path_s));
    sjs::binary::RelationFileInfo infoR, infoS;
    sjs::binary::BinaryReadOptions opt;
    opt.generate_ids_if_missing = true;
    opt.drop_empty = false;
    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(path_r.string(), &R, &infoR, opt, &local_err)) {
      SetErr(err, "Failed reading binary R from " + path_r.string() + ": " + local_err);
      return false;
    }
    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(path_s.string(), &S, &infoS, opt, &local_err)) {
      SetErr(err, "Failed reading binary S from " + path_s.string() + ": " + local_err);
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

    const fs::path path_r = ResolvePathByWalkingUp(fs::path(cfg.dataset.path_r));
    const fs::path path_s = ResolvePathByWalkingUp(fs::path(cfg.dataset.path_s));

    char sep = ',';
    if (auto v = cfg.run.extra.find("csv_sep"); v != cfg.run.extra.end()) {
      if (!v->second.empty()) {
        if (v->second == "tab" || v->second == "\\t") sep = '\t';
        else sep = v->second[0];
      }
    }

    if (!sjs::csv::ReadBoxesSimple<Dim, sjs::Scalar>(path_r.string(), &R, sep, /*has_header=*/true, &local_err)) {
      SetErr(err, "Failed reading CSV R from " + path_r.string() + ": " + local_err);
      return false;
    }
    if (!sjs::csv::ReadBoxesSimple<Dim, sjs::Scalar>(path_s.string(), &S, sep, /*has_header=*/true, &local_err)) {
      SetErr(err, "Failed reading CSV S from " + path_s.string() + ": " + local_err);
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

template <int Dim>
inline bool WriteResultRow(sjs::csv::Writer& w,
                           const sjs::Config& cfg,
                           const sjs::Dataset<Dim, sjs::Scalar>& ds,
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
      << "  --method=<ours|kd_tree>\n"
      << "  --variant=<sampling|enum_sampling>\n"
      << "  --t=<num_samples>\n"
      << "  --seed=<seed>           (base seed; repeats add +rep)\n"
      << "  --repeats=<k>\n"
      << "  --out_dir=<dir>\n"
      << "  --results_file=<path>   (default: <out_dir>/run.csv)\n"
      << "  --write_samples=0|1\n"
      << "  --enum_cap=<cap>        (EnumSampling; 0 means no cap)\n"
      << "\nDataset flags:\n"
      << "  --dataset_source=<synthetic|binary|csv>\n"
      << "  --dataset=<name>\n"
      << "  --dim=<2|3|4|5>\n"
      << "  --path_r=<file> --path_s=<file>    (for binary/csv)\n"
      << "\nSynthetic flags:\n"
      << "  --gen=<alacarte_rectgen|alacarte-rectgen|rectgen|alacarte>\n"
      << "  --n_r=<N> --n_s=<N>\n"
      << "  --alpha=<float>\n"
      << "  --gen_seed=<seed>\n"
      << "  plus generator-specific params, e.g. --rectgen_script=... --audit_pairs=2000000\n"
      << "\nBaselines supported in this build (Dim=2/3/4/5):\n"
      << sjs::baselines::BaselineHelpSupportedDims()
      << "\n";
}

template <int Dim>
int RunMainForDim(const sjs::Config& cfg, const sjs::ArgMap& args) {
  std::string err;

  sjs::Dataset<Dim, sjs::Scalar> ds;
  sjs::synthetic::Report gen_report;
  sjs::synthetic::Report* gen_report_ptr = nullptr;

  {
    std::string local_err;
    if (!LoadOrGenerateDataset<Dim>(cfg, &ds, &gen_report, &local_err)) {
      SJS_LOG_ERROR(local_err);
      return 3;
    }
    if (cfg.dataset.source == sjs::DataSource::Synthetic) {
      gen_report_ptr = &gen_report;
      SJS_LOG_INFO("Generated dataset:", ds.name,
                   "dim=", Dim,
                   "R=", static_cast<unsigned long long>(ds.R.Size()),
                   "S=", static_cast<unsigned long long>(ds.S.Size()),
                   "gen=", cfg.dataset.synthetic.generator,
                   "report=", gen_report.ToJsonLite());
    } else {
      SJS_LOG_INFO("Loaded dataset:", ds.name,
                   "dim=", Dim,
                   "R=", static_cast<unsigned long long>(ds.R.Size()),
                   "S=", static_cast<unsigned long long>(ds.S.Size()));
    }
  }

  auto baseline = sjs::baselines::CreateBaseline<Dim>(cfg, &err);
  if (!baseline) {
    SJS_LOG_ERROR("CreateBaseline failed:", err);
    PrintUsage();
    return 4;
  }
  SJS_LOG_INFO("Baseline:", baseline->Name(),
               "method=", sjs::ToString(baseline->method()),
               "variant=", sjs::ToString(baseline->variant()),
               "dim=", Dim);

  const fs::path out_dir(cfg.output.out_dir);
  if (!EnsureDir(out_dir, &err)) {
    SJS_LOG_ERROR("Cannot create out_dir:", err);
    return 5;
  }

  std::string results_path = cfg.output.out_dir + "/run.csv";
  if (auto v = args.Get("results_file")) results_path = std::string(*v);

  if (!EnsureDir(fs::path(results_path).parent_path(), &err)) {
    SJS_LOG_ERROR("Cannot create results parent dir:", err);
    return 5;
  }

  sjs::csv::Writer writer(results_path, sjs::csv::Dialect{','}, &err);
  if (!writer.Ok()) {
    SJS_LOG_ERROR("Cannot open results file:", results_path, "err=", err);
    return 5;
  }

  if (!writer.WriteHeader(ResultHeader(), &err)) {
    SJS_LOG_ERROR("Failed writing CSV header:", err);
    return 5;
  }

  fs::path samples_dir = out_dir / "samples";
  const bool write_samples = cfg.run.write_samples;
  if (write_samples) {
    if (!EnsureDir(samples_dir, &err)) {
      SJS_LOG_ERROR("Cannot create samples_dir:", err);
      return 5;
    }
  }

  std::vector<double> wall_ms_all;
  std::vector<double> wall_ms_ok;
  sjs::u64 ok_reps = 0;
  wall_ms_all.reserve(static_cast<sjs::usize>(cfg.run.repeats));
  wall_ms_ok.reserve(static_cast<sjs::usize>(cfg.run.repeats));

  for (sjs::u64 rep = 0; rep < cfg.run.repeats; ++rep) {
    const sjs::u64 seed = cfg.run.seed + rep;

    sjs::baselines::RunReport report;
    std::string local_err;

    sjs::Stopwatch sw;
    bool ok = false;

    switch (cfg.run.variant) {
      case sjs::Variant::Sampling:
        ok = sjs::baselines::RunSamplingOnce<Dim, sjs::Scalar>(
            baseline.get(), ds, cfg, seed, &report, &local_err);
        break;
      case sjs::Variant::EnumSampling:
        ok = sjs::baselines::RunEnumSamplingOnce<Dim, sjs::Scalar>(
            baseline.get(), ds, cfg, seed, &report, &local_err);
        break;
      case sjs::Variant::Adaptive:
        ok = false;
        local_err = "Adaptive variant has been removed; use --variant=sampling or --variant=enum_sampling.";
        break;
    }

    const double wall_ms = sw.ElapsedMillis();
    wall_ms_all.push_back(wall_ms);
    if (ok) {
      wall_ms_ok.push_back(wall_ms);
      ++ok_reps;
    }

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

    if (!WriteResultRow(writer, cfg, ds, report, static_cast<int>(rep), wall_ms, gen_report_ptr)) {
      SJS_LOG_ERROR("Failed writing result row to CSV:", results_path);
      return 6;
    }

    if (write_samples && report.ok && !report.samples.pairs.empty()) {
      const std::string base =
          SanitizeFilename(ds.name) + "__" +
          SanitizeFilename(sjs::ToString(report.method)) + "__" +
          SanitizeFilename(sjs::ToString(report.variant)) + "__t" +
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

  const sjs::Summary sum_all = sjs::Summarize(wall_ms_all);
  SJS_LOG_INFO("Wall-time summary (ms) [all]:", sum_all.ToJson());

  if (!wall_ms_ok.empty()) {
    const sjs::Summary sum_ok = sjs::Summarize(wall_ms_ok);
    const double ok_rate = (cfg.run.repeats > 0)
                               ? (static_cast<double>(ok_reps) / static_cast<double>(cfg.run.repeats))
                               : 0.0;
    SJS_LOG_INFO("Wall-time summary (ms) [ok-only]:", sum_ok.ToJson(), "ok_rate=", ok_rate);
  } else {
    SJS_LOG_WARN("No successful repetitions; ok_rate=0.0");
  }

  SJS_LOG_INFO("Wrote results to:", results_path);
  if (write_samples) {
    SJS_LOG_INFO("Wrote samples to:", samples_dir.string());
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

  return sjs::DispatchSupportedDim(cfg.dataset.dim, [&](auto tag) {
    constexpr int Dim = decltype(tag)::value;
    return sjs::apps::RunMainForDim<Dim>(cfg, args);
  });
}
