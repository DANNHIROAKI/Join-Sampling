// src/apps/sjs_sweep.cpp
//
// Parameter sweep harness.
//
// This app repeatedly calls the single-run protocol (same as sjs_run),
// but over a grid/list of parameters, and writes:
//
//   (1) sweep_raw.csv     : one row per run (repeat)
//   (2) sweep_summary.csv : aggregated stats per (dataset params, method, variant, t)
//
// The sweep interface is intentionally lightweight: lists are comma-separated.
//
// Examples:
//   ./sjs_sweep --dataset_source=synthetic --gen=stripe --dataset=alpha_sweep \
//               --n_r=200000 --n_s=200000 --gen_seed=1 \
//               --sweep_alpha=1e-8,3e-8,1e-7,3e-7,1e-6 \
//               --sweep_methods=ours,aabb,interval_tree \
//               --sweep_variants=sampling,enum_sampling,adaptive \
//               --sweep_t=10000,100000 --repeats=3 --seed=1 \
//               --out_dir=out/alpha_sweep
//
// Notes:
//   - For synthetic datasets, dataset is re-generated for each (n_r,n_s,alpha,gen_seed).
//   - For binary/csv datasets, dataset is loaded once; only method/variant/t/seed vary.

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
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
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

inline std::string Trim(std::string_view s) {
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.remove_suffix(1);
  return std::string(s);
}

inline std::vector<std::string> SplitList(std::string_view s, char sep = ',') {
  std::vector<std::string> out;
  std::string_view rest = s;
  while (!rest.empty()) {
    const usize pos = rest.find(sep);
    std::string_view tok = (pos == std::string_view::npos) ? rest : rest.substr(0, pos);
    tok = std::string_view(Trim(tok));
    if (!tok.empty()) out.emplace_back(tok);
    if (pos == std::string_view::npos) break;
    rest.remove_prefix(pos + 1);
  }
  return out;
}

inline std::vector<double> ParseDoubleList(std::string_view s) {
  std::vector<double> out;
  for (const auto& tok : SplitList(s)) {
    try {
      out.push_back(std::stod(tok));
    } catch (...) {
      // ignore invalid tokens; caller can validate lengths if needed
    }
  }
  return out;
}

inline std::vector<u64> ParseU64List(std::string_view s) {
  std::vector<u64> out;
  for (const auto& tok : SplitList(s)) {
    try {
      std::size_t idx = 0;
      unsigned long long v = std::stoull(tok, &idx, 10);
      if (idx == tok.size()) out.push_back(static_cast<u64>(v));
    } catch (...) {
      // ignore
    }
  }
  return out;
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

inline std::string FormatSci(double x, int prec = 3) {
  std::ostringstream oss;
  oss << std::scientific << std::setprecision(prec) << x;
  return oss.str();
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
    sjs::binary::BinaryReadOptions opt;
    opt.generate_ids_if_missing = true;
    opt.drop_empty = false;
    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(cfg.dataset.path_r, &R, nullptr, opt, &local_err)) {
      SetErr(err, "Failed reading binary R: " + local_err);
      return false;
    }
    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(cfg.dataset.path_s, &S, nullptr, opt, &local_err)) {
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

inline std::vector<std::string> RawHeader() {
  return {
      "dataset",
      "generator",
      "alpha",
      "n_r",
      "n_s",
      "method",
      "variant",
      "t",
      "rep",
      "seed",
      "ok",
      "wall_ms",
      "count_value",
      "count_exact",
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

inline std::vector<std::string> SummaryHeader() {
  return {
      "dataset",
      "generator",
      "alpha",
      "n_r",
      "n_s",
      "method",
      "variant",
      "t",
      "repeats",
      "ok_rate",
      "wall_mean_ms",
      "wall_stdev_ms",
      "wall_p50_ms",
      "wall_p95_ms",
      "count_mean",
      "count_stdev",
      "count_exact_fraction",
      "note",
  };
}

inline void PrintUsage() {
  std::cerr
      << "sjs_sweep: parameter sweep runner\n\n"
      << "This app understands all flags of sjs_run, plus sweep lists:\n"
      << "  --sweep_methods=m1,m2,...\n"
      << "  --sweep_variants=v1,v2,...             (sampling|enum_sampling|adaptive)\n"
      << "  --sweep_t=t1,t2,...\n"
      << "  --sweep_seeds=s1,s2,...                (overrides --repeats)\n"
      << "\nSynthetic-only sweeps:\n"
      << "  --sweep_alpha=a1,a2,...\n"
      << "  --sweep_n_r=n1,n2,...  --sweep_n_s=m1,m2,...\n"
      << "\nOutput:\n"
      << "  --raw_file=<path>       (default: <out_dir>/sweep_raw.csv)\n"
      << "  --summary_file=<path>   (default: <out_dir>/sweep_summary.csv)\n"
      << "\nBaselines supported in this build (Dim=2):\n"
      << sjs::baselines::BaselineHelp2D()
      << "\n";
}

}  // namespace

}  // namespace apps
}  // namespace sjs

int main(int argc, char** argv) {
  using sjs::u64;
  using sjs::usize;
  sjs::ArgMap args = sjs::ArgMap::FromArgv(argc, argv);
  if (sjs::apps::IsHelpRequested(args)) {
    sjs::apps::PrintUsage();
    return 0;
  }

  sjs::Config cfg0 = sjs::Config::FromArgs(argc, argv);
  std::string err;
  if (!cfg0.Validate(&err)) {
    SJS_LOG_ERROR("Config validation failed:", err);
    sjs::apps::PrintUsage();
    return 2;
  }
  sjs::Logger::Instance().SetConfig(cfg0.logging);

  if (cfg0.dataset.dim != 2) {
    SJS_LOG_ERROR("This build currently supports only Dim=2. Got --dim=", cfg0.dataset.dim);
    return 2;
  }

  // --------------------------
  // Parse sweep lists (or defaults)
  // --------------------------
  std::vector<double> alpha_list;
  if (auto v = args.Get("sweep_alpha")) alpha_list = sjs::apps::ParseDoubleList(*v);
  if (alpha_list.empty()) alpha_list.push_back(cfg0.dataset.synthetic.alpha);

  std::vector<u64> n_r_list;
  if (auto v = args.Get("sweep_n_r")) n_r_list = sjs::apps::ParseU64List(*v);
  if (n_r_list.empty()) n_r_list.push_back(cfg0.dataset.synthetic.n_r);

  std::vector<u64> n_s_list;
  if (auto v = args.Get("sweep_n_s")) n_s_list = sjs::apps::ParseU64List(*v);
  if (n_s_list.empty()) n_s_list.push_back(cfg0.dataset.synthetic.n_s);

  std::vector<u64> t_list;
  if (auto v = args.Get("sweep_t")) t_list = sjs::apps::ParseU64List(*v);
  if (t_list.empty()) t_list.push_back(cfg0.run.t);

  // Methods
  std::vector<sjs::Method> methods;
  if (auto v = args.Get("sweep_methods")) {
    for (const auto& tok : sjs::apps::SplitList(*v)) {
      sjs::Method m = sjs::Method::Unknown;
      if (!sjs::ParseMethod(tok, &m) || m == sjs::Method::Unknown) {
        SJS_LOG_ERROR("Unknown method in --sweep_methods:", tok);
        return 2;
      }
      methods.push_back(m);
    }
  }
  if (methods.empty()) {
    methods.push_back(cfg0.run.method);
  }

  // Variants
  std::vector<sjs::Variant> variants;
  if (auto v = args.Get("sweep_variants")) {
    for (const auto& tok : sjs::apps::SplitList(*v)) {
      sjs::Variant vv = sjs::Variant::Sampling;
      if (!sjs::ParseVariant(tok, &vv)) {
        SJS_LOG_ERROR("Unknown variant in --sweep_variants:", tok);
        return 2;
      }
      variants.push_back(vv);
    }
  }
  if (variants.empty()) {
    variants.push_back(cfg0.run.variant);
  }

  // Seeds: if sweep_seeds provided, it overrides repeats.
  std::vector<u64> seed_list;
  if (auto v = args.Get("sweep_seeds")) seed_list = sjs::apps::ParseU64List(*v);

  // Output
  const fs::path out_dir(cfg0.output.out_dir);
  if (!sjs::apps::EnsureDir(out_dir, &err)) {
    SJS_LOG_ERROR("Cannot create out_dir:", err);
    return 5;
  }

  std::string raw_path = (out_dir / "sweep_raw.csv").string();
  if (auto v = args.Get("raw_file")) raw_path = std::string(*v);

  std::string summary_path = (out_dir / "sweep_summary.csv").string();
  if (auto v = args.Get("summary_file")) summary_path = std::string(*v);

  if (!sjs::apps::EnsureDir(fs::path(raw_path).parent_path(), &err)) {
    SJS_LOG_ERROR("Cannot create raw parent dir:", err);
    return 5;
  }
  if (!sjs::apps::EnsureDir(fs::path(summary_path).parent_path(), &err)) {
    SJS_LOG_ERROR("Cannot create summary parent dir:", err);
    return 5;
  }

  sjs::csv::Writer raw_writer(raw_path, sjs::csv::Dialect{','}, &err);
  if (!raw_writer.Ok()) {
    SJS_LOG_ERROR("Cannot open raw CSV:", raw_path, "err=", err);
    return 5;
  }
  sjs::csv::Writer sum_writer(summary_path, sjs::csv::Dialect{','}, &err);
  if (!sum_writer.Ok()) {
    SJS_LOG_ERROR("Cannot open summary CSV:", summary_path, "err=", err);
    return 5;
  }

  raw_writer.WriteHeader(sjs::apps::RawHeader(), &err);
  sum_writer.WriteHeader(sjs::apps::SummaryHeader(), &err);

  // For binary/csv, load dataset once (no dataset-param sweep).
  sjs::Dataset<2, sjs::Scalar> fixed_ds;
  sjs::synthetic::Report fixed_gen_report;
  bool has_fixed_ds = false;

  if (cfg0.dataset.source != sjs::DataSource::Synthetic) {
    std::string local_err;
    if (!sjs::apps::LoadOrGenerateDataset<2>(cfg0, &fixed_ds, &fixed_gen_report, &local_err)) {
      SJS_LOG_ERROR(local_err);
      return 3;
    }
    has_fixed_ds = true;
    SJS_LOG_INFO("Loaded fixed dataset:", fixed_ds.name,
                 "R=", static_cast<unsigned long long>(fixed_ds.R.Size()),
                 "S=", static_cast<unsigned long long>(fixed_ds.S.Size()));
  }

  // --------------------------
  // Sweep loops
  // --------------------------
  u64 total_runs = 0;
  u64 ok_runs = 0;

  for (u64 n_r : n_r_list) {
    for (u64 n_s : n_s_list) {
      for (double alpha : alpha_list) {
        // Generate dataset for this (n_r,n_s,alpha) if synthetic; otherwise reuse.
        sjs::Dataset<2, sjs::Scalar> ds;
        sjs::synthetic::Report gen_report;
        const sjs::synthetic::Report* gen_report_ptr = nullptr;

        if (cfg0.dataset.source == sjs::DataSource::Synthetic) {
          sjs::Config cfg_ds = cfg0;
          cfg_ds.dataset.synthetic.n_r = n_r;
          cfg_ds.dataset.synthetic.n_s = n_s;
          cfg_ds.dataset.synthetic.alpha = alpha;

          // Auto-name dataset so rows remain self-describing.
          {
            std::ostringstream nm;
            nm << cfg0.dataset.name
               << "__nr" << static_cast<unsigned long long>(n_r)
               << "__ns" << static_cast<unsigned long long>(n_s)
               << "__a" << sjs::apps::FormatSci(alpha, 3);
            cfg_ds.dataset.name = sjs::apps::SanitizeFilename(nm.str());
          }

          std::string local_err;
          if (!sjs::apps::LoadOrGenerateDataset<2>(cfg_ds, &ds, &gen_report, &local_err)) {
            SJS_LOG_ERROR("Dataset generation failed:", local_err);
            return 3;
          }
          gen_report_ptr = &gen_report;
        } else {
          ds = fixed_ds;
          gen_report_ptr = nullptr;
        }

        for (u64 t : t_list) {
          for (sjs::Method method : methods) {
            for (sjs::Variant variant : variants) {
              // Prepare cfg for this run group.
              sjs::Config cfg = cfg0;
              cfg.run.t = t;
              cfg.run.method = method;
              cfg.run.variant = variant;
              if (cfg.dataset.source == sjs::DataSource::Synthetic) {
                cfg.dataset.synthetic.n_r = n_r;
                cfg.dataset.synthetic.n_s = n_s;
                cfg.dataset.synthetic.alpha = alpha;
                cfg.dataset.name = ds.name;  // keep consistent with generated ds
              }

              // Determine repeats + seeds.
              std::vector<u64> seeds;
              if (!seed_list.empty()) {
                seeds = seed_list;
              } else {
                seeds.resize(static_cast<usize>(cfg.run.repeats));
                for (u64 rep = 0; rep < cfg.run.repeats; ++rep) seeds[static_cast<usize>(rep)] = cfg.run.seed + rep;
              }

              // Create baseline.
              std::string b_err;
              auto baseline = sjs::baselines::CreateBaseline2D(method, variant, &b_err);
              if (!baseline) {
                SJS_LOG_ERROR("CreateBaseline2D failed:", b_err);
                return 4;
              }

              // Collect per-repeat stats.
              std::vector<double> wall_ms;
              wall_ms.reserve(seeds.size());
              std::vector<double> count_vals;
              count_vals.reserve(seeds.size());
              u64 exact_count = 0;
              u64 ok_in_group = 0;

              for (usize i = 0; i < seeds.size(); ++i) {
                const u64 seed = seeds[i];
                sjs::baselines::RunReport rep_out;
                std::string local_err;

                sjs::Stopwatch sw;
                bool ok = false;
                switch (variant) {
                  case sjs::Variant::Sampling:
                    ok = sjs::baselines::RunSamplingOnce<2, sjs::Scalar>(
                        baseline.get(), ds, cfg, seed, &rep_out, &local_err);
                    break;
                  case sjs::Variant::EnumSampling:
                    ok = sjs::baselines::RunEnumSamplingOnce<2, sjs::Scalar>(
                        baseline.get(), ds, cfg, seed, &rep_out, &local_err);
                    break;
                  case sjs::Variant::Adaptive:
                    ok = sjs::baselines::RunAdaptiveOnce<2, sjs::Scalar>(
                        baseline.get(), ds, cfg, seed, &rep_out, &local_err);
                    break;
                }

                const double wms = sw.ElapsedMillis();
                wall_ms.push_back(wms);
                total_runs++;

                if (!ok) {
                  rep_out.ok = false;
                  rep_out.error = local_err;
                } else {
                  ok_runs++;
                  ok_in_group++;
                  if (rep_out.count.exact) exact_count++;
                }
                count_vals.push_back(static_cast<double>(rep_out.count.value));

                // Raw row.
                raw_writer.WriteRowV(
                    ds.name,
                    (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.generator : ""),
                    (cfg.dataset.source == sjs::DataSource::Synthetic ? alpha : std::numeric_limits<double>::quiet_NaN()),
                    static_cast<unsigned long long>(ds.R.Size()),
                    static_cast<unsigned long long>(ds.S.Size()),
                    sjs::ToString(method),
                    sjs::ToString(variant),
                    static_cast<unsigned long long>(t),
                    static_cast<unsigned long long>(i),
                    static_cast<unsigned long long>(seed),
                    (rep_out.ok ? 1 : 0),
                    wms,
                    static_cast<double>(rep_out.count.value),
                    (rep_out.count.exact ? 1 : 0),
                    (rep_out.used_enumeration ? 1 : 0),
                    (rep_out.enumeration_truncated ? 1 : 0),
                    static_cast<unsigned long long>(rep_out.enumeration_cap),
                    rep_out.adaptive_branch,
                    static_cast<unsigned long long>(rep_out.adaptive_pilot_pairs),
                    rep_out.phases.ToJsonMillis(),
                    rep_out.note,
                    cfg.ToJsonLite(),
                    (gen_report_ptr ? gen_report_ptr->ToJsonLite() : "{}"));
              }

              // Summary row.
              const sjs::Summary wall_sum = sjs::Summarize(wall_ms);
              const sjs::Summary cnt_sum = sjs::Summarize(count_vals);
              const double ok_rate = (seeds.empty() ? 0.0 : static_cast<double>(ok_in_group) /
                                                        static_cast<double>(seeds.size()));
              const double exact_frac = (seeds.empty() ? 0.0 : static_cast<double>(exact_count) /
                                                          static_cast<double>(seeds.size()));

              sum_writer.WriteRowV(
                  ds.name,
                  (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.generator : ""),
                  (cfg.dataset.source == sjs::DataSource::Synthetic ? alpha : std::numeric_limits<double>::quiet_NaN()),
                  static_cast<unsigned long long>(ds.R.Size()),
                  static_cast<unsigned long long>(ds.S.Size()),
                  sjs::ToString(method),
                  sjs::ToString(variant),
                  static_cast<unsigned long long>(t),
                  static_cast<unsigned long long>(seeds.size()),
                  ok_rate,
                  wall_sum.mean,
                  wall_sum.stdev,
                  wall_sum.median,
                  wall_sum.p95,
                  cnt_sum.mean,
                  cnt_sum.stdev,
                  exact_frac,
                  "");
            }
          }
        }
      }
    }
  }

  SJS_LOG_INFO("Sweep finished. total_runs=", static_cast<unsigned long long>(total_runs),
               " ok_runs=", static_cast<unsigned long long>(ok_runs));
  SJS_LOG_INFO("Raw CSV:", raw_path);
  SJS_LOG_INFO("Summary CSV:", summary_path);
  return 0;
}
