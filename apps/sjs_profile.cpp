// apps/sjs_profile.cpp
//
// Profiling helper: run a single (dataset, method, variant, t) configuration
// and print a human-readable phase breakdown + key counters.
//
// This is not a micro-benchmark; it is designed to quickly answer:
//   - Where is the time going? (build / count / enumerate / sample)
//   - Is the run exact or estimated?
//   - Did adaptive branch enumerate-all or fallback?
//
// Example:
//   ./sjs_profile --dataset_source=synthetic --gen=stripe --dataset=demo
//     --n_r=200000 --n_s=200000 --alpha=1e-6 --gen_seed=1
//     --method=ours --variant=adaptive --t=100000 --seed=1
//     --out_dir=out/profile_demo

#include "baselines/baseline_factory_2d.h"

#include "sjs/baselines/runners/adaptive_runner.h"
#include "sjs/baselines/runners/enum_sampling_runner.h"
#include "sjs/baselines/runners/sampling_runner.h"

#include "sjs/core/config.h"
#include "sjs/core/logging.h"
#include "sjs/core/timer.h"

#include "sjs/data/synthetic/generator.h"
#include "sjs/io/binary_io.h"
#include "sjs/io/csv_io.h"
#include "sjs/io/dataset.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
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

inline bool ParseBool(const std::optional<std::string_view>& v, bool def) {
  if (!v) return def;
  bool b = def;
  if (sjs::detail::ParseBool(*v, &b)) return b;
  return def;
}

inline sjs::usize ParseUsize(const std::optional<std::string_view>& v, sjs::usize def) {
  if (!v) return def;
  sjs::u64 x = 0;
  if (sjs::detail::ParseU64(*v, &x)) return static_cast<sjs::usize>(x);
  return def;
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

inline void PrintUsage() {
  std::cerr
      << "sjs_profile: run once and print phase breakdown / counters\n\n"
      << "Common flags (same as sjs_run):\n"
      << "  --method=<ours|aabb|interval_tree|kd_tree|r_tree|range_tree|pbsm|tlsop|sirs|rejection|tsunami>\n"
      << "  --variant=<sampling|enum_sampling|adaptive>\n"
      << "  --t=<num_samples>\n"
      << "  --seed=<seed>\n"
      << "  --out_dir=<dir>\n"
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
      << "  plus generator-specific params\n"
      << "\nProfile-specific flags:\n"
      << "  --topk=<k>              (default 20) show top-k phases by time\n"
      << "  --out_file=<path>       (optional) write one JSONL profile record\n"
      << "  --print_json=0|1        (default 0) also print JSON to stdout\n"
      << "\nBaselines supported in this build (Dim=2):\n"
      << sjs::baselines::BaselineHelp2D()
      << "\n";
}

struct PhaseLine {
  std::string name;
  double ms{0.0};
};

inline void PrintPhases(const sjs::PhaseRecorder& pr, sjs::usize topk) {
  auto snap = pr.SnapshotNanosSorted();
  std::vector<PhaseLine> lines;
  lines.reserve(snap.size());
  double total_ms = 0.0;
  for (const auto& kv : snap) {
    const double ms = static_cast<double>(kv.second) / 1e6;
    total_ms += ms;
    lines.push_back(PhaseLine{kv.first, ms});
  }

  std::sort(lines.begin(), lines.end(), [](const PhaseLine& a, const PhaseLine& b) {
    return a.ms > b.ms;
  });

  std::cerr << "\nPhase breakdown (top " << static_cast<unsigned long long>(topk)
            << ", total=" << std::fixed << std::setprecision(3) << total_ms << " ms):\n";
  std::cerr << "  ms\t%\tphase\n";

  const sjs::usize k = std::min(topk, lines.size());
  for (sjs::usize i = 0; i < k; ++i) {
    const double pct = (total_ms > 0.0) ? (100.0 * lines[i].ms / total_ms) : 0.0;
    std::cerr << "  " << std::fixed << std::setprecision(3) << lines[i].ms
              << "\t" << std::setprecision(1) << pct
              << "\t" << lines[i].name << "\n";
  }

  if (lines.size() > k) {
    double rest_ms = 0.0;
    for (sjs::usize i = k; i < lines.size(); ++i) rest_ms += lines[i].ms;
    const double pct = (total_ms > 0.0) ? (100.0 * rest_ms / total_ms) : 0.0;
    std::cerr << "  " << std::fixed << std::setprecision(3) << rest_ms
              << "\t" << std::setprecision(1) << pct
              << "\t" << "(rest " << static_cast<unsigned long long>(lines.size() - k) << " phases)\n";
  }
}

inline std::string JsonEscape(std::string_view s) {
  std::string out;
  out.reserve(s.size() + 8);
  for (char c : s) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          // control char -> skip
        } else {
          out.push_back(c);
        }
    }
  }
  return out;
}

inline std::string ProfileRecordJson(const sjs::Config& cfg,
                                    const sjs::Dataset<2, sjs::Scalar>& ds,
                                    const sjs::baselines::RunReport& r,
                                    double wall_ms) {
  std::ostringstream oss;
  oss << "{";
  oss << "\"dataset\":\"" << JsonEscape(ds.name) << "\",";
  oss << "\"method\":\"" << sjs::ToString(r.method) << "\",";
  oss << "\"variant\":\"" << sjs::ToString(r.variant) << "\",";
  oss << "\"t\":" << static_cast<unsigned long long>(r.t) << ',';
  oss << "\"seed\":" << static_cast<unsigned long long>(r.seed) << ',';
  oss << "\"ok\":" << (r.ok ? 1 : 0) << ',';
  oss << "\"wall_ms\":" << std::fixed << std::setprecision(6) << wall_ms << ',';
  oss << "\"count\":" << static_cast<double>(r.count.value) << ',';
  oss << "\"count_exact\":" << (r.count.exact ? 1 : 0) << ',';
  oss << "\"count_stderr\":" << static_cast<double>(r.count.stderr) << ',';
  oss << "\"used_enumeration\":" << (r.used_enumeration ? 1 : 0) << ',';
  oss << "\"enum_truncated\":" << (r.enumeration_truncated ? 1 : 0) << ',';
  oss << "\"enum_cap\":" << static_cast<unsigned long long>(r.enumeration_cap) << ',';
  oss << "\"adaptive_branch\":\"" << JsonEscape(r.adaptive_branch) << "\",";
  oss << "\"adaptive_pilot_pairs\":" << static_cast<unsigned long long>(r.adaptive_pilot_pairs) << ',';
  oss << "\"enum_stats_pass1\":" << r.enum_stats_pass1.ToJsonLite() << ',';
  oss << "\"enum_stats_pass2\":" << r.enum_stats_pass2.ToJsonLite() << ',';
  oss << "\"phases\":" << r.phases.ToJsonMillis() << ',';
  oss << "\"config\":" << cfg.ToJsonLite();
  oss << "}";
  return oss.str();
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

  const sjs::usize topk = sjs::apps::ParseUsize(args.Get("topk"), 20);
  const bool print_json = sjs::apps::ParseBool(args.Get("print_json"), false);

  // Load/generate dataset.
  sjs::Dataset<2, sjs::Scalar> ds;
  sjs::synthetic::Report gen_report;
  if (!sjs::apps::LoadOrGenerateDataset<2>(cfg, &ds, &gen_report, &err)) {
    SJS_LOG_ERROR("Dataset load/generate failed:", err);
    return 3;
  }

  if (cfg.dataset.source == sjs::DataSource::Synthetic) {
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

  // Run once.
  sjs::baselines::RunReport report;
  std::string run_err;
  sjs::Stopwatch sw;

  bool ok = false;
  switch (cfg.run.variant) {
    case sjs::Variant::Sampling:
      ok = sjs::baselines::RunSamplingOnce<2, sjs::Scalar>(baseline.get(), ds, cfg, cfg.run.seed, &report, &run_err);
      break;
    case sjs::Variant::EnumSampling:
      ok = sjs::baselines::RunEnumSamplingOnce<2, sjs::Scalar>(baseline.get(), ds, cfg, cfg.run.seed, &report, &run_err);
      break;
    case sjs::Variant::Adaptive:
      ok = sjs::baselines::RunAdaptiveOnce<2, sjs::Scalar>(baseline.get(), ds, cfg, cfg.run.seed, &report, &run_err);
      break;
  }

  const double wall_ms = sw.ElapsedMillis();

  if (!ok) {
    SJS_LOG_ERROR("Run failed:", run_err);
    return 5;
  }

  // Print headline.
  std::cerr << "\n=== Profile summary ===\n";
  std::cerr << "dataset:  " << ds.name << " (|R|=" << static_cast<unsigned long long>(ds.R.Size())
            << ", |S|=" << static_cast<unsigned long long>(ds.S.Size()) << ")\n";
  std::cerr << "method:   " << sjs::ToString(report.method) << "\n";
  std::cerr << "variant:  " << sjs::ToString(report.variant) << "\n";
  std::cerr << "t:        " << static_cast<unsigned long long>(report.t) << "\n";
  std::cerr << "seed:     " << static_cast<unsigned long long>(report.seed) << "\n";
  std::cerr << "wall_ms:  " << std::fixed << std::setprecision(3) << wall_ms << "\n";
  std::cerr << "count:    " << static_cast<double>(report.count.value)
            << (report.count.exact ? " (exact)" : " (est)")
            << " stderr=" << static_cast<double>(report.count.stderr)
            << "\n";
  if (!report.adaptive_branch.empty()) {
    std::cerr << "adaptive: " << report.adaptive_branch
              << " pilot_pairs=" << static_cast<unsigned long long>(report.adaptive_pilot_pairs)
              << "\n";
  }
  if (report.used_enumeration) {
    std::cerr << "enum:     cap=" << static_cast<unsigned long long>(report.enumeration_cap)
              << " pass1=" << static_cast<unsigned long long>(report.enumeration_pairs_pass1)
              << " pass2=" << static_cast<unsigned long long>(report.enumeration_pairs_pass2)
              << " truncated=" << (report.enumeration_truncated ? 1 : 0)
              << "\n";
  }
  std::cerr << "samples:  " << static_cast<unsigned long long>(report.samples.Size())
            << (report.samples.with_replacement ? " (with_replacement)" : "")
            << (report.samples.weighted ? " (weighted)" : "")
            << "\n";

  // Print join stats if available.
  if (report.used_enumeration) {
    std::cerr << "enum_stats_pass1: " << report.enum_stats_pass1.ToJsonLite() << "\n";
    if (report.enumeration_pairs_pass2 > 0) {
      std::cerr << "enum_stats_pass2: " << report.enum_stats_pass2.ToJsonLite() << "\n";
    }
  }

  // Print phases.
  sjs::apps::PrintPhases(report.phases, topk);

  // Optional JSONL output.
  std::string out_file;
  if (auto v = args.Get("out_file")) {
    out_file = std::string(*v);
  } else if (args.Has("out_dir")) {
    out_file = (fs::path(cfg.output.out_dir) / "profile.jsonl").string();
  }

  if (!out_file.empty()) {
    if (!sjs::apps::EnsureDir(fs::path(out_file).parent_path(), &err)) {
      SJS_LOG_WARN("Cannot create out_file parent dir:", err);
    } else {
      std::ofstream ofs(out_file, std::ios::app);
      if (!ofs) {
        SJS_LOG_WARN("Cannot open out_file:", out_file);
      } else {
        const std::string rec = sjs::apps::ProfileRecordJson(cfg, ds, report, wall_ms);
        ofs << rec << "\n";
        SJS_LOG_INFO("Wrote profile record:", out_file);
        if (print_json) {
          std::cout << rec << "\n";
        }
      }
    }
  } else if (print_json) {
    std::cout << sjs::apps::ProfileRecordJson(cfg, ds, report, wall_ms) << "\n";
  }

  return 0;
}
