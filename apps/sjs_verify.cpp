// apps/sjs_verify.cpp
//
// Small-scale correctness + sampling-quality verification.
//
// What it does:
//   1) Load/generate a (small) dataset
//   2) Compute oracle join size |J| via brute-force O(|R||S|)
//   3) Run one baseline (method+variant) to get:
//        - count result (exact/estimate)
//        - sample pairs
//   4) If the oracle join universe is small enough, evaluate sample uniformity:
//        - missing_in_universe (should be 0)
//        - chi-square goodness-of-fit vs uniform over join pairs
//        - uniqueness/duplicates, L1/Linf deviations
//        - autocorrelation sanity check (hashed to [0,1))
//        - (optional) KS test on hashed pairs
//
// Intended usage: small n (e.g., n_r,n_s <= 5k) to keep oracle feasible.
//
// Example:
//   ./sjs_verify --dataset_source=synthetic --gen=stripe --dataset=verify
//                --n_r=2000 --n_s=2000 --alpha=1e-3 --gen_seed=1
//                --method=ours --variant=sampling --t=20000 --seed=1 --repeats=3

#include "baselines/baseline_factory_2d.h"

#include "sjs/baselines/runners/adaptive_runner.h"
#include "sjs/baselines/runners/enum_sampling_runner.h"
#include "sjs/baselines/runners/sampling_runner.h"

#include "sjs/core/config.h"
#include "sjs/core/logging.h"
#include "sjs/core/types.h"

#include "sjs/data/synthetic/generator.h"
#include "sjs/io/binary_io.h"
#include "sjs/io/csv_io.h"
#include "sjs/io/dataset.h"

#include "sjs/join/join_oracle.h"
#include "sjs/sampling/sample_quality.h"

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
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

inline bool ParseU64(const std::optional<std::string_view>& v, u64* out) {
  if (!v || !out) return false;
  try {
    std::size_t idx = 0;
    const std::string sv(*v);
    unsigned long long x = std::stoull(sv, &idx, 10);
    if (idx != sv.size()) return false;
    
    *out = static_cast<u64>(x);
    return true;
  } catch (...) {
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
      << "sjs_verify: small-scale correctness & sampling-quality checks\n\n"
      << "This app understands the sjs_run flags, plus oracle controls:\n"
      << "  --oracle_max_checks=<u64>     (default 50000000)  limit on |R|*|S|\n"
      << "  --oracle_collect_limit=<u64>  (default 1000000)   max |J| to fully collect pairs\n"
      << "  --oracle_cap=<u64>            (default 0)         collect at most this many pairs\n"
      << "\nTip: use small n_r,n_s (<= 5000) for feasible oracle.\n\n"
      << "Baselines supported in this build (Dim=2):\n"
      << sjs::baselines::BaselineHelp2D()
      << "\n";
}

inline double RelErr(double est, double truth) {
  if (truth == 0.0) return std::numeric_limits<double>::quiet_NaN();
  return (est - truth) / truth;
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

  u64 oracle_max_checks = 50'000'000ULL;
  u64 oracle_collect_limit = 1'000'000ULL;
  u64 oracle_cap = 0;

  sjs::apps::ParseU64(args.Get("oracle_max_checks"), &oracle_max_checks);
  sjs::apps::ParseU64(args.Get("oracle_collect_limit"), &oracle_collect_limit);
  sjs::apps::ParseU64(args.Get("oracle_cap"), &oracle_cap);

  // Load/generate dataset.
  sjs::Dataset<2, sjs::Scalar> ds;
  sjs::synthetic::Report gen_report;
  if (!sjs::apps::LoadOrGenerateDataset<2>(cfg, &ds, &gen_report, &err)) {
    SJS_LOG_ERROR("Dataset load/generate failed:", err);
    return 3;
  }

  const u64 n_r = static_cast<u64>(ds.R.Size());
  const u64 n_s = static_cast<u64>(ds.S.Size());
  const __uint128_t checks128 = static_cast<__uint128_t>(n_r) * static_cast<__uint128_t>(n_s);
  if (checks128 > static_cast<__uint128_t>(oracle_max_checks)) {
    SJS_LOG_ERROR("Oracle would require |R|*|S|=", static_cast<unsigned long long>(checks128),
                  " checks > oracle_max_checks=", static_cast<unsigned long long>(oracle_max_checks),
                  ". Reduce n or raise --oracle_max_checks.");
    return 4;
  }

  SJS_LOG_INFO("Dataset:", ds.name,
               "R=", static_cast<unsigned long long>(n_r),
               "S=", static_cast<unsigned long long>(n_s));

  // Oracle count.
  sjs::join::JoinStats oracle_stats;
  const u64 oracle_count = sjs::join::CountNaive<2, sjs::Scalar>(ds.R, ds.S, &oracle_stats);
  SJS_LOG_INFO("Oracle |J| =", static_cast<unsigned long long>(oracle_count),
               "candidate_checks=", static_cast<unsigned long long>(oracle_stats.candidate_checks));

  // Oracle universe (optional).
  std::vector<sjs::PairId> universe;
  bool have_full_universe = false;

  if (oracle_count <= oracle_collect_limit) {
    universe = sjs::join::CollectNaivePairs<2, sjs::Scalar>(ds.R, ds.S, /*cap=*/0, nullptr);
    have_full_universe = (universe.size() == static_cast<sjs::usize>(oracle_count));
    SJS_LOG_INFO("Collected full universe pairs:", static_cast<unsigned long long>(universe.size()));
  } else if (oracle_cap > 0) {
    universe = sjs::join::CollectNaivePairs<2, sjs::Scalar>(ds.R, ds.S, oracle_cap, nullptr);
    have_full_universe = false;
    SJS_LOG_WARN("Join too large to fully collect (|J|=", static_cast<unsigned long long>(oracle_count),
                 "). Collected cap=", static_cast<unsigned long long>(oracle_cap),
                 " pairs; will skip uniformity tests.");
  } else {
    SJS_LOG_WARN("Join too large to collect (|J|=", static_cast<unsigned long long>(oracle_count),
                 "). Will skip uniformity tests.");
  }

  // Baseline
  auto baseline = sjs::baselines::CreateBaseline2D(cfg, &err);
  if (!baseline) {
    SJS_LOG_ERROR("CreateBaseline2D failed:", err);
    return 5;
  }

  // Repeat runs.
  const u64 repeats = cfg.run.repeats;
  for (u64 rep = 0; rep < repeats; ++rep) {
    const u64 seed = cfg.run.seed + rep;

    sjs::baselines::RunReport out;
    std::string local_err;

    bool ok = false;
    switch (cfg.run.variant) {
      case sjs::Variant::Sampling:
        ok = sjs::baselines::RunSamplingOnce<2, sjs::Scalar>(baseline.get(), ds, cfg, seed, &out, &local_err);
        break;
      case sjs::Variant::EnumSampling:
        ok = sjs::baselines::RunEnumSamplingOnce<2, sjs::Scalar>(baseline.get(), ds, cfg, seed, &out, &local_err);
        break;
      case sjs::Variant::Adaptive:
        ok = sjs::baselines::RunAdaptiveOnce<2, sjs::Scalar>(baseline.get(), ds, cfg, seed, &out, &local_err);
        break;
    }

    std::cout << "---- run rep=" << rep << " seed=" << seed << " ----\n";
    if (!ok) {
      std::cout << "FAILED: " << local_err << "\n";
      continue;
    }

    const double est = static_cast<double>(out.count.value);
    const double truth = static_cast<double>(oracle_count);
    const double rel = sjs::apps::RelErr(est, truth);

    std::cout << "method=" << sjs::ToString(out.method)
              << " variant=" << sjs::ToString(out.variant)
              << " t=" << out.t << "\n";
    std::cout << "count=" << std::setprecision(17) << est
              << (out.count.exact ? " (exact)" : " (est)")
              << "  oracle=" << truth
              << "  rel_err=" << rel << "\n";
    std::cout << "samples=" << out.samples.Size() << "\n";

    if (have_full_universe && !universe.empty() && !out.samples.pairs.empty()) {
      const auto uni = sjs::sampling::quality::EvaluatePairUniformity(
          sjs::Span<const sjs::PairId>(universe.data(), universe.size()),
          sjs::Span<const sjs::PairId>(out.samples.pairs.data(), out.samples.pairs.size()));

      const double ac1 = sjs::sampling::quality::AutocorrelationHashedPairs(
          sjs::Span<const sjs::PairId>(out.samples.pairs.data(), out.samples.pairs.size()), /*lag=*/1);

      // KS sanity check: map samples into [0,1) using a universe-rank + jitter scheme
      // to avoid systematic rejection on small discrete universes (see sample_quality.h).
      const auto ks = sjs::sampling::quality::KSPairsHashUniform01RankJitter(
          sjs::Span<const sjs::PairId>(universe.data(), universe.size()),
          sjs::Span<const sjs::PairId>(out.samples.pairs.data(), out.samples.pairs.size()));

      std::cout << "quality:\n";
      std::cout << "  missing_in_universe=" << uni.missing_in_universe << "\n";
      std::cout << "  unique_fraction=" << uni.unique_fraction
                << " duplicate_fraction=" << uni.duplicate_fraction << "\n";
      std::cout << "  l1=" << uni.l1 << " l_inf=" << uni.l_inf
                << " max_rel_error=" << uni.max_rel_error << "\n";
      std::cout << "  chi2_stat=" << uni.chi2.statistic
                << " df=" << uni.chi2.df
                << " p_value=" << uni.chi2.p_value << "\n";
      std::cout << "  autocorr_hash_lag1=" << ac1 << "\n";
      std::cout << "  ks_hash_uniform01 D=" << ks.D << " p=" << ks.p_value << "\n";
    } else {
      std::cout << "quality: skipped (universe not collected)\n";
    }
  }

  return 0;
}