#pragma once
// sjs/baselines/runners/adaptive_runner.h
//
// Runner for Variant::Adaptive.
//
// SJS v3 alignment
// ---------------
// In SJS v3, Adaptive corresponds to "Framework III" (budgeted caching / prefetch)
// to reduce (or eliminate) second-pass work without changing i.i.d. uniform targets.
//
// This runner enforces the same top-level protocol as other runners and delegates
// the internal caching strategy to the baseline:
//   Reset -> Build -> Count -> Sample
//
// Budget/caching knobs are baseline-specific and read from cfg.run,
// e.g., cfg.run.j_star as budget B, cfg.run.extra[...] for thresholds.

#include "sjs/baselines/baseline_api.h"
#include "sjs/core/assert.h"
#include "sjs/core/logging.h"

#include <string>

namespace sjs {
namespace baselines {

template <int Dim, class T = Scalar>
bool RunAdaptiveOnce(IBaseline<Dim, T>* baseline,
                     const Dataset<Dim, T>& dataset,
                     const Config& cfg,
                     u64 seed,
                     RunReport* out,
                     std::string* err = nullptr) {
  if (!baseline) {
    if (err) *err = "RunAdaptiveOnce: baseline is null";
    return false;
  }
  if (!out) {
    if (err) *err = "RunAdaptiveOnce: out is null";
    return false;
  }

  out->ok = false;
  out->error.clear();
  out->method = baseline->method();
  out->variant = Variant::Adaptive;
  out->baseline_name = std::string(baseline->Name());
  out->dataset_name = dataset.name;
  out->seed = seed;
  out->t = cfg.run.t;

  out->count = CountResult{};
  out->samples.Clear();

  out->used_enumeration = false;
  out->enumeration_truncated = false;
  out->enumeration_cap = cfg.run.enum_cap;

  out->enumeration_pairs_pass1 = 0;
  out->enumeration_pairs_pass2 = 0;
  out->enum_stats_pass1 = join::JoinStats{};
  out->enum_stats_pass2 = join::JoinStats{};

  out->adaptive_branch = "framework_III";
  out->adaptive_pilot_pairs = 0;
  out->note.clear();
  out->phases.Clear();

  std::string local_err;

  // SJS v3: separate RNG streams across phases.
  Rng rng_count(DeriveSeed(seed, 1));
  Rng rng_sample(DeriveSeed(seed, 2));

  baseline->Reset();

  {
    auto _ = out->phases.Scoped("run_build");
    if (!baseline->Build(dataset, cfg, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  {
    auto _ = out->phases.Scoped("run_count");
    if (!baseline->Count(cfg, &rng_count, &out->count, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  {
    auto _ = out->phases.Scoped("run_sample");
    if (!baseline->Sample(cfg, &rng_sample, &out->samples, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  // Defensive validation.
  {
    std::string v_err;
    if (!out->samples.Validate(&v_err)) {
      out->error = "SampleSet validation failed: " + v_err;
      if (err) *err = out->error;
      return false;
    }
  }

  out->ok = true;
  return true;
}

}  // namespace baselines
}  // namespace sjs
