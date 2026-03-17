#pragma once
// baselines/runners/sampling_runner.h
//
// Runner for Variant::Sampling.
//
// This enforces a consistent protocol across baselines:
//   Reset -> Build -> Count -> Sample
// and records top-level phase timings.
//
// Baselines may additionally record sub-phases into the same PhaseRecorder.

#include "baselines/baseline_api.h"
#include "core/assert.h"
#include "core/logging.h"

#include <string>

namespace sjs {
namespace baselines {

template <int Dim, class T = Scalar>
bool RunSamplingOnce(IBaseline<Dim, T>* baseline,
                     const Dataset<Dim, T>& dataset,
                     const Config& cfg,
                     u64 seed,
                     RunReport* out,
                     std::string* err = nullptr) {
  if (!baseline) {
    if (err) *err = "RunSamplingOnce: baseline is null";
    return false;
  }
  if (!out) {
    if (err) *err = "RunSamplingOnce: out is null";
    return false;
  }

  out->ok = false;
  out->error.clear();
  out->method = baseline->method();
  out->variant = Variant::Sampling;
  out->baseline_name = std::string(baseline->Name());
  out->dataset_name = dataset.name;
  out->seed = seed;
  out->t = cfg.run.t;
  out->count = CountResult{};
  out->samples.Clear();
  out->used_enumeration = false;
  out->enumeration_truncated = false;
  out->enumeration_cap = 0;
  out->enumeration_pairs_pass1 = 0;
  out->enumeration_pairs_pass2 = 0;
  out->enum_stats_pass1.Reset();
  out->enum_stats_pass2.Reset();
  out->adaptive_branch.clear();
  out->adaptive_pilot_pairs = 0;
  out->note.clear();
  out->phases.Clear();

  std::string local_err;

  // SJS v3: use independent randomness streams across phases.
  // (This improves reproducibility and avoids unintended coupling between
  // Count() and Sample() for adaptive baselines that use RNG in both.)
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

  // Defensive validation (weights etc.).
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
