#pragma once
// sjs/baselines/runners/adaptive_runner.h
//
// Runner for Variant::Adaptive.
//
// Goal: provide a *baseline-agnostic* adaptive protocol so that every
// baseline's adaptive version can be evaluated consistently.
//
// Strategy implemented here (deterministic, simple, and extensible):
//   1) Reset -> Build
//   2) Pilot enumeration: enumerate up to (j_star + 1) join pairs (and also
//      respecting enum_cap if configured). This tells us whether |J| <= j_star.
//   3a) If join finishes within the pilot budget:
//         - We have enumerated all join pairs (|J| = N <= j_star) and stored them.
//         - Report exact count N.
//         - Sample t pairs uniformly WITH replacement from the stored list.
//      (No second enumeration pass needed.)
//   3b) Otherwise (join is large or unknown due to caps):
//         - Fall back to baseline-provided Count + Sample (sampling-based).
//
// Why this is reasonable:
//   - This matches common SIGMOD-style adaptive designs: "enumerate when output is small,
//     sample when output is huge".
//   - It makes the trigger explicit via cfg.run.j_star.
//   - It avoids bias by NOT using the pilot prefix as samples in the large-join branch.
//
// Notes:
//   - If cfg.run.j_star == 0, the runner skips pilot enumeration and directly falls back.
//   - If cfg.run.enum_cap is set and smaller than j_star, the pilot may be truncated;
//     in that case, the runner will likely fall back.

#include "sjs/baselines/baseline_api.h"
#include "sjs/core/assert.h"
#include "sjs/core/logging.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

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
  out->used_enumeration = true;  // pilot always tries to enumerate unless j_star==0
  out->enumeration_truncated = false;
  out->enumeration_cap = cfg.run.enum_cap;
  out->enumeration_pairs_pass1 = 0;
  out->enumeration_pairs_pass2 = 0;
  out->enum_stats_pass1.Reset();
  out->enum_stats_pass2.Reset();
  out->adaptive_branch.clear();
  out->adaptive_pilot_pairs = 0;
  out->note.clear();
  out->phases.Clear();

  std::string local_err;
  Rng rng(seed);

  baseline->Reset();

  {
    auto _ = out->phases.Scoped("run_build");
    if (!baseline->Build(dataset, cfg, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  const u64 t = cfg.run.t;
  const u64 j_star = cfg.run.j_star;
  const u64 enum_cap = cfg.run.enum_cap;

  if (j_star == 0) {
    // Degenerate config: never enumerate.
    out->used_enumeration = false;
    out->adaptive_branch = "fallback_sampling_no_pilot";

    {
      auto _ = out->phases.Scoped("run_fallback_count");
      if (!baseline->Count(cfg, &rng, &out->count, &out->phases, &local_err)) {
        out->error = local_err;
        if (err) *err = local_err;
        return false;
      }
    }

    {
      auto _ = out->phases.Scoped("run_fallback_sample");
      if (!baseline->Sample(cfg, &rng, &out->samples, &out->phases, &local_err)) {
        out->error = local_err;
        if (err) *err = local_err;
        return false;
      }
    }

    // Validation.
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

  // --------------------------
  // Pilot enumeration up to (j_star + 1) (and also enum_cap if set)
  // --------------------------
  std::unique_ptr<IJoinEnumerator> stream;
  {
    auto _ = out->phases.Scoped("run_pilot_enum_prepare");
    stream = baseline->Enumerate(cfg, &out->phases, &local_err);
    if (!stream) {
      if (local_err.empty()) local_err = "baseline->Enumerate returned null stream";
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  // Compute pilot limit safely (avoid overflow).
  u64 pilot_limit = j_star;
  if (pilot_limit != std::numeric_limits<u64>::max()) {
    pilot_limit = pilot_limit + 1;
  }
  if (enum_cap > 0) {
    pilot_limit = std::min(pilot_limit, enum_cap);
  }

  std::vector<PairId> pilot_pairs;
  pilot_pairs.reserve(static_cast<usize>(std::min<u64>(pilot_limit, 1'000'000ULL)));

  bool hit_limit = false;
  {
    auto _ = out->phases.Scoped("run_pilot_enum_scan");
    stream->Reset();
    PairId p;
    while (stream->Next(&p)) {
      pilot_pairs.push_back(p);
      if (pilot_pairs.size() >= static_cast<usize>(pilot_limit)) {
        hit_limit = true;
        break;
      }
    }
  }

  out->adaptive_pilot_pairs = static_cast<u64>(pilot_pairs.size());
  out->enumeration_pairs_pass1 = static_cast<u64>(pilot_pairs.size());
  out->enum_stats_pass1 = stream->Stats();

  if (!hit_limit) {
    // Join exhausted within the pilot budget -> treat as small join.
    const u64 N = static_cast<u64>(pilot_pairs.size());
    out->adaptive_branch = "enumerate_all";
    out->count = MakeExactCount(N);

    out->samples.Clear();
    out->samples.with_replacement = true;
    out->samples.weighted = false;
    out->samples.weights.clear();

    if (N == 0 || t == 0) {
      out->ok = true;
      return true;
    }

    // Sample uniformly with replacement from the stored join list.
    {
      auto _ = out->phases.Scoped("run_small_join_sample_from_list");
      out->samples.pairs.resize(static_cast<usize>(t));
      for (u64 i = 0; i < t; ++i) {
        const u64 idx = rng.UniformU64(N);
        out->samples.pairs[static_cast<usize>(i)] = pilot_pairs[static_cast<usize>(idx)];
      }
    }

    // Validation.
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

  // If we hit the limit, we do not know |J| exactly (it might be > j_star or just > enum_cap).
  out->adaptive_branch = "fallback_sampling";
  if (enum_cap > 0 && pilot_limit == enum_cap) {
    out->enumeration_truncated = true;
    out->note = "pilot enumeration truncated by enum_cap; adaptive decision may be conservative";
  } else {
    out->enumeration_truncated = true;  // truncated relative to full join
  }

  // --------------------------
  // Fallback: baseline-provided Count + Sample (sampling-based)
  // --------------------------
  {
    auto _ = out->phases.Scoped("run_fallback_count");
    if (!baseline->Count(cfg, &rng, &out->count, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  {
    auto _ = out->phases.Scoped("run_fallback_sample");
    if (!baseline->Sample(cfg, &rng, &out->samples, &out->phases, &local_err)) {
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  // Validation.
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
