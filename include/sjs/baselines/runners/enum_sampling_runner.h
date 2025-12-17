#pragma once
// sjs/baselines/runners/enum_sampling_runner.h
//
// Runner for Variant::EnumSampling (Enumerate + Sampling).
//
// Protocol (consistent across baselines):
//   Reset -> Build -> Enumerate(pass1 count) -> Draw ranks -> Enumerate(pass2 select)
//
// Sampling scheme:
//   - Uniform sampling WITH replacement from the enumerated join stream.
//   - Two-pass rank sampling:
//       pass1: count N (or cap)
//       pass2: re-scan stream and pick elements at sampled ranks
//
// Notes:
//   - Enumeration order must be deterministic across Reset() calls.
//   - enum_cap (cfg.run.enum_cap) truncates enumeration if >0.

#include "sjs/baselines/baseline_api.h"
#include "sjs/core/assert.h"
#include "sjs/core/logging.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {

namespace detail {

struct RankReq {
  u64 rank = 0;  // in [0, N)
  u64 slot = 0;  // which output position to fill
};

inline bool RankReqLess(const RankReq& a, const RankReq& b) noexcept {
  if (a.rank < b.rank) return true;
  if (b.rank < a.rank) return false;
  return a.slot < b.slot;
}

}  // namespace detail

template <int Dim, class T = Scalar>
bool RunEnumSamplingOnce(IBaseline<Dim, T>* baseline,
                         const Dataset<Dim, T>& dataset,
                         const Config& cfg,
                         u64 seed,
                         RunReport* out,
                         std::string* err = nullptr) {
  if (!baseline) {
    if (err) *err = "RunEnumSamplingOnce: baseline is null";
    return false;
  }
  if (!out) {
    if (err) *err = "RunEnumSamplingOnce: out is null";
    return false;
  }

  out->ok = false;
  out->error.clear();
  out->method = baseline->method();
  out->variant = Variant::EnumSampling;
  out->baseline_name = std::string(baseline->Name());
  out->dataset_name = dataset.name;
  out->seed = seed;
  out->t = cfg.run.t;
  out->count = CountResult{};
  out->samples.Clear();
  out->used_enumeration = true;
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

  std::unique_ptr<IJoinEnumerator> stream;
  {
    auto _ = out->phases.Scoped("run_enum_prepare");
    stream = baseline->Enumerate(cfg, &out->phases, &local_err);
    if (!stream) {
      if (local_err.empty()) local_err = "baseline->Enumerate returned null stream";
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }
  }

  // --------------------------
  // Pass 1: count stream length (possibly capped)
  // --------------------------
  u64 N = 0;
  const u64 cap = cfg.run.enum_cap;
  {
    auto _ = out->phases.Scoped("run_enum_pass1_count");
    stream->Reset();
    PairId p;
    while (stream->Next(&p)) {
      ++N;
      if (cap > 0 && N >= cap) {
        out->enumeration_truncated = true;
        break;
      }
    }
  }
  out->enumeration_pairs_pass1 = N;
  out->enum_stats_pass1 = stream->Stats();

  // For Enum+Sampling, the count is exact if not truncated.
  if (!out->enumeration_truncated) {
    out->count = MakeExactCount(N);
  } else {
    out->count = MakeEstimateCount(static_cast<long double>(N));
    out->note = "enumeration truncated by enum_cap; count/sample are for prefix only";
  }

  // If N==0, there is nothing to sample.
  const u64 t = cfg.run.t;
  if (N == 0 || t == 0) {
    out->ok = true;
    return true;
  }

  // --------------------------
  // Draw ranks with replacement from [0, N)
  // --------------------------
  std::vector<detail::RankReq> req;
  req.resize(static_cast<usize>(t));
  {
    auto _ = out->phases.Scoped("run_draw_ranks");
    for (u64 i = 0; i < t; ++i) {
      req[static_cast<usize>(i)].rank = rng.UniformU64(N);
      req[static_cast<usize>(i)].slot = i;
    }
    std::sort(req.begin(), req.end(), detail::RankReqLess);
  }

  // --------------------------
  // Pass 2: scan stream and pick requested ranks
  // --------------------------
  out->samples.pairs.assign(static_cast<usize>(t), PairId{});
  out->samples.weights.clear();
  out->samples.weighted = false;
  out->samples.with_replacement = true;

  {
    auto _ = out->phases.Scoped("run_enum_pass2_select");
    stream->Reset();
    PairId p;
    u64 idx = 0;
    usize j = 0;

    // We only consider the first N pairs (N may be capped).
    while (idx < N && stream->Next(&p)) {
      while (j < req.size() && req[j].rank == idx) {
        out->samples.pairs[static_cast<usize>(req[j].slot)] = p;
        ++j;
      }
      if (j >= req.size()) {
        // We've collected all ranks; can stop early.
        break;
      }
      ++idx;
    }

    // If we stopped because stream ended early before reaching N, that's a correctness bug.
    if (idx < N && j < req.size()) {
      local_err = "EnumSampling pass2 ended early: enumerator shorter than pass1 count or nondeterministic Reset()";
      out->error = local_err;
      if (err) *err = local_err;
      return false;
    }

    // If last requested rank is near N-1, j will reach req.size(); else early stop.
  }
  out->enumeration_pairs_pass2 = stream->Stats().output_pairs;  // how many pairs were scanned in pass2
  out->enum_stats_pass2 = stream->Stats();

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
