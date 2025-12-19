#pragma once
// sjs/baselines/aabb/enum_sampling.h
//
// Plane Sweep + Dynamic AABB-Tree baseline (Variant::EnumSampling).
//
// Baseline v2.0 design (Enumerate+Sampling):
//   1) One plane sweep to enumerate *all* join pairs J into a vector Pairs.
//   2) Uniform i.i.d. sampling with replacement by sampling indices in Pairs.
//
// This is the "materialize then sample" baseline described in
// "AABB-Tree Baseline v2.0.md". It is intentionally memory-heavy
// (space Θ(|J|)) and serves as a correctness/quality reference.

#include "sjs/baselines/aabb/sampling.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace aabb {

template <int Dim, class T>
class AABBEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "AABBEnumSamplingBaseline requires Dim >= 2");
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::AABB; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "aabb_enum+sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)err;
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = true;
    // Empty relations are fine; enumeration will be empty.
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic
    if (!built_ || !ds_) {
      if (err) *err = "AABBEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "AABBEnumSamplingBaseline::Count: out is null";
      return false;
    }

    // Efficient exact count using the Sampling variant's Phase-1 scheme.
    auto scoped = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");
    AABBSamplingBaseline<Dim, T> counter;
    if (!counter.Build(*ds_, /*cfg=*/cfg, phases, err)) return false;
    return counter.Count(cfg, /*rng=*/nullptr, out, phases, err);
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "AABBEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "AABBEnumSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "AABBEnumSamplingBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u32 t = static_cast<u32>(cfg.run.t);
    if (t == 0) return true;

    // --------------------------
    // Phase 1: enumerate all join pairs into Pairs
    // --------------------------
    std::vector<PairId> pairs;
    {
      auto scoped = phases ? phases->Scoped("enumerate_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");

      auto stream = Enumerate(cfg, phases, err);
      if (!stream) return false;

      PairId p;
      while (stream->Next(&p)) {
        pairs.push_back(p);
      }
    }

    const u64 W = static_cast<u64>(pairs.size());
    if (W == 0) {
      // Empty join.
      return true;
    }

    // --------------------------
    // Phase 2: uniform i.i.d. sampling from the materialized array
    // --------------------------
    {
      auto scoped = phases ? phases->Scoped("sample_from_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");
      out->pairs.resize(static_cast<usize>(t));
      for (u32 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(W);
        out->pairs[static_cast<usize>(i)] = pairs[static_cast<usize>(idx)];
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !ds_) {
      if (err) *err = "AABBEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::AABBJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                join::SideTieBreak::RBeforeS);
  }

 private:
  const DatasetT* ds_ = nullptr;
  bool built_ = false;
};

}  // namespace aabb
}  // namespace baselines
}  // namespace sjs
