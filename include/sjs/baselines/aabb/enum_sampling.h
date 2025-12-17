#pragma once
// sjs/baselines/aabb/enum_sampling.h
//
// Plane Sweep + Dynamic AABB-Tree baseline (Variant::EnumSampling).
//
// This variant exposes a deterministic enumerator (sweep + AABB-tree) and
// implements Sample() by two-pass rank sampling over that enumerator.
//
// In the experimental harness, EnumSampling variants are typically executed via
// baselines::runners::RunEnumSamplingOnce(), which directly calls Enumerate()
// and performs rank sampling itself. Still, we provide Sample() here for
// completeness.

#include "sjs/baselines/aabb/sampling.h"
#include "sjs/sampling/rank_sampling.h"

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
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = true;
    if (ds.R.Size() == 0 || ds.S.Size() == 0) {
      // Empty relation is fine; enumeration will be empty.
      return true;
    }
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;
    if (!built_ || !ds_) {
      if (err) *err = "AABBEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "AABBEnumSamplingBaseline::Count: out is null";
      return false;
    }

    // Efficient exact count using the same Phase-1 scheme as the Sampling variant.
    // (No randomness, no materialization of pairs.)
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

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    auto scoped = phases ? phases->Scoped("rank_sampling") : PhaseRecorder::ScopedPhase(nullptr, "");

    auto stream = Enumerate(cfg, phases, err);
    if (!stream) return false;

    std::vector<PairId> samples;
    sampling::RankSamplingInfo info;
    if (!sampling::RankSampleWithReplacement<IJoinEnumerator, PairId>(stream.get(), t, rng, &samples, &info, err)) {
      return false;
    }

    // If join is empty, RankSampleWithReplacement returns 0 samples.
    out->pairs = std::move(samples);
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
