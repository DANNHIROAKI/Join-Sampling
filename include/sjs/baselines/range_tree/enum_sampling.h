#pragma once
// sjs/baselines/range_tree/enum_sampling.h
//
// Plane Sweep + Dynamic Range-Tree baseline (Variant::EnumSampling).
//
// This variant exposes a deterministic enumerator (sweep + ActiveRangeTree)
// and implements Sample() via two-pass rank sampling over that enumerator.
//
// In the experimental harness, EnumSampling variants are typically executed via
// baselines::runners::RunEnumSamplingOnce(), which directly calls Enumerate()
// and performs rank sampling itself. Still, we provide Sample() here for
// completeness.

#include "sjs/baselines/range_tree/sampling.h"
#include "sjs/sampling/rank_sampling.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace range_tree {

template <int Dim, class T>
class RangeTreeEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim == 2, "RangeTreeEnumSamplingBaseline is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::RangeTree; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "range_tree_enum+sampling"; }

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
    (void)rng;
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Count: out is null";
      return false;
    }

    // Exact count using the Phase-1 scheme from the Sampling variant.
    auto scoped = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");
    RangeTreeSamplingBaseline<Dim, T> counter;
    if (!counter.Build(*ds_, cfg, phases, err)) return false;
    return counter.Count(cfg, /*rng=*/nullptr, out, phases, err);
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: out is null";
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
    out->pairs = std::move(samples);
    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RangeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                 join::SideTieBreak::RBeforeS);
  }

 private:
  const DatasetT* ds_ = nullptr;
  bool built_ = false;
};

}  // namespace range_tree
}  // namespace baselines
}  // namespace sjs
