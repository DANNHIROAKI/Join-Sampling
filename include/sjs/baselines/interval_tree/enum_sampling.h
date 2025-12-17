#pragma once
// sjs/baselines/interval_tree/enum_sampling.h
//
// Plane Sweep + IntervalTreeSampler baseline (Variant::EnumSampling).
//
// This variant exposes a deterministic enumerator based on the same sweep +
// interval-tree (A/B) machinery. Sampling is performed via generic two-pass
// rank sampling over the enumeration stream (see sjs/sampling/rank_sampling.h).
//
// In the experimental harness, EnumSampling baselines are usually executed via
// baselines::runners::RunEnumSamplingOnce(), which directly calls Enumerate().
// We still provide Count() and Sample() here for completeness.

#include "sjs/baselines/interval_tree/sampling.h"
#include "sjs/sampling/rank_sampling.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace interval_tree {

template <int Dim, class T = Scalar>
class IntervalTreeEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::IntervalTree; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "interval_tree_enum+sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = true;

    if constexpr (Dim != 2) {
      if (err) *err = "IntervalTreeEnumSamplingBaseline: currently only Dim=2 is supported";
      built_ = false;
      return false;
    }
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;
    if (!built_ || !ds_) {
      if (err) *err = "IntervalTreeEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "IntervalTreeEnumSamplingBaseline::Count: out is null";
      return false;
    }

    // Use the O(n log n) exact counting scheme from the Sampling variant.
    auto scoped = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");
    IntervalTreeSamplingBaseline<Dim, T> counter;
    if (!counter.Build(*ds_, cfg, phases, err)) return false;
    return counter.Count(cfg, /*rng=*/nullptr, out, phases, err);
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "IntervalTreeEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "IntervalTreeEnumSamplingBaseline::Sample: null rng/out";
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
      if (err) *err = "IntervalTreeEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    // Self-contained enumerator using interval-tree A/B reporting.
    return std::make_unique<detail::IntervalTreeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                        join::SideTieBreak::RBeforeS);
  }

 private:
  const DatasetT* ds_{nullptr};
  bool built_{false};
};

}  // namespace interval_tree
}  // namespace baselines
}  // namespace sjs
