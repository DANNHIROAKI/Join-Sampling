#pragma once
// sjs/baselines/ours/enum_sampling.h
//
// Our method (Enumerate+Sampling variant).
//
// In the unified experiment protocol, Enum+Sampling baselines expose a
// deterministic join enumerator; sampling is then performed by a generic
// two-pass rank-sampling routine (see sjs/sampling/rank_sampling.h) to obtain
// i.i.d. uniform samples with replacement.
//
// This baseline intentionally keeps the enumerator simple and robust by using
// the existing deterministic plane-sweep join stream.

#include "sjs/baselines/ours/sampling.h"  // reuse PlaneSweepEnumeratorWrapper
#include "sjs/sampling/rank_sampling.h"

#include <memory>
#include <string>
#include <string_view>

namespace sjs {
namespace baselines {
namespace ours {

template <int Dim, class T = Scalar>
class OursEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "ours_enum+sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
  }

  bool Build(const DatasetT& ds,
             const Config& cfg,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    auto ph = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = true;

    if constexpr (Dim != 2) {
      // The project currently only runs 2D. We keep a friendly error rather than static_assert.
      if (err) *err = "OursEnumSamplingBaseline: currently only Dim=2 is supported";
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
    if (!built_ || ds_ == nullptr) {
      if (err) *err = "OursEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    auto ph = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");

    join::PlaneSweepOptions opt;
    opt.axis = 0;
    opt.side_order = join::SideTieBreak::RBeforeS;
    opt.skip_axis_check = true;

    const u64 n = join::CountPlaneSweep<Dim, T>(ds_->R, ds_->S, opt);
    if (out) *out = MakeExactCount(n);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || ds_ == nullptr) {
      if (err) *err = "OursEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "OursEnumSamplingBaseline::Sample: null rng/out";
      return false;
    }

    auto ph = phases ? phases->Scoped("sample") : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 t = cfg.run.t;
    out->Clear();

    // Build a deterministic enumerator and perform two-pass rank sampling.
    std::string enum_err;
    auto enumerator = Enumerate(cfg, nullptr, &enum_err);
    if (!enumerator) {
      if (err) *err = "OursEnumSamplingBaseline::Sample: Enumerate failed: " + enum_err;
      return false;
    }

    sampling::RankSamplingInfo info;
    std::vector<PairId> samples;
    if (!sampling::RankSampleWithReplacement<IJoinEnumerator, PairId>(enumerator.get(), t, rng, &samples, &info, err)) {
      return false;
    }

    out->pairs = std::move(samples);
    out->weighted = false;
    out->with_replacement = true;
    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || ds_ == nullptr) {
      if (err) *err = "OursEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    join::PlaneSweepOptions opt;
    opt.axis = 0;
    opt.side_order = join::SideTieBreak::RBeforeS;
    opt.skip_axis_check = true;

    return std::make_unique<detail::PlaneSweepEnumeratorWrapper<Dim, T>>(ds_->R, ds_->S, opt);
  }

 private:
  const DatasetT* ds_{nullptr};
  bool built_{false};
};

}  // namespace ours
}  // namespace baselines
}  // namespace sjs
