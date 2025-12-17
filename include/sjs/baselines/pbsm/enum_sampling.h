#pragma once
// sjs/baselines/pbsm/enum_sampling.h
//
// PBSM (Partition-Based Spatial-Merge) baseline (Variant::EnumSampling).
//
// This variant exposes the deterministic PBSM join enumerator and performs
// uniform sampling WITH replacement by two-pass rank sampling over that
// enumerator.
//
// In the experiment harness, EnumSampling variants are typically executed via
// baselines::runners::RunEnumSamplingOnce(), which directly calls Enumerate()
// and applies rank sampling externally. We still provide Sample() here for
// completeness and for standalone use.

#include "sjs/baselines/pbsm/sampling.h"  // provides PBSMIndex + enumerator

#include "sjs/sampling/rank_sampling.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace pbsm {

template <int Dim, class T = Scalar>
class PBSMEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "PBSMEnumSamplingBaseline requires Dim>=2");
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::PBSM; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "pbsm_enum+sampling"; }

  void Reset() override {
    index_.Reset();
    built_ = false;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    built_ = false;
    if (!index_.Build(ds, cfg, phases, err)) return false;
    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;
    if (!built_ || !index_.Built()) {
      if (err) *err = "PBSMEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "PBSMEnumSamplingBaseline::Count: out is null";
      return false;
    }
    auto scoped = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");

    detail::PBSMJoinEnumerator<Dim, T> stream(&index_);
    u64 W = 0;
    PairId tmp;
    while (stream.Next(&tmp)) ++W;
    *out = MakeExactCount(W);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !index_.Built()) {
      if (err) *err = "PBSMEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "PBSMEnumSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "PBSMEnumSamplingBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    auto scoped = phases ? phases->Scoped("rank_sampling") : PhaseRecorder::ScopedPhase(nullptr, "");
    detail::PBSMJoinEnumerator<Dim, T> stream(&index_);

    std::vector<PairId> samples;
    sampling::RankSamplingInfo info;
    if (!sampling::RankSampleWithReplacement<detail::PBSMJoinEnumerator<Dim, T>, PairId>(
            &stream, t, rng, &samples, &info, err)) {
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
    if (!built_ || !index_.Built()) {
      if (err) *err = "PBSMEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::PBSMJoinEnumerator<Dim, T>>(&index_);
  }

 private:
  detail::PBSMIndex<Dim, T> index_;
  bool built_ = false;
};

}  // namespace pbsm
}  // namespace baselines
}  // namespace sjs
