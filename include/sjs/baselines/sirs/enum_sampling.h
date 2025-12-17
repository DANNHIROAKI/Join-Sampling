#pragma once
// sjs/baselines/sirs/enum_sampling.h
//
// SIRS baseline (Variant::EnumSampling).
//
// This variant exposes a deterministic join enumerator and performs uniform
// sampling WITH replacement by two-pass rank sampling over the stream.
//
// Enumerate+Sampling is intentionally simple and robust:
//   Pass1: count N=|J| by full enumeration.
//   Pass2: rescan and pick the t requested ranks.
//
// See sjs/sampling/rank_sampling.h for details.

#include "sjs/baselines/sirs/sampling.h"  // provides SIRSState + enumerator
#include "sjs/sampling/rank_sampling.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace sirs {

template <int Dim, class T = Scalar>
class SIRSEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::SIRS; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "sirs_enum+sampling"; }

  void Reset() override { st_.Reset(); }

  bool Build(const DatasetT& ds,
             const Config& cfg,
             PhaseRecorder* phases,
             std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    st_.Reset();
    return st_.BuildIndex(ds, cfg, phases, err);
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic
    if (!st_.built) {
      if (err) *err = "SIRSEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "SIRSEnumSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");

    // For convenience, we reuse the same per-outer range counting used by the
    // Sampling variant. (This avoids enumerating the full join just to count.)
    if (!st_.ComputeDegreesAndAlias(cfg, phases, err)) return false;
    *out = MakeExactCount(st_.W);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!st_.built) {
      if (err) *err = "SIRSEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "SIRSEnumSamplingBaseline::Sample: null rng/out";
      return false;
    }

    auto scoped = phases ? phases->Scoped("sample") : PhaseRecorder::ScopedPhase(nullptr, "");

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    // Two-pass rank sampling over a deterministic enumerator.
    std::string enum_err;
    auto stream = Enumerate(cfg, nullptr, &enum_err);
    if (!stream) {
      if (err) *err = "SIRSEnumSamplingBaseline::Sample: Enumerate failed: " + enum_err;
      return false;
    }

    sampling::RankSamplingInfo info;
    std::vector<PairId> samples;
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
    (void)phases;
    if (!st_.built) {
      if (err) *err = "SIRSEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::SIRSJoinEnumerator<Dim, T>>(&st_);
  }

 private:
  detail::SIRSState<Dim, T> st_;
};

}  // namespace sirs
}  // namespace baselines
}  // namespace sjs
