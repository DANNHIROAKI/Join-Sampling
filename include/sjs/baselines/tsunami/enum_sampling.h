#pragma once
// sjs/baselines/tsunami/enum_sampling.h
//
// Tsunami baseline (Variant::EnumSampling).
//
// This variant exposes a deterministic enumerator (see sampling.h) and relies
// on two-pass rank sampling over that stream (implemented by the runner).
//
// The baseline also provides a convenience Sample() method using the generic
// rank-sampling helper.

#include "sjs/baselines/tsunami/sampling.h"

#include "sjs/sampling/rank_sampling.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace tsunami {

template <int Dim, class T = Scalar>
class TsunamiEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 1, "TsunamiEnumSamplingBaseline requires Dim >= 1");
  static_assert(2 * Dim <= kMaxSupportedDim,
                "TsunamiEnumSamplingBaseline uses Point<2*Dim>; increase kMaxSupportedDim or reduce Dim");

  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::Tsunami; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "tsunami_enum_sampling"; }

  void Reset() override {
    prep_.Reset();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    Reset();
    return prep_.Build(ds, cfg, phases, err);
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;
    if (!out) {
      if (err) *err = "TsunamiEnumSamplingBaseline::Count: out is null";
      return false;
    }
    *out = MakeExactCount(0);

    if (!prep_.built()) {
      if (err) *err = "TsunamiEnumSamplingBaseline::Count: call Build() first";
      return false;
    }

    // Exact count by enumerating (stream length).
    auto scoped = phases ? phases->Scoped("count_by_enumeration") : PhaseRecorder::ScopedPhase(nullptr, "");
    auto stream = Enumerate(cfg, phases, err);
    if (!stream) return false;

    stream->Reset();
    PairId tmp{};
    u64 n = 0;
    while (stream->Next(&tmp)) {
      ++n;
    }
    *out = MakeExactCount(n);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!prep_.built()) {
      if (err) *err = "TsunamiEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "TsunamiEnumSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "TsunamiEnumSamplingBaseline::Sample: out is null";
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
    if (!prep_.built()) {
      if (err) *err = "TsunamiEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::TsunamiJoinEnumerator<Dim, T>>(&prep_);
  }

 private:
  detail::TsunamiPreproc<Dim, T> prep_;
};

}  // namespace tsunami
}  // namespace baselines
}  // namespace sjs
