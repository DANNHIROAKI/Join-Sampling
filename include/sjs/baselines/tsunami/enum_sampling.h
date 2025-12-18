#pragma once
// sjs/baselines/tsunami/enum_sampling.h
//
// Tsunami baseline (Variant::EnumSampling) — strictly aligned with docs/Baseline/Tsunami’20 Baseline.md.
//
// This variant implements:
//   TSUNAMI-Enumerate-ArraySample
//
// Algorithm (exact, i.i.d., uniform WITH replacement):
//   - Enumerate the full join result J into an array AllPairs (Θ(|J|) memory).
//   - Sample t indices i.i.d. uniformly from [0..|J|-1] and output AllPairs[idx].
//
// This is the most direct/naïve baseline and is intended to demonstrate the baseline limitation:
// worst-case memory/time still depends on |J|.

#include "sjs/baselines/tsunami/sampling.h"

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
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

  void Reset() override { prep_.Reset(); }

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

    // Exact count by streaming enumeration (no materialization needed for Count()).
    auto scoped = phases ? phases->Scoped("count_by_stream_enumeration") : PhaseRecorder::ScopedPhase(nullptr, "");
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

    // Enumerate full join result into AllPairs (Θ(|J|) memory).
    auto scoped = phases ? phases->Scoped("enumerate_all_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");
    std::vector<PairId> all_pairs;
    {
      auto stream = Enumerate(cfg, phases, err);
      if (!stream) return false;

      stream->Reset();
      PairId tmp{};
      while (stream->Next(&tmp)) {
        all_pairs.push_back(tmp);
      }
    }

    const u64 W = static_cast<u64>(all_pairs.size());
    if (W == 0) {
      return true;  // empty join
    }

    // Array-sample with replacement.
    out->pairs.resize(static_cast<usize>(t));
    for (u64 i = 0; i < t; ++i) {
      const u64 idx = rng->UniformU64(W);  // [0..W-1]
      out->pairs[static_cast<usize>(i)] = all_pairs[static_cast<usize>(idx)];
    }
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
