#pragma once
// sjs/baselines/tsunami/enum_sampling.h
//
// Tsunami baseline (Variant::EnumSampling).
//
// This variant implements TSUNAMI-Enumerate-ArraySample per the baseline writeup:
//  1. Enumerate all join pairs into AllPairs array
//  2. Sample t pairs uniformly with replacement from AllPairs
//
// This is the naive but straightforward approach: enumerate first, then sample.
// It has good constant factors for small joins but requires O(|J|) memory.

#include "sjs/baselines/tsunami/sampling.h"

#include <algorithm>
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
    counted_ = false;
    all_pairs_.clear();
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

    // Enumerate all pairs to get exact count and populate all_pairs_.
    auto scoped = phases ? phases->Scoped("enumerate_all_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");
    
    all_pairs_.clear();
    auto stream = Enumerate(cfg, phases, err);
    if (!stream) return false;

    stream->Reset();
    PairId tmp{};
    while (stream->Next(&tmp)) {
      all_pairs_.push_back(tmp);
    }
    
    counted_ = true;
    *out = MakeExactCount(static_cast<u64>(all_pairs_.size()));
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

    // Ensure all pairs are enumerated.
    if (!counted_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    if (all_pairs_.empty()) {
      // Empty join.
      return true;
    }

    // Sample t pairs uniformly with replacement from AllPairs.
    auto scoped = phases ? phases->Scoped("sample_from_all_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");
    const u64 W = static_cast<u64>(all_pairs_.size());
    
    out->pairs.resize(static_cast<usize>(t));
    for (u64 i = 0; i < t; ++i) {
      const u64 idx = rng->UniformU64(W);
      out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(idx)];
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
  
  bool counted_ = false;
  std::vector<PairId> all_pairs_;  // Explicit storage of all join pairs
};

}  // namespace tsunami
}  // namespace baselines
}  // namespace sjs
