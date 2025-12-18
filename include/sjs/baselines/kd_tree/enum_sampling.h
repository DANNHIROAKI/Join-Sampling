#pragma once
// sjs/baselines/kd_tree/enum_sampling.h
//
// KD-tree baseline (Variant::EnumSampling): enumerate join pairs deterministically
// using sweep + KD-tree reporting, then draw i.i.d. uniform samples via two-pass
// rank sampling.
//
// This is a faithful implementation of the Enum+Sampling protocol used for baselines
// in the project: the runner owns the two-pass logic, but this baseline provides:
//   - Enumerate(): a deterministic join stream (KDJoinEnumerator)
//   - Sample(): convenience wrapper around sampling::RankSampleWithReplacement
//   - Count(): convenience wrapper (exact count) via a sampling baseline pass

#include "sjs/baselines/kd_tree/sampling.h"
#include "sjs/sampling/rank_sampling.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace kd_tree {

template <int Dim, class T = Scalar>
class KDTreeEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "KDTreeEnumSamplingBaseline requires Dim >= 2");

  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::KDTree; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "kd_tree_enum_sampling"; }

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
    return true;
  }

  // Convenience exact count by reusing the sampling baseline's Phase-1.
  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;
    if (!built_ || !ds_) {
      if (err) *err = "KDTreeEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    auto scoped = phases ? phases->Scoped("count_via_sampling_baseline") : PhaseRecorder::ScopedPhase(nullptr, "");

    KDTreeSamplingBaseline<Dim, T> counter;
    counter.Reset();
    if (!counter.Build(*ds_, cfg, phases, err)) return false;

    // Count is deterministic; rng is unused.
    return counter.Count(cfg, /*rng=*/nullptr, out, phases, err);
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "KDTreeEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "KDTreeEnumSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "KDTreeEnumSamplingBaseline::Sample: out is null";
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

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");

    if (!built_ || !ds_) {
      if (err) *err = "KDTreeEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    return std::make_unique<detail::KDJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                              join::SideTieBreak::RBeforeS);
  }

 private:
  const DatasetT* ds_ = nullptr;
  bool built_ = false;
};

}  // namespace kd_tree
}  // namespace baselines
}  // namespace sjs

