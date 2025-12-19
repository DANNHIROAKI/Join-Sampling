#pragma once
// sjs/baselines/interval_tree/enum_sampling.h
//
// IntervalTree Baseline (v2.0) — Variant::EnumSampling
// ----------------------------------------------------
// Matches the "Enumerate+Sampling" baseline in:
//   docs/Baseline/IntervalTree Baseline v2.0.md
//
// Algorithm:
//   1) One x-sweep enumerates all join pairs J into an array Pairs.
//   2) Output t i.i.d. samples by drawing uniform indices in [0,|Pairs|).
//
// Notes:
//   - This is intended as a classic baseline: O(n log n + |J|) time and O(|J|) memory.
//   - The sweep/overlap logic is shared with sampling.h via IntervalTreeJoinEnumerator.

#include "sjs/baselines/interval_tree/sampling.h"

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
  std::string_view Name() const noexcept override { return "interval_tree_enum_sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;

    pairs_.clear();
    pairs_valid_ = false;
    W_ = 0;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)err;  // unused parameter (interface requirement)
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    ds_ = &ds;
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

    if (!built_ || !ds_) {
      if (err) *err = "IntervalTreeEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "IntervalTreeEnumSamplingBaseline::Count: out is null";
      return false;
    }

    // Enumerate and cache Pairs (this is the baseline's main cost).
    if (!pairs_valid_) {
      auto scoped = phases ? phases->Scoped("enumerate_all_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");

      pairs_.clear();
      std::unique_ptr<IJoinEnumerator> it = Enumerate(cfg, phases, err);
      if (!it) return false;

      PairId p;
      while (it->Next(&p)) pairs_.push_back(p);

      W_ = static_cast<u64>(pairs_.size());
      pairs_valid_ = true;
    }

    *out = MakeExactCount(W_);
    return true;
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
    if (!rng) {
      if (err) *err = "IntervalTreeEnumSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "IntervalTreeEnumSamplingBaseline::Sample: out is null";
      return false;
    }

    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    if (t == 0) return true;

    // Ensure Pairs exist.
    if (!pairs_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }
    if (W_ == 0) return true;

    out->pairs.resize(static_cast<usize>(t));

    auto scoped = phases ? phases->Scoped("sample_from_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");

    for (u32 i = 0; i < t; ++i) {
      const u64 idx = rng->UniformU64(W_);
      out->pairs[static_cast<usize>(i)] = pairs_[static_cast<usize>(idx)];
    }
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

    // Streaming enumerator (one sweep) — used for materializing Pairs and also for external callers.
    return std::make_unique<detail::IntervalTreeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S);
  }

 private:
  const DatasetT* ds_{nullptr};
  bool built_{false};

  std::vector<PairId> pairs_;
  bool pairs_valid_{false};
  u64 W_{0};
};

}  // namespace interval_tree
}  // namespace baselines
}  // namespace sjs