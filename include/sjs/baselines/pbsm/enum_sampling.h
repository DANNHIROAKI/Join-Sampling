#pragma once
// sjs/baselines/pbsm/enum_sampling.h
//
// PBSM (Partition-Based Spatial-Merge) baseline (Variant::EnumSampling).
//
// Baseline: "Enumerate+Sampling" (see Tsitsigkos’19 Baseline.md §4.3 / Theorem 1)
//
//   1) Build PBSM partitions once (multi-assignment + per-partition sorting).
//   2) EnumerateUniquePairs() once and MATERIALIZE the entire join result J into an array.
//   3) Draw t i.i.d. unbiased indices in [0, |J|) WITH replacement and return Pairs[idx].
//
// Notes
// -----
// - This variant intentionally uses O(|J|) memory.
// - Random integers in [0,W) MUST be unbiased (Baseline Appendix A1). We rely on
//   Rng::UniformU64(W) to implement unbiased generation (e.g., rejection sampling).

#include "sjs/baselines/pbsm/sampling.h"  // provides PBSMIndex + enumerator

#include "sjs/core/assert.h"

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
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
    have_enum_ = false;
    W_ = 0;
    std::vector<PairId>().swap(all_pairs_);  // free memory
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    built_ = false;
    have_enum_ = false;
    W_ = 0;
    std::vector<PairId>().swap(all_pairs_);  // free memory from previous runs

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
    (void)rng;  // deterministic; RNG not needed for counting
    if (!built_ || !index_.Built()) {
      if (err) *err = "PBSMEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "PBSMEnumSamplingBaseline::Count: out is null";
      return false;
    }

    // Enumerate+Sampling naturally materializes all pairs; we do that here so a
    // subsequent Sample() can reuse the same one-pass enumeration result.
    if (!EnsureEnumerated_(phases, err)) return false;

    *out = MakeExactCount(W_);
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

    if (!EnsureEnumerated_(phases, err)) return false;

    const u64 W = W_;
    if (W == 0) return true;

    auto scoped = phases ? phases->Scoped("pbsm_enum_sampling_draw") : PhaseRecorder::ScopedPhase(nullptr, "");
    out->pairs.resize(static_cast<usize>(t));
    for (u64 i = 0; i < t; ++i) {
      const u64 idx = rng->UniformU64(W);  // must be unbiased (Baseline Appendix A1)
      out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(idx)];
    }
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
  bool EnsureEnumerated_(PhaseRecorder* phases, std::string* err) {
    (void)err;  // unused parameter (interface requirement)
    if (have_enum_) return true;

    auto scoped = phases ? phases->Scoped("pbsm_enum_sampling_enumerate") : PhaseRecorder::ScopedPhase(nullptr, "");

    std::vector<PairId>().swap(all_pairs_);
    all_pairs_.clear();

    detail::PBSMJoinEnumerator<Dim, T> stream(&index_);
    stream.Reset();

    PairId p;
    while (stream.Next(&p)) {
      all_pairs_.push_back(p);
    }

    W_ = static_cast<u64>(all_pairs_.size());
    have_enum_ = true;
    return true;
  }

  detail::PBSMIndex<Dim, T> index_;
  bool built_ = false;

  // Cached one-pass enumeration result (Pairs) for Enumerate+Sampling.
  bool have_enum_ = false;
  u64 W_ = 0;
  std::vector<PairId> all_pairs_;
};

}  // namespace pbsm
}  // namespace baselines
}  // namespace sjs
