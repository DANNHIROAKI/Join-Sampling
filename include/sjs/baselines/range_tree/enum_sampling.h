#pragma once
// sjs/baselines/range_tree/enum_sampling.h
//
// RangeTree baseline (Variant::EnumSampling) — RT-Enumerate (Enumerate+Sampling)
//
// This variant is the *materialize-all* baseline described in
//   docs/Baseline/RangeTree Baseline v2.0.md  (RT-Enumerate)
//
//
// Algorithm (v2.0 §3.1):
//   1) One plane sweep on axis 0.
//      For each START(q): Report all opposite active rectangles that intersect q
//      in the remaining dimensions (here only y), and append the corresponding
//      join pairs into an array Pairs.
//   2) After the sweep, |Pairs| = |J|. Produce t i.i.d. uniform samples from J
//      by sampling array indices uniformly with replacement.
//
// Notes
// -----
// * This baseline intentionally has time/space dependency on |J|, and may be
//   infeasible when the join result is huge. It serves as a correctness and
//   constant-factor reference.
//
// Implementation detail:
// * We reuse the deterministic RangeJoinEnumerator (same sweep + ActiveRangeTree
//   machinery as the Sampling baseline) to materialize all pairs.

#include "sjs/baselines/range_tree/sampling.h"

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace range_tree {

template <int Dim, class T = Scalar>
class RangeTreeEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim == 2, "RangeTreeEnumSamplingBaseline is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::RangeTree; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "range_tree_enum+sampling"; }

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

  // Exact |J| by full enumeration (consistent with RT-Enumerate).
  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;  // deterministic
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("enumerate_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    auto stream = std::make_unique<detail::RangeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                       join::SideTieBreak::RBeforeS);

    PairId p;
    u64 W = 0;
    while (stream->Next(&p)) {
      if (W == std::numeric_limits<u64>::max()) {
        if (err) *err = "RangeTreeEnumSamplingBaseline::Count: |J| overflowed u64";
        return false;
      }
      ++W;
    }

    *out = MakeExactCount(W);
    return true;
  }

  // Materialize all pairs then sample indices uniformly with replacement.
  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t64 = cfg.run.t;
    if (t64 == 0) return true;
    if (t64 > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: t too large for u32 slots";
      return false;
    }
    const u32 t = static_cast<u32>(t64);

    auto scoped = phases ? phases->Scoped("enumerate_materialize") : PhaseRecorder::ScopedPhase(nullptr, "");

    auto stream = Enumerate(cfg, phases, err);
    if (!stream) return false;

    std::vector<PairId> pairs;
    PairId p;
    while (stream->Next(&p)) {
      pairs.push_back(p);
    }

    const u64 W = static_cast<u64>(pairs.size());
    if (W == 0) {
      return true;  // empty join
    }

    out->pairs.resize(static_cast<usize>(t));
    for (u32 i = 0; i < t; ++i) {
      const u64 j = rng->UniformU64(W);
      out->pairs[static_cast<usize>(i)] = pairs[static_cast<usize>(j)];
    }
    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RangeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                 join::SideTieBreak::RBeforeS);
  }

 private:
  const DatasetT* ds_ = nullptr;
  bool built_ = false;
};

}  // namespace range_tree
}  // namespace baselines
}  // namespace sjs
