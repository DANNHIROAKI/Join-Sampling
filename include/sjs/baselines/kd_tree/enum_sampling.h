#pragma once
// sjs/baselines/kd_tree/enum_sampling.h
//
// KD-tree baseline (Variant::EnumSampling).
//
// Baseline v2.0 design: Enumerate+Sampling
//   1) One sweep enumerates all join pairs J deterministically (using the same
//      sweep + global rank embedding + static KD-tree + active on/off machinery
//      as the Sampling baseline).
//   2) Draw t i.i.d. uniform samples with replacement by random indexing into
//      the materialized pair array.
//
// Notes:
//   - This variant is intended for regimes where |J| is small enough to fit
//     in memory. When |J| is large, use Variant::Sampling or Variant::Adaptive.
//

#include "sjs/baselines/kd_tree/sampling.h"

#include <cstdint>
#include <limits>
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

  // Exact count by full enumeration (|J| = number of enumerated pairs).
  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;  // deterministic
    if (!built_ || !ds_) {
      if (err) *err = "KDTreeEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "KDTreeEnumSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("enumerate_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    auto stream = Enumerate(cfg, phases, err);
    if (!stream) return false;

    u64 cnt = 0;
    PairId tmp;
    while (stream->Next(&tmp)) {
      if (cnt == std::numeric_limits<u64>::max()) {
        if (err) *err = "KDTreeEnumSamplingBaseline::Count: |J| overflowed u64";
        return false;
      }
      ++cnt;
    }

    *out = MakeExactCount(cnt);
    return true;
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

    const u64 t64 = cfg.run.t;
    if (t64 == 0) return true;
    if (t64 > static_cast<u64>(std::numeric_limits<usize>::max())) {
      if (err) *err = "KDTreeEnumSamplingBaseline::Sample: t too large";
      return false;
    }
    const usize t = static_cast<usize>(t64);

    // Step 1: enumerate all join pairs.
    auto scoped = phases ? phases->Scoped("enumerate_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");

    auto stream = Enumerate(cfg, phases, err);
    if (!stream) return false;

    std::vector<PairId> pairs;
    pairs.reserve(1024);  // heuristic; grows as needed

    PairId p;
    while (stream->Next(&p)) {
      pairs.push_back(p);
    }

    const u64 N = static_cast<u64>(pairs.size());
    if (N == 0) {
      // Empty join -> empty sample set.
      return true;
    }

    // Step 2: i.i.d. uniform sampling with replacement via random indexing.
    out->pairs.resize(t);
    for (usize i = 0; i < t; ++i) {
      const u64 idx = rng->UniformU64(N);  // in [0, N)
      out->pairs[i] = pairs[static_cast<usize>(idx)];
    }

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
