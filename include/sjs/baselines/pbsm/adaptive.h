#pragma once
// sjs/baselines/pbsm/adaptive.h
//
// PBSM (Partition-Based Spatial-Merge) baseline (Variant::Adaptive).
//
// This implements Baseline "Adaptive+Sampling" (see Tsitsigkos’19 Baseline.md §4.5):
//   - Phase 1: enumerate the unique PBSM join stream once, counting |J| exactly.
//     While the running count W stays <= J*, also materialize the pairs.
//     If W exceeds J*, discard the materialized pairs and continue counting only.
//   - If |J| <= J*: sample directly from the materialized list (one-pass total).
//   - Else: perform a second PBSM enumeration pass and select pairs by sorted
//     i.i.d. ranks in [0,|J|) (two-pass total).
//
// This header is compatible with the project-wide adaptive runner as well:
// the runner may do pilot enumeration and then call Count()/Sample() for the
// fallback branch. The logic here still produces correct uniform samples.

#include "sjs/baselines/pbsm/sampling.h"  // provides PBSMIndex + enumerator

#include "sjs/core/assert.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace pbsm {

template <int Dim, class T = Scalar>
class PBSMAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "PBSMAdaptiveBaseline requires Dim>=2");
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::PBSM; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "pbsm_adaptive"; }

  void Reset() override {
    index_.Reset();
    built_ = false;
    have_count_ = false;
    W_ = 0;
    mode_ = Mode::Unknown;
    all_pairs_.clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    built_ = false;
    have_count_ = false;
    W_ = 0;
    mode_ = Mode::Unknown;
    all_pairs_.clear();
    if (!index_.Build(ds, cfg, phases, err)) return false;
    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic
    if (!built_ || !index_.Built()) {
      if (err) *err = "PBSMAdaptiveBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "PBSMAdaptiveBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("pbsm_adaptive_pass1") : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 J_star_cfg = cfg.run.j_star;
    const u64 cap_cfg = cfg.run.enum_cap;
    const u64 J_star = (cap_cfg == 0) ? J_star_cfg : std::min(J_star_cfg, cap_cfg);

    detail::PBSMJoinEnumerator<Dim, T> stream(&index_);
    stream.Reset();

    W_ = 0;
    have_count_ = false;
    mode_ = Mode::Enumerate;
    all_pairs_.clear();
    if (J_star > 0) {
      // Reserve a bit to reduce reallocs but don't go crazy.
      const u64 cap = std::min<u64>(J_star, 1'000'000ULL);
      all_pairs_.reserve(static_cast<usize>(cap));
    }

    PairId p;
    while (stream.Next(&p)) {
      ++W_;
      if (mode_ == Mode::Enumerate) {
        if (W_ <= J_star) {
          all_pairs_.push_back(p);
        } else {
          // Switch to count-only mode; discard materialization.
          mode_ = Mode::CountOnly;
          // Baseline: discard partial enumeration and free memory (Pairs = empty_and_free()).
          std::vector<PairId>().swap(all_pairs_);
        }
      }
    }

    // If we never exceeded J*, we enumerated all join pairs.
    if (mode_ == Mode::Enumerate) {
      // Sanity: all_pairs_ should contain exactly W_ pairs.
      // (We don't assert hard; keep experiment robust.)
    }

    have_count_ = true;
    *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !index_.Built()) {
      if (err) *err = "PBSMAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "PBSMAdaptiveBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "PBSMAdaptiveBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->weighted = false;
    out->with_replacement = true;

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    if (!have_count_) {
      CountResult cr;
      if (!Count(cfg, /*rng=*/nullptr, &cr, phases, err)) return false;
    }
    const u64 W = W_;
    if (W == 0) return true;

    // Small-join branch: direct sampling from materialized pairs.
    if (mode_ == Mode::Enumerate && all_pairs_.size() == static_cast<usize>(W)) {
      auto scoped = phases ? phases->Scoped("pbsm_adaptive_sample_small") : PhaseRecorder::ScopedPhase(nullptr, "");
      out->pairs.resize(static_cast<usize>(t));
      for (u64 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(W);
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(idx)];
      }
      return true;
    }

    // Large-join branch: pass2 rank selection over the deterministic PBSM stream.
    auto scoped = phases ? phases->Scoped("pbsm_adaptive_sample_large") : PhaseRecorder::ScopedPhase(nullptr, "");

    struct RankReq {
      u64 rank;
      u64 slot;
    };
    std::vector<RankReq> req;
    req.resize(static_cast<usize>(t));
    for (u64 i = 0; i < t; ++i) {
      req[static_cast<usize>(i)] = RankReq{rng->UniformU64(W), i};
    }
    std::sort(req.begin(), req.end(), [](const RankReq& a, const RankReq& b) {
      if (a.rank < b.rank) return true;
      if (b.rank < a.rank) return false;
      return a.slot < b.slot;
    });

    out->pairs.assign(static_cast<usize>(t), PairId{});

    detail::PBSMJoinEnumerator<Dim, T> stream(&index_);
    stream.Reset();
    PairId p;
    u64 cur = 0;
    usize k = 0;
    while (stream.Next(&p)) {
      while (k < req.size() && req[k].rank == cur) {
        out->pairs[static_cast<usize>(req[k].slot)] = p;
        ++k;
      }
      if (k == req.size()) break;
      ++cur;
    }
    if (k != req.size()) {
      if (err) *err = "PBSMAdaptiveBaseline::Sample: stream ended early in pass2 (cardinality mismatch: expected W pairs from EnumerateUniquePairs)";
      return false;
    }
    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !index_.Built()) {
      if (err) *err = "PBSMAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::PBSMJoinEnumerator<Dim, T>>(&index_);
  }

 private:
  enum class Mode : u8 {
    Unknown = 0,
    Enumerate = 1,
    CountOnly = 2,
  };

  detail::PBSMIndex<Dim, T> index_;
  bool built_ = false;

  bool have_count_ = false;
  u64 W_ = 0;
  Mode mode_ = Mode::Unknown;
  std::vector<PairId> all_pairs_;
};

}  // namespace pbsm
}  // namespace baselines
}  // namespace sjs
