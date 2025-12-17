#pragma once
// sjs/baselines/tlsop/adaptive.h
//
// TLSOP baseline (Two-layer SOP) (Variant::Adaptive).
//
// This follows the adaptive strategy described in the uploaded "Tsitsigkos’23.md":
//   - Pass 1: enumerate the TLSOP JoinStream once, counting |J| exactly.
//     While the running count W stays <= J*, also materialize the pairs.
//     Once W exceeds J*, immediately release the materialized buffer and continue
//     counting only.
//   - If |J| <= J*: sample directly from the materialized list (1 pass total).
//   - Else: perform a second enumeration pass and select by sorted i.i.d. ranks in
//     [0,|J|) (2 passes total).
//
// Compatibility: This baseline also provides Enumerate() so it can be used with the
// project-wide adaptive runner (which may do its own pilot enumeration). Even if the
// runner is used, the Count/Sample logic here remains correct (may incur extra passes).

#include "sjs/baselines/tlsop/sampling.h"  // TLSOPIndex + TLSOPJoinEnumerator

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
namespace tlsop {

template <int Dim, class T = Scalar>
class TLSOPAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim == 2, "TLSOPAdaptiveBaseline is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::TLSOP; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "tlsop_adaptive"; }

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
      if (err) *err = "TLSOPAdaptiveBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "TLSOPAdaptiveBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("tlsop_adaptive_pass1") : PhaseRecorder::ScopedPhase(nullptr, "");

    // As in pbsm adaptive, clamp J* by enum_cap if user configured a cap for enumeration-based runs.
    const u64 J_star_cfg = cfg.run.j_star;
    const u64 cap_cfg = cfg.run.enum_cap;
    const u64 J_star = (cap_cfg == 0) ? J_star_cfg : std::min(J_star_cfg, cap_cfg);

    detail::TLSOPJoinEnumerator<Dim, T> stream(&index_);
    stream.Reset();

    W_ = 0;
    have_count_ = false;
    mode_ = Mode::Enumerate;
    all_pairs_.clear();

    if (J_star > 0) {
      const u64 reserve_cap = std::min<u64>(J_star, 1'000'000ULL);
      all_pairs_.reserve(static_cast<usize>(reserve_cap));
    }

    PairId p;
    while (stream.Next(&p)) {
      ++W_;
      if (mode_ == Mode::Enumerate) {
        if (W_ <= J_star) {
          all_pairs_.push_back(p);
        } else {
          mode_ = Mode::CountOnly;
          all_pairs_.clear();
          all_pairs_.shrink_to_fit();
        }
      }
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
      if (err) *err = "TLSOPAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "TLSOPAdaptiveBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "TLSOPAdaptiveBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->weighted = false;
    out->with_replacement = true;

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    if (!have_count_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    const u64 W = W_;
    if (W == 0) return true;

    // Small-join branch: sample directly from the materialized list.
    if (mode_ == Mode::Enumerate && all_pairs_.size() == static_cast<usize>(W)) {
      auto scoped = phases ? phases->Scoped("tlsop_adaptive_sample_small") : PhaseRecorder::ScopedPhase(nullptr, "");
      out->pairs.resize(static_cast<usize>(t));
      for (u64 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(W);
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(idx)];
      }
      return true;
    }

    // Large-join branch: second pass rank selection over the deterministic join stream.
    auto scoped = phases ? phases->Scoped("tlsop_adaptive_sample_large") : PhaseRecorder::ScopedPhase(nullptr, "");

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

    detail::TLSOPJoinEnumerator<Dim, T> stream(&index_);
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
      if (err) *err = "TLSOPAdaptiveBaseline::Sample: stream ended early in pass2 (non-deterministic enumerate?)";
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
      if (err) *err = "TLSOPAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::TLSOPJoinEnumerator<Dim, T>>(&index_);
  }

 private:
  enum class Mode : u8 {
    Unknown = 0,
    Enumerate = 1,
    CountOnly = 2,
  };

  detail::TLSOPIndex<Dim, T> index_;
  bool built_ = false;

  bool have_count_ = false;
  u64 W_ = 0;
  Mode mode_ = Mode::Unknown;
  std::vector<PairId> all_pairs_;
};

}  // namespace tlsop
}  // namespace baselines
}  // namespace sjs
