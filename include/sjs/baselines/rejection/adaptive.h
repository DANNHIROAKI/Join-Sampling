#pragma once
// sjs/baselines/rejection/adaptive.h
//
// AGR-BoxJoin baseline — Variant::Adaptive (AGR-AS)
//
// Doc-aligned adaptive strategy (docs/Baseline/Amagata’25 Baseline.md v2):
//  - Preprocess once (same as AGR-S): build M_i, CellMap, (optional) SlabIndex, μ/alias.
//  - Enumerate true join pairs until either:
//      * stream ends (|J| <= J*): keep AllPairs and return exact |J|
//      * more than J* pairs are seen (|J| > J*): clear AllPairs (avoid prefix bias)
//        and switch to AGR-S for sampling.
//  - Sample():
//      * small join: array uniform sampling with replacement from AllPairs
//      * large join: AGR-S rejection sampling using the already-built alias tables
//
// Empty join handling:
//  - If enumeration finds 0 pairs, return empty and allow Sample() to return empty.

#include "sjs/baselines/rejection/sampling.h"

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
namespace rejection {

template <int Dim, class T = Scalar>
class RejectionAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "RejectionAdaptiveBaseline expects Dim>=2");
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Rejection; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "rejection_adaptive"; }

  void Reset() override {
    st_.Reset();
    mode_ = Mode::Unknown;
    all_pairs_.clear();
    all_pairs_.shrink_to_fit();
    exact_n_ = 0;
  }

  bool Build(const DatasetT& ds,
             const Config& cfg,
             PhaseRecorder* phases,
             std::string* err) override {
    Reset();
    if (!st_.BuildIndex(ds, cfg, phases, err)) return false;
    mode_ = Mode::Unknown;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    if (!st_.built) {
      if (err) *err = "RejectionAdaptiveBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RejectionAdaptiveBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("adaptive_phase1") : PhaseRecorder::ScopedPhase(nullptr, "");

    // Effective J* bound: obey enum_cap if provided.
    u64 J_star = cfg.run.j_star;
    if (cfg.run.enum_cap > 0) J_star = std::min(J_star, cfg.run.enum_cap);

    // Special-case: J_star==0 still enumerates minimally to distinguish empty/non-empty,
    // preventing the "AGR-S on empty join" non-termination scenario.
    const u64 limit = (J_star == 0) ? 1ULL : (J_star + 1ULL);

    all_pairs_.clear();
    exact_n_ = 0;

    auto stream = std::make_unique<detail::RejectionJoinEnumerator<Dim, T>>(&st_);
    PairId p;
    u64 seen = 0;
    bool switch_to_sampling = false;

    while (stream->Next(&p)) {
      ++seen;
      if (seen <= J_star) {
        all_pairs_.push_back(p);
      } else {
        switch_to_sampling = true;
        break;
      }
      if (seen >= limit) {
        if (J_star == 0) switch_to_sampling = true;  // found at least one pair
        break;
      }
    }

    if (!switch_to_sampling) {
      // Fully enumerated (|J| <= J*).
      mode_ = Mode::Enumerate;
      exact_n_ = seen;
      *out = MakeExactCount(static_cast<u64>(exact_n_));
      return true;
    }

    // Switch branch: discard prefix to avoid bias (doc §4.3).
    mode_ = Mode::Sampling;
    all_pairs_.clear();

    // Provide a pilot estimate of |J| using RNG-isolated pilot (doc §5.4).
    if (!rng) {
      *out = MakeEstimateCount(std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              0);
      return true;
    }
    Rng tmp = *rng;
    return st_.EstimateCountByPilot(&tmp, out, phases, err);
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!st_.built) {
      if (err) *err = "RejectionAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "RejectionAdaptiveBaseline::Sample: null rng/out";
      return false;
    }

    auto scoped = phases ? phases->Scoped("adaptive_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    if (mode_ == Mode::Unknown) {
      if (err) *err = "RejectionAdaptiveBaseline::Sample: Count() must be called before Sample()";
      return false;
    }

    if (mode_ == Mode::Enumerate) {
      if (all_pairs_.empty()) return true;  // empty join
      out->pairs.resize(static_cast<usize>(t));
      const u64 N = static_cast<u64>(all_pairs_.size());
      for (u64 i = 0; i < t; ++i) {
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(rng->UniformU64(N))];
      }
      return true;
    }

    // mode_ == Sampling: run AGR-S (join guaranteed non-empty due to switch condition).
    if (st_.mu_sum == 0) return true;  // defensive; implies join empty

    out->pairs.reserve(static_cast<usize>(t));
    while (out->pairs.size() < static_cast<usize>(t)) {
      u32 ai = 0, bi = 0;
      if (!st_.DrawProposal(rng, &ai, &bi)) return true;
      if (!st_.ABox(ai).Intersects(st_.BBox(bi))) continue;
      out->pairs.push_back(st_.MakeOutputPair(ai, bi));
    }
    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!st_.built) {
      if (err) *err = "RejectionAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RejectionJoinEnumerator<Dim, T>>(&st_);
  }

 private:
  enum class Mode : u8 { Unknown = 0, Enumerate = 1, Sampling = 2 };

  detail::RejectionState<Dim, T> st_;
  Mode mode_{Mode::Unknown};

  std::vector<PairId> all_pairs_;  // used when |J| <= J*
  u64 exact_n_{0};
};

}  // namespace rejection
}  // namespace baselines
}  // namespace sjs
