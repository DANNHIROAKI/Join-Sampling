#pragma once
// sjs/baselines/tsunami/adaptive.h
//
// Tsunami baseline (Variant::Adaptive).
//
// Baseline-level adaptive strategy (per Tsunami’20 baseline writeup):
//   - Phase 1: always enumerate Tsunami query results to compute exact deg[q]
//     and W=|J|. While W <= J* (cfg.run.j_star), also materialize join pairs
//     into an array AllPairs.
//   - If join finishes with W <= J*: sample directly from AllPairs.
//   - Otherwise: discard AllPairs and fall back to the two-pass rank sampling
//     protocol using the precomputed deg[] (no need to re-count).

#include "sjs/baselines/tsunami/sampling.h"

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
namespace tsunami {

template <int Dim, class T = Scalar>
class TsunamiAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 1, "TsunamiAdaptiveBaseline requires Dim >= 1");
  static_assert(2 * Dim <= kMaxSupportedDim,
                "TsunamiAdaptiveBaseline uses Point<2*Dim>; increase kMaxSupportedDim or reduce Dim");

  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;

  Method method() const noexcept override { return Method::Tsunami; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "tsunami_adaptive"; }

  void Reset() override {
    prep_.Reset();
    counted_ = false;
    W_ = 0;
    deg_order_.clear();
    mode_ = Mode::Unknown;
    all_pairs_.clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    Reset();
    if (!prep_.Build(ds, cfg, phases, err)) return false;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic
    if (!out) {
      if (err) *err = "TsunamiAdaptiveBaseline::Count: out is null";
      return false;
    }
    *out = MakeExactCount(0);

    if (!prep_.built()) {
      if (err) *err = "TsunamiAdaptiveBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count_and_maybe_enumerate")
                         : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 JSTAR = cfg.run.j_star;

    counted_ = false;
    W_ = 0;
    deg_order_.assign(prep_.QueryOrder().size(), 0ULL);
    mode_ = (JSTAR > 0) ? Mode::EnumerateAll : Mode::CountOnly;
    all_pairs_.clear();
    if (mode_ == Mode::EnumerateAll) {
      // Reserve conservatively to avoid huge upfront allocations.
      all_pairs_.reserve(static_cast<usize>(std::min<u64>(JSTAR, 1'000'000ULL)));
    }

    const auto& rel_q = prep_.QueryRel();
    const auto& order = prep_.QueryOrder();

    for (usize qi = 0; qi < order.size(); ++qi) {
      const u32 q_idx = order[qi];
      const BoxT& q = rel_q.boxes[static_cast<usize>(q_idx)];

      typename detail::TsunamiEmbedding<Dim, T>::RangeBox range;
      if (!prep_.MakeRangeForQuery(q, &range)) {
        deg_order_[qi] = 0;
        continue;
      }

      typename detail::TsunamiPreproc<Dim, T>::KD::RangeIter it(&prep_.Index(), range);

      u64 deg = 0;
      usize data_pt = 0;
      while (it.Next(&data_pt)) {
        ++deg;
        if (W_ == std::numeric_limits<u64>::max()) {
          if (err) *err = "TsunamiAdaptiveBaseline::Count: |J| overflowed u64";
          return false;
        }
        ++W_;

        if (mode_ == Mode::EnumerateAll) {
          if (all_pairs_.size() < static_cast<usize>(JSTAR)) {
            all_pairs_.push_back(prep_.MakePair(q_idx, static_cast<u32>(data_pt)));
          } else {
            // Switch to count-only and free memory.
            mode_ = Mode::CountOnly;
            std::vector<PairId>().swap(all_pairs_);
          }
        }
      }

      deg_order_[qi] = deg;
    }

    counted_ = true;
    *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!prep_.built()) {
      if (err) *err = "TsunamiAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "TsunamiAdaptiveBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "TsunamiAdaptiveBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    if (!counted_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    if (W_ == 0) {
      return true;
    }

    const u64 JSTAR = cfg.run.j_star;
    if (mode_ == Mode::EnumerateAll) {
      // Small join branch: sample from stored list.
      if (JSTAR == 0) {
        if (err) *err = "TsunamiAdaptiveBaseline::Sample: internal error (EnumerateAll but JSTAR==0)";
        return false;
      }
      if (all_pairs_.size() != static_cast<usize>(W_)) {
        if (err) *err = "TsunamiAdaptiveBaseline::Sample: internal error (AllPairs not complete)";
        return false;
      }

      auto scoped = phases ? phases->Scoped("sample_from_all_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");
      out->pairs.resize(static_cast<usize>(t));
      for (u64 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(W_);
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(idx)];
      }
      return true;
    }

    // Large join branch: reuse deg[] computed in Phase 1 and run Pass 2 only.
    auto scoped = phases ? phases->Scoped("fallback_two_pass_rank_sample") : PhaseRecorder::ScopedPhase(nullptr, "");
    std::vector<PairId> pairs;
    if (!detail::TwoPassRankSampleUsingDeg(prep_, deg_order_, W_, t, rng, &pairs, phases, err)) {
      return false;
    }
    out->pairs = std::move(pairs);
    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!prep_.built()) {
      if (err) *err = "TsunamiAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::TsunamiJoinEnumerator<Dim, T>>(&prep_);
  }

 private:
  enum class Mode : u8 {
    Unknown = 0,
    EnumerateAll = 1,
    CountOnly = 2,
  };

  detail::TsunamiPreproc<Dim, T> prep_;

  bool counted_ = false;
  u64 W_ = 0;
  std::vector<u64> deg_order_;      // aligned with prep_.QueryOrder()
  Mode mode_ = Mode::Unknown;
  std::vector<PairId> all_pairs_;   // only when mode_==EnumerateAll
};

}  // namespace tsunami
}  // namespace baselines
}  // namespace sjs
