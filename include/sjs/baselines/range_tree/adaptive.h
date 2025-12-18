#pragma once
// sjs/baselines/range_tree/adaptive.h
//
// Plane Sweep + Dynamic Range-Tree baseline (Variant::Adaptive).
//
// Design reference: docs/Baseline/RangeTree Baseline v2.0.md
//
// Adaptive strategy (mirrors the project's common adaptive pattern):
//   - During Count(), we run the Phase-1 sweep to compute exact |J| and w_e.
//   - If |J| stays <= j_star, we simultaneously enumerate and store all join
//     pairs in memory (so sampling later is O(t)).
//   - If |J| exceeds j_star, we switch to "count-only" mode and do not store pairs.
//   - Sample():
//       * If we stored all pairs (small join): sample with replacement from them.
//       * Otherwise: perform Phase-2 alias + Phase-3 conditional range sampling
//         (same as the Sampling variant), reusing the w_e computed in Count().
//
// This variant is useful to avoid doing two full sweeps when joins are tiny.

#include "sjs/baselines/range_tree/sampling.h"
#include "sjs/sampling/alias_table.h"

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
namespace range_tree {

template <int Dim, class T = Scalar>
class RangeTreeAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim == 2, "RangeTreeAdaptiveBaseline is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;
  static constexpr int K = 2 * (Dim - 1);

  Method method() const noexcept override { return Method::RangeTree; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "range_tree_adaptive"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;
    mode_small_ = false;
    std::vector<PairId>().swap(all_pairs_);

    events_.clear();
    start_id_of_event_.clear();
    w_total_.clear();

    embed_.Reset();
    rt_r_.Clear();
    rt_s_.Clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    ds_ = &ds;

    // Events.
    {
      auto _ = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    // START ids.
    {
      auto _ = phases ? phases->Scoped("build_start_index") : PhaseRecorder::ScopedPhase(nullptr, "");
      start_id_of_event_.assign(events_.size(), kInvalidStartId);
      u32 sid = 0;
      for (usize i = 0; i < events_.size(); ++i) {
        if (events_[i].kind == join::EventKind::Start) start_id_of_event_[i] = sid++;
      }
      w_total_.assign(static_cast<usize>(sid), 0ULL);
    }

    // Embedding + range trees.
    {
      std::vector<detail::RankPoint<K>> pts_r;
      std::vector<detail::RankPoint<K>> pts_s;
      if (!embed_.Build(ds.R, ds.S, &pts_r, &pts_s, phases, err)) return false;
      auto _ = phases ? phases->Scoped("build_range_trees") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!rt_r_.Build(pts_r, err)) return false;
      if (!rt_s_.Build(pts_s, err)) return false;
    }

    built_ = true;
    weights_valid_ = false;
    W_ = 0;
    mode_small_ = false;
    std::vector<PairId>().swap(all_pairs_);
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("adaptive_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 j_star = cfg.run.j_star;
    const bool try_materialize = (j_star > 0);

    if (ds_->R.Size() == 0 || ds_->S.Size() == 0) {
      W_ = 0;
      weights_valid_ = true;
      mode_small_ = true;
      all_pairs_.clear();
      std::fill(w_total_.begin(), w_total_.end(), 0ULL);
      *out = MakeExactCount(0);
      return true;
    }

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    rt_r_.ResetToEmpty();
    rt_s_.ResetToEmpty();
    all_pairs_.clear();
    mode_small_ = try_materialize;

    u64 W = 0;

    // Sweep.
    for (usize ev_i = 0; ev_i < events_.size(); ++ev_i) {
      const join::Event& e = events_[ev_i];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) rt_r_.Deactivate(static_cast<u32>(e.index));
        else rt_s_.Deactivate(static_cast<u32>(e.index));
        continue;
      }

      // START
      const u32 sid = start_id_of_event_[ev_i];
      SJS_DASSERT(sid != kInvalidStartId);

      const BoxT& q = (e.side == join::Side::R)
                          ? ds_->R.boxes[static_cast<usize>(e.index)]
                          : ds_->S.boxes[static_cast<usize>(e.index)];

      detail::RankBox<K> range;
      u64 w = 0;
      if (embed_.MakeQueryRange(q, &range)) {
        if (e.side == join::Side::R) {
          w = rt_s_.Count(range);
        } else {
          w = rt_r_.Count(range);
        }
      }
      w_total_[static_cast<usize>(sid)] = w;

      if (w > 0) {
        if (W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "RangeTreeAdaptiveBaseline::Count: |J| overflowed u64";
          return false;
        }
        W += w;
      }

      // Optional materialization while join remains small.
      if (mode_small_ && w > 0) {
        if (W > j_star) {
          mode_small_ = false;
          std::vector<PairId>().swap(all_pairs_);
        } else {
          // Enumerate actual pairs for this START.
          std::vector<u32> hits;
          if (embed_.MakeQueryRange(q, &range)) {
            if (e.side == join::Side::R) {
              rt_s_.Report(range, &hits, /*st=*/nullptr);
              // Emit (R start, S hit)
              const Id rid = ds_->R.GetId(static_cast<usize>(e.index));
              for (u32 s_idx : hits) {
                all_pairs_.push_back(PairId{rid, ds_->S.GetId(static_cast<usize>(s_idx))});
              }
            } else {
              rt_r_.Report(range, &hits, /*st=*/nullptr);
              const Id sid_id = ds_->S.GetId(static_cast<usize>(e.index));
              for (u32 r_idx : hits) {
                all_pairs_.push_back(PairId{ds_->R.GetId(static_cast<usize>(r_idx)), sid_id});
              }
            }
          }
        }
      }

      // Activate current.
      if (e.side == join::Side::R) rt_r_.Activate(static_cast<u32>(e.index));
      else rt_s_.Activate(static_cast<u32>(e.index));
    }

    W_ = W;
    weights_valid_ = true;
    *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t64 = cfg.run.t;
    if (t64 == 0) return true;
    if (t64 > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Sample: t too large for u32 slots";
      return false;
    }
    const u32 t = static_cast<u32>(t64);

    // Ensure Count() ran so we have w_e and potentially all_pairs_.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }
    if (W_ == 0) return true;

    // If join is small and we materialized it, sampling is trivial.
    if (mode_small_) {
      if (all_pairs_.empty()) {
        // This can happen when W_==0 (handled above) or due to defensive clears.
        return true;
      }
      auto scoped = phases ? phases->Scoped("sample_from_materialized") : PhaseRecorder::ScopedPhase(nullptr, "");
      out->pairs.resize(static_cast<usize>(t));
      const u64 n = static_cast<u64>(all_pairs_.size());
      for (u32 i = 0; i < t; ++i) {
        const u64 j = rng->UniformU64(n);
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(j)];
      }
      return true;
    }

    // Otherwise fall back to 3-phase sampling (reuse w_total_ computed in Count()).
    {
      auto scoped = phases ? phases->Scoped("phase2_alias") : PhaseRecorder::ScopedPhase(nullptr, "");

      detail::NonZeroStartAlias ev_alias;
      if (!ev_alias.Build(Span<const u64>(w_total_.data(), w_total_.size()), err)) return false;

      struct SlotAssign {
        u32 sid;
        u32 slot;
      };
      std::vector<SlotAssign> asg;
      asg.reserve(static_cast<usize>(t));
      for (u32 j = 0; j < t; ++j) {
        const u32 sid = ev_alias.SampleSid(rng);
        asg.push_back(SlotAssign{sid, j});
      }
      std::sort(asg.begin(), asg.end(), [](const SlotAssign& a, const SlotAssign& b) {
        if (a.sid < b.sid) return true;
        if (b.sid < a.sid) return false;
        return a.slot < b.slot;
      });

      auto __ = phases ? phases->Scoped("phase3_sweep") : PhaseRecorder::ScopedPhase(nullptr, "");
      out->pairs.resize(static_cast<usize>(t));

      rt_r_.ResetToEmpty();
      rt_s_.ResetToEmpty();

      usize ptr = 0;
      for (usize ev_i = 0; ev_i < events_.size(); ++ev_i) {
        const join::Event& e = events_[ev_i];

        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) rt_r_.Deactivate(static_cast<u32>(e.index));
          else rt_s_.Deactivate(static_cast<u32>(e.index));
          continue;
        }

        const u32 sid = start_id_of_event_[ev_i];
        SJS_DASSERT(sid != kInvalidStartId);

        usize end = ptr;
        while (end < asg.size() && asg[end].sid == sid) ++end;
        const u32 k = static_cast<u32>(end - ptr);

        if (k > 0) {
          if (w_total_[static_cast<usize>(sid)] == 0) {
            if (err) *err = "RangeTreeAdaptiveBaseline::Sample: sampled a START with w_e==0 (unexpected)";
            return false;
          }

          const BoxT& q = (e.side == join::Side::R)
                              ? ds_->R.boxes[static_cast<usize>(e.index)]
                              : ds_->S.boxes[static_cast<usize>(e.index)];
          detail::RankBox<K> range;
          if (!embed_.MakeQueryRange(q, &range)) {
            if (err) *err = "RangeTreeAdaptiveBaseline::Sample: empty query range for START that has slots";
            return false;
          }

          std::vector<u32> picked;
          picked.reserve(static_cast<usize>(k));

          if (e.side == join::Side::R) {
            if (!rt_s_.Sample(range, k, rng, &picked, err)) return false;
            for (u32 j = 0; j < k; ++j) {
              const u32 slot = asg[ptr + j].slot;
              const u32 s_idx = picked[static_cast<usize>(j)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(e.index)),
                                                          ds_->S.GetId(static_cast<usize>(s_idx))};
            }
          } else {
            if (!rt_r_.Sample(range, k, rng, &picked, err)) return false;
            for (u32 j = 0; j < k; ++j) {
              const u32 slot = asg[ptr + j].slot;
              const u32 r_idx = picked[static_cast<usize>(j)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(r_idx)),
                                                          ds_->S.GetId(static_cast<usize>(e.index))};
            }
          }
          ptr = end;
        }

        if (e.side == join::Side::R) rt_r_.Activate(static_cast<u32>(e.index));
        else rt_s_.Activate(static_cast<u32>(e.index));
      }

      if (ptr != asg.size()) {
        if (err) *err = "RangeTreeAdaptiveBaseline::Sample: internal error (ptr != asg.size())";
        return false;
      }

      return true;
    }
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RangeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                 join::SideTieBreak::RBeforeS);
  }

 private:
  static constexpr u32 kInvalidStartId = std::numeric_limits<u32>::max();

  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  // Phase-1 derived.
  bool weights_valid_ = false;
  u64 W_ = 0;

  // Small-join mode.
  bool mode_small_ = false;
  std::vector<PairId> all_pairs_;

  // Preprocessing.
  std::vector<join::Event> events_;
  std::vector<u32> start_id_of_event_;
  std::vector<u64> w_total_;

  detail::GlobalRankEmbedding<Dim, T> embed_;
  detail::ActiveRangeTree<K> rt_r_;
  detail::ActiveRangeTree<K> rt_s_;
};

}  // namespace range_tree
}  // namespace baselines
}  // namespace sjs
