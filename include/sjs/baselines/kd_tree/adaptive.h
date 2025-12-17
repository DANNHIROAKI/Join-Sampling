#pragma once
// sjs/baselines/kd_tree/adaptive.h
//
// KD-tree baseline (Variant::Adaptive).
//
// This implements the adaptive strategy described in "KD-Tree Baseline.md":
//   - Phase 1 sweep always computes w_e and W = sum_e w_e.
//   - While W <= J* (cfg.run.j_star), also enumerates and stores all join pairs.
//   - If W exceeds J*, discard stored pairs and continue in count-only mode.
//   - Sampling:
//       * If |J| <= J*, sample directly from the stored join list.
//       * Else, fall back to the three-phase sampling protocol (alias + second sweep)
//         exactly like Variant::Sampling.
//
// Note: The project also provides an adaptive runner that performs a pilot enumeration.
// This baseline-level adaptive implementation remains useful when you want the baseline
// to be self-contained.

#include "sjs/baselines/kd_tree/sampling.h"

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
namespace kd_tree {

template <int Dim, class T = Scalar>
class KDTreeAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "KDTreeAdaptiveBaseline requires Dim >= 2");

  using BoxT = Box<Dim, T>;
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::KDTree; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "kd_tree_adaptive"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;

    events_.clear();
    start_id_of_event_.clear();

    embed_.Reset();
    kd_r_.Clear();
    kd_s_.Clear();

    weights_valid_ = false;
    W_ = 0;
    w_total_.clear();

    mode_ = Mode::Unknown;
    all_pairs_.clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    ds_ = &ds;

    // Events.
    {
      auto _ = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = join::BuildSweepEvents(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    // Start-id map.
    start_id_of_event_.resize(events_.size());
    u32 start_cnt = 0;
    for (usize i = 0; i < events_.size(); ++i) {
      if (events_[i].kind == join::EventKind::Start) {
        start_id_of_event_[i] = start_cnt++;
      } else {
        start_id_of_event_[i] = kInvalidStartId;
      }
    }
    w_total_.assign(static_cast<usize>(start_cnt), 0ULL);

    // Rank embedding + KD-trees.
    std::vector<typename detail::ActiveKDTree<K>::PointT> pts_r;
    std::vector<typename detail::ActiveKDTree<K>::PointT> pts_s;
    {
      auto _ = phases ? phases->Scoped("build_embedding") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!embed_.Build(ds.R, ds.S, &pts_r, &pts_s, phases, err)) return false;
    }
    {
      auto _ = phases ? phases->Scoped("build_kdtrees") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!kd_r_.Build(pts_r, err)) return false;
      if (!kd_s_.Build(pts_s, err)) return false;
    }

    // Free temporary point arrays; nodes contain copied points.
    std::vector<typename detail::ActiveKDTree<K>::PointT>().swap(pts_r);
    std::vector<typename detail::ActiveKDTree<K>::PointT>().swap(pts_s);

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic
    if (!out) {
      if (err) *err = "KDTreeAdaptiveBaseline::Count: out is null";
      return false;
    }
    *out = MakeExactCount(0);

    if (!built_ || !ds_) {
      if (err) *err = "KDTreeAdaptiveBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1") : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 JSTAR = cfg.run.j_star;

    // Reset phase-1 caches.
    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    W_ = 0;
    weights_valid_ = false;

    mode_ = Mode::EnumerateAll;
    all_pairs_.clear();
    if (JSTAR > 0) {
      all_pairs_.reserve(static_cast<usize>(std::min<u64>(JSTAR, 1'000'000ULL)));
    }

    kd_r_.ResetToEmpty();
    kd_s_.ResetToEmpty();

    // Sweep.
    for (usize pos = 0; pos < events_.size(); ++pos) {
      const join::Event& e = events_[pos];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) {
          kd_r_.Deactivate(static_cast<u32>(e.index));
        } else {
          kd_s_.Deactivate(static_cast<u32>(e.index));
        }
        continue;
      }

      // START
      const u32 sid = start_id_of_event_[pos];
      SJS_DASSERT(sid != kInvalidStartId);

      const BoxT& q = (e.side == join::Side::R) ? ds_->R.boxes[e.index] : ds_->S.boxes[e.index];

      u64 w = 0;
      detail::RankBox<K> range;
      if (embed_.MakeQueryRange(q, &range)) {
        if (e.side == join::Side::R) {
          w = kd_s_.Count(range);
        } else {
          w = kd_r_.Count(range);
        }
      }

      w_total_[static_cast<usize>(sid)] = w;
      W_ += w;

      // Optional enumeration while W <= JSTAR.
      if (mode_ == Mode::EnumerateAll && (JSTAR > 0) && (W_ <= JSTAR)) {
        if (w > 0 && embed_.MakeQueryRange(q, &range)) {
          std::vector<u32> hits;
          if (e.side == join::Side::R) {
            kd_s_.Report(range, &hits);
            for (u32 h : hits) {
              all_pairs_.push_back(PairId{ds_->R.GetId(static_cast<usize>(e.index)),
                                          ds_->S.GetId(static_cast<usize>(h))});
            }
          } else {
            kd_r_.Report(range, &hits);
            for (u32 h : hits) {
              all_pairs_.push_back(PairId{ds_->R.GetId(static_cast<usize>(h)),
                                          ds_->S.GetId(static_cast<usize>(e.index))});
            }
          }
        }
      }

      // If W crosses JSTAR, switch to count-only and free memory.
      if (mode_ == Mode::EnumerateAll && (JSTAR > 0) && (W_ > JSTAR)) {
        mode_ = Mode::CountOnly;
        std::vector<PairId>().swap(all_pairs_);
      }

      // Activate current.
      if (e.side == join::Side::R) {
        kd_r_.Activate(static_cast<u32>(e.index));
      } else {
        kd_s_.Activate(static_cast<u32>(e.index));
      }
    }

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
      if (err) *err = "KDTreeAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "KDTreeAdaptiveBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "KDTreeAdaptiveBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    // Ensure phase-1 is available.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, rng, &tmp, phases, err)) return false;
    }

    if (W_ == 0) {
      // Empty join.
      return true;
    }

    // Branch A: join small and fully stored.
    if (mode_ == Mode::EnumerateAll) {
      // Sanity: in this branch we should have all pairs stored.
      // If JSTAR==0, EnumerateAll is meaningless; but Count never stores.
      const u64 N = static_cast<u64>(all_pairs_.size());
      if (N != W_) {
        // Non-fatal, but indicates a logic error.
        if (err) *err = "KDTreeAdaptiveBaseline::Sample: inconsistent all_pairs_ size vs W_";
        return false;
      }

      out->pairs.resize(static_cast<usize>(t));
      for (u64 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(N);
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(idx)];
      }
      return true;
    }

    // Branch B: large join -> fall back to alias + second sweep (same as Sampling variant).
    auto scoped = phases ? phases->Scoped("adaptive_fallback_sampling") : PhaseRecorder::ScopedPhase(nullptr, "");

    sampling::AliasTable alias;
    {
      auto _ = phases ? phases->Scoped("phase2_alias") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!alias.BuildFromU64(Span<const u64>(w_total_.data(), w_total_.size()), err)) {
        return false;
      }
    }

    struct Slot {
      u32 sid;
      u32 slot;
    };

    std::vector<Slot> slots;
    slots.reserve(static_cast<usize>(t));
    {
      auto _ = phases ? phases->Scoped("phase2_assign") : PhaseRecorder::ScopedPhase(nullptr, "");
      for (u32 j = 0; j < static_cast<u32>(t); ++j) {
        const u32 sid = alias.Sample(rng);
        slots.push_back(Slot{sid, j});
      }
      std::sort(slots.begin(), slots.end(), [](const Slot& a, const Slot& b) {
        if (a.sid < b.sid) return true;
        if (b.sid < a.sid) return false;
        return a.slot < b.slot;
      });
    }

    out->pairs.resize(static_cast<usize>(t));

    // Re-sweep and sample.
    {
      auto _ = phases ? phases->Scoped("phase3_sweep") : PhaseRecorder::ScopedPhase(nullptr, "");

      kd_r_.ResetToEmpty();
      kd_s_.ResetToEmpty();

      usize ptr = 0;
      for (usize pos = 0; pos < events_.size(); ++pos) {
        const join::Event& e = events_[pos];

        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) {
            kd_r_.Deactivate(static_cast<u32>(e.index));
          } else {
            kd_s_.Deactivate(static_cast<u32>(e.index));
          }
          continue;
        }

        const u32 sid = start_id_of_event_[pos];
        SJS_DASSERT(sid != kInvalidStartId);

        // Group slots for this start.
        usize end = ptr;
        while (end < slots.size() && slots[end].sid == sid) ++end;
        const u32 need = static_cast<u32>(end - ptr);

        if (need > 0) {
          const BoxT& q = (e.side == join::Side::R) ? ds_->R.boxes[e.index] : ds_->S.boxes[e.index];

          detail::RankBox<K> range;
          if (!embed_.MakeQueryRange(q, &range)) {
            if (err) *err = "KDTreeAdaptiveBaseline::Sample: picked a START with empty query range";
            return false;
          }

          std::vector<u32> picked;
          picked.reserve(static_cast<usize>(need));

          bool ok = false;
          if (e.side == join::Side::R) {
            ok = kd_s_.Sample(range, need, rng, &picked, err);
          } else {
            ok = kd_r_.Sample(range, need, rng, &picked, err);
          }
          if (!ok) {
            if (err && err->empty()) {
              *err = "KDTreeAdaptiveBaseline::Sample: KD sampling failed (empty?)";
            }
            return false;
          }

          SJS_DASSERT(picked.size() == static_cast<usize>(need));
          for (u32 i = 0; i < need; ++i) {
            const u32 slot_id = slots[ptr + i].slot;
            const u32 other = picked[static_cast<usize>(i)];
            if (e.side == join::Side::R) {
              (*out).pairs[static_cast<usize>(slot_id)] =
                  PairId{ds_->R.GetId(static_cast<usize>(e.index)), ds_->S.GetId(static_cast<usize>(other))};
            } else {
              (*out).pairs[static_cast<usize>(slot_id)] =
                  PairId{ds_->R.GetId(static_cast<usize>(other)), ds_->S.GetId(static_cast<usize>(e.index))};
            }
          }

          ptr = end;
        }

        // Activate current.
        if (e.side == join::Side::R) {
          kd_r_.Activate(static_cast<u32>(e.index));
        } else {
          kd_s_.Activate(static_cast<u32>(e.index));
        }

        if (ptr == slots.size()) {
          // We can still keep sweeping to properly deactivate, but since we reset at the
          // start of Phase 3 we can safely break here (all remaining slots are filled).
          // Keeping the sweep improves internal consistency if someone inspects the trees.
          // We'll just continue for simplicity.
        }
      }

      if (ptr != slots.size()) {
        if (err) *err = "KDTreeAdaptiveBaseline::Sample: some slots were not filled";
        return false;
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !ds_) {
      if (err) *err = "KDTreeAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::KDJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                              join::SideTieBreak::RBeforeS);
  }

 private:
  static constexpr int K = 2 * (Dim - 1);
  static constexpr u32 kInvalidStartId = std::numeric_limits<u32>::max();

  enum class Mode : u8 {
    Unknown = 0,
    EnumerateAll = 1,
    CountOnly = 2,
  };

  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  std::vector<join::Event> events_;
  std::vector<u32> start_id_of_event_;

  detail::GlobalRankEmbedding<Dim, T> embed_;
  detail::ActiveKDTree<K> kd_r_;
  detail::ActiveKDTree<K> kd_s_;

  bool weights_valid_ = false;
  u64 W_ = 0;
  std::vector<u64> w_total_;

  Mode mode_ = Mode::Unknown;
  std::vector<PairId> all_pairs_;
};

}  // namespace kd_tree
}  // namespace baselines
}  // namespace sjs
