#pragma once
// sjs/baselines/kd_tree/adaptive.h
//
// KD-tree baseline (Variant::Adaptive).
//
// Baseline v2.0 design: Adaptive+Sampling
//   Phase 1 (one sweep):
//     - Always compute per-START weights w_e and W = sum_e w_e = |J|.
//     - While W <= J* (cfg.run.j_star), also enumerate and store all join pairs.
//     - If W exceeds J*, discard stored pairs and continue in COUNT_ONLY mode.
//   Sampling:
//     - If |J| <= J*, sample directly from the stored pair array (uniform, i.i.d., with replacement).
//     - Otherwise, reuse w_e and W from Phase 1 to do event-level alias+slot allocation,
//       then perform a second sweep and call KD.Sample locally per START.
//
// Strict alignment with Baseline v2.0 "implementation checklist":
//   #4: filter out events with w_e==0 before building alias.
//   #7: when switching, discard AllPairs but keep w_e and W.
//

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

  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;
  static constexpr int K = 2 * (Dim - 1);

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

    // 1) Build sweep events (must satisfy END-before-START and stable START tie-break).
    {
      auto _ = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = join::BuildSweepEvents(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    // 2) Build a mapping: event index -> START id (sid), for START events only.
    {
      auto _ = phases ? phases->Scoped("build_start_index") : PhaseRecorder::ScopedPhase(nullptr, "");
      start_id_of_event_.assign(events_.size(), kInvalidStartId);
      u32 sid = 0;
      for (usize i = 0; i < events_.size(); ++i) {
        if (events_[i].kind == join::EventKind::Start) {
          start_id_of_event_[i] = sid++;
        }
      }
      w_total_.assign(static_cast<usize>(sid), 0ULL);
    }

    // 3) Global rank embedding + 2 KD-trees (static structure, dynamic active counts).
    std::vector<detail::RankPoint<K>> pts_r;
    std::vector<detail::RankPoint<K>> pts_s;
    {
      auto _ = phases ? phases->Scoped("build_embedding") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!embed_.Build(ds.R, ds.S, &pts_r, &pts_s, phases, err)) return false;
    }
    {
      auto _ = phases ? phases->Scoped("build_kdtrees") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!kd_r_.Build(pts_r, err)) return false;
      if (!kd_s_.Build(pts_s, err)) return false;
    }

    // Free temporary point arrays; KD nodes contain copied points.
    std::vector<detail::RankPoint<K>>().swap(pts_r);
    std::vector<detail::RankPoint<K>>().swap(pts_s);

    built_ = true;
    weights_valid_ = false;
    W_ = 0;
    mode_ = Mode::Unknown;
    all_pairs_.clear();
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic
    if (!built_ || !ds_) {
      if (err) *err = "KDTreeAdaptiveBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "KDTreeAdaptiveBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count_and_maybe_enumerate")
                         : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 JSTAR = cfg.run.j_star;

    // Reset caches.
    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    W_ = 0;
    weights_valid_ = false;

    mode_ = Mode::EnumerateAll;
    all_pairs_.clear();
    if (JSTAR > 0) {
      all_pairs_.reserve(static_cast<usize>(std::min<u64>(JSTAR, 1'000'000ULL)));
    }

    // Start with empty active sets.
    kd_r_.ResetToEmpty();
    kd_s_.ResetToEmpty();

    // Sweep.
    for (usize ev_i = 0; ev_i < events_.size(); ++ev_i) {
      const join::Event& e = events_[ev_i];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) kd_r_.Deactivate(e.index);
        else kd_s_.Deactivate(e.index);
        continue;
      }

      // START
      const u32 sid = start_id_of_event_[ev_i];
      SJS_DASSERT(sid != kInvalidStartId);

      const BoxT& q = (e.side == join::Side::R) ? ds_->R.boxes[static_cast<usize>(e.index)]
                                                : ds_->S.boxes[static_cast<usize>(e.index)];

      detail::RankBox<K> range;
      const bool has_range = embed_.MakeQueryRange(q, &range);

      u64 w = 0;
      if (has_range) {
        w = (e.side == join::Side::R) ? kd_s_.Count(range) : kd_r_.Count(range);
      }
      w_total_[static_cast<usize>(sid)] = w;

      // Accumulate W with overflow guard.
      if (w > 0) {
        if (W_ > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "KDTreeAdaptiveBaseline::Count: |J| overflowed u64";
          return false;
        }
        W_ += w;
      }

      // While in EnumerateAll mode, decide whether to enumerate this START's local join block.
      if (mode_ == Mode::EnumerateAll) {
        if (W_ <= JSTAR) {
          // Enumerate this START's block if non-empty.
          if (w > 0 && has_range) {
            std::vector<u32> hits;
            if (e.side == join::Side::R) {
              kd_s_.Report(range, &hits);
              for (u32 h : hits) {
                all_pairs_.push_back(
                    PairId{ds_->R.GetId(static_cast<usize>(e.index)), ds_->S.GetId(static_cast<usize>(h))});
              }
            } else {
              kd_r_.Report(range, &hits);
              for (u32 h : hits) {
                all_pairs_.push_back(
                    PairId{ds_->R.GetId(static_cast<usize>(h)), ds_->S.GetId(static_cast<usize>(e.index))});
              }
            }
          }
        } else {
          // Switch to COUNT_ONLY and discard stored pairs (Baseline v2.0 checklist #7).
          mode_ = Mode::CountOnly;
          std::vector<PairId>().swap(all_pairs_);
        }
      }

      // Activate current box in its own KD-tree.
      if (e.side == join::Side::R) kd_r_.Activate(e.index);
      else kd_s_.Activate(e.index);
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

    const u64 t64 = cfg.run.t;
    if (t64 == 0) return true;
    if (t64 > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "KDTreeAdaptiveBaseline::Sample: t too large for u32 slots";
      return false;
    }
    const u32 t = static_cast<u32>(t64);

    // Ensure Phase 1 results exist (w_e and W).
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    if (W_ == 0) {
      // Empty join.
      return true;
    }

    // -----------------
    // Branch A: no switch happened -> all pairs stored (|J| <= JSTAR).
    // -----------------
    if (mode_ == Mode::EnumerateAll) {
      const u64 N = static_cast<u64>(all_pairs_.size());
      if (N != W_) {
        if (err) *err = "KDTreeAdaptiveBaseline::Sample: inconsistent all_pairs_ size vs W_";
        return false;
      }

      out->pairs.resize(static_cast<usize>(t));
      for (u32 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(N);
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(idx)];
      }
      return true;
    }

    // -----------------
    // Branch B: switched -> fallback to Sampling protocol (Phase2+Phase3).
    // -----------------
    auto scoped = phases ? phases->Scoped("adaptive_fallback_sampling")
                         : PhaseRecorder::ScopedPhase(nullptr, "");

    // Phase 2: build alias table over START events with w_e>0, then assign slots.
    std::vector<u32> nz_sids;
    std::vector<u64> nz_w;
    nz_sids.reserve(w_total_.size());
    nz_w.reserve(w_total_.size());
    for (u32 sid = 0; sid < static_cast<u32>(w_total_.size()); ++sid) {
      const u64 w = w_total_[static_cast<usize>(sid)];
      if (w == 0) continue;
      nz_sids.push_back(sid);
      nz_w.push_back(w);
    }
    if (nz_sids.empty()) {
      if (err) *err = "KDTreeAdaptiveBaseline::Sample: internal error (W_>0 but no positive w_e)";
      return false;
    }

    sampling::AliasTable alias;
    {
      auto _ = phases ? phases->Scoped("phase2_alias") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!alias.BuildFromU64(Span<const u64>(nz_w.data(), nz_w.size()), err)) {
        return false;
      }
    }

    struct SlotAssign {
      u32 sid;   // START id (original sid)
      u32 slot;  // output position
    };

    std::vector<SlotAssign> slots;
    slots.reserve(static_cast<usize>(t));
    {
      auto _ = phases ? phases->Scoped("phase2_assign") : PhaseRecorder::ScopedPhase(nullptr, "");
      for (u32 j = 0; j < t; ++j) {
        const u32 p = alias.Sample(rng);  // index into nz_sids/nz_w
        const u32 sid = nz_sids[static_cast<usize>(p)];
        slots.push_back(SlotAssign{sid, j});
      }
      std::sort(slots.begin(), slots.end(), [](const SlotAssign& a, const SlotAssign& b) {
        if (a.sid < b.sid) return true;
        if (b.sid < a.sid) return false;
        return a.slot < b.slot;
      });
    }

    // Phase 3: second sweep + local KD.Sample at STARTs with slots.
    out->pairs.resize(static_cast<usize>(t));
    {
      auto _ = phases ? phases->Scoped("phase3_sweep") : PhaseRecorder::ScopedPhase(nullptr, "");

      kd_r_.ResetToEmpty();
      kd_s_.ResetToEmpty();

      usize ptr = 0;
      for (usize ev_i = 0; ev_i < events_.size(); ++ev_i) {
        const join::Event& e = events_[ev_i];

        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) kd_r_.Deactivate(e.index);
          else kd_s_.Deactivate(e.index);
          continue;
        }

        const u32 sid = start_id_of_event_[ev_i];
        SJS_DASSERT(sid != kInvalidStartId);

        // Slots for this START.
        usize end = ptr;
        while (end < slots.size() && slots[end].sid == sid) ++end;
        const u32 need = static_cast<u32>(end - ptr);

        if (need > 0) {
          const BoxT& q = (e.side == join::Side::R) ? ds_->R.boxes[static_cast<usize>(e.index)]
                                                    : ds_->S.boxes[static_cast<usize>(e.index)];

          detail::RankBox<K> range;
          if (!embed_.MakeQueryRange(q, &range)) {
            if (err) *err = "KDTreeAdaptiveBaseline::Sample: empty query range for START that has slots";
            return false;
          }

          std::vector<u32> picked;
          picked.reserve(static_cast<usize>(need));

          bool ok = false;
          if (e.side == join::Side::R) ok = kd_s_.Sample(range, need, rng, &picked, err);
          else ok = kd_r_.Sample(range, need, rng, &picked, err);

          if (!ok) {
            if (err && err->empty()) {
              *err = "KDTreeAdaptiveBaseline::Sample: KD sampling failed";
            }
            return false;
          }

          SJS_DASSERT(picked.size() == static_cast<usize>(need));
          for (u32 j = 0; j < need; ++j) {
            const u32 slot = slots[ptr + j].slot;
            const u32 other = picked[static_cast<usize>(j)];
            if (e.side == join::Side::R) {
              out->pairs[static_cast<usize>(slot)] =
                  PairId{ds_->R.GetId(static_cast<usize>(e.index)), ds_->S.GetId(static_cast<usize>(other))};
            } else {
              out->pairs[static_cast<usize>(slot)] =
                  PairId{ds_->R.GetId(static_cast<usize>(other)), ds_->S.GetId(static_cast<usize>(e.index))};
            }
          }

          ptr = end;
          if (ptr == slots.size()) {
            // All requested samples have been filled; we can break early.
            // (Active sets are local to Phase 3 and will be reset next call anyway.)
            // break;
          }
        }

        // Activate current.
        if (e.side == join::Side::R) kd_r_.Activate(e.index);
        else kd_s_.Activate(e.index);
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
  std::vector<u64> w_total_;  // per-START weights w_e

  Mode mode_ = Mode::Unknown;
  std::vector<PairId> all_pairs_;  // materialized join pairs when |J| <= JSTAR
};

}  // namespace kd_tree
}  // namespace baselines
}  // namespace sjs
