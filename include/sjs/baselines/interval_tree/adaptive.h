#pragma once
// sjs/baselines/interval_tree/adaptive.h
//
// IntervalTree Baseline (v2.0) — Variant::Adaptive
// ------------------------------------------------
// Matches the "Adaptive-IT" baseline in:
//   docs/Baseline/IntervalTree Baseline v2.0.md
//
// Idea:
//   - Phase1 performs the usual counting sweep, while *optionally* enumerating pairs
//     into AllPairs as long as W stays <= J_star.
//   - If W exceeds J_star at any point, switch to CountOnly mode (discard AllPairs).
//   - Sample():
//       * Enumerate mode  -> uniform index sampling from AllPairs.
//       * CountOnly mode  -> Sampling-IT's Phase2+Phase3 (slot plan + resweep).
//
// This header reuses the same data structures and slot-plan builder in sampling.h.

#include "sjs/baselines/interval_tree/sampling.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace interval_tree {

template <int Dim, class T = Scalar>
class IntervalTreeAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::IntervalTree; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "interval_tree_adaptive"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;

    events_.clear();
    start_id_of_event_.clear();

    y_endpoints_.clear();
    ylo_idx_r_.clear();
    yhi_idx_r_.clear();
    ylo_idx_s_.clear();
    yhi_idx_s_.clear();
    m_a_ = 0;

    active_r_.Clear();
    active_s_.Clear();

    w_total_.clear();
    w_a_.clear();
    w_b_.clear();
    W_ = 0;
    weights_valid_ = false;

    mode_ = Mode::CountOnly;
    all_pairs_.clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;

    if constexpr (Dim != 2) {
      if (err) *err = "IntervalTreeAdaptiveBaseline: only Dim==2 is supported";
      return false;
    }

    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    ds_ = &ds;

    if (ds.R.Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        ds.S.Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "IntervalTreeAdaptiveBaseline::Build: relation size exceeds u32 limit";
      return false;
    }

    // Events.
    {
      auto ph = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = detail::BuildSweepEvents2D<Dim, T>(ds.R, ds.S);
    }

    // START id mapping.
    {
      auto ph = phases ? phases->Scoped("build_start_ids") : PhaseRecorder::ScopedPhase(nullptr, "");
      start_id_of_event_.assign(events_.size(), -1);
      i32 sid = 0;
      for (usize i = 0; i < events_.size(); ++i) {
        if (events_[i].kind == join::EventKind::Start) {
          start_id_of_event_[i] = sid++;
        }
      }
      const usize E = static_cast<usize>(sid);
      w_total_.assign(E, 0ULL);
      w_a_.assign(E, 0ULL);
      w_b_.assign(E, 0ULL);
    }

    // y-endpoints {Ly,Ry}.
    {
      auto ph = phases ? phases->Scoped("build_y_endpoints") : PhaseRecorder::ScopedPhase(nullptr, "");
      y_endpoints_.clear();
      y_endpoints_.reserve(2 * (ds.R.Size() + ds.S.Size()));
      for (const auto& b : ds.R.boxes) {
        y_endpoints_.push_back(b.lo.v[1]);
        y_endpoints_.push_back(b.hi.v[1]);
      }
      for (const auto& b : ds.S.boxes) {
        y_endpoints_.push_back(b.lo.v[1]);
        y_endpoints_.push_back(b.hi.v[1]);
      }
      std::sort(y_endpoints_.begin(), y_endpoints_.end());
      y_endpoints_.erase(std::unique(y_endpoints_.begin(), y_endpoints_.end()), y_endpoints_.end());
      const u32 m_end = static_cast<u32>(y_endpoints_.size());
      m_a_ = (m_end >= 2) ? (m_end - 1) : 0;
    }

    // Precompute IT-A indices.
    {
      auto ph = phases ? phases->Scoped("build_y_indices") : PhaseRecorder::ScopedPhase(nullptr, "");
      auto lb_end = [&](T v) -> u32 {
        const auto it = std::lower_bound(y_endpoints_.begin(), y_endpoints_.end(), v);
        return static_cast<u32>(it - y_endpoints_.begin());
      };

      ylo_idx_r_.resize(ds.R.Size());
      yhi_idx_r_.resize(ds.R.Size());
      for (usize i = 0; i < ds.R.Size(); ++i) {
        const auto& b = ds.R.boxes[i];
        const u32 lo = lb_end(b.lo.v[1]);
        const u32 hi = lb_end(b.hi.v[1]);
        if (!(lo < hi)) {
          if (err) *err = "IntervalTreeAdaptiveBaseline::Build: invalid R y-interval (lo>=hi)";
          return false;
        }
        if (m_a_ > 0 && lo >= m_a_) {
          if (err) *err = "IntervalTreeAdaptiveBaseline::Build: R ylo maps to last endpoint (unexpected)";
          return false;
        }
        if (hi > m_a_) {
          if (err) *err = "IntervalTreeAdaptiveBaseline::Build: R yhi index out of range";
          return false;
        }
        ylo_idx_r_[i] = lo;
        yhi_idx_r_[i] = hi;
      }

      ylo_idx_s_.resize(ds.S.Size());
      yhi_idx_s_.resize(ds.S.Size());
      for (usize i = 0; i < ds.S.Size(); ++i) {
        const auto& b = ds.S.boxes[i];
        const u32 lo = lb_end(b.lo.v[1]);
        const u32 hi = lb_end(b.hi.v[1]);
        if (!(lo < hi)) {
          if (err) *err = "IntervalTreeAdaptiveBaseline::Build: invalid S y-interval (lo>=hi)";
          return false;
        }
        if (m_a_ > 0 && lo >= m_a_) {
          if (err) *err = "IntervalTreeAdaptiveBaseline::Build: S ylo maps to last endpoint (unexpected)";
          return false;
        }
        if (hi > m_a_) {
          if (err) *err = "IntervalTreeAdaptiveBaseline::Build: S yhi index out of range";
          return false;
        }
        ylo_idx_s_[i] = lo;
        yhi_idx_s_[i] = hi;
      }
    }

    // Init active samplers.
    {
      auto ph = phases ? phases->Scoped("build_active") : PhaseRecorder::ScopedPhase(nullptr, "");

      std::vector<detail::Key2<T>> keys_r;
      keys_r.reserve(ds.R.Size());
      for (usize i = 0; i < ds.R.Size(); ++i) keys_r.push_back(detail::Key2<T>{ds.R.boxes[i].lo.v[1], static_cast<u32>(i)});

      std::vector<detail::Key2<T>> keys_s;
      keys_s.reserve(ds.S.Size());
      for (usize i = 0; i < ds.S.Size(); ++i) keys_s.push_back(detail::Key2<T>{ds.S.boxes[i].lo.v[1], static_cast<u32>(i)});

      active_r_.Init(static_cast<u32>(ds.R.Size()), m_a_, std::move(keys_r));
      active_s_.Init(static_cast<u32>(ds.S.Size()), m_a_, std::move(keys_s));
    }

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic for counting

    if (!built_ || !ds_) {
      if (err) *err = "IntervalTreeAdaptiveBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "IntervalTreeAdaptiveBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count_adaptive") : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 J_star = cfg.run.j_star;

    // Start in Enumerate mode if we have a positive threshold; otherwise CountOnly.
    mode_ = (J_star > 0) ? Mode::Enumerate : Mode::CountOnly;
    all_pairs_.clear();

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    std::fill(w_a_.begin(), w_a_.end(), 0ULL);
    std::fill(w_b_.begin(), w_b_.end(), 0ULL);

    active_r_.ResetToEmpty();
    active_s_.ResetToEmpty();

    u64 W = 0;

    // Helper for appending all pairs for a START event into all_pairs_ (only when mode==Enumerate).
    auto append_pairs_for_start = [&](join::Side side, u32 idx) {
      std::vector<u32> cand;
      cand.reserve(256);

      if (side == join::Side::R) {
        const u32 a_idx = ylo_idx_r_[static_cast<usize>(idx)];
        const T a = ds_->R.boxes[static_cast<usize>(idx)].lo.v[1];
        const T b = ds_->R.boxes[static_cast<usize>(idx)].hi.v[1];
        active_s_.ReportA(a_idx, &cand);
        active_s_.ReportB(a, b, &cand);

        const Id qid = ds_->R.GetId(static_cast<usize>(idx));
        for (u32 other : cand) {
          all_pairs_.push_back(PairId{qid, ds_->S.GetId(static_cast<usize>(other))});
        }
      } else {
        const u32 a_idx = ylo_idx_s_[static_cast<usize>(idx)];
        const T a = ds_->S.boxes[static_cast<usize>(idx)].lo.v[1];
        const T b = ds_->S.boxes[static_cast<usize>(idx)].hi.v[1];
        active_r_.ReportA(a_idx, &cand);
        active_r_.ReportB(a, b, &cand);

        const Id qid = ds_->S.GetId(static_cast<usize>(idx));
        for (u32 other : cand) {
          all_pairs_.push_back(PairId{ds_->R.GetId(static_cast<usize>(other)), qid});
        }
      }
    };

    for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
      const auto& e = events_[ev_pos];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) active_r_.Erase(e.index);
        else active_s_.Erase(e.index);
        continue;
      }

      // START
      const i32 sid_i32 = start_id_of_event_[ev_pos];
      SJS_DASSERT(sid_i32 >= 0);
      const usize sid = static_cast<usize>(sid_i32);

      u64 wa = 0, wb = 0, w = 0;

      if (e.side == join::Side::R) {
        const u32 a_idx = ylo_idx_r_[static_cast<usize>(e.index)];
        const u32 hi_idx = yhi_idx_r_[static_cast<usize>(e.index)];
        const T a = ds_->R.boxes[static_cast<usize>(e.index)].lo.v[1];
        const T b = ds_->R.boxes[static_cast<usize>(e.index)].hi.v[1];

        wa = active_s_.CountA(a_idx);
        wb = active_s_.CountB(a, b);
        w = wa + wb;

        w_a_[sid] = wa;
        w_b_[sid] = wb;
        w_total_[sid] = w;

        if (w != 0 && W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "IntervalTreeAdaptiveBaseline::Count: |J| overflow (u64)";
          return false;
        }
        W += w;

        // If still in Enumerate mode, append pairs for this START as long as W <= J_star.
        if (mode_ == Mode::Enumerate) {
          if (W <= J_star) {
            append_pairs_for_start(join::Side::R, e.index);
          } else {
            mode_ = Mode::CountOnly;
            all_pairs_.clear();
          }
        }

        active_r_.Insert(e.index, a_idx, hi_idx);

      } else {  // S START
        const u32 a_idx = ylo_idx_s_[static_cast<usize>(e.index)];
        const u32 hi_idx = yhi_idx_s_[static_cast<usize>(e.index)];
        const T a = ds_->S.boxes[static_cast<usize>(e.index)].lo.v[1];
        const T b = ds_->S.boxes[static_cast<usize>(e.index)].hi.v[1];

        wa = active_r_.CountA(a_idx);
        wb = active_r_.CountB(a, b);
        w = wa + wb;

        w_a_[sid] = wa;
        w_b_[sid] = wb;
        w_total_[sid] = w;

        if (w != 0 && W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "IntervalTreeAdaptiveBaseline::Count: |J| overflow (u64)";
          return false;
        }
        W += w;

        if (mode_ == Mode::Enumerate) {
          if (W <= J_star) {
            append_pairs_for_start(join::Side::S, e.index);
          } else {
            mode_ = Mode::CountOnly;
            all_pairs_.clear();
          }
        }

        active_s_.Insert(e.index, a_idx, hi_idx);
      }
    }

    W_ = W;
    weights_valid_ = true;

    // Sanity: if we stayed in Enumerate mode, we must have materialized exactly W pairs.
    if (mode_ == Mode::Enumerate) {
      SJS_DASSERT(all_pairs_.size() == static_cast<usize>(W_));
    }

    *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "IntervalTreeAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "IntervalTreeAdaptiveBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "IntervalTreeAdaptiveBaseline::Sample: out is null";
      return false;
    }

    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    if (t == 0) return true;

    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }
    if (W_ == 0) return true;

    out->pairs.resize(static_cast<usize>(t));

    if (mode_ == Mode::Enumerate) {
      auto scoped = phases ? phases->Scoped("sample_small_join") : PhaseRecorder::ScopedPhase(nullptr, "");
      for (u32 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(W_);
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(idx)];
      }
      return true;
    }

    // CountOnly mode: reuse Sampling-IT Phase2+Phase3 (slot plan + resweep).
    detail::SlotPlan plan;
    {
      auto scoped2 = phases ? phases->Scoped("phase2_slots") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!detail::BuildSlotPlan(w_total_, w_a_, w_b_, t, rng, &plan, err)) return false;
    }

    {
      auto scoped3 = phases ? phases->Scoped("phase3_resweep") : PhaseRecorder::ScopedPhase(nullptr, "");
      active_r_.ResetToEmpty();
      active_s_.ResetToEmpty();

      std::vector<u32> tmp;
      tmp.reserve(1024);

      const usize E = w_total_.size();
      SJS_DASSERT(plan.a_off.size() == E + 1);
      SJS_DASSERT(plan.b_off.size() == E + 1);

      for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
        const auto& e = events_[ev_pos];

        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) active_r_.Erase(e.index);
          else active_s_.Erase(e.index);
          continue;
        }

        const i32 sid_i32 = start_id_of_event_[ev_pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 sid = static_cast<u32>(sid_i32);

        const u32 a_begin = plan.a_off[static_cast<usize>(sid)];
        const u32 a_end = plan.a_off[static_cast<usize>(sid + 1)];
        const u32 b_begin = plan.b_off[static_cast<usize>(sid)];
        const u32 b_end = plan.b_off[static_cast<usize>(sid + 1)];
        const u32 kA = a_end - a_begin;
        const u32 kB = b_end - b_begin;

        if (e.side == join::Side::R) {
          const u32 a_idx = ylo_idx_r_[static_cast<usize>(e.index)];
          const u32 hi_idx = yhi_idx_r_[static_cast<usize>(e.index)];
          const T a = ds_->R.boxes[static_cast<usize>(e.index)].lo.v[1];
          const T b = ds_->R.boxes[static_cast<usize>(e.index)].hi.v[1];

          if (kA > 0) {
            if (!active_s_.SampleA(a_idx, kA, rng, &tmp)) {
              if (err) *err = "IntervalTreeAdaptiveBaseline::Sample: empty A-set in phase3 (R-start)";
              return false;
            }
            for (u32 i = 0; i < kA; ++i) {
              const u32 slot = plan.a_slots[static_cast<usize>(a_begin + i)];
              const u32 other = tmp[static_cast<usize>(i)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(e.index)),
                                                           ds_->S.GetId(static_cast<usize>(other))};
            }
          }

          if (kB > 0) {
            if (!active_s_.SampleB(a, b, kB, rng, &tmp)) {
              if (err) *err = "IntervalTreeAdaptiveBaseline::Sample: empty B-set in phase3 (R-start)";
              return false;
            }
            for (u32 i = 0; i < kB; ++i) {
              const u32 slot = plan.b_slots[static_cast<usize>(b_begin + i)];
              const u32 other = tmp[static_cast<usize>(i)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(e.index)),
                                                           ds_->S.GetId(static_cast<usize>(other))};
            }
          }

          active_r_.Insert(e.index, a_idx, hi_idx);

        } else {
          const u32 a_idx = ylo_idx_s_[static_cast<usize>(e.index)];
          const u32 hi_idx = yhi_idx_s_[static_cast<usize>(e.index)];
          const T a = ds_->S.boxes[static_cast<usize>(e.index)].lo.v[1];
          const T b = ds_->S.boxes[static_cast<usize>(e.index)].hi.v[1];

          if (kA > 0) {
            if (!active_r_.SampleA(a_idx, kA, rng, &tmp)) {
              if (err) *err = "IntervalTreeAdaptiveBaseline::Sample: empty A-set in phase3 (S-start)";
              return false;
            }
            for (u32 i = 0; i < kA; ++i) {
              const u32 slot = plan.a_slots[static_cast<usize>(a_begin + i)];
              const u32 other = tmp[static_cast<usize>(i)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(other)),
                                                           ds_->S.GetId(static_cast<usize>(e.index))};
            }
          }

          if (kB > 0) {
            if (!active_r_.SampleB(a, b, kB, rng, &tmp)) {
              if (err) *err = "IntervalTreeAdaptiveBaseline::Sample: empty B-set in phase3 (S-start)";
              return false;
            }
            for (u32 i = 0; i < kB; ++i) {
              const u32 slot = plan.b_slots[static_cast<usize>(b_begin + i)];
              const u32 other = tmp[static_cast<usize>(i)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(other)),
                                                           ds_->S.GetId(static_cast<usize>(e.index))};
            }
          }

          active_s_.Insert(e.index, a_idx, hi_idx);
        }
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");

    if (!built_ || !ds_) {
      if (err) *err = "IntervalTreeAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    return std::make_unique<detail::IntervalTreeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S);
  }

 private:
  enum class Mode { Enumerate, CountOnly };

  const DatasetT* ds_{nullptr};
  bool built_{false};

  // Preprocessed events and START-id mapping.
  std::vector<detail::SweepEvent<T>> events_;
  std::vector<i32> start_id_of_event_;

  // y endpoints and IT-A indices per rectangle.
  std::vector<T> y_endpoints_;
  u32 m_a_{0};
  std::vector<u32> ylo_idx_r_;
  std::vector<u32> yhi_idx_r_;
  std::vector<u32> ylo_idx_s_;
  std::vector<u32> yhi_idx_s_;

  // Active samplers.
  detail::IntervalTreeSampler2D<T> active_r_;
  detail::IntervalTreeSampler2D<T> active_s_;

  // Phase1 weights.
  std::vector<u64> w_total_;
  std::vector<u64> w_a_;
  std::vector<u64> w_b_;
  u64 W_{0};
  bool weights_valid_{false};

  // Adaptive state.
  Mode mode_{Mode::CountOnly};
  std::vector<PairId> all_pairs_;
};

}  // namespace interval_tree
}  // namespace baselines
}  // namespace sjs