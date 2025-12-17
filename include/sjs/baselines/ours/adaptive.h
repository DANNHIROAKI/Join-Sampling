#pragma once
// sjs/baselines/ours/adaptive.h
//
// Our method (Adaptive variant).
//
// The idea follows "Our Method.md":
//   - If the join size |J| is small (<= j*), materialize all join pairs and
//     sample from the materialized list (Enum+Sampling style).
//   - Otherwise, fall back to the three-phase sampling scheme (Phase1 weights,
//     alias on START events, Phase3 second sweep sampling) without materializing
//     J.
//
// In this codebase the adaptive runner may also implement its own pilot scan;
// nevertheless, this baseline implements the full adaptive logic internally so
// you can run it either way.

#include "sjs/baselines/ours/sampling.h"  // reuse 2D indices + enumerator wrapper

#include "sjs/core/assert.h"
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
namespace ours {

// Adaptive baseline.
template <int Dim, class T = Scalar>
class OursAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "ours_adaptive"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;
    small_pairs_.clear();
    mode_ = Mode::Unknown;

    events_.clear();
    start_id_of_event_.clear();
    start_event_pos_.clear();

    y_coords_.clear();
    ylo_rank_r_.clear();
    yhi_lb_rank_r_.clear();
    ylo_rank_s_.clear();
    yhi_lb_rank_s_.clear();

    active_r_.Clear();
    active_s_.Clear();

    w_total_.clear();
    w_a_.clear();
    w_b_.clear();
  }

  bool Build(const DatasetT& ds,
             const Config& cfg,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;

    if constexpr (Dim != 2) {
      if (err) *err = "OursAdaptiveBaseline: currently only Dim=2 is supported";
      return false;
    }

    auto ph = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;

    // Safety: we store indices as u32.
    if (ds.R.Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        ds.S.Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "OursAdaptiveBaseline: dataset too large for u32 handles";
      return false;
    }

    // Build sweep events.
    {
      auto ph2 = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      // We build the exact same event ordering as the sampling variant.
      join::PlaneSweepOptions opt;
      opt.axis = 0;
      opt.side_order = join::SideTieBreak::RBeforeS;
      opt.skip_axis_check = true;
      events_ = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, opt.axis, opt.side_order);
    }

    // Map START events to compact ids.
    start_id_of_event_.assign(events_.size(), -1);
    start_event_pos_.clear();
    start_event_pos_.reserve(ds.R.Size() + ds.S.Size());
    for (usize i = 0; i < events_.size(); ++i) {
      if (events_[i].kind == join::EventKind::Start) {
        const i32 sid = static_cast<i32>(start_event_pos_.size());
        start_id_of_event_[i] = sid;
        start_event_pos_.push_back(i);
      }
    }

    // Build y-domain = unique y-lower values across both relations.
    {
      auto ph2 = phases ? phases->Scoped("build_y_domain") : PhaseRecorder::ScopedPhase(nullptr, "");
      y_coords_.clear();
      y_coords_.reserve(ds.R.Size() + ds.S.Size());
      for (const auto& b : ds.R.boxes) y_coords_.push_back(b.lo.v[1]);
      for (const auto& b : ds.S.boxes) y_coords_.push_back(b.lo.v[1]);
      std::sort(y_coords_.begin(), y_coords_.end());
      y_coords_.erase(std::unique(y_coords_.begin(), y_coords_.end()), y_coords_.end());
      if (y_coords_.empty()) {
        if (err) *err = "OursAdaptiveBaseline: empty y-domain";
        return false;
      }
    }

    // Precompute ranks.
    {
      auto ph2 = phases ? phases->Scoped("build_y_ranks") : PhaseRecorder::ScopedPhase(nullptr, "");
      auto lb = [&](T x) -> u32 {
        const auto it = std::lower_bound(y_coords_.begin(), y_coords_.end(), x);
        return static_cast<u32>(std::distance(y_coords_.begin(), it));
      };

      ylo_rank_r_.resize(ds.R.Size());
      yhi_lb_rank_r_.resize(ds.R.Size());
      for (usize i = 0; i < ds.R.Size(); ++i) {
        const auto& b = ds.R.boxes[i];
        const u32 lo = lb(b.lo.v[1]);
        // y-lo must exist in domain.
        SJS_DASSERT(lo < y_coords_.size() && y_coords_[lo] == b.lo.v[1]);
        ylo_rank_r_[i] = lo;
        yhi_lb_rank_r_[i] = lb(b.hi.v[1]);  // may equal |domain|
      }

      ylo_rank_s_.resize(ds.S.Size());
      yhi_lb_rank_s_.resize(ds.S.Size());
      for (usize i = 0; i < ds.S.Size(); ++i) {
        const auto& b = ds.S.boxes[i];
        const u32 lo = lb(b.lo.v[1]);
        SJS_DASSERT(lo < y_coords_.size() && y_coords_[lo] == b.lo.v[1]);
        ylo_rank_s_[i] = lo;
        yhi_lb_rank_s_[i] = lb(b.hi.v[1]);
      }
    }

    // Init active indices per side.
    {
      auto ph2 = phases ? phases->Scoped("build_active_indices") : PhaseRecorder::ScopedPhase(nullptr, "");
      active_r_.Init(static_cast<u32>(ds.R.Size()), static_cast<u32>(y_coords_.size()));
      active_s_.Init(static_cast<u32>(ds.S.Size()), static_cast<u32>(y_coords_.size()));
    }

    // Allocate weights.
    const usize E = start_event_pos_.size();
    w_total_.assign(E, 0ULL);
    w_a_.assign(E, 0ULL);
    w_b_.assign(E, 0ULL);

    weights_valid_ = false;
    W_ = 0;
    small_pairs_.clear();
    mode_ = Mode::Unknown;
    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;
    if (!built_ || ds_ == nullptr) {
      if (err) *err = "OursAdaptiveBaseline::Count: call Build() first";
      return false;
    }

    auto ph = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 j_star = cfg.run.j_star;

    // Reset weights.
    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    std::fill(w_a_.begin(), w_a_.end(), 0ULL);
    std::fill(w_b_.begin(), w_b_.end(), 0ULL);

    // Materialize join pairs only if we stay within the threshold.
    small_pairs_.clear();
    if (j_star > 0) {
      // Reserve modestly to avoid huge allocations; exact cap is j_star.
      small_pairs_.reserve(static_cast<usize>(std::min<u64>(j_star, 1'000'000ULL)));
    }

    bool keep_pairs = (j_star > 0);

    // Active sets start empty. (They should be empty again at the end of the scan
    // because each rectangle is inserted exactly once and erased exactly once.)

    u64 W = 0;

    // Scratch buffers for reporting (only used if keep_pairs).
    std::vector<u32> tmp;
    tmp.reserve(1024);

    for (usize i = 0; i < events_.size(); ++i) {
      const auto& ev = events_[i];
      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        if (ev.side == join::Side::R) {
          active_r_.Erase(handle, ylo_rank_r_[ev.index]);
        } else {
          active_s_.Erase(handle, ylo_rank_s_[ev.index]);
        }
        continue;
      }

      // START event.
      const i32 sid = start_id_of_event_[i];
      SJS_DASSERT(sid >= 0);

      const bool q_is_r = (ev.side == join::Side::R);
      const u32 q_ylo = q_is_r ? ylo_rank_r_[ev.index] : ylo_rank_s_[ev.index];
      const u32 q_yhi = q_is_r ? yhi_lb_rank_r_[ev.index] : yhi_lb_rank_s_[ev.index];

      const detail::ActiveIndex2D& other = q_is_r ? active_s_ : active_r_;

      const u64 wa = other.CountA(q_ylo);
      const u64 wb = other.CountB(q_ylo, q_yhi);
      const u64 w = wa + wb;

      w_a_[static_cast<usize>(sid)] = wa;
      w_b_[static_cast<usize>(sid)] = wb;
      w_total_[static_cast<usize>(sid)] = w;

      if (w > 0) {
        if (keep_pairs) {
          // If adding these pairs would exceed the threshold, switch to large mode.
          // At this point small_pairs_.size() equals W (pairs collected so far).
          if (static_cast<u64>(small_pairs_.size()) + w > j_star) {
            keep_pairs = false;
            small_pairs_.clear();
            small_pairs_.shrink_to_fit();  // release memory; join is large anyway.
          } else {
            // Enumerate pairs for this START event via ReportA/ReportB.
            u64 reported = 0;

            // Pattern A.
            tmp.clear();
            other.ReportA(q_ylo, &tmp);
            reported += static_cast<u64>(tmp.size());
            for (u32 oh : tmp) {
              if (q_is_r) {
                small_pairs_.push_back(PairId{ds_->R.GetId(ev.index), ds_->S.GetId(oh)});
              } else {
                small_pairs_.push_back(PairId{ds_->R.GetId(oh), ds_->S.GetId(ev.index)});
              }
            }

            // Pattern B.
            tmp.clear();
            other.ReportB(q_ylo, q_yhi, &tmp);
            reported += static_cast<u64>(tmp.size());
            for (u32 oh : tmp) {
              if (q_is_r) {
                small_pairs_.push_back(PairId{ds_->R.GetId(ev.index), ds_->S.GetId(oh)});
              } else {
                small_pairs_.push_back(PairId{ds_->R.GetId(oh), ds_->S.GetId(ev.index)});
              }
            }

            // Optional sanity: we should have reported exactly w pairs.
            SJS_DASSERT(reported == w);
          }
        }

        W += w;
      }

      // Insert q into its side's active set.
      if (q_is_r) {
        active_r_.Insert(handle, q_ylo, q_yhi);
      } else {
        active_s_.Insert(handle, q_ylo, q_yhi);
      }
    }

    // After a full sweep, active sets should be empty again.

    W_ = W;
    weights_valid_ = true;

    if (keep_pairs) {
      // If we kept pairs, we must have materialized all W pairs.
      SJS_DASSERT(static_cast<u64>(small_pairs_.size()) == W_);
      mode_ = Mode::SmallMaterialized;
    } else {
      mode_ = Mode::LargeSampling;
    }

    if (out) *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || ds_ == nullptr) {
      if (err) *err = "OursAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "OursAdaptiveBaseline::Sample: null rng/out";
      return false;
    }

    auto ph = phases ? phases->Scoped("sample") : PhaseRecorder::ScopedPhase(nullptr, "");

    // Ensure weights (and small materialization decision) are ready.
    if (!weights_valid_ || mode_ == Mode::Unknown) {
      CountResult dummy;
      if (!Count(cfg, nullptr, &dummy, phases, err)) return false;
    }

    const u64 t = cfg.run.t;
    out->Clear();
    out->with_replacement = true;

    if (t == 0 || W_ == 0) {
      return true;
    }

    if (mode_ == Mode::SmallMaterialized) {
      // Sample from the stored list.
      out->pairs.resize(static_cast<usize>(t));
      for (u64 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(W_);
        out->pairs[static_cast<usize>(i)] = small_pairs_[static_cast<usize>(idx)];
      }
      out->weighted = false;
      return true;
    }

    // Large mode: run the same Phase2+Phase3 sampling as the sampling baseline.
    // Phase 2: alias on START events.
    sampling::AliasTable alias;
    if (!alias.BuildFromU64(Span<const u64>(w_total_))) {
      if (err) *err = "OursAdaptiveBaseline::Sample: failed to build alias table";
      return false;
    }

    struct Assignment {
      u32 start_id;
      u8 pat;      // 0=A, 1=B
      u32 slot;    // output slot
    };

    std::vector<Assignment> assign;
    assign.reserve(static_cast<usize>(t));

    // Draw assignments.
    for (u32 slot = 0; slot < static_cast<u32>(t); ++slot) {
      // alias.Sample returns a start-id.
      const u32 sid = static_cast<u32>(alias.Sample(rng));
      const u64 wa = w_a_[sid];
      const u64 wb = w_b_[sid];
      const u64 w = wa + wb;
      if (w == 0) {
        // Should not happen when W_>0; re-draw defensively.
        u32 tries = 0;
        u32 sid2 = sid;
        u64 w2 = w;
        while (tries++ < 16 && w2 == 0) {
          sid2 = static_cast<u32>(alias.Sample(rng));
          w2 = w_a_[sid2] + w_b_[sid2];
        }
        if (w2 == 0) {
          if (err) *err = "OursAdaptiveBaseline::Sample: alias produced zero-weight event repeatedly";
          return false;
        }
        // Use sid2.
        const u64 wa2 = w_a_[sid2];
        const u64 wb2 = w_b_[sid2];
        const u64 tot2 = wa2 + wb2;
        u8 pat = 0;
        if (wa2 == 0) pat = 1;
        else if (wb2 == 0) pat = 0;
        else pat = (rng->UniformU64(tot2) < wa2) ? 0 : 1;
        assign.push_back(Assignment{sid2, pat, slot});
        continue;
      }

      u8 pat = 0;
      if (wa == 0) pat = 1;
      else if (wb == 0) pat = 0;
      else pat = (rng->UniformU64(w) < wa) ? 0 : 1;

      assign.push_back(Assignment{sid, pat, slot});
    }

    std::sort(assign.begin(), assign.end(), [](const Assignment& a, const Assignment& b) {
      if (a.start_id != b.start_id) return a.start_id < b.start_id;
      if (a.pat != b.pat) return a.pat < b.pat;
      return a.slot < b.slot;
    });

    out->pairs.resize(static_cast<usize>(t));
    out->weighted = false;

    // Phase 3: second sweep.

    usize ptr = 0;
    std::vector<u32> sampled;
    sampled.reserve(1024);

    for (usize i = 0; i < events_.size(); ++i) {
      const auto& ev = events_[i];
      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        if (ev.side == join::Side::R) {
          active_r_.Erase(handle, ylo_rank_r_[ev.index]);
        } else {
          active_s_.Erase(handle, ylo_rank_s_[ev.index]);
        }
        continue;
      }

      const i32 sid_i = start_id_of_event_[i];
      SJS_DASSERT(sid_i >= 0);
      const u32 sid = static_cast<u32>(sid_i);

      const bool q_is_r = (ev.side == join::Side::R);
      const u32 q_ylo = q_is_r ? ylo_rank_r_[ev.index] : ylo_rank_s_[ev.index];
      const u32 q_yhi = q_is_r ? yhi_lb_rank_r_[ev.index] : yhi_lb_rank_s_[ev.index];

      const detail::ActiveIndex2D& other = q_is_r ? active_s_ : active_r_;

      // Process assignments for this START event.
      while (ptr < assign.size() && assign[ptr].start_id == sid) {
        const u8 pat = assign[ptr].pat;
        const usize begin = ptr;
        while (ptr < assign.size() && assign[ptr].start_id == sid && assign[ptr].pat == pat) {
          ++ptr;
        }
        const u32 k = static_cast<u32>(ptr - begin);

        sampled.clear();
        bool ok = false;
        if (pat == 0) {
          ok = other.SampleA(q_ylo, k, rng, &sampled);
        } else {
          ok = other.SampleB(q_ylo, q_yhi, k, rng, &sampled);
        }
        if (!ok || sampled.size() != k) {
          if (err) *err = "OursAdaptiveBaseline::Sample: failed to sample from K_e";
          return false;
        }

        for (u32 j = 0; j < k; ++j) {
          const u32 slot = assign[begin + j].slot;
          const u32 oh = sampled[j];
          if (q_is_r) {
            out->pairs[slot] = PairId{ds_->R.GetId(ev.index), ds_->S.GetId(oh)};
          } else {
            out->pairs[slot] = PairId{ds_->R.GetId(oh), ds_->S.GetId(ev.index)};
          }
        }
      }

      // Insert q.
      if (q_is_r) {
        active_r_.Insert(handle, q_ylo, q_yhi);
      } else {
        active_s_.Insert(handle, q_ylo, q_yhi);
      }
    }

    if (ptr != assign.size()) {
      if (err) *err = "OursAdaptiveBaseline::Sample: internal error (not all assignments consumed)";
      return false;
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || ds_ == nullptr) {
      if (err) *err = "OursAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    join::PlaneSweepOptions opt;
    opt.axis = 0;
    opt.side_order = join::SideTieBreak::RBeforeS;
    opt.skip_axis_check = true;

    return std::make_unique<detail::PlaneSweepEnumeratorWrapper<Dim, T>>(ds_->R, ds_->S, opt);
  }

 private:
  enum class Mode {
    Unknown,
    SmallMaterialized,
    LargeSampling,
  };

  const DatasetT* ds_{nullptr};
  bool built_{false};

  std::vector<join::Event> events_;
  std::vector<i32> start_id_of_event_;
  std::vector<usize> start_event_pos_;

  std::vector<T> y_coords_;
  std::vector<u32> ylo_rank_r_;
  std::vector<u32> yhi_lb_rank_r_;
  std::vector<u32> ylo_rank_s_;
  std::vector<u32> yhi_lb_rank_s_;

  detail::ActiveIndex2D active_r_;
  detail::ActiveIndex2D active_s_;

  std::vector<u64> w_total_;
  std::vector<u64> w_a_;
  std::vector<u64> w_b_;

  u64 W_{0};
  bool weights_valid_{false};

  std::vector<PairId> small_pairs_;
  Mode mode_{Mode::Unknown};
};

}  // namespace ours
}  // namespace baselines
}  // namespace sjs
