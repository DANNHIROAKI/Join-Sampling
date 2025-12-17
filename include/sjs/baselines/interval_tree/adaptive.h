#pragma once
// sjs/baselines/interval_tree/adaptive.h
//
// Plane Sweep + IntervalTreeSampler baseline (Variant::Adaptive).
//
// This follows the "Adaptive+Sampling" variant in "IntervalTree Baseline.md":
//   - Phase1: sweep and compute exact weights w_e = w_e^A + w_e^B for every START.
//             While the running total W <= J*, also enumerate join pairs into
//             AllPairs (exactly |J| pairs in the small-join regime).
//             Once W exceeds J*, discard AllPairs and switch to COUNT_ONLY.
//   - If W <= J*: sample directly from AllPairs (i.i.d. uniform, with replacement).
//   - Else: run Phase2+Phase3 identical to Variant::Sampling.
//
// NOTE: The project also provides a baseline-agnostic adaptive runner
// (sjs/baselines/runners/adaptive_runner.h) that does pilot enumeration and can
// fall back to Sampling. This class is still useful standalone and matches the
// per-baseline description.

#include "sjs/baselines/interval_tree/sampling.h"

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
    weights_valid_ = false;
    W_ = 0;
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

    all_pairs_.clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;
    mode_ = Mode::Unknown;
    all_pairs_.clear();

    if constexpr (Dim != 2) {
      if (err) *err = "IntervalTreeAdaptiveBaseline: currently only Dim=2 is supported";
      return false;
    }

    if (ds.R.Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        ds.S.Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "IntervalTreeAdaptiveBaseline::Build: relation size exceeds u32 limit";
      return false;
    }

    // Build and order sweep events.
    events_ = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);

    // Assign START ids.
    start_id_of_event_.assign(events_.size(), -1);
    usize start_cnt = 0;
    for (usize i = 0; i < events_.size(); ++i) {
      if (events_[i].kind == join::EventKind::Start) {
        start_id_of_event_[i] = static_cast<i32>(start_cnt);
        ++start_cnt;
      }
    }
    start_event_pos_.assign(start_cnt, 0);
    for (usize i = 0; i < events_.size(); ++i) {
      const i32 sid = start_id_of_event_[i];
      if (sid >= 0) start_event_pos_[static_cast<usize>(sid)] = i;
    }

    // Build y-domain from unique y-lower values.
    y_coords_.clear();
    y_coords_.reserve(ds.R.Size() + ds.S.Size());
    for (usize i = 0; i < ds.R.Size(); ++i) y_coords_.push_back(ds.R.boxes[i].lo.v[1]);
    for (usize i = 0; i < ds.S.Size(); ++i) y_coords_.push_back(ds.S.boxes[i].lo.v[1]);
    std::sort(y_coords_.begin(), y_coords_.end());
    y_coords_.erase(std::unique(y_coords_.begin(), y_coords_.end()), y_coords_.end());

    const u32 m = static_cast<u32>(y_coords_.size());

    // Precompute ranks.
    auto lb_rank = [&](T y) -> u32 {
      const auto it = std::lower_bound(y_coords_.begin(), y_coords_.end(), y);
      return static_cast<u32>(it - y_coords_.begin());
    };

    ylo_rank_r_.resize(ds.R.Size());
    yhi_lb_rank_r_.resize(ds.R.Size());
    for (usize i = 0; i < ds.R.Size(); ++i) {
      const auto& b = ds.R.boxes[i];
      ylo_rank_r_[i] = lb_rank(b.lo.v[1]);
      yhi_lb_rank_r_[i] = lb_rank(b.hi.v[1]);
    }
    ylo_rank_s_.resize(ds.S.Size());
    yhi_lb_rank_s_.resize(ds.S.Size());
    for (usize i = 0; i < ds.S.Size(); ++i) {
      const auto& b = ds.S.boxes[i];
      ylo_rank_s_[i] = lb_rank(b.lo.v[1]);
      yhi_lb_rank_s_[i] = lb_rank(b.hi.v[1]);
    }

    active_r_.Init(static_cast<u32>(ds.R.Size()), m);
    active_s_.Init(static_cast<u32>(ds.S.Size()), m);

    w_total_.assign(start_cnt, 0ULL);
    w_a_.assign(start_cnt, 0ULL);
    w_b_.assign(start_cnt, 0ULL);

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic
    if (!built_ || !ds_) {
      if (err) *err = "IntervalTreeAdaptiveBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "IntervalTreeAdaptiveBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    // Effective enumeration threshold. If enum_cap>0, never enumerate more than that.
    const u64 J_star_cfg = cfg.run.j_star;
    const u64 cap_cfg = cfg.run.enum_cap;
    const u64 J_star = (cap_cfg == 0) ? J_star_cfg : std::min(J_star_cfg, cap_cfg);

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    std::fill(w_a_.begin(), w_a_.end(), 0ULL);
    std::fill(w_b_.begin(), w_b_.end(), 0ULL);

    W_ = 0;
    weights_valid_ = false;
    all_pairs_.clear();

    if (J_star == 0) {
      mode_ = Mode::CountOnly;
    } else {
      mode_ = Mode::Enumerate;
      // A heuristic reserve to reduce reallocs; bounded.
      all_pairs_.reserve(static_cast<usize>(std::min<u64>(J_star, 1'000'000ULL)));
    }

    // Ensure active indices are empty.
    active_r_.ResetToEmpty();
    active_s_.ResetToEmpty();

    u64 W = 0;
    std::vector<u32> tmp;

    for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
      const auto& e = events_[ev_pos];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) {
          const u32 h = static_cast<u32>(e.index);
          active_r_.Erase(h);
        } else {
          const u32 h = static_cast<u32>(e.index);
          active_s_.Erase(h);
        }
        continue;
      }

      // START.
      const i32 sid_i32 = start_id_of_event_[ev_pos];
      SJS_DASSERT(sid_i32 >= 0);
      const usize sid = static_cast<usize>(sid_i32);

      if (e.side == join::Side::R) {
        const u32 q_ylo = ylo_rank_r_[e.index];
        const u32 q_yhi = yhi_lb_rank_r_[e.index];

        const u64 wa = active_s_.CountA(q_ylo);
        const u64 wb = active_s_.CountB(q_ylo, q_yhi);
        const u64 w = wa + wb;

        w_a_[sid] = wa;
        w_b_[sid] = wb;
        w_total_[sid] = w;

        if (w != 0 && W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "IntervalTreeAdaptiveBaseline::Count: join size overflow (u64)";
          return false;
        }
        W += w;

        // "Count first, then threshold" (as emphasized in the markdown).
        if (mode_ == Mode::Enumerate) {
          if (W <= J_star) {
            tmp.clear();
            active_s_.ReportA(q_ylo, &tmp);
            active_s_.ReportB(q_ylo, q_yhi, &tmp);
            const Id rid = ds_->R.GetId(e.index);
            for (u32 s_idx : tmp) {
              all_pairs_.push_back(PairId{rid, ds_->S.GetId(static_cast<usize>(s_idx))});
            }
          } else {
            mode_ = Mode::CountOnly;
            all_pairs_.clear();
            all_pairs_.shrink_to_fit();
          }
        }

        active_r_.Insert(static_cast<u32>(e.index), q_ylo, q_yhi);

      } else {
        const u32 q_ylo = ylo_rank_s_[e.index];
        const u32 q_yhi = yhi_lb_rank_s_[e.index];

        const u64 wa = active_r_.CountA(q_ylo);
        const u64 wb = active_r_.CountB(q_ylo, q_yhi);
        const u64 w = wa + wb;

        w_a_[sid] = wa;
        w_b_[sid] = wb;
        w_total_[sid] = w;

        if (w != 0 && W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "IntervalTreeAdaptiveBaseline::Count: join size overflow (u64)";
          return false;
        }
        W += w;

        if (mode_ == Mode::Enumerate) {
          if (W <= J_star) {
            tmp.clear();
            active_r_.ReportA(q_ylo, &tmp);
            active_r_.ReportB(q_ylo, q_yhi, &tmp);
            const Id sid_id = ds_->S.GetId(e.index);
            for (u32 r_idx : tmp) {
              all_pairs_.push_back(PairId{ds_->R.GetId(static_cast<usize>(r_idx)), sid_id});
            }
          } else {
            mode_ = Mode::CountOnly;
            all_pairs_.clear();
            all_pairs_.shrink_to_fit();
          }
        }

        active_s_.Insert(static_cast<u32>(e.index), q_ylo, q_yhi);
      }
    }

    // If we stayed in Enumerate mode, we must have AllPairs == W (unless J_star==0).
    if (mode_ == Mode::Enumerate && J_star > 0) {
      SJS_DASSERT(all_pairs_.size() == static_cast<usize>(W));
    }

    W_ = W;
    weights_valid_ = true;
    *out = MakeExactCount(W);
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

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u32 t = cfg.run.t;
    if (t == 0) return true;

    // Ensure Phase1 has been done.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    if (W_ == 0) return true;  // empty join

    // Branch A: small join, we have materialized all pairs.
    if (mode_ == Mode::Enumerate) {
      auto scoped = phases ? phases->Scoped("sample_from_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");

      const u64 W = static_cast<u64>(all_pairs_.size());
      if (W == 0) return true;

      out->pairs.resize(t);
      for (u32 j = 0; j < t; ++j) {
        const u64 idx = rng->UniformU64(W);
        out->pairs[j] = all_pairs_[static_cast<usize>(idx)];
      }
      return true;
    }

    // Branch B: large join, run Phase2+Phase3 identical to the Sampling variant.

    struct Assignment {
      u32 sid;
      u8 pat;   // 0=A, 1=B
      u32 slot; // output slot
    };

    std::vector<Assignment> assign;
    {
      auto scoped = phases ? phases->Scoped("phase2_assign") : PhaseRecorder::ScopedPhase(nullptr, "");

      sampling::AliasTable alias;
      if (!alias.BuildFromU64(Span<const u64>(w_total_), err)) {
        if (err && err->empty()) *err = "IntervalTreeAdaptiveBaseline::Sample: failed to build alias table";
        return false;
      }

      assign.reserve(t);
      for (u32 j = 0; j < t; ++j) {
        u32 sid = 0;
        u64 w = 0;
        for (int tries = 0; tries < 16; ++tries) {
          sid = static_cast<u32>(alias.Sample(rng));
          w = w_total_[sid];
          if (w > 0) break;
        }
        if (w == 0) {
          if (err) *err = "IntervalTreeAdaptiveBaseline::Sample: alias produced only zero-weight events (unexpected)";
          return false;
        }

        const u64 wa = w_a_[sid];
        const u64 wb = w_b_[sid];
        SJS_DASSERT(wa + wb == w);

        u8 pat = 0;
        if (wa == 0) pat = 1;
        else if (wb == 0) pat = 0;
        else {
          const u64 r = rng->UniformU64(w);
          pat = (r < wa) ? 0 : 1;
        }

        assign.push_back(Assignment{sid, pat, j});
      }

      std::sort(assign.begin(), assign.end(), [](const Assignment& a, const Assignment& b) {
        if (a.sid != b.sid) return a.sid < b.sid;
        if (a.pat != b.pat) return a.pat < b.pat;
        return a.slot < b.slot;
      });
    }

    out->pairs.resize(t);

    {
      auto scoped = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

      // Ensure active indices are empty.
      active_r_.ResetToEmpty();
      active_s_.ResetToEmpty();

      usize ptr = 0;
      std::vector<u32> tmp_handles;

      for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
        const auto& e = events_[ev_pos];

        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) {
            active_r_.Erase(static_cast<u32>(e.index));
          } else {
            active_s_.Erase(static_cast<u32>(e.index));
          }
          continue;
        }

        const i32 sid_i32 = start_id_of_event_[ev_pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 sid = static_cast<u32>(sid_i32);

        // Consume all assignments for this START.
        while (ptr < assign.size() && assign[ptr].sid == sid) {
          const u8 pat = assign[ptr].pat;
          const usize begin = ptr;
          while (ptr < assign.size() && assign[ptr].sid == sid && assign[ptr].pat == pat) ++ptr;
          const u32 k = static_cast<u32>(ptr - begin);
          if (k == 0) continue;

          bool ok = false;

          if (e.side == join::Side::R) {
            const u32 q_ylo = ylo_rank_r_[e.index];
            const u32 q_yhi = yhi_lb_rank_r_[e.index];

            if (pat == 0) ok = active_s_.SampleA(q_ylo, k, rng, &tmp_handles);
            else ok = active_s_.SampleB(q_ylo, q_yhi, k, rng, &tmp_handles);

            if (!ok) {
              if (err) *err = "IntervalTreeAdaptiveBaseline::Sample: empty candidate set in phase3 (R-start)";
              return false;
            }

            for (u32 i = 0; i < k; ++i) {
              const u32 slot = assign[begin + i].slot;
              const u32 other_h = tmp_handles[i];
              out->pairs[slot] = PairId{ds_->R.GetId(e.index), ds_->S.GetId(static_cast<usize>(other_h))};
            }

          } else {
            const u32 q_ylo = ylo_rank_s_[e.index];
            const u32 q_yhi = yhi_lb_rank_s_[e.index];

            if (pat == 0) ok = active_r_.SampleA(q_ylo, k, rng, &tmp_handles);
            else ok = active_r_.SampleB(q_ylo, q_yhi, k, rng, &tmp_handles);

            if (!ok) {
              if (err) *err = "IntervalTreeAdaptiveBaseline::Sample: empty candidate set in phase3 (S-start)";
              return false;
            }

            for (u32 i = 0; i < k; ++i) {
              const u32 slot = assign[begin + i].slot;
              const u32 other_h = tmp_handles[i];
              out->pairs[slot] = PairId{ds_->R.GetId(static_cast<usize>(other_h)), ds_->S.GetId(e.index)};
            }
          }
        }

        // Insert current rectangle after fulfilling its slots.
        if (e.side == join::Side::R) {
          active_r_.Insert(static_cast<u32>(e.index), ylo_rank_r_[e.index], yhi_lb_rank_r_[e.index]);
        } else {
          active_s_.Insert(static_cast<u32>(e.index), ylo_rank_s_[e.index], yhi_lb_rank_s_[e.index]);
        }
      }

      if (ptr != assign.size()) {
        if (err) *err = "IntervalTreeAdaptiveBaseline::Sample: internal error (not all assignments consumed)";
        return false;
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || !ds_) {
      if (err) *err = "IntervalTreeAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    return std::make_unique<detail::IntervalTreeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                        join::SideTieBreak::RBeforeS);
  }

 private:
  enum class Mode : u8 {
    Unknown = 0,
    Enumerate = 1,
    CountOnly = 2,
  };

  const DatasetT* ds_{nullptr};
  bool built_{false};

  std::vector<join::Event> events_;
  std::vector<i32> start_id_of_event_;  // -1 for END
  std::vector<usize> start_event_pos_;

  std::vector<T> y_coords_;
  std::vector<u32> ylo_rank_r_;
  std::vector<u32> yhi_lb_rank_r_;
  std::vector<u32> ylo_rank_s_;
  std::vector<u32> yhi_lb_rank_s_;

  detail::IntervalTreeSampler2D active_r_;
  detail::IntervalTreeSampler2D active_s_;

  std::vector<u64> w_total_;
  std::vector<u64> w_a_;
  std::vector<u64> w_b_;

  u64 W_{0};
  bool weights_valid_{false};
  Mode mode_{Mode::Unknown};

  std::vector<PairId> all_pairs_;
};

}  // namespace interval_tree
}  // namespace baselines
}  // namespace sjs
