#pragma once
// sjs/baselines/ours/adaptive.h
//
// "Our Method" — Adaptive variant.
//
// Matches "Our Method.md":
//   Phase1: single sweep counts exact |J| and per-event weights w_e^A/w_e^B.
//           If |J| remains <= j*, also materialize all pairs into AllPairs.
//           Once |J| exceeds j*, discard AllPairs and keep counting only.
//   If |J| <= j*: sample by uniform indexing into AllPairs.
//   Else: run Phase2+Phase3 of the sampling variant (alias on w_e then second sweep).

#include "sjs/baselines/ours/sampling.h"  // reuse Ours2DContext + SlotPlan2D + ActiveIndex2D

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
namespace ours {

template <int Dim, class T = Scalar>
class OursAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "ours_adaptive"; }

  void Reset() override {
    built_ = false;
    ctx_.Reset();

    w_total_.clear();
    w_a_.clear();
    w_b_.clear();

    W_ = 0;
    weights_valid_ = false;

    small_pairs_.clear();
    mode_ = Mode::Unknown;
  }

  bool Build(const DatasetT& ds,
             const Config& cfg,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    Reset();

    if (!ctx_.Build(ds, phases, err)) return false;

    const usize E = ctx_.num_start_events();
    w_total_.assign(E, 0ULL);
    w_a_.assign(E, 0ULL);
    w_b_.assign(E, 0ULL);

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;

    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursAdaptiveBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count_adaptive") : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 j_star = cfg.run.j_star;

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    std::fill(w_a_.begin(), w_a_.end(), 0ULL);
    std::fill(w_b_.begin(), w_b_.end(), 0ULL);

    small_pairs_.clear();
    if (j_star > 0) {
      // reserve a modest amount; cap is j_star but that can be huge.
      small_pairs_.reserve(static_cast<usize>(std::min<u64>(j_star, 1'000'000ULL)));
    }

    bool keep_pairs = (j_star > 0);
    mode_ = Mode::Unknown;

    ctx_.ResetActive();

    auto& ar = ctx_.active_r();
    auto& as = ctx_.active_s();

    const auto& events = ctx_.events();
    const auto& sid_of_pos = ctx_.start_id_of_event();

    const auto& ylo_r = ctx_.ylo_rank_r();
    const auto& yhi_r = ctx_.yhi_lb_rank_r();
    const auto& ylo_s = ctx_.ylo_rank_s();
    const auto& yhi_s = ctx_.yhi_lb_rank_s();

    const auto* ds = ctx_.dataset();
    SJS_DASSERT(ds != nullptr);

    u64 W = 0;
    std::vector<u32> tmp;
    tmp.reserve(1024);

    for (usize pos = 0; pos < events.size(); ++pos) {
      const auto& ev = events[pos];
      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        if (ev.side == join::Side::R) ar.Erase(handle);
        else as.Erase(handle);
        continue;
      }

      // START
      const i32 sid_i32 = sid_of_pos[pos];
      SJS_DASSERT(sid_i32 >= 0);
      const u32 sid = static_cast<u32>(sid_i32);

      const bool q_is_r = (ev.side == join::Side::R);
      const u32 q_ylo = q_is_r ? ylo_r[ev.index] : ylo_s[ev.index];
      const u32 q_yhi = q_is_r ? yhi_r[ev.index] : yhi_s[ev.index];

      const detail::ActiveIndex2D& other = q_is_r ? as : ar;

      const u64 wa = other.CountA(q_ylo);
      const u64 wb = other.CountB(q_ylo, q_yhi);
      const u64 w = wa + wb;

      w_a_[sid] = wa;
      w_b_[sid] = wb;
      w_total_[sid] = w;

      const u64 W_new = W + w;

      // Adaptive decision (exactly as in the doc): decide using the updated W.
      if (keep_pairs && w > 0) {
        if (W_new > j_star) {
          // Switch to COUNT_ONLY: discard all previously materialized pairs.
          keep_pairs = false;
          small_pairs_.clear();
          small_pairs_.shrink_to_fit();
        } else {
          // Still within threshold: enumerate this event's pairs and append.
          tmp.clear();
          other.ReportA(q_ylo, &tmp);
          for (u32 oh : tmp) {
            if (q_is_r) small_pairs_.push_back(PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)});
            else small_pairs_.push_back(PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)});
          }

          tmp.clear();
          other.ReportB(q_ylo, q_yhi, &tmp);
          for (u32 oh : tmp) {
            if (q_is_r) small_pairs_.push_back(PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)});
            else small_pairs_.push_back(PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)});
          }
        }
      }

      W = W_new;

      // Insert q
      if (q_is_r) ar.Insert(handle, q_ylo, q_yhi);
      else as.Insert(handle, q_ylo, q_yhi);
    }

    ctx_.ResetActive();

    W_ = W;
    weights_valid_ = true;

    if (keep_pairs) {
      mode_ = Mode::SmallMaterialized;
      SJS_DASSERT(static_cast<u64>(small_pairs_.size()) == W_);
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
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "OursAdaptiveBaseline::Sample: null rng/out";
      return false;
    }

    auto scoped = phases ? phases->Scoped("sample") : PhaseRecorder::ScopedPhase(nullptr, "");

    // Ensure Phase1 is done (counts + mode).
    if (!weights_valid_ || mode_ == Mode::Unknown) {
      CountResult dummy;
      if (!Count(cfg, nullptr, &dummy, phases, err)) return false;
    }

    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    if (t == 0 || W_ == 0) return true;

    if (mode_ == Mode::SmallMaterialized) {
      out->pairs.resize(t);
      for (u32 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(W_);
        out->pairs[i] = small_pairs_[static_cast<usize>(idx)];
      }
      return true;
    }

    // Large branch: Phase2+Phase3 (same as Sampling variant).
    detail::SlotPlan2D plan;
    {
      auto scoped2 = phases ? phases->Scoped("phase2_plan") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!detail::BuildSlotPlan2D(t, rng, w_total_, w_a_, w_b_, &plan, err)) {
        if (err && err->empty()) *err = "OursAdaptiveBaseline::Sample: failed to build slot plan";
        return false;
      }
    }

    out->pairs.resize(t);

    {
      auto scoped3 = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

      ctx_.ResetActive();

      auto& ar = ctx_.active_r();
      auto& as = ctx_.active_s();

      const auto& events = ctx_.events();
      const auto& sid_of_pos = ctx_.start_id_of_event();

      const auto& ylo_r = ctx_.ylo_rank_r();
      const auto& yhi_r = ctx_.yhi_lb_rank_r();
      const auto& ylo_s = ctx_.ylo_rank_s();
      const auto& yhi_s = ctx_.yhi_lb_rank_s();

      const auto* ds = ctx_.dataset();
      SJS_DASSERT(ds != nullptr);

      std::vector<u32> sampled;

      for (usize pos = 0; pos < events.size(); ++pos) {
        const auto& ev = events[pos];
        const u32 handle = static_cast<u32>(ev.index);

        if (ev.kind == join::EventKind::End) {
          if (ev.side == join::Side::R) ar.Erase(handle);
          else as.Erase(handle);
          continue;
        }

        const i32 sid_i32 = sid_of_pos[pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 sid = static_cast<u32>(sid_i32);

        const bool q_is_r = (ev.side == join::Side::R);
        const u32 q_ylo = q_is_r ? ylo_r[ev.index] : ylo_s[ev.index];
        const u32 q_yhi = q_is_r ? yhi_r[ev.index] : yhi_s[ev.index];

        const detail::ActiveIndex2D& other = q_is_r ? as : ar;

        // Pattern A slots
        {
          const u32 begin = plan.offset_a[sid];
          const u32 end = plan.offset_a[sid + 1];
          const u32 k = end - begin;
          if (k > 0) {
            sampled.clear();
            const bool ok = other.SampleA(q_ylo, k, rng, &sampled);
            if (!ok || sampled.size() != k) {
              if (err) *err = "OursAdaptiveBaseline::Sample: SampleA failed (inconsistent weights)";
              return false;
            }
            for (u32 i = 0; i < k; ++i) {
              const u32 slot = plan.slots_a[begin + i];
              const u32 oh = sampled[i];
              if (q_is_r) out->pairs[slot] = PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)};
              else out->pairs[slot] = PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)};
            }
          }
        }

        // Pattern B slots
        {
          const u32 begin = plan.offset_b[sid];
          const u32 end = plan.offset_b[sid + 1];
          const u32 k = end - begin;
          if (k > 0) {
            sampled.clear();
            const bool ok = other.SampleB(q_ylo, q_yhi, k, rng, &sampled);
            if (!ok || sampled.size() != k) {
              if (err) *err = "OursAdaptiveBaseline::Sample: SampleB failed (inconsistent weights)";
              return false;
            }
            for (u32 i = 0; i < k; ++i) {
              const u32 slot = plan.slots_b[begin + i];
              const u32 oh = sampled[i];
              if (q_is_r) out->pairs[slot] = PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)};
              else out->pairs[slot] = PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)};
            }
          }
        }

        // Insert q
        if (q_is_r) ar.Insert(handle, q_ylo, q_yhi);
        else as.Insert(handle, q_ylo, q_yhi);
      }

      ctx_.ResetActive();
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    join::PlaneSweepOptions opt;
    opt.axis = 0;
    opt.side_order = join::SideTieBreak::RBeforeS;
    opt.skip_axis_check = true;

    const auto* ds = ctx_.dataset();
    return std::make_unique<detail::PlaneSweepEnumeratorWrapper<Dim, T>>(ds->R, ds->S, opt);
  }

 private:
  enum class Mode {
    Unknown,
    SmallMaterialized,
    LargeSampling,
  };

  bool built_{false};
  detail::Ours2DContext<Dim, T> ctx_;

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
