#pragma once
// sjs/baselines/sirs/adaptive.h
//
// SIRS baseline (Variant::Adaptive).
//
// Adaptive strategy (per "SIRS’21.md"):
//   - Try to enumerate the join pairs while the running total |J| stays <= J*.
//     (We store all pairs in memory when join is small.)
//   - If the running total would exceed J*, stop enumerating, clear stored pairs,
//     and compute remaining degrees via range COUNT to obtain exact |J|.
//   - Sampling:
//       * if join was fully enumerated (|J|<=J*): sample uniformly from stored pairs.
//       * otherwise: sample OUTER by degree (alias) and INNER by range sampling.
//
// This thresholded switch keeps small joins fast/accurate while avoiding
// materializing huge joins.

#include "sjs/baselines/sirs/sampling.h"  // SIRSState + enumerator

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
namespace sirs {

template <int Dim, class T = Scalar>
class SIRSAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::SIRS; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "sirs_adaptive"; }

  void Reset() override {
    st_.Reset();
    enumerated_all_ = false;
    threshold_ = 0;
    all_pairs_.clear();
  }

  bool Build(const DatasetT& ds,
             const Config& cfg,
             PhaseRecorder* phases,
             std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    if (!st_.BuildIndex(ds, cfg, phases, err)) return false;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic
    if (!st_.built) {
      if (err) *err = "SIRSAdaptiveBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "SIRSAdaptiveBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");

    // Threshold selection.
    u64 J_star = cfg.run.j_star;
    if (cfg.run.enum_cap > 0 && cfg.run.enum_cap < J_star) J_star = cfg.run.enum_cap;
    threshold_ = J_star;

    // Clear previous state.
    all_pairs_.clear();
    enumerated_all_ = true;

    st_.weights_valid = false;
    st_.outer_alias.Clear();
    st_.W = 0;
    st_.deg.assign(st_.outer ? st_.outer->Size() : 0, 0ULL);

    if (!st_.outer || !st_.inner) {
      if (err) *err = "SIRSAdaptiveBaseline::Count: internal state missing relations";
      return false;
    }

    const u32 n_outer = static_cast<u32>(st_.outer->Size());
    const u32 n_inner = static_cast<u32>(st_.inner->Size());

    if (n_outer == 0 || n_inner == 0) {
      st_.W = 0;
      enumerated_all_ = true;
      *out = MakeExactCount(0);
      return true;
    }

    // We keep W_total as the exact total join size accumulated per outer.
    u64 W_total = 0;

    // While enumerated_all_ is true, we enumerate pairs and store them, but we
    // stop as soon as we would exceed the threshold.
    for (u32 outer_idx = 0; outer_idx < n_outer; ++outer_idx) {
      const auto& a = st_.outer->boxes[static_cast<usize>(outer_idx)];
      const auto q = st_.QueryForOuter(a);

      if (q.IsEmpty()) {
        st_.deg[static_cast<usize>(outer_idx)] = 0;
        continue;
      }

      if (enumerated_all_) {
        // Enumerate matches for this outer until we hit the threshold.
        u64 local = 0;
        auto cur = st_.kd.MakeCursor(q);
        u32 inner_handle = 0;

        while (cur.Next(&inner_handle, nullptr)) {
          // If adding this pair would exceed the threshold, switch to counting-only.
          if (W_total + local + 1 > J_star) {
            enumerated_all_ = false;
            all_pairs_.clear();

            // Compute the exact degree for this outer without enumerating more.
            const u64 total_local = st_.kd.Count(q);
            st_.deg[static_cast<usize>(outer_idx)] = total_local;

            if (W_total > std::numeric_limits<u64>::max() - total_local) {
              if (err) *err = "SIRSAdaptiveBaseline::Count: |J| overflowed u64";
              return false;
            }
            W_total += total_local;

            // From now on, do count-only for remaining outers.
            break;
          }

          // Still within threshold: store the pair.
          all_pairs_.push_back(st_.MakePair(outer_idx, inner_handle));
          ++local;
        }

        if (enumerated_all_) {
          // We enumerated all matches for this outer.
          st_.deg[static_cast<usize>(outer_idx)] = local;
          if (W_total > std::numeric_limits<u64>::max() - local) {
            if (err) *err = "SIRSAdaptiveBaseline::Count: |J| overflowed u64";
            return false;
          }
          W_total += local;
        } else {
          // Switched in the middle of this outer: remaining outers will be count-only.
          continue;
        }
      }

      if (!enumerated_all_) {
        // Count-only mode.
        const u64 w = st_.kd.Count(q);
        st_.deg[static_cast<usize>(outer_idx)] = w;
        if (W_total > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "SIRSAdaptiveBaseline::Count: |J| overflowed u64";
          return false;
        }
        W_total += w;
      }
    }

    st_.W = W_total;

    if (enumerated_all_) {
      // Sanity: if we enumerated all pairs, we should have stored exactly W pairs.
      if (all_pairs_.size() != static_cast<usize>(st_.W)) {
        // Not fatal (could happen with empty query / numeric oddities), but it
        // indicates a bug in the cursor or query mapping.
        if (err) *err = "SIRSAdaptiveBaseline::Count: stored pairs != computed |J| (bug)";
        return false;
      }
      st_.weights_valid = false;
      st_.outer_alias.Clear();
    } else {
      // Build alias for sampling.
      if (!st_.outer_alias.BuildFromU64(Span<const u64>(st_.deg), err)) return false;
      st_.weights_valid = true;
    }

    *out = MakeExactCount(st_.W);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!st_.built) {
      if (err) *err = "SIRSAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "SIRSAdaptiveBaseline::Sample: null rng/out";
      return false;
    }

    auto scoped = phases ? phases->Scoped("sample") : PhaseRecorder::ScopedPhase(nullptr, "");

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    const u64 t64 = cfg.run.t;
    if (t64 == 0) return true;
    if (t64 > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "SIRSAdaptiveBaseline::Sample: t too large for u32 slots";
      return false;
    }
    const u32 t = static_cast<u32>(t64);

    if (st_.W == 0) return true;

    if (enumerated_all_) {
      if (all_pairs_.empty()) {
        // |J|>0 but all_pairs_ empty should not happen.
        if (err) *err = "SIRSAdaptiveBaseline::Sample: enumerated_all_=true but all_pairs_ is empty";
        return false;
      }

      out->pairs.resize(static_cast<usize>(t));
      const u64 N = static_cast<u64>(all_pairs_.size());
      for (u32 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(N);
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(idx)];
      }
      return true;
    }

    // Sampling branch: require alias/degrees.
    if (!st_.weights_valid) {
      if (err) *err = "SIRSAdaptiveBaseline::Sample: call Count() first";
      return false;
    }

    // Assign sample slots to outers.
    struct Assign {
      u32 outer;
      u32 slot;
    };
    std::vector<Assign> asg;
    asg.reserve(static_cast<usize>(t));
    for (u32 i = 0; i < t; ++i) {
      const u32 outer_idx = static_cast<u32>(st_.outer_alias.Sample(rng));
      asg.push_back(Assign{outer_idx, i});
    }
    std::sort(asg.begin(), asg.end(), [](const Assign& a, const Assign& b) {
      if (a.outer < b.outer) return true;
      if (b.outer < a.outer) return false;
      return a.slot < b.slot;
    });

    out->pairs.resize(static_cast<usize>(t));

    std::vector<u32> inner_handles;
    usize i = 0;
    while (i < asg.size()) {
      const u32 outer_idx = asg[i].outer;
      usize j = i + 1;
      while (j < asg.size() && asg[j].outer == outer_idx) ++j;
      const u32 k = static_cast<u32>(j - i);

      if (outer_idx >= st_.outer->Size()) {
        if (err) *err = "SIRSAdaptiveBaseline::Sample: outer index out of range";
        return false;
      }
      if (st_.deg[static_cast<usize>(outer_idx)] == 0) {
        // Should not happen if alias table is built from deg.
        if (err) *err = "SIRSAdaptiveBaseline::Sample: sampled outer with deg=0";
        return false;
      }

      const auto& a = st_.outer->boxes[static_cast<usize>(outer_idx)];
      const auto q = st_.QueryForOuter(a);

      inner_handles.clear();
      if (!st_.kd.Sample(q, k, rng, &inner_handles, err)) return false;
      if (inner_handles.size() != k) {
        if (err) *err = "SIRSAdaptiveBaseline::Sample: kd.Sample produced wrong number of samples";
        return false;
      }

      for (u32 tslot = 0; tslot < k; ++tslot) {
        const u32 slot = asg[i + tslot].slot;
        out->pairs[static_cast<usize>(slot)] = st_.MakePair(outer_idx, inner_handles[static_cast<usize>(tslot)]);
      }

      i = j;
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!st_.built) {
      if (err) *err = "SIRSAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::SIRSJoinEnumerator<Dim, T>>(&st_);
  }

 private:
  detail::SIRSState<Dim, T> st_;

  bool enumerated_all_{false};
  u64 threshold_{0};
  std::vector<PairId> all_pairs_;
};

}  // namespace sirs
}  // namespace baselines
}  // namespace sjs
