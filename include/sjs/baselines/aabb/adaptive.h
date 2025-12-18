#pragma once
// sjs/baselines/aabb/adaptive.h
//
// Plane Sweep + Dynamic AABB-Tree baseline (Variant::Adaptive).
//
// This follows the "Adaptive+Sampling" variant in "AABB-Tree Baseline v2.0.md":
//   - Phase1: sweep and compute weights w_e exactly using CountIntersect.
//            While W <= J*, also enumerate pairs (via ReportIntersect) into
//            AllPairs; if W exceeds J*, discard and switch to COUNT_ONLY.
//   - If W <= J*: sample directly from AllPairs (uniform i.i.d., with replacement).
//   - Else: run Phase2+Phase3 (alias + two-pass sampling) identical to the
//           Sampling variant.
//
// NOTE: The project also provides a baseline-agnostic adaptive runner that can
// do pilot enumeration before falling back to Sampling. This class is still
// useful standalone and for matching the paper's description.

#include "sjs/baselines/aabb/sampling.h"

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
namespace aabb {

template <int Dim, class T>
class AABBAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "AABBAdaptiveBaseline requires Dim >= 2");
  using DatasetT = Dataset<Dim, T>;
  using ProjBoxT = Box<Dim - 1, T>;

  Method method() const noexcept override { return Method::AABB; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "aabb_adaptive"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;
    mode_ = Mode::Unknown;

    events_.clear();
    start_id_of_event_.clear();
    w_total_.clear();

    proj_r_.clear();
    proj_s_.clear();
    handle_r_.clear();
    handle_s_.clear();
    tree_r_.Clear();
    tree_s_.Clear();

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

    if (ds.R.Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        ds.S.Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "AABBAdaptiveBaseline::Build: relation size exceeds u32 limit";
      return false;
    }

    events_ = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);
    start_id_of_event_.assign(events_.size(), -1);
    usize start_cnt = 0;
    for (usize i = 0; i < events_.size(); ++i) {
      if (events_[i].kind == join::EventKind::Start) {
        start_id_of_event_[i] = static_cast<i32>(start_cnt);
        ++start_cnt;
      }
    }
    w_total_.assign(start_cnt, 0ULL);

    proj_r_.resize(ds.R.Size());
    for (usize i = 0; i < ds.R.Size(); ++i) {
      proj_r_[i] = detail::ProjectDropFirst<Dim, T>(ds.R.boxes[i]);
    }
    proj_s_.resize(ds.S.Size());
    for (usize i = 0; i < ds.S.Size(); ++i) {
      proj_s_[i] = detail::ProjectDropFirst<Dim, T>(ds.S.boxes[i]);
    }

    handle_r_.assign(ds.R.Size(), kNull);
    handle_s_.assign(ds.S.Size(), kNull);

    tree_r_.Clear();
    tree_s_.Clear();
    tree_r_.ReserveLeaves(ds.R.Size());
    tree_s_.ReserveLeaves(ds.S.Size());

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
      if (err) *err = "AABBAdaptiveBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "AABBAdaptiveBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 J_star_cfg = cfg.run.j_star;
    const u64 cap_cfg = cfg.run.enum_cap;
    const u64 J_star = (cap_cfg == 0) ? J_star_cfg : std::min(J_star_cfg, cap_cfg);

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    W_ = 0;
    weights_valid_ = false;
    mode_ = Mode::Enumerate;  // optimistic; may switch to CountOnly
    all_pairs_.clear();
    if (J_star > 0) {
      all_pairs_.reserve(static_cast<usize>(std::min<u64>(J_star, 1'000'000ULL)));
    }

    tree_r_.Clear();
    tree_s_.Clear();
    std::fill(handle_r_.begin(), handle_r_.end(), kNull);
    std::fill(handle_s_.begin(), handle_s_.end(), kNull);

    u64 W = 0;
    std::vector<u32> tmp;

    for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
      const join::Event& e = events_[ev_pos];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) {
          const int h = handle_r_[e.index];
          if (h != kNull) tree_r_.Remove(h);
          handle_r_[e.index] = kNull;
        } else {
          const int h = handle_s_[e.index];
          if (h != kNull) tree_s_.Remove(h);
          handle_s_[e.index] = kNull;
        }
        continue;
      }

      const i32 sid_i32 = start_id_of_event_[ev_pos];
      SJS_DASSERT(sid_i32 >= 0);
      const usize sid = static_cast<usize>(sid_i32);

      if (e.side == join::Side::R) {
        const ProjBoxT& q = proj_r_[e.index];
        const u64 w = tree_s_.CountIntersect(q);
        w_total_[sid] = w;
        W += w;

        // "Count first, then threshold" (as emphasized in the markdown).
        if (mode_ == Mode::Enumerate && J_star > 0) {
          if (W <= J_star) {
            tmp.clear();
            tree_s_.ReportIntersect(q, &tmp);
            const Id rid = ds_->R.GetId(e.index);
            for (u32 s_idx : tmp) {
              all_pairs_.push_back(PairId{rid, ds_->S.GetId(static_cast<usize>(s_idx))});
            }
          } else {
            mode_ = Mode::CountOnly;
            all_pairs_.clear();
            all_pairs_.shrink_to_fit();
          }
        } else if (mode_ == Mode::Enumerate && J_star == 0) {
          // J_star==0 means "never enumerate".
          mode_ = Mode::CountOnly;
        }

        handle_r_[e.index] = tree_r_.Insert(static_cast<u32>(e.index), q);
      } else {
        const ProjBoxT& q = proj_s_[e.index];
        const u64 w = tree_r_.CountIntersect(q);
        w_total_[sid] = w;
        W += w;

        if (mode_ == Mode::Enumerate && J_star > 0) {
          if (W <= J_star) {
            tmp.clear();
            tree_r_.ReportIntersect(q, &tmp);
            const Id sid_id = ds_->S.GetId(e.index);
            for (u32 r_idx : tmp) {
              all_pairs_.push_back(PairId{ds_->R.GetId(static_cast<usize>(r_idx)), sid_id});
            }
          } else {
            mode_ = Mode::CountOnly;
            all_pairs_.clear();
            all_pairs_.shrink_to_fit();
          }
        } else if (mode_ == Mode::Enumerate && J_star == 0) {
          mode_ = Mode::CountOnly;
        }

        handle_s_[e.index] = tree_s_.Insert(static_cast<u32>(e.index), q);
      }
    }

    // If we never exceeded threshold and enumerated, we must have AllPairs == W.
    if (mode_ == Mode::Enumerate) {
      // Note: in rare cases we may have skipped enumeration due to J_star==0.
      if (J_star > 0 && all_pairs_.size() != static_cast<usize>(W)) {
        // Deterministic sanity check; do not fail hard in release.
        SJS_DASSERT(all_pairs_.size() == static_cast<usize>(W));
      }
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
      if (err) *err = "AABBAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "AABBAdaptiveBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "AABBAdaptiveBaseline::Sample: out is null";
      return false;
    }

    const u32 t = cfg.run.t;
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

    // Small join branch: sample from explicit AllPairs.
    if (mode_ == Mode::Enumerate && !all_pairs_.empty()) {
      auto scoped = phases ? phases->Scoped("sample_from_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");
      out->pairs.resize(static_cast<usize>(t));
      for (u32 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(static_cast<u64>(all_pairs_.size()));
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(idx)];
      }
      return true;
    }

    // Large join branch: Phase2 + Phase3 (same as Sampling variant).
    struct Assignment {
      u32 eid;
      u32 slot;
    };
    auto assign_less = [](const Assignment& a, const Assignment& b) noexcept {
      if (a.eid < b.eid) return true;
      if (b.eid < a.eid) return false;
      return a.slot < b.slot;
    };

    std::vector<Assignment> assign;
    {
      auto scoped = phases ? phases->Scoped("phase2_assign") : PhaseRecorder::ScopedPhase(nullptr, "");
      // Build alias only on positive-weight START events (w_e > 0), as required by the v2.0 design.
      std::vector<u32> pos_eids;
      std::vector<u64> pos_w;
      pos_eids.reserve(w_total_.size());
      pos_w.reserve(w_total_.size());
      for (u32 eid = 0; eid < static_cast<u32>(w_total_.size()); ++eid) {
        const u64 w = w_total_[static_cast<usize>(eid)];
        if (w == 0) continue;
        pos_eids.push_back(eid);
        pos_w.push_back(w);
      }
      if (pos_w.empty()) {
        if (err) *err = "AABBAdaptiveBaseline::Sample: internal error (W>0 but no positive-weight events)";
        return false;
      }

      sampling::AliasTable alias;
      if (!alias.BuildFromU64(Span<const u64>(pos_w), err)) {
        if (err && err->empty()) *err = "AABBAdaptiveBaseline::Sample: failed to build alias table";
        return false;
      }

      assign.reserve(t);
      for (u32 j = 0; j < t; ++j) {
        const u32 idx = static_cast<u32>(alias.Sample(rng));
        const u32 eid = pos_eids[static_cast<usize>(idx)];
        assign.push_back(Assignment{eid, j});
      }

      std::sort(assign.begin(), assign.end(), assign_less);
    }

    out->pairs.assign(static_cast<usize>(t), PairId{});

    {
      auto scoped = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");
      tree_r_.Clear();
      tree_s_.Clear();
      std::fill(handle_r_.begin(), handle_r_.end(), kNull);
      std::fill(handle_s_.begin(), handle_s_.end(), kNull);

      usize a_ptr = 0;
      std::vector<u32> picked;

      for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
        const join::Event& e = events_[ev_pos];
        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) {
            const int h = handle_r_[e.index];
            if (h != kNull) tree_r_.Remove(h);
            handle_r_[e.index] = kNull;
          } else {
            const int h = handle_s_[e.index];
            if (h != kNull) tree_s_.Remove(h);
            handle_s_[e.index] = kNull;
          }
          continue;
        }

        const i32 sid_i32 = start_id_of_event_[ev_pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 eid = static_cast<u32>(sid_i32);

        const usize begin = a_ptr;
        while (a_ptr < assign.size() && assign[a_ptr].eid == eid) {
          ++a_ptr;
        }
        const u32 t_e = static_cast<u32>(a_ptr - begin);

        if (t_e > 0) {
          picked.clear();
          if (e.side == join::Side::R) {
            const ProjBoxT& q = proj_r_[e.index];
            if (!tree_s_.SampleIntersect(q, t_e, rng, &picked)) {
              if (err) *err = "AABBAdaptiveBaseline::Sample: SampleIntersect returned empty unexpectedly";
              return false;
            }
            const Id rid = ds_->R.GetId(e.index);
            for (u32 u = 0; u < t_e; ++u) {
              const u32 s_idx = picked[static_cast<usize>(u)];
              out->pairs[static_cast<usize>(assign[begin + u].slot)] = PairId{rid, ds_->S.GetId(s_idx)};
            }
          } else {
            const ProjBoxT& q = proj_s_[e.index];
            if (!tree_r_.SampleIntersect(q, t_e, rng, &picked)) {
              if (err) *err = "AABBAdaptiveBaseline::Sample: SampleIntersect returned empty unexpectedly";
              return false;
            }
            const Id sid_id = ds_->S.GetId(e.index);
            for (u32 u = 0; u < t_e; ++u) {
              const u32 r_idx = picked[static_cast<usize>(u)];
              out->pairs[static_cast<usize>(assign[begin + u].slot)] = PairId{ds_->R.GetId(r_idx), sid_id};
            }
          }
        }

        // Insert START.
        if (e.side == join::Side::R) {
          const ProjBoxT& q = proj_r_[e.index];
          handle_r_[e.index] = tree_r_.Insert(static_cast<u32>(e.index), q);
        } else {
          const ProjBoxT& q = proj_s_[e.index];
          handle_s_[e.index] = tree_s_.Insert(static_cast<u32>(e.index), q);
        }
      }

      if (a_ptr != assign.size()) {
        if (err) *err = "AABBAdaptiveBaseline::Sample: internal error (did not consume all assignments)";
        return false;
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
      if (err) *err = "AABBAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::AABBJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                join::SideTieBreak::RBeforeS);
  }

 private:
  static constexpr int kNull = -1;

  enum class Mode : u8 {
    Unknown = 0,
    Enumerate = 1,
    CountOnly = 2,
  };

  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  std::vector<join::Event> events_;
  std::vector<i32> start_id_of_event_;
  std::vector<u64> w_total_;

  std::vector<ProjBoxT> proj_r_;
  std::vector<ProjBoxT> proj_s_;

  detail::DynamicAABBTree<Dim - 1, T> tree_r_;
  detail::DynamicAABBTree<Dim - 1, T> tree_s_;

  std::vector<int> handle_r_;
  std::vector<int> handle_s_;

  // Phase-1 outputs/state.
  u64 W_ = 0;
  bool weights_valid_ = false;
  Mode mode_ = Mode::Unknown;

  // Only used when mode_==Enumerate and W <= J*.
  std::vector<PairId> all_pairs_;
};

}  // namespace aabb
}  // namespace baselines
}  // namespace sjs
