#pragma once
// sjs/baselines/r_tree/adaptive.h
//
// Plane Sweep + Dynamic R-Tree baseline (Variant::Adaptive).
//
// Adaptive strategy (from "R-Tree Baseline v2.0.md"):
//   - Always perform Phase-1 counting via a sweep + R-tree.
//   - While the running total W <= J_star, also enumerate and materialize all
//     pairs into AllPairs (exact join).
//   - If W surpasses J_star, drop the materialized pairs and switch to the
//     Sampling-style Phase-2/3 procedure.
//
// When materialization succeeds (W <= J_star), Sample() simply draws uniform
// i.i.d. samples (with replacement) from AllPairs.

#include "sjs/baselines/r_tree/sampling.h"
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
namespace r_tree {

template <int Dim, class T>
class RTreeAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "RTreeAdaptiveBaseline requires Dim >= 2");
  using DatasetT = Dataset<Dim, T>;
  using ProjBoxT = Box<Dim - 1, T>;

  Method method() const noexcept override { return Method::RTree; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "rtree_adaptive"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;
    mode_ = Mode::CountOnly;

    events_.clear();
    start_id_of_event_.clear();
    w_total_.clear();

    proj_r_.clear();
    proj_s_.clear();

    all_pairs_.clear();

    tree_r_.Clear();
    tree_s_.Clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;
    mode_ = Mode::CountOnly;
    all_pairs_.clear();

    if (ds.R.Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        ds.S.Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RTreeAdaptiveBaseline::Build: relation size exceeds u32 limit";
      return false;
    }

    // Parse R-tree options from cfg.run.extra.
    typename detail::DynamicRTree<Dim - 1, T>::Options opt;
    {
      u32 m = 32;
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_max_children", m);
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_M", m);
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_m", m);
      if (m < 2) m = 2;
      opt.max_children = m;

      u32 minc = 0;
      minc = detail::ExtraU32Or(cfg.run.extra, "rtree_min_children", minc);
      if (minc > 0) opt.min_children = minc;
      opt.ignore_duplicate_insert = true;
    }
    opt_ = opt;

    // Events and start-id mapping.
    events_ = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);
    start_id_of_event_.assign(events_.size(), -1);
    usize start_cnt = 0;
    for (usize i = 0; i < events_.size(); ++i) {
      if (events_[i].kind == join::EventKind::Start) start_id_of_event_[i] = static_cast<i32>(start_cnt++);
    }
    w_total_.assign(start_cnt, 0ULL);

    // Projection caches.
    proj_r_.resize(ds.R.Size());
    for (usize i = 0; i < ds.R.Size(); ++i) proj_r_[i] = detail::ProjectDropFirst<Dim, T>(ds.R.boxes[i]);
    proj_s_.resize(ds.S.Size());
    for (usize i = 0; i < ds.S.Size(); ++i) proj_s_[i] = detail::ProjectDropFirst<Dim, T>(ds.S.boxes[i]);

    // Init trees (capacity only); cleared per sweep.
    tree_r_.Init(static_cast<u32>(ds.R.Size()), opt_);
    tree_s_.Init(static_cast<u32>(ds.S.Size()), opt_);

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;
    if (!built_ || !ds_) {
      if (err) *err = "RTreeAdaptiveBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RTreeAdaptiveBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_adaptive_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    // Determine enumeration threshold J_star.
    u64 J_star = cfg.run.j_star;
    if (cfg.run.enum_cap > 0) J_star = std::min<u64>(J_star, cfg.run.enum_cap);

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    W_ = 0;
    weights_valid_ = false;
    mode_ = Mode::Enumerate;
    all_pairs_.clear();

    tree_r_.Clear();
    tree_s_.Clear();

    std::vector<u32> hits;
    hits.reserve(256);

    u64 W = 0;
    for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
      const join::Event& e = events_[ev_pos];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) {
          (void)tree_r_.Remove(static_cast<u32>(e.index));
        } else {
          (void)tree_s_.Remove(static_cast<u32>(e.index));
        }
        continue;
      }

      // START
      const i32 sid_i32 = start_id_of_event_[ev_pos];
      SJS_DASSERT(sid_i32 >= 0);
      const u32 sid = static_cast<u32>(sid_i32);

      if (e.side == join::Side::R) {
        const ProjBoxT& q = proj_r_[e.index];
        const u64 w = tree_s_.CountIntersect(q);
        w_total_[sid] = w;
        W += w;

        if (mode_ == Mode::Enumerate) {
          if (W <= J_star) {
            hits.clear();
            tree_s_.ReportIntersect(q, &hits);
            // Safety: ReportIntersect should match CountIntersect.
            // (We don't abort on mismatch to avoid killing long runs; still deterministic.)
            const Id rid = ds_->R.GetId(e.index);
            for (u32 sidx : hits) {
              all_pairs_.push_back(PairId{rid, ds_->S.GetId(static_cast<usize>(sidx))});
            }
          } else {
            // Trigger switch.
            mode_ = Mode::CountOnly;
            all_pairs_.clear();
          }
        }

        (void)tree_r_.Insert(static_cast<u32>(e.index), q);
      } else {
        const ProjBoxT& q = proj_s_[e.index];
        const u64 w = tree_r_.CountIntersect(q);
        w_total_[sid] = w;
        W += w;

        if (mode_ == Mode::Enumerate) {
          if (W <= J_star) {
            hits.clear();
            tree_r_.ReportIntersect(q, &hits);
            const Id sid_id = ds_->S.GetId(e.index);
            for (u32 ridx : hits) {
              all_pairs_.push_back(PairId{ds_->R.GetId(static_cast<usize>(ridx)), sid_id});
            }
          } else {
            mode_ = Mode::CountOnly;
            all_pairs_.clear();
          }
        }

        (void)tree_s_.Insert(static_cast<u32>(e.index), q);
      }
    }

    W_ = W;
    weights_valid_ = true;

    if (mode_ == Mode::Enumerate) {
      // If we stayed in ENUMERATE mode, AllPairs should contain the whole join.
      // We do not hard-assert equality because W could be 0, and some queries may have
      // empty boxes ignored. Still, this is a helpful sanity check in debug builds.
      SJS_DASSERT(static_cast<u64>(all_pairs_.size()) == W_);
    }

    *out = baselines::MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "RTreeAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "RTreeAdaptiveBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "RTreeAdaptiveBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    // Ensure weights exist.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, rng, &tmp, phases, err)) return false;
    }

    if (W_ == 0) {
      // Empty join.
      return true;
    }

    if (mode_ == Mode::Enumerate) {
      auto scoped = phases ? phases->Scoped("sample_from_materialized") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (all_pairs_.empty()) {
        // Should not happen when W_>0.
        if (err) *err = "RTreeAdaptiveBaseline::Sample: materialized mode but AllPairs is empty";
        return false;
      }
      out->pairs.resize(static_cast<usize>(t));
      const u64 N = static_cast<u64>(all_pairs_.size());
      for (u64 i = 0; i < t; ++i) {
        const u64 j = rng->UniformU64(N);
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(j)];
      }
      return true;
    }

    // COUNT_ONLY mode: execute Baseline v2.0 "Sampling" Phase2+Phase3.
    auto scoped = phases ? phases->Scoped("phase2+3_sampling") : PhaseRecorder::ScopedPhase(nullptr, "");

    // This baseline stores slot indices in u32 (and DynamicRTree::SampleIntersect takes u32),
    // so we require t <= 2^32-1 (practically always true).
    if (t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RTreeAdaptiveBaseline::Sample: t exceeds u32 limit";
      return false;
    }
    const u32 t_u32 = static_cast<u32>(t);

    const u32 num_starts = static_cast<u32>(w_total_.size());

    // Build an alias distribution over START events with positive weight only.
    std::vector<u32> pos_eids;
    std::vector<u64> pos_w;
    pos_eids.reserve(num_starts);
    pos_w.reserve(num_starts);
    for (u32 eid = 0; eid < num_starts; ++eid) {
      const u64 w = w_total_[static_cast<usize>(eid)];
      if (w > 0) {
        pos_eids.push_back(eid);
        pos_w.push_back(w);
      }
    }
    if (pos_eids.empty()) {
      if (err) *err = "RTreeAdaptiveBaseline::Sample: no positive-weight START events (W>0)";
      return false;
    }

    // --------------------------
    // Phase 2: event alias + slots[start_id].append(j) (Baseline v2.0 §3.2 Phase2)
    // --------------------------
    // Compact CSR layout:
    //   - slot_offsets[eid]..slot_offsets[eid+1] are the output positions assigned to START eid.
    //   - slots_data stores all output positions contiguously.
    std::vector<u32> slot_offsets;
    std::vector<u32> slots_data;
    {
      sampling::AliasTable alias;
      if (!alias.BuildFromU64(Span<const u64>(pos_w), err)) {
        if (err && err->empty()) *err = "RTreeAdaptiveBaseline::Sample: alias build failed";
        return false;
      }

      std::vector<u32> slot_counts(static_cast<usize>(num_starts), 0U);
      std::vector<u32> chosen_eid(static_cast<usize>(t_u32), 0U);

      for (u32 j = 0; j < t_u32; ++j) {
        const u32 idx = static_cast<u32>(alias.Sample(rng));  // index into pos_eids
        const u32 eid = pos_eids[static_cast<usize>(idx)];
        chosen_eid[static_cast<usize>(j)] = eid;
        slot_counts[static_cast<usize>(eid)]++;
      }

      slot_offsets.assign(static_cast<usize>(num_starts + 1), 0U);
      for (u32 eid = 0; eid < num_starts; ++eid) {
        slot_offsets[static_cast<usize>(eid + 1)] =
            slot_offsets[static_cast<usize>(eid)] + slot_counts[static_cast<usize>(eid)];
      }
      SJS_DASSERT(slot_offsets.back() == t_u32);

      slots_data.assign(static_cast<usize>(t_u32), 0U);
      std::vector<u32> write_pos(slot_offsets.begin(), slot_offsets.end() - 1);
      for (u32 j = 0; j < t_u32; ++j) {
        const u32 eid = chosen_eid[static_cast<usize>(j)];
        const u32 pos = write_pos[static_cast<usize>(eid)]++;
        slots_data[static_cast<usize>(pos)] = j;
      }
    }

    out->pairs.assign(static_cast<usize>(t_u32), PairId{});

    // --------------------------
    // Phase 3: second sweep + local sampling (Baseline v2.0 §3.2 Phase3)
    // --------------------------
    tree_r_.Clear();
    tree_s_.Clear();

    std::vector<u32> picked;
    picked.reserve(256);

    for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
      const join::Event& e = events_[ev_pos];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) {
          (void)tree_r_.Remove(static_cast<u32>(e.index));
        } else {
          (void)tree_s_.Remove(static_cast<u32>(e.index));
        }
        continue;
      }

      const i32 sid_i32 = start_id_of_event_[ev_pos];
      SJS_DASSERT(sid_i32 >= 0);
      const u32 eid = static_cast<u32>(sid_i32);

      const u32 begin = slot_offsets[static_cast<usize>(eid)];
      const u32 end = slot_offsets[static_cast<usize>(eid + 1)];
      const u32 need = end - begin;

      if (e.side == join::Side::R) {
        const ProjBoxT& q = proj_r_[e.index];
        if (need > 0) {
          picked.clear();
          if (!tree_s_.SampleIntersect(q, need, rng, &picked)) {
            if (err) *err = "RTreeAdaptiveBaseline::Sample: SampleIntersect failed (unexpected empty K_e)";
            return false;
          }
          const Id rid = ds_->R.GetId(e.index);
          for (u32 u = 0; u < need; ++u) {
            const u32 sidx = picked[static_cast<usize>(u)];
            const u32 slot = slots_data[static_cast<usize>(begin + u)];
            out->pairs[static_cast<usize>(slot)] = PairId{rid, ds_->S.GetId(static_cast<usize>(sidx))};
          }
        }
        (void)tree_r_.Insert(static_cast<u32>(e.index), q);
      } else {
        const ProjBoxT& q = proj_s_[e.index];
        if (need > 0) {
          picked.clear();
          if (!tree_r_.SampleIntersect(q, need, rng, &picked)) {
            if (err) *err = "RTreeAdaptiveBaseline::Sample: SampleIntersect failed (unexpected empty K_e)";
            return false;
          }
          const Id sid_id = ds_->S.GetId(e.index);
          for (u32 u = 0; u < need; ++u) {
            const u32 ridx = picked[static_cast<usize>(u)];
            const u32 slot = slots_data[static_cast<usize>(begin + u)];
            out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(ridx)), sid_id};
          }
        }
        (void)tree_s_.Insert(static_cast<u32>(e.index), q);
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
      if (err) *err = "RTreeAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RTreeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, opt_, /*axis=*/0,
                                                                 join::SideTieBreak::RBeforeS);
  }

 private:
  enum class Mode : u8 { Enumerate = 0, CountOnly = 1 };

  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  std::vector<join::Event> events_;
  std::vector<i32> start_id_of_event_;

  std::vector<u64> w_total_;
  u64 W_ = 0;
  bool weights_valid_ = false;

  // Projection caches.
  std::vector<ProjBoxT> proj_r_;
  std::vector<ProjBoxT> proj_s_;

  // Active-set indexes.
  detail::DynamicRTree<Dim - 1, T> tree_r_;
  detail::DynamicRTree<Dim - 1, T> tree_s_;

  typename detail::DynamicRTree<Dim - 1, T>::Options opt_{};

  // Enumerate mode storage.
  Mode mode_ = Mode::CountOnly;
  std::vector<PairId> all_pairs_;
};

}  // namespace r_tree
}  // namespace baselines
}  // namespace sjs
