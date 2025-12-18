#pragma once
// sjs/baselines/r_tree/enum_sampling.h
//
// Plane Sweep + Dynamic R-Tree baseline (Variant::EnumSampling).
//
// Baseline v2.0 alignment:
//   Version A: Enumerate+Sampling (explicitly enumerate the full join J,
//   materialize all intersecting pairs, then sample i.i.d. uniformly with
//   replacement from that array).
//
// This variant is only appropriate when |J| is small enough to fit in memory.
// For large joins, use Variant::Sampling (two-pass) or Variant::Adaptive.

#include "sjs/baselines/r_tree/sampling.h"

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
class RTreeEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "RTreeEnumSamplingBaseline requires Dim >= 2");
  using DatasetT = Dataset<Dim, T>;
  using ProjBoxT = Box<Dim - 1, T>;

  Method method() const noexcept override { return Method::RTree; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "rtree_enum+sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    pairs_valid_ = false;

    events_.clear();
    proj_r_.clear();
    proj_s_.clear();
    pairs_.clear();

    tree_r_.Clear();
    tree_s_.Clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = false;
    pairs_valid_ = false;
    pairs_.clear();

    if (ds.R.Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        ds.S.Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RTreeEnumSamplingBaseline::Build: relation size exceeds u32 limit";
      return false;
    }

    // Parse R-tree options (same keys as other R-tree baselines).
    typename detail::DynamicRTree<Dim - 1, T>::Options opt;
    {
      u32 m = 32;
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_max_children", m);
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_M", m);
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_m", m);  // legacy alias
      if (m < 2) m = 2;
      opt.max_children = m;

      u32 minc = 0;
      minc = detail::ExtraU32Or(cfg.run.extra, "rtree_min_children", minc);
      if (minc > 0) opt.min_children = minc;

      opt.ignore_duplicate_insert = true;
    }
    opt_ = opt;

    // Build sweep events (axis 0 = x1) with deterministic tie-breaking.
    events_ = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);

    // Projection caches (drop sweep axis).
    proj_r_.resize(ds.R.Size());
    for (usize i = 0; i < ds.R.Size(); ++i) proj_r_[i] = detail::ProjectDropFirst<Dim, T>(ds.R.boxes[i]);
    proj_s_.resize(ds.S.Size());
    for (usize i = 0; i < ds.S.Size(); ++i) proj_s_[i] = detail::ProjectDropFirst<Dim, T>(ds.S.boxes[i]);

    // Init trees (capacity only); cleared before enumeration.
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
    (void)cfg;
    (void)rng;
    if (!built_ || !ds_) {
      if (err) *err = "RTreeEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RTreeEnumSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("enumerate_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    if (!EnsureEnumerated(phases, err)) return false;
    *out = baselines::MakeExactCount(static_cast<u64>(pairs_.size()));
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "RTreeEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "RTreeEnumSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "RTreeEnumSamplingBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    auto scoped = phases ? phases->Scoped("enumerate+array_sampling") : PhaseRecorder::ScopedPhase(nullptr, "");

    if (!EnsureEnumerated(phases, err)) return false;
    const u64 W = static_cast<u64>(pairs_.size());
    if (W == 0) {
      // Empty join.
      return true;
    }

    out->pairs.resize(static_cast<usize>(t));
    for (u64 i = 0; i < t; ++i) {
      const u64 j = rng->UniformU64(W);
      out->pairs[static_cast<usize>(i)] = pairs_[static_cast<usize>(j)];
    }
    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !ds_) {
      if (err) *err = "RTreeEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RTreeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, opt_, /*axis=*/0,
                                                                 join::SideTieBreak::RBeforeS);
  }

 private:
  bool EnsureEnumerated(PhaseRecorder* phases, std::string* err) {
    (void)err;  // unused parameter (interface requirement)
    if (pairs_valid_) return true;
    if (!built_ || !ds_) return false;

    auto scoped = phases ? phases->Scoped("enumerate_pairs") : PhaseRecorder::ScopedPhase(nullptr, "");

    pairs_.clear();
    tree_r_.Clear();
    tree_s_.Clear();

    std::vector<u32> hits;
    hits.reserve(256);

    // Version A (Baseline v2.0 §3.1): one sweep, report all pairs at each START.
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
      if (e.side == join::Side::R) {
        const ProjBoxT& q = proj_r_[e.index];
        hits.clear();
        tree_s_.ReportIntersect(q, &hits);

        const Id rid = ds_->R.GetId(e.index);
        for (u32 sidx : hits) {
          pairs_.push_back(PairId{rid, ds_->S.GetId(static_cast<usize>(sidx))});
        }

        (void)tree_r_.Insert(static_cast<u32>(e.index), q);
      } else {
        const ProjBoxT& q = proj_s_[e.index];
        hits.clear();
        tree_r_.ReportIntersect(q, &hits);

        const Id sid = ds_->S.GetId(e.index);
        for (u32 ridx : hits) {
          pairs_.push_back(PairId{ds_->R.GetId(static_cast<usize>(ridx)), sid});
        }

        (void)tree_s_.Insert(static_cast<u32>(e.index), q);
      }
    }

    pairs_valid_ = true;
    return true;
  }

  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  // Deterministic sweep events (axis 0).
  std::vector<join::Event> events_;

  // Projection caches.
  std::vector<ProjBoxT> proj_r_;
  std::vector<ProjBoxT> proj_s_;

  // Active-set R-trees.
  detail::DynamicRTree<Dim - 1, T> tree_r_;
  detail::DynamicRTree<Dim - 1, T> tree_s_;

  // R-tree options.
  typename detail::DynamicRTree<Dim - 1, T>::Options opt_{};

  // Materialized join pairs (Version A).
  bool pairs_valid_ = false;
  std::vector<PairId> pairs_;
};

}  // namespace r_tree
}  // namespace baselines
}  // namespace sjs
