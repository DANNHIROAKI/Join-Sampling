#pragma once
// sjs/baselines/kd_tree/sampling.h
//
// KdTree baseline (HighDims) — Sampling variant (KDS-Sampling).
//
// This implements the Chapter-5 baseline from SJS-HighDims.md:
//   - Embed each box s in S into a 2d-dimensional point p(s) = (lo(s), hi(s)).
//   - For each r in R, the match set is the orthogonal range query
//       Q(r) = Π_i (-∞, r.hi[i]) × Π_i (r.lo[i], +∞)
//     on p(S), using strict inequalities (half-open boxes => strict overlap).
//   - Count: w(r) = COUNT(Q(r)), and |J| = Σ_r w(r).
//   - Sample t i.i.d pairs:
//       1) sample r with Pr(r)=w(r)/|J|
//       2) sample s uniformly from Q(r)
//
// This is *not* a sweep-based event-block method.

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/embedding.h"
#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <array>
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

namespace detail {

template <class T>
struct OpenInterval {
  bool has_lo = false;
  T lo{};
  bool has_hi = false;
  T hi{};
};

// 2d-dimensional query for intersection in all dimensions.
template <int Dim, class T>
struct Query2D {
  static constexpr int D = 2 * Dim;
  std::array<OpenInterval<T>, D> I{};
};

// Build Q(r): (-inf, r.hi[i]) for lo-dims, (r.lo[i], inf) for hi-dims.
template <int Dim, class T>
inline Query2D<Dim, T> MakeQueryForLeftBox(const Box<Dim, T>& r) {
  Query2D<Dim, T> q;
  for (int i = 0; i < Dim; ++i) {
    // lo(s)[i] < hi(r)[i]
    q.I[i].has_lo = false;
    q.I[i].has_hi = true;
    q.I[i].hi = r.hi.v[i];

    // hi(s)[i] > lo(r)[i]
    q.I[Dim + i].has_lo = true;
    q.I[Dim + i].lo = r.lo.v[i];
    q.I[Dim + i].has_hi = false;
  }
  return q;
}

template <int Dim, class T>
class KdTreeIndex {
 public:
  static_assert(Dim >= 1, "KdTreeIndex requires Dim >= 1");
  static constexpr int D = 2 * Dim;

  using BoxT = Box<Dim, T>;
  using QryT = Query2D<Dim, T>;
  using PointT = Point<D, T>;

  struct Pt {
    PointT p{};     // embedded coordinates
    u32 handle = 0; // S-handle (index in relation)
  };

  struct Node {
    u32 l = 0;  // inclusive
    u32 r = 0;  // exclusive
    int left = -1;
    int right = -1;
    bool leaf = false;
    std::array<T, D> mn{};
    std::array<T, D> mx{};
  };

  void Reset() {
    pts_.clear();
    nodes_.clear();
    root_ = -1;
  }

  bool BuildFromRelation(const Relation<Dim, T>& rel, u32 leaf_size = 32U) {
    Reset();
    leaf_size_ = (leaf_size == 0U) ? 1U : leaf_size;

    const usize n = rel.Size();
    pts_.resize(n);
    for (usize i = 0; i < n; ++i) {
      const BoxT& b = rel.boxes[i];
      pts_[i].p = EmbedLowerUpper<Dim, T>(b);
      pts_[i].handle = static_cast<u32>(i);
    }

    if (n == 0) {
      root_ = -1;
      return true;
    }

    nodes_.reserve(2 * n);
    root_ = BuildNode(0U, static_cast<u32>(n), 0U);
    return true;
  }

  u64 Count(const QryT& q) const {
    if (root_ < 0) return 0ULL;
    return CountNode(root_, q);
  }

  void Report(const QryT& q, std::vector<u32>* out) const {
    if (!out) return;
    out->clear();
    if (root_ < 0) return;
    ReportNode(root_, q, out);
  }

  bool Sample(const QryT& q, u64 k, Rng* rng, std::vector<u32>* out) const {
    if (!out) return false;
    out->clear();
    if (k == 0) return true;
    if (!rng) return false;
    if (root_ < 0) return true;

    // Collect canonical blocks for sampling.
    struct Block {
      u64 w = 0;
      bool contained = false;
      u32 l = 0;
      u32 r = 0;
      u32 hit_off = 0;
      u32 hit_len = 0;
    };

    std::vector<Block> blocks;
    blocks.reserve(64);
    std::vector<u32> hits;
    hits.reserve(256);

    CollectBlocks(root_, q, &blocks, &hits);

    u64 total = 0;
    for (const auto& b : blocks) {
      if (b.w == 0) continue;
      if (total > (std::numeric_limits<u64>::max() - b.w)) {
        return false;
      }
      total += b.w;
    }
    if (total == 0) return true;

    std::vector<u64> pref;
    pref.resize(blocks.size());
    u64 run = 0;
    for (usize i = 0; i < blocks.size(); ++i) {
      run += blocks[i].w;
      pref[i] = run;
    }

    out->reserve(static_cast<usize>(k));
    for (u64 i = 0; i < k; ++i) {
      const u64 x = rng->UniformU64(total) + 1ULL;
      const auto it = std::lower_bound(pref.begin(), pref.end(), x);
      if (it == pref.end()) return false;
      const usize bi = static_cast<usize>(it - pref.begin());
      const Block& b = blocks[bi];
      if (b.w == 0) return false;

      u32 pick = 0;
      if (b.contained) {
        const u32 len = b.r - b.l;
        if (len == 0) return false;
        const u32 off = static_cast<u32>(rng->UniformU64(len));
        pick = b.l + off;
      } else {
        if (b.hit_len == 0) return false;
        const u32 off = static_cast<u32>(rng->UniformU64(b.hit_len));
        pick = hits[static_cast<usize>(b.hit_off + off)];
      }
      if (pick >= pts_.size()) return false;
      out->push_back(pts_[pick].handle);
    }

    return true;
  }

 private:
  // Query predicate for a point.
  static bool PointInQuery(const Pt& pt, const QryT& q) {
    for (int d = 0; d < D; ++d) {
      const T x = pt.p.v[d];
      const auto& I = q.I[d];
      if (I.has_hi && !(x < I.hi)) return false;  // strict
      if (I.has_lo && !(x > I.lo)) return false;  // strict
    }
    return true;
  }

  enum class Rel { Disjoint, Partial, Contained };

  static Rel RelateMBR(const Node& n, const QryT& q) {
    bool all_contained = true;

    for (int d = 0; d < D; ++d) {
      const auto& I = q.I[d];
      const T mn = n.mn[d];
      const T mx = n.mx[d];

      if (I.has_lo && I.has_hi) {
        // (lo, hi)
        if (!(mn < I.hi) || !(mx > I.lo)) {
          // disjoint if mn >= hi OR mx <= lo
          if (mn >= I.hi || mx <= I.lo) return Rel::Disjoint;
        }
        // contained if mn > lo AND mx < hi
        if (!(mn > I.lo && mx < I.hi)) all_contained = false;
      } else if (I.has_hi) {
        // (-inf, hi)
        if (mn >= I.hi) return Rel::Disjoint;
        if (!(mx < I.hi)) all_contained = false;
      } else if (I.has_lo) {
        // (lo, +inf)
        if (mx <= I.lo) return Rel::Disjoint;
        if (!(mn > I.lo)) all_contained = false;
      } else {
        // (-inf, +inf)
        continue;
      }
    }

    return all_contained ? Rel::Contained : Rel::Partial;
  }

  void ComputeMBR(u32 l, u32 r, std::array<T, D>* mn, std::array<T, D>* mx) const {
    SJS_DASSERT(mn && mx);
    SJS_DASSERT(l < r);
    *mn = pts_[l].p.v;
    *mx = pts_[l].p.v;
    for (u32 i = l + 1; i < r; ++i) {
      for (int d = 0; d < D; ++d) {
        (*mn)[d] = std::min((*mn)[d], pts_[i].p.v[d]);
        (*mx)[d] = std::max((*mx)[d], pts_[i].p.v[d]);
      }
    }
  }

  int BuildNode(u32 l, u32 r, u32 depth) {
    const u32 n = r - l;
    SJS_DASSERT(n > 0);

    const int id = static_cast<int>(nodes_.size());
    nodes_.push_back(Node{});
    Node& node = nodes_.back();
    node.l = l;
    node.r = r;

    ComputeMBR(l, r, &node.mn, &node.mx);

    if (n <= leaf_size_) {
      node.leaf = true;
      node.left = -1;
      node.right = -1;
      return id;
    }

    const int split_dim = static_cast<int>(depth % static_cast<u32>(D));
    const u32 mid = l + (n / 2U);

    auto cmp = [split_dim](const Pt& a, const Pt& b) {
      const T ax = a.p.v[split_dim];
      const T bx = b.p.v[split_dim];
      if (ax < bx) return true;
      if (bx < ax) return false;
      return a.handle < b.handle;
    };

    std::nth_element(pts_.begin() + l, pts_.begin() + mid, pts_.begin() + r, cmp);

    node.leaf = false;
    node.left = BuildNode(l, mid, depth + 1U);
    node.right = BuildNode(mid, r, depth + 1U);
    return id;
  }

  u64 CountNode(int id, const QryT& q) const {
    if (id < 0) return 0ULL;
    const Node& n = nodes_[static_cast<usize>(id)];
    const Rel rel = RelateMBR(n, q);
    if (rel == Rel::Disjoint) return 0ULL;
    if (rel == Rel::Contained) return static_cast<u64>(n.r - n.l);

    if (n.leaf) {
      u64 cnt = 0;
      for (u32 i = n.l; i < n.r; ++i) {
        if (PointInQuery(pts_[i], q)) ++cnt;
      }
      return cnt;
    }

    return CountNode(n.left, q) + CountNode(n.right, q);
  }

  void ReportNode(int id, const QryT& q, std::vector<u32>* out) const {
    if (id < 0) return;
    const Node& n = nodes_[static_cast<usize>(id)];
    const Rel rel = RelateMBR(n, q);
    if (rel == Rel::Disjoint) return;

    if (rel == Rel::Contained) {
      for (u32 i = n.l; i < n.r; ++i) out->push_back(pts_[i].handle);
      return;
    }

    if (n.leaf) {
      for (u32 i = n.l; i < n.r; ++i) {
        if (PointInQuery(pts_[i], q)) out->push_back(pts_[i].handle);
      }
      return;
    }

    ReportNode(n.left, q, out);
    ReportNode(n.right, q, out);
  }

  // Internal helper: uses Block definition local to Sample().
  template <class BlockT>
  void CollectBlocks(int id, const QryT& q, std::vector<BlockT>* blocks, std::vector<u32>* hits) const {
    if (id < 0) return;
    const Node& n = nodes_[static_cast<usize>(id)];
    const Rel rel = RelateMBR(n, q);
    if (rel == Rel::Disjoint) return;

    if (rel == Rel::Contained) {
      BlockT b;
      b.contained = true;
      b.l = n.l;
      b.r = n.r;
      b.w = static_cast<u64>(n.r - n.l);
      blocks->push_back(b);
      return;
    }

    if (n.leaf) {
      const u32 off = static_cast<u32>(hits->size());
      for (u32 i = n.l; i < n.r; ++i) {
        if (PointInQuery(pts_[i], q)) hits->push_back(i);
      }
      const u32 len = static_cast<u32>(hits->size()) - off;
      if (len > 0) {
        BlockT b;
        b.contained = false;
        b.hit_off = off;
        b.hit_len = len;
        b.w = static_cast<u64>(len);
        blocks->push_back(b);
      }
      return;
    }

    CollectBlocks<BlockT>(n.left, q, blocks, hits);
    CollectBlocks<BlockT>(n.right, q, blocks, hits);
  }

  std::vector<Pt> pts_;
  std::vector<Node> nodes_;
  int root_ = -1;
  u32 leaf_size_ = 32U;
};

// Slot-plan helper: group sampled indices by key.
struct SlotPlan {
  std::vector<u32> offset;  // size E+1
  std::vector<u32> slots;   // length t
};

inline bool BuildSlotPlan(const std::vector<u64>& weights, u32 t, Rng* rng, SlotPlan* out) {
  if (!out || !rng) return false;
  const usize E = weights.size();
  if (E == 0) return false;

  out->offset.assign(E + 1, 0U);
  out->slots.clear();
  out->slots.resize(t);

  if (t == 0) return true;

  sampling::AliasTable alias;
  if (!alias.BuildFromU64(weights)) return false;

  std::vector<u32> idx_of_slot;
  idx_of_slot.resize(t);
  std::vector<u32> cnt(E, 0U);

  for (u32 j = 0; j < t; ++j) {
    const usize eid = alias.Sample(rng);
    if (eid >= E) return false;
    idx_of_slot[j] = static_cast<u32>(eid);
    cnt[eid]++;
  }

  for (usize i = 0; i < E; ++i) out->offset[i + 1] = out->offset[i] + cnt[i];

  std::vector<u32> cur = out->offset;
  for (u32 slot = 0; slot < t; ++slot) {
    const u32 eid = idx_of_slot[slot];
    const u32 pos = cur[eid]++;
    out->slots[pos] = slot;
  }

  return true;
}

// Streaming enumerator using kd-tree REPORT grouped by left box.
template <int Dim, class T>
class KdTreeJoinEnumerator : public IJoinEnumerator {
 public:
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;

  KdTreeJoinEnumerator(const DatasetT* ds, const KdTreeIndex<Dim, T>* idx)
      : ds_(ds), idx_(idx) {
    Reset();
  }

  bool Next(PairId* out) override {
    if (!out || !ds_ || !idx_) return false;

    while (r_pos_ < ds_->R.Size()) {
      if (match_pos_ < matches_.size()) {
        const u32 s_handle = matches_[match_pos_++];
        *out = PairId{ds_->R.GetId(r_pos_), ds_->S.GetId(s_handle)};
        stats_.output_pairs++;
        return true;
      }

      // Load matches for next r.
      matches_.clear();
      match_pos_ = 0;
      const BoxT& rbox = ds_->R.boxes[r_pos_];
      const auto q = MakeQueryForLeftBox<Dim, T>(rbox);
      idx_->Report(q, &matches_);
      r_pos_++;
    }

    return false;
  }

  void Reset() override {
    r_pos_ = 0;
    match_pos_ = 0;
    matches_.clear();
    stats_ = join::JoinStats{};
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  const DatasetT* ds_ = nullptr;
  const KdTreeIndex<Dim, T>* idx_ = nullptr;
  usize r_pos_ = 0;
  usize match_pos_ = 0;
  std::vector<u32> matches_;
  join::JoinStats stats_;
};

}  // namespace detail

// --------------------------
// Baseline: KdTreeSamplingBaseline
// --------------------------

template <int Dim, class T = Scalar>
class KdTreeSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 1, "KdTreeSamplingBaseline requires Dim >= 1");

  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;
  using BoxT = Box<Dim, T>;

  Method method() const noexcept override { return Method::KDTree; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "kd_tree_sampling"; }

  void Reset() override {
    built_ = false;
    ds_ = nullptr;
    idx_s_.Reset();
    w_r_.clear();
    W_ = 0;
    weights_valid_ = false;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    Reset();

    // Require proper (non-empty) half-open boxes.
    {
      std::string tmp;
      if (!ds.Validate(/*require_proper=*/true, &tmp)) {
        if (err) *err = "KdTreeSamplingBaseline::Build: invalid dataset: " + tmp;
        return false;
      }
    }

    ds_ = &ds;

    auto scoped = phases ? phases->Scoped("build_kd_tree") : PhaseRecorder::ScopedPhase(nullptr, "");

    if (!idx_s_.BuildFromRelation(ds.S)) {
      if (err) *err = "KdTreeSamplingBaseline::Build: failed to build kd-tree";
      return false;
    }

    w_r_.assign(ds.R.Size(), 0ULL);
    built_ = true;
    return true;
  }

  bool Count(const Config& cfg, Rng* rng, CountResult* out, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)rng;

    if (!built_ || !ds_) {
      if (err) *err = "KdTreeSamplingBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("count_join_via_kd_tree") : PhaseRecorder::ScopedPhase(nullptr, "");

    std::fill(w_r_.begin(), w_r_.end(), 0ULL);
    u64 W = 0;

    for (usize i = 0; i < ds_->R.Size(); ++i) {
      const BoxT& rbox = ds_->R.boxes[i];
      const auto q = detail::MakeQueryForLeftBox<Dim, T>(rbox);
      const u64 w = idx_s_.Count(q);
      w_r_[i] = w;
      if (w > 0) {
        if (W > (std::numeric_limits<u64>::max() - w)) {
          if (err) *err = "KdTreeSamplingBaseline::Count: |J| overflowed u64";
          return false;
        }
        W += w;
      }
    }

    W_ = W;
    weights_valid_ = true;

    if (out) *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "KdTreeSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "KdTreeSamplingBaseline::Sample: null rng/out";
      return false;
    }
    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "KdTreeSamplingBaseline::Sample: run.t too large (must fit in u32)";
      return false;
    }

    const u32 t = static_cast<u32>(cfg.run.t);
    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    if (t == 0) return true;

    // Ensure weights are ready.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/rng, &tmp, phases, err)) return false;
    }
    if (W_ == 0) return true;

    auto scoped = phases ? phases->Scoped("sample_join_via_kd_tree") : PhaseRecorder::ScopedPhase(nullptr, "");

    detail::SlotPlan plan;
    if (!detail::BuildSlotPlan(w_r_, t, rng, &plan)) {
      if (err) *err = "KdTreeSamplingBaseline::Sample: failed to build slot plan";
      return false;
    }

    out->pairs.resize(t);

    std::vector<u32> samp;
    samp.reserve(128);

    for (u32 r = 0; r < static_cast<u32>(ds_->R.Size()); ++r) {
      const u32 off = plan.offset[static_cast<usize>(r)];
      const u32 end = plan.offset[static_cast<usize>(r) + 1];
      const u32 k = end - off;
      if (k == 0) continue;

      // Sample k right partners for this r.
      const BoxT& rbox = ds_->R.boxes[static_cast<usize>(r)];
      const auto q = detail::MakeQueryForLeftBox<Dim, T>(rbox);

      samp.clear();
      if (!idx_s_.Sample(q, static_cast<u64>(k), rng, &samp) || samp.size() != k) {
        if (err) *err = "KdTreeSamplingBaseline::Sample: kd-tree SAMPLE failed";
        return false;
      }

      for (u32 j = 0; j < k; ++j) {
        const u32 slot = plan.slots[static_cast<usize>(off + j)];
        const u32 s_handle = samp[static_cast<usize>(j)];
        out->pairs[slot] = PairId{ds_->R.GetId(static_cast<usize>(r)), ds_->S.GetId(static_cast<usize>(s_handle))};
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || !ds_) {
      if (err) *err = "KdTreeSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::KdTreeJoinEnumerator<Dim, T>>(ds_, &idx_s_);
  }

 private:
  bool built_ = false;
  const DatasetT* ds_ = nullptr;

  detail::KdTreeIndex<Dim, T> idx_s_;

  std::vector<u64> w_r_;
  u64 W_ = 0;
  bool weights_valid_ = false;
};

}  // namespace kd_tree
}  // namespace baselines
}  // namespace sjs
