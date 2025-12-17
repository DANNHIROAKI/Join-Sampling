#pragma once
// sjs/baselines/interval_tree/sampling.h
//
// Plane Sweep + IntervalTreeSampler baseline (Variant::Sampling).
//
// This implements the algorithm in "IntervalTree Baseline.md":
//   - Sweep on axis 0 (x).
//   - Maintain two dynamic y-interval samplers (one for R, one for S).
//   - At each START event, x-overlap is guaranteed by sweep semantics, so the
//     join condition reduces to 1D overlap on y.
//   - Use the A/B decomposition for y-interval overlap:
//        Overlap([a,b)) = A([a,b)) \uplus B([a,b))
//        A:  l <= a < r     (stabbing at point a)
//        B:  a < l < b      (left endpoint in (a,b))
//   - Two-pass sampling:
//        Phase1: count wA,wB per START and W=|J|
//        Phase2: alias on w and assign sample slots to (event, A/B)
//        Phase3: resweep and locally sample i.i.d. from A or B to fill slots
//
// Notes on coordinate compression
// ------------------------------
// In 2D, the query point for A is always a rectangle's y-lower value (a = Ly(q)).
// To reduce memory overhead, we build the y-domain from *unique y-lower values*
// across both relations. This is sufficient for correctness of A/B in our setting.
// If you later want to support arbitrary query points, extend the domain to all
// y-endpoints (both Ly and Ry).

#include "sjs/baselines/baseline_api.h"
#include "sjs/core/assert.h"
#include "sjs/join/sweep_events.h"
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

namespace detail {

// --------------------------
// Dense bucket helper
// --------------------------
// A bucket stores handles densely and supports:
//   - O(1) uniform sampling
//   - O(1) delete by (node, position)
// We achieve O(1) delete by storing, for each inserted element, a "backref"
// index into a per-handle placement list.

struct BucketItem {
  u32 handle = 0;
  u32 backref = 0;  // index within this handle's placement list
};

struct Placement {
  u32 node = 0;     // segment tree node id
  u32 pos = 0;      // position inside node bucket
  u32 backref = 0;  // redundant but handy for debugging
};

// --------------------------
// Segment tree for A-type queries (stabbing at a point)
// --------------------------
// Domain: discrete points [0..m). Query is a point q.
// Insert interval [lo, hi) (in rank space) by canonical cover.
// Then for a query point q, all intervals that contain q appear exactly once
// across buckets on the root-to-leaf path of q.
class StabbingSegTree {
 public:
  StabbingSegTree() = default;

  void Clear() {
    n_handles_ = 0;
    m_ = 0;
    p_ = 0;
    max_refs_ = 0;
    nodes_.clear();
    placements_.clear();
    placement_size_.clear();
    touched_nodes_.clear();
    node_touched_.clear();
  }

  void Init(u32 n_handles, u32 m_points) {
    Clear();
    n_handles_ = n_handles;
    m_ = m_points;
    if (m_ == 0) return;

    // p_ = smallest power of two >= m_
    p_ = 1;
    while (p_ < m_) p_ <<= 1;

    // max refs per interval for canonical cover is <= 2*log2(p_) + 2.
    u32 logp = 0;
    for (u32 x = p_; x > 1; x >>= 1) ++logp;
    max_refs_ = 2 * logp + 4;

    nodes_.resize(static_cast<usize>(2 * p_));

    placements_.resize(static_cast<usize>(n_handles_) * static_cast<usize>(max_refs_));
    placement_size_.assign(n_handles_, 0);

    touched_nodes_.clear();
    node_touched_.assign(nodes_.size(), 0);
  }

  bool Empty() const noexcept { return m_ == 0; }

  void ResetToEmpty() {
    // Clear only nodes that were ever non-empty.
    for (u32 nid : touched_nodes_) {
      nodes_[static_cast<usize>(nid)].items.clear();
      node_touched_[static_cast<usize>(nid)] = 0;
    }
    touched_nodes_.clear();
    // Reset placements.
    std::fill(placement_size_.begin(), placement_size_.end(), 0);
  }

  // Insert interval [lo, hi) in rank space.
  // Precondition: 0 <= lo < hi <= m_. (hi may equal m_)
  void Insert(u32 handle, u32 lo, u32 hi) {
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(m_ > 0);
    SJS_DASSERT(lo < hi);
    SJS_DASSERT(hi <= m_);

    const u32 base = handle * max_refs_;
    u32& psz = placement_size_[handle];
    SJS_DASSERT(psz == 0);

    u32 l = lo + p_;
    u32 r = hi + p_;
    while (l < r) {
      if (l & 1U) {
        AddToNode(handle, l, base, psz);
        ++l;
      }
      if (r & 1U) {
        --r;
        AddToNode(handle, r, base, psz);
      }
      l >>= 1;
      r >>= 1;
    }

    // Defensive: should not exceed max_refs_.
    SJS_DASSERT(psz <= max_refs_);
  }

  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_handles_);
    const u32 base = handle * max_refs_;
    const u32 psz = placement_size_[handle];
    // It is valid to call Erase on a non-active handle only if psz==0.
    if (psz == 0) return;

    // Remove from each referenced node (swap-delete inside bucket).
    for (u32 i = 0; i < psz; ++i) {
      const Placement pl = placements_[static_cast<usize>(base + i)];
      RemoveFromNode(pl.node, pl.pos);
    }
    placement_size_[handle] = 0;
  }

  // Count intervals containing query point q.
  u64 Count(u32 q) const {
    SJS_DASSERT(q < m_);
    u64 cnt = 0;
    u32 idx = q + p_;
    while (idx > 0) {
      cnt += static_cast<u64>(nodes_[static_cast<usize>(idx)].items.size());
      idx >>= 1;
    }
    return cnt;
  }

  // Report intervals containing query point q.
  void Report(u32 q, std::vector<u32>* out) const {
    if (!out) return;
    SJS_DASSERT(q < m_);
    u32 idx = q + p_;
    while (idx > 0) {
      const auto& vec = nodes_[static_cast<usize>(idx)].items;
      for (const auto& it : vec) out->push_back(it.handle);
      idx >>= 1;
    }
  }

  // Sample k i.i.d. handles uniformly from the stabbing set at q.
  // Returns false if the stabbing set is empty.
  bool Sample(u32 q, u32 k, Rng* rng, std::vector<u32>* out) const {
    if (!rng || !out) return false;
    out->clear();
    if (k == 0) return true;
    if (m_ == 0) return false;
    SJS_DASSERT(q < m_);

    // Collect non-empty buckets along the path.
    u32 idx = q + p_;
    std::vector<u32> nodes;
    nodes.reserve(64);
    std::vector<u64> weights;
    weights.reserve(64);

    u64 total = 0;
    while (idx > 0) {
      const u64 w = static_cast<u64>(nodes_[static_cast<usize>(idx)].items.size());
      if (w > 0) {
        nodes.push_back(idx);
        weights.push_back(w);
        total += w;
      }
      idx >>= 1;
    }
    if (total == 0) return false;

    sampling::AliasTable alias;
    std::string tmp_err;
    if (!alias.BuildFromU64(Span<const u64>(weights), &tmp_err)) {
      // Should not fail for non-negative weights.
      return false;
    }

    out->reserve(static_cast<usize>(k));
    for (u32 i = 0; i < k; ++i) {
      const usize bi = alias.Sample(rng);
      const u32 node_id = nodes[bi];
      const auto& bucket = nodes_[static_cast<usize>(node_id)].items;
      SJS_DASSERT(!bucket.empty());
      const u32 pos = static_cast<u32>(rng->UniformU64(static_cast<u64>(bucket.size())));
      out->push_back(bucket[static_cast<usize>(pos)].handle);
    }
    return true;
  }

 private:
  struct Node {
    std::vector<BucketItem> items;
  };

  void TouchNode(u32 node) {
    const usize i = static_cast<usize>(node);
    if (!node_touched_[i]) {
      node_touched_[i] = 1;
      touched_nodes_.push_back(node);
    }
  }

  inline usize PlacementIndex(u32 handle, u32 backref) const noexcept {
    return static_cast<usize>(handle) * static_cast<usize>(max_refs_) + static_cast<usize>(backref);
  }

  void AddToNode(u32 handle, u32 node, u32 base, u32& psz) {
    SJS_DASSERT(psz < max_refs_);
    TouchNode(node);

    Node& nd = nodes_[static_cast<usize>(node)];
    const u32 pos = static_cast<u32>(nd.items.size());
    const u32 backref = psz;
    nd.items.push_back(BucketItem{handle, backref});

    placements_[static_cast<usize>(base + psz)] = Placement{node, pos, backref};
    ++psz;
  }

  void RemoveFromNode(u32 node, u32 pos) {
    Node& nd = nodes_[static_cast<usize>(node)];
    auto& vec = nd.items;
    SJS_DASSERT(pos < vec.size());

    const BucketItem moved = vec.back();
    vec[pos] = moved;
    vec.pop_back();

    // Update moved element's recorded position.
    const usize pidx = PlacementIndex(moved.handle, moved.backref);
    placements_[pidx].pos = pos;
  }

  u32 n_handles_{0};
  u32 m_{0};
  u32 p_{0};
  u32 max_refs_{0};

  std::vector<Node> nodes_;
  std::vector<Placement> placements_;
  std::vector<u32> placement_size_;

  // Sparse reset support.
  std::vector<u32> touched_nodes_;
  std::vector<u8> node_touched_;
};

// --------------------------
// Segment tree for B-type queries (start-in-range)
// --------------------------
// Domain: discrete ranks of left endpoints in [0..m).
// Each active interval contributes a point at rank = rank(Ly).
// We insert each point into all nodes on its leaf-to-root path.
// Query a rank range [l,r) by canonical decomposition; buckets are disjoint.
class StartRangeSegTree {
 public:
  StartRangeSegTree() = default;

  void Clear() {
    n_handles_ = 0;
    m_ = 0;
    p_ = 0;
    max_refs_ = 0;
    nodes_.clear();
    placements_.clear();
    placement_size_.clear();
    touched_nodes_.clear();
    node_touched_.clear();
  }

  void Init(u32 n_handles, u32 m_points) {
    Clear();
    n_handles_ = n_handles;
    m_ = m_points;
    if (m_ == 0) return;

    p_ = 1;
    while (p_ < m_) p_ <<= 1;

    // Each point is stored along a leaf-to-root path: log2(p_)+1 nodes.
    u32 logp = 0;
    for (u32 x = p_; x > 1; x >>= 1) ++logp;
    max_refs_ = logp + 2;

    nodes_.resize(static_cast<usize>(2 * p_));
    placements_.resize(static_cast<usize>(n_handles_) * static_cast<usize>(max_refs_));
    placement_size_.assign(n_handles_, 0);

    touched_nodes_.clear();
    node_touched_.assign(nodes_.size(), 0);
  }

  void ResetToEmpty() {
    for (u32 nid : touched_nodes_) {
      nodes_[static_cast<usize>(nid)].items.clear();
      node_touched_[static_cast<usize>(nid)] = 0;
    }
    touched_nodes_.clear();
    std::fill(placement_size_.begin(), placement_size_.end(), 0);
  }

  // Insert a point (handle) at rank.
  void Insert(u32 handle, u32 rank) {
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(rank < m_);

    const u32 base = handle * max_refs_;
    u32& psz = placement_size_[handle];
    SJS_DASSERT(psz == 0);

    u32 node = rank + p_;
    while (node > 0) {
      AddToNode(handle, node, base, psz);
      node >>= 1;
    }
    SJS_DASSERT(psz <= max_refs_);
  }

  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_handles_);
    const u32 base = handle * max_refs_;
    const u32 psz = placement_size_[handle];
    if (psz == 0) return;
    for (u32 i = 0; i < psz; ++i) {
      const Placement pl = placements_[static_cast<usize>(base + i)];
      RemoveFromNode(pl.node, pl.pos);
    }
    placement_size_[handle] = 0;
  }

  // Count points with rank in [l,r).
  u64 CountRange(u32 l, u32 r) const {
    if (l >= r) return 0;
    SJS_DASSERT(r <= m_);

    u64 cnt = 0;
    u32 L = l + p_;
    u32 R = r + p_;
    while (L < R) {
      if (L & 1U) {
        cnt += static_cast<u64>(nodes_[static_cast<usize>(L)].items.size());
        ++L;
      }
      if (R & 1U) {
        --R;
        cnt += static_cast<u64>(nodes_[static_cast<usize>(R)].items.size());
      }
      L >>= 1;
      R >>= 1;
    }
    return cnt;
  }

  // Report points with rank in [l,r).
  void ReportRange(u32 l, u32 r, std::vector<u32>* out) const {
    if (!out) return;
    if (l >= r) return;
    SJS_DASSERT(r <= m_);

    u32 L = l + p_;
    u32 R = r + p_;
    while (L < R) {
      if (L & 1U) {
        const auto& vec = nodes_[static_cast<usize>(L)].items;
        for (const auto& it : vec) out->push_back(it.handle);
        ++L;
      }
      if (R & 1U) {
        --R;
        const auto& vec = nodes_[static_cast<usize>(R)].items;
        for (const auto& it : vec) out->push_back(it.handle);
      }
      L >>= 1;
      R >>= 1;
    }
  }

  // Sample k i.i.d. points with rank in [l,r).
  // Returns false if the set is empty.
  bool SampleRange(u32 l, u32 r, u32 k, Rng* rng, std::vector<u32>* out) const {
    if (!rng || !out) return false;
    out->clear();
    if (k == 0) return true;
    if (l >= r) return false;
    SJS_DASSERT(r <= m_);

    // Collect non-empty canonical buckets for [l,r).
    std::vector<u32> nodes;
    nodes.reserve(64);
    std::vector<u64> weights;
    weights.reserve(64);

    u32 L = l + p_;
    u32 R = r + p_;
    while (L < R) {
      if (L & 1U) {
        const u64 w = static_cast<u64>(nodes_[static_cast<usize>(L)].items.size());
        if (w > 0) {
          nodes.push_back(L);
          weights.push_back(w);
        }
        ++L;
      }
      if (R & 1U) {
        --R;
        const u64 w = static_cast<u64>(nodes_[static_cast<usize>(R)].items.size());
        if (w > 0) {
          nodes.push_back(R);
          weights.push_back(w);
        }
      }
      L >>= 1;
      R >>= 1;
    }

    if (nodes.empty()) return false;

    sampling::AliasTable alias;
    std::string tmp_err;
    if (!alias.BuildFromU64(Span<const u64>(weights), &tmp_err)) {
      return false;
    }

    out->reserve(static_cast<usize>(k));
    for (u32 i = 0; i < k; ++i) {
      const usize bi = alias.Sample(rng);
      const u32 node_id = nodes[bi];
      const auto& bucket = nodes_[static_cast<usize>(node_id)].items;
      SJS_DASSERT(!bucket.empty());
      const u32 pos = static_cast<u32>(rng->UniformU64(static_cast<u64>(bucket.size())));
      out->push_back(bucket[static_cast<usize>(pos)].handle);
    }
    return true;
  }

 private:
  struct Node {
    std::vector<BucketItem> items;
  };

  void TouchNode(u32 node) {
    const usize i = static_cast<usize>(node);
    if (!node_touched_[i]) {
      node_touched_[i] = 1;
      touched_nodes_.push_back(node);
    }
  }

  inline usize PlacementIndex(u32 handle, u32 backref) const noexcept {
    return static_cast<usize>(handle) * static_cast<usize>(max_refs_) + static_cast<usize>(backref);
  }

  void AddToNode(u32 handle, u32 node, u32 base, u32& psz) {
    SJS_DASSERT(psz < max_refs_);
    TouchNode(node);

    Node& nd = nodes_[static_cast<usize>(node)];
    const u32 pos = static_cast<u32>(nd.items.size());
    const u32 backref = psz;
    nd.items.push_back(BucketItem{handle, backref});

    placements_[static_cast<usize>(base + psz)] = Placement{node, pos, backref};
    ++psz;
  }

  void RemoveFromNode(u32 node, u32 pos) {
    Node& nd = nodes_[static_cast<usize>(node)];
    auto& vec = nd.items;
    SJS_DASSERT(pos < vec.size());

    const BucketItem moved = vec.back();
    vec[pos] = moved;
    vec.pop_back();

    const usize pidx = PlacementIndex(moved.handle, moved.backref);
    placements_[pidx].pos = pos;
  }

  u32 n_handles_{0};
  u32 m_{0};
  u32 p_{0};
  u32 max_refs_{0};

  std::vector<Node> nodes_;
  std::vector<Placement> placements_;
  std::vector<u32> placement_size_;

  std::vector<u32> touched_nodes_;
  std::vector<u8> node_touched_;
};

// --------------------------
// IntervalTreeSampler2D: IT-A + IT-B
// --------------------------
class IntervalTreeSampler2D {
 public:
  void Clear() {
    n_handles_ = 0;
    m_ = 0;
    active_size_ = 0;
    a_.Clear();
    b_.Clear();
  }

  void Init(u32 n_handles, u32 m_domain) {
    Clear();
    n_handles_ = n_handles;
    m_ = m_domain;
    active_size_ = 0;
    a_.Init(n_handles_, m_);
    b_.Init(n_handles_, m_);
  }

  void ResetToEmpty() {
    active_size_ = 0;
    a_.ResetToEmpty();
    b_.ResetToEmpty();
  }

  u32 ActiveSize() const noexcept { return active_size_; }

  void Insert(u32 handle, u32 ylo_rank, u32 yhi_lb_rank) {
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(ylo_rank < m_);
    SJS_DASSERT(yhi_lb_rank <= m_);
    SJS_DASSERT(ylo_rank < yhi_lb_rank);

    a_.Insert(handle, ylo_rank, yhi_lb_rank);
    b_.Insert(handle, ylo_rank);
    ++active_size_;
  }

  void Erase(u32 handle) {
    // It is ok to call Erase on a handle that is not active (no-op), but
    // in our sweep it should be active exactly once.
    a_.Erase(handle);
    b_.Erase(handle);
    if (active_size_ > 0) --active_size_;
  }

  // A-type count: intervals with l <= a < r.
  u64 CountA(u32 a_rank) const { return a_.Count(a_rank); }

  // B-type count: intervals with a < l < b.
  u64 CountB(u32 a_rank, u32 b_lb_rank) const {
    // open interval (a,b) -> ranks (a_rank, b_lb_rank)
    const u32 l = a_rank + 1;
    const u32 r = b_lb_rank;
    if (l >= r) return 0;
    return b_.CountRange(l, r);
  }

  void ReportA(u32 a_rank, std::vector<u32>* out) const { a_.Report(a_rank, out); }

  void ReportB(u32 a_rank, u32 b_lb_rank, std::vector<u32>* out) const {
    const u32 l = a_rank + 1;
    const u32 r = b_lb_rank;
    if (l >= r) return;
    b_.ReportRange(l, r, out);
  }

  bool SampleA(u32 a_rank, u32 k, Rng* rng, std::vector<u32>* out) const {
    return a_.Sample(a_rank, k, rng, out);
  }

  bool SampleB(u32 a_rank, u32 b_lb_rank, u32 k, Rng* rng, std::vector<u32>* out) const {
    const u32 l = a_rank + 1;
    const u32 r = b_lb_rank;
    if (l >= r) {
      if (out) out->clear();
      return false;
    }
    return b_.SampleRange(l, r, k, rng, out);
  }

 private:
  u32 n_handles_{0};
  u32 m_{0};
  u32 active_size_{0};

  StabbingSegTree a_;
  StartRangeSegTree b_;
};

// --------------------------
// Deterministic join enumerator: sweep + IntervalTreeSampler2D
// --------------------------
// Enumerates every intersecting (R,S) pair exactly once by the same decomposition:
// each pair is emitted at the later START event under the sweep ordering.

template <int Dim, class T>
class IntervalTreeJoinEnumerator final : public IJoinEnumerator {
 public:
  using RelationT = Relation<Dim, T>;

  IntervalTreeJoinEnumerator(const RelationT* R,
                            const RelationT* S,
                            int axis,
                            join::SideTieBreak side_order)
      : R_(R), S_(S), axis_(axis), side_order_(side_order) {
    SJS_DASSERT(R_ != nullptr && S_ != nullptr);
    BuildPreproc();
    Reset();
  }

  void Reset() override {
    // Reset streaming state.
    ev_pos_ = 0;
    in_emit_ = false;
    cand_.clear();
    cand_pos_ = 0;
    cur_side_ = join::Side::R;
    cur_index_ = 0;

    // Reset sweep state.
    active_r_.ResetToEmpty();
    active_s_.ResetToEmpty();
    active_r_size_ = 0;
    active_s_size_ = 0;

    stats_.Reset();
    stats_.num_events = static_cast<u64>(events_.size());
  }

  bool Next(PairId* out) override {
    if (!out) return false;

    // If we finished emitting candidates for a START event, we still need to
    // insert that START rectangle into its active structure before continuing.
    if (in_emit_ && cand_pos_ >= cand_.size()) {
      FinishStartEvent();
    }

    // If we are mid-emission, output next candidate.
    if (in_emit_) {
      SJS_DASSERT(cand_pos_ < cand_.size());
      EmitOne(out);
      return true;
    }

    // Otherwise advance events until we find a START with at least one match,
    // or we reach end.
    while (ev_pos_ < events_.size()) {
      const join::Event& e = events_[ev_pos_++];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) {
          active_r_.Erase(static_cast<u32>(e.index));
          if (active_r_size_ > 0) --active_r_size_;
        } else {
          active_s_.Erase(static_cast<u32>(e.index));
          if (active_s_size_ > 0) --active_s_size_;
        }
        stats_.active_max_r = std::max<u64>(stats_.active_max_r, active_r_size_);
        stats_.active_max_s = std::max<u64>(stats_.active_max_s, active_s_size_);
        continue;
      }

      // START event.
      cur_side_ = e.side;
      cur_index_ = e.index;
      cand_.clear();
      cand_pos_ = 0;

      if (cur_side_ == join::Side::R) {
        const u32 a = ylo_rank_r_[cur_index_];
        const u32 b = yhi_lb_rank_r_[cur_index_];
        active_s_.ReportA(a, &cand_);
        active_s_.ReportB(a, b, &cand_);
      } else {
        const u32 a = ylo_rank_s_[cur_index_];
        const u32 b = yhi_lb_rank_s_[cur_index_];
        active_r_.ReportA(a, &cand_);
        active_r_.ReportB(a, b, &cand_);
      }

      if (cand_.empty()) {
        // No output for this START; still need to insert it and continue.
        if (cur_side_ == join::Side::R) {
          active_r_.Insert(static_cast<u32>(cur_index_), ylo_rank_r_[cur_index_], yhi_lb_rank_r_[cur_index_]);
          ++active_r_size_;
        } else {
          active_s_.Insert(static_cast<u32>(cur_index_), ylo_rank_s_[cur_index_], yhi_lb_rank_s_[cur_index_]);
          ++active_s_size_;
        }
        stats_.active_max_r = std::max<u64>(stats_.active_max_r, active_r_size_);
        stats_.active_max_s = std::max<u64>(stats_.active_max_s, active_s_size_);
        continue;
      }

      // Begin emitting pairs for this START.
      in_emit_ = true;
      EmitOne(out);
      return true;
    }

    // End of stream.
    return false;
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  void BuildPreproc() {
    // Build events.
    events_ = join::BuildSweepEvents<Dim, T>(*R_, *S_, axis_, side_order_);

    // Build y-domain from unique y-lower values (axis 1).
    std::vector<T> y;
    y.reserve(R_->Size() + S_->Size());
    for (const auto& b : R_->boxes) y.push_back(b.lo.v[1]);
    for (const auto& b : S_->boxes) y.push_back(b.lo.v[1]);
    std::sort(y.begin(), y.end());
    y.erase(std::unique(y.begin(), y.end()), y.end());
    y_coords_ = std::move(y);

    const u32 m = static_cast<u32>(y_coords_.size());
    if (m == 0) {
      // Degenerate; keep empty structures.
      active_r_.Init(static_cast<u32>(R_->Size()), 0);
      active_s_.Init(static_cast<u32>(S_->Size()), 0);
      return;
    }

    auto lb = [&](T v) -> u32 {
      const auto it = std::lower_bound(y_coords_.begin(), y_coords_.end(), v);
      return static_cast<u32>(std::distance(y_coords_.begin(), it));
    };

    ylo_rank_r_.resize(R_->Size());
    yhi_lb_rank_r_.resize(R_->Size());
    for (usize i = 0; i < R_->Size(); ++i) {
      const auto& b = R_->boxes[i];
      const u32 lo = lb(b.lo.v[1]);
      ylo_rank_r_[i] = lo;
      yhi_lb_rank_r_[i] = lb(b.hi.v[1]);
    }

    ylo_rank_s_.resize(S_->Size());
    yhi_lb_rank_s_.resize(S_->Size());
    for (usize i = 0; i < S_->Size(); ++i) {
      const auto& b = S_->boxes[i];
      const u32 lo = lb(b.lo.v[1]);
      ylo_rank_s_[i] = lo;
      yhi_lb_rank_s_[i] = lb(b.hi.v[1]);
    }

    active_r_.Init(static_cast<u32>(R_->Size()), m);
    active_s_.Init(static_cast<u32>(S_->Size()), m);
  }

  void EmitOne(PairId* out) {
    const u32 other_h = cand_[cand_pos_++];

    if (cur_side_ == join::Side::R) {
      *out = PairId{R_->GetId(cur_index_), S_->GetId(static_cast<usize>(other_h))};
    } else {
      *out = PairId{R_->GetId(static_cast<usize>(other_h)), S_->GetId(cur_index_)};
    }

    ++stats_.output_pairs;
    ++stats_.candidate_checks;

    // If this was the last candidate, we keep in_emit_=true; the next call will
    // finish the event by inserting cur into the active structure.
    if (cand_pos_ >= cand_.size()) {
      // no-op here
    }
  }

  void FinishStartEvent() {
    // Insert the START rectangle into its side's active structure.
    if (cur_side_ == join::Side::R) {
      active_r_.Insert(static_cast<u32>(cur_index_), ylo_rank_r_[cur_index_], yhi_lb_rank_r_[cur_index_]);
      ++active_r_size_;
    } else {
      active_s_.Insert(static_cast<u32>(cur_index_), ylo_rank_s_[cur_index_], yhi_lb_rank_s_[cur_index_]);
      ++active_s_size_;
    }
    stats_.active_max_r = std::max<u64>(stats_.active_max_r, active_r_size_);
    stats_.active_max_s = std::max<u64>(stats_.active_max_s, active_s_size_);

    in_emit_ = false;
    cand_.clear();
    cand_pos_ = 0;
  }

  const RelationT* R_;
  const RelationT* S_;
  int axis_{0};
  join::SideTieBreak side_order_{join::SideTieBreak::RBeforeS};

  // Preprocessed sweep.
  std::vector<join::Event> events_;

  std::vector<T> y_coords_;
  std::vector<u32> ylo_rank_r_;
  std::vector<u32> yhi_lb_rank_r_;
  std::vector<u32> ylo_rank_s_;
  std::vector<u32> yhi_lb_rank_s_;

  IntervalTreeSampler2D active_r_;
  IntervalTreeSampler2D active_s_;

  // Streaming state.
  usize ev_pos_{0};
  bool in_emit_{false};
  join::Side cur_side_{join::Side::R};
  usize cur_index_{0};
  std::vector<u32> cand_;
  usize cand_pos_{0};

  // Simple active sizes for stats.
  u64 active_r_size_{0};
  u64 active_s_size_{0};

  join::JoinStats stats_;
};

}  // namespace detail

// --------------------------
// Baseline: IntervalTree + Sampling
// --------------------------

template <int Dim, class T = Scalar>
class IntervalTreeSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::IntervalTree; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "interval_tree_sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;

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

    W_ = 0;
    weights_valid_ = false;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;

    if constexpr (Dim != 2) {
      if (err) *err = "IntervalTreeSamplingBaseline: currently only Dim=2 is supported";
      return false;
    }

    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;

    if (ds.R.Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        ds.S.Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "IntervalTreeSamplingBaseline::Build: relation size exceeds u32 limit";
      return false;
    }

    // Events.
    {
      auto ph = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    // START id mapping.
    {
      auto ph = phases ? phases->Scoped("build_start_ids") : PhaseRecorder::ScopedPhase(nullptr, "");
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
    }

    // y-domain: unique y-lower values.
    {
      auto ph = phases ? phases->Scoped("build_y_domain") : PhaseRecorder::ScopedPhase(nullptr, "");
      y_coords_.clear();
      y_coords_.reserve(ds.R.Size() + ds.S.Size());
      for (const auto& b : ds.R.boxes) y_coords_.push_back(b.lo.v[1]);
      for (const auto& b : ds.S.boxes) y_coords_.push_back(b.lo.v[1]);
      std::sort(y_coords_.begin(), y_coords_.end());
      y_coords_.erase(std::unique(y_coords_.begin(), y_coords_.end()), y_coords_.end());
      if (y_coords_.empty()) {
        if (err) *err = "IntervalTreeSamplingBaseline::Build: empty y domain";
        return false;
      }
    }

    // Precompute ranks.
    {
      auto ph = phases ? phases->Scoped("build_y_ranks") : PhaseRecorder::ScopedPhase(nullptr, "");
      auto lb = [&](T v) -> u32 {
        const auto it = std::lower_bound(y_coords_.begin(), y_coords_.end(), v);
        return static_cast<u32>(std::distance(y_coords_.begin(), it));
      };

      ylo_rank_r_.resize(ds.R.Size());
      yhi_lb_rank_r_.resize(ds.R.Size());
      for (usize i = 0; i < ds.R.Size(); ++i) {
        const auto& b = ds.R.boxes[i];
        const u32 lo = lb(b.lo.v[1]);
        SJS_DASSERT(lo < y_coords_.size() && y_coords_[lo] == b.lo.v[1]);
        ylo_rank_r_[i] = lo;
        yhi_lb_rank_r_[i] = lb(b.hi.v[1]);
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

    // Init active samplers.
    {
      auto ph = phases ? phases->Scoped("build_active") : PhaseRecorder::ScopedPhase(nullptr, "");
      const u32 m = static_cast<u32>(y_coords_.size());
      active_r_.Init(static_cast<u32>(ds.R.Size()), m);
      active_s_.Init(static_cast<u32>(ds.S.Size()), m);
    }

    // Allocate weights per START event.
    {
      auto ph = phases ? phases->Scoped("build_weights") : PhaseRecorder::ScopedPhase(nullptr, "");
      const usize E = start_event_pos_.size();
      w_total_.assign(E, 0ULL);
      w_a_.assign(E, 0ULL);
      w_b_.assign(E, 0ULL);
    }

    built_ = true;
    weights_valid_ = false;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;  // deterministic

    if (!built_ || !ds_) {
      if (err) *err = "IntervalTreeSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "IntervalTreeSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    std::fill(w_a_.begin(), w_a_.end(), 0ULL);
    std::fill(w_b_.begin(), w_b_.end(), 0ULL);

    // Ensure active sets are empty.
    active_r_.ResetToEmpty();
    active_s_.ResetToEmpty();

    u64 W = 0;

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

      // START event.
      const i32 sid_i32 = start_id_of_event_[ev_pos];
      SJS_DASSERT(sid_i32 >= 0);
      const usize sid = static_cast<usize>(sid_i32);

      if (e.side == join::Side::R) {
        const u32 a = ylo_rank_r_[e.index];
        const u32 b = yhi_lb_rank_r_[e.index];

        const u64 wa = active_s_.CountA(a);
        const u64 wb = active_s_.CountB(a, b);
        const u64 w = wa + wb;

        w_a_[sid] = wa;
        w_b_[sid] = wb;
        w_total_[sid] = w;

        if (W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "IntervalTreeSamplingBaseline::Count: |J| overflow (u64)";
          return false;
        }
        W += w;

        active_r_.Insert(static_cast<u32>(e.index), a, b);
      } else {
        const u32 a = ylo_rank_s_[e.index];
        const u32 b = yhi_lb_rank_s_[e.index];

        const u64 wa = active_r_.CountA(a);
        const u64 wb = active_r_.CountB(a, b);
        const u64 w = wa + wb;

        w_a_[sid] = wa;
        w_b_[sid] = wb;
        w_total_[sid] = w;

        if (W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "IntervalTreeSamplingBaseline::Count: |J| overflow (u64)";
          return false;
        }
        W += w;

        active_s_.Insert(static_cast<u32>(e.index), a, b);
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
      if (err) *err = "IntervalTreeSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "IntervalTreeSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "IntervalTreeSamplingBaseline::Sample: out is null";
      return false;
    }

    const u32 t = cfg.run.t;

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    if (t == 0) return true;

    // Ensure weights are available.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    if (W_ == 0) {
      // Empty join.
      return true;
    }

    // --------------------------
    // Phase2: assign sample slots to (event, A/B)
    // --------------------------
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
        if (err && err->empty()) *err = "IntervalTreeSamplingBaseline::Sample: failed to build alias table";
        return false;
      }

      assign.reserve(static_cast<usize>(t));
      for (u32 j = 0; j < t; ++j) {
        // Robustly avoid any accidental zero-weight event (should not happen).
        u32 sid = 0;
        u64 w = 0;
        for (int tries = 0; tries < 16; ++tries) {
          sid = static_cast<u32>(alias.Sample(rng));
          w = w_total_[sid];
          if (w > 0) break;
        }
        if (w == 0) {
          if (err) *err = "IntervalTreeSamplingBaseline::Sample: alias produced only zero-weight events";
          return false;
        }

        const u64 wa = w_a_[sid];
        const u64 wb = w_b_[sid];
        SJS_DASSERT(wa + wb == w);

        u8 pat = 0;
        if (wa == 0) {
          pat = 1;
        } else if (wb == 0) {
          pat = 0;
        } else {
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

    out->pairs.resize(static_cast<usize>(t));

    // --------------------------
    // Phase3: second sweep; for each START event, fulfill its assigned slots
    // --------------------------
    {
      auto scoped = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

      active_r_.ResetToEmpty();
      active_s_.ResetToEmpty();

      usize ptr = 0;
      std::vector<u32> tmp;
      tmp.reserve(1024);

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
          while (ptr < assign.size() && assign[ptr].sid == sid && assign[ptr].pat == pat) {
            ++ptr;
          }
          const u32 k = static_cast<u32>(ptr - begin);
          if (k == 0) continue;

          bool ok = false;

          if (e.side == join::Side::R) {
            const u32 a = ylo_rank_r_[e.index];
            const u32 b = yhi_lb_rank_r_[e.index];

            if (pat == 0) ok = active_s_.SampleA(a, k, rng, &tmp);
            else ok = active_s_.SampleB(a, b, k, rng, &tmp);

            if (!ok) {
              if (err) *err = "IntervalTreeSamplingBaseline::Sample: empty candidate set in phase3 (R-start)";
              return false;
            }

            for (u32 i = 0; i < k; ++i) {
              const u32 slot = assign[begin + i].slot;
              const u32 other_h = tmp[static_cast<usize>(i)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(e.index), ds_->S.GetId(other_h)};
            }

          } else {
            const u32 a = ylo_rank_s_[e.index];
            const u32 b = yhi_lb_rank_s_[e.index];

            if (pat == 0) ok = active_r_.SampleA(a, k, rng, &tmp);
            else ok = active_r_.SampleB(a, b, k, rng, &tmp);

            if (!ok) {
              if (err) *err = "IntervalTreeSamplingBaseline::Sample: empty candidate set in phase3 (S-start)";
              return false;
            }

            for (u32 i = 0; i < k; ++i) {
              const u32 slot = assign[begin + i].slot;
              const u32 other_h = tmp[static_cast<usize>(i)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(other_h), ds_->S.GetId(e.index)};
            }
          }
        }

        // Insert the START rectangle into its side.
        if (e.side == join::Side::R) {
          const u32 a = ylo_rank_r_[e.index];
          const u32 b = yhi_lb_rank_r_[e.index];
          active_r_.Insert(static_cast<u32>(e.index), a, b);
        } else {
          const u32 a = ylo_rank_s_[e.index];
          const u32 b = yhi_lb_rank_s_[e.index];
          active_s_.Insert(static_cast<u32>(e.index), a, b);
        }
      }

      if (ptr != assign.size()) {
        if (err) *err = "IntervalTreeSamplingBaseline::Sample: internal error (not all assignments consumed)";
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
      if (err) *err = "IntervalTreeSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    // Self-contained enumerator (its own active state).
    return std::make_unique<detail::IntervalTreeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                        join::SideTieBreak::RBeforeS);
  }

 private:
  const DatasetT* ds_{nullptr};
  bool built_{false};

  // Preprocessed sweep events.
  std::vector<join::Event> events_;
  std::vector<i32> start_id_of_event_;  // -1 for END events
  std::vector<usize> start_event_pos_;  // start-id -> event position

  // y-domain (unique y-lower values).
  std::vector<T> y_coords_;

  // Ranks per rectangle.
  std::vector<u32> ylo_rank_r_;
  std::vector<u32> yhi_lb_rank_r_;
  std::vector<u32> ylo_rank_s_;
  std::vector<u32> yhi_lb_rank_s_;

  // Active interval samplers (per side).
  detail::IntervalTreeSampler2D active_r_;
  detail::IntervalTreeSampler2D active_s_;

  // Phase1 weights per START event.
  std::vector<u64> w_total_;
  std::vector<u64> w_a_;
  std::vector<u64> w_b_;

  u64 W_{0};
  bool weights_valid_{false};
};

}  // namespace interval_tree
}  // namespace baselines
}  // namespace sjs
