#pragma once
// sjs/baselines/interval_tree/sampling.h
//
// IntervalTree Baseline (v2.0) — Variant::Sampling
// -------------------------------------------------
// This implementation is written to match the design in:
//   docs/Baseline/IntervalTree Baseline v2.0.md
//
// Problem: Given two sets of axis-aligned half-open rectangles in R^2,
// sample t i.i.d. uniform intersecting cross-pairs without materializing J.
//
// High-level algorithm (Sampling-IT):
//   Phase1: one x-sweep (axis 0) that *only counts* for every START event e:
//           w_e^A (A-type), w_e^B (B-type), w_e = w_e^A + w_e^B, and W=|J|.
//   Phase2: build an event-level alias table over START events with weights w_e,
//           then for each output slot independently sample (event e, type A/B)
//           and store it into per-event slot lists S_e^A / S_e^B.
//   Phase3: a second x-sweep; at each START event e, batch-sample exactly
//           |S_e^A| objects from A and |S_e^B| objects from B, and fill outputs.
//
// Key design details enforced here:
//   * Half-open semantics => END events must be processed before START at the same x.
//   * START tie-break must be deterministic (side order R before S, then id).
//   * Overlap([a,b)) is decomposed as a disjoint union A ⊎ B:
//         A: l <= a < r           (stabbing query at point a)
//         B: a < l < b            (left endpoint in open interval (a,b))
//   * IT-A compresses *all y endpoints* {Ly,Ry} and uses an atomic-cell segment tree.
//   * IT-B uses per-relation keys (Ly, id) to correctly realize strict inequalities.
//
// Notes:
//   - This baseline currently supports Dim==2 only.
//   - The code below is header-only and is used by enum_sampling.h and adaptive.h.

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/join/sweep_events.h"   // join::Side, join::EventKind, join::JoinStats
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
// Sweep event (internal): strict order matching the v2.0 design.
// --------------------------

template <class T>
struct SweepEvent {
  T x;
  join::EventKind kind;  // Start / End
  join::Side side;       // R / S
  u32 index;             // index within the side (fits in u32 by construction)
};

inline constexpr int KindOrder(join::EventKind k) noexcept {
  // Half-open semantics: END before START at the same x.
  return (k == join::EventKind::End) ? 0 : 1;
}

inline constexpr int SideOrder(join::Side s) noexcept {
  // START tie-break: R before S (matches "class priority" in the markdown).
  return (s == join::Side::R) ? 0 : 1;
}

template <class T>
inline bool SweepEventLess(const SweepEvent<T>& a, const SweepEvent<T>& b) noexcept {
  if (a.x < b.x) return true;
  if (b.x < a.x) return false;

  const int ka = KindOrder(a.kind);
  const int kb = KindOrder(b.kind);
  if (ka != kb) return ka < kb;

  // Only START events need the side-order tie-break to make the "later START"
  // assignment deterministic for pairs whose Lx ties.
  if (a.kind == join::EventKind::Start) {
    const int sa = SideOrder(a.side);
    const int sb = SideOrder(b.side);
    if (sa != sb) return sa < sb;
  }

  // Final deterministic tie-break by (side,index).
  const int sa = SideOrder(a.side);
  const int sb = SideOrder(b.side);
  if (sa != sb) return sa < sb;
  return a.index < b.index;
}

template <int Dim, class T>
std::vector<SweepEvent<T>> BuildSweepEvents2D(const Relation<Dim, T>& R, const Relation<Dim, T>& S) {
  static_assert(Dim == 2, "BuildSweepEvents2D only supports Dim==2");

  const usize nR = R.Size();
  const usize nS = S.Size();

  std::vector<SweepEvent<T>> ev;
  ev.reserve(2 * (nR + nS));

  for (usize i = 0; i < nR; ++i) {
    const auto& b = R.boxes[i];
    ev.push_back(SweepEvent<T>{b.lo.v[0], join::EventKind::Start, join::Side::R, static_cast<u32>(i)});
    ev.push_back(SweepEvent<T>{b.hi.v[0], join::EventKind::End, join::Side::R, static_cast<u32>(i)});
  }
  for (usize i = 0; i < nS; ++i) {
    const auto& b = S.boxes[i];
    ev.push_back(SweepEvent<T>{b.lo.v[0], join::EventKind::Start, join::Side::S, static_cast<u32>(i)});
    ev.push_back(SweepEvent<T>{b.hi.v[0], join::EventKind::End, join::Side::S, static_cast<u32>(i)});
  }

  std::sort(ev.begin(), ev.end(), SweepEventLess<T>);
  return ev;
}

// --------------------------
// IT-B key: (Ly, id) for strict open-interval queries a < Ly < b.
// --------------------------

template <class T>
struct Key2 {
  T y;
  u32 id;
};

template <class T>
inline bool Key2Less(const Key2<T>& a, const Key2<T>& b) noexcept {
  if (a.y < b.y) return true;
  if (b.y < a.y) return false;
  return a.id < b.id;
}

// --------------------------
// Bucket primitives: dense vector + O(1) swap-delete using per-item placements.
// --------------------------

struct BucketItem {
  u32 handle;   // rectangle/interval handle (0..n_handles)
  u32 backref;  // which placement record to use for this (handle,node) membership
};

struct Placement {
  u32 node = 0;
  u32 pos = 0;
  u32 backref = 0;
};

// --------------------------
// IT-A: stabbing segment tree (canonical cover buckets)
// --------------------------

class StabbingSegTree {
 public:
  struct Node {
    std::vector<BucketItem> items;
  };

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

    if (m_ == 0) {
      // Degenerate (empty y-domain). All operations should effectively be no-ops.
      return;
    }

    p_ = 1;
    while (p_ < m_) p_ <<= 1;

    // Canonical cover size is O(log p). Keep a small safety margin.
    const u32 logp = static_cast<u32>(std::max<int>(1, 32 - __builtin_clz(p_)));
    max_refs_ = 2 * logp + 4;

    nodes_.resize(static_cast<usize>(2 * p_));
    placements_.resize(static_cast<usize>(n_handles_) * max_refs_);
    placement_size_.assign(static_cast<usize>(n_handles_), 0);

    node_touched_.assign(static_cast<usize>(2 * p_), 0);
    touched_nodes_.reserve(2 * logp + 8);
  }

  void ResetToEmpty() {
    // Clear all buckets that were touched since last reset.
    for (u32 node : touched_nodes_) {
      nodes_[static_cast<usize>(node)].items.clear();
      node_touched_[static_cast<usize>(node)] = 0;
    }
    touched_nodes_.clear();

    std::fill(placement_size_.begin(), placement_size_.end(), 0u);
  }

  void Insert(u32 handle, u32 lo, u32 hi) {
    if (m_ == 0) return;
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(lo < hi);
    SJS_DASSERT(hi <= m_);
    SJS_DASSERT(placement_size_[static_cast<usize>(handle)] == 0);

    u32& ps = placement_size_[static_cast<usize>(handle)];
    InsertRec(/*node=*/1, /*l=*/0, /*r=*/p_, lo, hi, handle, &ps);
    SJS_DASSERT(ps <= max_refs_);
  }

  void Erase(u32 handle) {
    if (m_ == 0) return;
    SJS_DASSERT(handle < n_handles_);
    u32& ps = placement_size_[static_cast<usize>(handle)];
    // Erase all recorded placements for this handle.
    for (u32 i = 0; i < ps; ++i) {
      const usize pidx = PlacementIndex(handle, i);
      const Placement pl = placements_[pidx];
      RemoveFromNode(pl.node, pl.pos);
    }
    ps = 0;
  }

  u64 Count(u32 q) const {
    if (m_ == 0) return 0;
    SJS_DASSERT(q < m_);
    u64 sum = 0;
    u32 node = 1;
    u32 l = 0, r = p_;
    while (true) {
      sum += nodes_[static_cast<usize>(node)].items.size();
      if (r - l == 1) break;
      const u32 mid = (l + r) >> 1;
      if (q < mid) {
        node = node << 1;
        r = mid;
      } else {
        node = (node << 1) | 1;
        l = mid;
      }
    }
    return sum;
  }

  void Report(u32 q, std::vector<u32>* out) const {
    if (!out) return;
    if (m_ == 0) return;
    SJS_DASSERT(q < m_);
    u32 node = 1;
    u32 l = 0, r = p_;
    while (true) {
      const auto& vec = nodes_[static_cast<usize>(node)].items;
      for (const auto& it : vec) out->push_back(it.handle);
      if (r - l == 1) break;
      const u32 mid = (l + r) >> 1;
      if (q < mid) {
        node = node << 1;
        r = mid;
      } else {
        node = (node << 1) | 1;
        l = mid;
      }
    }
  }

  bool Sample(u32 q, u32 k, Rng* rng, std::vector<u32>* out) const {
    if (!out) return false;
    out->clear();
    if (k == 0) return true;
    if (m_ == 0) return false;
    if (!rng) return false;
    SJS_DASSERT(q < m_);

    // Collect buckets along root-to-leaf path and their sizes.
    std::vector<u32> nodes;
    std::vector<u64> weights;

    u32 node = 1;
    u32 l = 0, r = p_;
    while (true) {
      const u64 w = static_cast<u64>(nodes_[static_cast<usize>(node)].items.size());
      if (w > 0) {
        nodes.push_back(node);
        weights.push_back(w);
      }
      if (r - l == 1) break;
      const u32 mid = (l + r) >> 1;
      if (q < mid) {
        node = node << 1;
        r = mid;
      } else {
        node = (node << 1) | 1;
        l = mid;
      }
    }

    if (weights.empty()) return false;

    sampling::AliasTable alias;
    std::string tmp_err;
    if (!alias.BuildFromU64(Span<const u64>(weights), &tmp_err)) {
      return false;
    }

    out->reserve(static_cast<usize>(k));
    for (u32 i = 0; i < k; ++i) {
      const usize idx = alias.Sample(rng);
      const u32 nd = nodes[idx];
      const auto& vec = nodes_[static_cast<usize>(nd)].items;
      const u32 pos = static_cast<u32>(rng->UniformU64(vec.size()));
      out->push_back(vec[static_cast<usize>(pos)].handle);
    }
    return true;
  }

 private:
  usize PlacementIndex(u32 handle, u32 backref) const noexcept { return static_cast<usize>(handle) * max_refs_ + backref; }

  void MarkTouched(u32 node) {
    if (node_touched_[static_cast<usize>(node)]) return;
    node_touched_[static_cast<usize>(node)] = 1;
    touched_nodes_.push_back(node);
  }

  void AddToNode(u32 node, u32 handle, u32 backref, u32* ps) {
    MarkTouched(node);

    Node& nd = nodes_[static_cast<usize>(node)];
    const u32 pos = static_cast<u32>(nd.items.size());
    nd.items.push_back(BucketItem{handle, backref});

    const usize pidx = PlacementIndex(handle, backref);
    placements_[pidx] = Placement{node, pos, backref};

    (*ps)++;
  }

  void InsertRec(u32 node,
                 u32 l,
                 u32 r,
                 u32 ql,
                 u32 qr,
                 u32 handle,
                 u32* ps) {
    if (ql <= l && r <= qr) {
      AddToNode(node, handle, *ps, ps);
      return;
    }
    const u32 mid = (l + r) >> 1;
    if (ql < mid) InsertRec(node << 1, l, mid, ql, qr, handle, ps);
    if (qr > mid) InsertRec((node << 1) | 1, mid, r, ql, qr, handle, ps);
  }

  void RemoveFromNode(u32 node, u32 pos) {
    Node& nd = nodes_[static_cast<usize>(node)];
    auto& vec = nd.items;
    SJS_DASSERT(pos < vec.size());

    const BucketItem moved = vec.back();
    vec[pos] = moved;
    vec.pop_back();

    // Update placement of the moved item.
    const usize pidx = PlacementIndex(moved.handle, moved.backref);
    placements_[pidx].pos = pos;
  }

  u32 n_handles_{0};
  u32 m_{0};       // number of atomic cells (leaf segments)
  u32 p_{0};       // power-of-two leaf base
  u32 max_refs_{0};

  std::vector<Node> nodes_;
  std::vector<Placement> placements_;
  std::vector<u32> placement_size_;

  std::vector<u32> touched_nodes_;
  std::vector<u8> node_touched_;
};

// --------------------------
// IT-B: start-in-range segment tree (root-to-leaf path buckets)
// --------------------------

class StartRangeSegTree {
 public:
  struct Node {
    std::vector<BucketItem> items;
  };

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

    const u32 logp = static_cast<u32>(std::max<int>(1, 32 - __builtin_clz(p_)));
    max_refs_ = logp + 2;  // root-to-leaf length <= logp+1

    nodes_.resize(static_cast<usize>(2 * p_));
    placements_.resize(static_cast<usize>(n_handles_) * max_refs_);
    placement_size_.assign(static_cast<usize>(n_handles_), 0);

    node_touched_.assign(static_cast<usize>(2 * p_), 0);
    touched_nodes_.reserve(max_refs_);
  }

  void ResetToEmpty() {
    for (u32 node : touched_nodes_) {
      nodes_[static_cast<usize>(node)].items.clear();
      node_touched_[static_cast<usize>(node)] = 0;
    }
    touched_nodes_.clear();

    std::fill(placement_size_.begin(), placement_size_.end(), 0u);
  }

  void Insert(u32 handle, u32 pos) {
    if (m_ == 0) return;
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(pos < m_);
    SJS_DASSERT(placement_size_[static_cast<usize>(handle)] == 0);

    u32& ps = placement_size_[static_cast<usize>(handle)];
    u32 node = 1;
    u32 l = 0, r = p_;
    while (true) {
      AddToNode(node, handle, ps, &ps);
      if (r - l == 1) break;
      const u32 mid = (l + r) >> 1;
      if (pos < mid) {
        node = node << 1;
        r = mid;
      } else {
        node = (node << 1) | 1;
        l = mid;
      }
    }
    SJS_DASSERT(ps <= max_refs_);
  }

  void Erase(u32 handle) {
    if (m_ == 0) return;
    SJS_DASSERT(handle < n_handles_);
    u32& ps = placement_size_[static_cast<usize>(handle)];
    for (u32 i = 0; i < ps; ++i) {
      const usize pidx = PlacementIndex(handle, i);
      const Placement pl = placements_[pidx];
      RemoveFromNode(pl.node, pl.pos);
    }
    ps = 0;
  }

  u64 CountRange(u32 ql, u32 qr) const {
    if (m_ == 0) return 0;
    if (ql >= qr) return 0;
    SJS_DASSERT(qr <= m_);

    u64 sum = 0;
    CountRangeRec(1, 0, p_, ql, qr, &sum);
    return sum;
  }

  void ReportRange(u32 ql, u32 qr, std::vector<u32>* out) const {
    if (!out) return;
    if (m_ == 0) return;
    if (ql >= qr) return;
    SJS_DASSERT(qr <= m_);
    ReportRangeRec(1, 0, p_, ql, qr, out);
  }

  bool SampleRange(u32 ql, u32 qr, u32 k, Rng* rng, std::vector<u32>* out) const {
    if (!out) return false;
    out->clear();
    if (k == 0) return true;
    if (m_ == 0) return false;
    if (!rng) return false;
    if (ql >= qr) return false;
    SJS_DASSERT(qr <= m_);

    // Collect canonical decomposition nodes and their bucket sizes.
    std::vector<u32> nodes;
    std::vector<u64> weights;
    CollectDecomp(1, 0, p_, ql, qr, &nodes, &weights);

    if (weights.empty()) return false;

    sampling::AliasTable alias;
    std::string tmp_err;
    if (!alias.BuildFromU64(Span<const u64>(weights), &tmp_err)) {
      return false;
    }

    out->reserve(static_cast<usize>(k));
    for (u32 i = 0; i < k; ++i) {
      const usize idx = alias.Sample(rng);
      const u32 nd = nodes[idx];
      const auto& vec = nodes_[static_cast<usize>(nd)].items;
      const u32 pos = static_cast<u32>(rng->UniformU64(vec.size()));
      out->push_back(vec[static_cast<usize>(pos)].handle);
    }
    return true;
  }

 private:
  usize PlacementIndex(u32 handle, u32 backref) const noexcept { return static_cast<usize>(handle) * max_refs_ + backref; }

  void MarkTouched(u32 node) {
    if (node_touched_[static_cast<usize>(node)]) return;
    node_touched_[static_cast<usize>(node)] = 1;
    touched_nodes_.push_back(node);
  }

  void AddToNode(u32 node, u32 handle, u32 backref, u32* ps) {
    MarkTouched(node);

    Node& nd = nodes_[static_cast<usize>(node)];
    const u32 pos = static_cast<u32>(nd.items.size());
    nd.items.push_back(BucketItem{handle, backref});

    const usize pidx = PlacementIndex(handle, backref);
    placements_[pidx] = Placement{node, pos, backref};

    (*ps)++;
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

  void CountRangeRec(u32 node, u32 l, u32 r, u32 ql, u32 qr, u64* sum) const {
    if (qr <= l || r <= ql) return;
    if (ql <= l && r <= qr) {
      *sum += nodes_[static_cast<usize>(node)].items.size();
      return;
    }
    const u32 mid = (l + r) >> 1;
    CountRangeRec(node << 1, l, mid, ql, qr, sum);
    CountRangeRec((node << 1) | 1, mid, r, ql, qr, sum);
  }

  void ReportRangeRec(u32 node, u32 l, u32 r, u32 ql, u32 qr, std::vector<u32>* out) const {
    if (qr <= l || r <= ql) return;
    if (ql <= l && r <= qr) {
      const auto& vec = nodes_[static_cast<usize>(node)].items;
      for (const auto& it : vec) out->push_back(it.handle);
      return;
    }
    const u32 mid = (l + r) >> 1;
    ReportRangeRec(node << 1, l, mid, ql, qr, out);
    ReportRangeRec((node << 1) | 1, mid, r, ql, qr, out);
  }

  void CollectDecomp(u32 node,
                     u32 l,
                     u32 r,
                     u32 ql,
                     u32 qr,
                     std::vector<u32>* out_nodes,
                     std::vector<u64>* out_weights) const {
    if (qr <= l || r <= ql) return;
    if (ql <= l && r <= qr) {
      const u64 w = static_cast<u64>(nodes_[static_cast<usize>(node)].items.size());
      if (w > 0) {
        out_nodes->push_back(node);
        out_weights->push_back(w);
      }
      return;
    }
    const u32 mid = (l + r) >> 1;
    CollectDecomp(node << 1, l, mid, ql, qr, out_nodes, out_weights);
    CollectDecomp((node << 1) | 1, mid, r, ql, qr, out_nodes, out_weights);
  }

  u32 n_handles_{0};
  u32 m_{0};       // number of leaf points in B-rank domain (== n_handles)
  u32 p_{0};
  u32 max_refs_{0};

  std::vector<Node> nodes_;
  std::vector<Placement> placements_;
  std::vector<u32> placement_size_;

  std::vector<u32> touched_nodes_;
  std::vector<u8> node_touched_;
};

// --------------------------
// IntervalTreeSampler2D: IT-A + IT-B (per side)
// --------------------------

template <class T>
class IntervalTreeSampler2D {
 public:
  void Clear() {
    n_handles_ = 0;
    m_a_ = 0;
    keys_sorted_.clear();
    b_rank_by_handle_.clear();
    active_.clear();
    active_size_ = 0;
    a_.Clear();
    b_.Clear();
  }

  // keys: vector of size n_handles, where keys[i] = (Ly(handle), id(handle)).
  // In this baseline we use id(handle)=handle index for uniqueness.
  void Init(u32 n_handles, u32 m_a_domain, std::vector<Key2<T>> keys) {
    Clear();
    n_handles_ = n_handles;
    m_a_ = m_a_domain;

    // Build IT-A.
    a_.Init(n_handles_, m_a_);

    // Build IT-B (rank domain size is n_handles_).
    b_.Init(n_handles_, n_handles_);

    // Store and sort keys for strict bounds.
    keys_sorted_ = std::move(keys);
    if (keys_sorted_.size() != static_cast<usize>(n_handles_)) {
      // Inconsistent input; keep an empty structure to avoid UB.
      keys_sorted_.clear();
      b_rank_by_handle_.assign(static_cast<usize>(n_handles_), 0u);
    } else {
      std::sort(keys_sorted_.begin(), keys_sorted_.end(), Key2Less<T>);
      b_rank_by_handle_.assign(static_cast<usize>(n_handles_), 0u);
      for (u32 rank = 0; rank < n_handles_; ++rank) {
        const u32 h = keys_sorted_[static_cast<usize>(rank)].id;
        SJS_DASSERT(h < n_handles_);
        b_rank_by_handle_[static_cast<usize>(h)] = rank;
      }
    }

    active_.assign(static_cast<usize>(n_handles_), 0u);
    active_size_ = 0;
  }

  void ResetToEmpty() {
    a_.ResetToEmpty();
    b_.ResetToEmpty();
    std::fill(active_.begin(), active_.end(), 0u);
    active_size_ = 0;
  }

  u32 ActiveSize() const noexcept { return active_size_; }

  void Insert(u32 handle, u32 ylo_a_idx, u32 yhi_a_idx) {
    SJS_DASSERT(handle < n_handles_);
    if (active_[static_cast<usize>(handle)]) return;  // already active (should not happen in a correct sweep)
    active_[static_cast<usize>(handle)] = 1;
    ++active_size_;

    a_.Insert(handle, ylo_a_idx, yhi_a_idx);

    // IT-B: insert at rank of (Ly, id).
    if (n_handles_ > 0) {
      const u32 pos = b_rank_by_handle_[static_cast<usize>(handle)];
      b_.Insert(handle, pos);
    }
  }

  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_handles_);
    if (!active_[static_cast<usize>(handle)]) return;  // not active (no-op)
    active_[static_cast<usize>(handle)] = 0;
    if (active_size_ > 0) --active_size_;

    a_.Erase(handle);
    b_.Erase(handle);
  }

  // A-type (stabbing): l <= a < r. Here a is represented as an atomic-cell index.
  u64 CountA(u32 a_idx) const { return a_.Count(a_idx); }

  void ReportA(u32 a_idx, std::vector<u32>* out) const { a_.Report(a_idx, out); }

  bool SampleA(u32 a_idx, u32 k, Rng* rng, std::vector<u32>* out) const { return a_.Sample(a_idx, k, rng, out); }

  // B-type: a < l < b, where a and b are y-coordinates.
  // Implemented via key ranks using upper/lower bounds on (Ly,id).
  u64 CountB(T a, T b) const {
    if (n_handles_ == 0) return 0;
    const auto [L, R] = BRange(a, b);
    return b_.CountRange(L, R);
  }

  void ReportB(T a, T b, std::vector<u32>* out) const {
    if (!out) return;
    if (n_handles_ == 0) return;
    const auto [L, R] = BRange(a, b);
    b_.ReportRange(L, R, out);
  }

  bool SampleB(T a, T b, u32 k, Rng* rng, std::vector<u32>* out) const {
    if (n_handles_ == 0) {
      if (out) out->clear();
      return false;
    }
    const auto [L, R] = BRange(a, b);
    if (L >= R) {
      if (out) out->clear();
      return false;
    }
    return b_.SampleRange(L, R, k, rng, out);
  }

 private:
  // Compute rank interval [L,R) s.t. a < Ly < b (open interval).
  std::pair<u32, u32> BRange(T a, T b) const {
    // If a>=b, (a,b) is empty.
    if (!(a < b)) return {0u, 0u};

    // L = upper_bound((a, +inf))  => first key with y > a.
    const Key2<T> kL{a, std::numeric_limits<u32>::max()};
    const auto itL = std::upper_bound(keys_sorted_.begin(), keys_sorted_.end(), kL,
                                      [](const Key2<T>& x, const Key2<T>& y) { return Key2Less<T>(x, y); });
    const u32 L = static_cast<u32>(itL - keys_sorted_.begin());

    // R = lower_bound((b, -inf)) => first key with y >= b.
    const Key2<T> kR{b, 0u};
    const auto itR = std::lower_bound(keys_sorted_.begin(), keys_sorted_.end(), kR,
                                      [](const Key2<T>& x, const Key2<T>& y) { return Key2Less<T>(x, y); });
    const u32 R = static_cast<u32>(itR - keys_sorted_.begin());

    // Clamp (should already be within [0,n_handles]).
    return {std::min<u32>(L, n_handles_), std::min<u32>(R, n_handles_)};
  }

  u32 n_handles_{0};
  u32 m_a_{0};

  std::vector<Key2<T>> keys_sorted_;        // sorted by (Ly, id)
  std::vector<u32> b_rank_by_handle_;       // handle -> rank in keys_sorted_
  std::vector<u8> active_;                  // handle -> active flag
  u32 active_size_{0};

  StabbingSegTree a_;
  StartRangeSegTree b_;
};

// --------------------------
// Slot plan for Phase2 (event alias + per-event slot lists)
// We store S_e^A and S_e^B as flat arrays with per-event offsets.
// --------------------------

struct SlotPlan {
  // Offsets are size (E+1); slots for event e are in [off[e], off[e+1]).
  std::vector<u32> a_off;
  std::vector<u32> b_off;
  std::vector<u32> a_slots;
  std::vector<u32> b_slots;

  void Clear() {
    a_off.clear();
    b_off.clear();
    a_slots.clear();
    b_slots.clear();
  }
};

// Build per-event slot lists for t samples.
// w_total/w_a/w_b are per-START-event weights with w_total[e]=w_a[e]+w_b[e].
//
// This matches the markdown Phase2 semantics:
//   - sample event e with prob w_total[e]/W
//   - then sample type A vs B with prob w_a[e]/w_total[e] and w_b[e]/w_total[e]
//   - record the output slot index into S_e^A or S_e^B
inline bool BuildSlotPlan(const std::vector<u64>& w_total,
                          const std::vector<u64>& w_a,
                          const std::vector<u64>& w_b,
                          u32 t,
                          Rng* rng,
                          SlotPlan* plan,
                          std::string* err) {
  if (!plan) return false;
  plan->Clear();
  if (t == 0) return true;
  if (!rng) {
    if (err) *err = "BuildSlotPlan: rng is null";
    return false;
  }

  const usize E = w_total.size();
  if (w_a.size() != E || w_b.size() != E) {
    if (err) *err = "BuildSlotPlan: inconsistent weight array sizes";
    return false;
  }

  // Filter out zero-weight events, as required by the v2.0 design.
  std::vector<u32> nz_sids;
  std::vector<u64> nz_w;
  nz_sids.reserve(E);
  nz_w.reserve(E);
  for (u32 sid = 0; sid < static_cast<u32>(E); ++sid) {
    const u64 w = w_total[static_cast<usize>(sid)];
    if (w > 0) {
      nz_sids.push_back(sid);
      nz_w.push_back(w);
    }
  }
  if (nz_w.empty()) {
    if (err) *err = "BuildSlotPlan: no positive-weight events (W==0?)";
    return false;
  }

  sampling::AliasTable alias;
  if (!alias.BuildFromU64(Span<const u64>(nz_w), err)) {
    if (err && err->empty()) *err = "BuildSlotPlan: failed to build event alias table";
    return false;
  }

  // First pass: sample (sid, pat) for each output slot and count per-event needs.
  std::vector<u32> sid_of_slot(static_cast<usize>(t));
  std::vector<u8> pat_of_slot(static_cast<usize>(t));  // 0=A, 1=B

  std::vector<u32> cnt_a(E, 0u);
  std::vector<u32> cnt_b(E, 0u);

  for (u32 j = 0; j < t; ++j) {
    const usize picked = alias.Sample(rng);
    const u32 sid = nz_sids[picked];

    const u64 w = w_total[static_cast<usize>(sid)];
    const u64 wa = w_a[static_cast<usize>(sid)];
    const u64 wb = w_b[static_cast<usize>(sid)];
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

    sid_of_slot[static_cast<usize>(j)] = sid;
    pat_of_slot[static_cast<usize>(j)] = pat;

    if (pat == 0) ++cnt_a[static_cast<usize>(sid)];
    else ++cnt_b[static_cast<usize>(sid)];
  }

  // Prefix sums -> offsets.
  plan->a_off.assign(E + 1, 0u);
  plan->b_off.assign(E + 1, 0u);
  for (usize sid = 0; sid < E; ++sid) {
    plan->a_off[sid + 1] = plan->a_off[sid] + cnt_a[sid];
    plan->b_off[sid + 1] = plan->b_off[sid] + cnt_b[sid];
  }
  const u32 total_a = plan->a_off[E];
  const u32 total_b = plan->b_off[E];
  plan->a_slots.assign(static_cast<usize>(total_a), 0u);
  plan->b_slots.assign(static_cast<usize>(total_b), 0u);

  // Second pass: materialize flat slot lists.
  std::vector<u32> write_a = plan->a_off;
  std::vector<u32> write_b = plan->b_off;

  for (u32 j = 0; j < t; ++j) {
    const u32 sid = sid_of_slot[static_cast<usize>(j)];
    const u8 pat = pat_of_slot[static_cast<usize>(j)];

    if (pat == 0) {
      const u32 pos = write_a[static_cast<usize>(sid)]++;
      plan->a_slots[static_cast<usize>(pos)] = j;
    } else {
      const u32 pos = write_b[static_cast<usize>(sid)]++;
      plan->b_slots[static_cast<usize>(pos)] = j;
    }
  }

  return true;
}

// --------------------------
// Deterministic join enumerator: sweep + IntervalTreeSampler2D
// Enumerates each intersecting (R,S) pair exactly once (later START).
// --------------------------

template <int Dim, class T>
class IntervalTreeJoinEnumerator final : public IJoinEnumerator {
 public:
  using RelationT = Relation<Dim, T>;

  IntervalTreeJoinEnumerator(const RelationT* R, const RelationT* S)
      : R_(R), S_(S) {
    SJS_DASSERT(R_ != nullptr && S_ != nullptr);
    BuildPreproc();
    Reset();
  }

  void Reset() override {
    ev_pos_ = 0;
    in_emit_ = false;
    cand_.clear();
    cand_pos_ = 0;
    cur_side_ = join::Side::R;
    cur_index_ = 0;

    active_r_.ResetToEmpty();
    active_s_.ResetToEmpty();
    active_r_size_ = 0;
    active_s_size_ = 0;

    stats_.Reset();
    stats_.num_events = static_cast<u64>(events_.size());
  }

  bool Next(PairId* out) override {
    if (!out) return false;

    // If we finished emitting candidates for a START event, insert that START rectangle now.
    if (in_emit_ && cand_pos_ >= cand_.size()) {
      FinishStartEvent();
    }

    if (in_emit_) {
      EmitOne(out);
      return true;
    }

    while (ev_pos_ < events_.size()) {
      const auto& e = events_[ev_pos_++];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) {
          active_r_.Erase(e.index);
          if (active_r_size_ > 0) --active_r_size_;
        } else {
          active_s_.Erase(e.index);
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
        const u32 a_idx = ylo_idx_r_[static_cast<usize>(cur_index_)];
        const T a = R_->boxes[static_cast<usize>(cur_index_)].lo.v[1];
        const T b = R_->boxes[static_cast<usize>(cur_index_)].hi.v[1];
        active_s_.ReportA(a_idx, &cand_);
        active_s_.ReportB(a, b, &cand_);
      } else {
        const u32 a_idx = ylo_idx_s_[static_cast<usize>(cur_index_)];
        const T a = S_->boxes[static_cast<usize>(cur_index_)].lo.v[1];
        const T b = S_->boxes[static_cast<usize>(cur_index_)].hi.v[1];
        active_r_.ReportA(a_idx, &cand_);
        active_r_.ReportB(a, b, &cand_);
      }

      if (cand_.empty()) {
        // No output for this START; still need to insert it and continue.
        InsertCurrentStart();
        continue;
      }

      in_emit_ = true;
      EmitOne(out);
      return true;
    }

    return false;
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  void BuildPreproc() {
    static_assert(Dim == 2, "IntervalTreeJoinEnumerator only supports Dim==2");

    // Build events with strict ordering.
    events_ = BuildSweepEvents2D<Dim, T>(*R_, *S_);

    // Build y-endpoints {Ly,Ry} across both relations.
    y_endpoints_.clear();
    y_endpoints_.reserve(2 * (R_->Size() + S_->Size()));
    for (const auto& b : R_->boxes) {
      y_endpoints_.push_back(b.lo.v[1]);
      y_endpoints_.push_back(b.hi.v[1]);
    }
    for (const auto& b : S_->boxes) {
      y_endpoints_.push_back(b.lo.v[1]);
      y_endpoints_.push_back(b.hi.v[1]);
    }
    std::sort(y_endpoints_.begin(), y_endpoints_.end());
    y_endpoints_.erase(std::unique(y_endpoints_.begin(), y_endpoints_.end()), y_endpoints_.end());

    const u32 m_end = static_cast<u32>(y_endpoints_.size());
    const u32 m_a = (m_end >= 2) ? (m_end - 1) : 0;

    auto lb_end = [&](T v) -> u32 {
      const auto it = std::lower_bound(y_endpoints_.begin(), y_endpoints_.end(), v);
      return static_cast<u32>(it - y_endpoints_.begin());
    };

    ylo_idx_r_.resize(R_->Size());
    yhi_idx_r_.resize(R_->Size());
    for (usize i = 0; i < R_->Size(); ++i) {
      const auto& b = R_->boxes[i];
      const u32 lo = lb_end(b.lo.v[1]);
      const u32 hi = lb_end(b.hi.v[1]);
      SJS_DASSERT(lo < hi);
      SJS_DASSERT(lo < m_a || m_a == 0);
      SJS_DASSERT(hi <= m_a);
      ylo_idx_r_[i] = lo;
      yhi_idx_r_[i] = hi;
    }

    ylo_idx_s_.resize(S_->Size());
    yhi_idx_s_.resize(S_->Size());
    for (usize i = 0; i < S_->Size(); ++i) {
      const auto& b = S_->boxes[i];
      const u32 lo = lb_end(b.lo.v[1]);
      const u32 hi = lb_end(b.hi.v[1]);
      SJS_DASSERT(lo < hi);
      SJS_DASSERT(lo < m_a || m_a == 0);
      SJS_DASSERT(hi <= m_a);
      ylo_idx_s_[i] = lo;
      yhi_idx_s_[i] = hi;
    }

    // Build IT-B keys per relation: (Ly, id) where id is the handle index.
    std::vector<Key2<T>> keys_r;
    keys_r.reserve(R_->Size());
    for (usize i = 0; i < R_->Size(); ++i) keys_r.push_back(Key2<T>{R_->boxes[i].lo.v[1], static_cast<u32>(i)});

    std::vector<Key2<T>> keys_s;
    keys_s.reserve(S_->Size());
    for (usize i = 0; i < S_->Size(); ++i) keys_s.push_back(Key2<T>{S_->boxes[i].lo.v[1], static_cast<u32>(i)});

    active_r_.Init(static_cast<u32>(R_->Size()), m_a, std::move(keys_r));
    active_s_.Init(static_cast<u32>(S_->Size()), m_a, std::move(keys_s));
  }

  void EmitOne(PairId* out) {
    const u32 other_h = cand_[cand_pos_++];

    if (cur_side_ == join::Side::R) {
      *out = PairId{R_->GetId(static_cast<usize>(cur_index_)), S_->GetId(static_cast<usize>(other_h))};
    } else {
      *out = PairId{R_->GetId(static_cast<usize>(other_h)), S_->GetId(static_cast<usize>(cur_index_))};
    }

    ++stats_.output_pairs;
    ++stats_.candidate_checks;
  }

  void InsertCurrentStart() {
    if (cur_side_ == join::Side::R) {
      active_r_.Insert(cur_index_, ylo_idx_r_[static_cast<usize>(cur_index_)], yhi_idx_r_[static_cast<usize>(cur_index_)]);
      ++active_r_size_;
    } else {
      active_s_.Insert(cur_index_, ylo_idx_s_[static_cast<usize>(cur_index_)], yhi_idx_s_[static_cast<usize>(cur_index_)]);
      ++active_s_size_;
    }
    stats_.active_max_r = std::max<u64>(stats_.active_max_r, active_r_size_);
    stats_.active_max_s = std::max<u64>(stats_.active_max_s, active_s_size_);
  }

  void FinishStartEvent() {
    InsertCurrentStart();
    in_emit_ = false;
    cand_.clear();
    cand_pos_ = 0;
  }

  const RelationT* R_;
  const RelationT* S_;

  // Preprocessed sweep.
  std::vector<SweepEvent<T>> events_;

  // Global y endpoints and per-rect indices for IT-A.
  std::vector<T> y_endpoints_;
  std::vector<u32> ylo_idx_r_;
  std::vector<u32> yhi_idx_r_;
  std::vector<u32> ylo_idx_s_;
  std::vector<u32> yhi_idx_s_;

  // Active interval samplers.
  IntervalTreeSampler2D<T> active_r_;
  IntervalTreeSampler2D<T> active_s_;

  // Streaming state.
  usize ev_pos_{0};
  bool in_emit_{false};
  join::Side cur_side_{join::Side::R};
  u32 cur_index_{0};
  std::vector<u32> cand_;
  usize cand_pos_{0};

  // Simple active sizes for stats.
  u64 active_r_size_{0};
  u64 active_s_size_{0};

  join::JoinStats stats_;
};

}  // namespace detail

// --------------------------
// Baseline: IntervalTree + Sampling (two-pass)
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
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;

    if constexpr (Dim != 2) {
      if (err) *err = "IntervalTreeSamplingBaseline: only Dim==2 is supported";
      return false;
    }

    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    ds_ = &ds;

    if (ds.R.Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        ds.S.Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "IntervalTreeSamplingBaseline::Build: relation size exceeds u32 limit";
      return false;
    }

    // Events (strict ordering).
    {
      auto ph = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = detail::BuildSweepEvents2D<Dim, T>(ds.R, ds.S);
    }

    // START id mapping (event position -> start-id, -1 for END).
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

    // y-endpoints {Ly,Ry} across both relations.
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

    // Precompute IT-A indices (lo,hi) in the atomic-cell domain.
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
          if (err) *err = "IntervalTreeSamplingBaseline::Build: invalid R y-interval (lo>=hi)";
          return false;
        }
        if (m_a_ > 0 && lo >= m_a_) {
          if (err) *err = "IntervalTreeSamplingBaseline::Build: R ylo maps to last endpoint (unexpected)";
          return false;
        }
        if (hi > m_a_) {
          if (err) *err = "IntervalTreeSamplingBaseline::Build: R yhi index out of range";
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
          if (err) *err = "IntervalTreeSamplingBaseline::Build: invalid S y-interval (lo>=hi)";
          return false;
        }
        if (m_a_ > 0 && lo >= m_a_) {
          if (err) *err = "IntervalTreeSamplingBaseline::Build: S ylo maps to last endpoint (unexpected)";
          return false;
        }
        if (hi > m_a_) {
          if (err) *err = "IntervalTreeSamplingBaseline::Build: S yhi index out of range";
          return false;
        }
        ylo_idx_s_[i] = lo;
        yhi_idx_s_[i] = hi;
      }
    }

    // Init active samplers (IT-A + IT-B) for each side.
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

    active_r_.ResetToEmpty();
    active_s_.ResetToEmpty();

    u64 W = 0;

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

      if (e.side == join::Side::R) {
        const u32 a_idx = ylo_idx_r_[static_cast<usize>(e.index)];
        const u32 hi_idx = yhi_idx_r_[static_cast<usize>(e.index)];
        const T a = ds_->R.boxes[static_cast<usize>(e.index)].lo.v[1];
        const T b = ds_->R.boxes[static_cast<usize>(e.index)].hi.v[1];

        const u64 wa = active_s_.CountA(a_idx);
        const u64 wb = active_s_.CountB(a, b);
        const u64 w = wa + wb;

        w_a_[sid] = wa;
        w_b_[sid] = wb;
        w_total_[sid] = w;

        if (w != 0 && W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "IntervalTreeSamplingBaseline::Count: |J| overflow (u64)";
          return false;
        }
        W += w;

        active_r_.Insert(e.index, a_idx, hi_idx);

      } else {  // S START
        const u32 a_idx = ylo_idx_s_[static_cast<usize>(e.index)];
        const u32 hi_idx = yhi_idx_s_[static_cast<usize>(e.index)];
        const T a = ds_->S.boxes[static_cast<usize>(e.index)].lo.v[1];
        const T b = ds_->S.boxes[static_cast<usize>(e.index)].hi.v[1];

        const u64 wa = active_r_.CountA(a_idx);
        const u64 wb = active_r_.CountB(a, b);
        const u64 w = wa + wb;

        w_a_[sid] = wa;
        w_b_[sid] = wb;
        w_total_[sid] = w;

        if (w != 0 && W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "IntervalTreeSamplingBaseline::Count: |J| overflow (u64)";
          return false;
        }
        W += w;

        active_s_.Insert(e.index, a_idx, hi_idx);
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

    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    if (t == 0) return true;

    // Ensure Phase1 weights are available.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    if (W_ == 0) return true;  // empty join

    // Phase2: event alias + per-event slot lists.
    detail::SlotPlan plan;
    {
      auto scoped2 = phases ? phases->Scoped("phase2_slots") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!detail::BuildSlotPlan(w_total_, w_a_, w_b_, t, rng, &plan, err)) return false;
    }

    out->pairs.resize(static_cast<usize>(t));

    // Phase3: second sweep; fulfill slots at each START and then insert the START rectangle.
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
              if (err) *err = "IntervalTreeSamplingBaseline::Sample: empty A-set in phase3 (R-start)";
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
              if (err) *err = "IntervalTreeSamplingBaseline::Sample: empty B-set in phase3 (R-start)";
              return false;
            }
            for (u32 i = 0; i < kB; ++i) {
              const u32 slot = plan.b_slots[static_cast<usize>(b_begin + i)];
              const u32 other = tmp[static_cast<usize>(i)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(e.index)),
                                                           ds_->S.GetId(static_cast<usize>(other))};
            }
          }

          // Insert current START after fulfilling its slots.
          active_r_.Insert(e.index, a_idx, hi_idx);

        } else {  // S START
          const u32 a_idx = ylo_idx_s_[static_cast<usize>(e.index)];
          const u32 hi_idx = yhi_idx_s_[static_cast<usize>(e.index)];
          const T a = ds_->S.boxes[static_cast<usize>(e.index)].lo.v[1];
          const T b = ds_->S.boxes[static_cast<usize>(e.index)].hi.v[1];

          if (kA > 0) {
            if (!active_r_.SampleA(a_idx, kA, rng, &tmp)) {
              if (err) *err = "IntervalTreeSamplingBaseline::Sample: empty A-set in phase3 (S-start)";
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
              if (err) *err = "IntervalTreeSamplingBaseline::Sample: empty B-set in phase3 (S-start)";
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
      if (err) *err = "IntervalTreeSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    return std::make_unique<detail::IntervalTreeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S);
  }

 private:
  const DatasetT* ds_{nullptr};
  bool built_{false};

  // Preprocessed events and START-id mapping.
  std::vector<detail::SweepEvent<T>> events_;
  std::vector<i32> start_id_of_event_;  // -1 for END

  // y endpoints (global) and IT-A indices per rectangle.
  std::vector<T> y_endpoints_;
  u32 m_a_{0};  // number of atomic cells (Y.size()-1)

  std::vector<u32> ylo_idx_r_;
  std::vector<u32> yhi_idx_r_;
  std::vector<u32> ylo_idx_s_;
  std::vector<u32> yhi_idx_s_;

  // Active interval samplers (per side).
  detail::IntervalTreeSampler2D<T> active_r_;
  detail::IntervalTreeSampler2D<T> active_s_;

  // Phase1 weights per START event (sid order).
  std::vector<u64> w_total_;
  std::vector<u64> w_a_;
  std::vector<u64> w_b_;

  u64 W_{0};
  bool weights_valid_{false};
};

}  // namespace interval_tree
}  // namespace baselines
}  // namespace sjs
