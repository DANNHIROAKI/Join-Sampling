#pragma once
// sjs/baselines/ours/sampling.h
//
// "Our Method" — Sampling variant (3-phase uniform sampler).
//
// This header is written to match the design in "Our Method.md":
//
//  - Common preprocessing (all variants share):
//      * Build plane-sweep events on axis 0 with half-open semantics:
//          - sort by coordinate
//          - END before START at the same coordinate
//          - START tie-break by SideTieBreak (R before S by default)
//      * For each side (R/S) build per-pattern indices (skeleton) and keep the
//        active set empty.
//
//  - Sampling variant (Ours):
//      Phase 1: single sweep, for each START event e with box q compute exact
//               w_e^A and w_e^B and w_e = w_e^A + w_e^B. Let W = sum_e w_e = |J|.
//      Phase 2: build event-level alias on (w_e). For each output slot j=1..t:
//               sample event E_j ~ w_e/W; sample pattern G_j in {A,B} with
//               Pr(G_j=A|E_j=e)=w_e^A/w_e; then append slot j into S_e^G.
//      Phase 3: second sweep. When visiting START event e for q:
//               for each pattern g, let k = |S_e^g|. If k>0, call
//                 Sample_g(q*, k)
//               to get k i.i.d. uniform "other-side" rectangles, then fill
//               Ans at those slots, keeping pair order as (R,S).
//
// Dim notes:
//  - The repository currently instantiates Dim=2. This implementation is
//    specialized for Dim=2 but keeps the template parameter so it integrates
//    cleanly with the baseline factory.

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/join/join_enumerator.h"  // for PlaneSweepJoinStream wrapper
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
namespace ours {

namespace detail {

// -----------------------------------------------------------------------------
// Pattern definitions for 2D (y-axis)
// -----------------------------------------------------------------------------
// For half-open intervals [lo,hi):
//   [a0,a1) intersects [b0,b1)  <=>  a0 < b1 AND b0 < a1.
// When sweep axis (x) overlap is guaranteed by active-set membership, a pair
// intersects iff y-intervals intersect.
//
// Partition by comparing the lower endpoints (ties go to Pattern A):
//   Pattern A: y_r_lo <= y_q_lo < y_r_hi
//   Pattern B: y_q_lo < y_r_lo < y_q_hi
// These sets are disjoint and cover all y-intersections.

// -----------------------------------------------------------------------------
// RangePointSegTree (Pattern B): dynamic point range count/sample/report
// -----------------------------------------------------------------------------
// We maintain active points keyed by y_r_lo (rank on the compressed ylo-domain).
//
// Key design goal (per Our Method.md):
//   SampleRange(l,r,k) should cost O(log m + k) for fixed query (l,r), by
//   building a small alias table once and drawing k samples from it.
//
// Technique:
//   - Segment tree over ranks [0, P), P = next power-of-two >= m.
//   - Insert a point into *all* nodes on its leaf-to-root path.
//   - For a query range [l,r), compute the canonical segment-tree cover nodes.
//     These cover nodes are disjoint and each leaf in the range belongs to
//     exactly one cover node; therefore each point in [l,r) appears in exactly
//     one cover node bucket (even though it is stored on multiple ancestors).
//   - To sample uniformly from all points in [l,r): choose a cover node with
//     probability proportional to its bucket size, then choose uniformly from
//     that bucket.
//
// Deletion is O(log m) using swap-delete + backrefs, same as the stabbing tree.

class RangePointSegTree {
 public:
  struct Placement {
    u32 node = 0;  // segment-tree node index
    u32 pos = 0;   // position inside nodes_[node].items
  };

  struct Item {
    u32 handle = 0;
    u32 backref = 0;  // index into handle's placement list
  };

  void Init(u32 num_handles, u32 num_ranks) {
    n_handles_ = num_handles;
    m_ = num_ranks;

    // Next power-of-two.
    p_ = 1;
    while (p_ < m_) p_ <<= 1;

    nodes_.assign(static_cast<usize>(2 * p_), Node{});

    // log2(p_)
    u32 logp = 0;
    for (u32 x = p_; x > 1; x >>= 1) ++logp;
    // Path length leaf->root inclusive.
    max_refs_ = static_cast<u32>(logp + 1);

    placement_size_.assign(static_cast<usize>(n_handles_), 0);
    placements_.assign(static_cast<usize>(n_handles_) * static_cast<usize>(max_refs_), Placement{});
  }

  void Clear() {
    n_handles_ = 0;
    m_ = 0;
    p_ = 0;
    max_refs_ = 0;
    nodes_.clear();
    placements_.clear();
    placement_size_.clear();
  }

  // Keep the skeleton but drop all active points.
  void ResetActive() {
    for (auto& n : nodes_) n.items.clear();
    std::fill(placement_size_.begin(), placement_size_.end(), 0U);
  }

  void Insert(u32 handle, u32 rank) {
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(rank < m_);
    SJS_DASSERT(placement_size_[handle] == 0);

    u32 backref = 0;
    u32 idx = rank + p_;
    while (idx > 0) {
      AddToNode(handle, idx, backref);
      ++backref;
      idx >>= 1;
    }
    SJS_DASSERT(backref <= max_refs_);
    placement_size_[handle] = backref;
  }

  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_handles_);
    const u32 sz = placement_size_[handle];
    for (u32 i = 0; i < sz; ++i) {
      const Placement p = placements_[PlacementIndex(handle, i)];
      RemoveFromNode(p.node, p.pos);
    }
    placement_size_[handle] = 0;
  }

  u64 CountRange(u32 l, u32 r) const {
    if (r <= l || m_ == 0) return 0ULL;
    l = std::min(l, m_);
    r = std::min(r, m_);
    if (r <= l) return 0ULL;

    std::vector<u32> cover;
    Decompose(l, r, &cover);
    u64 total = 0;
    for (u32 node : cover) total += static_cast<u64>(nodes_[node].items.size());
    return total;
  }

  // Sample k handles uniformly from ranks in [l, r).
  // Returns false iff the range is empty.
  bool SampleRange(u32 l, u32 r, u32 k, Rng* rng, std::vector<u32>* out) const {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(out != nullptr);
    out->clear();
    out->reserve(k);

    if (k == 0) return true;

    if (r <= l || m_ == 0) return false;
    l = std::min(l, m_);
    r = std::min(r, m_);
    if (r <= l) return false;

    std::vector<u32> cover;
    std::vector<u64> weights;
    Decompose(l, r, &cover);

    u64 total = 0;
    weights.reserve(cover.size());
    for (u32 node : cover) {
      const u64 w = static_cast<u64>(nodes_[node].items.size());
      weights.push_back(w);
      total += w;
    }

    if (total == 0) return false;

    sampling::AliasTable alias;
    if (!alias.BuildFromU64(Span<const u64>(weights))) {
      return false;
    }

    for (u32 i = 0; i < k; ++i) {
      const usize bi = alias.Sample(rng);
      const u32 node = cover[bi];
      const auto& bucket = nodes_[node].items;
      SJS_DASSERT(!bucket.empty());
      const u32 pos = rng->UniformU32(static_cast<u32>(bucket.size()));
      out->push_back(bucket[pos].handle);
    }
    return true;
  }

  // Report all handles in [l, r) in a deterministic left-to-right cover order.
  void ReportRange(u32 l, u32 r, std::vector<u32>* out) const {
    SJS_DASSERT(out != nullptr);
    if (r <= l || m_ == 0) return;
    l = std::min(l, m_);
    r = std::min(r, m_);
    if (r <= l) return;

    std::vector<u32> cover;
    Decompose(l, r, &cover);
    for (u32 node : cover) {
      const auto& bucket = nodes_[node].items;
      for (const auto& it : bucket) out->push_back(it.handle);
    }
  }

 private:
  struct Node {
    std::vector<Item> items;
  };

  usize PlacementIndex(u32 handle, u32 backref) const {
    return static_cast<usize>(handle) * static_cast<usize>(max_refs_) + static_cast<usize>(backref);
  }

  void AddToNode(u32 handle, u32 node, u32 backref) {
    SJS_DASSERT(node < nodes_.size());
    SJS_DASSERT(backref < max_refs_);

    auto& bucket = nodes_[node].items;
    const u32 pos = static_cast<u32>(bucket.size());

    placements_[PlacementIndex(handle, backref)] = Placement{node, pos};
    bucket.push_back(Item{handle, backref});
  }

  void RemoveFromNode(u32 node, u32 pos) {
    SJS_DASSERT(node < nodes_.size());
    auto& bucket = nodes_[node].items;
    SJS_DASSERT(!bucket.empty());
    SJS_DASSERT(pos < bucket.size());

    const usize last_pos = bucket.size() - 1;
    if (static_cast<usize>(pos) != last_pos) {
      const Item swapped = bucket[last_pos];
      bucket[pos] = swapped;
      placements_[PlacementIndex(swapped.handle, swapped.backref)].pos = pos;
    }
    bucket.pop_back();
  }

  // Canonical cover decomposition of [l,r) (0<=l<=r<=m_).
  // The returned nodes are disjoint and ordered left-to-right.
  void Decompose(u32 l, u32 r, std::vector<u32>* out) const {
    out->clear();
    if (r <= l) return;

    u32 L = l + p_;
    u32 R = r + p_;

    std::vector<u32> left;
    std::vector<u32> right;
    left.reserve(32);
    right.reserve(32);

    while (L < R) {
      if (L & 1U) left.push_back(L++);
      if (R & 1U) right.push_back(--R);
      L >>= 1;
      R >>= 1;
    }

    out->reserve(left.size() + right.size());
    out->insert(out->end(), left.begin(), left.end());
    for (auto it = right.rbegin(); it != right.rend(); ++it) out->push_back(*it);
  }

  u32 n_handles_{0};
  u32 m_{0};
  u32 p_{0};
  u32 max_refs_{0};

  std::vector<Node> nodes_;               // size 2*p_
  std::vector<Placement> placements_;     // size n_handles_*max_refs_
  std::vector<u32> placement_size_;       // per-handle placement count
};

// -----------------------------------------------------------------------------
// StabbingSegTree (Pattern A): dynamic interval stabbing count/sample/report
// -----------------------------------------------------------------------------
// Identical to the previous implementation, with an explicit ResetActive().

class StabbingSegTree {
 public:
  struct Placement {
    u32 node = 0;
    u32 pos = 0;
  };

  struct Item {
    u32 handle = 0;
    u32 backref = 0;
  };

  void Init(u32 num_handles, u32 num_points) {
    n_handles_ = num_handles;
    m_ = num_points;

    p_ = 1;
    while (p_ < m_) p_ <<= 1;

    nodes_.assign(static_cast<usize>(2 * p_), Node{});

    // log2(p_)
    u32 logp = 0;
    for (u32 x = p_; x > 1; x >>= 1) ++logp;

    // Interval range decomposition touches <= 2*log2(p_) nodes.
    max_refs_ = static_cast<u32>(2 * logp + 4);

    placement_size_.assign(static_cast<usize>(n_handles_), 0);
    placements_.assign(static_cast<usize>(n_handles_) * static_cast<usize>(max_refs_), Placement{});
  }

  void Clear() {
    n_handles_ = 0;
    m_ = 0;
    p_ = 0;
    max_refs_ = 0;
    nodes_.clear();
    placements_.clear();
    placement_size_.clear();
  }

  // Keep the skeleton but drop all active intervals.
  void ResetActive() {
    for (auto& n : nodes_) n.items.clear();
    std::fill(placement_size_.begin(), placement_size_.end(), 0U);
  }

  // Insert interval [L, R) on ranks.
  void Insert(u32 handle, u32 L, u32 R) {
    SJS_DASSERT(handle < n_handles_);
    if (R <= L || m_ == 0) return;
    L = std::min(L, m_);
    R = std::min(R, m_);
    if (R <= L) return;

    u32 l = L + p_;
    u32 r = R + p_;

    while (l < r) {
      if (l & 1U) {
        AddToNode(handle, l);
        ++l;
      }
      if (r & 1U) {
        --r;
        AddToNode(handle, r);
      }
      l >>= 1;
      r >>= 1;
    }
  }

  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_handles_);
    const u32 sz = placement_size_[handle];
    for (u32 i = 0; i < sz; ++i) {
      const Placement p = placements_[PlacementIndex(handle, i)];
      RemoveFromNode(p.node, p.pos);
    }
    placement_size_[handle] = 0;
  }

  u64 Count(u32 q) const {
    if (m_ == 0 || q >= m_) return 0ULL;
    u64 total = 0;
    u32 idx = q + p_;
    while (idx > 0) {
      total += static_cast<u64>(nodes_[idx].items.size());
      idx >>= 1;
    }
    return total;
  }

  bool Sample(u32 q, u32 k, Rng* rng, std::vector<u32>* out) const {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(out != nullptr);
    out->clear();
    out->reserve(k);

    if (k == 0) return true;
    if (m_ == 0 || q >= m_) return false;

    // Collect non-empty buckets on root-to-leaf path.
    std::vector<u32> path_nodes;
    std::vector<u64> weights;
    path_nodes.reserve(32);
    weights.reserve(32);

    u32 idx = q + p_;
    u64 total = 0;
    while (idx > 0) {
      const u64 w = static_cast<u64>(nodes_[idx].items.size());
      if (w > 0) {
        path_nodes.push_back(idx);
        weights.push_back(w);
        total += w;
      }
      idx >>= 1;
    }

    if (total == 0) return false;

    sampling::AliasTable alias;
    if (!alias.BuildFromU64(Span<const u64>(weights))) {
      return false;
    }

    for (u32 i = 0; i < k; ++i) {
      const usize bi = alias.Sample(rng);
      const u32 node = path_nodes[bi];
      const auto& bucket = nodes_[node].items;
      SJS_DASSERT(!bucket.empty());
      const u32 pos = rng->UniformU32(static_cast<u32>(bucket.size()));
      out->push_back(bucket[pos].handle);
    }
    return true;
  }

  void Report(u32 q, std::vector<u32>* out) const {
    SJS_DASSERT(out != nullptr);
    if (m_ == 0 || q >= m_) return;

    u32 idx = q + p_;
    while (idx > 0) {
      const auto& bucket = nodes_[idx].items;
      for (const auto& it : bucket) out->push_back(it.handle);
      idx >>= 1;
    }
  }

 private:
  struct Node {
    std::vector<Item> items;
  };

  usize PlacementIndex(u32 handle, u32 backref) const {
    return static_cast<usize>(handle) * static_cast<usize>(max_refs_) + static_cast<usize>(backref);
  }

  void AddToNode(u32 handle, u32 node) {
    SJS_DASSERT(node < nodes_.size());

    u32& sz = placement_size_[handle];
    SJS_DASSERT(sz < max_refs_);

    auto& bucket = nodes_[node].items;
    const u32 pos = static_cast<u32>(bucket.size());

    placements_[PlacementIndex(handle, sz)] = Placement{node, pos};
    bucket.push_back(Item{handle, sz});
    ++sz;
  }

  void RemoveFromNode(u32 node, u32 pos) {
    SJS_DASSERT(node < nodes_.size());
    auto& bucket = nodes_[node].items;
    SJS_DASSERT(!bucket.empty());
    SJS_DASSERT(pos < bucket.size());

    const usize last_pos = bucket.size() - 1;
    if (static_cast<usize>(pos) != last_pos) {
      const Item swapped = bucket[last_pos];
      bucket[pos] = swapped;
      placements_[PlacementIndex(swapped.handle, swapped.backref)].pos = pos;
    }
    bucket.pop_back();
  }

  u32 n_handles_{0};
  u32 m_{0};
  u32 p_{0};
  u32 max_refs_{0};

  std::vector<Node> nodes_;
  std::vector<Placement> placements_;
  std::vector<u32> placement_size_;
};

// -----------------------------------------------------------------------------
// ActiveIndex2D: per-side wrapper for both patterns
// -----------------------------------------------------------------------------

class ActiveIndex2D {
 public:
  void Init(u32 num_handles, u32 num_ylo_ranks) {
    n_ = num_handles;
    m_ = num_ylo_ranks;
    stab_.Init(num_handles, num_ylo_ranks);
    pts_.Init(num_handles, num_ylo_ranks);
  }

  void Clear() {
    n_ = 0;
    m_ = 0;
    stab_.Clear();
    pts_.Clear();
  }

  void ResetActive() {
    stab_.ResetActive();
    pts_.ResetActive();
  }

  void Insert(u32 handle, u32 ylo_rank, u32 yhi_lb_rank) {
    SJS_DASSERT(handle < n_);
    SJS_DASSERT(ylo_rank < m_);
    // yhi_lb_rank is in [0, m_].
    stab_.Insert(handle, ylo_rank, yhi_lb_rank);
    pts_.Insert(handle, ylo_rank);
  }

  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_);
    stab_.Erase(handle);
    pts_.Erase(handle);
  }

  // Pattern A: y_r contains y_q_lo.
  u64 CountA(u32 y_q_lo_rank) const { return stab_.Count(y_q_lo_rank); }

  // Pattern B: y_r_lo in (y_q_lo, y_q_hi).
  u64 CountB(u32 y_q_lo_rank, u32 y_q_hi_lb_rank) const {
    if (m_ == 0) return 0ULL;
    const u32 l = (y_q_lo_rank + 1 <= m_) ? (y_q_lo_rank + 1) : m_;
    const u32 r = std::min(y_q_hi_lb_rank, m_);
    return pts_.CountRange(l, r);
  }

  bool SampleA(u32 y_q_lo_rank, u32 k, Rng* rng, std::vector<u32>* out) const {
    return stab_.Sample(y_q_lo_rank, k, rng, out);
  }

  bool SampleB(u32 y_q_lo_rank, u32 y_q_hi_lb_rank, u32 k, Rng* rng, std::vector<u32>* out) const {
    const u32 l = (y_q_lo_rank + 1 <= m_) ? (y_q_lo_rank + 1) : m_;
    const u32 r = std::min(y_q_hi_lb_rank, m_);
    return pts_.SampleRange(l, r, k, rng, out);
  }

  void ReportA(u32 y_q_lo_rank, std::vector<u32>* out) const { stab_.Report(y_q_lo_rank, out); }

  void ReportB(u32 y_q_lo_rank, u32 y_q_hi_lb_rank, std::vector<u32>* out) const {
    const u32 l = (y_q_lo_rank + 1 <= m_) ? (y_q_lo_rank + 1) : m_;
    const u32 r = std::min(y_q_hi_lb_rank, m_);
    pts_.ReportRange(l, r, out);
  }

 private:
  u32 n_{0};
  u32 m_{0};
  StabbingSegTree stab_;
  RangePointSegTree pts_;
};

// -----------------------------------------------------------------------------
// Shared 2D preprocessing context (events + ranks + active indices)
// -----------------------------------------------------------------------------

template <int Dim, class T>
class Ours2DContext {
 public:
  using DatasetT = Dataset<Dim, T>;

  void Reset() {
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
  }

  bool Build(const DatasetT& ds, PhaseRecorder* phases, std::string* err) {
    if constexpr (Dim != 2) {
      if (err) *err = "Ours2DContext: only Dim=2 is implemented";
      return false;
    }

    ds_ = &ds;

    const usize nR = ds.R.Size();
    const usize nS = ds.S.Size();
    if (nR > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        nS > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "Ours2DContext: relation size exceeds u32";
      return false;
    }

    // 1) Events
    {
      auto scoped = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = join::BuildSweepEvents<Dim, T>(ds, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    // 2) START id mapping (dense 0..|E|-1)
    start_id_of_event_.assign(events_.size(), -1);
    start_event_pos_.clear();
    start_event_pos_.reserve(nR + nS);
    for (usize i = 0; i < events_.size(); ++i) {
      if (events_[i].kind == join::EventKind::Start) {
        start_id_of_event_[i] = static_cast<i32>(start_event_pos_.size());
        start_event_pos_.push_back(i);
      }
    }

    // 3) y-domain from all y-lower endpoints
    {
      auto scoped = phases ? phases->Scoped("build_y_domain") : PhaseRecorder::ScopedPhase(nullptr, "");
      y_coords_.clear();
      y_coords_.reserve(nR + nS);
      for (const auto& b : ds.R.boxes) y_coords_.push_back(b.lo.v[1]);
      for (const auto& b : ds.S.boxes) y_coords_.push_back(b.lo.v[1]);
      std::sort(y_coords_.begin(), y_coords_.end());
      y_coords_.erase(std::unique(y_coords_.begin(), y_coords_.end()), y_coords_.end());
    }

    if (y_coords_.empty()) {
      if (err) *err = "Ours2DContext: empty y-domain";
      return false;
    }

    const u32 m = static_cast<u32>(y_coords_.size());

    auto lb_rank = [&](T v) -> u32 {
      const auto it = std::lower_bound(y_coords_.begin(), y_coords_.end(), v);
      return static_cast<u32>(std::distance(y_coords_.begin(), it));
    };

    // 4) Precompute ranks per box
    {
      auto scoped = phases ? phases->Scoped("build_ranks") : PhaseRecorder::ScopedPhase(nullptr, "");

      ylo_rank_r_.resize(nR);
      yhi_lb_rank_r_.resize(nR);
      for (usize i = 0; i < nR; ++i) {
        const auto& b = ds.R.boxes[i];
        const u32 lo = lb_rank(b.lo.v[1]);
        SJS_DASSERT(lo < m && y_coords_[lo] == b.lo.v[1]);
        ylo_rank_r_[i] = lo;
        yhi_lb_rank_r_[i] = lb_rank(b.hi.v[1]);
      }

      ylo_rank_s_.resize(nS);
      yhi_lb_rank_s_.resize(nS);
      for (usize i = 0; i < nS; ++i) {
        const auto& b = ds.S.boxes[i];
        const u32 lo = lb_rank(b.lo.v[1]);
        SJS_DASSERT(lo < m && y_coords_[lo] == b.lo.v[1]);
        ylo_rank_s_[i] = lo;
        yhi_lb_rank_s_[i] = lb_rank(b.hi.v[1]);
      }
    }

    // 5) Build active indices skeleton
    {
      auto scoped = phases ? phases->Scoped("build_active_indices") : PhaseRecorder::ScopedPhase(nullptr, "");
      active_r_.Init(static_cast<u32>(nR), m);
      active_s_.Init(static_cast<u32>(nS), m);
      active_r_.ResetActive();
      active_s_.ResetActive();
    }

    built_ = true;
    return true;
  }

  bool built() const noexcept { return built_; }
  const DatasetT* dataset() const noexcept { return ds_; }

  const std::vector<join::Event>& events() const noexcept { return events_; }
  const std::vector<i32>& start_id_of_event() const noexcept { return start_id_of_event_; }
  usize num_start_events() const noexcept { return start_event_pos_.size(); }

  const std::vector<u32>& ylo_rank_r() const noexcept { return ylo_rank_r_; }
  const std::vector<u32>& yhi_lb_rank_r() const noexcept { return yhi_lb_rank_r_; }
  const std::vector<u32>& ylo_rank_s() const noexcept { return ylo_rank_s_; }
  const std::vector<u32>& yhi_lb_rank_s() const noexcept { return yhi_lb_rank_s_; }

  ActiveIndex2D& active_r() noexcept { return active_r_; }
  ActiveIndex2D& active_s() noexcept { return active_s_; }
  const ActiveIndex2D& active_r() const noexcept { return active_r_; }
  const ActiveIndex2D& active_s() const noexcept { return active_s_; }

  void ResetActive() {
    active_r_.ResetActive();
    active_s_.ResetActive();
  }

 private:
  const DatasetT* ds_{nullptr};
  bool built_{false};

  std::vector<join::Event> events_;
  std::vector<i32> start_id_of_event_;   // per event position (size events_.size())
  std::vector<usize> start_event_pos_;   // positions of START events in events_

  std::vector<T> y_coords_;              // unique y-lower values
  std::vector<u32> ylo_rank_r_;
  std::vector<u32> yhi_lb_rank_r_;
  std::vector<u32> ylo_rank_s_;
  std::vector<u32> yhi_lb_rank_s_;

  ActiveIndex2D active_r_;
  ActiveIndex2D active_s_;
};

// -----------------------------------------------------------------------------
// Phase2 slot plan (shared by Sampling + Adaptive large branch)
// -----------------------------------------------------------------------------

struct SlotPlan2D {
  // Offsets are size E+1; slots arrays are size t.
  std::vector<u32> offset_a;
  std::vector<u32> offset_b;
  std::vector<u32> slots_a;
  std::vector<u32> slots_b;

  void Clear() {
    offset_a.clear();
    offset_b.clear();
    slots_a.clear();
    slots_b.clear();
  }
};

inline bool BuildSlotPlan2D(u32 t,
                           Rng* rng,
                           const std::vector<u64>& w_total,
                           const std::vector<u64>& w_a,
                           const std::vector<u64>& w_b,
                           SlotPlan2D* plan,
                           std::string* err) {
  SJS_DASSERT(rng != nullptr);
  SJS_DASSERT(plan != nullptr);
  plan->Clear();

  const usize E = w_total.size();
  if (E == 0) {
    if (err) *err = "BuildSlotPlan2D: empty event set";
    return false;
  }

  sampling::AliasTable alias;
  if (!alias.BuildFromU64(Span<const u64>(w_total), err)) {
    if (err && err->empty()) *err = "BuildSlotPlan2D: failed to build alias";
    return false;
  }

  // First pass: per-slot assignment + per-event counts.
  std::vector<u32> eid_of_slot;
  std::vector<u8> pat_of_slot;
  eid_of_slot.resize(t);
  pat_of_slot.resize(t);

  std::vector<u32> cnt_a(E, 0U);
  std::vector<u32> cnt_b(E, 0U);

  for (u32 j = 0; j < t; ++j) {
    const u32 eid = static_cast<u32>(alias.Sample(rng));
    const u64 wa = w_a[eid];
    const u64 wb = w_b[eid];
    const u64 w = wa + wb;
    SJS_DASSERT(w == w_total[eid]);

    // Choose pattern conditional on the event.
    u8 pat = 0;  // 0=A, 1=B
    if (wa == 0) {
      pat = 1;
    } else if (wb == 0) {
      pat = 0;
    } else {
      const u64 r = rng->UniformU64(w);
      pat = (r < wa) ? 0 : 1;
    }

    eid_of_slot[j] = eid;
    pat_of_slot[j] = pat;

    if (pat == 0) {
      ++cnt_a[eid];
    } else {
      ++cnt_b[eid];
    }
  }

  // Prefix sums -> offsets.
  plan->offset_a.resize(E + 1);
  plan->offset_b.resize(E + 1);
  plan->offset_a[0] = 0;
  plan->offset_b[0] = 0;
  for (usize e = 0; e < E; ++e) {
    plan->offset_a[e + 1] = plan->offset_a[e] + cnt_a[e];
    plan->offset_b[e + 1] = plan->offset_b[e] + cnt_b[e];
  }

  const u32 total_a = plan->offset_a[E];
  const u32 total_b = plan->offset_b[E];
  SJS_DASSERT(total_a + total_b == t);

  plan->slots_a.resize(total_a);
  plan->slots_b.resize(total_b);

  // Second pass: stable-fill the slot indices into the flat arrays.
  std::vector<u32> cur_a = plan->offset_a;
  std::vector<u32> cur_b = plan->offset_b;

  for (u32 j = 0; j < t; ++j) {
    const u32 eid = eid_of_slot[j];
    const u8 pat = pat_of_slot[j];
    if (pat == 0) {
      plan->slots_a[cur_a[eid]++] = j;
    } else {
      plan->slots_b[cur_b[eid]++] = j;
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
// Wrapper enumerator around join::PlaneSweepJoinStream so it satisfies
// baselines::IJoinEnumerator (used by the experimental framework).
// -----------------------------------------------------------------------------

template <int Dim, class T>
class PlaneSweepEnumeratorWrapper final : public baselines::IJoinEnumerator {
 public:
  PlaneSweepEnumeratorWrapper(const Relation<Dim, T>& R,
                              const Relation<Dim, T>& S,
                              join::PlaneSweepOptions opt)
      : stream_(R, S, opt) {}

  void Reset() override { stream_.Reset(); }
  bool Next(PairId* out) override { return stream_.Next(out); }
  const join::JoinStats& Stats() const noexcept override { return stream_.Stats(); }

 private:
  join::PlaneSweepJoinStream<Dim, T> stream_;
};

}  // namespace detail

// -----------------------------------------------------------------------------
// OursSamplingBaseline (Sampling variant, Dim=2)
// -----------------------------------------------------------------------------

template <int Dim, class T = Scalar>
class OursSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "ours_sampling"; }

  void Reset() override {
    ctx_.Reset();
    built_ = false;
    weights_valid_ = false;
    W_ = 0;
    w_total_.clear();
    w_a_.clear();
    w_b_.clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    Reset();

    if (!ctx_.Build(ds, phases, err)) {
      return false;
    }

    const usize E = ctx_.num_start_events();
    w_total_.assign(E, 0ULL);
    w_a_.assign(E, 0ULL);
    w_b_.assign(E, 0ULL);

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg, Rng* rng, CountResult* out, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)rng;

    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursSamplingBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    std::fill(w_a_.begin(), w_a_.end(), 0ULL);
    std::fill(w_b_.begin(), w_b_.end(), 0ULL);

    ctx_.ResetActive();

    u64 W = 0;

    auto& ar = ctx_.active_r();
    auto& as = ctx_.active_s();

    const auto& events = ctx_.events();
    const auto& sid_of_pos = ctx_.start_id_of_event();

    const auto& ylo_r = ctx_.ylo_rank_r();
    const auto& yhi_r = ctx_.yhi_lb_rank_r();
    const auto& ylo_s = ctx_.ylo_rank_s();
    const auto& yhi_s = ctx_.yhi_lb_rank_s();

    for (usize pos = 0; pos < events.size(); ++pos) {
      const auto& ev = events[pos];
      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        if (ev.side == join::Side::R) {
          ar.Erase(handle);
        } else {
          as.Erase(handle);
        }
        continue;
      }

      // START event.
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

      W += w;

      // Insert q into its side.
      if (q_is_r) {
        ar.Insert(handle, q_ylo, q_yhi);
      } else {
        as.Insert(handle, q_ylo, q_yhi);
      }
    }

    // Sweep ends with empty active sets, but we keep it explicit for safety.
    ctx_.ResetActive();

    W_ = W;
    weights_valid_ = true;

    if (out) *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "OursSamplingBaseline::Sample: null rng/out";
      return false;
    }

    const u32 t = cfg.run.t;

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    if (t == 0) return true;

    // Ensure Phase1 weights exist.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, nullptr, &tmp, phases, err)) return false;
    }

    if (W_ == 0) {
      // Empty join.
      return true;
    }

    // Phase2: slot plan.
    detail::SlotPlan2D plan;
    {
      auto scoped = phases ? phases->Scoped("phase2_plan") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!detail::BuildSlotPlan2D(t, rng, w_total_, w_a_, w_b_, &plan, err)) {
        if (err && err->empty()) *err = "OursSamplingBaseline::Sample: failed to build slot plan";
        return false;
      }
    }

    out->pairs.resize(t);

    // Phase3: second sweep.
    {
      auto scoped = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

      ctx_.ResetActive();

      auto& ar = ctx_.active_r();
      auto& as = ctx_.active_s();

      const auto& events = ctx_.events();
      const auto& sid_of_pos = ctx_.start_id_of_event();

      const auto& ylo_r = ctx_.ylo_rank_r();
      const auto& yhi_r = ctx_.yhi_lb_rank_r();
      const auto& ylo_s = ctx_.ylo_rank_s();
      const auto& yhi_s = ctx_.yhi_lb_rank_s();

      std::vector<u32> sampled;

      const auto* ds = ctx_.dataset();
      SJS_DASSERT(ds != nullptr);

      for (usize pos = 0; pos < events.size(); ++pos) {
        const auto& ev = events[pos];
        const u32 handle = static_cast<u32>(ev.index);

        if (ev.kind == join::EventKind::End) {
          if (ev.side == join::Side::R) {
            ar.Erase(handle);
          } else {
            as.Erase(handle);
          }
          continue;
        }

        const i32 sid_i32 = sid_of_pos[pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 sid = static_cast<u32>(sid_i32);

        const bool q_is_r = (ev.side == join::Side::R);
        const u32 q_ylo = q_is_r ? ylo_r[ev.index] : ylo_s[ev.index];
        const u32 q_yhi = q_is_r ? yhi_r[ev.index] : yhi_s[ev.index];

        const detail::ActiveIndex2D& other = q_is_r ? as : ar;

        // Pattern A slots.
        {
          const u32 begin = plan.offset_a[sid];
          const u32 end = plan.offset_a[sid + 1];
          const u32 k = end - begin;
          if (k > 0) {
            sampled.clear();
            const bool ok = other.SampleA(q_ylo, k, rng, &sampled);
            if (!ok || sampled.size() != k) {
              if (err) *err = "OursSamplingBaseline::Sample: SampleA failed (inconsistent weights)";
              return false;
            }
            for (u32 i = 0; i < k; ++i) {
              const u32 slot = plan.slots_a[begin + i];
              const u32 oh = sampled[i];
              if (q_is_r) {
                out->pairs[slot] = PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)};
              } else {
                out->pairs[slot] = PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)};
              }
            }
          }
        }

        // Pattern B slots.
        {
          const u32 begin = plan.offset_b[sid];
          const u32 end = plan.offset_b[sid + 1];
          const u32 k = end - begin;
          if (k > 0) {
            sampled.clear();
            const bool ok = other.SampleB(q_ylo, q_yhi, k, rng, &sampled);
            if (!ok || sampled.size() != k) {
              if (err) *err = "OursSamplingBaseline::Sample: SampleB failed (inconsistent weights)";
              return false;
            }
            for (u32 i = 0; i < k; ++i) {
              const u32 slot = plan.slots_b[begin + i];
              const u32 oh = sampled[i];
              if (q_is_r) {
                out->pairs[slot] = PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)};
              } else {
                out->pairs[slot] = PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)};
              }
            }
          }
        }

        // Insert q into its side.
        if (q_is_r) {
          ar.Insert(handle, q_ylo, q_yhi);
        } else {
          as.Insert(handle, q_ylo, q_yhi);
        }
      }

      ctx_.ResetActive();
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursSamplingBaseline::Enumerate: call Build() first";
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
  bool built_{false};

  detail::Ours2DContext<Dim, T> ctx_;

  std::vector<u64> w_total_;
  std::vector<u64> w_a_;
  std::vector<u64> w_b_;

  u64 W_{0};
  bool weights_valid_{false};
};

}  // namespace ours
}  // namespace baselines
}  // namespace sjs
