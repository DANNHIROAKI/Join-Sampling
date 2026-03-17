#pragma once
// baselines/ours/sampling.h
//
// "Our Method" — exact i.i.d. sampling baseline aligned with 研究大纲 SJS-nD.
//
// Supported dimensions in this repository build: Dim = 2,3,4,5.
//
// High-level structure:
//   - Dim=2 keeps the original optimized 2D fast path.
//   - Dim>=3 uses the fixed-mode recursive structures D_{m,g} from Chapter 2.
//
// Sampling variant (Framework II in Chapter 3):
//   Phase 1: sweep once, compute each START-event block weight w_e and its
//            per-mode decomposition. Summing over events gives |J|.
//   Phase 2: sample output slots by event weight, then by mode weight inside
//            the chosen event; this materializes a slot plan grouped by
//            (event, mode).
//   Phase 3: sweep again and fill each (event, mode) group by exact uniform
//            sampling from the corresponding local structure, preserving the
//            global (R,S) PairId orientation.
//
// Correctness-critical conventions:
//   - Half-open boxes [lo, hi) on every axis.
//   - Sweep events on axis 0 with END before START at equal coordinates.
//   - Public PairId values use Relation::GetId(), which must be unique within
//     each relation (enforced by dataset validation).

#include "baselines/baseline_api.h"
#include "baselines/ours/detail/context_nd.h"
#include "baselines/ours/detail/local_mode_count_index.h"
#include "baselines/ours/detail/local_mode_index.h"
#include "baselines/ours/detail/report_enumerator_nd.h"
#include "baselines/ours/detail/slot_plan_nd.h"

#include "core/assert.h"
#include "join/join_enumerator.h"  // optional (kept for debugging / comparison)
#include "join/sweep_events.h"
#include "sampling/alias_table.h"

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
// Key design goal (per SJS v3):
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
  // Fanout for Stage3 k-ary tree.
  // k=4 is a practical sweet spot (depth halves vs binary while keeping cover sizes reasonable).
  static constexpr u32 kFanout = 4;

  void Init(u32 num_handles, u32 num_ranks) {
    n_handles_ = num_handles;
    m_ = num_ranks;

    // Leaves: smallest power of kFanout >= m_.
    p_ = 1;
    depth_ = 0;
    while (p_ < m_) {
      p_ *= kFanout;
      ++depth_;
    }

    // Level layout: level 0 = leaves, level depth_ = root.
    level_offs_.assign(static_cast<usize>(depth_ + 1U), 0U);
    usize total = 0;
    u32 cnt = p_;
    for (u32 lvl = 0; lvl <= depth_; ++lvl) {
      level_offs_[lvl] = static_cast<u32>(total);
      total += static_cast<usize>(cnt);
      cnt /= kFanout;
      if (cnt == 0) break;
    }

    nodes_.assign(total, Node{});

    max_refs_ = depth_ + 1U;  // one ancestor per level

    pos_in_node_.assign(static_cast<usize>(n_handles_) * static_cast<usize>(max_refs_), 0U);
    rank_of_handle_.assign(static_cast<usize>(n_handles_), kInvalidRank);

    touched_flag_.assign(total, 0);
    touched_nodes_.clear();

    active_count_ = 0;
  }

  void Clear() {
    n_handles_ = 0;
    m_ = 0;
    p_ = 0;
    depth_ = 0;
    max_refs_ = 0;
    level_offs_.clear();
    nodes_.clear();
    pos_in_node_.clear();
    rank_of_handle_.clear();
    touched_flag_.clear();
    touched_nodes_.clear();
    active_count_ = 0;
  }

  // Keep the skeleton but drop all active points.
  void ResetActive() {
    if (active_count_ == 0) return;

    for (const u32 node : touched_nodes_) {
      nodes_[node].items.clear();
      touched_flag_[node] = 0;
    }
    touched_nodes_.clear();

    std::fill(rank_of_handle_.begin(), rank_of_handle_.end(), kInvalidRank);
    active_count_ = 0;
  }

  // Optional performance hint: reserve bucket capacities for hot nodes.
  // Reserve the top 'levels' levels near the root.
  void ReserveHotBuckets(u32 max_active, u32 levels = 4) {
    if (max_active == 0 || nodes_.empty()) return;
    if (depth_ == 0) return;

    for (u32 d = 0; d < levels; ++d) {
      if (d > depth_) break;
      const u32 lvl = depth_ - d;
      const u32 first = level_offs_[lvl];
      const u32 last = (lvl == depth_) ? static_cast<u32>(nodes_.size()) : level_offs_[lvl + 1U];
      const u32 nodes_at = last - first;
      const u32 cap_per = std::max<u32>(1U, max_active / std::max<u32>(1U, nodes_at));
      for (u32 node = first; node < last; ++node) {
        auto& b = nodes_[node].items;
        if (b.capacity() < cap_per) b.reserve(cap_per);
      }
    }
  }

  void Insert(u32 handle, u32 rank) {
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(rank < m_);
    SJS_DASSERT(rank_of_handle_[handle] == kInvalidRank);

    rank_of_handle_[handle] = rank;
    ++active_count_;

    const usize base = static_cast<usize>(handle) * static_cast<usize>(max_refs_);

    u32 idx = rank;
    for (u32 lvl = 0; lvl <= depth_; ++lvl) {
      const u32 node = NodeId(lvl, idx);
      AddToNode(handle, node, base, lvl);
      idx /= kFanout;
    }
  }

  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_handles_);
    const u32 rank = rank_of_handle_[handle];
    if (rank == kInvalidRank) return;

    const usize base = static_cast<usize>(handle) * static_cast<usize>(max_refs_);

    u32 idx = rank;
    for (u32 lvl = 0; lvl <= depth_; ++lvl) {
      const u32 node = NodeId(lvl, idx);
      const u32 pos = pos_in_node_[base + static_cast<usize>(lvl)];
      RemoveFromNode(node, pos, lvl);
      idx /= kFanout;
    }

    rank_of_handle_[handle] = kInvalidRank;
    SJS_DASSERT(active_count_ > 0);
    --active_count_;
  }

  // Allocation-free exact count of active handles in ranks [l, r).
  u64 CountRange(u32 l, u32 r) const {
    if (r <= l || m_ == 0) return 0ULL;
    l = std::min(l, m_);
    r = std::min(r, m_);
    if (r <= l) return 0ULL;

    std::array<u32, kMaxCoverNodes> cover;
    const u32 cover_sz = DecomposeOrdered(l, r, &cover);

    u64 total = 0;
    for (u32 i = 0; i < cover_sz; ++i) {
      total += static_cast<u64>(nodes_[cover[i]].items.size());
    }
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

    std::array<u32, kMaxCoverNodes> cover;
    const u32 cover_sz = DecomposeOrdered(l, r, &cover);

    std::array<u32, kMaxCoverNodes> nodes;
    std::array<u64, kMaxCoverNodes> weights;
    u32 nz = 0;
    u64 total = 0;

    for (u32 i = 0; i < cover_sz; ++i) {
      const u32 node = cover[i];
      const u64 w = static_cast<u64>(nodes_[node].items.size());
      if (w == 0) continue;
      nodes[nz] = node;
      weights[nz] = w;
      total += w;
      ++nz;
    }

    if (total == 0 || nz == 0) return false;

    static constexpr u32 kAliasMinK = 64;
    if (nz > 1 && k >= kAliasMinK) {
      sampling::AliasSmall<kMaxCoverNodes> alias;
      (void)alias.BuildFromU64(Span<const u64>(weights.data(), nz), /*err=*/nullptr);

      for (u32 i = 0; i < k; ++i) {
        const u32 bi = static_cast<u32>(alias.Sample(rng));
        const u32 node = nodes[bi];
        const auto& bucket = nodes_[node].items;
        const u32 pos = rng->UniformU32(static_cast<u32>(bucket.size()));
        out->push_back(bucket[pos]);
      }
      return true;
    }

    for (u32 i = 0; i < k; ++i) {
      const u64 x = rng->UniformU64(total);
      u64 cum = 0;
      u32 bi = 0;
      for (; bi < nz; ++bi) {
        cum += weights[bi];
        if (x < cum) break;
      }
      if (bi >= nz) bi = nz - 1;

      const u32 node = nodes[bi];
      const auto& bucket = nodes_[node].items;
      const u32 pos = rng->UniformU32(static_cast<u32>(bucket.size()));
      out->push_back(bucket[pos]);
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

    std::array<u32, kMaxCoverNodes> cover;
    const u32 cover_sz = DecomposeOrdered(l, r, &cover);

    for (u32 i = 0; i < cover_sz; ++i) {
      const auto& bucket = nodes_[cover[i]].items;
      for (const u32 h : bucket) out->push_back(h);
    }
  }

 private:
  struct Node {
    std::vector<u32> items;
  };

  static constexpr u32 kInvalidRank = std::numeric_limits<u32>::max();
  static constexpr usize kMaxCoverNodes = 256;

  u32 NodeId(u32 lvl, u32 idx) const {
    SJS_DASSERT(lvl <= depth_);
    const u32 off = level_offs_[lvl];
    return off + idx;
  }

  void Touch(u32 node) {
    if (touched_flag_[node] == 0) {
      touched_flag_[node] = 1;
      touched_nodes_.push_back(node);
    }
  }

  void AddToNode(u32 handle, u32 node, usize base, u32 lvl) {
    auto& bucket = nodes_[node].items;
    if (bucket.empty()) Touch(node);

    const u32 pos = static_cast<u32>(bucket.size());
    pos_in_node_[base + static_cast<usize>(lvl)] = pos;
    bucket.push_back(handle);
  }

  void RemoveFromNode(u32 node, u32 pos, u32 lvl) {
    auto& bucket = nodes_[node].items;
    const usize last_pos = bucket.size() - 1;
    if (static_cast<usize>(pos) != last_pos) {
      const u32 swapped_handle = bucket[last_pos];
      bucket[pos] = swapped_handle;
      pos_in_node_[static_cast<usize>(swapped_handle) * static_cast<usize>(max_refs_) + static_cast<usize>(lvl)] = pos;
    }
    bucket.pop_back();
  }

  // Canonical cover decomposition of [l,r) (0<=l<=r<=m_).
  // Returns disjoint nodes in left-to-right order.
  u32 DecomposeOrdered(u32 l, u32 r, std::array<u32, kMaxCoverNodes>* out) const {
    SJS_DASSERT(out != nullptr);
    if (r <= l) return 0;

    std::array<u32, kMaxCoverNodes> left;
    std::array<u32, kMaxCoverNodes> right;
    u32 lsz = 0;
    u32 rsz = 0;

    u32 lvl = 0;
    while (l < r) {
      SJS_DASSERT(lvl <= depth_);

      while (l < r && (l % kFanout) != 0U) {
        SJS_DASSERT(lsz < kMaxCoverNodes);
        left[lsz++] = NodeId(lvl, l);
        ++l;
      }
      while (l < r && (r % kFanout) != 0U) {
        --r;
        SJS_DASSERT(rsz < kMaxCoverNodes);
        right[rsz++] = NodeId(lvl, r);
      }

      l /= kFanout;
      r /= kFanout;
      ++lvl;
    }

    const u32 out_sz = lsz + rsz;
    SJS_DASSERT(out_sz <= kMaxCoverNodes);

    u32 out_i = 0;
    for (u32 i = 0; i < lsz; ++i) (*out)[out_i++] = left[i];
    for (u32 i = 0; i < rsz; ++i) (*out)[out_i++] = right[rsz - 1U - i];
    return out_i;
  }

  u32 n_handles_{0};
  u32 m_{0};
  u32 p_{0};
  u32 depth_{0};
  u32 max_refs_{0};

  std::vector<u32> level_offs_;  // per-level base offset into nodes_
  std::vector<Node> nodes_;
  std::vector<u32> pos_in_node_;    // handle*max_refs_ + level -> position
  std::vector<u32> rank_of_handle_; // per-handle rank (kInvalidRank => inactive)

  std::vector<u8> touched_flag_;
  std::vector<u32> touched_nodes_;

  u32 active_count_{0};
};



class StabbingSegTree {
 public:
  static constexpr u32 kFanout = 4;

  void Init(u32 num_handles, u32 num_points) {
    n_handles_ = num_handles;
    m_ = num_points;

    // Leaves: smallest power of kFanout >= m_.
    p_ = 1;
    depth_ = 0;
    while (p_ < m_) {
      p_ *= kFanout;
      ++depth_;
    }

    level_offs_.assign(static_cast<usize>(depth_ + 1U), 0U);
    usize total = 0;
    u32 cnt = p_;
    for (u32 lvl = 0; lvl <= depth_; ++lvl) {
      level_offs_[lvl] = static_cast<u32>(total);
      total += static_cast<usize>(cnt);
      cnt /= kFanout;
      if (cnt == 0) break;
    }

    nodes_.assign(total, Node{});

    touched_flag_.assign(total, 0);
    touched_nodes_.clear();

    refs_.clear();
    ref_off_.assign(static_cast<usize>(n_handles_), 0U);
    ref_len_.assign(static_cast<usize>(n_handles_), 0U);

    active_count_ = 0;
  }

  void Clear() {
    n_handles_ = 0;
    m_ = 0;
    p_ = 0;
    depth_ = 0;
    level_offs_.clear();
    nodes_.clear();
    touched_flag_.clear();
    touched_nodes_.clear();
    refs_.clear();
    ref_off_.clear();
    ref_len_.clear();
    active_count_ = 0;
  }

  // Keep the skeleton but drop all active intervals.
  void ResetActive() {
    if (active_count_ == 0) return;

    for (const u32 node : touched_nodes_) {
      nodes_[node].items.clear();
      touched_flag_[node] = 0;
    }
    touched_nodes_.clear();

    refs_.clear();
    std::fill(ref_len_.begin(), ref_len_.end(), 0U);
    std::fill(ref_off_.begin(), ref_off_.end(), 0U);

    active_count_ = 0;
  }

  void ReserveHotBuckets(u32 max_active, u32 levels = 4) {
    if (max_active == 0 || nodes_.empty()) return;
    if (depth_ == 0) return;

    for (u32 d = 0; d < levels; ++d) {
      if (d > depth_) break;
      const u32 lvl = depth_ - d;
      const u32 first = level_offs_[lvl];
      const u32 last = (lvl == depth_) ? static_cast<u32>(nodes_.size()) : level_offs_[lvl + 1U];
      const u32 nodes_at = last - first;
      const u32 cap_per = std::max<u32>(1U, max_active / std::max<u32>(1U, nodes_at));
      for (u32 node = first; node < last; ++node) {
        auto& b = nodes_[node].items;
        if (b.capacity() < cap_per) b.reserve(cap_per);
      }
    }
  }

  // Insert interval [L, R) on ranks.
  void Insert(u32 handle, u32 L, u32 R) {
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(ref_len_[handle] == 0U);

    if (R <= L || m_ == 0) return;
    L = std::min(L, m_);
    R = std::min(R, m_);
    if (R <= L) return;

    std::array<u32, kMaxCoverNodes> cover;
    const u32 cover_sz = DecomposeOrdered(L, R, &cover);

    const u32 off = static_cast<u32>(refs_.size());

    for (u32 i = 0; i < cover_sz; ++i) {
      const u32 node = cover[i];
      auto& bucket = nodes_[node].items;
      if (bucket.empty()) Touch(node);

      const u32 pos = static_cast<u32>(bucket.size());
      const u32 ref_idx = static_cast<u32>(refs_.size());
      bucket.push_back(Item{handle, ref_idx});
      refs_.push_back(Ref{node, pos});
    }

    ref_off_[handle] = off;
    ref_len_[handle] = cover_sz;
    ++active_count_;
  }

  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_handles_);
    const u32 len = ref_len_[handle];
    if (len == 0U) return;
    const u32 off = ref_off_[handle];

    for (u32 i = 0; i < len; ++i) {
      const u32 ref_idx = off + i;
      const Ref ref = refs_[ref_idx];
      auto& bucket = nodes_[ref.node].items;

      const u32 pos = ref.pos;
      const u32 last_pos = static_cast<u32>(bucket.size() - 1U);
      if (pos != last_pos) {
        const Item swapped = bucket[last_pos];
        bucket[pos] = swapped;
        refs_[swapped.ref_idx].pos = pos;
      }
      bucket.pop_back();
    }

    ref_len_[handle] = 0U;
    ref_off_[handle] = 0U;
    SJS_DASSERT(active_count_ > 0);
    --active_count_;
  }

  u64 Count(u32 q) const {
    if (m_ == 0 || q >= m_) return 0ULL;
    u64 total = 0;

    u32 idx = q;
    for (u32 lvl = 0; lvl <= depth_; ++lvl) {
      const u32 node = NodeId(lvl, idx);
      total += static_cast<u64>(nodes_[node].items.size());
      idx /= kFanout;
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

    std::array<u32, kMaxPathNodes> path_nodes;
    std::array<u64, kMaxPathNodes> weights;
    u32 n = 0;
    u64 total = 0;

    u32 idx = q;
    for (u32 lvl = 0; lvl <= depth_; ++lvl) {
      const u32 node = NodeId(lvl, idx);
      const u64 w = static_cast<u64>(nodes_[node].items.size());
      if (w > 0) {
        path_nodes[n] = node;
        weights[n] = w;
        total += w;
        ++n;
      }
      idx /= kFanout;
    }

    if (total == 0 || n == 0) return false;

    static constexpr u32 kAliasMinK = 64;
    if (n > 1 && k >= kAliasMinK) {
      sampling::AliasSmall<kMaxPathNodes> alias;
      (void)alias.BuildFromU64(Span<const u64>(weights.data(), n), /*err=*/nullptr);

      for (u32 i = 0; i < k; ++i) {
        const u32 bi = static_cast<u32>(alias.Sample(rng));
        const u32 node = path_nodes[bi];
        const auto& bucket = nodes_[node].items;
        const u32 pos = rng->UniformU32(static_cast<u32>(bucket.size()));
        out->push_back(bucket[pos].handle);
      }
      return true;
    }

    for (u32 i = 0; i < k; ++i) {
      const u64 x = rng->UniformU64(total);
      u64 cum = 0;
      u32 bi = 0;
      for (; bi < n; ++bi) {
        cum += weights[bi];
        if (x < cum) break;
      }
      if (bi >= n) bi = n - 1;

      const u32 node = path_nodes[bi];
      const auto& bucket = nodes_[node].items;
      const u32 pos = rng->UniformU32(static_cast<u32>(bucket.size()));
      out->push_back(bucket[pos].handle);
    }
    return true;
  }

  void Report(u32 q, std::vector<u32>* out) const {
    SJS_DASSERT(out != nullptr);
    if (m_ == 0 || q >= m_) return;

    u32 idx = q;
    for (u32 lvl = 0; lvl <= depth_; ++lvl) {
      const u32 node = NodeId(lvl, idx);
      const auto& bucket = nodes_[node].items;
      for (const auto& it : bucket) out->push_back(it.handle);
      idx /= kFanout;
    }
  }

 private:
  struct Item {
    u32 handle;
    u32 ref_idx;
  };

  struct Ref {
    u32 node;
    u32 pos;
  };

  struct Node {
    std::vector<Item> items;
  };

  static constexpr usize kMaxCoverNodes = 256;
  static constexpr usize kMaxPathNodes = 64;

  u32 NodeId(u32 lvl, u32 idx) const {
    SJS_DASSERT(lvl <= depth_);
    return level_offs_[lvl] + idx;
  }

  void Touch(u32 node) {
    if (touched_flag_[node] == 0) {
      touched_flag_[node] = 1;
      touched_nodes_.push_back(node);
    }
  }

  u32 DecomposeOrdered(u32 l, u32 r, std::array<u32, kMaxCoverNodes>* out) const {
    SJS_DASSERT(out != nullptr);
    if (r <= l) return 0;

    std::array<u32, kMaxCoverNodes> left;
    std::array<u32, kMaxCoverNodes> right;
    u32 lsz = 0;
    u32 rsz = 0;

    u32 lvl = 0;
    while (l < r) {
      SJS_DASSERT(lvl <= depth_);

      while (l < r && (l % kFanout) != 0U) {
        SJS_DASSERT(lsz < kMaxCoverNodes);
        left[lsz++] = NodeId(lvl, l);
        ++l;
      }
      while (l < r && (r % kFanout) != 0U) {
        --r;
        SJS_DASSERT(rsz < kMaxCoverNodes);
        right[rsz++] = NodeId(lvl, r);
      }

      l /= kFanout;
      r /= kFanout;
      ++lvl;
    }

    const u32 out_sz = lsz + rsz;
    SJS_DASSERT(out_sz <= kMaxCoverNodes);

    u32 out_i = 0;
    for (u32 i = 0; i < lsz; ++i) (*out)[out_i++] = left[i];
    for (u32 i = 0; i < rsz; ++i) (*out)[out_i++] = right[rsz - 1U - i];
    return out_i;
  }

  u32 n_handles_{0};
  u32 m_{0};
  u32 p_{0};
  u32 depth_{0};

  std::vector<u32> level_offs_;
  std::vector<Node> nodes_;

  std::vector<u8> touched_flag_;
  std::vector<u32> touched_nodes_;

  std::vector<Ref> refs_;
  std::vector<u32> ref_off_;
  std::vector<u32> ref_len_;

  u32 active_count_{0};
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

  void ReserveHotBuckets(u32 max_active) {
    // Both trees can have very hot buckets near the root/top levels.
    // Reserving a small number of top levels reduces reallocations in the sweep.
    pts_.ReserveHotBuckets(max_active);
    stab_.ReserveHotBuckets(max_active);
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
// CountIndex2D: count-only active-set index (two Fenwick trees)
//
// Used only in Phase 1 (Count) to compute exact wa/wb without maintaining the
// heavier sampling/report structures.
//
//   - Pattern A count: stabbing count via difference Fenwick
//       Insert interval [ylo, yhi_lb):  diff[ylo]+=1, diff[yhi_lb]-=1
//       CountA(q_ylo) = prefix_sum(diff, q_ylo)
//
//   - Pattern B count: point range count via point Fenwick
//       Insert point at ylo: pts[ylo]+=1
//       CountB(q_ylo, q_yhi_lb) = sum(pts, (q_ylo+1 .. q_yhi_lb))
// -----------------------------------------------------------------------------

class Fenwick {
 public:
  void Init(u32 n) {
    n_ = n;
    // BIT values represent counts of active rectangles.
    // In our benchmark scale (n up to a few million), i32 is sufficient and
    // reduces memory bandwidth vs i64.
    bit_.assign(static_cast<usize>(n_ + 1), 0);
  }

  void Clear() {
    n_ = 0;
    bit_.clear();
  }

  void Reset() { std::fill(bit_.begin(), bit_.end(), 0); }

  // Add delta at index i (0-based). Requires i < n_.
  void Add(u32 i, i32 delta) {
    SJS_DASSERT(i < n_);
    // BIT uses 1-based internally.
    for (u32 x = i + 1; x <= n_; x += x & (~x + 1)) {
      bit_[x] += delta;
    }
  }

  // Prefix sum over [0..i] inclusive. If i >= n_, returns SumAll().
  i64 PrefixSumI64(u32 i) const {
    if (n_ == 0) return 0;
    if (i >= n_) i = n_ - 1;
    i64 res = 0;
    for (u32 x = i + 1; x > 0; x &= (x - 1)) {
      res += bit_[x];
    }
    return res;
  }

  i64 RangeSumI64(u32 l, u32 r) const {
    if (r <= l || n_ == 0) return 0;
    if (l >= n_) return 0;
    if (r > n_) r = n_;
    if (r <= l) return 0;
    const i64 a = PrefixSumI64(r - 1);
    const i64 b = (l == 0) ? 0 : PrefixSumI64(l - 1);
    return a - b;
  }

 private:
  u32 n_{0};
  std::vector<i32> bit_;  // index 0 unused
};

class CountIndex2D {
 public:
  void Init(u32 num_ylo_ranks) {
    m_ = num_ylo_ranks;
    // diff has length m_+1 (we update at yhi_lb which can be == m_).
    diff_.Init(m_ + 1);
    pts_.Init(m_);
    ResetActive();
  }

  void Clear() {
    m_ = 0;
    diff_.Clear();
    pts_.Clear();
  }

  void ResetActive() {
    diff_.Reset();
    pts_.Reset();
  }

  // Insert rectangle with (ylo_rank, yhi_lb_rank) into active set.
  void Insert(u32 ylo_rank, u32 yhi_lb_rank) {
    if (m_ == 0) return;
    SJS_DASSERT(ylo_rank < m_);
    SJS_DASSERT(yhi_lb_rank <= m_);

    pts_.Add(ylo_rank, +1);
    diff_.Add(ylo_rank, +1);
    diff_.Add(yhi_lb_rank, -1);
  }

  void Erase(u32 ylo_rank, u32 yhi_lb_rank) {
    if (m_ == 0) return;
    SJS_DASSERT(ylo_rank < m_);
    SJS_DASSERT(yhi_lb_rank <= m_);

    pts_.Add(ylo_rank, -1);
    diff_.Add(ylo_rank, -1);
    diff_.Add(yhi_lb_rank, +1);
  }

  // Pattern A: y_r contains y_q_lo.
  u64 CountA(u32 y_q_lo_rank) const {
    if (m_ == 0 || y_q_lo_rank >= m_) return 0ULL;
    const i64 v = diff_.PrefixSumI64(y_q_lo_rank);
    return (v <= 0) ? 0ULL : static_cast<u64>(v);
  }

  // Pattern B: y_r_lo in (y_q_lo, y_q_hi).
  u64 CountB(u32 y_q_lo_rank, u32 y_q_hi_lb_rank) const {
    if (m_ == 0) return 0ULL;
    u32 l = y_q_lo_rank + 1;
    if (l > m_) l = m_;
    const u32 r = std::min(y_q_hi_lb_rank, m_);
    const i64 v = pts_.RangeSumI64(l, r);
    return (v <= 0) ? 0ULL : static_cast<u64>(v);
  }

 private:
  u32 m_{0};
  Fenwick diff_;  // size m_+1
  Fenwick pts_;   // size m_
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

    // Event SoA.
    start_id_of_event_.clear();
    ev_kind_side_.clear();
    ev_handle_.clear();
    ev_ylo_rank_.clear();
    ev_yhi_lb_rank_.clear();
    num_start_events_ = 0;

    // y-domain.
    y_coords_.clear();
    m_ = 0;

    // Handle -> id mapping.
    id_of_r_handle_.clear();
    id_of_s_handle_.clear();

    active_r_.Clear();
    active_s_.Clear();
  }

  bool Build(const DatasetT& ds, PhaseRecorder* phases, std::string* err) {
    if constexpr (Dim != 2) {
      if (err) *err = "Ours2DContext: only Dim=2 is implemented";
      return false;
    }

    Reset();
    ds_ = &ds;

    const usize nR64 = ds.R.Size();
    const usize nS64 = ds.S.Size();
    if (nR64 > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        nS64 > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "Ours2DContext: relation size exceeds u32";
      return false;
    }

    const u32 nR = static_cast<u32>(nR64);
    const u32 nS = static_cast<u32>(nS64);

    n_r_ = nR;
    n_s_ = nS;

    // 1) Build plane-sweep events (AoS from join::BuildSweepEvents).
    std::vector<join::Event> events;
    {
      auto scoped = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events = join::BuildSweepEvents<Dim, T>(ds, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    const usize num_events = events.size();

    // 2) START id mapping + handle renumbering + kind/side cache.
    //
    // Handle renumbering (Stage 2): assign dense handles 0..n-1 in START order
    // within each side. This improves locality of handle-indexed arrays inside
    // the active indices during phase3.
    std::vector<u32> handle_of_r_index;
    std::vector<u32> handle_of_s_index;
    {
      auto scoped = phases ? phases->Scoped("build_event_maps") : PhaseRecorder::ScopedPhase(nullptr, "");

      static constexpr u32 kInvalidHandle = std::numeric_limits<u32>::max();

      handle_of_r_index.assign(nR, kInvalidHandle);
      handle_of_s_index.assign(nS, kInvalidHandle);

      id_of_r_handle_.resize(nR);
      id_of_s_handle_.resize(nS);

      start_id_of_event_.assign(num_events, -1);
      ev_kind_side_.resize(num_events);
      num_start_events_ = 0;

      u32 next_r = 0;
      u32 next_s = 0;

      for (usize pos = 0; pos < num_events; ++pos) {
        const auto& ev = events[pos];
        ev_kind_side_[pos] = (static_cast<u8>(ev.kind) << 1) | static_cast<u8>(ev.side);

        if (ev.kind != join::EventKind::Start) continue;

        // Dense START id in [0, num_start_events_).
        start_id_of_event_[pos] = static_cast<i32>(num_start_events_++);

        if (ev.side == join::Side::R) {
          const u32 idx = static_cast<u32>(ev.index);
          SJS_DASSERT(idx < nR);
          SJS_DASSERT(handle_of_r_index[idx] == kInvalidHandle);
          const u32 h = next_r++;
          handle_of_r_index[idx] = h;
          id_of_r_handle_[h] = ds.R.GetId(idx);
        } else {
          const u32 idx = static_cast<u32>(ev.index);
          SJS_DASSERT(idx < nS);
          SJS_DASSERT(handle_of_s_index[idx] == kInvalidHandle);
          const u32 h = next_s++;
          handle_of_s_index[idx] = h;
          id_of_s_handle_[h] = ds.S.GetId(idx);
        }
      }

      if (next_r != nR || next_s != nS) {
        if (err) *err = "Ours2DContext: handle renumbering failed (missing START events?)";
        return false;
      }
    }

    // 3) y-domain (unique y-lower endpoints) + ylo-ranks via sort+merge.
    std::vector<u32> ylo_rank_r;
    std::vector<u32> yhi_lb_rank_r;
    std::vector<u32> ylo_rank_s;
    std::vector<u32> yhi_lb_rank_s;

    {
      auto scoped = phases ? phases->Scoped("build_y_domain") : PhaseRecorder::ScopedPhase(nullptr, "");

      y_coords_.clear();
      y_coords_.reserve(static_cast<usize>(nR) + static_cast<usize>(nS));

      ylo_rank_r.resize(nR);
      ylo_rank_s.resize(nS);

      struct YEntry {
        T v;
        u32 g;  // global index: 0..nR-1 => R, nR..nR+nS-1 => S
      };

      std::vector<YEntry> ylo;
      ylo.reserve(static_cast<usize>(nR) + static_cast<usize>(nS));

      for (u32 i = 0; i < nR; ++i) ylo.push_back(YEntry{ds.R.boxes[i].lo.v[1], i});
      for (u32 i = 0; i < nS; ++i) ylo.push_back(YEntry{ds.S.boxes[i].lo.v[1], static_cast<u32>(nR + i)});

      std::sort(ylo.begin(), ylo.end(),
                [](const YEntry& a, const YEntry& b) { return a.v < b.v; });

      if (!ylo.empty()) {
        u32 rank = 0;
        T prev = ylo[0].v;
        y_coords_.push_back(prev);

        for (const auto& e : ylo) {
          if (e.v != prev) {
            prev = e.v;
            y_coords_.push_back(prev);
            ++rank;
          }
          if (e.g < nR) {
            ylo_rank_r[e.g] = rank;
          } else {
            ylo_rank_s[e.g - nR] = rank;
          }
        }
      }
    }

    if (y_coords_.empty()) {
      if (err) *err = "Ours2DContext: empty y-domain";
      return false;
    }
    m_ = static_cast<u32>(y_coords_.size());

    // 4) yhi lower-bound ranks via sort+merge against y_coords_.
    {
      auto scoped = phases ? phases->Scoped("build_ranks") : PhaseRecorder::ScopedPhase(nullptr, "");

      yhi_lb_rank_r.resize(nR);
      yhi_lb_rank_s.resize(nS);

      struct YEntry {
        T v;
        u32 g;  // global index
      };

      std::vector<YEntry> yhi;
      yhi.reserve(static_cast<usize>(nR) + static_cast<usize>(nS));

      for (u32 i = 0; i < nR; ++i) yhi.push_back(YEntry{ds.R.boxes[i].hi.v[1], i});
      for (u32 i = 0; i < nS; ++i) yhi.push_back(YEntry{ds.S.boxes[i].hi.v[1], static_cast<u32>(nR + i)});

      std::sort(yhi.begin(), yhi.end(),
                [](const YEntry& a, const YEntry& b) { return a.v < b.v; });

      u32 r = 0;
      for (const auto& e : yhi) {
        while (r < m_ && y_coords_[static_cast<usize>(r)] < e.v) ++r;
        const u32 lb = r;  // in [0, m_]
        if (e.g < nR) {
          yhi_lb_rank_r[e.g] = lb;
        } else {
          yhi_lb_rank_s[e.g - nR] = lb;
        }
      }
    }

    // 5) Cache per-event endpoints + handle (SoA) for cache-friendly sweeps.
    {
      auto scoped = phases ? phases->Scoped("build_event_cache") : PhaseRecorder::ScopedPhase(nullptr, "");

      ev_handle_.resize(num_events);
      ev_ylo_rank_.resize(num_events);
      ev_yhi_lb_rank_.resize(num_events);

      for (usize pos = 0; pos < num_events; ++pos) {
        const auto& ev = events[pos];

        const u32 idx = static_cast<u32>(ev.index);
        if (ev.side == join::Side::R) {
          SJS_DASSERT(idx < nR);
          const u32 h = handle_of_r_index[idx];
          SJS_DASSERT(h < nR);
          ev_handle_[pos] = h;
          ev_ylo_rank_[pos] = ylo_rank_r[idx];
          ev_yhi_lb_rank_[pos] = yhi_lb_rank_r[idx];
        } else {
          SJS_DASSERT(idx < nS);
          const u32 h = handle_of_s_index[idx];
          SJS_DASSERT(h < nS);
          ev_handle_[pos] = h;
          ev_ylo_rank_[pos] = ylo_rank_s[idx];
          ev_yhi_lb_rank_[pos] = yhi_lb_rank_s[idx];
        }
      }
    }

    // Free large build-time temporaries before allocating the active indices.
    std::vector<join::Event>().swap(events);
    std::vector<u32>().swap(handle_of_r_index);
    std::vector<u32>().swap(handle_of_s_index);
    std::vector<u32>().swap(ylo_rank_r);
    std::vector<u32>().swap(ylo_rank_s);
    std::vector<u32>().swap(yhi_lb_rank_r);
    std::vector<u32>().swap(yhi_lb_rank_s);

    // 6) Build active indices skeleton.
    {
      auto scoped = phases ? phases->Scoped("build_active_indices") : PhaseRecorder::ScopedPhase(nullptr, "");
      active_r_.Init(nR, m_);
      active_s_.Init(nS, m_);
      active_r_.ResetActive();
      active_s_.ResetActive();
    }

    built_ = true;
    return true;
  }

  bool built() const noexcept { return built_; }
  const DatasetT* dataset() const noexcept { return ds_; }

  // Event SoA
  usize num_events() const noexcept { return ev_kind_side_.size(); }
  usize num_start_events() const noexcept { return static_cast<usize>(num_start_events_); }
  u32 num_ylo_ranks() const noexcept { return m_; }

  const std::vector<i32>& start_id_of_event() const noexcept { return start_id_of_event_; }
  const std::vector<u8>& ev_kind_side() const noexcept { return ev_kind_side_; }
  const std::vector<u32>& ev_handle() const noexcept { return ev_handle_; }
  const std::vector<u32>& ev_ylo_rank() const noexcept { return ev_ylo_rank_; }
  const std::vector<u32>& ev_yhi_lb_rank() const noexcept { return ev_yhi_lb_rank_; }

  // Handle -> Id (stable external ids from Dataset)
  Id IdR(u32 h) const noexcept {
    SJS_DASSERT(h < id_of_r_handle_.size());
    return id_of_r_handle_[static_cast<usize>(h)];
  }
  Id IdS(u32 h) const noexcept {
    SJS_DASSERT(h < id_of_s_handle_.size());
    return id_of_s_handle_[static_cast<usize>(h)];
  }

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

  u32 n_r_{0};
  u32 n_s_{0};

  // Event SoA (length = 2*(nR+nS)).
  std::vector<i32> start_id_of_event_;  // -1 for END, else dense START id
  std::vector<u8> ev_kind_side_;        // (kind<<1)|side, kind:0=end,1=start; side:0=R,1=S
  std::vector<u32> ev_handle_;          // per event: side-local handle (renumbered)
  std::vector<u32> ev_ylo_rank_;        // per event: ylo rank in [0,m_)
  std::vector<u32> ev_yhi_lb_rank_;     // per event: lower_bound(yhi) rank in [0,m_]

  u32 num_start_events_{0};

  // y-domain
  std::vector<T> y_coords_;
  u32 m_{0};

  // Handle -> external id.
  std::vector<Id> id_of_r_handle_;
  std::vector<Id> id_of_s_handle_;

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
                           const std::vector<u64>& w_a,
                           const std::vector<u64>& w_b,
                           SlotPlan2D* plan,
                           std::string* err) {
  SJS_DASSERT(rng != nullptr);
  SJS_DASSERT(plan != nullptr);
  plan->Clear();

  const usize E = w_a.size();
  if (w_b.size() != E) {
    if (err) *err = "BuildSlotPlan2D: w_a/w_b size mismatch";
    return false;
  }
  if (E == 0) {
    if (err) *err = "BuildSlotPlan2D: empty event set";
    return false;
  }

  // Pack per-slot (eid,pattern) into one u32: info = (eid<<1) | pat.
  // This requires eid to fit in 31 bits.
  if (E > (static_cast<usize>(std::numeric_limits<u32>::max()) >> 1)) {
    if (err) *err = "BuildSlotPlan2D: too many events to pack slot_info";
    return false;
  }

  // Build an alias table over events with weight w_e = w_a[e] + w_b[e].
  // Then, conditioned on the sampled event eid, choose Pattern A with
  // probability w_a[eid]/(w_a[eid] + w_b[eid]).
  std::vector<u64> w_event;
  w_event.resize(E);

  // Per-event threshold for Pattern A:
  //   pat = A  iff  u < thr_a[eid]
  // where u is uniform on [0, 2^32). We use AliasTable's threshold encoding.
  std::vector<u32> thr_a;
  thr_a.resize(E);

  for (usize e = 0; e < E; ++e) {
    const u64 wa = w_a[e];
    const u64 wb = w_b[e];
    const u64 w = wa + wb;
    w_event[e] = w;

    if (w == 0) {
      thr_a[e] = 0;
    } else if (wa == 0) {
      thr_a[e] = 0;
    } else if (wa == w) {
      thr_a[e] = sampling::AliasTable::kProbOne;
    } else {
#if defined(__SIZEOF_INT128__)
      const __uint128_t num = static_cast<__uint128_t>(wa) << 32;  // wa * 2^32
      u64 thr = static_cast<u64>(num / static_cast<__uint128_t>(w));
#else
      const long double x =
          (static_cast<long double>(wa) * sampling::AliasTable::kProbScaleLd) /
          static_cast<long double>(w);
      u64 thr = static_cast<u64>(x);
#endif
      if (thr >= sampling::AliasTable::kProbScaleU64) thr = sampling::AliasTable::kProbScaleU64 - 1;
      thr_a[e] = static_cast<u32>(thr);
    }
  }

  sampling::AliasTable alias;
  if (!alias.BuildFromU64(Span<const u64>(w_event), err)) {
    if (err && err->empty()) *err = "BuildSlotPlan2D: failed to build event alias";
    return false;
  }

  // Release large temporary early (saves peak memory for large E).
  std::vector<u64>().swap(w_event);

  // First pass: per-slot assignment + per-event counts.
  std::vector<u32> slot_info;
  slot_info.resize(t);

  std::vector<u32> cnt_a(E, 0U);
  std::vector<u32> cnt_b(E, 0U);

  for (u32 j = 0; j < t; ++j) {
    const u32 eid = static_cast<u32>(alias.Sample(rng));  // in [0, E)

    const u32 thr = thr_a[static_cast<usize>(eid)];
    const u32 pat = (thr == sampling::AliasTable::kProbOne) ? 0U : ((rng->NextU32() < thr) ? 0U : 1U);

    const u32 info = (eid << 1) | pat;
    slot_info[j] = info;

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
    const u32 info = slot_info[j];
    const u32 eid = info >> 1;
    const u32 pat = info & 1U;
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

// -----------------------------------------------------------------------------
// OursReportJoinEnumerator2D
// -----------------------------------------------------------------------------
// A deterministic join enumerator that follows the same sweep + pattern
// decomposition used by "Our Method".
//
// Why this exists:
//   The project-wide EnumSampling and Adaptive runners rely on baseline->Enumerate()
//   and may perform one or two passes over the join stream. Returning a generic
//   PlaneSweep enumerator here can be catastrophically slow on adversarial cases
//   (e.g., large active sets), even if the baseline's Count/Sample logic is fast.
//
// This enumerator uses the same ActiveIndex2D (stabbing + range-point segtree)
// to *directly* report true intersections, avoiding O(|active|) scans.

template <int Dim, class T>
class OursReportJoinEnumerator2D final : public baselines::IJoinEnumerator {
 public:
  using Ctx = Ours2DContext<Dim, T>;

  explicit OursReportJoinEnumerator2D(Ctx* ctx) : ctx_(ctx) {
    static_assert(Dim == 2, "OursReportJoinEnumerator2D is only implemented for Dim=2");
    SJS_DASSERT(ctx_ != nullptr);
    Reset();
  }

  void Reset() override {
    SJS_DASSERT(ctx_ != nullptr);
    ctx_->ResetActive();

    pos_ = 0;
    stage_ = Stage::Scan;
    tmp_.clear();
    tmp_i_ = 0;

    active_r_ = 0;
    active_s_ = 0;

    stats_.Reset();
    stats_.num_events = static_cast<u64>(ctx_->num_events());
  }

  bool Next(PairId* out) override {
    SJS_DASSERT(out != nullptr);
    SJS_DASSERT(ctx_ != nullptr);

    auto& ar = ctx_->active_r();
    auto& as = ctx_->active_s();

    const auto& ks = ctx_->ev_kind_side();
    const auto& h = ctx_->ev_handle();
    const auto& ev_ylo = ctx_->ev_ylo_rank();
    const auto& ev_yhi = ctx_->ev_yhi_lb_rank();

    while (true) {
      // If we're in the middle of emitting pairs for the current START event.
      if (stage_ == Stage::Emit) {
        if (tmp_i_ < tmp_.size()) {
          const u32 oh = tmp_[tmp_i_++];
          if (q_is_r_) {
            *out = PairId{ctx_->IdR(q_handle_), ctx_->IdS(oh)};
          } else {
            *out = PairId{ctx_->IdR(oh), ctx_->IdS(q_handle_)};
          }
          ++stats_.output_pairs;
          ++stats_.candidate_checks;
          return true;
        }

        // Finished current buffer; advance to the next pattern or finalize START.
        if (pat_ == Pattern::A) {
          // Switch to Pattern B.
          pat_ = Pattern::B;
          tmp_.clear();
          const detail::ActiveIndex2D& other = q_is_r_ ? as : ar;
          other.ReportB(q_ylo_, q_yhi_, &tmp_);
          tmp_i_ = 0;
          continue;
        }

        // Finished Pattern B: insert q and move on.
        if (q_is_r_) {
          ar.Insert(q_handle_, q_ylo_, q_yhi_);
          ++active_r_;
          if (active_r_ > stats_.active_max_r) stats_.active_max_r = active_r_;
        } else {
          as.Insert(q_handle_, q_ylo_, q_yhi_);
          ++active_s_;
          if (active_s_ > stats_.active_max_s) stats_.active_max_s = active_s_;
        }

        stage_ = Stage::Scan;
        ++pos_;
        continue;
      }

      // Scan events until we find the next START (or reach end).
      if (pos_ >= ks.size()) {
        // Leave ctx_ in a clean state for safety.
        ctx_->ResetActive();
        return false;
      }

      const u8 kind_side = ks[pos_];
      const bool is_start = ((kind_side >> 1) != 0);
      const bool is_r = ((kind_side & 1U) == 0);

      const u32 handle = h[pos_];

      if (!is_start) {
        // END event.
        if (is_r) {
          ar.Erase(handle);
          SJS_DASSERT(active_r_ > 0);
          --active_r_;
        } else {
          as.Erase(handle);
          SJS_DASSERT(active_s_ > 0);
          --active_s_;
        }
        ++pos_;
        continue;
      }

      // START event: prepare Pattern A buffer, but do NOT insert q yet.
      q_is_r_ = is_r;
      q_handle_ = handle;
      q_ylo_ = ev_ylo[pos_];
      q_yhi_ = ev_yhi[pos_];

      pat_ = Pattern::A;
      tmp_.clear();
      const detail::ActiveIndex2D& other = q_is_r_ ? as : ar;
      other.ReportA(q_ylo_, &tmp_);
      tmp_i_ = 0;
      stage_ = Stage::Emit;
      // Loop back to emit from tmp_.
    }
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  enum class Stage : u8 {
    Scan = 0,
    Emit = 1,
  };

  enum class Pattern : u8 {
    A = 0,
    B = 1,
  };

  Ctx* ctx_{nullptr};

  // Iteration state.
  usize pos_{0};
  Stage stage_{Stage::Scan};
  Pattern pat_{Pattern::A};

  // Current START event (q).
  bool q_is_r_{true};
  u32 q_handle_{0};
  u32 q_ylo_{0};
  u32 q_yhi_{0};

  // Buffer of opposite-side handles for the current pattern.
  std::vector<u32> tmp_;
  usize tmp_i_{0};

  // Lightweight stats.
  join::JoinStats stats_;
  u64 active_r_{0};
  u64 active_s_{0};
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
  static constexpr u32 kModeCount = detail::kModeCountV<Dim>;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "ours_sampling"; }

  void Reset() override {
    ctx_.Reset();
    nd_ctx_.Reset();
    built_ = false;
    weights_valid_ = false;
    W_ = 0;
    max_active_r_ = 0;
    max_active_s_ = 0;
    w_a_.clear();
    w_b_.clear();
    nd_w_modes_.clear();
    count_r_.Clear();
    count_s_.Clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    Reset();

    if constexpr (Dim == 2) {
      if (!ctx_.Build(ds, phases, err)) {
        return false;
      }

      const usize E = ctx_.num_start_events();
      w_a_.assign(E, 0ULL);
      w_b_.assign(E, 0ULL);

      // Count-only BIT indices for Phase 1.
      const u32 m = ctx_.num_ylo_ranks();
      count_r_.Init(m);
      count_s_.Init(m);

      built_ = true;
      return true;
    } else {
      if (!nd_ctx_.Build(ds, phases, err)) {
        return false;
      }
      nd_w_modes_.assign(nd_ctx_.num_start_events(), detail::ModeWeights<Dim>{});
      built_ = true;
      return true;
    }
  }

  bool Count(const Config& cfg, Rng* rng, CountResult* out, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)rng;

    if constexpr (Dim == 2) {
      if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
        if (err) *err = "OursSamplingBaseline::Count: call Build() first";
        return false;
      }

      auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

      std::fill(w_a_.begin(), w_a_.end(), 0ULL);
      std::fill(w_b_.begin(), w_b_.end(), 0ULL);

      // Count-only active sets (two BITs per side).
      count_r_.ResetActive();
      count_s_.ResetActive();

      u64 W = 0;

      // Track maximum active-set sizes (for optional capacity reservations in phase3).
      u32 active_r_sz = 0;
      u32 active_s_sz = 0;
      u32 max_active_r = 0;
      u32 max_active_s = 0;

      const auto& ks = ctx_.ev_kind_side();
      const auto& sid_of_pos = ctx_.start_id_of_event();
      const auto& ev_ylo = ctx_.ev_ylo_rank();
      const auto& ev_yhi = ctx_.ev_yhi_lb_rank();

      for (usize pos = 0; pos < ks.size(); ++pos) {
        const u8 kind_side = ks[pos];
        const bool is_start = ((kind_side >> 1) != 0);
        const bool is_r = ((kind_side & 1U) == 0);

        if (!is_start) {
          // END event.
          if (is_r) {
            count_r_.Erase(ev_ylo[pos], ev_yhi[pos]);
            SJS_DASSERT(active_r_sz > 0);
            --active_r_sz;
          } else {
            count_s_.Erase(ev_ylo[pos], ev_yhi[pos]);
            SJS_DASSERT(active_s_sz > 0);
            --active_s_sz;
          }
          continue;
        }

        // START event.
        const i32 sid_i32 = sid_of_pos[pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 sid = static_cast<u32>(sid_i32);

        const u32 q_ylo = ev_ylo[pos];
        const u32 q_yhi = ev_yhi[pos];

        const detail::CountIndex2D& other = is_r ? count_s_ : count_r_;

        const u64 wa = other.CountA(q_ylo);
        const u64 wb = other.CountB(q_ylo, q_yhi);
        const u64 w = wa + wb;

        w_a_[sid] = wa;
        w_b_[sid] = wb;
        W += w;

        // Insert q into its side (and update active size tracking).
        if (is_r) {
          count_r_.Insert(q_ylo, q_yhi);
          ++active_r_sz;
          max_active_r = std::max(max_active_r, active_r_sz);
        } else {
          count_s_.Insert(q_ylo, q_yhi);
          ++active_s_sz;
          max_active_s = std::max(max_active_s, active_s_sz);
        }
      }

      max_active_r_ = max_active_r;
      max_active_s_ = max_active_s;

      W_ = W;
      weights_valid_ = true;

      if (out) *out = MakeExactCount(W_);
      return true;
    } else {
      if (!built_ || !nd_ctx_.built() || nd_ctx_.dataset() == nullptr) {
        if (err) *err = "OursSamplingBaseline::Count: call Build() first";
        return false;
      }

      auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");
      for (auto& arr : nd_w_modes_) arr.fill(0ULL);
      nd_ctx_.ResetActive();

      u64 W = 0;
      u32 active_r_sz = 0;
      u32 active_s_sz = 0;
      u32 max_active_r = 0;
      u32 max_active_s = 0;

      const auto& ks = nd_ctx_.ev_kind_side();
      const auto& sid_of_pos = nd_ctx_.start_id_of_event();
      const auto& ev_handle = nd_ctx_.ev_handle();

      for (usize pos = 0; pos < ks.size(); ++pos) {
        const u8 kind_side = ks[pos];
        const bool is_start = ((kind_side >> 1U) != 0U);
        const bool is_r = ((kind_side & 1U) == 0U);
        const u32 handle = ev_handle[pos];

        if (!is_start) {
          if (is_r) {
            nd_ctx_.EraseR(handle);
            SJS_DASSERT(active_r_sz > 0);
            --active_r_sz;
          } else {
            nd_ctx_.EraseS(handle);
            SJS_DASSERT(active_s_sz > 0);
            --active_s_sz;
          }
          continue;
        }

        const i32 sid_i32 = sid_of_pos[pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 sid = static_cast<u32>(sid_i32);
        const auto& q = is_r ? nd_ctx_.BoxR(handle) : nd_ctx_.BoxS(handle);
        const auto& other = is_r ? nd_ctx_.modes_s() : nd_ctx_.modes_r();

        u64 w = 0ULL;
        for (u32 m = 0; m < kModeCount; ++m) {
          const u64 c = other.CountMode(m, q);
          nd_w_modes_[sid][m] = c;
          w += c;
        }
        W += w;

        if (is_r) {
          nd_ctx_.InsertR(handle);
          ++active_r_sz;
          max_active_r = std::max(max_active_r, active_r_sz);
        } else {
          nd_ctx_.InsertS(handle);
          ++active_s_sz;
          max_active_s = std::max(max_active_s, active_s_sz);
        }
      }

      nd_ctx_.ResetActive();
      max_active_r_ = max_active_r;
      max_active_s_ = max_active_s;
      W_ = W;
      weights_valid_ = true;
      if (out) *out = MakeExactCount(W_);
      return true;
    }
  }

  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    if (!rng || !out) {
      if (err) *err = "OursSamplingBaseline::Sample: null rng/out";
      return false;
    }

    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "OursSamplingBaseline::Sample: run.t too large (must fit in u32)";
      return false;
    }
    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    if (t == 0) return true;

    if constexpr (Dim == 2) {
      if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
        if (err) *err = "OursSamplingBaseline::Sample: call Build() first";
        return false;
      }

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
        if (!detail::BuildSlotPlan2D(t, rng, w_a_, w_b_, &plan, err)) {
          if (err && err->empty()) *err = "OursSamplingBaseline::Sample: failed to build slot plan";
          return false;
        }
      }

      out->pairs.resize(t);

      // Phase3: second sweep (sampling).
      {
        auto scoped = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

        ctx_.ResetActive();

        // Optional: reserve capacities for hot buckets using max active sizes observed in phase1.
        ctx_.active_r().ReserveHotBuckets(max_active_r_);
        ctx_.active_s().ReserveHotBuckets(max_active_s_);

        auto& ar = ctx_.active_r();
        auto& as = ctx_.active_s();

        const auto& ks = ctx_.ev_kind_side();
        const auto& sid_of_pos = ctx_.start_id_of_event();
        const auto& ev_handle = ctx_.ev_handle();
        const auto& ev_ylo = ctx_.ev_ylo_rank();
        const auto& ev_yhi = ctx_.ev_yhi_lb_rank();

        std::vector<u32> sampled;

        for (usize pos = 0; pos < ks.size(); ++pos) {
          const u8 kind_side = ks[pos];
          const bool is_start = ((kind_side >> 1) != 0);
          const bool is_r = ((kind_side & 1U) == 0);

          const u32 q_handle = ev_handle[pos];
          const u32 q_ylo = ev_ylo[pos];
          const u32 q_yhi = ev_yhi[pos];

          if (!is_start) {
            // END event.
            if (is_r) {
              ar.Erase(q_handle);
            } else {
              as.Erase(q_handle);
            }
            continue;
          }

          const i32 sid_i32 = sid_of_pos[pos];
          SJS_DASSERT(sid_i32 >= 0);
          const u32 sid = static_cast<u32>(sid_i32);

          const detail::ActiveIndex2D& other = is_r ? as : ar;

          // Cache q's external id once per START (micro-optimization).
          const Id q_id = is_r ? ctx_.IdR(q_handle) : ctx_.IdS(q_handle);

          // Pattern A slots.
          {
            const u32 begin = plan.offset_a[sid];
            const u32 end = plan.offset_a[sid + 1];
            const u32 k = end - begin;
            if (k > 0) {
              const bool ok = other.SampleA(q_ylo, k, rng, &sampled);
              if (!ok || sampled.size() != k) {
                if (err) *err = "OursSamplingBaseline::Sample: SampleA failed (inconsistent weights)";
                return false;
              }
              for (u32 i = 0; i < k; ++i) {
                const u32 slot = plan.slots_a[begin + i];
                const u32 oh = sampled[i];
                if (is_r) {
                  out->pairs[slot] = PairId{q_id, ctx_.IdS(oh)};
                } else {
                  out->pairs[slot] = PairId{ctx_.IdR(oh), q_id};
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
              const bool ok = other.SampleB(q_ylo, q_yhi, k, rng, &sampled);
              if (!ok || sampled.size() != k) {
                if (err) *err = "OursSamplingBaseline::Sample: SampleB failed (inconsistent weights)";
                return false;
              }
              for (u32 i = 0; i < k; ++i) {
                const u32 slot = plan.slots_b[begin + i];
                const u32 oh = sampled[i];
                if (is_r) {
                  out->pairs[slot] = PairId{q_id, ctx_.IdS(oh)};
                } else {
                  out->pairs[slot] = PairId{ctx_.IdR(oh), q_id};
                }
              }
            }
          }

          // Insert q into its side.
          if (is_r) {
            ar.Insert(q_handle, q_ylo, q_yhi);
          } else {
            as.Insert(q_handle, q_ylo, q_yhi);
          }
        }
      }

      return true;
    } else {
      if (!built_ || !nd_ctx_.built() || nd_ctx_.dataset() == nullptr) {
        if (err) *err = "OursSamplingBaseline::Sample: call Build() first";
        return false;
      }

      if (!weights_valid_) {
        CountResult tmp;
        if (!Count(cfg, nullptr, &tmp, phases, err)) return false;
      }
      if (W_ == 0) return true;

      detail::SlotPlanND<Dim> plan;
      {
        auto scoped = phases ? phases->Scoped("phase2_plan") : PhaseRecorder::ScopedPhase(nullptr, "");
        if (!detail::BuildSlotPlanND<Dim>(t, rng, nd_w_modes_, &plan, err)) {
          if (err && err->empty()) *err = "OursSamplingBaseline::Sample: failed to build slot plan";
          return false;
        }
      }

      out->pairs.resize(t);

      {
        auto scoped = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");
        nd_ctx_.ResetActive();
        std::vector<u32> sampled;

        const auto& ks = nd_ctx_.ev_kind_side();
        const auto& sid_of_pos = nd_ctx_.start_id_of_event();
        const auto& ev_handle = nd_ctx_.ev_handle();

        for (usize pos = 0; pos < ks.size(); ++pos) {
          const u8 kind_side = ks[pos];
          const bool is_start = ((kind_side >> 1U) != 0U);
          const bool is_r = ((kind_side & 1U) == 0U);
          const u32 q_handle = ev_handle[pos];

          if (!is_start) {
            if (is_r) nd_ctx_.EraseR(q_handle);
            else nd_ctx_.EraseS(q_handle);
            continue;
          }

          const i32 sid_i32 = sid_of_pos[pos];
          SJS_DASSERT(sid_i32 >= 0);
          const u32 sid = static_cast<u32>(sid_i32);
          bool need_any = false;
          for (u32 m = 0; m < kModeCount; ++m) {
            if (plan.offsets_by_mode[m][static_cast<usize>(sid)] !=
                plan.offsets_by_mode[m][static_cast<usize>(sid + 1)]) {
              need_any = true;
              break;
            }
          }

          const auto& q = is_r ? nd_ctx_.BoxR(q_handle) : nd_ctx_.BoxS(q_handle);
          const auto& other = is_r ? nd_ctx_.modes_s() : nd_ctx_.modes_r();
          const Id q_id = is_r ? nd_ctx_.IdR(q_handle) : nd_ctx_.IdS(q_handle);

          if (need_any) {
            for (u32 m = 0; m < kModeCount; ++m) {
              const u32 begin = plan.offsets_by_mode[m][static_cast<usize>(sid)];
              const u32 end = plan.offsets_by_mode[m][static_cast<usize>(sid + 1)];
              const u32 k = end - begin;
              if (k == 0) continue;
              if (!other.SampleMode(m, q, k, rng, &sampled) || sampled.size() != k) {
                if (err) *err = "OursSamplingBaseline::Sample: fixed-mode sampling failed (inconsistent weights)";
                return false;
              }
              for (u32 i = 0; i < k; ++i) {
                const u32 slot = plan.slots_by_mode[m][static_cast<usize>(begin + i)];
                const u32 oh = sampled[i];
                if (is_r) {
                  out->pairs[slot] = PairId{q_id, nd_ctx_.IdS(oh)};
                } else {
                  out->pairs[slot] = PairId{nd_ctx_.IdR(oh), q_id};
                }
              }
            }
          }

          if (is_r) nd_ctx_.InsertR(q_handle);
          else nd_ctx_.InsertS(q_handle);
        }
        nd_ctx_.ResetActive();
      }
      return true;
    }
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)phases;
    if constexpr (Dim == 2) {
      if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
        if (err) *err = "OursSamplingBaseline::Enumerate: call Build() first";
        return nullptr;
      }
      return std::make_unique<detail::OursReportJoinEnumerator2D<Dim, T>>(&ctx_);
    } else {
      if (!built_ || !nd_ctx_.built() || nd_ctx_.dataset() == nullptr) {
        if (err) *err = "OursSamplingBaseline::Enumerate: call Build() first";
        return nullptr;
      }
      return std::make_unique<detail::OursReportJoinEnumeratorND<Dim, T>>(&nd_ctx_);
    }
  }

 private:
  bool built_{false};

  detail::Ours2DContext<Dim, T> ctx_;
  detail::OursNDContext<Dim, T> nd_ctx_;

  std::vector<u64> w_a_;
  std::vector<u64> w_b_;
  std::vector<detail::ModeWeights<Dim>> nd_w_modes_;

  detail::CountIndex2D count_r_;
  detail::CountIndex2D count_s_;

  u32 max_active_r_{0};
  u32 max_active_s_{0};

  u64 W_{0};
  bool weights_valid_{false};
};


extern template class OursSamplingBaseline<2, sjs::Scalar>;
extern template class OursSamplingBaseline<3, sjs::Scalar>;
extern template class OursSamplingBaseline<4, sjs::Scalar>;
extern template class OursSamplingBaseline<5, sjs::Scalar>;

}  // namespace ours
}  // namespace baselines
}  // namespace sjs
