#pragma once
// baselines/ours/detail/local_mode_index.h
//
// Recursive fixed-mode local structures for the high-dimensional SJS extension.
//
// For Dim >= 3, after sweeping axis 0 we maintain one fixed-mode structure
// \mathcal D_{Dim-1, g} per mode g in {A,B}^{Dim-1}. Each structure is built on
// the full static universe of one relation, and recursively routes objects and
// queries through local static universes U(v) exactly as in the research outline.

#include "baselines/ours/detail/mode_mask.h"
#include "core/assert.h"
#include "core/rng.h"
#include "core/types.h"
#include "io/dataset.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace ours {
namespace detail {

namespace nd_detail {

inline u32 SampleWeightedIndex(const std::vector<u64>& weights, u64 total, Rng* rng) {
  SJS_DASSERT(rng != nullptr);
  SJS_DASSERT(total > 0);
  const u64 x = rng->UniformU64(total);
  u64 cum = 0;
  for (u32 i = 0; i < static_cast<u32>(weights.size()); ++i) {
    cum += weights[static_cast<usize>(i)];
    if (x < cum) return i;
  }
  return static_cast<u32>(weights.size() - 1);
}

inline void BuildGroupedSlots(const std::vector<u64>& weights,
                              u32 k,
                              Rng* rng,
                              std::vector<u32>* offsets,
                              std::vector<u32>* slots,
                              std::vector<u32>* counts,
                              std::string* err) {
  SJS_DASSERT(offsets != nullptr);
  SJS_DASSERT(slots != nullptr);
  SJS_DASSERT(counts != nullptr);
  const usize n = weights.size();
  offsets->assign(n + 1U, 0U);
  slots->assign(static_cast<usize>(k), 0U);
  counts->assign(n, 0U);
  if (k == 0 || n == 0) return;

  u64 total = 0;
  for (u64 w : weights) total += w;
  if (total == 0) {
    if (err) *err = "BuildGroupedSlots: zero total weight";
    return;
  }

  std::vector<u32> slot_node(static_cast<usize>(k), 0U);
  for (u32 slot = 0; slot < k; ++slot) {
    const u32 idx = SampleWeightedIndex(weights, total, rng);
    slot_node[static_cast<usize>(slot)] = idx;
    ++(*counts)[static_cast<usize>(idx)];
  }

  for (usize i = 0; i < n; ++i) {
    (*offsets)[i + 1U] = (*offsets)[i] + (*counts)[i];
  }
  std::vector<u32> cur = *offsets;
  for (u32 slot = 0; slot < k; ++slot) {
    const u32 idx = slot_node[static_cast<usize>(slot)];
    const u32 pos = cur[static_cast<usize>(idx)]++;
    (*slots)[static_cast<usize>(pos)] = slot;
  }
}

}  // namespace nd_detail


template <int Dim, class T>
class FixedModeLocalIndexND {
 public:
  using RelationT = Relation<Dim, T>;
  using BoxT = Box<Dim, T>;

  FixedModeLocalIndexND() = default;

  void Clear() {
    rel_ = nullptr;
    axis_offset_ = 1;
    local_dim_ = 0;
    mode_mask_ = 0;
    head_is_a_ = true;
    root_ = kInvalidNode;
    handles_.clear();
    active_.clear();
    assign_nodes_.clear();
    bucket_pos_.clear();
    nodes_.clear();
    coords_.clear();
    rank_lo_.clear();
    local_rank_.clear();
  }

  bool Build(const RelationT& rel,
             const std::vector<u32>& universe_handles,
             int axis_offset,
             u32 mode_mask,
             std::string* err) {
    if (err) err->clear();
    Clear();
    rel_ = &rel;
    axis_offset_ = axis_offset;
    local_dim_ = Dim - axis_offset;
    if (local_dim_ <= 0) {
      if (err) *err = "FixedModeLocalIndexND::Build: local_dim must be >= 1";
      return false;
    }

    const u32 mode_cap = (local_dim_ >= 31) ? std::numeric_limits<u32>::max()
                                            : ((1u << static_cast<u32>(local_dim_)) - 1u);
    mode_mask_ = mode_mask & mode_cap;
    head_is_a_ = ((mode_mask_ & 1u) == 0u);

    handles_ = universe_handles;
    std::sort(handles_.begin(), handles_.end());
    handles_.erase(std::unique(handles_.begin(), handles_.end()), handles_.end());
    active_.assign(handles_.size(), 0);
    assign_nodes_.assign(handles_.size(), {});
    bucket_pos_.assign(handles_.size(), {});
    local_rank_.assign(handles_.size(), 0U);

    if (handles_.empty()) {
      root_ = kInvalidNode;
      return true;
    }

    if (head_is_a_) {
      BuildA(err);
    } else {
      BuildB(err);
    }
    if (err && !err->empty()) return false;
    return true;
  }

  bool Empty() const noexcept { return handles_.empty(); }

  void ResetActive() {
    for (usize i = 0; i < handles_.size(); ++i) {
      if (active_[i] != 0) Erase(handles_[i]);
    }
  }

  void Insert(u32 handle) {
    if (!rel_ || root_ == kInvalidNode) return;
    const u32 local = LocalIndexOf(handle);
    SJS_DASSERT(local < handles_.size());
    if (active_[local] != 0) return;
    active_[local] = 1;

    const auto& nodes = assign_nodes_[local];
    if (local_dim_ == 1) {
      auto& posv = bucket_pos_[local];
      posv.resize(nodes.size());
      for (usize slot = 0; slot < nodes.size(); ++slot) {
        const u32 node_id = nodes[slot];
        auto& bucket = nodes_[node_id].bucket;
        posv[slot] = static_cast<u32>(bucket.size());
        bucket.push_back(BucketItem{handle, static_cast<u32>(slot)});
      }
    } else {
      for (u32 node_id : nodes) {
        auto* child = nodes_[node_id].child.get();
        SJS_DASSERT(child != nullptr);
        child->Insert(handle);
      }
    }
  }

  void Erase(u32 handle) {
    if (!rel_ || root_ == kInvalidNode) return;
    const u32 local = LocalIndexOf(handle);
    SJS_DASSERT(local < handles_.size());
    if (active_[local] == 0) return;

    const auto& nodes = assign_nodes_[local];
    if (local_dim_ == 1) {
      auto& posv = bucket_pos_[local];
      SJS_DASSERT(posv.size() == nodes.size());
      for (usize slot = 0; slot < nodes.size(); ++slot) {
        const u32 node_id = nodes[slot];
        auto& bucket = nodes_[node_id].bucket;
        const u32 pos = posv[slot];
        SJS_DASSERT(pos < bucket.size());
        const usize last_pos = bucket.size() - 1U;
        if (static_cast<usize>(pos) != last_pos) {
          const BucketItem swapped = bucket[last_pos];
          bucket[static_cast<usize>(pos)] = swapped;
          const u32 swapped_local = LocalIndexOf(swapped.handle);
          bucket_pos_[swapped_local][swapped.slot] = pos;
        }
        bucket.pop_back();
      }
      posv.clear();
    } else {
      for (u32 node_id : nodes) {
        auto* child = nodes_[node_id].child.get();
        SJS_DASSERT(child != nullptr);
        child->Erase(handle);
      }
    }

    active_[local] = 0;
  }

  u64 Count(const BoxT& q) const {
    if (!rel_ || root_ == kInvalidNode) return 0ULL;
    std::vector<u32> qnodes;
    CollectQueryNodes(q, &qnodes);
    if (qnodes.empty()) return 0ULL;

    u64 total = 0ULL;
    if (local_dim_ == 1) {
      for (u32 node_id : qnodes) {
        total += static_cast<u64>(nodes_[node_id].bucket.size());
      }
      return total;
    }

    for (u32 node_id : qnodes) {
      const auto* child = nodes_[node_id].child.get();
      if (!child) continue;
      total += child->Count(q);
    }
    return total;
  }

  void Report(const BoxT& q, std::vector<u32>* out) const {
    SJS_DASSERT(out != nullptr);
    if (!rel_ || root_ == kInvalidNode) return;
    std::vector<u32> qnodes;
    CollectQueryNodes(q, &qnodes);
    if (qnodes.empty()) return;

    if (local_dim_ == 1) {
      for (u32 node_id : qnodes) {
        const auto& bucket = nodes_[node_id].bucket;
        for (const auto& item : bucket) out->push_back(item.handle);
      }
      return;
    }

    for (u32 node_id : qnodes) {
      const auto* child = nodes_[node_id].child.get();
      if (child) child->Report(q, out);
    }
  }

  bool Sample(const BoxT& q, u32 k, Rng* rng, std::vector<u32>* out) const {
    SJS_DASSERT(out != nullptr);
    SJS_DASSERT(rng != nullptr);
    out->clear();
    out->reserve(k);
    if (k == 0) return true;
    if (!rel_ || root_ == kInvalidNode) return false;

    std::vector<u32> qnodes;
    CollectQueryNodes(q, &qnodes);
    if (qnodes.empty()) return false;

    std::vector<u64> weights;
    std::vector<u32> live_nodes;
    weights.reserve(qnodes.size());
    live_nodes.reserve(qnodes.size());

    if (local_dim_ == 1) {
      u64 total = 0ULL;
      for (u32 node_id : qnodes) {
        const u64 w = static_cast<u64>(nodes_[node_id].bucket.size());
        if (w == 0) continue;
        live_nodes.push_back(node_id);
        weights.push_back(w);
        total += w;
      }
      if (total == 0ULL) return false;

      for (u32 i = 0; i < k; ++i) {
        const u32 idx = nd_detail::SampleWeightedIndex(weights, total, rng);
        const auto& bucket = nodes_[live_nodes[idx]].bucket;
        const u32 pos = rng->UniformU32(static_cast<u32>(bucket.size()));
        out->push_back(bucket[static_cast<usize>(pos)].handle);
      }
      return true;
    }

    u64 total = 0ULL;
    for (u32 node_id : qnodes) {
      const auto* child = nodes_[node_id].child.get();
      if (!child) continue;
      const u64 w = child->Count(q);
      if (w == 0) continue;
      live_nodes.push_back(node_id);
      weights.push_back(w);
      total += w;
    }
    if (total == 0ULL) return false;

    std::vector<u32> offsets;
    std::vector<u32> slots;
    std::vector<u32> counts;
    std::string group_err;
    nd_detail::BuildGroupedSlots(weights, k, rng, &offsets, &slots, &counts, &group_err);
    if (!group_err.empty()) return false;

    out->assign(static_cast<usize>(k), 0U);
    std::vector<u32> sampled;
    for (usize i = 0; i < live_nodes.size(); ++i) {
      const u32 cnt = counts[i];
      if (cnt == 0) continue;
      const auto* child = nodes_[live_nodes[i]].child.get();
      SJS_DASSERT(child != nullptr);
      if (!child->Sample(q, cnt, rng, &sampled) || sampled.size() != cnt) return false;
      const u32 begin = offsets[i];
      for (u32 j = 0; j < cnt; ++j) {
        const u32 slot = slots[static_cast<usize>(begin + j)];
        (*out)[static_cast<usize>(slot)] = sampled[static_cast<usize>(j)];
      }
    }
    return true;
  }

 private:
  struct BucketItem {
    u32 handle = 0;
    u32 slot = 0;  // position in assign_nodes_[local]
  };

  struct Node {
    u32 l = 0;
    u32 r = 0;
    i32 left = -1;
    i32 right = -1;

    // Used only for local_dim_ == 1.
    std::vector<BucketItem> bucket;

    // Used only for local_dim_ > 1.
    std::vector<u32> member_handles;
    std::unique_ptr<FixedModeLocalIndexND> child;
  };

  static constexpr u32 kInvalidNode = std::numeric_limits<u32>::max();

  const RelationT* rel_{nullptr};
  int axis_offset_{1};
  int local_dim_{0};
  u32 mode_mask_{0};
  bool head_is_a_{true};
  u32 root_{kInvalidNode};

  std::vector<u32> handles_;  // sorted global handles of the local universe U
  std::vector<u8> active_;
  std::vector<std::vector<u32>> assign_nodes_;
  std::vector<std::vector<u32>> bucket_pos_;  // valid only for local_dim_ == 1

  std::vector<Node> nodes_;

  // A-mode outer data.
  std::vector<T> coords_;  // unique sorted endpoints, size = num_leaves + 1

  // B-mode outer data.
  std::vector<T> rank_lo_;    // lower endpoints in local key order, size = |U|
  std::vector<u32> local_rank_;  // local handle -> rank in [0, |U|)

  u32 LocalIndexOf(u32 handle) const {
    const auto it = std::lower_bound(handles_.begin(), handles_.end(), handle);
    SJS_DASSERT(it != handles_.end() && *it == handle);
    return static_cast<u32>(it - handles_.begin());
  }

  T LoOf(u32 local) const { return rel_->boxes[static_cast<usize>(handles_[local])].lo.v[static_cast<usize>(axis_offset_)]; }
  T HiOf(u32 local) const { return rel_->boxes[static_cast<usize>(handles_[local])].hi.v[static_cast<usize>(axis_offset_)]; }
  Id IdOf(u32 local) const { return rel_->GetId(static_cast<usize>(handles_[local])); }

  u32 BuildTree(u32 l, u32 r) {
    const u32 id = static_cast<u32>(nodes_.size());
    nodes_.push_back(Node{});
    nodes_[id].l = l;
    nodes_[id].r = r;
    if (r - l > 1) {
      const u32 mid = l + (r - l) / 2U;
      nodes_[id].left = static_cast<i32>(BuildTree(l, mid));
      nodes_[id].right = static_cast<i32>(BuildTree(mid, r));
    }
    return id;
  }

  void CollectCover(u32 node_id, u32 ql, u32 qr, std::vector<u32>* out) const {
    const Node& node = nodes_[node_id];
    if (qr <= node.l || node.r <= ql) return;
    if (ql <= node.l && node.r <= qr) {
      out->push_back(node_id);
      return;
    }
    if (node.left >= 0) CollectCover(static_cast<u32>(node.left), ql, qr, out);
    if (node.right >= 0) CollectCover(static_cast<u32>(node.right), ql, qr, out);
  }

  void CollectPath(u32 node_id, u32 leaf, std::vector<u32>* out) const {
    const Node& node = nodes_[node_id];
    SJS_DASSERT(node.l <= leaf && leaf < node.r);
    out->push_back(node_id);
    if (node.r - node.l == 1) return;
    const u32 mid = node.l + (node.r - node.l) / 2U;
    if (leaf < mid) {
      CollectPath(static_cast<u32>(node.left), leaf, out);
    } else {
      CollectPath(static_cast<u32>(node.right), leaf, out);
    }
  }

  void BuildA(std::string* err) {
    coords_.clear();
    coords_.reserve(handles_.size() * 2U);
    for (u32 h : handles_) {
      const auto& b = rel_->boxes[static_cast<usize>(h)];
      coords_.push_back(b.lo.v[static_cast<usize>(axis_offset_)]);
      coords_.push_back(b.hi.v[static_cast<usize>(axis_offset_)]);
    }
    std::sort(coords_.begin(), coords_.end());
    coords_.erase(std::unique(coords_.begin(), coords_.end()), coords_.end());
    if (coords_.size() < 2U) {
      if (err) *err = "FixedModeLocalIndexND::BuildA: coordinate universe too small";
      return;
    }

    nodes_.clear();
    root_ = BuildTree(0U, static_cast<u32>(coords_.size() - 1U));
    bucket_pos_.assign(handles_.size(), {});

    std::vector<u32> cover;
    cover.reserve(64);
    for (u32 local = 0; local < static_cast<u32>(handles_.size()); ++local) {
      const T lo = LoOf(local);
      const T hi = HiOf(local);
      const u32 l = static_cast<u32>(std::lower_bound(coords_.begin(), coords_.end(), lo) - coords_.begin());
      const u32 r = static_cast<u32>(std::lower_bound(coords_.begin(), coords_.end(), hi) - coords_.begin());
      SJS_DASSERT(l < r);
      cover.clear();
      CollectCover(root_, l, r, &cover);
      assign_nodes_[local] = cover;
      if (local_dim_ > 1) {
        const u32 h = handles_[local];
        for (u32 node_id : cover) nodes_[node_id].member_handles.push_back(h);
      }
    }

    if (local_dim_ > 1) BuildChildren(err);
  }

  void BuildB(std::string* err) {
    struct Key {
      T lo{};
      Id id{};
      u32 local{};
    };

    std::vector<Key> keys;
    keys.reserve(handles_.size());
    for (u32 local = 0; local < static_cast<u32>(handles_.size()); ++local) {
      keys.push_back(Key{LoOf(local), IdOf(local), local});
    }
    std::sort(keys.begin(), keys.end(), [](const Key& a, const Key& b) {
      if (a.lo < b.lo) return true;
      if (b.lo < a.lo) return false;
      return a.id < b.id;
    });

    rank_lo_.resize(keys.size());
    for (u32 rank = 0; rank < static_cast<u32>(keys.size()); ++rank) {
      rank_lo_[static_cast<usize>(rank)] = keys[static_cast<usize>(rank)].lo;
      local_rank_[keys[static_cast<usize>(rank)].local] = rank;
    }

    nodes_.clear();
    root_ = BuildTree(0U, static_cast<u32>(keys.size()));
    bucket_pos_.assign(handles_.size(), {});

    std::vector<u32> path;
    path.reserve(64);
    for (u32 local = 0; local < static_cast<u32>(handles_.size()); ++local) {
      path.clear();
      CollectPath(root_, local_rank_[local], &path);
      assign_nodes_[local] = path;
      if (local_dim_ > 1) {
        const u32 h = handles_[local];
        for (u32 node_id : path) nodes_[node_id].member_handles.push_back(h);
      }
    }

    if (local_dim_ > 1) BuildChildren(err);
  }

  void BuildChildren(std::string* err) {
    const u32 tail_mode = (mode_mask_ >> 1U);
    for (auto& node : nodes_) {
      if (node.member_handles.empty()) continue;
      node.child = std::make_unique<FixedModeLocalIndexND>();
      if (!node.child->Build(*rel_, node.member_handles, axis_offset_ + 1, tail_mode, err)) return;
      std::vector<u32>().swap(node.member_handles);
    }
  }

  void CollectQueryNodes(const BoxT& q, std::vector<u32>* out) const {
    SJS_DASSERT(out != nullptr);
    out->clear();
    if (root_ == kInvalidNode) return;
    const T qlo = q.lo.v[static_cast<usize>(axis_offset_)];
    const T qhi = q.hi.v[static_cast<usize>(axis_offset_)];

    if (head_is_a_) {
      if (coords_.empty()) return;
      if (qlo < coords_.front()) return;
      if (!(qlo < coords_.back())) return;
      const auto it = std::upper_bound(coords_.begin(), coords_.end(), qlo);
      const u32 leaf = static_cast<u32>((it - coords_.begin()) - 1);
      CollectPath(root_, leaf, out);
      return;
    }

    if (rank_lo_.empty()) return;
    const u32 l = static_cast<u32>(std::upper_bound(rank_lo_.begin(), rank_lo_.end(), qlo) - rank_lo_.begin());
    const u32 r = static_cast<u32>(std::lower_bound(rank_lo_.begin(), rank_lo_.end(), qhi) - rank_lo_.begin());
    if (l >= r) return;
    CollectCover(root_, l, r, out);
  }
};


template <int Dim, class T>
class ModeFamilyND {
 public:
  using RelationT = Relation<Dim, T>;
  using BoxT = Box<Dim, T>;
  static constexpr u32 kModeCount = kModeCountV<Dim>;

  void Clear() {
    rel_ = nullptr;
    for (auto& idx : indices_) idx.reset();
  }

  bool Build(const RelationT& rel, std::string* err) {
    if (err) err->clear();
    Clear();
    rel_ = &rel;
    std::vector<u32> universe(rel.Size());
    for (u32 i = 0; i < static_cast<u32>(rel.Size()); ++i) universe[static_cast<usize>(i)] = i;
    for (u32 mode = 0; mode < kModeCount; ++mode) {
      indices_[static_cast<usize>(mode)] = std::make_unique<FixedModeLocalIndexND<Dim, T>>();
      if (!indices_[static_cast<usize>(mode)]->Build(rel, universe, /*axis_offset=*/1, mode, err)) {
        return false;
      }
    }
    return true;
  }

  void ResetActive() {
    for (auto& idx : indices_) {
      if (idx) idx->ResetActive();
    }
  }

  void InsertAll(u32 handle) {
    for (auto& idx : indices_) {
      if (idx) idx->Insert(handle);
    }
  }

  void EraseAll(u32 handle) {
    for (auto& idx : indices_) {
      if (idx) idx->Erase(handle);
    }
  }

  u64 CountMode(u32 mode, const BoxT& q) const {
    SJS_DASSERT(mode < kModeCount);
    const auto& idx = indices_[static_cast<usize>(mode)];
    return idx ? idx->Count(q) : 0ULL;
  }

  void ReportMode(u32 mode, const BoxT& q, std::vector<u32>* out) const {
    SJS_DASSERT(mode < kModeCount);
    const auto& idx = indices_[static_cast<usize>(mode)];
    if (idx) idx->Report(q, out);
  }

  bool SampleMode(u32 mode, const BoxT& q, u32 k, Rng* rng, std::vector<u32>* out) const {
    SJS_DASSERT(mode < kModeCount);
    const auto& idx = indices_[static_cast<usize>(mode)];
    return idx ? idx->Sample(q, k, rng, out) : false;
  }

 private:
  const RelationT* rel_{nullptr};
  std::array<std::unique_ptr<FixedModeLocalIndexND<Dim, T>>, kModeCount> indices_{};
};

}  // namespace detail
}  // namespace ours
}  // namespace baselines
}  // namespace sjs
