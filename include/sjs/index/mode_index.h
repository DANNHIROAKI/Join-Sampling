#pragma once
// sjs/index/mode_index.h
//
// §2.4 Fixed-mode multi-dimensional recursive structure D_{m,g}.
//
// We fix scan dimension = axis 0 (as in SJS-HighDims.md).
// Remaining projected dimensions are axes 1..Dim-1, count m = Dim-1.
//
// A mode g in {A,B}^m is encoded as a bitmask over remaining dims:
//   bit=0 -> A (stabbing):   L(s) <= L(q) < R(s)
//   bit=1 -> B (rank-range): L(q) < L(s) < R(q)
//
// D_{m,g} is implemented recursively:
// - Base m=1: use 1D stabbing segtree (A) or rank-range segtree (B) with DenseBucket buckets.
// - Recursive m>1: outer segtree nodes each hold an inner D_{m-1,g_{2..m}} for remaining dims.
//   Outer query produces a disjoint node set, and we aggregate inner results without duplication.
//   SAMPLE uses alias + assignment + grouped inner sampling (Theorem 2.8).

#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/core/types.h"

#include "sjs/geometry/box.h"

#include "sjs/index/dense_bucket.h"
#include "sjs/index/occ_pool.h"
#include "sjs/index/rank_space.h"
#include "sjs/index/segtree_common.h"

#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace sjs {
namespace index {

// Per-projected-axis shared universe (used by all modes and all recursion nodes).
template <class ValueT = Scalar, class HandleT = Id>
struct ModeAxisData {
  using Value = ValueT;
  using Handle = HandleT;
  using RankSpaceT = RankSpace<Value, Handle>;

  // For stabbing (A): sorted unique endpoints X, leaves represent [X[i],X[i+1})
  std::shared_ptr<const std::vector<Value>> coords;
  u32 stab_leaves = 0;
  u32 stab_base = 1;

  // For rank-range (B): rank space on a(u)=L(u)
  std::shared_ptr<const RankSpaceT> rank;
  u32 rank_leaves = 0;
  u32 rank_base = 1;
};

// Runtime fixed-mode index D_{m,g}, where g is encoded by mode_mask (low bits).
// axis_offset tells which projected axis this node corresponds to:
//   axis_offset=0 -> full axis 1
//   axis_offset=1 -> full axis 2
//   ...
template <int Dim, class ValueT = Scalar, class HandleT = Id>
class ModeIndexND {
 public:
  static_assert(Dim >= 2, "ModeIndexND requires Dim >= 2 (otherwise m=0).");

  using Value = ValueT;
  using Handle = HandleT;
  using BoxT = Box<Dim, Value>;
  using AxisData = ModeAxisData<Value, Handle>;
  using Bucket = DenseBucket<Handle>;

  ModeIndexND() = default;

  ModeIndexND(const ModeIndexND&) = delete;
  ModeIndexND& operator=(const ModeIndexND&) = delete;
  ModeIndexND(ModeIndexND&&) noexcept = default;
  ModeIndexND& operator=(ModeIndexND&&) noexcept = default;

  // mode_mask encodes the remaining bits starting at this axis_offset.
  // remaining dims = (Dim-1 - axis_offset).
  void Init(const std::vector<BoxT>* boxes,
            const std::vector<AxisData>* axes,
            OccPool* pool,
            int axis_offset,
            u32 mode_mask) {
    boxes_ = boxes;
    axes_ = axes;
    pool_ = pool;
    axis_offset_ = axis_offset;
    mask_ = mode_mask;

    SJS_ASSERT(boxes_ != nullptr);
    SJS_ASSERT(axes_ != nullptr);
    SJS_ASSERT(pool_ != nullptr);

    const int m_total = Dim - 1;
    SJS_ASSERT(axis_offset_ >= 0 && axis_offset_ < m_total);
    remaining_ = m_total - axis_offset_;

    // Leaf iff remaining_ == 1
    if (remaining_ == 1) {
      store_.template emplace<BucketMap>();
    } else {
      store_.template emplace<ChildMap>();
    }
  }

  bool IsInitialized() const noexcept { return boxes_ != nullptr && axes_ != nullptr && pool_ != nullptr; }
  int AxisOffset() const noexcept { return axis_offset_; }
  int RemainingDims() const noexcept { return remaining_; }
  u32 ModeMask() const noexcept { return mask_; }

  // INSERT(s): insert handle into this fixed-mode structure.
  // Deletion is global via shared OccPool::EraseAll(handle).
  void Insert(Handle h) {
    if (!IsInitialized()) return;
    if (remaining_ <= 0) return;

    std::vector<u32> nodes;
    OuterNodesForInsert_(h, &nodes);
    if (nodes.empty()) return;

    if (IsLeaf_()) {
      for (u32 v : nodes) {
        Bucket* bkt = GetOrCreateBucket_(v);
        bkt->Insert(h);
      }
    } else {
      for (u32 v : nodes) {
        ModeIndexND* child = GetOrCreateChild_(v);
        child->Insert(h);
      }
    }
  }

  // COUNT_g(q)
  u64 Count(const BoxT& q) const {
    if (!IsInitialized()) return 0;
    if (remaining_ <= 0) return 0;

    std::vector<u32> nodes;
    OuterNodesForQuery_(q, &nodes);
    if (nodes.empty()) return 0;

    u64 sum = 0;
    if (IsLeaf_()) {
      for (u32 v : nodes) {
        const Bucket* bkt = FindBucket_(v);
        if (bkt) sum += static_cast<u64>(bkt->Size());
      }
    } else {
      for (u32 v : nodes) {
        const ModeIndexND* child = FindChild_(v);
        if (child) sum += child->Count(q);
      }
    }
    return sum;
  }

  // REPORT_g(q): append handles (no duplicates).
  void Report(const BoxT& q, std::vector<Handle>* out) const {
    SJS_ASSERT(out != nullptr);
    if (!IsInitialized()) return;
    if (remaining_ <= 0) return;

    std::vector<u32> nodes;
    OuterNodesForQuery_(q, &nodes);

    if (IsLeaf_()) {
      for (u32 v : nodes) {
        const Bucket* bkt = FindBucket_(v);
        if (bkt) bkt->Report(out);
      }
    } else {
      for (u32 v : nodes) {
        const ModeIndexND* child = FindChild_(v);
        if (child) child->Report(q, out);
      }
    }
  }

  // SAMPLE_g(q,k): return k i.i.d. uniform samples (with replacement).
  // Implements Theorem 2.8 via alias over disjoint outer nodes + grouped inner sampling.
  void Sample(const BoxT& q, u64 k, Rng* rng, std::vector<Handle>* out) const {
    SJS_ASSERT(rng != nullptr);
    SJS_ASSERT(out != nullptr);
    out->clear();
    if (k == 0) return;
    if (!IsInitialized()) return;
    if (remaining_ <= 0) return;

    std::vector<u32> nodes;
    OuterNodesForQuery_(q, &nodes);
    if (nodes.empty()) return;

    // c_v for each outer node v
    std::vector<u64> c(nodes.size(), 0);
    u64 C = 0;

    if (IsLeaf_()) {
      for (usize i = 0; i < nodes.size(); ++i) {
        const Bucket* bkt = FindBucket_(nodes[i]);
        const u64 sz = bkt ? static_cast<u64>(bkt->Size()) : 0ULL;
        c[i] = sz;
        C += sz;
      }
    } else {
      for (usize i = 0; i < nodes.size(); ++i) {
        const ModeIndexND* child = FindChild_(nodes[i]);
        const u64 sz = child ? child->Count(q) : 0ULL;
        c[i] = sz;
        C += sz;
      }
    }

    if (C == 0) return;

    // Alias over nodes with weights c[i]
    sampling::AliasTable alias;
    (void)alias.BuildFromU64(Span<const u64>(c.data(), c.size()));

    out->assign(static_cast<usize>(k), Handle{});

    // Step 2: assign each sample position to a node index i
    std::vector<std::vector<usize>> pos(nodes.size());
    for (u64 t = 0; t < k; ++t) {
      const usize i = alias.Sample(rng);
      pos[i].push_back(static_cast<usize>(t));
    }

    // Step 3/4: per node i, sample k_i times inside its block and write back.
    if (IsLeaf_()) {
      for (usize i = 0; i < nodes.size(); ++i) {
        const usize ki = pos[i].size();
        if (ki == 0) continue;

        const Bucket* bkt = FindBucket_(nodes[i]);
        // c[i]>0 => bucket must exist and non-empty
        SJS_DASSERT(bkt != nullptr && !bkt->Empty());

        for (usize j = 0; j < ki; ++j) {
          Handle pick{};
          const bool ok = bkt->SampleOne(rng, &pick);
          SJS_DASSERT(ok);
          (*out)[pos[i][j]] = pick;
        }
      }
    } else {
      for (usize i = 0; i < nodes.size(); ++i) {
        const usize ki = pos[i].size();
        if (ki == 0) continue;

        const ModeIndexND* child = FindChild_(nodes[i]);
        SJS_DASSERT(child != nullptr);

        std::vector<Handle> tmp;
        child->Sample(q, static_cast<u64>(ki), rng, &tmp);
        SJS_DASSERT(tmp.size() == ki);

        for (usize j = 0; j < ki; ++j) {
          (*out)[pos[i][j]] = tmp[j];
        }
      }
    }
  }

 private:
  using BucketMap = std::unordered_map<u32, std::unique_ptr<Bucket>>;
  using ChildMap = std::unordered_map<u32, std::unique_ptr<ModeIndexND>>;

  bool IsLeaf_() const noexcept { return remaining_ == 1; }
  bool IsRange_() const noexcept { return (mask_ & 1u) != 0u; }
  u32 ChildMask_() const noexcept { return mask_ >> 1; }

  // Full axis index in Box: scan axis is 0, projected axes start at 1.
  int FullAxis_() const noexcept { return axis_offset_ + 1; }

  const AxisData& CurAxis_() const {
    SJS_ASSERT(axes_ != nullptr);
    return (*axes_)[static_cast<usize>(axis_offset_)];
  }

  // For stabbing endpoints (insertion): map endpoint value to coords index.
  u32 EndpointIndex_(Value v) const {
    const auto& ax = CurAxis_();
    SJS_DASSERT(ax.coords != nullptr);
    const auto& X = *ax.coords;
    auto it = std::lower_bound(X.begin(), X.end(), v);
    SJS_DASSERT_MSG(it != X.end(), "EndpointIndex_: endpoint should be within coords universe");
    const u32 idx = static_cast<u32>(it - X.begin());
    // debug: endpoints from the stored relation should be present exactly
    SJS_DASSERT_MSG(idx < X.size() && (*it == v), "EndpointIndex_: endpoint not found exactly in coords");
    return idx;
  }

  // For stabbing query point: find segment i s.t. X[i] <= x < X[i+1]
  // Return stab_leaves as invalid if out-of-range.
  u32 SegmentIndex_(Value x) const {
    const auto& ax = CurAxis_();
    if (ax.coords == nullptr || ax.coords->size() < 2) return ax.stab_leaves;

    const auto& X = *ax.coords;
    auto it = std::upper_bound(X.begin(), X.end(), x);
    if (it == X.begin()) return ax.stab_leaves;
    const u32 i = static_cast<u32>((it - X.begin()) - 1);
    if (i >= ax.stab_leaves) return ax.stab_leaves;
    return i;
  }

  void OuterNodesForInsert_(Handle h, std::vector<u32>* out) const {
    SJS_ASSERT(out != nullptr);
    out->clear();

    const auto& ax = CurAxis_();
    const int axis = FullAxis_();
    const auto& b = (*boxes_)[static_cast<usize>(h)];

    if (IsRange_()) {
      if (ax.rank_leaves == 0) return;
      SJS_DASSERT(ax.rank != nullptr);
      const u32 r = ax.rank->RankOf(h);
      if (r >= ax.rank_leaves) return;
      PathToRoot(r, ax.rank_base, out);
      return;
    }

    // Stabbing (A): canonical cover of [L,R)
    if (ax.stab_leaves == 0) return;
    SJS_DASSERT(ax.coords != nullptr);

    const Value L = b.lo[axis];
    const Value R = b.hi[axis];
    if (!(L < R)) return;

    const u32 l = EndpointIndex_(L);
    const u32 r = EndpointIndex_(R);
    if (l >= r) return;

    DecomposeCover(l, r, ax.stab_base, out);
  }

  void OuterNodesForQuery_(const BoxT& q, std::vector<u32>* out) const {
    SJS_ASSERT(out != nullptr);
    out->clear();

    const auto& ax = CurAxis_();
    const int axis = FullAxis_();

    if (IsRange_()) {
      // (ell, r) -> [L,R) -> cover nodes
      if (ax.rank_leaves == 0) return;
      SJS_DASSERT(ax.rank != nullptr);

      const Value ell = q.lo[axis];
      const Value rr = q.hi[axis];
      if (!(ell < rr)) return;

      const auto lr = ax.rank->OpenToRankRange(ell, rr);
      const u32 L = lr.first;
      const u32 R = lr.second;
      if (L >= R) return;

      DecomposeCoverOrdered(L, R, ax.rank_base, out);
      return;
    }

    // Stabbing (A): query path for point x = L(q)
    if (ax.stab_leaves == 0) return;
    const Value x = q.lo[axis];
    const u32 seg = SegmentIndex_(x);
    if (seg >= ax.stab_leaves) return;

    PathFromRoot(seg, ax.stab_base, out);
  }

  Bucket* GetOrCreateBucket_(u32 node) {
    BucketMap& m = std::get<BucketMap>(store_);
    auto& ptr = m[node];
    if (!ptr) ptr = std::make_unique<Bucket>(pool_);
    return ptr.get();
  }

  const Bucket* FindBucket_(u32 node) const {
    const BucketMap& m = std::get<BucketMap>(store_);
    auto it = m.find(node);
    if (it == m.end()) return nullptr;
    return it->second.get();
  }

  ModeIndexND* GetOrCreateChild_(u32 node) {
    ChildMap& m = std::get<ChildMap>(store_);
    auto& ptr = m[node];
    if (!ptr) {
      ptr = std::make_unique<ModeIndexND>();
      ptr->Init(boxes_, axes_, pool_, axis_offset_ + 1, ChildMask_());
    }
    return ptr.get();
  }

  const ModeIndexND* FindChild_(u32 node) const {
    const ChildMap& m = std::get<ChildMap>(store_);
    auto it = m.find(node);
    if (it == m.end()) return nullptr;
    return it->second.get();
  }

  const std::vector<BoxT>* boxes_{nullptr};
  const std::vector<AxisData>* axes_{nullptr};
  OccPool* pool_{nullptr};

  int axis_offset_{0};
  int remaining_{0};
  u32 mask_{0};

  // Leaf: BucketMap. Non-leaf: ChildMap.
  std::variant<BucketMap, ChildMap> store_{BucketMap{}};
};

}  // namespace index
}  // namespace sjs
