#pragma once
// sjs/index/segtree_stabbing_dyn.h
//
// §2.2 Dynamic stabbing segment tree (half-open intervals).
//
// Maintains a dynamic subset of intervals [L,R) over a fixed coordinate universe X.
// Query is a point x; answer set = { u | L(u) <= x < R(u) }.
//
// Canonical cover insertion + query root-to-leaf path yields a disjoint decomposition
// (Lemma 2.3 / Corollary 2.4 in SJS-HighDims.md), enabling exact COUNT/REPORT/SAMPLE.
//
// Notes:
//  - This class stores DenseBucket<Handle> at segment-tree nodes created on-demand.
//  - Deletion is handled globally via OccPool::EraseAll(handle) (shared pool).
//  - Coordinates universe X is shared via std::shared_ptr to avoid duplication.

#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/core/types.h"

#include "sjs/index/dense_bucket.h"
#include "sjs/index/occ_pool.h"
#include "sjs/index/segtree_common.h"

#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <vector>

namespace sjs {
namespace index {

template <class ValueT = Scalar, class HandleT = Id>
class SegTreeStabbingDyn {
 public:
  using Value = ValueT;
  using Handle = HandleT;
  using Bucket = DenseBucket<Handle>;

  SegTreeStabbingDyn() = default;

  // coords: sorted unique endpoint coordinates X = {x0 < x1 < ... < xM}.
  // Leaves represent elementary segments [x_i, x_{i+1}) for i in [0, M-1].
  void Init(std::shared_ptr<const std::vector<Value>> coords, OccPool* pool) {
    coords_ = std::move(coords);
    pool_ = pool;
    buckets_.clear();

    SJS_ASSERT(coords_ != nullptr);
    SJS_ASSERT_MSG(coords_->size() >= 2, "SegTreeStabbingDyn: coords must have >= 2 values");
    SJS_ASSERT(pool_ != nullptr);

    leaves_ = static_cast<u32>(coords_->size() - 1);   // number of elementary segments
    base_ = NextPow2(leaves_);
  }

  bool IsInitialized() const noexcept { return coords_ != nullptr && pool_ != nullptr; }

  u32 Leaves() const noexcept { return leaves_; }
  u32 Base() const noexcept { return base_; }

  const std::vector<Value>& Coords() const {
    SJS_ASSERT(coords_ != nullptr);
    return *coords_;
  }

  // Insert interval [a,b) for handle h.
  // Requires a < b, and a,b are endpoints from the universe used to build coords.
  void InsertInterval(Value a, Value b, Handle h) {
    if (!IsInitialized()) return;
    if (!(a < b)) return;

    const u32 l = EndpointIndex_(a);
    const u32 r = EndpointIndex_(b);
    if (l >= r) return;

    std::vector<u32> cover;
    DecomposeCover(l, r, base_, &cover);

    for (u32 v : cover) {
      Bucket* bkt = GetOrCreateBucket_(v);
      bkt->Insert(h);
    }
  }

  // COUNT(x): number of intervals containing point x.
  u64 Count(Value x) const {
    if (!IsInitialized()) return 0;

    std::vector<u32> path;
    QueryPathNodes_(x, &path);
    if (path.empty()) return 0;

    u64 sum = 0;
    for (u32 v : path) {
      const Bucket* bkt = FindBucket_(v);
      if (bkt) sum += static_cast<u64>(bkt->Size());
    }
    return sum;
  }

  // REPORT(x): append all handles of intervals containing x (no duplicates).
  void Report(Value x, std::vector<Handle>* out) const {
    SJS_ASSERT(out != nullptr);
    if (!IsInitialized()) return;

    std::vector<u32> path;
    QueryPathNodes_(x, &path);

    for (u32 v : path) {
      const Bucket* bkt = FindBucket_(v);
      if (bkt) bkt->Report(out);
    }
  }

  // SAMPLE(x,k): return k i.i.d. uniform samples (with replacement) from intervals containing x.
  // Implements Theorem 2.5 by alias sampling over disjoint buckets on the query path.
  void Sample(Value x, u64 k, Rng* rng, std::vector<Handle>* out) const {
    SJS_ASSERT(rng != nullptr);
    SJS_ASSERT(out != nullptr);
    out->clear();
    if (k == 0) return;
    if (!IsInitialized()) return;

    std::vector<u32> path;
    QueryPathNodes_(x, &path);
    if (path.empty()) return;

    std::vector<u64> w(path.size(), 0);
    u64 total = 0;
    for (usize i = 0; i < path.size(); ++i) {
      const Bucket* bkt = FindBucket_(path[i]);
      const u64 sz = bkt ? static_cast<u64>(bkt->Size()) : 0ULL;
      w[i] = sz;
      total += sz;
    }
    if (total == 0) return;

    sampling::AliasTable alias;
    (void)alias.BuildFromU64(Span<const u64>(w.data(), w.size()));

    out->assign(static_cast<usize>(k), Handle{});
    // For base 1D stabbing, grouping is optional; we implement direct i.i.d. sampling.
    for (u64 t = 0; t < k; ++t) {
      const usize idx = alias.Sample(rng);
      const Bucket* bkt = FindBucket_(path[idx]);
      SJS_DASSERT(bkt != nullptr && !bkt->Empty());
      Handle pick{};
      const bool ok = bkt->SampleOne(rng, &pick);
      SJS_DASSERT(ok);
      (*out)[static_cast<usize>(t)] = pick;
    }
  }

 private:
  Bucket* GetOrCreateBucket_(u32 node) {
    auto& ptr = buckets_[node];
    if (!ptr) ptr = std::make_unique<Bucket>(pool_);
    return ptr.get();
  }

  const Bucket* FindBucket_(u32 node) const {
    auto it = buckets_.find(node);
    if (it == buckets_.end()) return nullptr;
    return it->second.get();
  }

  // Endpoint index in coords (must exist for interval endpoints from the universe).
  u32 EndpointIndex_(Value v) const {
    const auto& X = *coords_;
    auto it = std::lower_bound(X.begin(), X.end(), v);
    SJS_DASSERT_MSG(it != X.end(), "EndpointIndex_: endpoint should be within universe");
    const u32 idx = static_cast<u32>(it - X.begin());
    // In debug: endpoint must be exactly present (interval endpoints are from universe).
    SJS_DASSERT_MSG(idx < X.size() && (*it == v), "EndpointIndex_: endpoint not found exactly in coords");
    return idx;
  }

  // Segment index i such that X[i] <= x < X[i+1].
  // Returns leaves_ as invalid if x is out of range.
  u32 SegmentIndex_(Value x) const {
    const auto& X = *coords_;
    if (X.size() < 2) return leaves_;
    auto it = std::upper_bound(X.begin(), X.end(), x);
    if (it == X.begin()) return leaves_;
    const u32 i = static_cast<u32>((it - X.begin()) - 1);
    if (i >= leaves_) return leaves_;
    return i;
  }

  void QueryPathNodes_(Value x, std::vector<u32>* out) const {
    SJS_ASSERT(out != nullptr);
    out->clear();

    const u32 seg = SegmentIndex_(x);
    if (seg >= leaves_) return;

    PathFromRoot(seg, base_, out);
  }

  std::shared_ptr<const std::vector<Value>> coords_{};
  OccPool* pool_{nullptr};

  u32 leaves_{0};
  u32 base_{1};

  // node-index -> bucket (heap allocated for stable pointer in OccPool)
  std::unordered_map<u32, std::unique_ptr<Bucket>> buckets_;
};

}  // namespace index
}  // namespace sjs
