#pragma once
// sjs/index/segtree_rankrange_dyn.h
//
// §2.3 Dynamic rank-range segment tree.
//
// Universe is the static set of handles U with keys a(u)=L(u) (left endpoints).
// We build a RankSpace that maps each handle u to a unique rank in [0,N).
//
// Updates:
//  - INSERT(u): put u into buckets on Path(rank(u)) (leaf->root), O(log N).
// Queries:
//  - Given open interval (ell, r), map to rank range [L,R) using upper/lower bound,
//    then decompose into canonical cover nodes U(L,R) and aggregate buckets.
//  - Disjointness is guaranteed by Lemma 2.6 in SJS-HighDims.md.
//
// Supports COUNT/REPORT/SAMPLE with alias sampling over cover nodes.

#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/core/types.h"

#include "sjs/index/dense_bucket.h"
#include "sjs/index/occ_pool.h"
#include "sjs/index/rank_space.h"
#include "sjs/index/segtree_common.h"

#include "sjs/sampling/alias_table.h"

#include <memory>
#include <unordered_map>
#include <vector>

namespace sjs {
namespace index {

template <class ValueT = Scalar, class HandleT = Id>
class SegTreeRankRangeDyn {
 public:
  using Value = ValueT;
  using Handle = HandleT;
  using Bucket = DenseBucket<Handle>;
  using RankSpaceT = RankSpace<Value, Handle>;

  SegTreeRankRangeDyn() = default;

  // rank: shared rank space for this axis (built from all handles' a(u)=L(u)).
  void Init(std::shared_ptr<const RankSpaceT> rank, OccPool* pool) {
    rank_ = std::move(rank);
    pool_ = pool;
    buckets_.clear();

    SJS_ASSERT(rank_ != nullptr);
    SJS_ASSERT(pool_ != nullptr);

    leaves_ = static_cast<u32>(rank_->Size());  // N
    base_ = NextPow2(leaves_);
  }

  bool IsInitialized() const noexcept { return rank_ != nullptr && pool_ != nullptr; }

  u32 Leaves() const noexcept { return leaves_; }
  u32 Base() const noexcept { return base_; }

  // INSERT(u): insert handle u into all buckets on its leaf->root path.
  void Insert(Handle u) {
    if (!IsInitialized()) return;
    if (leaves_ == 0) return;

    const u32 r = rank_->RankOf(u);
    if (r >= leaves_) return;

    std::vector<u32> path;
    PathToRoot(r, base_, &path);

    for (u32 v : path) {
      Bucket* bkt = GetOrCreateBucket_(v);
      bkt->Insert(u);
    }
  }

  // COUNT(ell,r): open interval (ell,r)
  u64 Count(Value ell, Value r) const {
    if (!IsInitialized()) return 0;
    if (!(ell < r)) return 0;

    std::vector<u32> cover;
    QueryCoverNodes_(ell, r, &cover);
    if (cover.empty()) return 0;

    u64 sum = 0;
    for (u32 v : cover) {
      const Bucket* bkt = FindBucket_(v);
      if (bkt) sum += static_cast<u64>(bkt->Size());
    }
    return sum;
  }

  void Report(Value ell, Value r, std::vector<Handle>* out) const {
    SJS_ASSERT(out != nullptr);
    if (!IsInitialized()) return;
    if (!(ell < r)) return;

    std::vector<u32> cover;
    QueryCoverNodes_(ell, r, &cover);

    for (u32 v : cover) {
      const Bucket* bkt = FindBucket_(v);
      if (bkt) bkt->Report(out);
    }
  }

  void Sample(Value ell, Value r, u64 k, Rng* rng, std::vector<Handle>* out) const {
    SJS_ASSERT(rng != nullptr);
    SJS_ASSERT(out != nullptr);
    out->clear();
    if (k == 0) return;
    if (!IsInitialized()) return;
    if (!(ell < r)) return;

    std::vector<u32> cover;
    QueryCoverNodes_(ell, r, &cover);
    if (cover.empty()) return;

    std::vector<u64> w(cover.size(), 0);
    u64 total = 0;
    for (usize i = 0; i < cover.size(); ++i) {
      const Bucket* bkt = FindBucket_(cover[i]);
      const u64 sz = bkt ? static_cast<u64>(bkt->Size()) : 0ULL;
      w[i] = sz;
      total += sz;
    }
    if (total == 0) return;

    sampling::AliasTable alias;
    (void)alias.BuildFromU64(Span<const u64>(w.data(), w.size()));

    out->assign(static_cast<usize>(k), Handle{});
    for (u64 t = 0; t < k; ++t) {
      const usize idx = alias.Sample(rng);
      const Bucket* bkt = FindBucket_(cover[idx]);
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

  void QueryCoverNodes_(Value ell, Value r, std::vector<u32>* out) const {
    SJS_ASSERT(out != nullptr);
    out->clear();

    if (leaves_ == 0) return;

    const auto lr = rank_->OpenToRankRange(ell, r);  // [L,R)
    const u32 L = lr.first;
    const u32 R = lr.second;
    if (L >= R) return;

    DecomposeCoverOrdered(L, R, base_, out);
  }

  std::shared_ptr<const RankSpaceT> rank_{};
  OccPool* pool_{nullptr};

  u32 leaves_{0};
  u32 base_{1};

  std::unordered_map<u32, std::unique_ptr<Bucket>> buckets_;
};

}  // namespace index
}  // namespace sjs
