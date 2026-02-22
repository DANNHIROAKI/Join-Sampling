#pragma once
// sjs/baselines/range_tree/sampling.h
//
// RangeTree baseline (HighDims) — Sampling variant (Framework II).
//
// This header implements the Chapter-4 baseline from SJS-HighDims.md:
//   - Sweep on axis 0 with END-before-START (half-open semantics).
//   - At each START(q), the partner set K_e consists of active boxes on the
//     opposite side whose projections on dims 1..Dim-1 intersect q.
//   - Represent each box b by an embedded point p(b) in D=2*(Dim-1) dims:
//       p(b) = (L_2..L_Dim, R_2..R_Dim)
//     where L_j = lo[j], R_j = hi[j].
//   - Intersection on dims 1..Dim-1 becomes an orthogonal (open) range query
//     in embedded space:
//       for each j in {2..Dim}:
//         L_j(b) < R_j(q)   and   R_j(b) > L_j(q)
//
// We implement a dynamic D-dim range tree with a Treap 1D engine (k=1 case)
// supporting open interval COUNT/REPORT/SAMPLE. This provides event-block
// primitives:
//   - COUNT(e)  = |K_e|
//   - REPORT(e) enumerates K_e once
//   - SAMPLE(e,k) returns k i.i.d. uniform (with replacement) samples from K_e
//
// NOTE: This is a correctness-focused baseline. Range trees have
// O(N log^{D-1} N) worst-case space and are expected to become impractical
// for larger Dim.

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/join/sweep_events.h"
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
namespace range_tree {

namespace detail {

// --------------------------
// Open interval helper
// --------------------------

template <class T>
struct OpenInterval {
  bool has_lo = false;
  bool has_hi = false;
  T lo{};
  T hi{};
};

// --------------------------
// 1D key with deterministic tie-break (x, id)
// --------------------------

template <class T>
struct Key1D {
  T x{};
  i64 id = 0;  // unique id for tie-break (use box handle)
};

template <class T>
inline bool operator<(const Key1D<T>& a, const Key1D<T>& b) noexcept {
  if (a.x < b.x) return true;
  if (b.x < a.x) return false;
  return a.id < b.id;
}

template <class T>
inline bool operator<=(const Key1D<T>& a, const Key1D<T>& b) noexcept {
  return !(b < a);
}

template <class T>
inline bool operator==(const Key1D<T>& a, const Key1D<T>& b) noexcept {
  return a.id == b.id && !(a.x < b.x) && !(b.x < a.x);
}

inline constexpr i64 kIdMinSentinel = std::numeric_limits<i64>::min();
inline constexpr i64 kIdMaxSentinel = std::numeric_limits<i64>::max();

// Deterministic "random" priority from id (splitmix64 finalizer).
inline u64 PriorityFromId(i64 id) noexcept {
  u64 x = static_cast<u64>(id);
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  x = x ^ (x >> 31);
  return x;
}

// --------------------------
// Treap (1D): dynamic set of Key1D<T>
// Supports open interval COUNT/REPORT/SAMPLE via split/merge.
// --------------------------

template <class T>
class Treap1D {
 public:
  Treap1D() = default;
  Treap1D(const Treap1D&) = delete;
  Treap1D& operator=(const Treap1D&) = delete;

  ~Treap1D() { Clear(); }

  void Clear() {
    Destroy(root_);
    root_ = nullptr;
  }

  bool Empty() const noexcept { return root_ == nullptr; }
  u64 Size() const noexcept { return static_cast<u64>(SizeOf(root_)); }

  void Insert(const Key1D<T>& key) {
    // Standard split-based insert.
    Node* a = nullptr;
    Node* b = nullptr;
    Split(root_, key, a, b);  // a: <= key, b: > key

    // Optional de-dup: if key already exists, it is in 'a'.
    // We do not remove it; insertion order is controlled by active flags.
    Node* n = new Node(key);
    root_ = Merge(Merge(a, n), b);
  }

  void Erase(const Key1D<T>& key) {
    // Split out <=key and >key
    Node* a = nullptr;
    Node* b = nullptr;
    Split(root_, key, a, b);

    // Split a into <=key_prev and (key_prev, key]
    Key1D<T> key_prev{key.x, key.id - 1};
    Node* a1 = nullptr;
    Node* mid = nullptr;
    Split(a, key_prev, a1, mid);

    // mid should contain only 'key' (unique). Remove it safely.
    if (mid) {
      Node* merged = Merge(mid->l, mid->r);
      mid->l = mid->r = nullptr;
      delete mid;
      mid = merged;
    }

    root_ = Merge(Merge(a1, mid), b);
  }

  // Count keys with x in (lo, hi) under open semantics.
  u64 Count(const OpenInterval<T>& I) {
    if (!root_) return 0ULL;

    Node* a = nullptr;
    Node* b = nullptr;
    Node* c = nullptr;

    if (I.has_lo) {
      const Key1D<T> kl{I.lo, kIdMaxSentinel};
      Split(root_, kl, a, b);  // b: > lo
    } else {
      b = root_;
    }

    if (I.has_hi) {
      const Key1D<T> kh{I.hi, kIdMinSentinel};
      Split(b, kh, b, c);  // b: < hi
    }

    const u64 ans = static_cast<u64>(SizeOf(b));

    // Merge back.
    if (I.has_hi) b = Merge(b, c);
    if (I.has_lo) root_ = Merge(a, b);

    return ans;
  }

  // Report all ids in the open interval; appends to 'out'.
  void Report(const OpenInterval<T>& I, std::vector<u32>* out) {
    if (!out) return;
    if (!root_) return;

    Node* a = nullptr;
    Node* b = nullptr;
    Node* c = nullptr;

    if (I.has_lo) {
      const Key1D<T> kl{I.lo, kIdMaxSentinel};
      Split(root_, kl, a, b);
    } else {
      b = root_;
    }

    if (I.has_hi) {
      const Key1D<T> kh{I.hi, kIdMinSentinel};
      Split(b, kh, b, c);
    }

    InorderAppend(b, out);

    // Merge back.
    if (I.has_hi) b = Merge(b, c);
    if (I.has_lo) root_ = Merge(a, b);
  }

  // Sample one id uniformly from the open interval.
  // Returns false if empty.
  bool Sample(const OpenInterval<T>& I, Rng* rng, u32* out_id) {
    if (!rng || !out_id) return false;
    if (!root_) return false;

    Node* a = nullptr;
    Node* b = nullptr;
    Node* c = nullptr;

    if (I.has_lo) {
      const Key1D<T> kl{I.lo, kIdMaxSentinel};
      Split(root_, kl, a, b);
    } else {
      b = root_;
    }

    if (I.has_hi) {
      const Key1D<T> kh{I.hi, kIdMinSentinel};
      Split(b, kh, b, c);
    }

    const u32 n = SizeOf(b);
    bool ok = false;
    if (n > 0) {
      const u32 k = static_cast<u32>(rng->UniformU64(static_cast<u64>(n)));
      const Node* kth = Kth(b, k);
      if (kth) {
        *out_id = static_cast<u32>(kth->key.id);
        ok = true;
      }
    }

    // Merge back.
    if (I.has_hi) b = Merge(b, c);
    if (I.has_lo) root_ = Merge(a, b);

    return ok;
  }

 private:
  struct Node {
    explicit Node(const Key1D<T>& k) : key(k), prio(PriorityFromId(k.id)) {}

    Key1D<T> key;
    u64 prio = 0;
    Node* l = nullptr;
    Node* r = nullptr;
    u32 sz = 1;
  };

  static u32 SizeOf(const Node* n) noexcept { return n ? n->sz : 0U; }

  static void Pull(Node* n) noexcept {
    if (!n) return;
    n->sz = 1U + SizeOf(n->l) + SizeOf(n->r);
  }

  static void Destroy(Node* n) {
    if (!n) return;
    Destroy(n->l);
    Destroy(n->r);
    delete n;
  }

  // Split by key: left gets keys <= key, right gets keys > key.
  static void Split(Node* cur, const Key1D<T>& key, Node*& left, Node*& right) {
    if (!cur) {
      left = nullptr;
      right = nullptr;
      return;
    }
    if (cur->key <= key) {
      Split(cur->r, key, cur->r, right);
      left = cur;
      Pull(left);
    } else {
      Split(cur->l, key, left, cur->l);
      right = cur;
      Pull(right);
    }
  }

  static Node* Merge(Node* a, Node* b) {
    if (!a) return b;
    if (!b) return a;
    if (a->prio < b->prio) {
      a->r = Merge(a->r, b);
      Pull(a);
      return a;
    } else {
      b->l = Merge(a, b->l);
      Pull(b);
      return b;
    }
  }

  static const Node* Kth(const Node* cur, u32 k) {
    // 0-based.
    if (!cur) return nullptr;
    const u32 left_sz = SizeOf(cur->l);
    if (k < left_sz) return Kth(cur->l, k);
    if (k == left_sz) return cur;
    return Kth(cur->r, k - left_sz - 1);
  }

  static void InorderAppend(const Node* cur, std::vector<u32>* out) {
    if (!cur || !out) return;
    // Iterative to avoid deep recursion in adversarial cases.
    std::vector<const Node*> st;
    st.reserve(64);
    const Node* x = cur;
    while (x || !st.empty()) {
      while (x) {
        st.push_back(x);
        x = x->l;
      }
      x = st.back();
      st.pop_back();
      out->push_back(static_cast<u32>(x->key.id));
      x = x->r;
    }
  }

  Node* root_ = nullptr;
};

// --------------------------
// RangeTree (k dims, Offset in a DTotal-dim point)
// Segment-tree representation of the primary dimension.
// Associated structure at each node is a (k-1)-dim range tree.
// Base case (k=1) is Treap1D.
// --------------------------

template <int DTotal, class T>
struct RTPoint {
  Point<DTotal, T> coords{};
  i64 id = 0;  // unique (handle)
};

template <int K, int Offset, int DTotal, class T>
class RangeTree;

// Base case: 1D treap.
template <int Offset, int DTotal, class T>
class RangeTree<1, Offset, DTotal, T> {
 public:
  using PointT = RTPoint<DTotal, T>;

  RangeTree() = default;
  RangeTree(const RangeTree&) = delete;
  RangeTree& operator=(const RangeTree&) = delete;

  void Reset() {
    pts_ = nullptr;
    treap_.Clear();
  }

  void Build(const std::vector<PointT>* pts) {
    pts_ = pts;
    treap_.Clear();
  }

  void Insert(u32 pidx) {
    const auto& p = (*pts_)[static_cast<usize>(pidx)];
    treap_.Insert(Key1D<T>{p.coords.v[static_cast<usize>(Offset)], p.id});
  }

  void Erase(u32 pidx) {
    const auto& p = (*pts_)[static_cast<usize>(pidx)];
    treap_.Erase(Key1D<T>{p.coords.v[static_cast<usize>(Offset)], p.id});
  }

  struct Block {
    Treap1D<T>* treap = nullptr;
    OpenInterval<T> I;
  };

  template <class IntervalArray>
  void Blocks(const IntervalArray& q, std::vector<Block>* out) {
    if (!out) return;
    Block b;
    b.treap = &treap_;
    b.I = q[static_cast<usize>(Offset)];
    out->push_back(b);
  }

 private:
  const std::vector<PointT>* pts_ = nullptr;
  Treap1D<T> treap_;
};

// General case: k>=2.
template <int K, int Offset, int DTotal, class T>
class RangeTree {
  static_assert(K >= 2, "RangeTree<K>: K must be >= 2 in this specialization");

 public:
  using PointT = RTPoint<DTotal, T>;
  using Child = RangeTree<K - 1, Offset + 1, DTotal, T>;
  using Block = typename RangeTree<1, Offset + (K - 1), DTotal, T>::Block;  // base Block type

  RangeTree() = default;
  RangeTree(const RangeTree&) = delete;
  RangeTree& operator=(const RangeTree&) = delete;

  void Reset() {
    pts_ = nullptr;
    order_.clear();
    nodes_.clear();
    root_ = -1;
  }

  void Build(const std::vector<PointT>* pts, std::vector<u32> indices) {
    pts_ = pts;
    order_ = std::move(indices);
    std::sort(order_.begin(), order_.end(), [&](u32 a, u32 b) {
      const auto& pa = (*pts_)[static_cast<usize>(a)];
      const auto& pb = (*pts_)[static_cast<usize>(b)];
      const Key1D<T> ka{pa.coords.v[static_cast<usize>(Offset)], pa.id};
      const Key1D<T> kb{pb.coords.v[static_cast<usize>(Offset)], pb.id};
      return ka < kb;
    });

    nodes_.clear();
    nodes_.reserve(order_.empty() ? 0 : (2 * order_.size()));
    root_ = BuildNode(0U, static_cast<u32>(order_.size()));
  }

  void Insert(u32 pidx) {
    if (order_.empty()) return;
    const u32 pos = FindPos(pidx);
    InsertRec(root_, pos, pidx);
  }

  void Erase(u32 pidx) {
    if (order_.empty()) return;
    const u32 pos = FindPos(pidx);
    EraseRec(root_, pos, pidx);
  }

  template <class IntervalArray>
  void Blocks(const IntervalArray& q, std::vector<typename Child::Block>* out) {
    if (!out) return;
    if (order_.empty()) return;

    const auto& I = q[static_cast<usize>(Offset)];
    u32 L = 0;
    u32 R = static_cast<u32>(order_.size());

    if (I.has_lo) {
      const Key1D<T> kl{I.lo, kIdMaxSentinel};
      L = UpperBoundKey(kl);
    }
    if (I.has_hi) {
      const Key1D<T> kh{I.hi, kIdMinSentinel};
      R = LowerBoundKey(kh);
    }

    if (L >= R) return;
    BlocksRec(root_, L, R, q, out);
  }

 private:
  struct Node {
    u32 l = 0;
    u32 r = 0;
    i32 left = -1;
    i32 right = -1;
    std::unique_ptr<Child> assoc;
  };

  const std::vector<PointT>* pts_ = nullptr;
  std::vector<u32> order_;
  std::vector<Node> nodes_;
  i32 root_ = -1;

  i32 BuildNode(u32 l, u32 r) {
    const i32 id = static_cast<i32>(nodes_.size());
    nodes_.push_back(Node{});
    Node& n = nodes_.back();
    n.l = l;
    n.r = r;

    // Build associated (K-1)-dim tree on this node's point set.
    n.assoc = std::make_unique<Child>();
    if constexpr (K - 1 == 1) {
      // Base case: no need to pass indices.
      n.assoc->Build(pts_);
    } else {
      std::vector<u32> sub;
      sub.reserve(static_cast<usize>(r - l));
      for (u32 i = l; i < r; ++i) sub.push_back(order_[static_cast<usize>(i)]);
      n.assoc->Build(pts_, std::move(sub));
    }

    if (r - l <= 1) {
      return id;
    }

    const u32 mid = l + ((r - l) >> 1);
    n.left = BuildNode(l, mid);
    n.right = BuildNode(mid, r);
    return id;
  }

  Key1D<T> KeyOf(u32 pidx) const {
    const auto& p = (*pts_)[static_cast<usize>(pidx)];
    return Key1D<T>{p.coords.v[static_cast<usize>(Offset)], p.id};
  }

  u32 FindPos(u32 pidx) const {
    const Key1D<T> k = KeyOf(pidx);
    const auto it = std::lower_bound(order_.begin(), order_.end(), k, [&](u32 idx, const Key1D<T>& key) {
      const auto& p = (*pts_)[static_cast<usize>(idx)];
      const Key1D<T> kk{p.coords.v[static_cast<usize>(Offset)], p.id};
      return kk < key;
    });
    // Key must exist in this tree.
    SJS_DASSERT(it != order_.end());
    return static_cast<u32>(it - order_.begin());
  }

  u32 LowerBoundKey(const Key1D<T>& k) const {
    const auto it = std::lower_bound(order_.begin(), order_.end(), k, [&](u32 idx, const Key1D<T>& key) {
      const auto& p = (*pts_)[static_cast<usize>(idx)];
      const Key1D<T> kk{p.coords.v[static_cast<usize>(Offset)], p.id};
      return kk < key;
    });
    return static_cast<u32>(it - order_.begin());
  }

  u32 UpperBoundKey(const Key1D<T>& k) const {
    const auto it = std::upper_bound(order_.begin(), order_.end(), k, [&](const Key1D<T>& key, u32 idx) {
      const auto& p = (*pts_)[static_cast<usize>(idx)];
      const Key1D<T> kk{p.coords.v[static_cast<usize>(Offset)], p.id};
      return key < kk;
    });
    return static_cast<u32>(it - order_.begin());
  }

  void InsertRec(i32 node_id, u32 pos, u32 pidx) {
    Node& n = nodes_[static_cast<usize>(node_id)];
    n.assoc->Insert(pidx);
    if (n.left < 0) return;  // leaf
    const u32 mid = n.l + ((n.r - n.l) >> 1);
    if (pos < mid) {
      InsertRec(n.left, pos, pidx);
    } else {
      InsertRec(n.right, pos, pidx);
    }
  }

  void EraseRec(i32 node_id, u32 pos, u32 pidx) {
    Node& n = nodes_[static_cast<usize>(node_id)];
    n.assoc->Erase(pidx);
    if (n.left < 0) return;
    const u32 mid = n.l + ((n.r - n.l) >> 1);
    if (pos < mid) {
      EraseRec(n.left, pos, pidx);
    } else {
      EraseRec(n.right, pos, pidx);
    }
  }

  template <class IntervalArray>
  void BlocksRec(i32 node_id,
                 u32 ql,
                 u32 qr,
                 const IntervalArray& q,
                 std::vector<typename Child::Block>* out) {
    Node& n = nodes_[static_cast<usize>(node_id)];
    if (qr <= n.l || n.r <= ql) return;
    if (ql <= n.l && n.r <= qr) {
      n.assoc->Blocks(q, out);
      return;
    }
    if (n.left < 0) {
      // leaf segment (single point) that overlaps
      n.assoc->Blocks(q, out);
      return;
    }
    BlocksRec(n.left, ql, qr, q, out);
    BlocksRec(n.right, ql, qr, q, out);
  }
};

// --------------------------
// RangeTreeIndex: wraps a D=2*(Dim-1) range tree for one relation.
// Dynamic active set is maintained via Insert/Erase.
// Query builds the Chapter-4 orthogonal constraints.
// --------------------------

template <int Dim, class T>
class RangeTreeIndex {
  static_assert(Dim >= 2, "RangeTreeIndex requires Dim>=2");

 public:
  static constexpr int m = Dim - 1;
  static constexpr int D = 2 * (Dim - 1);

  using BoxT = Box<Dim, T>;
  using PointT = RTPoint<D, T>;
  using TreeT = RangeTree<D, 0, D, T>;
  using BlockT = typename RangeTree<1, D - 1, D, T>::Block;  // same Block type for all

  RangeTreeIndex() = default;

  void Reset() {
    pts_.clear();
    tree_.Reset();
  }

  bool BuildFromRelation(const Relation<Dim, T>& rel) {
    const usize n = rel.Size();
    pts_.resize(n);
    for (usize i = 0; i < n; ++i) {
      const auto& b = rel.boxes[i];
      Point<D, T> p;
      const usize mm = static_cast<usize>(m);
      // coords[0..m-1] = lo[1..Dim-1]
      // coords[m..2m-1] = hi[1..Dim-1]
      for (int j = 1; j < Dim; ++j) {
        const usize idx = static_cast<usize>(j - 1);
        const usize jj = static_cast<usize>(j);
        p.v[idx] = b.lo.v[jj];
        p.v[idx + mm] = b.hi.v[jj];
      }
      pts_[i].coords = p;
      pts_[i].id = static_cast<i64>(i);  // unique id = handle
    }

    std::vector<u32> idx;
    idx.reserve(n);
    for (usize i = 0; i < n; ++i) idx.push_back(static_cast<u32>(i));
    tree_.Build(&pts_, std::move(idx));
    return true;
  }

  void Insert(u32 handle) { tree_.Insert(handle); }
  void Erase(u32 handle) { tree_.Erase(handle); }

  // COUNT: number of active points whose embedded coords satisfy the query.
  u64 Count(const BoxT& q) {
    std::array<OpenInterval<T>, D> Q{};
    BuildQuery(q, &Q);

    using Block = typename RangeTree<1, D - 1, D, T>::Block;
    std::vector<Block> blocks;
    blocks.reserve(64);
    tree_.Blocks(Q, &blocks);

    u64 total = 0;
    for (auto& b : blocks) {
      total += b.treap->Count(b.I);
    }
    return total;
  }

  // REPORT: enumerate all active handles satisfying the query (order arbitrary).
  void Report(const BoxT& q, std::vector<u32>* out) {
    if (!out) return;
    std::array<OpenInterval<T>, D> Q{};
    BuildQuery(q, &Q);

    using Block = typename RangeTree<1, D - 1, D, T>::Block;
    std::vector<Block> blocks;
    blocks.reserve(64);
    tree_.Blocks(Q, &blocks);

    for (auto& b : blocks) {
      b.treap->Report(b.I, out);
    }
  }

  // SAMPLE: k i.i.d. uniform (with replacement) handles satisfying the query.
  void Sample(const BoxT& q, u32 k, Rng* rng, std::vector<u32>* out) {
    if (!out || k == 0) return;
    out->clear();
    if (!rng) return;

    std::array<OpenInterval<T>, D> Q{};
    BuildQuery(q, &Q);

    using Block = typename RangeTree<1, D - 1, D, T>::Block;
    std::vector<Block> blocks;
    blocks.reserve(64);
    tree_.Blocks(Q, &blocks);
    if (blocks.empty()) return;

    std::vector<u64> w;
    w.reserve(blocks.size());
    u64 total = 0;
    for (auto& b : blocks) {
      const u64 c = b.treap->Count(b.I);
      w.push_back(c);
      total += c;
    }
    if (total == 0) return;

    // Prefix sums for weighted block selection.
    std::vector<u64> pref;
    pref.reserve(w.size());
    u64 run = 0;
    for (u64 x : w) {
      run += x;
      pref.push_back(run);
    }

    out->reserve(k);
    for (u32 i = 0; i < k; ++i) {
      const u64 t = rng->UniformU64(total) + 1ULL;  // 1..total
      const auto it = std::lower_bound(pref.begin(), pref.end(), t);
      const usize bi = static_cast<usize>(it - pref.begin());
      SJS_DASSERT(bi < blocks.size());
      u32 hid = 0;
      const bool ok = blocks[bi].treap->Sample(blocks[bi].I, rng, &hid);
      if (!ok) {
        // Should not happen if weights are consistent, but be defensive.
        continue;
      }
      out->push_back(hid);
    }
  }

 private:
  static void BuildQuery(const BoxT& q, std::array<OpenInterval<T>, D>* out) {
    SJS_DASSERT(out != nullptr);
    const usize mm = static_cast<usize>(m);

    // First m dims: L_j < R_j(q)  => (-inf, hi)
    for (int j = 1; j < Dim; ++j) {
      const usize idx = static_cast<usize>(j - 1);
      (*out)[idx].has_lo = false;
      (*out)[idx].has_hi = true;
      (*out)[idx].hi = q.hi.v[static_cast<usize>(j)];
    }

    // Last m dims: R_j > L_j(q) => (lo, +inf)
    for (int j = 1; j < Dim; ++j) {
      const usize idx = static_cast<usize>(j - 1);
      (*out)[idx + mm].has_lo = true;
      (*out)[idx + mm].lo = q.lo.v[static_cast<usize>(j)];
      (*out)[idx + mm].has_hi = false;
    }
  }

  std::vector<PointT> pts_;
  TreeT tree_;
};

// --------------------------
// Slot plan: group t output slots by event id (sid) using exact weights.
// --------------------------

struct SlotPlanND {
  // offsets size = E+1; slots size = t.
  // For event i, its assigned slots are slots[offset[i]..offset[i+1]).
  std::vector<u32> offset;
  std::vector<u32> slots;
};

inline bool BuildSlotPlanND(const std::vector<u64>& weights, u64 t, Rng* rng, SlotPlanND* out) {
  if (!out || !rng) return false;
  out->offset.clear();
  out->slots.clear();

  const usize E = weights.size();
  if (E == 0) return false;
  if (t == 0) {
    out->offset.assign(E + 1, 0U);
    return true;
  }
  if (t > static_cast<u64>(std::numeric_limits<u32>::max())) {
    return false;  // slot indices are stored as u32
  }

  sampling::AliasTable alias;
  if (!alias.BuildFromU64(weights)) {
    // All weights are zero.
    return false;
  }

  // 1) Draw event id per slot once.
  const u32 t32 = static_cast<u32>(t);
  std::vector<u32> eid_of_slot;
  eid_of_slot.resize(static_cast<usize>(t32));

  std::vector<u32> cnt(E, 0U);
  for (u32 slot = 0; slot < t32; ++slot) {
    const u32 sid = alias.Sample(rng);
    SJS_DASSERT(static_cast<usize>(sid) < E);
    eid_of_slot[static_cast<usize>(slot)] = sid;
    ++cnt[static_cast<usize>(sid)];
  }

  // 2) Prefix sums -> offsets.
  out->offset.resize(E + 1);
  out->offset[0] = 0U;
  for (usize i = 0; i < E; ++i) {
    out->offset[i + 1] = out->offset[i] + cnt[i];
  }

  // 3) Stable fill slots grouped by sid.
  out->slots.resize(static_cast<usize>(t32));
  std::vector<u32> cursor = out->offset;
  for (u32 slot = 0; slot < t32; ++slot) {
    const u32 sid = eid_of_slot[static_cast<usize>(slot)];
    const u32 pos = cursor[static_cast<usize>(sid)]++;
    out->slots[static_cast<usize>(pos)] = slot;
  }

  return true;
}

// --------------------------
// Context: sweep events + two dynamic RangeTreeIndex structures (R and S).
// --------------------------

template <int Dim, class T>
class RangeTreeContext {
 public:
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;

  bool built() const noexcept { return built_; }
  const DatasetT* dataset() const noexcept { return ds_; }

  const std::vector<join::Event>& events() const noexcept { return events_; }
  const std::vector<i32>& start_id_of_event() const noexcept { return start_id_of_event_; }
  usize num_start_events() const noexcept { return num_start_events_; }

  const RangeTreeIndex<Dim, T>& index_r() const noexcept { return idx_r_; }
  const RangeTreeIndex<Dim, T>& index_s() const noexcept { return idx_s_; }

  RangeTreeIndex<Dim, T>& active_r() noexcept { return idx_r_; }
  RangeTreeIndex<Dim, T>& active_s() noexcept { return idx_s_; }

  void Reset() {
    built_ = false;
    ds_ = nullptr;
    events_.clear();
    start_id_of_event_.clear();
    num_start_events_ = 0;

    idx_r_.Reset();
    idx_s_.Reset();

    active_r_flag_.clear();
    active_s_flag_.clear();
    active_r_list_.clear();
    active_s_list_.clear();
    active_r_sz_ = 0;
    active_s_sz_ = 0;
  }

  bool Build(const DatasetT& ds, PhaseRecorder* phases, std::string* err) {
    Reset();

    // Validate dataset (proper half-open boxes).
    {
      std::string tmp;
      if (!ds.Validate(/*require_proper=*/true, &tmp)) {
        if (err) *err = "RangeTreeContext::Build: invalid dataset: " + tmp;
        return false;
      }
    }

    // Basic size guards (we use u32 handles and u32 START ids).
    const usize nR = ds.R.Size();
    const usize nS = ds.S.Size();
    if (nR > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        nS > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RangeTreeContext::Build: relation size exceeds u32 handle range";
      return false;
    }
    if ((nR + nS) > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RangeTreeContext::Build: number of START events exceeds u32 range";
      return false;
    }

    ds_ = &ds;

    {
      auto _ = phases ? phases->Scoped("build_sweep_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = join::BuildSweepEvents<Dim, T>(ds, /*axis=*/0, join::SideTieBreak::RBeforeS);

      // Build dense START id mapping.
      start_id_of_event_.assign(events_.size(), -1);
      num_start_events_ = 0;
      for (usize pos = 0; pos < events_.size(); ++pos) {
        if (events_[pos].kind == join::EventKind::Start) {
          // num_start_events_ is at most nR+nS which we already bounded to u32.
          start_id_of_event_[pos] = static_cast<i32>(num_start_events_);
          ++num_start_events_;
        }
      }
    }

    {
      auto _ = phases ? phases->Scoped("build_range_trees") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!idx_r_.BuildFromRelation(ds.R)) {
        if (err) *err = "RangeTreeContext::Build: failed to build range tree for R";
        return false;
      }
      if (!idx_s_.BuildFromRelation(ds.S)) {
        if (err) *err = "RangeTreeContext::Build: failed to build range tree for S";
        return false;
      }
    }

    active_r_flag_.assign(ds.R.Size(), 0U);
    active_s_flag_.assign(ds.S.Size(), 0U);
    active_r_list_.clear();
    active_s_list_.clear();
    active_r_sz_ = 0;
    active_s_sz_ = 0;

    built_ = true;
    return true;
  }

  // Clear active sets (needed if a sweep is interrupted and restarted).
  void ResetActive() {
    // Delete any still-active handles.
    for (u32 h : active_r_list_) {
      const usize hu = static_cast<usize>(h);
      if (hu < active_r_flag_.size() && active_r_flag_[hu]) {
        idx_r_.Erase(h);
        active_r_flag_[hu] = 0U;
      }
    }
    for (u32 h : active_s_list_) {
      const usize hu = static_cast<usize>(h);
      if (hu < active_s_flag_.size() && active_s_flag_[hu]) {
        idx_s_.Erase(h);
        active_s_flag_[hu] = 0U;
      }
    }
    active_r_list_.clear();
    active_s_list_.clear();
    active_r_sz_ = 0;
    active_s_sz_ = 0;
  }

  // Access a box by event.
  const BoxT& BoxOf(const join::Event& ev) const {
    SJS_DASSERT(ds_ != nullptr);
    if (ev.side == join::Side::R) {
      return ds_->R.boxes[static_cast<usize>(ev.index)];
    }
    return ds_->S.boxes[static_cast<usize>(ev.index)];
  }

  // Insert/delete helpers.
  void Insert(join::Side side, u32 handle) {
    if (side == join::Side::R) {
      if (handle >= active_r_flag_.size()) return;
      if (active_r_flag_[handle]) return;
      active_r_flag_[handle] = 1U;
      active_r_list_.push_back(handle);
      idx_r_.Insert(handle);
      ++active_r_sz_;
    } else {
      if (handle >= active_s_flag_.size()) return;
      if (active_s_flag_[handle]) return;
      active_s_flag_[handle] = 1U;
      active_s_list_.push_back(handle);
      idx_s_.Insert(handle);
      ++active_s_sz_;
    }
  }

  void Delete(join::Side side, u32 handle) {
    if (side == join::Side::R) {
      if (handle >= active_r_flag_.size()) return;
      if (!active_r_flag_[handle]) return;
      active_r_flag_[handle] = 0U;
      idx_r_.Erase(handle);
      if (active_r_sz_ > 0) --active_r_sz_;
    } else {
      if (handle >= active_s_flag_.size()) return;
      if (!active_s_flag_[handle]) return;
      active_s_flag_[handle] = 0U;
      idx_s_.Erase(handle);
      if (active_s_sz_ > 0) --active_s_sz_;
    }
  }

  u32 active_size_r() const noexcept { return active_r_sz_; }
  u32 active_size_s() const noexcept { return active_s_sz_; }

 private:
  bool built_ = false;
  const DatasetT* ds_ = nullptr;

  std::vector<join::Event> events_;
  std::vector<i32> start_id_of_event_;
  usize num_start_events_ = 0;

  RangeTreeIndex<Dim, T> idx_r_;
  RangeTreeIndex<Dim, T> idx_s_;

  std::vector<u8> active_r_flag_;
  std::vector<u8> active_s_flag_;
  std::vector<u32> active_r_list_;
  std::vector<u32> active_s_list_;
  u32 active_r_sz_ = 0;
  u32 active_s_sz_ = 0;
};

// --------------------------
// Deterministic REPORT-based join enumerator (sweep + range tree REPORT).
// --------------------------

template <int Dim, class T>
class RangeTreeReportJoinEnumeratorND final : public IJoinEnumerator {
 public:
  using Ctx = RangeTreeContext<Dim, T>;
  using BoxT = Box<Dim, T>;

  explicit RangeTreeReportJoinEnumeratorND(Ctx* ctx) : ctx_(ctx) {
    SJS_DASSERT(ctx_ != nullptr);
    Reset();
  }

  void Reset() override {
    done_ = false;
    buf_.clear();
    buf_pos_ = 0;
    ev_pos_ = 0;
    stats_ = join::JoinStats{};

    if (ctx_) ctx_->ResetActive();
  }

  bool Next(PairId* out) override {
    if (!out) return false;
    if (done_ || !ctx_ || !ctx_->built() || ctx_->dataset() == nullptr) return false;

    const auto* ds = ctx_->dataset();

    while (true) {
      // Drain buffered partners for the current START event.
      if (buf_pos_ < buf_.size()) {
        const u32 oh = buf_[buf_pos_++];
        stats_.output_pairs++;

        if (cur_side_ == join::Side::R) {
          *out = PairId{ds->R.GetId(static_cast<usize>(cur_handle_)),
                         ds->S.GetId(static_cast<usize>(oh))};
        } else {
          *out = PairId{ds->R.GetId(static_cast<usize>(oh)),
                         ds->S.GetId(static_cast<usize>(cur_handle_))};
        }
        return true;
      }

      // Advance to next START event.
      buf_.clear();
      buf_pos_ = 0;

      if (ev_pos_ >= ctx_->events().size()) {
        done_ = true;
        ctx_->ResetActive();
        return false;
      }

      const auto& ev = ctx_->events()[ev_pos_++];
      stats_.num_events++;
      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        ctx_->Delete(ev.side, handle);
        continue;
      }

      // START: query before inserting.
      cur_side_ = ev.side;
      cur_handle_ = handle;

      const BoxT& q = ctx_->BoxOf(ev);
      if (ev.side == join::Side::R) {
        ctx_->active_s().Report(q, &buf_);
      } else {
        ctx_->active_r().Report(q, &buf_);
      }

      stats_.candidate_checks += static_cast<u64>(buf_.size());
      stats_.active_max_r = std::max(stats_.active_max_r, static_cast<u64>(ctx_->active_size_r()));
      stats_.active_max_s = std::max(stats_.active_max_s, static_cast<u64>(ctx_->active_size_s()));

      // Insert q after reporting.
      ctx_->Insert(ev.side, handle);

      // Loop; if buf_ empty we'll continue scanning.
    }
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  Ctx* ctx_ = nullptr;
  bool done_ = false;

  std::vector<u32> buf_;
  usize buf_pos_ = 0;
  usize ev_pos_ = 0;

  join::Side cur_side_ = join::Side::R;
  u32 cur_handle_ = 0;

  join::JoinStats stats_;
};

}  // namespace detail

// --------------------------
// Baseline: RangeTreeSamplingBaseline (Framework II)
// --------------------------

template <int Dim, class T = Scalar>
class RangeTreeSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "RangeTreeSamplingBaseline requires Dim>=2");

  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;
  using BoxT = Box<Dim, T>;

  Method method() const noexcept override { return Method::RangeTree; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "range_tree_sampling"; }

  void Reset() override {
    built_ = false;
    ctx_.Reset();

    w_total_.clear();
    W_ = 0;
    weights_valid_ = false;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    Reset();
    if (!ctx_.Build(ds, phases, err)) return false;

    const usize E = ctx_.num_start_events();
    w_total_.assign(E, 0ULL);
    built_ = true;
    return true;
  }

  bool Count(const Config& cfg, Rng* rng, CountResult* out, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)rng;
    if (!out) return false;
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "RangeTreeSamplingBaseline::Count: call Build() first";
      return false;
    }

    auto _ = phases ? phases->Scoped("pass1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    W_ = 0;

    ctx_.ResetActive();

    const auto& events = ctx_.events();
    const auto& sid_of_pos = ctx_.start_id_of_event();

    for (usize pos = 0; pos < events.size(); ++pos) {
      const auto& ev = events[pos];
      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        ctx_.Delete(ev.side, handle);
        continue;
      }

      // START
      const i32 sid_i32 = sid_of_pos[pos];
      SJS_DASSERT(sid_i32 >= 0);
      const usize sid = static_cast<usize>(sid_i32);
      SJS_DASSERT(sid < w_total_.size());

      const BoxT& q = ctx_.BoxOf(ev);
      u64 w = 0;
      if (ev.side == join::Side::R) {
        w = ctx_.active_s().Count(q);
      } else {
        w = ctx_.active_r().Count(q);
      }

      w_total_[sid] = w;
      W_ += w;

      // Insert after counting.
      ctx_.Insert(ev.side, handle);
    }

    ctx_.ResetActive();

    weights_valid_ = true;
    *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    if (!out || !rng) return false;
    out->Clear();

    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "RangeTreeSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!weights_valid_) {
      if (err) *err = "RangeTreeSamplingBaseline::Sample: call Count() first";
      return false;
    }

    const auto* ds = ctx_.dataset();
    const u64 t = cfg.run.t;

    out->pairs.resize(static_cast<usize>(t));
    out->weighted = false;
    out->with_replacement = true;

    if (W_ == 0 || t == 0) {
      out->pairs.clear();
      return true;
    }

    // ---------------------------------------------------------------------
    // Randomness protocol (SJS-HighDims §3.1.4):
    //   seed(i, phase_id, ctr) = Hash(master_seed, i, phase_id, ctr)
    //
    // IMPORTANT: do NOT use cfg.run.seed here.
    // The experiment runner offsets the per-repeat seed (cfg.run.seed + rep)
    // and passes phase-separated RNG streams into Count() and Sample().
    // Using cfg.run.seed would ignore that offset, causing repeats to share
    // within-event randomness.
    // ---------------------------------------------------------------------
    const u64 master_seed = rng->NextU64();
    static constexpr u64 kPhaseAssign = 1;  // slot->event assignment
    static constexpr u64 kPhaseSweep  = 2;  // within-event partner sampling
    static constexpr u64 kCtr0 = 0;

    // -----------------
    // Planning: assign t output slots to events by weights w_i/W.
    // -----------------
    detail::SlotPlanND plan;
    {
      auto _ = phases ? phases->Scoped("planning") : PhaseRecorder::ScopedPhase(nullptr, "");
      // Independent randomness stream for assignment.
      Rng rng_assign(DeriveSeed(master_seed, /*i=*/0, kPhaseAssign, kCtr0));
      if (!detail::BuildSlotPlanND(w_total_, t, &rng_assign, &plan)) {
        out->pairs.clear();
        return true;
      }
    }

    // -----------------
    // Pass 2: second sweep, sample only assigned blocks.
    // -----------------
    {
      auto _ = phases ? phases->Scoped("pass2_fill") : PhaseRecorder::ScopedPhase(nullptr, "");

      ctx_.ResetActive();
      const auto& events = ctx_.events();
      const auto& sid_of_pos = ctx_.start_id_of_event();

      std::vector<u32> sampled;
      sampled.reserve(1024);

      for (usize pos = 0; pos < events.size(); ++pos) {
        const auto& ev = events[pos];
        const u32 handle = static_cast<u32>(ev.index);

        if (ev.kind == join::EventKind::End) {
          ctx_.Delete(ev.side, handle);
          continue;
        }

        const i32 sid_i32 = sid_of_pos[pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 sid = static_cast<u32>(sid_i32);
        const usize sid_u = static_cast<usize>(sid);

        const u32 begin = plan.offset[sid_u];
        const u32 end = plan.offset[sid_u + 1];
        const u32 k = end - begin;

        const bool q_is_r = (ev.side == join::Side::R);
        const BoxT& q = ctx_.BoxOf(ev);

        if (k > 0) {
          // Per-event RNG stream for reproducibility.
          Rng ev_rng(DeriveSeed(master_seed, static_cast<u64>(sid), kPhaseSweep, kCtr0));

          sampled.clear();
          if (q_is_r) {
            ctx_.active_s().Sample(q, k, &ev_rng, &sampled);
          } else {
            ctx_.active_r().Sample(q, k, &ev_rng, &sampled);
          }

          if (sampled.size() != static_cast<usize>(k)) {
            if (err) *err = "RangeTreeSamplingBaseline::Sample: SAMPLE returned wrong size";
            ctx_.ResetActive();
            return false;
          }

          for (u32 i = 0; i < k; ++i) {
            const u32 slot = plan.slots[static_cast<usize>(begin + i)];
            const u32 oh = sampled[static_cast<usize>(i)];
            if (q_is_r) {
              out->pairs[static_cast<usize>(slot)] =
                  PairId{ds->R.GetId(static_cast<usize>(handle)), ds->S.GetId(static_cast<usize>(oh))};
            } else {
              out->pairs[static_cast<usize>(slot)] =
                  PairId{ds->R.GetId(static_cast<usize>(oh)), ds->S.GetId(static_cast<usize>(handle))};
            }
          }
        }

        // Insert q after sampling.
        ctx_.Insert(ev.side, handle);
      }

      ctx_.ResetActive();
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "RangeTreeSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RangeTreeReportJoinEnumeratorND<Dim, T>>(&ctx_);
  }

 private:
  bool built_{false};
  detail::RangeTreeContext<Dim, T> ctx_;

  std::vector<u64> w_total_;
  u64 W_{0};
  bool weights_valid_{false};
};

}  // namespace range_tree
}  // namespace baselines
}  // namespace sjs