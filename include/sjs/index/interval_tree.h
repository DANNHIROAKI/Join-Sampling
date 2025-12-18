#pragma once
// sjs/index/interval_tree.h
//
// Dynamic 1D interval tree for half-open intervals [lo, hi).
//
// Primary intended use (2D join baseline):
//   - Sweep on x, maintain active set in an IntervalTree on y.
//   - When a START event arrives, query overlaps on y to get candidates,
//     then do full box-box test to confirm intersection.
//
// Implementation:
//   - Randomized treap keyed by (lo, handle) with subtree augmentation max_hi.
//   - Supports insert/erase/query in expected O(log n + out).
//   - Preallocated node array: each interval is addressed by a stable handle (usize).
//     This avoids per-operation heap allocations.
//
// Notes:
//   - Half-open overlap: [a0,a1) intersects [b0,b1) iff a0 < b1 && b0 < a1.
//   - If an interval is empty (lo >= hi), it is treated as non-insertable.
//
// Extensibility:
//   - Works for any axis; wrapper IntervalTreeOnAxis<Dim> is provided for Box<Dim>.
//
// Disclaimer:
//   - This is a baseline-quality data structure for experiments, not a fully
//     tuned production interval tree.

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/geometry/box.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace index {

namespace detail {

// SplitMix64 finalizer (deterministic pseudo-random).
inline u64 SplitMix64(u64 x) noexcept {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  x = x ^ (x >> 31);
  return x;
}

inline u32 PriorityFromHandle(u64 seed, usize handle) noexcept {
  const u64 h = SplitMix64(seed ^ static_cast<u64>(handle));
  return static_cast<u32>(h >> 32);  // upper 32 bits
}

template <class T>
constexpr bool Intersects1DHalfOpen(T a0, T a1, T b0, T b1) noexcept {
  return (a0 < b1) && (b0 < a1);
}

}  // namespace detail

class IntervalTree {
 public:
  struct Options {
    // Deterministic seed used to derive treap priorities from handles.
    u64 priority_seed = 0xC001D00DULL;

    // Defensive: if true, Insert() will check for duplicate insert and ignore.
    // If false, duplicate insert is an assertion failure in debug.
    bool ignore_duplicate_insert = false;
  };

  IntervalTree() = default;

  explicit IntervalTree(usize capacity) { Init(capacity, Options{}); }
  explicit IntervalTree(usize capacity, Options opt) { Init(capacity, opt); }

  void Init(usize capacity) { Init(capacity, Options{}); }

  void Init(usize capacity, Options opt) {
    opt_ = opt;
    capacity_ = capacity;
    nodes_.assign(capacity_, Node{});
    Reset();
  }

  void Reset() {
    root_ = kNull;
    size_ = 0;
    for (auto& n : nodes_) {
      n.left = kNull;
      n.right = kNull;
      n.in_tree = 0;
      n.lo = 0;
      n.hi = 0;
      n.max_hi = -std::numeric_limits<Scalar>::infinity();
      n.prio = 0;
    }
  }

  usize Capacity() const noexcept { return capacity_; }
  usize Size() const noexcept { return size_; }
  bool Empty() const noexcept { return size_ == 0; }

  bool Contains(usize handle) const noexcept {
    if (handle >= capacity_) return false;
    return nodes_[handle].in_tree != 0;
  }

  // Insert an interval for `handle`.
  // If lo>=hi, the interval is empty and will not be inserted (returns false).
  bool Insert(usize handle, Scalar lo, Scalar hi, std::string* err = nullptr) {
    if (handle >= capacity_) {
      if (err) *err = "IntervalTree::Insert: handle out of range";
      return false;
    }
    if (!(lo < hi)) {
      // empty interval under half-open semantics
      return false;
    }

    Node& n = nodes_[handle];

    if (n.in_tree) {
      if (opt_.ignore_duplicate_insert) return true;
      SJS_DASSERT_MSG(false, "IntervalTree::Insert: duplicate insert");
      if (err) *err = "IntervalTree::Insert: duplicate insert";
      return false;
    }

    n.lo = lo;
    n.hi = hi;
    n.max_hi = hi;
    n.left = kNull;
    n.right = kNull;
    n.prio = detail::PriorityFromHandle(opt_.priority_seed, handle);
    n.in_tree = 1;

    root_ = InsertRec(root_, static_cast<u32>(handle));
    ++size_;
    return true;
  }

  // Erase an interval by handle; returns true if removed, false if handle not present.
  bool Erase(usize handle) {
    if (handle >= capacity_) return false;
    Node& n = nodes_[handle];
    if (!n.in_tree) return false;
    root_ = EraseRec(root_, static_cast<u32>(handle));
    n.in_tree = 0;
    // Keep fields as-is; they may be useful for debugging.
    --size_;
    return true;
  }

  // Query all active intervals that overlap [qlo, qhi).
  // Calls emit(handle) for each matching handle.
  //
  // If you return false from emit, the query stops early.
  template <class EmitFn>
  void Query(Scalar qlo, Scalar qhi, EmitFn&& emit) const {
    if (!(qlo < qhi)) return;     // empty query
    if (root_ == kNull) return;
    QueryRec(root_, qlo, qhi, emit);
  }

 private:
  static constexpr u32 kNull = std::numeric_limits<u32>::max();

  struct Node {
    Scalar lo{0};
    Scalar hi{0};
    Scalar max_hi{-std::numeric_limits<Scalar>::infinity()};
    u32 left{kNull};
    u32 right{kNull};
    u32 prio{0};
    u8 in_tree{0};
    // padding intentionally omitted; struct is POD
  };

  Options opt_{};
  usize capacity_{0};
  usize size_{0};
  u32 root_{kNull};
  std::vector<Node> nodes_;

  // Ordering in BST: (lo, handle)
  bool Less(u32 a, u32 b) const noexcept {
    const Node& A = nodes_[a];
    const Node& B = nodes_[b];
    if (A.lo < B.lo) return true;
    if (B.lo < A.lo) return false;
    return a < b;
  }

  void Pull(u32 idx) noexcept {
    Node& n = nodes_[idx];
    Scalar m = n.hi;
    if (n.left != kNull) m = std::max(m, nodes_[n.left].max_hi);
    if (n.right != kNull) m = std::max(m, nodes_[n.right].max_hi);
    n.max_hi = m;
  }

  u32 RotateRight(u32 y) noexcept {
    const u32 x = nodes_[y].left;
    SJS_DASSERT(x != kNull);
    nodes_[y].left = nodes_[x].right;
    nodes_[x].right = y;
    Pull(y);
    Pull(x);
    return x;
  }

  u32 RotateLeft(u32 x) noexcept {
    const u32 y = nodes_[x].right;
    SJS_DASSERT(y != kNull);
    nodes_[x].right = nodes_[y].left;
    nodes_[y].left = x;
    Pull(x);
    Pull(y);
    return y;
  }

  u32 InsertRec(u32 root, u32 h) noexcept {
    if (root == kNull) {
      return h;
    }

    if (Less(h, root)) {
      nodes_[root].left = InsertRec(nodes_[root].left, h);
      if (nodes_[nodes_[root].left].prio < nodes_[root].prio) {
        root = RotateRight(root);
      }
    } else {
      nodes_[root].right = InsertRec(nodes_[root].right, h);
      if (nodes_[nodes_[root].right].prio < nodes_[root].prio) {
        root = RotateLeft(root);
      }
    }

    Pull(root);
    return root;
  }

  u32 EraseRec(u32 root, u32 h) noexcept {
    if (root == kNull) return kNull;

    if (root == h) {
      // Merge children by rotating the smaller-priority child up.
      if (nodes_[root].left == kNull) return nodes_[root].right;
      if (nodes_[root].right == kNull) return nodes_[root].left;

      // Rotate the child with smaller priority (min-heap).
      if (nodes_[nodes_[root].left].prio < nodes_[nodes_[root].right].prio) {
        root = RotateRight(root);
        nodes_[root].right = EraseRec(nodes_[root].right, h);
      } else {
        root = RotateLeft(root);
        nodes_[root].left = EraseRec(nodes_[root].left, h);
      }
      Pull(root);
      return root;
    }

    if (Less(h, root)) {
      nodes_[root].left = EraseRec(nodes_[root].left, h);
    } else {
      nodes_[root].right = EraseRec(nodes_[root].right, h);
    }
    Pull(root);
    return root;
  }

  template <class EmitFn>
  void QueryRec(u32 node, Scalar qlo, Scalar qhi, EmitFn& emit) const {
    if (node == kNull) return;

    // Prune: if qlo >= max_hi in this subtree, no interval has hi > qlo.
    if (!(qlo < nodes_[node].max_hi)) return;

    const Node& n = nodes_[node];

    // Explore left first.
    if (n.left != kNull) {
      QueryRec(n.left, qlo, qhi, emit);
    }

    // Check this node.
    if (detail::Intersects1DHalfOpen<Scalar>(n.lo, n.hi, qlo, qhi)) {
      if (!emit(static_cast<usize>(node))) return;
    }

    // Explore right: if this node's lo >= qhi, then all rights have lo >= qhi => cannot overlap.
    if (n.right != kNull && n.lo < qhi) {
      QueryRec(n.right, qlo, qhi, emit);
    }
  }
};

// Convenience wrapper: maintain intervals from Box<Dim> projected on one axis.
template <int Dim, class T = Scalar>
class IntervalTreeOnAxis {
 public:
  using BoxT = Box<Dim, T>;

  struct Options {
    int axis = 1;                       // default y-axis for 2D
    IntervalTree::Options tree_opt{};
  };

  IntervalTreeOnAxis() = default;

  IntervalTreeOnAxis(usize capacity) { Init(capacity, Options{}); }
  IntervalTreeOnAxis(usize capacity, Options opt) { Init(capacity, opt); }

  void Init(usize capacity) { Init(capacity, Options{}); }

  void Init(usize capacity, Options opt) {
    opt_ = opt;
    SJS_ASSERT(opt_.axis >= 0 && opt_.axis < Dim);
    tree_.Init(capacity, opt_.tree_opt);
  }

  void Reset() { tree_.Reset(); }

  usize Size() const noexcept { return tree_.Size(); }
  bool Contains(usize handle) const noexcept { return tree_.Contains(handle); }

  bool InsertBox(usize handle, const BoxT& b, std::string* err = nullptr) {
    const int a = opt_.axis;
    const Scalar lo = static_cast<Scalar>(b.lo[a]);
    const Scalar hi = static_cast<Scalar>(b.hi[a]);
    return tree_.Insert(handle, lo, hi, err);
  }

  bool Erase(usize handle) { return tree_.Erase(handle); }

  template <class EmitFn>
  void QueryBox(const BoxT& q, EmitFn&& emit) const {
    const int a = opt_.axis;
    const Scalar qlo = static_cast<Scalar>(q.lo[a]);
    const Scalar qhi = static_cast<Scalar>(q.hi[a]);
    tree_.Query(qlo, qhi, emit);
  }

 private:
  Options opt_{};
  IntervalTree tree_{};
};

}  // namespace index
}  // namespace sjs
