#pragma once
// sjs/index/aabb_tree.h
//
// Static AABB tree (BVH) for axis-aligned boxes.
//
// Intended use:
//  - Baseline spatial join: build an AABBTree on S, then for each r in R,
//    query candidates that intersect r.
//
// Design goals:
//  - Small, dependency-free implementation suitable for experiments.
//  - Generic in dimension (Dim), works for 2D now and higher-D later.
//  - Stable and deterministic (no RNG required).
//
// Data layout:
//  - The tree stores a permutation of object indices in `indices_`.
//  - Each node refers to a contiguous range [begin,end) in `indices_`.
//  - Internal nodes split by median along the widest axis of the node bbox.
//  - Leaves scan their range and call the user callback.
//
// Complexity:
//  - Build: O(n log n) expected (due to nth_element at each recursion).
//  - Query: O(log n + out) typical; worst-case O(n).

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/predicates.h"

#include <algorithm>
#include <array>
#include <limits>
#include <utility>
#include <vector>

namespace sjs {
namespace index {

template <int Dim, class T = Scalar>
class AABBTree {
 public:
  using BoxT = Box<Dim, T>;

  struct BuildOptions {
    u32 leaf_size = 8;  // max objects per leaf
  };

  AABBTree() = default;

  void Build(Span<const BoxT> boxes) { Build(boxes, BuildOptions{}); }

  // Build over externally-owned boxes. The caller must keep `boxes` alive.
  void Build(Span<const BoxT> boxes, BuildOptions opt) {
    boxes_ = boxes.data();
    n_ = boxes.size();
    opt_ = opt;
    nodes_.clear();
    indices_.clear();
    root_ = kNull;

    if (n_ == 0) return;

    indices_.resize(n_);
    for (usize i = 0; i < n_; ++i) indices_[i] = static_cast<u32>(i);

    nodes_.reserve(n_ * 2);
    root_ = BuildRec(0, static_cast<u32>(n_));
  }

  bool Empty() const noexcept { return n_ == 0; }
  usize Size() const noexcept { return n_; }
  usize NodeCount() const noexcept { return nodes_.size(); }

  BoxT Bounds() const noexcept {
    if (root_ == kNull) return BoxT::Empty();
    return nodes_[root_].bbox;
  }

  // Query boxes that intersect `q`. Calls emit(obj_index) for each candidate.
  // Stop early if emit returns false.
  template <class EmitFn>
  void QueryIntersections(const BoxT& q, EmitFn&& emit) const {
    if (root_ == kNull) return;
    if (q.IsEmpty()) return;

    // Iterative DFS to avoid recursion for deep trees.
    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const u32 ni = stack.back();
      stack.pop_back();
      const Node& node = nodes_[ni];
      if (!IntersectsHalfOpen<Dim, T>(node.bbox, q)) continue;

      if (node.IsLeaf()) {
        for (u32 k = node.begin; k < node.end; ++k) {
          const u32 obj = indices_[k];
          if (IntersectsHalfOpen<Dim, T>(boxes_[obj], q)) {
            if (!emit(static_cast<usize>(obj))) return;
          }
        }
      } else {
        stack.push_back(node.left);
        stack.push_back(node.right);
      }
    }
  }

  // Count intersections (no callback).
  u64 CountIntersections(const BoxT& q) const {
    u64 c = 0;
    QueryIntersections(q, [&](usize) {
      ++c;
      return true;
    });
    return c;
  }

 private:
  static constexpr u32 kNull = std::numeric_limits<u32>::max();

  struct Node {
    BoxT bbox = BoxT::Empty();
    u32 left = kNull;
    u32 right = kNull;
    u32 begin = 0;
    u32 end = 0;

    bool IsLeaf() const noexcept { return left == kNull && right == kNull; }
  };

  const BoxT* boxes_{nullptr};
  usize n_{0};
  BuildOptions opt_{};

  std::vector<Node> nodes_;
  std::vector<u32> indices_;
  u32 root_{kNull};

  // Compute bbox union of boxes in indices_[begin:end).
  BoxT ComputeBounds(u32 begin, u32 end) const {
    BoxT b = BoxT::Empty();
    for (u32 i = begin; i < end; ++i) {
      b.ExpandToIncludeBox(boxes_[indices_[i]]);
    }
    return b;
  }

  int ChooseSplitAxis(const BoxT& b) const noexcept {
    // Split along widest axis.
    int axis = 0;
    T best = b.Width(0);
    for (int d = 1; d < Dim; ++d) {
      const T w = b.Width(d);
      if (w > best) {
        best = w;
        axis = d;
      }
    }
    return axis;
  }

  u32 BuildRec(u32 begin, u32 end) {
    SJS_DASSERT(begin < end);
    const u32 count = end - begin;

    Node node;
    node.begin = begin;
    node.end = end;
    node.bbox = ComputeBounds(begin, end);

    const u32 my_index = static_cast<u32>(nodes_.size());
    nodes_.push_back(node);

    if (count <= opt_.leaf_size) {
      // Leaf.
      return my_index;
    }

    const int axis = ChooseSplitAxis(node.bbox);
    const u32 mid = begin + count / 2;

    // Partition indices by center along chosen axis.
    std::nth_element(indices_.begin() + begin,
                     indices_.begin() + mid,
                     indices_.begin() + end,
                    [&](u32 ia, u32 ib) {
                      const T ca = boxes_[ia].Center()[axis];
                      const T cb = boxes_[ib].Center()[axis];
                       if (ca < cb) return true;
                       if (cb < ca) return false;
                       return ia < ib;
                     });

    // Handle degenerate split (shouldn't happen often).
    u32 left_begin = begin;
    u32 left_end = mid;
    u32 right_begin = mid;
    u32 right_end = end;
    if (left_begin == left_end || right_begin == right_end) {
      left_end = begin + count / 2;
      right_begin = left_end;
      // If still degenerate, force leaf.
      if (left_begin == left_end || right_begin == right_end) {
        return my_index;
      }
    }

    const u32 left = BuildRec(left_begin, left_end);
    const u32 right = BuildRec(right_begin, right_end);

    nodes_[my_index].left = left;
    nodes_[my_index].right = right;
    // bbox already computed; keep.
    return my_index;
  }
};

}  // namespace index
}  // namespace sjs
