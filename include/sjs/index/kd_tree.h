#pragma once
// sjs/index/kd_tree.h
//
// Static KD-tree for points with orthogonal range queries.
//
// Intended use:
//  - Baselines that embed rectangles into points in 2*Dim dimensions,
//    then query an orthogonal range (KD-tree or RangeTree style).
//
// Design goals:
//  - Deterministic build.
//  - Works for Dim>=1, but note that KD-tree performance degrades as Dim grows.
//  - Supports streaming range reporting and fast counting with bbox pruning.
//
// Geometry conventions:
//  - Range queries are half-open boxes: p in [lo,hi) for each dimension.
//
// Implementation:
//  - Balanced KD-tree by recursive median split.
//  - Node stores bounding box (half-open) of its points to allow pruning.
//  - Leaves contain up to leaf_size points scanned linearly.
//
// Complexity:
//  - Build: O(n log n) expected.
//  - Query: typical O(log n + out) in low dims; can degrade in high dims.

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/geometry/point.h"
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
class KDTree {
 public:
  using PointT = Point<Dim, T>;
  using BoxT = Box<Dim, T>;

  struct BuildOptions {
    u32 leaf_size = 16;
  };

  KDTree() = default;

  void Build(Span<const PointT> points) { Build(points, BuildOptions{}); }

  void Build(Span<const PointT> points, BuildOptions opt) {
    pts_ = points.data();
    n_ = points.size();
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

  template <class EmitFn>
  void RangeQuery(const BoxT& range, EmitFn&& emit) const {
    if (root_ == kNull) return;
    if (range.IsEmpty()) return;

    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const u32 ni = stack.back();
      stack.pop_back();
      const Node& node = nodes_[ni];

      if (!IntersectsHalfOpen<Dim, T>(node.bbox, range)) continue;

      if (node.IsLeaf()) {
        for (u32 k = node.begin; k < node.end; ++k) {
          const u32 idx = indices_[k];
          if (ContainsHalfOpen<Dim, T>(range, pts_[idx])) {
            if (!emit(static_cast<usize>(idx))) return;
          }
        }
      } else {
        stack.push_back(node.left);
        stack.push_back(node.right);
      }
    }
  }

  u64 RangeCount(const BoxT& range) const {
    u64 c = 0;
    RangeQuery(range, [&](usize) {
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
    u8 axis = 0;  // split axis for debugging (not required for query)

    bool IsLeaf() const noexcept { return left == kNull && right == kNull; }
  };

  const PointT* pts_{nullptr};
  usize n_{0};
  BuildOptions opt_{};

  std::vector<Node> nodes_;
  std::vector<u32> indices_;
  u32 root_{kNull};

  BoxT ComputeBounds(u32 begin, u32 end) const {
    BoxT b = BoxT::Empty();
    for (u32 i = begin; i < end; ++i) {
      b.ExpandToIncludePoint(pts_[indices_[i]]);
    }
    return b;
  }

  int ChooseSplitAxis(const BoxT& bbox) const noexcept {
    // Choose the axis with the largest width in the bbox.
    int axis = 0;
    T best = bbox.Width(0);
    for (int d = 1; d < Dim; ++d) {
      const T w = bbox.Width(d);
      if (w > best) {
        best = w;
        axis = d;
      }
    }
    return axis;
  }

  u32 BuildRec(u32 begin, u32 end) {
    const u32 count = end - begin;
    SJS_DASSERT(count > 0);

    Node node;
    node.begin = begin;
    node.end = end;
    node.bbox = ComputeBounds(begin, end);

    const u32 my = static_cast<u32>(nodes_.size());
    nodes_.push_back(node);

    if (count <= opt_.leaf_size) {
      return my;
    }

    const int axis = ChooseSplitAxis(node.bbox);
    nodes_[my].axis = static_cast<u8>(axis);

    const u32 mid = begin + count / 2;

    std::nth_element(indices_.begin() + begin,
                     indices_.begin() + mid,
                     indices_.begin() + end,
                    [&](u32 ia, u32 ib) {
                      const T a = pts_[ia][axis];
                      const T b = pts_[ib][axis];
                       if (a < b) return true;
                       if (b < a) return false;
                       return ia < ib;
                     });

    // Degenerate splits: if all points have identical coord on axis, just split by index.
    u32 left_begin = begin;
    u32 left_end = mid;
    u32 right_begin = mid;
    u32 right_end = end;

    if (left_begin == left_end || right_begin == right_end) {
      left_end = begin + count / 2;
      right_begin = left_end;
      if (left_begin == left_end || right_begin == right_end) {
        return my;  // leaf fallback
      }
    }

    const u32 left = BuildRec(left_begin, left_end);
    const u32 right = BuildRec(right_begin, right_end);

    nodes_[my].left = left;
    nodes_[my].right = right;
    return my;
  }
};

}  // namespace index
}  // namespace sjs
