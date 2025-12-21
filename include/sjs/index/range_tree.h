#pragma once
// sjs/index/range_tree.h
//
// RangeTree baseline for orthogonal range queries over points.
//
// In 2D, this implements a classic range tree / segment tree on x with
// y-sorted lists at nodes:
//
//   - Build: sort points by x, build balanced tree by index ranges.
//   - Each node stores its points' indices sorted by y.
//   - Query: decompose x-range into O(log n) nodes, for each node binary search
//            y-range, report points.
//
// Complexity (2D):
//   - Build: O(n log n) time, O(n log n) memory for y-lists.
//   - Query: O(log^2 n + out).
//
// For Dim==1, this degrades to a sorted array with binary search.
// For Dim>2, this implementation still builds an x-tree with y-sorted lists,
// and then checks remaining dimensions (2..Dim-1) by filtering candidates.
// This is not a full multi-dimensional range tree (which would be very memory
// intensive), but provides a reasonable, scalable baseline.
//
// Geometry convention:
//   - Queries are half-open: p in [lo,hi) per dimension.

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/geometry/point.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/predicates.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace sjs {
namespace index {

template <int Dim, class T = Scalar>
class RangeTree {
 public:
  static_assert(Dim >= 1, "RangeTree requires Dim >= 1");

  using PointT = Point<Dim, T>;
  using BoxT = Box<Dim, T>;

  struct BuildOptions {
    // Leaf segment size: larger reduces memory/overhead but increases leaf scan time.
    u32 leaf_size = 32;
  };

  RangeTree() = default;

  void Build(Span<const PointT> points) { Build(points, BuildOptions{}); }

  void Build(Span<const PointT> points, BuildOptions opt) {
    pts_ = points.data();
    n_ = points.size();
    opt_ = opt;

    nodes_.clear();
    x_order_.clear();
    y_pool_.clear();
    root_ = kNull;

    if (n_ == 0) return;
    SJS_ASSERT(opt_.leaf_size >= 1);

    if (n_ > static_cast<usize>(std::numeric_limits<u32>::max())) {
      // For simplicity we use u32 indices internally.
      SJS_ASSERT_MSG(false, "RangeTree: n too large for u32 internal indices");
    }

    x_order_.resize(n_);
    for (u32 i = 0; i < static_cast<u32>(n_); ++i) x_order_[i] = i;

    // Sort by x (dim0), tie-break by index for determinism.
    std::sort(x_order_.begin(), x_order_.end(), [&](u32 a, u32 b) {
      const T ax = pts_[a].v[0];
      const T bx = pts_[b].v[0];
      if (ax < bx) return true;
      if (bx < ax) return false;
      return a < b;
    });

    nodes_.reserve(n_ * 2);
    root_ = BuildRec(0, static_cast<u32>(n_));
  }

  bool Empty() const noexcept { return n_ == 0; }
  usize Size() const noexcept { return n_; }
  usize NodeCount() const noexcept { return nodes_.size(); }

  // Report all point indices inside `range`.
  template <class EmitFn>
  void RangeQuery(const BoxT& range, EmitFn&& emit) const {
    if (root_ == kNull) return;
    if (range.IsEmpty()) return;

    if constexpr (Dim == 1) {
      if (!RangeQuery1D(range, emit)) return;
    } else {
      // x-range -> index range [L,R) in x_order_
      const T xlo = range.lo.v[0];
      const T xhi = range.hi.v[0];
      const u32 L = LowerBoundX(xlo);
      const u32 R = LowerBoundX(xhi);
      if (L >= R) return;

      const T ylo = range.lo.v[1];
      const T yhi = range.hi.v[1];
      if (!(ylo < yhi)) return;

      if (!RangeQueryRec(root_, L, R, ylo, yhi, range, emit)) return;
    }
  }

  // Count points inside `range`.
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
    u32 left = kNull;
    u32 right = kNull;
    u32 l = 0;  // index range in x_order_ [l,r)
    u32 r = 0;
    u32 y_off = 0;
    u32 y_len = 0;

    bool IsLeaf() const noexcept { return left == kNull && right == kNull; }
  };

  const PointT* pts_{nullptr};
  usize n_{0};
  BuildOptions opt_{};

  std::vector<u32> x_order_;   // indices sorted by x
  std::vector<u32> y_pool_;    // concatenated y-sorted lists for nodes
  std::vector<Node> nodes_;
  u32 root_{kNull};

  // Compare point indices by y (dim1), tie by index.
  bool LessY(u32 a, u32 b) const noexcept {
    const T ay = pts_[a].v[1];
    const T by = pts_[b].v[1];
    if (ay < by) return true;
    if (by < ay) return false;
    return a < b;
  }

  // Build segment tree node for x_order_[l:r).
  u32 BuildRec(u32 l, u32 r) {
    SJS_DASSERT(l < r);
    Node node;
    node.l = l;
    node.r = r;

    const u32 my = static_cast<u32>(nodes_.size());
    nodes_.push_back(node);

    const u32 count = r - l;
    if (count <= opt_.leaf_size || Dim == 1) {
      // Leaf: y-list is just the points in this segment, sorted by y (if Dim>=2).
      std::vector<u32> tmp;
      tmp.reserve(count);
      for (u32 i = l; i < r; ++i) tmp.push_back(x_order_[i]);
      if constexpr (Dim >= 2) {
        std::sort(tmp.begin(), tmp.end(), [&](u32 a, u32 b) { return LessY(a, b); });
      }
      AppendList(my, tmp);
      return my;
    }

    const u32 mid = l + count / 2;
    const u32 left = BuildRec(l, mid);
    const u32 right = BuildRec(mid, r);
    nodes_[my].left = left;
    nodes_[my].right = right;

    // Merge children's y-lists.
    if constexpr (Dim >= 2) {
      const Node& A = nodes_[left];
      const Node& B = nodes_[right];
      const u32* a0 = y_pool_.data() + A.y_off;
      const u32* a1 = a0 + A.y_len;
      const u32* b0 = y_pool_.data() + B.y_off;
      const u32* b1 = b0 + B.y_len;

      std::vector<u32> merged;
      merged.reserve(static_cast<usize>(A.y_len) + static_cast<usize>(B.y_len));
      const u32* pa = a0;
      const u32* pb = b0;
      while (pa != a1 && pb != b1) {
        const u32 ia = *pa;
        const u32 ib = *pb;
        if (LessY(ia, ib)) {
          merged.push_back(ia);
          ++pa;
        } else {
          merged.push_back(ib);
          ++pb;
        }
      }
      while (pa != a1) { merged.push_back(*pa++); }
      while (pb != b1) { merged.push_back(*pb++); }

      AppendList(my, merged);
    } else {
      // Dim==1: not reached (handled above).
      AppendList(my, std::vector<u32>{});
    }

    return my;
  }

  void AppendList(u32 node_index, const std::vector<u32>& list) {
    Node& n = nodes_[node_index];
    n.y_off = static_cast<u32>(y_pool_.size());
    n.y_len = static_cast<u32>(list.size());
    y_pool_.insert(y_pool_.end(), list.begin(), list.end());
  }

  // Lower_bound in x_order_ for x >= value (dim0).
  u32 LowerBoundX(T value) const {
    auto it = std::lower_bound(
        x_order_.begin(), x_order_.end(), value,
        [&](u32 idx, const T& v) { return pts_[idx].v[0] < v; });
    return static_cast<u32>(it - x_order_.begin());
  }

  // Lower_bound within a node's y-list for y >= value.
  u32 LowerBoundY(const Node& node, T value) const {
    const u32* b = y_pool_.data() + node.y_off;
    const u32* e = b + node.y_len;
    auto it = std::lower_bound(
        b, e, value,
        [&](u32 idx, const T& v) { return pts_[idx].v[1] < v; });
    return static_cast<u32>(it - b);
  }

  // Filter remaining dimensions (d>=2) for Dim>2.
  bool PassOtherDims(u32 idx, const BoxT& range) const noexcept {
    if constexpr (Dim <= 2) {
      (void)idx;
      (void)range;
      return true;
    } else {
      for (int d = 2; d < Dim; ++d) {
        const T x = pts_[idx].v[d];
        if (x < range.lo.v[d]) return false;
        if (!(x < range.hi.v[d])) return false;
      }
      return true;
    }
  }

  template <class EmitFn>
bool RangeQueryRec(u32 node_index,
                   u32 ql, u32 qr,
                   T ylo, T yhi,
                   const BoxT& full_range,
                   EmitFn& emit) const {
  const Node& node = nodes_[node_index];
  if (qr <= node.l || node.r <= ql) return true;  // disjoint
  if (ql <= node.l && node.r <= qr) {
    // fully covered in x -> use y-list
    const u32 a = LowerBoundY(node, ylo);
    const u32 b = LowerBoundY(node, yhi);
    const u32* base = y_pool_.data() + node.y_off;
    for (u32 i = a; i < b; ++i) {
      const u32 idx = base[i];
      // Fast x check is implicit due to decomposition; y check due to bounds; check other dims.
      if (PassOtherDims(idx, full_range)) {
        if (!emit(static_cast<usize>(idx))) return false;
      }
    }
    return true;
  }

  if (node.IsLeaf()) {
    // Leaf with partial overlap: scan its x-range and check full range.
    for (u32 i = std::max(node.l, ql); i < std::min(node.r, qr); ++i) {
      const u32 idx = x_order_[i];
      if (ContainsHalfOpen<Dim, T>(full_range, pts_[idx])) {
        if (!emit(static_cast<usize>(idx))) return false;
      }
    }
    return true;
  }

  if (!RangeQueryRec(node.left, ql, qr, ylo, yhi, full_range, emit)) return false;
  if (!RangeQueryRec(node.right, ql, qr, ylo, yhi, full_range, emit)) return false;
  return true;
}

  template <class EmitFn>
bool RangeQuery1D(const BoxT& range, EmitFn& emit) const {
  const T xlo = range.lo.v[0];
  const T xhi = range.hi.v[0];
  if (!(xlo < xhi)) return true;
  const u32 L = LowerBoundX(xlo);
  const u32 R = LowerBoundX(xhi);
  for (u32 i = L; i < R; ++i) {
    const u32 idx = x_order_[i];
    // In 1D, x check is enough.
    if (!emit(static_cast<usize>(idx))) return false;
  }
  return true;
}
};

}  // namespace index
}  // namespace sjs
