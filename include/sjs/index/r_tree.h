#pragma once
// sjs/index/r_tree.h
//
// Packed R-tree for axis-aligned boxes (static bulk-loaded).
//
// Intended use:
//  - Baseline join: build RTree on S, then query each r in R.
//
// This implementation focuses on simplicity and determinism:
//  - Bulk-loading only (no dynamic insert/delete).
//  - Node fanout up to M (default 32).
//  - Uses STR-style packing in 2D; for Dim>2 uses axis-cycling sort-and-pack.
//
// Representation (packed):
//  - A Node stores a range [begin,end) into a global `children_` vector.
//  - Each child reference contains a bbox and a ref:
//      - if node.leaf=true: ref is object index.
//      - else: ref is child node index.
//
// Query:
//  - Depth-first traversal; prune by bbox intersection.
//
// Notes:
//  - Geometry uses half-open boxes, and intersection tests follow that convention.
//  - This is a baseline-quality R-tree meant for experiments, not a fully optimized
//    implementation with sophisticated splits/reinsertion.

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/predicates.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace sjs {
namespace index {

template <int Dim, class T = Scalar>
class RTree {
 public:
  using BoxT = Box<Dim, T>;

  struct BuildOptions {
    u32 max_children = 32;   // M
    bool use_str_2d = true;  // use STR packing when Dim==2
  };

  RTree() = default;

  void Build(Span<const BoxT> boxes) { Build(boxes, BuildOptions{}); }

  void Build(Span<const BoxT> boxes, BuildOptions opt) {
    boxes_ = boxes.data();
    n_ = boxes.size();
    opt_ = opt;

    nodes_.clear();
    children_.clear();
    root_ = kNull;
    height_ = 0;

    if (n_ == 0) return;
    SJS_ASSERT(opt_.max_children >= 2);

    // Initial refs correspond to objects.
    std::vector<Ref> level;
    level.reserve(n_);
    for (u32 i = 0; i < static_cast<u32>(n_); ++i) {
      level.push_back(Ref{boxes_[i], i});
    }

    bool leaf_level = true;
    u32 level_id = 0;
    while (true) {
      // Build one level of nodes from refs.
      std::vector<Ref> parent = BuildLevel(std::move(level), leaf_level, level_id);
      ++height_;
      if (parent.size() == 1) {
        // parent[0].ref is the root node index
        root_ = parent[0].ref;
        break;
      }
      leaf_level = false;
      level = std::move(parent);
      ++level_id;
    }
  }

  bool Empty() const noexcept { return n_ == 0; }
  usize Size() const noexcept { return n_; }
  usize NodeCount() const noexcept { return nodes_.size(); }
  u32 Height() const noexcept { return height_; }

  BoxT Bounds() const noexcept {
    if (root_ == kNull) return BoxT::Empty();
    return nodes_[root_].bbox;
  }

  template <class EmitFn>
  void QueryIntersections(const BoxT& q, EmitFn&& emit) const {
    if (root_ == kNull) return;
    if (q.IsEmpty()) return;

    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const u32 ni = stack.back();
      stack.pop_back();
      const Node& node = nodes_[ni];

      if (!IntersectsHalfOpen<Dim, T>(node.bbox, q)) continue;

      if (node.leaf) {
        for (u32 i = node.begin; i < node.end; ++i) {
          const Ref& r = children_[i];
          if (IntersectsHalfOpen<Dim, T>(r.bbox, q)) {
            if (!emit(static_cast<usize>(r.ref))) return;
          }
        }
      } else {
        for (u32 i = node.begin; i < node.end; ++i) {
          const Ref& r = children_[i];
          if (IntersectsHalfOpen<Dim, T>(r.bbox, q)) {
            stack.push_back(r.ref);
          }
        }
      }
    }
  }

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

  struct Ref {
    BoxT bbox;
    u32 ref;  // object index (leaf) or node index (internal)
  };

  struct Node {
    BoxT bbox = BoxT::Empty();
    u32 begin = 0;
    u32 end = 0;
    bool leaf = false;
  };

  const BoxT* boxes_{nullptr};
  usize n_{0};
  BuildOptions opt_{};

  std::vector<Node> nodes_;
  std::vector<Ref> children_;
  u32 root_{kNull};
  u32 height_{0};

  static usize CeilDiv(usize a, usize b) noexcept {
    return (a + b - 1) / b;
  }

  static T CenterCoord(const BoxT& b, int axis) noexcept {
    return (b.lo[axis] + b.hi[axis]) / T(2);
  }

  static BoxT UnionAll(Span<const Ref> refs) {
    BoxT b = BoxT::Empty();
    for (const auto& r : refs) b.ExpandToIncludeBox(r.bbox);
    return b;
  }

  // Build one tree level: pack `refs` into nodes with up to M children.
  std::vector<Ref> BuildLevel(std::vector<Ref>&& refs, bool leaf_level, u32 level_id) {
    const u32 M = opt_.max_children;
    if (refs.empty()) return {};

    // Reorder refs for packing.
    if (opt_.use_str_2d && Dim == 2) {
      STRSort2D(&refs, M);
    } else {
      AxisCycleSort(&refs, static_cast<int>(level_id % Dim));
    }

    std::vector<Ref> parent;
    parent.reserve(CeilDiv(refs.size(), M));

    usize pos = 0;
    while (pos < refs.size()) {
      const usize end = std::min(refs.size(), pos + static_cast<usize>(M));
      const u32 begin_off = static_cast<u32>(children_.size());

      // Append children for this node.
      for (usize i = pos; i < end; ++i) children_.push_back(refs[i]);

      const u32 end_off = static_cast<u32>(children_.size());
      Node node;
      node.begin = begin_off;
      node.end = end_off;
      node.leaf = leaf_level;
      node.bbox = UnionAll(Span<const Ref>(children_.data() + begin_off, end_off - begin_off));

      const u32 node_index = static_cast<u32>(nodes_.size());
      nodes_.push_back(node);

      parent.push_back(Ref{node.bbox, node_index});
      pos = end;
    }

    return parent;
  }

  // STR packing for 2D (Sort-Tile-Recursive).
  // Steps:
  //  1) sort by center x
  //  2) slice into groups of size S*M, where S = ceil(sqrt(num_nodes))
  //  3) within each slice sort by center y
  //  4) pack sequentially into nodes of size M
  static void STRSort2D(std::vector<Ref>* refs, u32 M) {
    SJS_ASSERT(refs != nullptr);
    auto& v = *refs;
    const usize n = v.size();
    if (n <= M) return;

    const usize num_nodes = CeilDiv(n, static_cast<usize>(M));
    const usize S = static_cast<usize>(std::ceil(std::sqrt(static_cast<double>(num_nodes))));
    const usize slice_cap = S * static_cast<usize>(M);

    std::sort(v.begin(), v.end(), [](const Ref& a, const Ref& b) {
      const T ax = CenterCoord(a.bbox, 0);
      const T bx = CenterCoord(b.bbox, 0);
      if (ax < bx) return true;
      if (bx < ax) return false;
      // tie-break for determinism
      const T ay = CenterCoord(a.bbox, 1);
      const T by = CenterCoord(b.bbox, 1);
      if (ay < by) return true;
      if (by < ay) return false;
      return a.ref < b.ref;
    });

    for (usize s0 = 0; s0 < n; s0 += slice_cap) {
      const usize s1 = std::min(n, s0 + slice_cap);
      std::sort(v.begin() + s0, v.begin() + s1, [](const Ref& a, const Ref& b) {
        const T ay = CenterCoord(a.bbox, 1);
        const T by = CenterCoord(b.bbox, 1);
        if (ay < by) return true;
        if (by < ay) return false;
        const T ax = CenterCoord(a.bbox, 0);
        const T bx = CenterCoord(b.bbox, 0);
        if (ax < bx) return true;
        if (bx < ax) return false;
        return a.ref < b.ref;
      });
    }
  }

  // Generic fallback: sort by center on a chosen axis.
  static void AxisCycleSort(std::vector<Ref>* refs, int axis) {
    SJS_ASSERT(refs != nullptr);
    auto& v = *refs;
    std::sort(v.begin(), v.end(), [axis](const Ref& a, const Ref& b) {
      const T aa = CenterCoord(a.bbox, axis);
      const T bb = CenterCoord(b.bbox, axis);
      if (aa < bb) return true;
      if (bb < aa) return false;
      return a.ref < b.ref;
    });
  }
};

}  // namespace index
}  // namespace sjs
