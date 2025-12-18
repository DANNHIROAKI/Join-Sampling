#pragma once
// sjs/baselines/aabb/sampling.h
//
// Plane Sweep + Dynamic AABB-Tree baseline (Variant::Sampling).
//
// This implements the algorithm described in "AABB-Tree Baseline v2.0.md":
//   - Sweep on axis 0 (x1).
//   - Maintain two *dynamic* AABB trees over the projection to axes 1..Dim-1.
//   - Phase 1 (first sweep): for each START event e, compute w_e by
//       CountIntersect(q*) on the opposite active tree; accumulate W=|J|.
//   - Phase 2: build alias table on {w_e} and assign t sample slots to events.
//   - Phase 3 (second sweep): for each START event with t_e slots, draw t_e
//       i.i.d. uniform samples from the local intersection set K_e using
//       BuildGuide + weighted descent, and fill the output slots.
//
// Notes:
//   - Correct for half-open boxes [lo,hi) on each axis.
//   - Complexity depends on BVH quality (worst-case may degrade).
//   - This header also provides a deterministic enumerator based on the same
//     sweep + AABB-tree machinery (used by other variants and runners).

#include "sjs/baselines/baseline_api.h"
#include "sjs/core/assert.h"
#include "sjs/join/sweep_events.h"
#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace aabb {

namespace detail {

// --------------------------
// Projection helpers
// --------------------------

// Project a Dim-dimensional box to (Dim-1) dimensions by dropping axis 0.
template <int Dim, class T>
inline Box<Dim - 1, T> ProjectDropFirst(const Box<Dim, T>& b) noexcept {
  static_assert(Dim >= 2, "ProjectDropFirst requires Dim >= 2");
  Box<Dim - 1, T> out;
  for (int i = 1; i < Dim; ++i) {
    out.lo.v[i - 1] = b.lo.v[i];
    out.hi.v[i - 1] = b.hi.v[i];
  }
  return out;
}

// --------------------------
// DynamicAABBTree (dynamic BVH)
// --------------------------
// A small, self-contained dynamic BVH implementation inspired by Box2D's
// dynamic tree. It supports insert/delete by handle and query-time operations
// needed by the AABB baseline:
//   - CountIntersect(Q)
//   - ReportIntersect(Q)
//   - SampleIntersect(Q, k): exact uniform sampling among intersecting leaves
//     via guide construction.

template <int Dim, class T>
class DynamicAABBTree {
 public:
  static_assert(Dim >= 1, "DynamicAABBTree requires Dim >= 1");
  using BoxT = Box<Dim, T>;

  struct QueryStats {
    u64 node_tests = 0;  // node AABB vs query tests
    u64 leaf_tests = 0;  // leaf AABB vs query tests

    void Reset() {
      node_tests = 0;
      leaf_tests = 0;
    }
  };

  DynamicAABBTree() = default;

  void Clear() {
    root_ = kNull;
    nodes_.clear();
  }

  void ReserveLeaves(usize leaves) {
    // A full binary tree with L leaves has at most 2L-1 nodes.
    const usize want = (leaves > 0) ? (2 * leaves + 8) : 0;
    if (want > nodes_.capacity()) nodes_.reserve(want);
  }

  bool Empty() const noexcept { return root_ == kNull; }

  u32 LeafCount() const noexcept {
    if (root_ == kNull) return 0;
    return nodes_[static_cast<usize>(root_)].size;
  }

  // Insert a leaf with object index `obj` and its AABB.
  // Returns a handle (node index) that can later be used for deletion.
  int Insert(u32 obj, const BoxT& aabb) {
    const int leaf = NewLeaf(obj, aabb);
    InsertLeaf(leaf);
    return leaf;
  }

  // Remove a leaf by handle (as returned by Insert).
  void Remove(int leaf) {
    if (leaf == kNull) return;
    RemoveLeaf(leaf);
  }

  // Count leaves whose AABB intersects query Q (half-open semantics).
  u64 CountIntersect(const BoxT& q, QueryStats* st = nullptr) const {
    if (root_ == kNull) return 0;
    if (q.IsEmpty()) return 0;

    u64 cnt = 0;
    std::vector<int> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const int idx = stack.back();
      stack.pop_back();
      const Node& n = nodes_[static_cast<usize>(idx)];
      if (st) ++st->node_tests;

      if (!IntersectsHalfOpen(n.aabb, q)) continue;

      // Full accept: if the subtree conservative AABB is contained in Q,
      // every leaf true AABB in this subtree is also contained in Q, hence intersects Q.
      if (ContainedHalfOpen(n.aabb, q)) {
        cnt += static_cast<u64>(n.size);
        continue;
      }

      if (n.IsLeaf()) {
        if (st) ++st->leaf_tests;
        if (IntersectsHalfOpen(n.true_aabb, q)) ++cnt;
      } else {
        stack.push_back(n.left);
        stack.push_back(n.right);
      }
    }
    return cnt;
  }

  // Report all intersecting leaves (returns their stored object indices).
  void ReportIntersect(const BoxT& q, std::vector<u32>* out, QueryStats* st = nullptr) const {
    if (!out) return;
    if (root_ == kNull) return;
    if (q.IsEmpty()) return;

    std::vector<int> stack;
    stack.reserve(64);
    stack.push_back(root_);
    while (!stack.empty()) {
      const int idx = stack.back();
      stack.pop_back();
      const Node& n = nodes_[static_cast<usize>(idx)];
      if (st) ++st->node_tests;
      if (!IntersectsHalfOpen(n.aabb, q)) continue;

      if (n.IsLeaf()) {
        if (st) ++st->leaf_tests;
        // Leaf must be tested against the exact (true) half-open box.
        if (IntersectsHalfOpen(n.true_aabb, q)) out->push_back(n.obj);
      } else {
        stack.push_back(n.left);
        stack.push_back(n.right);
      }
    }
  }

  // Sample k i.i.d. leaves uniformly among those intersecting Q.
  // Returns false if there are no intersecting leaves (i.e., |K(Q)|==0).
  bool SampleIntersect(const BoxT& q,
                       u32 k,
                       Rng* rng,
                       std::vector<u32>* out,
                       QueryStats* st = nullptr) const {
    if (!rng || !out) return false;
    out->clear();
    if (k == 0) return true;
    if (root_ == kNull) return false;
    if (q.IsEmpty()) return false;

    guide_.clear();
    int g_root = BuildGuide(q, st);
    if (g_root < 0) return false;
    const u64 total = guide_[static_cast<usize>(g_root)].cnt;
    if (total == 0) return false;

    out->reserve(static_cast<usize>(k));
    for (u32 i = 0; i < k; ++i) {
      const u32 obj = SampleFromGuide(g_root, rng);
      out->push_back(obj);
    }
    return true;
  }

 private:
  static constexpr int kNull = -1;

  struct Node {
    // Conservative bounding box used for internal-node pruning and full-accept.
    // For leaves this may be a "fat AABB"; for exact semantics always use true_aabb.
    BoxT aabb;

    // Exact leaf AABB (half-open). Meaningful for leaves; for internal nodes may equal aabb.
    BoxT true_aabb;

    int parent = kNull;
    int left = kNull;
    int right = kNull;

    // Height (leaf=0). Used for balancing.
    int height = 0;

    // Number of leaves in subtree.
    u32 size = 1;

    // User payload (object index). Meaningful only for leaves.
    u32 obj = 0;

    bool IsLeaf() const noexcept { return left == kNull; }
  };

  // Guide nodes (per query, stored in mutable arena to reduce allocations).
  struct GuideNode {
    enum class Type : u8 {
      Full = 0,    // subtree fully contained in query: sample uniformly from subtree
      Leaf = 1,    // single intersecting leaf
      Partial = 2  // partially intersecting internal node: recurse into guide children
    };

    Type type = Type::Leaf;
    u64 cnt = 0;          // number of intersecting leaves in this guide subtree
    int tree_node = kNull;  // for Full: index into BVH nodes_
    u32 obj = 0;            // for Leaf: object index
    int child0 = kNull;     // for Partial
    int child1 = kNull;     // for Partial (may be kNull if only one child intersects)
  };

  // BVH storage.
  std::vector<Node> nodes_;
  int root_ = kNull;

  // Mutable guide arena reused between SampleIntersect calls.
  mutable std::vector<GuideNode> guide_;

  static BoxT UnionBox(const BoxT& a, const BoxT& b) {
    BoxT u = a;
    u.ExpandToIncludeBox(b);
    return u;
  }

  // Half-open intersection test: [lo, hi) semantics in every dimension.
  static bool IntersectsHalfOpen(const BoxT& a, const BoxT& b) noexcept {
    for (int i = 0; i < Dim; ++i) {
      const T lo = (a.lo.v[i] > b.lo.v[i]) ? a.lo.v[i] : b.lo.v[i];
      const T hi = (a.hi.v[i] < b.hi.v[i]) ? a.hi.v[i] : b.hi.v[i];
      if (!(lo < hi)) return false;
    }
    return true;
  }

  // Return true iff inner ⊆ outer (both treated as half-open boxes).
  static bool ContainedHalfOpen(const BoxT& inner, const BoxT& outer) noexcept {
    for (int i = 0; i < Dim; ++i) {
      if (!(inner.lo.v[i] >= outer.lo.v[i])) return false;
      if (!(inner.hi.v[i] <= outer.hi.v[i])) return false;
    }
    return true;
  }

  // A cheap, monotone "perimeter-like" measure used by the insertion heuristic.
  // This is more overflow-resistant than full volume in higher dimensions.
  static double Measure(const BoxT& b) noexcept {
    double s = 0.0;
    for (int i = 0; i < Dim; ++i) {
      s += static_cast<double>(b.hi.v[i] - b.lo.v[i]);
    }
    return s;
  }

  int NewLeaf(u32 obj, const BoxT& aabb) {
    Node n;
    n.aabb = aabb;
    n.true_aabb = aabb;
    n.parent = kNull;
    n.left = kNull;
    n.right = kNull;
    n.height = 0;
    n.size = 1;
    n.obj = obj;
    nodes_.push_back(n);
    return static_cast<int>(nodes_.size() - 1);
  }

  int NewInternal() {
    Node n;
    n.parent = kNull;
    n.left = kNull;
    n.right = kNull;
    n.height = 1;
    n.size = 0;
    n.obj = 0;
    nodes_.push_back(n);
    return static_cast<int>(nodes_.size() - 1);
  }

  void InsertLeaf(int leaf) {
    if (root_ == kNull) {
      root_ = leaf;
      nodes_[static_cast<usize>(root_)].parent = kNull;
      return;
    }

    // Find the best sibling for this leaf.
    const BoxT leafAABB = nodes_[static_cast<usize>(leaf)].aabb;
    int sibling = FindBestSibling(leafAABB);

    // Create a new parent.
    const int oldParent = nodes_[static_cast<usize>(sibling)].parent;
    const int newParent = NewInternal();
    Node& p = nodes_[static_cast<usize>(newParent)];
    p.parent = oldParent;
    p.left = sibling;
    p.right = leaf;
    p.aabb = UnionBox(nodes_[static_cast<usize>(sibling)].aabb, leafAABB);
    p.height = nodes_[static_cast<usize>(sibling)].height + 1;
    p.size = nodes_[static_cast<usize>(sibling)].size + 1;

    nodes_[static_cast<usize>(sibling)].parent = newParent;
    nodes_[static_cast<usize>(leaf)].parent = newParent;

    if (oldParent == kNull) {
      root_ = newParent;
    } else {
      Node& op = nodes_[static_cast<usize>(oldParent)];
      if (op.left == sibling) {
        op.left = newParent;
      } else {
        op.right = newParent;
      }
    }

    // Walk back up the tree fixing AABBs/heights/sizes and balancing.
    FixUpwards(nodes_[static_cast<usize>(leaf)].parent);
  }

  void RemoveLeaf(int leaf) {
    if (leaf == root_) {
      root_ = kNull;
      return;
    }

    const int parent = nodes_[static_cast<usize>(leaf)].parent;
    const int grandParent = nodes_[static_cast<usize>(parent)].parent;

    const int sibling = (nodes_[static_cast<usize>(parent)].left == leaf)
                            ? nodes_[static_cast<usize>(parent)].right
                            : nodes_[static_cast<usize>(parent)].left;

    if (grandParent != kNull) {
      // Destroy parent and connect sibling to grandParent.
      Node& gp = nodes_[static_cast<usize>(grandParent)];
      if (gp.left == parent) {
        gp.left = sibling;
      } else {
        gp.right = sibling;
      }
      nodes_[static_cast<usize>(sibling)].parent = grandParent;

      // Fix the tree upwards.
      FixUpwards(grandParent);
    } else {
      // Parent was the root.
      root_ = sibling;
      nodes_[static_cast<usize>(sibling)].parent = kNull;
    }

    // Optional: detach removed nodes (debug friendliness).
    nodes_[static_cast<usize>(parent)].parent = kNull;
    nodes_[static_cast<usize>(parent)].left = kNull;
    nodes_[static_cast<usize>(parent)].right = kNull;
    nodes_[static_cast<usize>(leaf)].parent = kNull;
  }

  int FindBestSibling(const BoxT& leafAABB) const {
    int index = root_;
    while (!nodes_[static_cast<usize>(index)].IsLeaf()) {
      const Node& n = nodes_[static_cast<usize>(index)];
      const int left = n.left;
      const int right = n.right;

      const double area = Measure(n.aabb);
      const BoxT combined = UnionBox(n.aabb, leafAABB);
      const double combinedArea = Measure(combined);

      // Cost of creating a new parent for this node and the new leaf.
      const double costParent = 2.0 * combinedArea;
      // Minimum cost of pushing the leaf down the tree.
      const double inheritanceCost = 2.0 * (combinedArea - area);

      // Cost for descending into left child.
      double costLeft;
      {
        const Node& L = nodes_[static_cast<usize>(left)];
        const BoxT aabb = UnionBox(L.aabb, leafAABB);
        const double newArea = Measure(aabb);
        if (L.IsLeaf()) {
          costLeft = newArea + inheritanceCost;
        } else {
          const double oldArea = Measure(L.aabb);
          costLeft = (newArea - oldArea) + inheritanceCost;
        }
      }

      // Cost for descending into right child.
      double costRight;
      {
        const Node& R = nodes_[static_cast<usize>(right)];
        const BoxT aabb = UnionBox(R.aabb, leafAABB);
        const double newArea = Measure(aabb);
        if (R.IsLeaf()) {
          costRight = newArea + inheritanceCost;
        } else {
          const double oldArea = Measure(R.aabb);
          costRight = (newArea - oldArea) + inheritanceCost;
        }
      }

      // Descend according to the smallest cost.
      if (costParent < costLeft && costParent < costRight) break;
      if (costLeft < costRight) {
        index = left;
      } else if (costRight < costLeft) {
        index = right;
      } else {
        // Deterministic tie-break.
        index = left;
      }
    }
    return index;
  }

  void FixUpwards(int start) {
    int index = start;
    while (index != kNull) {
      index = Balance(index);
      Node& n = nodes_[static_cast<usize>(index)];
      if (!n.IsLeaf()) {
        const Node& L = nodes_[static_cast<usize>(n.left)];
        const Node& R = nodes_[static_cast<usize>(n.right)];
        n.size = L.size + R.size;
        n.height = 1 + std::max(L.height, R.height);
        n.aabb = UnionBox(L.aabb, R.aabb);
      } else {
        n.size = 1;
        n.height = 0;
      }
      index = n.parent;
    }
  }

  // Rotate and balance the subtree rooted at iA.
  // Returns the new root index of the subtree.
  int Balance(int iA) {
    Node& A = nodes_[static_cast<usize>(iA)];
    if (A.IsLeaf() || A.height < 2) return iA;

    const int iB = A.left;
    const int iC = A.right;
    Node& B = nodes_[static_cast<usize>(iB)];
    Node& C = nodes_[static_cast<usize>(iC)];

    const int balance = C.height - B.height;

    // Rotate C up.
    if (balance > 1) {
      const int iF = C.left;
      const int iG = C.right;
      Node& F = nodes_[static_cast<usize>(iF)];
      Node& G = nodes_[static_cast<usize>(iG)];

      // Swap A and C.
      C.left = iA;
      const int oldParent = A.parent;
      C.parent = oldParent;
      A.parent = iC;

      // A's old parent should point to C.
      if (C.parent != kNull) {
        Node& P = nodes_[static_cast<usize>(C.parent)];
        if (P.left == iA) {
          P.left = iC;
        } else {
          P.right = iC;
        }
      } else {
        root_ = iC;
      }

      // Rotate.
      if (F.height > G.height) {
        C.right = iF;
        A.right = iG;
        G.parent = iA;

        // Update A.
        A.aabb = UnionBox(B.aabb, G.aabb);
        A.size = B.size + G.size;
        A.height = 1 + std::max(B.height, G.height);

        // Update C.
        C.aabb = UnionBox(A.aabb, F.aabb);
        C.size = A.size + F.size;
        C.height = 1 + std::max(A.height, F.height);
      } else {
        C.right = iG;
        A.right = iF;
        F.parent = iA;

        // Update A.
        A.aabb = UnionBox(B.aabb, F.aabb);
        A.size = B.size + F.size;
        A.height = 1 + std::max(B.height, F.height);

        // Update C.
        C.aabb = UnionBox(A.aabb, G.aabb);
        C.size = A.size + G.size;
        C.height = 1 + std::max(A.height, G.height);
      }

      return iC;
    }

    // Rotate B up.
    if (balance < -1) {
      const int iD = B.left;
      const int iE = B.right;
      Node& D = nodes_[static_cast<usize>(iD)];
      Node& E = nodes_[static_cast<usize>(iE)];

      // Swap A and B.
      B.left = iA;
      const int oldParent = A.parent;
      B.parent = oldParent;
      A.parent = iB;

      // A's old parent should point to B.
      if (B.parent != kNull) {
        Node& P = nodes_[static_cast<usize>(B.parent)];
        if (P.left == iA) {
          P.left = iB;
        } else {
          P.right = iB;
        }
      } else {
        root_ = iB;
      }

      // Rotate.
      if (D.height > E.height) {
        B.right = iD;
        A.left = iE;
        E.parent = iA;

        // Update A.
        A.aabb = UnionBox(C.aabb, E.aabb);
        A.size = C.size + E.size;
        A.height = 1 + std::max(C.height, E.height);

        // Update B.
        B.aabb = UnionBox(A.aabb, D.aabb);
        B.size = A.size + D.size;
        B.height = 1 + std::max(A.height, D.height);
      } else {
        B.right = iE;
        A.left = iD;
        D.parent = iA;

        // Update A.
        A.aabb = UnionBox(C.aabb, D.aabb);
        A.size = C.size + D.size;
        A.height = 1 + std::max(C.height, D.height);

        // Update B.
        B.aabb = UnionBox(A.aabb, E.aabb);
        B.size = A.size + E.size;
        B.height = 1 + std::max(A.height, E.height);
      }

      return iB;
    }

    return iA;
  }

  // --------------------------
  // Guide construction & sampling
  // --------------------------

  // Build a guide for query q.
  // Returns the guide root index, or -1 if the intersection set is empty.
  int BuildGuide(const BoxT& q, QueryStats* st) const {
    guide_.clear();
    if (root_ == kNull) return -1;

    struct Frame {
      int node = kNull;
      u8 state = 0;   // 0=enter, 1=left done, 2=right done
      int left_g = kNull;
      int right_g = kNull;
    };

    std::vector<Frame> stack;
    stack.reserve(64);
    stack.push_back(Frame{root_, 0, kNull, kNull});

    int root_g = kNull;
    bool done = false;

    auto finish = [&](int res_g) {
      stack.pop_back();
      if (stack.empty()) {
        root_g = res_g;
        done = true;
        return;
      }
      Frame& p = stack.back();
      if (p.state == 1) {
        p.left_g = res_g;
      } else {
        // state==2
        p.right_g = res_g;
      }
    };

    while (!done) {
      Frame& f = stack.back();
      const Node& n = nodes_[static_cast<usize>(f.node)];

      if (f.state == 0) {
        if (st) ++st->node_tests;
        if (!IntersectsHalfOpen(n.aabb, q)) {
          finish(kNull);
          continue;
        }

        if (ContainedHalfOpen(n.aabb, q)) {
          GuideNode g;
          g.type = GuideNode::Type::Full;
          g.cnt = static_cast<u64>(n.size);
          g.tree_node = f.node;
          const int idx = static_cast<int>(guide_.size());
          guide_.push_back(g);
          finish(idx);
          continue;
        }

        if (n.IsLeaf()) {
          if (st) ++st->leaf_tests;
          // Leaf must be tested against the exact (true) half-open box.
          if (!IntersectsHalfOpen(n.true_aabb, q)) {
            finish(kNull);
            continue;
          }
          GuideNode g;
          g.type = GuideNode::Type::Leaf;
          g.cnt = 1;
          g.obj = n.obj;
          const int idx = static_cast<int>(guide_.size());
          guide_.push_back(g);
          finish(idx);
          continue;
        }

        // Internal: descend left.
        f.state = 1;
        stack.push_back(Frame{n.left, 0, kNull, kNull});
        continue;
      }

      if (f.state == 1) {
        // Left done -> descend right.
        f.state = 2;
        stack.push_back(Frame{n.right, 0, kNull, kNull});
        continue;
      }

      // Both children done: build partial/forward.
      const int lg = f.left_g;
      const int rg = f.right_g;
      if (lg == kNull && rg == kNull) {
        finish(kNull);
        continue;
      }
      if (rg == kNull) {
        finish(lg);
        continue;
      }
      if (lg == kNull) {
        finish(rg);
        continue;
      }

      GuideNode g;
      g.type = GuideNode::Type::Partial;
      g.child0 = lg;
      g.child1 = rg;
      g.cnt = guide_[static_cast<usize>(lg)].cnt + guide_[static_cast<usize>(rg)].cnt;
      const int idx = static_cast<int>(guide_.size());
      guide_.push_back(g);
      finish(idx);
    }

    return root_g;
  }

  u32 SampleSubtreeUniform(int node, Rng* rng) const {
    int cur = node;
    while (!nodes_[static_cast<usize>(cur)].IsLeaf()) {
      const Node& n = nodes_[static_cast<usize>(cur)];
      const Node& L = nodes_[static_cast<usize>(n.left)];
      const Node& R = nodes_[static_cast<usize>(n.right)];
      const u64 sl = static_cast<u64>(L.size);
      const u64 sr = static_cast<u64>(R.size);
      const u64 pick = rng->UniformU64(sl + sr);
      cur = (pick < sl) ? n.left : n.right;
    }
    return nodes_[static_cast<usize>(cur)].obj;
  }

  u32 SampleFromGuide(int g_root, Rng* rng) const {
    int g = g_root;
    while (true) {
      const GuideNode& gn = guide_[static_cast<usize>(g)];
      switch (gn.type) {
        case GuideNode::Type::Full: {
          return SampleSubtreeUniform(gn.tree_node, rng);
        }
        case GuideNode::Type::Leaf: {
          return gn.obj;
        }
        case GuideNode::Type::Partial: {
          const int c0 = gn.child0;
          const int c1 = gn.child1;
          if (c1 == kNull) {
            g = c0;
            continue;
          }
          const u64 w0 = guide_[static_cast<usize>(c0)].cnt;
          const u64 w1 = guide_[static_cast<usize>(c1)].cnt;
          const u64 pick = rng->UniformU64(w0 + w1);
          g = (pick < w0) ? c0 : c1;
          continue;
        }
      }
    }
  }
};

// --------------------------
// Deterministic enumerator: Plane Sweep + AABB Tree
// --------------------------

template <int Dim, class T>
class AABBJoinEnumerator final : public baselines::IJoinEnumerator {
 public:
  static_assert(Dim >= 2, "AABBJoinEnumerator requires Dim >= 2");
  using BoxT = Box<Dim, T>;
  using ProjBoxT = Box<Dim - 1, T>;

  AABBJoinEnumerator(const Relation<Dim, T>* rel_r,
                     const Relation<Dim, T>* rel_s,
                     int axis = 0,
                     join::SideTieBreak tb = join::SideTieBreak::RBeforeS)
      : R_(rel_r), S_(rel_s), axis_(axis), side_order_(tb) {
    if (R_) {
      proj_r_.resize(R_->Size());
      for (usize i = 0; i < R_->Size(); ++i) {
        proj_r_[i] = ProjectDropFirst<Dim, T>((*R_).boxes[i]);
      }
      handle_r_.assign(R_->Size(), kNull);
      tree_r_.ReserveLeaves(R_->Size());
    }
    if (S_) {
      proj_s_.resize(S_->Size());
      for (usize i = 0; i < S_->Size(); ++i) {
        proj_s_[i] = ProjectDropFirst<Dim, T>((*S_).boxes[i]);
      }
      handle_s_.assign(S_->Size(), kNull);
      tree_s_.ReserveLeaves(S_->Size());
    }
    Reset();
  }

  void Reset() override {
    stats_.Reset();
    events_.clear();
    cur_candidates_.clear();
    cur_pos_ = 0;
    scanning_ = false;
    cur_side_ = join::Side::R;
    cur_index_ = 0;
    cur_id_ = kInvalidId;
    ev_pos_ = 0;

    tree_r_.Clear();
    tree_s_.Clear();
    std::fill(handle_r_.begin(), handle_r_.end(), kNull);
    std::fill(handle_s_.begin(), handle_s_.end(), kNull);

    if (!R_ || !S_) return;

    events_ = join::BuildSweepEvents<Dim, T>(*R_, *S_, axis_, side_order_);
    stats_.num_events = static_cast<u64>(events_.size());
  }

  bool Next(PairId* out) override {
    if (!out) return false;
    if (!R_ || !S_) return false;

    while (true) {
      if (scanning_) {
        if (cur_pos_ < cur_candidates_.size()) {
          const u32 other_idx = cur_candidates_[cur_pos_++];
          if (cur_side_ == join::Side::R) {
            *out = PairId{cur_id_, S_->GetId(static_cast<usize>(other_idx))};
          } else {
            *out = PairId{R_->GetId(static_cast<usize>(other_idx)), cur_id_};
          }
          ++stats_.output_pairs;
          return true;
        }

        // Finish this START: insert current into its own tree.
        if (cur_side_ == join::Side::R) {
          handle_r_[cur_index_] = tree_r_.Insert(static_cast<u32>(cur_index_), proj_r_[cur_index_]);
          stats_.active_max_r = std::max(stats_.active_max_r, static_cast<u64>(tree_r_.LeafCount()));
        } else {
          handle_s_[cur_index_] = tree_s_.Insert(static_cast<u32>(cur_index_), proj_s_[cur_index_]);
          stats_.active_max_s = std::max(stats_.active_max_s, static_cast<u64>(tree_s_.LeafCount()));
        }
        scanning_ = false;
        cur_candidates_.clear();
        cur_pos_ = 0;
        continue;
      }

      if (ev_pos_ >= events_.size()) return false;
      const join::Event& e = events_[ev_pos_++];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) {
          const int h = handle_r_[e.index];
          if (h != kNull) tree_r_.Remove(h);
          handle_r_[e.index] = kNull;
        } else {
          const int h = handle_s_[e.index];
          if (h != kNull) tree_s_.Remove(h);
          handle_s_[e.index] = kNull;
        }
        continue;
      }

      // START event: produce candidate list by querying opposite tree.
      cur_side_ = e.side;
      cur_index_ = e.index;
      cur_id_ = e.id;

      cur_candidates_.clear();
      typename DynamicAABBTree<Dim - 1, T>::QueryStats qst;
      if (cur_side_ == join::Side::R) {
        tree_s_.ReportIntersect(proj_r_[cur_index_], &cur_candidates_, &qst);
      } else {
        tree_r_.ReportIntersect(proj_s_[cur_index_], &cur_candidates_, &qst);
      }
      stats_.candidate_checks += qst.leaf_tests;

      cur_pos_ = 0;
      scanning_ = true;
      // Loop: either output first candidate, or (if none) fall through to insert.
    }
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  static constexpr int kNull = -1;

  const Relation<Dim, T>* R_ = nullptr;
  const Relation<Dim, T>* S_ = nullptr;
  int axis_ = 0;
  join::SideTieBreak side_order_ = join::SideTieBreak::RBeforeS;

  std::vector<join::Event> events_;

  // Projection caches.
  std::vector<ProjBoxT> proj_r_;
  std::vector<ProjBoxT> proj_s_;

  // Dynamic indexes over the active sets.
  DynamicAABBTree<Dim - 1, T> tree_r_;
  DynamicAABBTree<Dim - 1, T> tree_s_;

  // Handle per object index for fast deletion.
  std::vector<int> handle_r_;
  std::vector<int> handle_s_;

  // Current START buffering.
  std::vector<u32> cur_candidates_;
  usize cur_pos_ = 0;
  bool scanning_ = false;
  join::Side cur_side_ = join::Side::R;
  usize cur_index_ = 0;
  Id cur_id_ = kInvalidId;

  // Event cursor.
  usize ev_pos_ = 0;

  join::JoinStats stats_;
};

}  // namespace detail

// --------------------------
// AABBSamplingBaseline
// --------------------------

template <int Dim, class T>
class AABBSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "AABBSamplingBaseline requires Dim >= 2");
  using DatasetT = Dataset<Dim, T>;
  using ProjBoxT = Box<Dim - 1, T>;

  Method method() const noexcept override { return Method::AABB; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "aabb_sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;

    events_.clear();
    start_id_of_event_.clear();
    w_total_.clear();

    proj_r_.clear();
    proj_s_.clear();
    handle_r_.clear();
    handle_s_.clear();

    tree_r_.Clear();
    tree_s_.Clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;

    if (ds.R.Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        ds.S.Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "AABBSamplingBaseline::Build: relation size exceeds u32 limit";
      return false;
    }

    // Build sweep events on axis 0.
    events_ = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);
    start_id_of_event_.assign(events_.size(), -1);
    usize start_cnt = 0;
    for (usize i = 0; i < events_.size(); ++i) {
      if (events_[i].kind == join::EventKind::Start) {
        start_id_of_event_[i] = static_cast<i32>(start_cnt);
        ++start_cnt;
      }
    }
    w_total_.assign(start_cnt, 0ULL);

    // Projection caches (drop axis 0).
    proj_r_.resize(ds.R.Size());
    for (usize i = 0; i < ds.R.Size(); ++i) {
      proj_r_[i] = detail::ProjectDropFirst<Dim, T>(ds.R.boxes[i]);
    }
    proj_s_.resize(ds.S.Size());
    for (usize i = 0; i < ds.S.Size(); ++i) {
      proj_s_[i] = detail::ProjectDropFirst<Dim, T>(ds.S.boxes[i]);
    }

    handle_r_.assign(ds.R.Size(), kNull);
    handle_s_.assign(ds.S.Size(), kNull);

    tree_r_.Clear();
    tree_s_.Clear();
    tree_r_.ReserveLeaves(ds.R.Size());
    tree_s_.ReserveLeaves(ds.S.Size());

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;  // deterministic
    if (!built_ || !ds_) {
      if (err) *err = "AABBSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "AABBSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    W_ = 0;
    weights_valid_ = false;

    tree_r_.Clear();
    tree_s_.Clear();
    std::fill(handle_r_.begin(), handle_r_.end(), kNull);
    std::fill(handle_s_.begin(), handle_s_.end(), kNull);

    u64 W = 0;
    for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
      const join::Event& e = events_[ev_pos];
      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) {
          const int h = handle_r_[e.index];
          if (h != kNull) tree_r_.Remove(h);
          handle_r_[e.index] = kNull;
        } else {
          const int h = handle_s_[e.index];
          if (h != kNull) tree_s_.Remove(h);
          handle_s_[e.index] = kNull;
        }
        continue;
      }

      const i32 sid_i32 = start_id_of_event_[ev_pos];
      SJS_DASSERT(sid_i32 >= 0);
      const usize sid = static_cast<usize>(sid_i32);

      if (e.side == join::Side::R) {
        const ProjBoxT& q = proj_r_[e.index];
        const u64 w = tree_s_.CountIntersect(q);
        w_total_[sid] = w;
        W += w;
        handle_r_[e.index] = tree_r_.Insert(static_cast<u32>(e.index), q);
      } else {
        const ProjBoxT& q = proj_s_[e.index];
        const u64 w = tree_r_.CountIntersect(q);
        w_total_[sid] = w;
        W += w;
        handle_s_[e.index] = tree_s_.Insert(static_cast<u32>(e.index), q);
      }
    }

    W_ = W;
    weights_valid_ = true;
    *out = MakeExactCount(W);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "AABBSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "AABBSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "AABBSamplingBaseline::Sample: out is null";
      return false;
    }

    const u32 t = cfg.run.t;

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();
    if (t == 0) return true;

    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }
    if (W_ == 0) {
      // Empty join.
      return true;
    }

    struct Assignment {
      u32 eid;
      u32 slot;
    };
    auto assign_less = [](const Assignment& a, const Assignment& b) noexcept {
      if (a.eid < b.eid) return true;
      if (b.eid < a.eid) return false;
      return a.slot < b.slot;
    };

    // --------------------------
    // Phase 2: alias on events + slot assignments
    // --------------------------
    std::vector<Assignment> assign;
    {
      auto scoped = phases ? phases->Scoped("phase2_assign") : PhaseRecorder::ScopedPhase(nullptr, "");

      // Build alias only on positive-weight START events (w_e > 0), as required by the v2.0 design.
      std::vector<u32> pos_eids;
      std::vector<u64> pos_w;
      pos_eids.reserve(w_total_.size());
      pos_w.reserve(w_total_.size());
      for (u32 eid = 0; eid < static_cast<u32>(w_total_.size()); ++eid) {
        const u64 w = w_total_[static_cast<usize>(eid)];
        if (w == 0) continue;
        pos_eids.push_back(eid);
        pos_w.push_back(w);
      }
      if (pos_w.empty()) {
        if (err) *err = "AABBSamplingBaseline::Sample: internal error (W>0 but no positive-weight events)";
        return false;
      }

      sampling::AliasTable alias;
      if (!alias.BuildFromU64(Span<const u64>(pos_w), err)) {
        if (err && err->empty()) *err = "AABBSamplingBaseline::Sample: failed to build alias table";
        return false;
      }

      assign.reserve(t);
      for (u32 j = 0; j < t; ++j) {
        const u32 idx = static_cast<u32>(alias.Sample(rng));
        const u32 eid = pos_eids[static_cast<usize>(idx)];
        assign.push_back(Assignment{eid, j});
      }

      std::sort(assign.begin(), assign.end(), assign_less);
    }

    // Prepare output.
    out->pairs.assign(static_cast<usize>(t), PairId{});

    // --------------------------
    // Phase 3: second sweep + local sampling + slot fill
    // --------------------------
    {
      auto scoped = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

      tree_r_.Clear();
      tree_s_.Clear();
      std::fill(handle_r_.begin(), handle_r_.end(), kNull);
      std::fill(handle_s_.begin(), handle_s_.end(), kNull);

      usize a_ptr = 0;
      std::vector<u32> picked;

      for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
        const join::Event& e = events_[ev_pos];
        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) {
            const int h = handle_r_[e.index];
            if (h != kNull) tree_r_.Remove(h);
            handle_r_[e.index] = kNull;
          } else {
            const int h = handle_s_[e.index];
            if (h != kNull) tree_s_.Remove(h);
            handle_s_[e.index] = kNull;
          }
          continue;
        }

        const i32 sid_i32 = start_id_of_event_[ev_pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 eid = static_cast<u32>(sid_i32);

        // Find assignment range for this event id.
        const usize begin = a_ptr;
        while (a_ptr < assign.size() && assign[a_ptr].eid == eid) {
          ++a_ptr;
        }
        const u32 t_e = static_cast<u32>(a_ptr - begin);

        if (t_e > 0) {
          picked.clear();
          if (e.side == join::Side::R) {
            const ProjBoxT& q = proj_r_[e.index];
            if (!tree_s_.SampleIntersect(q, t_e, rng, &picked)) {
              if (err) *err = "AABBSamplingBaseline::Sample: SampleIntersect returned empty unexpectedly";
              return false;
            }
            const Id rid = ds_->R.GetId(e.index);
            for (u32 u = 0; u < t_e; ++u) {
              const u32 s_idx = picked[static_cast<usize>(u)];
              out->pairs[static_cast<usize>(assign[begin + u].slot)] = PairId{rid, ds_->S.GetId(s_idx)};
            }
          } else {
            const ProjBoxT& q = proj_s_[e.index];
            if (!tree_r_.SampleIntersect(q, t_e, rng, &picked)) {
              if (err) *err = "AABBSamplingBaseline::Sample: SampleIntersect returned empty unexpectedly";
              return false;
            }
            const Id sid = ds_->S.GetId(e.index);
            for (u32 u = 0; u < t_e; ++u) {
              const u32 r_idx = picked[static_cast<usize>(u)];
              out->pairs[static_cast<usize>(assign[begin + u].slot)] = PairId{ds_->R.GetId(r_idx), sid};
            }
          }
        }

        // Insert this START into its own active tree.
        if (e.side == join::Side::R) {
          const ProjBoxT& q = proj_r_[e.index];
          handle_r_[e.index] = tree_r_.Insert(static_cast<u32>(e.index), q);
        } else {
          const ProjBoxT& q = proj_s_[e.index];
          handle_s_[e.index] = tree_s_.Insert(static_cast<u32>(e.index), q);
        }
      }

      if (a_ptr != assign.size()) {
        if (err) *err = "AABBSamplingBaseline::Sample: internal error (did not consume all assignments)";
        return false;
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !ds_) {
      if (err) *err = "AABBSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::AABBJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                join::SideTieBreak::RBeforeS);
  }

 private:
  static constexpr int kNull = -1;

  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  // Precomputed sweep events and mapping from event position to START-id.
  std::vector<join::Event> events_;
  std::vector<i32> start_id_of_event_;  // -1 for END; [0,|E|) for START

  // Weight per START event.
  std::vector<u64> w_total_;
  u64 W_ = 0;
  bool weights_valid_ = false;

  // Projection caches.
  std::vector<ProjBoxT> proj_r_;
  std::vector<ProjBoxT> proj_s_;

  // Dynamic trees for active sets.
  detail::DynamicAABBTree<Dim - 1, T> tree_r_;
  detail::DynamicAABBTree<Dim - 1, T> tree_s_;

  // Handle arrays (indexed by relation index) for deletion during sweeps.
  std::vector<int> handle_r_;
  std::vector<int> handle_s_;
};

}  // namespace aabb
}  // namespace baselines
}  // namespace sjs
