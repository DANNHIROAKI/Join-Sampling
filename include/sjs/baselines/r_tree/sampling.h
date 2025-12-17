#pragma once
// sjs/baselines/r_tree/sampling.h
//
// Plane Sweep + Dynamic R-Tree baseline (Variant::Sampling).
//
// Implements the algorithm described in "R-Tree Baseline.md":
//   - Sweep on axis 0 (x1).
//   - Maintain two dynamic R-trees over the projected boxes on dims 1..Dim-1
//     (drop the sweep axis) for the active sets.
//   - Phase 1: for each START event e (query q*), compute w_e = CountIntersect(q*)
//     on the opposite active R-tree.
//   - Phase 2: build an alias table over {w_e} and assign t sample slots to events.
//   - Phase 3: second sweep; for each START with t_e>0, draw t_e i.i.d. uniform
//     samples from the local intersection set K_e using R-tree sampling.
//
// Notes
// -----
// * Half-open box semantics [lo,hi) are used everywhere (consistent with sjs::Box).
// * This header also provides a deterministic join enumerator based on the same
//   sweep + R-tree machinery (used by other variants and runners).
// * The dynamic R-tree implementation is baseline-grade (correct, deterministic,
//   reasonably fast) but not a full-blown R*-tree.

#include "sjs/baselines/baseline_api.h"
#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/join/sweep_events.h"
#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace r_tree {

namespace detail {

// --------------------------
// Small helpers
// --------------------------

inline bool ParseU32(std::string_view s, u32* out) {
  if (!out) return false;
  if (s.empty()) return false;
  try {
    std::size_t idx = 0;
    const unsigned long long v = std::stoull(std::string(s), &idx, 10);
    if (idx != s.size()) return false;
    if (v > static_cast<unsigned long long>(std::numeric_limits<u32>::max())) return false;
    *out = static_cast<u32>(v);
    return true;
  } catch (...) {
    return false;
  }
}

inline u32 ExtraU32Or(const std::unordered_map<std::string, std::string>& extra,
                      std::string_view key,
                      u32 default_v) {
  auto it = extra.find(std::string(key));
  if (it == extra.end()) return default_v;
  u32 v = default_v;
  if (!ParseU32(it->second, &v)) return default_v;
  return v;
}

// Project a Dim-dimensional box to Dim-1 by dropping axis 0.
// For Dim==2 this is just the y-interval.

template <int Dim, class T>
inline Box<Dim - 1, T> ProjectDropFirst(const Box<Dim, T>& b) noexcept {
  static_assert(Dim >= 2, "ProjectDropFirst requires Dim >= 2");
  Box<Dim - 1, T> out;
  for (int d = 1; d < Dim; ++d) {
    out.lo.v[d - 1] = b.lo.v[d];
    out.hi.v[d - 1] = b.hi.v[d];
  }
  return out;
}

// A monotone, overflow-resistant measure used for R-tree heuristics.
// We use L1 "perimeter-like" measure (sum of side lengths) instead of volume.
// This behaves better in higher dimensions and for very small rectangles.

template <int Dim, class T>
inline long double Measure(const Box<Dim, T>& b) noexcept {
  if (b.IsEmpty()) return 0.0L;
  long double s = 0.0L;
  for (int i = 0; i < Dim; ++i) {
    const long double w = static_cast<long double>(b.hi.v[i] - b.lo.v[i]);
    if (!(w > 0.0L)) return 0.0L;
    s += w;
  }
  return s;
}

template <int Dim, class T>
inline Box<Dim, T> UnionBox(const Box<Dim, T>& a, const Box<Dim, T>& b) noexcept {
  Box<Dim, T> u = a;
  u.ExpandToIncludeBox(b);
  return u;
}

// --------------------------
// Dynamic R-tree (baseline-grade)
// --------------------------

// This is a classic Guttman-style R-tree with:
//   - ChooseLeaf: minimal enlargement, then minimal area.
//   - Split: quadratic split.
//   - Delete: direct removal by (leaf,slot) handle and upward recomputation.
//
// We intentionally do NOT implement full CondenseTree underflow reinsertion.
// For our plane sweep workload (insert once, delete once) this is sufficient
// for correctness, and keeps the implementation compact.
//
// The tree stores only currently-active objects (insert at START, remove at END).

template <int Dim, class T>
class DynamicRTree {
 public:
  static_assert(Dim >= 1, "DynamicRTree requires Dim >= 1");
  using BoxT = Box<Dim, T>;

  struct Options {
    // Maximum children per node (M). Typical values: 16, 32, 64.
    u32 max_children = 32;
    // Minimum children per node after split (m). If 0, we set to ceil(M/2).
    u32 min_children = 0;

    // If true, ignore duplicate inserts of the same object id.
    // If false, duplicate insert is a debug assertion.
    bool ignore_duplicate_insert = true;
  };

  struct QueryStats {
    u64 node_tests = 0;   // number of node-level bbox tests
    u64 entry_tests = 0;  // number of leaf entry bbox tests

    void Reset() {
      node_tests = 0;
      entry_tests = 0;
    }
  };

  struct Handle {
    u32 node = kNull;
    u32 slot = 0;
    bool valid = false;
  };

  DynamicRTree() = default;

  void Init(u32 capacity, Options opt = {}) {
    opt_ = opt;
    if (opt_.max_children < 2) opt_.max_children = 2;
    if (opt_.min_children == 0) {
      opt_.min_children = (opt_.max_children + 1) / 2;  // ceil(M/2)
    }
    if (opt_.min_children < 1) opt_.min_children = 1;
    if (opt_.min_children > opt_.max_children) opt_.min_children = opt_.max_children;

    capacity_ = capacity;
    handles_.assign(static_cast<usize>(capacity_), Handle{});
    Clear();
  }

  void Clear() {
    nodes_.clear();
    root_ = kNull;
    // Keep handle storage but invalidate.
    for (auto& h : handles_) {
      h.node = kNull;
      h.slot = 0;
      h.valid = false;
    }
  }

  bool Empty() const noexcept { return root_ == kNull || nodes_[static_cast<usize>(root_)].size == 0; }

  u32 Size() const noexcept {
    if (root_ == kNull) return 0;
    return nodes_[static_cast<usize>(root_)].size;
  }

  bool Contains(u32 obj) const noexcept {
    if (obj >= capacity_) return false;
    return handles_[static_cast<usize>(obj)].valid;
  }

  // Insert an active object with its bounding box.
  // Returns false only on invalid input (obj out of range, empty box, etc.).
  bool Insert(u32 obj, const BoxT& b, std::string* err = nullptr) {
    if (obj >= capacity_) {
      if (err) *err = "DynamicRTree::Insert: obj out of range";
      return false;
    }
    if (b.IsEmpty()) {
      // Empty boxes contribute no join pairs; ignore insertion.
      return false;
    }

    Handle& h = handles_[static_cast<usize>(obj)];
    if (h.valid) {
      if (opt_.ignore_duplicate_insert) return true;
      SJS_DASSERT_MSG(false, "DynamicRTree::Insert: duplicate insert");
      if (err) *err = "DynamicRTree::Insert: duplicate insert";
      return false;
    }

    if (root_ == kNull) {
      const u32 r = NewNode(/*leaf=*/true, /*parent=*/kNull);
      root_ = r;
    }

    const u32 leaf = ChooseLeaf(root_, b);
    Node& ln = nodes_[static_cast<usize>(leaf)];
    SJS_DASSERT(ln.leaf);

    const u32 slot = static_cast<u32>(ln.children.size());
    ln.children.push_back(Entry{b, obj});
    h = Handle{leaf, slot, true};

    // Fix up bboxes/sizes and handle overflow.
    AdjustAfterInsert(leaf);
    return true;
  }

  // Remove an active object. Returns true if removed, false if obj not active.
  bool Remove(u32 obj) {
    if (obj >= capacity_) return false;
    Handle& h = handles_[static_cast<usize>(obj)];
    if (!h.valid) return false;
    const u32 leaf = h.node;
    const u32 slot = h.slot;
    if (leaf == kNull) return false;

    Node& ln = nodes_[static_cast<usize>(leaf)];
    SJS_DASSERT(ln.leaf);
    SJS_DASSERT(slot < ln.children.size());
    SJS_DASSERT(ln.children[static_cast<usize>(slot)].ref == obj);

    const u32 last = static_cast<u32>(ln.children.size() - 1);
    if (slot != last) {
      const Entry moved = ln.children[static_cast<usize>(last)];
      ln.children[static_cast<usize>(slot)] = moved;
      // Update handle of moved object.
      Handle& hm = handles_[static_cast<usize>(moved.ref)];
      SJS_DASSERT(hm.valid);
      hm.node = leaf;
      hm.slot = slot;
    }
    ln.children.pop_back();

    // Invalidate.
    h.valid = false;
    h.node = kNull;
    h.slot = 0;

    // Recompute upwards.
    AdjustUpwards(leaf);
    return true;
  }

  // Count active objects whose boxes intersect q.
  u64 CountIntersect(const BoxT& q, QueryStats* st = nullptr) const {
    if (root_ == kNull) return 0;
    if (q.IsEmpty()) return 0;

    if (st) st->Reset();

    u64 cnt = 0;
    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const u32 nid = stack.back();
      stack.pop_back();
      const Node& n = nodes_[static_cast<usize>(nid)];
      if (n.size == 0) continue;

      if (st) st->node_tests++;
      if (!n.bbox.Intersects(q)) continue;

      if (q.ContainsBox(n.bbox)) {
        cnt += static_cast<u64>(n.size);
        continue;
      }

      if (n.leaf) {
        for (const auto& e : n.children) {
          if (st) st->entry_tests++;
          if (e.bbox.Intersects(q)) cnt++;
        }
      } else {
        for (const auto& e : n.children) {
          stack.push_back(e.ref);
        }
      }
    }

    return cnt;
  }

  // Report active objects whose boxes intersect q.
  void ReportIntersect(const BoxT& q, std::vector<u32>* out, QueryStats* st = nullptr) const {
    if (!out) return;
    out->clear();
    if (root_ == kNull) return;
    if (q.IsEmpty()) return;

    if (st) st->Reset();

    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const u32 nid = stack.back();
      stack.pop_back();
      const Node& n = nodes_[static_cast<usize>(nid)];
      if (n.size == 0) continue;

      if (st) st->node_tests++;
      if (!n.bbox.Intersects(q)) continue;

      if (n.leaf) {
        for (const auto& e : n.children) {
          if (st) st->entry_tests++;
          if (e.bbox.Intersects(q)) out->push_back(e.ref);
        }
      } else {
        for (const auto& e : n.children) {
          stack.push_back(e.ref);
        }
      }
    }
  }

  // Draw k i.i.d. uniform samples WITH replacement from the set
  // { obj : bbox(obj) intersects q }.
  // Returns false if the intersection set is empty.
  bool SampleIntersect(const BoxT& q, u32 k, Rng* rng, std::vector<u32>* out, QueryStats* st = nullptr) const {
    if (!rng || !out) return false;
    out->clear();
    if (k == 0) return true;
    if (root_ == kNull) return false;
    if (q.IsEmpty()) return false;

    if (st) st->Reset();

    // Step A: build frontier components.
    struct Component {
      enum class Type : u8 { FullNode = 0, LeafHits = 1 };
      Type type = Type::LeafHits;
      u32 node = kNull;              // for FullNode
      std::vector<u32> hits;         // for LeafHits
      u32 weight = 0;                // #objects in this component
    };

    std::vector<Component> comps;
    comps.reserve(64);

    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const u32 nid = stack.back();
      stack.pop_back();
      const Node& n = nodes_[static_cast<usize>(nid)];
      if (n.size == 0) continue;

      if (st) st->node_tests++;
      if (!n.bbox.Intersects(q)) continue;

      if (q.ContainsBox(n.bbox)) {
        Component c;
        c.type = Component::Type::FullNode;
        c.node = nid;
        c.weight = n.size;
        comps.push_back(std::move(c));
        continue;
      }

      if (n.leaf) {
        Component c;
        c.type = Component::Type::LeafHits;
        c.node = nid;
        c.hits.clear();
        for (const auto& e : n.children) {
          if (st) st->entry_tests++;
          if (e.bbox.Intersects(q)) c.hits.push_back(e.ref);
        }
        c.weight = static_cast<u32>(c.hits.size());
        if (c.weight > 0) comps.push_back(std::move(c));
      } else {
        for (const auto& e : n.children) stack.push_back(e.ref);
      }
    }

    u64 total = 0;
    for (const auto& c : comps) total += static_cast<u64>(c.weight);
    if (total == 0) return false;

    // Step B: alias over components.
    std::vector<u64> w;
    w.reserve(comps.size());
    for (const auto& c : comps) w.push_back(static_cast<u64>(c.weight));

    sampling::AliasTable alias;
    if (!alias.BuildFromU64(Span<const u64>(w))) {
      return false;
    }

    out->reserve(static_cast<usize>(k));
    for (u32 i = 0; i < k; ++i) {
      const u32 ci = static_cast<u32>(alias.Sample(rng));
      const Component& c = comps[static_cast<usize>(ci)];
      if (c.type == Component::Type::LeafHits) {
        const u32 j = rng->UniformU32(c.weight);
        out->push_back(c.hits[static_cast<usize>(j)]);
      } else {
        const u32 obj = SampleFromSubtree(c.node, rng);
        out->push_back(obj);
      }
    }

    return true;
  }

 private:
  static constexpr u32 kNull = std::numeric_limits<u32>::max();

  struct Entry {
    BoxT bbox;
    u32 ref = 0;  // leaf: obj id, internal: child node id
  };

  struct Node {
    bool leaf = true;
    u32 parent = kNull;
    BoxT bbox = BoxT::Empty();
    u32 size = 0;  // number of leaf entries in subtree (active objects)
    std::vector<Entry> children;
  };

  Options opt_{};
  u32 capacity_ = 0;

  std::vector<Node> nodes_;
  u32 root_ = kNull;

  std::vector<Handle> handles_;

  u32 NewNode(bool leaf, u32 parent) {
    Node n;
    n.leaf = leaf;
    n.parent = parent;
    n.bbox = BoxT::Empty();
    n.size = 0;
    n.children.clear();
    n.children.reserve(static_cast<usize>(opt_.max_children + 1));
    nodes_.push_back(std::move(n));
    return static_cast<u32>(nodes_.size() - 1);
  }

  // Recompute bbox and size of a node from its children.
  void RecomputeNode(u32 nid) {
    Node& n = nodes_[static_cast<usize>(nid)];
    n.bbox = BoxT::Empty();
    if (n.leaf) {
      n.size = static_cast<u32>(n.children.size());
      for (const auto& e : n.children) {
        n.bbox.ExpandToIncludeBox(e.bbox);
      }
    } else {
      u32 sz = 0;
      for (auto& e : n.children) {
        const Node& ch = nodes_[static_cast<usize>(e.ref)];
        e.bbox = ch.bbox;  // keep cached bbox consistent
        sz += ch.size;
        n.bbox.ExpandToIncludeBox(e.bbox);
      }
      n.size = sz;
    }
  }

  // Find the index of child entry in parent that points to child_node.
  u32 FindChildSlot(u32 parent, u32 child_node) const {
    const Node& p = nodes_[static_cast<usize>(parent)];
    for (u32 i = 0; i < static_cast<u32>(p.children.size()); ++i) {
      if (p.children[static_cast<usize>(i)].ref == child_node) return i;
    }
    return kNull;
  }

  // Choose leaf by minimal enlargement.
  u32 ChooseLeaf(u32 start, const BoxT& b) const {
    u32 nid = start;
    while (true) {
      const Node& n = nodes_[static_cast<usize>(nid)];
      if (n.leaf) return nid;

      // Choose the child whose bbox needs least enlargement to include b.
      u32 best_i = 0;
      long double best_enl = std::numeric_limits<long double>::infinity();
      long double best_meas = std::numeric_limits<long double>::infinity();
      u32 best_ref = 0;

      for (u32 i = 0; i < static_cast<u32>(n.children.size()); ++i) {
        const Entry& e = n.children[static_cast<usize>(i)];
        const BoxT u = UnionBox<Dim, T>(e.bbox, b);
        const long double m0 = Measure<Dim, T>(e.bbox);
        const long double m1 = Measure<Dim, T>(u);
        const long double enl = m1 - m0;

        const u32 ref = e.ref;

        if (enl < best_enl) {
          best_enl = enl;
          best_meas = m0;
          best_i = i;
          best_ref = ref;
        } else if (enl == best_enl) {
          if (m0 < best_meas) {
            best_meas = m0;
            best_i = i;
            best_ref = ref;
          } else if (m0 == best_meas) {
            // Deterministic tie-break.
            if (ref < best_ref) {
              best_i = i;
              best_ref = ref;
            }
          }
        }
      }

      nid = n.children[static_cast<usize>(best_i)].ref;
    }
  }

  void AdjustUpwards(u32 from) {
    u32 nid = from;
    while (nid != kNull) {
      RecomputeNode(nid);
      nid = nodes_[static_cast<usize>(nid)].parent;
    }
  }

  // Quadratic split of node `nid`.
  // The node is split into (nid) and (new_nid). Returns new_nid.
  u32 SplitNode(u32 nid) {
    Node& n = nodes_[static_cast<usize>(nid)];
    const bool leaf = n.leaf;

    // Take ownership of current children.
    std::vector<Entry> entries = std::move(n.children);
    n.children.clear();
    n.children.reserve(static_cast<usize>(opt_.max_children + 1));

    const u32 new_nid = NewNode(leaf, n.parent);
    Node& nn = nodes_[static_cast<usize>(new_nid)];

    const u32 N = static_cast<u32>(entries.size());
    SJS_DASSERT(N == opt_.max_children + 1);

    // Pick seeds with max "waste".
    u32 seed_a = 0;
    u32 seed_b = 1;
    long double best_waste = -1.0L;
    for (u32 i = 0; i < N; ++i) {
      for (u32 j = i + 1; j < N; ++j) {
        const BoxT u = UnionBox<Dim, T>(entries[static_cast<usize>(i)].bbox,
                                       entries[static_cast<usize>(j)].bbox);
        const long double waste = Measure<Dim, T>(u) - Measure<Dim, T>(entries[static_cast<usize>(i)].bbox) -
                                  Measure<Dim, T>(entries[static_cast<usize>(j)].bbox);
        if (waste > best_waste || (waste == best_waste && (i < seed_a || (i == seed_a && j < seed_b)))) {
          best_waste = waste;
          seed_a = i;
          seed_b = j;
        }
      }
    }

    std::vector<u8> used(static_cast<usize>(N), 0);
    used[static_cast<usize>(seed_a)] = 1;
    used[static_cast<usize>(seed_b)] = 1;

    std::vector<Entry> group_a;
    std::vector<Entry> group_b;
    group_a.reserve(static_cast<usize>(opt_.max_children + 1));
    group_b.reserve(static_cast<usize>(opt_.max_children + 1));

    group_a.push_back(entries[static_cast<usize>(seed_a)]);
    group_b.push_back(entries[static_cast<usize>(seed_b)]);

    BoxT bbox_a = group_a[0].bbox;
    BoxT bbox_b = group_b[0].bbox;

    const u32 minc = opt_.min_children;

    auto remaining_count = [&]() -> u32 {
      u32 r = 0;
      for (u32 i = 0; i < N; ++i) if (!used[static_cast<usize>(i)]) ++r;
      return r;
    };

    while (true) {
      const u32 rem = remaining_count();
      if (rem == 0) break;

      // Force assignment to satisfy min-children.
      if (static_cast<u32>(group_a.size()) + rem == minc) {
        for (u32 i = 0; i < N; ++i) {
          if (used[static_cast<usize>(i)]) continue;
          used[static_cast<usize>(i)] = 1;
          group_a.push_back(entries[static_cast<usize>(i)]);
          bbox_a = UnionBox<Dim, T>(bbox_a, entries[static_cast<usize>(i)].bbox);
        }
        break;
      }
      if (static_cast<u32>(group_b.size()) + rem == minc) {
        for (u32 i = 0; i < N; ++i) {
          if (used[static_cast<usize>(i)]) continue;
          used[static_cast<usize>(i)] = 1;
          group_b.push_back(entries[static_cast<usize>(i)]);
          bbox_b = UnionBox<Dim, T>(bbox_b, entries[static_cast<usize>(i)].bbox);
        }
        break;
      }

      // Pick next entry with maximum |dA - dB|.
      u32 best_k = kNull;
      long double best_diff = -1.0L;
      long double best_dA = 0.0L;
      long double best_dB = 0.0L;

      for (u32 i = 0; i < N; ++i) {
        if (used[static_cast<usize>(i)]) continue;
        const BoxT& eb = entries[static_cast<usize>(i)].bbox;
        const long double mA0 = Measure<Dim, T>(bbox_a);
        const long double mB0 = Measure<Dim, T>(bbox_b);
        const long double mA1 = Measure<Dim, T>(UnionBox<Dim, T>(bbox_a, eb));
        const long double mB1 = Measure<Dim, T>(UnionBox<Dim, T>(bbox_b, eb));
        const long double dA = mA1 - mA0;
        const long double dB = mB1 - mB0;
        const long double diff = (dA > dB) ? (dA - dB) : (dB - dA);

        if (diff > best_diff || (diff == best_diff && i < best_k)) {
          best_diff = diff;
          best_k = i;
          best_dA = dA;
          best_dB = dB;
        }
      }

      SJS_DASSERT(best_k != kNull);
      used[static_cast<usize>(best_k)] = 1;
      const Entry e = entries[static_cast<usize>(best_k)];

      // Assign to group with smaller enlargement.
      bool to_a = false;
      if (best_dA < best_dB) {
        to_a = true;
      } else if (best_dB < best_dA) {
        to_a = false;
      } else {
        // Tie: smaller bbox measure.
        const long double mA = Measure<Dim, T>(bbox_a);
        const long double mB = Measure<Dim, T>(bbox_b);
        if (mA < mB) {
          to_a = true;
        } else if (mB < mA) {
          to_a = false;
        } else {
          // Tie: smaller group.
          if (group_a.size() <= group_b.size()) to_a = true;
          else to_a = false;
        }
      }

      if (to_a) {
        group_a.push_back(e);
        bbox_a = UnionBox<Dim, T>(bbox_a, e.bbox);
      } else {
        group_b.push_back(e);
        bbox_b = UnionBox<Dim, T>(bbox_b, e.bbox);
      }

      // Defensive: ensure we never violate min children feasibility.
      // (This should be guaranteed by the forced-assignment checks above.)
    }

    // Assign groups.
    n.children = std::move(group_a);
    nn.children = std::move(group_b);

    // Fix parent pointers for moved child nodes (internal) OR handles (leaf).
    if (leaf) {
      for (u32 i = 0; i < static_cast<u32>(n.children.size()); ++i) {
        const u32 obj = n.children[static_cast<usize>(i)].ref;
        Handle& h = handles_[static_cast<usize>(obj)];
        SJS_DASSERT(h.valid);
        h.node = nid;
        h.slot = i;
      }
      for (u32 i = 0; i < static_cast<u32>(nn.children.size()); ++i) {
        const u32 obj = nn.children[static_cast<usize>(i)].ref;
        Handle& h = handles_[static_cast<usize>(obj)];
        SJS_DASSERT(h.valid);
        h.node = new_nid;
        h.slot = i;
      }
    } else {
      for (auto& e : n.children) {
        nodes_[static_cast<usize>(e.ref)].parent = nid;
      }
      for (auto& e : nn.children) {
        nodes_[static_cast<usize>(e.ref)].parent = new_nid;
      }
    }

    // Recompute both nodes.
    RecomputeNode(nid);
    RecomputeNode(new_nid);

    return new_nid;
  }

  void AdjustAfterInsert(u32 leaf) {
    u32 nid = leaf;

    // Walk up, splitting any node that overflows.
    u32 new_nid = kNull;
    if (nodes_[static_cast<usize>(nid)].children.size() > opt_.max_children) {
      new_nid = SplitNode(nid);
    } else {
      AdjustUpwards(nid);
      return;
    }

    while (true) {
      // nid and new_nid are siblings after split.
      if (nodes_[static_cast<usize>(nid)].parent == kNull) {
        // Create a new root.
        const u32 new_root = NewNode(/*leaf=*/false, /*parent=*/kNull);
        nodes_[static_cast<usize>(new_root)].children.push_back(Entry{nodes_[static_cast<usize>(nid)].bbox, nid});
        nodes_[static_cast<usize>(new_root)].children.push_back(
            Entry{nodes_[static_cast<usize>(new_nid)].bbox, new_nid});
        nodes_[static_cast<usize>(nid)].parent = new_root;
        nodes_[static_cast<usize>(new_nid)].parent = new_root;
        root_ = new_root;
        RecomputeNode(new_root);
        break;
      }

      const u32 parent = nodes_[static_cast<usize>(nid)].parent;

      // Update parent's cached bbox for nid.
      {
        const u32 slot = FindChildSlot(parent, nid);
        SJS_DASSERT(slot != kNull);
        nodes_[static_cast<usize>(parent)].children[static_cast<usize>(slot)].bbox =
            nodes_[static_cast<usize>(nid)].bbox;
      }

      // Add new_nid to parent.
      nodes_[static_cast<usize>(parent)].children.push_back(
          Entry{nodes_[static_cast<usize>(new_nid)].bbox, new_nid});
      nodes_[static_cast<usize>(new_nid)].parent = parent;

      // Recompute parent.
      RecomputeNode(parent);

      if (nodes_[static_cast<usize>(parent)].children.size() > opt_.max_children) {
        // Split parent and continue upwards.
        new_nid = SplitNode(parent);
        nid = parent;
        continue;
      }

      // No overflow: adjust remaining ancestors and finish.
      AdjustUpwards(parent);
      break;
    }
  }

  u32 SampleFromSubtree(u32 nid, Rng* rng) const {
    SJS_DASSERT(rng);
    SJS_DASSERT(nid != kNull);

    u32 cur = nid;
    while (true) {
      const Node& n = nodes_[static_cast<usize>(cur)];
      SJS_DASSERT(n.size > 0);
      if (n.leaf) {
        const u32 j = rng->UniformU32(static_cast<u32>(n.children.size()));
        return n.children[static_cast<usize>(j)].ref;
      }

      // Choose child proportional to child.size.
      const u64 total = static_cast<u64>(n.size);
      u64 r = rng->UniformU64(total);
      for (const auto& e : n.children) {
        const Node& ch = nodes_[static_cast<usize>(e.ref)];
        const u32 cs = ch.size;
        if (cs == 0) continue;
        if (r < static_cast<u64>(cs)) {
          cur = e.ref;
          break;
        }
        r -= static_cast<u64>(cs);
      }
    }
  }
};

// --------------------------
// Deterministic enumerator: Plane Sweep + R-tree
// --------------------------

template <int Dim, class T>
class RTreeJoinEnumerator final : public baselines::IJoinEnumerator {
 public:
  static_assert(Dim >= 2, "RTreeJoinEnumerator requires Dim >= 2");
  using BoxT = Box<Dim, T>;
  using ProjBoxT = Box<Dim - 1, T>;

  explicit RTreeJoinEnumerator(const Relation<Dim, T>* rel_r,
                              const Relation<Dim, T>* rel_s,
                              typename DynamicRTree<Dim - 1, T>::Options opt = {},
                              int axis = 0,
                              join::SideTieBreak tb = join::SideTieBreak::RBeforeS)
      : R_(rel_r), S_(rel_s), axis_(axis), side_order_(tb), opt_(opt) {
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

    proj_r_.clear();
    proj_s_.clear();

    tree_r_.Clear();
    tree_s_.Clear();

    if (!R_ || !S_) return;

    // Projection caches.
    proj_r_.resize(R_->Size());
    for (usize i = 0; i < R_->Size(); ++i) proj_r_[i] = ProjectDropFirst<Dim, T>(R_->boxes[i]);
    proj_s_.resize(S_->Size());
    for (usize i = 0; i < S_->Size(); ++i) proj_s_[i] = ProjectDropFirst<Dim, T>(S_->boxes[i]);

    tree_r_.Init(static_cast<u32>(R_->Size()), opt_);
    tree_s_.Init(static_cast<u32>(S_->Size()), opt_);

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
          (void)tree_r_.Insert(static_cast<u32>(cur_index_), proj_r_[cur_index_]);
          stats_.active_max_r = std::max(stats_.active_max_r, static_cast<u64>(tree_r_.Size()));
        } else {
          (void)tree_s_.Insert(static_cast<u32>(cur_index_), proj_s_[cur_index_]);
          stats_.active_max_s = std::max(stats_.active_max_s, static_cast<u64>(tree_s_.Size()));
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
          (void)tree_r_.Remove(static_cast<u32>(e.index));
        } else {
          (void)tree_s_.Remove(static_cast<u32>(e.index));
        }
        continue;
      }

      // START event: query opposite tree.
      cur_side_ = e.side;
      cur_index_ = e.index;
      cur_id_ = e.id;

      cur_candidates_.clear();
      typename DynamicRTree<Dim - 1, T>::QueryStats qst;
      if (cur_side_ == join::Side::R) {
        tree_s_.ReportIntersect(proj_r_[cur_index_], &cur_candidates_, &qst);
      } else {
        tree_r_.ReportIntersect(proj_s_[cur_index_], &cur_candidates_, &qst);
      }
      stats_.candidate_checks += qst.entry_tests;

      cur_pos_ = 0;
      scanning_ = true;
    }
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  const Relation<Dim, T>* R_ = nullptr;
  const Relation<Dim, T>* S_ = nullptr;
  int axis_ = 0;
  join::SideTieBreak side_order_ = join::SideTieBreak::RBeforeS;
  typename DynamicRTree<Dim - 1, T>::Options opt_{};

  std::vector<join::Event> events_;
  std::vector<ProjBoxT> proj_r_;
  std::vector<ProjBoxT> proj_s_;

  DynamicRTree<Dim - 1, T> tree_r_;
  DynamicRTree<Dim - 1, T> tree_s_;

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
// RTreeSamplingBaseline
// --------------------------

template <int Dim, class T>
class RTreeSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "RTreeSamplingBaseline requires Dim >= 2");
  using DatasetT = Dataset<Dim, T>;
  using ProjBoxT = Box<Dim - 1, T>;

  Method method() const noexcept override { return Method::RTree; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "rtree_sampling"; }

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

    tree_r_.Clear();
    tree_s_.Clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;

    if (ds.R.Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        ds.S.Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RTreeSamplingBaseline::Build: relation size exceeds u32 limit";
      return false;
    }

    // Parse R-tree options (max children) from cfg.run.extra.
    typename detail::DynamicRTree<Dim - 1, T>::Options opt;
    {
      // Accept several aliases for convenience.
      u32 m = 32;
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_max_children", m);
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_M", m);
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_m", m);
      if (m < 2) m = 2;
      opt.max_children = m;

      u32 minc = 0;
      minc = detail::ExtraU32Or(cfg.run.extra, "rtree_min_children", minc);
      if (minc > 0) opt.min_children = minc;

      // Duplicate insert shouldn't happen in sweep; keep it strict.
      opt.ignore_duplicate_insert = true;
    }
    opt_ = opt;

    // Build sweep events on axis 0.
    events_ = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);
    start_id_of_event_.assign(events_.size(), -1);
    usize start_cnt = 0;
    for (usize i = 0; i < events_.size(); ++i) {
      if (events_[i].kind == join::EventKind::Start) {
        start_id_of_event_[i] = static_cast<i32>(start_cnt++);
      }
    }
    w_total_.assign(start_cnt, 0ULL);

    // Projection caches.
    proj_r_.resize(ds.R.Size());
    for (usize i = 0; i < ds.R.Size(); ++i) proj_r_[i] = detail::ProjectDropFirst<Dim, T>(ds.R.boxes[i]);
    proj_s_.resize(ds.S.Size());
    for (usize i = 0; i < ds.S.Size(); ++i) proj_s_[i] = detail::ProjectDropFirst<Dim, T>(ds.S.boxes[i]);

    // Init trees (capacity only); they are cleared per sweep.
    tree_r_.Init(static_cast<u32>(ds.R.Size()), opt_);
    tree_s_.Init(static_cast<u32>(ds.S.Size()), opt_);

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
      if (err) *err = "RTreeSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RTreeSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    W_ = 0;
    weights_valid_ = false;

    tree_r_.Clear();
    tree_s_.Clear();

    u64 W = 0;
    for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
      const join::Event& e = events_[ev_pos];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) {
          (void)tree_r_.Remove(static_cast<u32>(e.index));
        } else {
          (void)tree_s_.Remove(static_cast<u32>(e.index));
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
        (void)tree_r_.Insert(static_cast<u32>(e.index), q);
      } else {
        const ProjBoxT& q = proj_s_[e.index];
        const u64 w = tree_r_.CountIntersect(q);
        w_total_[sid] = w;
        W += w;
        (void)tree_s_.Insert(static_cast<u32>(e.index), q);
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
      if (err) *err = "RTreeSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "RTreeSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "RTreeSamplingBaseline::Sample: out is null";
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
      auto scoped2 = phases ? phases->Scoped("phase2_assign") : PhaseRecorder::ScopedPhase(nullptr, "");

      sampling::AliasTable alias;
      if (!alias.BuildFromU64(Span<const u64>(w_total_), err)) {
        if (err && err->empty()) *err = "RTreeSamplingBaseline::Sample: failed to build alias table";
        return false;
      }

      assign.reserve(t);
      for (u32 j = 0; j < t; ++j) {
        u32 eid = 0;
        u64 w = 0;
        for (int tries = 0; tries < 16; ++tries) {
          eid = static_cast<u32>(alias.Sample(rng));
          w = w_total_[eid];
          if (w > 0) break;
        }
        if (w == 0) {
          if (err) *err = "RTreeSamplingBaseline::Sample: alias produced only zero-weight events (unexpected)";
          return false;
        }
        assign.push_back(Assignment{eid, j});
      }
      std::sort(assign.begin(), assign.end(), assign_less);
    }

    out->pairs.assign(static_cast<usize>(t), PairId{});

    // --------------------------
    // Phase 3: second sweep + local sampling + slot fill
    // --------------------------
    {
      auto scoped3 = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

      tree_r_.Clear();
      tree_s_.Clear();

      usize a_ptr = 0;
      std::vector<u32> picked;

      for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
        const join::Event& e = events_[ev_pos];

        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) {
            (void)tree_r_.Remove(static_cast<u32>(e.index));
          } else {
            (void)tree_s_.Remove(static_cast<u32>(e.index));
          }
          continue;
        }

        const i32 sid_i32 = start_id_of_event_[ev_pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 eid = static_cast<u32>(sid_i32);

        // Find assignment range for this event id.
        const usize begin = a_ptr;
        while (a_ptr < assign.size() && assign[a_ptr].eid == eid) ++a_ptr;
        const u32 t_e = static_cast<u32>(a_ptr - begin);

        if (t_e > 0) {
          picked.clear();
          if (e.side == join::Side::R) {
            const ProjBoxT& q = proj_r_[e.index];
            if (!tree_s_.SampleIntersect(q, t_e, rng, &picked)) {
              if (err) *err = "RTreeSamplingBaseline::Sample: SampleIntersect returned empty unexpectedly";
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
              if (err) *err = "RTreeSamplingBaseline::Sample: SampleIntersect returned empty unexpectedly";
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
          (void)tree_r_.Insert(static_cast<u32>(e.index), q);
        } else {
          const ProjBoxT& q = proj_s_[e.index];
          (void)tree_s_.Insert(static_cast<u32>(e.index), q);
        }
      }

      if (a_ptr != assign.size()) {
        if (err) *err = "RTreeSamplingBaseline::Sample: internal error (did not consume all assignments)";
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
      if (err) *err = "RTreeSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RTreeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, opt_, /*axis=*/0,
                                                                 join::SideTieBreak::RBeforeS);
  }

 private:
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
  detail::DynamicRTree<Dim - 1, T> tree_r_;
  detail::DynamicRTree<Dim - 1, T> tree_s_;

  // R-tree options.
  typename detail::DynamicRTree<Dim - 1, T>::Options opt_{};
};

}  // namespace r_tree
}  // namespace baselines
}  // namespace sjs
