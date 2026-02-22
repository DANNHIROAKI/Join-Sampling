#pragma once
// sjs/index/segtree_common.h
//
// Common utilities for array-based segment trees used in SJS-HighDims (§2).
//
// We use the standard implicit binary tree layout:
//   - Choose p = next power-of-two >= n.
//   - Node indices are in [1, 2*p).
//   - Leaves for positions i in [0, p) are at index (p + i).
//   - Root is at index 1.
//
// This header provides:
//   - next power-of-two helper
//   - canonical cover decomposition of an interval [l, r)
//   - root-to-leaf / leaf-to-root paths

#include "sjs/core/assert.h"
#include "sjs/core/types.h"

#include <vector>

namespace sjs {
namespace index {

inline u32 NextPow2(u32 n) {
  if (n <= 1) return 1;
  --n;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  return n + 1;
}

inline u32 CeilLog2(u32 x) {
  SJS_DASSERT(x >= 1);
  u32 lg = 0;
  u32 p = 1;
  while (p < x) {
    p <<= 1;
    ++lg;
  }
  return lg;
}

// Canonical cover nodes for interval [l, r) in a tree with leaf base p.
// The returned nodes are disjoint and cover exactly [l, r).
inline void DecomposeCover(u32 l, u32 r, u32 p, std::vector<u32>* out) {
  SJS_DASSERT(out != nullptr);
  out->clear();
  if (r <= l) return;

  u32 L = l + p;
  u32 R = r + p;
  while (L < R) {
    if (L & 1U) out->push_back(L++);
    if (R & 1U) out->push_back(--R);
    L >>= 1;
    R >>= 1;
  }
}

// Canonical cover nodes in deterministic left-to-right order.
// This is useful for REPORT to be reproducible.
inline void DecomposeCoverOrdered(u32 l, u32 r, u32 p, std::vector<u32>* out) {
  SJS_DASSERT(out != nullptr);
  out->clear();
  if (r <= l) return;

  std::vector<u32> left;
  std::vector<u32> right;
  left.reserve(64);
  right.reserve(64);

  u32 L = l + p;
  u32 R = r + p;
  while (L < R) {
    if (L & 1U) left.push_back(L++);
    if (R & 1U) right.push_back(--R);
    L >>= 1;
    R >>= 1;
  }

  out->reserve(left.size() + right.size());
  out->insert(out->end(), left.begin(), left.end());
  for (usize i = 0; i < right.size(); ++i) {
    out->push_back(right[right.size() - 1 - i]);
  }
}

// Leaf-to-root path for position `pos` (0-based) in a tree with base p.
// Output order is leaf, parent, ..., root.
inline void PathToRoot(u32 pos, u32 p, std::vector<u32>* out) {
  SJS_DASSERT(out != nullptr);
  out->clear();
  u32 idx = pos + p;
  while (idx > 0) {
    out->push_back(idx);
    idx >>= 1;
  }
}

// Root-to-leaf path for position `pos` (0-based) in a tree with base p.
// Output order is root, ..., leaf.
inline void PathFromRoot(u32 pos, u32 p, std::vector<u32>* out) {
  SJS_DASSERT(out != nullptr);
  out->clear();
  std::vector<u32> tmp;
  PathToRoot(pos, p, &tmp);
  out->reserve(tmp.size());
  for (usize i = 0; i < tmp.size(); ++i) {
    out->push_back(tmp[tmp.size() - 1 - i]);
  }
}

}  // namespace index
}  // namespace sjs
