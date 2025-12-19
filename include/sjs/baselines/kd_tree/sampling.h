#pragma once
// sjs/baselines/kd_tree/sampling.h
//
// Plane Sweep + Global Rank Embedding + Static KD-Tree w/ dynamic active counts
// baseline (Variant::Sampling).
//
// This implements the algorithm described in "KD-Tree Baseline v2.0.md":
//   - Sweep on axis 0 (x1).
//   - Embed rectangles on the remaining axes (2..d) into points in 2(d-1) dims
//     using global ranks of (L_i, id) and (R_i, id) to make strict inequalities
//     deterministic for half-open boxes.
//   - Build a static KD-tree on each side's embedded points.
//   - Maintain dynamic active set by toggling points on START/END events and
//     updating subtree active counts.
//   - Phase 1: for each START event e (query box q), compute w_e = Count(Q(q))
//     on the opposite tree.
//   - Phase 2: build an alias table on {w_e} after filtering out w_e==0 (v2.0 checklist #4),
//     then assign t sample slots to events.
//   - Phase 3: second sweep; for each START with t_e>0, draw t_e i.i.d. uniform
//     samples from the local intersection set K_e using KD-tree range sampling.
//
// Notes:
//   - Correct for half-open boxes [lo,hi) in each dimension.
//   - This header also provides a deterministic join enumerator based on the
//     same sweep + KD machinery (used by other variants and runners).

#include "sjs/baselines/baseline_api.h"
#include "sjs/core/assert.h"
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
namespace kd_tree {

namespace detail {

// --------------------------
// Rank-space point/box (half-open)
// --------------------------

template <int K>
struct RankPoint {
  static_assert(K >= 1, "RankPoint<K>: K must be >= 1");
  std::array<u32, K> v{};

  u32& operator[](int i) noexcept { return v[static_cast<usize>(i)]; }
  const u32& operator[](int i) const noexcept { return v[static_cast<usize>(i)]; }
};

template <int K>
struct RankBox {
  static_assert(K >= 1, "RankBox<K>: K must be >= 1");
  RankPoint<K> lo{};
  RankPoint<K> hi{};

  bool IsEmpty() const noexcept {
    for (int d = 0; d < K; ++d) {
      if (!(lo[d] < hi[d])) return true;
    }
    return false;
  }
};

template <int K>
inline bool ContainsPointHalfOpen(const RankBox<K>& b, const RankPoint<K>& p) noexcept {
  for (int d = 0; d < K; ++d) {
    const u32 x = p[d];
    if (x < b.lo[d] || !(x < b.hi[d])) return false;
  }
  return true;
}

template <int K>
inline bool IntersectsHalfOpen(const RankBox<K>& a, const RankBox<K>& b) noexcept {
  for (int d = 0; d < K; ++d) {
    if (!(a.lo[d] < b.hi[d] && b.lo[d] < a.hi[d])) return false;
  }
  return true;
}

template <int K>
inline bool ContainsBoxHalfOpen(const RankBox<K>& outer, const RankBox<K>& inner) noexcept {
  for (int d = 0; d < K; ++d) {
    if (outer.lo[d] > inner.lo[d]) return false;
    if (inner.hi[d] > outer.hi[d]) return false;
  }
  return true;
}

template <int K>
inline RankBox<K> BoxFromPoint(const RankPoint<K>& p) noexcept {
  RankBox<K> b;
  for (int d = 0; d < K; ++d) {
    b.lo[d] = p[d];
    b.hi[d] = p[d] + 1;  // half-open singleton
  }
  return b;
}

template <int K>
inline void UnionInto(RankBox<K>* a, const RankBox<K>& b) noexcept {
  SJS_DASSERT(a);
  for (int d = 0; d < K; ++d) {
    if (b.lo[d] < a->lo[d]) a->lo[d] = b.lo[d];
    if (a->hi[d] < b.hi[d]) a->hi[d] = b.hi[d];
  }
}

// --------------------------
// Global rank embedding (shared coordinate system for both relations)
// --------------------------

template <int Dim, class T>
struct GlobalRankEmbedding {
  static_assert(Dim >= 2, "GlobalRankEmbedding requires Dim >= 2");
  static constexpr int M = Dim - 1;          // number of non-sweep dimensions
  static constexpr int K = 2 * (Dim - 1);    // embedding dimension

  struct Key {
    T value{};
    u32 gid{};  // global id for tie-break (must be unique across R and S)
  };

  struct KeyLess {
    bool operator()(const Key& a, const Key& b) const noexcept {
      if (a.value < b.value) return true;
      if (b.value < a.value) return false;
      return a.gid < b.gid;
    }
  };

  void Reset() {
    for (int i = 0; i < M; ++i) {
      A_[i].clear();
      B_[i].clear();
    }
    n_r_ = 0;
    n_s_ = 0;
    n_total_ = 0;
    built_ = false;
  }

  bool Built() const noexcept { return built_; }
  u32 Total() const noexcept { return n_total_; }
  u32 NR() const noexcept { return n_r_; }
  u32 NS() const noexcept { return n_s_; }

  // Build global rank arrays and output embedded points for each relation.
  // Global ids are assigned deterministically:
  //   gid(R[i]) = 1 + i
  //   gid(S[j]) = 1 + n_r + j
  bool Build(const Relation<Dim, T>& R,
             const Relation<Dim, T>& S,
             std::vector<RankPoint<K>>* out_pts_r,
             std::vector<RankPoint<K>>* out_pts_s,
             PhaseRecorder* phases,
             std::string* err) {
    Reset();

    if (!out_pts_r || !out_pts_s) {
      if (err) *err = "GlobalRankEmbedding::Build: out_pts_r/out_pts_s is null";
      return false;
    }

    const usize nr = R.Size();
    const usize ns = S.Size();
    const usize nt = nr + ns;

    if (nt > static_cast<usize>(std::numeric_limits<u32>::max()) - 2ULL) {
      if (err) *err = "GlobalRankEmbedding::Build: dataset too large for 32-bit ranks";
      return false;
    }

    n_r_ = static_cast<u32>(nr);
    n_s_ = static_cast<u32>(ns);
    n_total_ = static_cast<u32>(nt);

    // Prepare output point vectors.
    out_pts_r->assign(nr, RankPoint<K>{});
    out_pts_s->assign(ns, RankPoint<K>{});

    // Build A_i (L ranks) and B_i (R ranks) for each non-sweep dimension.
    {
      auto _ = phases ? phases->Scoped("build_rank_arrays") : PhaseRecorder::ScopedPhase(nullptr, "");
      KeyLess less;
      for (int i = 0; i < M; ++i) {
        A_[i].reserve(nt);
        B_[i].reserve(nt);
      }

      // Fill from R.
      for (u32 idx = 0; idx < n_r_; ++idx) {
        const u32 gid = 1U + idx;
        const auto& b = R.boxes[static_cast<usize>(idx)];
        for (int i = 0; i < M; ++i) {
          const int axis = i + 1;  // skip sweep axis 0
          A_[i].push_back(Key{b.lo.v[axis], gid});
          B_[i].push_back(Key{b.hi.v[axis], gid});
        }
      }

      // Fill from S.
      for (u32 idx = 0; idx < n_s_; ++idx) {
        const u32 gid = 1U + n_r_ + idx;
        const auto& b = S.boxes[static_cast<usize>(idx)];
        for (int i = 0; i < M; ++i) {
          const int axis = i + 1;
          A_[i].push_back(Key{b.lo.v[axis], gid});
          B_[i].push_back(Key{b.hi.v[axis], gid});
        }
      }

      // Sort and assign ranks to embedded points.
      for (int i = 0; i < M; ++i) {
        std::sort(A_[i].begin(), A_[i].end(), less);
        std::sort(B_[i].begin(), B_[i].end(), less);

        // rankA_i -> coordinate i
        for (u32 r = 0; r < n_total_; ++r) {
          const u32 gid = A_[i][static_cast<usize>(r)].gid;
          if (gid >= 1U && gid <= n_r_) {
            const u32 h = gid - 1U;
            (*out_pts_r)[static_cast<usize>(h)].v[static_cast<usize>(i)] = r;
          } else {
            const u32 h = gid - 1U - n_r_;
            if (h < n_s_) {
              (*out_pts_s)[static_cast<usize>(h)].v[static_cast<usize>(i)] = r;
            }
          }
        }

        // rankB_i -> coordinate (i + M)
        for (u32 r = 0; r < n_total_; ++r) {
          const u32 gid = B_[i][static_cast<usize>(r)].gid;
          if (gid >= 1U && gid <= n_r_) {
            const u32 h = gid - 1U;
            (*out_pts_r)[static_cast<usize>(h)].v[static_cast<usize>(i + M)] = r;
          } else {
            const u32 h = gid - 1U - n_r_;
            if (h < n_s_) {
              (*out_pts_s)[static_cast<usize>(h)].v[static_cast<usize>(i + M)] = r;
            }
          }
        }
      }
    }

    built_ = true;
    return true;
  }

  // Build the orthogonal query range Q(q) in rank space (half-open).
  // Returns false if the range is empty.
  bool MakeQueryRange(const Box<Dim, T>& q, RankBox<K>* out) const noexcept {
    if (!out) return false;
    if (!built_) return false;

    RankBox<K> qb;

    // For each non-sweep axis a=1..Dim-1 (mapped to i=0..M-1):
    //   L_a(r) < R_a(q)  -> rankA_i(r) in [0, lower_bound(A_i, (R_a(q), LOW)))
    //   R_a(r) > L_a(q)  -> rankB_i(r) in [upper_bound(B_i, (L_a(q), HIGH)), n_total)
    const KeyLess less;
    const u32 LOW = 0;
    const u32 HIGH = std::numeric_limits<u32>::max();

    for (int i = 0; i < M; ++i) {
      const int axis = i + 1;
      const T rhs = q.hi.v[axis];  // R_a(q)
      const T lhs = q.lo.v[axis];  // L_a(q)

      // ub = first position with (value, gid) >= (rhs, LOW)
      const Key k_hi{rhs, LOW};
      const auto it_hi = std::lower_bound(A_[i].begin(), A_[i].end(), k_hi, less);
      const u32 ub = static_cast<u32>(std::distance(A_[i].begin(), it_hi));

      // lb = first position with (value, gid) > (lhs, HIGH)
      const Key k_lo{lhs, HIGH};
      const auto it_lo = std::upper_bound(B_[i].begin(), B_[i].end(), k_lo, less);
      const u32 lb = static_cast<u32>(std::distance(B_[i].begin(), it_lo));

      qb.lo[i] = 0;
      qb.hi[i] = ub;

      qb.lo[i + M] = lb;
      qb.hi[i + M] = n_total_;
    }

    if (qb.IsEmpty()) return false;
    *out = qb;
    return true;
  }

 private:
  std::array<std::vector<Key>, M> A_{};  // L ranks
  std::array<std::vector<Key>, M> B_{};  // R ranks

  u32 n_r_ = 0;
  u32 n_s_ = 0;
  u32 n_total_ = 0;
  bool built_ = false;
};

// --------------------------
// Static KD-tree with dynamic active counts (range count/report/sample)
// --------------------------

template <int K>
class ActiveKDTree {
 public:
  static_assert(K >= 1, "ActiveKDTree<K>: K must be >= 1");

  using PointT = RankPoint<K>;
  using BoxT = RankBox<K>;

  struct QueryStats {
    u64 nodes_visited = 0;
    u64 point_tests = 0;
    u64 reported = 0;
  };

  ActiveKDTree() = default;

  void Clear() {
    nodes_.clear();
    node_of_handle_.clear();
    indices_.clear();
    root_ = kNull;
  }

  usize Size() const noexcept { return node_of_handle_.size(); }

  // Force the active set to empty (O(n)).
  void ResetToEmpty() {
    for (auto& n : nodes_) {
      n.active = 0;
      n.cnt = 0;
    }
  }

  bool Build(const std::vector<PointT>& pts, std::string* err = nullptr) {
    Clear();
    const usize n = pts.size();
    if (n == 0) {
      root_ = kNull;
      return true;
    }
    if (n > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "ActiveKDTree::Build: too many points for u32 handles";
      return false;
    }

    node_of_handle_.assign(n, kNull);
    indices_.resize(n);
    for (u32 i = 0; i < static_cast<u32>(n); ++i) indices_[static_cast<usize>(i)] = i;

    nodes_.reserve(n);
    root_ = BuildRec(pts, /*begin=*/0, /*end=*/n, /*depth=*/0, /*parent=*/kNull);

    // Ensure empty active set.
    ResetToEmpty();
    return true;
  }

  // Activate/deactivate by handle (0..n-1).
  void Activate(u32 handle) {
    if (handle >= node_of_handle_.size()) return;
    const u32 v = node_of_handle_[static_cast<usize>(handle)];
    if (v == kNull) return;
    if (nodes_[static_cast<usize>(v)].active) return;

    nodes_[static_cast<usize>(v)].active = 1;
    u32 cur = v;
    while (cur != kNull) {
      ++nodes_[static_cast<usize>(cur)].cnt;
      cur = nodes_[static_cast<usize>(cur)].parent;
    }
  }

  void Deactivate(u32 handle) {
    if (handle >= node_of_handle_.size()) return;
    const u32 v = node_of_handle_[static_cast<usize>(handle)];
    if (v == kNull) return;
    if (!nodes_[static_cast<usize>(v)].active) return;

    nodes_[static_cast<usize>(v)].active = 0;
    u32 cur = v;
    while (cur != kNull) {
      Node& n = nodes_[static_cast<usize>(cur)];
      SJS_DASSERT(n.cnt > 0);
      --n.cnt;
      cur = n.parent;
    }
  }

  u32 ActiveCount() const noexcept {
    if (root_ == kNull) return 0;
    return nodes_[static_cast<usize>(root_)].cnt;
  }

  // Count active points in range.
  u64 Count(const BoxT& q, QueryStats* st = nullptr) const {
    if (root_ == kNull) return 0;
    if (q.IsEmpty()) return 0;

    u64 ans = 0;
    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const u32 v = stack.back();
      stack.pop_back();

      const Node& n = nodes_[static_cast<usize>(v)];
      if (st) ++st->nodes_visited;

      if (n.cnt == 0) continue;
      if (!IntersectsHalfOpen<K>(n.bbox, q)) continue;

      if (ContainsBoxHalfOpen<K>(q, n.bbox)) {
        ans += static_cast<u64>(n.cnt);
        continue;
      }

      if (n.active) {
        if (st) ++st->point_tests;
        if (ContainsPointHalfOpen<K>(q, n.pt)) {
          ans += 1;
        }
      }

      if (n.left != kNull) stack.push_back(n.left);
      if (n.right != kNull) stack.push_back(n.right);
    }

    if (st) st->reported += ans;
    return ans;
  }

  // Report active point handles in range (deterministic order).
  void Report(const BoxT& q, std::vector<u32>* out, QueryStats* st = nullptr) const {
    if (!out) return;
    out->clear();
    if (root_ == kNull) return;
    if (q.IsEmpty()) return;

    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const u32 v = stack.back();
      stack.pop_back();

      const Node& n = nodes_[static_cast<usize>(v)];
      if (st) ++st->nodes_visited;

      if (n.cnt == 0) continue;
      if (!IntersectsHalfOpen<K>(n.bbox, q)) continue;

      if (ContainsBoxHalfOpen<K>(q, n.bbox)) {
        CollectAllActive(v, out);
        continue;
      }

      if (n.active) {
        if (st) ++st->point_tests;
        if (ContainsPointHalfOpen<K>(q, n.pt)) {
          out->push_back(n.handle);
        }
      }

      // Deterministic traversal: left before right.
      if (n.right != kNull) stack.push_back(n.right);
      if (n.left != kNull) stack.push_back(n.left);
    }

    if (st) st->reported += static_cast<u64>(out->size());
  }

  // Sample t handles uniformly with replacement from active points in range.
  // Returns false if the range is empty (no candidates).
  bool Sample(const BoxT& q, u32 t, Rng* rng, std::vector<u32>* out, std::string* err = nullptr) const {
    if (!out) {
      if (err) *err = "ActiveKDTree::Sample: out is null";
      return false;
    }
    out->clear();
    if (t == 0) return true;
    if (!rng) {
      if (err) *err = "ActiveKDTree::Sample: rng is null";
      return false;
    }
    if (root_ == kNull) return false;
    if (q.IsEmpty()) return false;

    std::vector<Piece> pieces;
    pieces.reserve(64);
    const u64 total = CollectPieces(q, &pieces);
    if (total == 0 || pieces.empty()) return false;

    std::vector<u64> weights;
    weights.resize(pieces.size());
    for (usize i = 0; i < pieces.size(); ++i) {
      weights[i] = static_cast<u64>(pieces[i].weight);
    }

    sampling::AliasTable alias;
    if (!alias.BuildFromU64(Span<const u64>(weights.data(), weights.size()), err)) {
      return false;
    }

    // Assign each output slot to a piece, then fill per piece.
    struct Assign {
      u32 piece;
      u32 slot;
    };

    std::vector<Assign> asg;
    asg.reserve(static_cast<usize>(t));
    for (u32 j = 0; j < t; ++j) {
      const u32 p = static_cast<u32>(alias.Sample(rng));
      asg.push_back(Assign{p, j});
    }
    std::sort(asg.begin(), asg.end(), [](const Assign& a, const Assign& b) {
      if (a.piece < b.piece) return true;
      if (b.piece < a.piece) return false;
      return a.slot < b.slot;
    });

    out->resize(static_cast<usize>(t));

    usize ptr = 0;
    while (ptr < asg.size()) {
      const u32 pidx = asg[ptr].piece;
      SJS_DASSERT(pidx < pieces.size());
      const Piece& pc = pieces[static_cast<usize>(pidx)];

      // Group [ptr, end)
      usize end = ptr + 1;
      while (end < asg.size() && asg[end].piece == pidx) ++end;

      if (pc.kind == PieceKind::Point) {
        const u32 h = nodes_[static_cast<usize>(pc.node)].handle;
        for (usize k = ptr; k < end; ++k) {
          (*out)[static_cast<usize>(asg[k].slot)] = h;
        }
      } else {
        for (usize k = ptr; k < end; ++k) {
          const u32 h = SampleUniformInSubtree(pc.node, rng);
          (*out)[static_cast<usize>(asg[k].slot)] = h;
        }
      }

      ptr = end;
    }

    return true;
  }

 private:
  static constexpr u32 kNull = std::numeric_limits<u32>::max();

  struct Node {
    PointT pt{};
    BoxT bbox{};
    u32 handle = 0;
    u32 left = kNull;
    u32 right = kNull;
    u32 parent = kNull;
    u32 cnt = 0;
    u8 active = 0;
  };

  enum class PieceKind : u8 { Subtree = 0, Point = 1 };
  struct Piece {
    u32 node;
    PieceKind kind;
    u32 weight;
  };

  // Build recursion.
  u32 BuildRec(const std::vector<PointT>& pts, usize begin, usize end, int depth, u32 parent) {
    if (!(begin < end)) return kNull;

    const int axis = (K == 0) ? 0 : (depth % K);
    const usize mid = begin + (end - begin) / 2;

    auto comp = [&](u32 a, u32 b) {
      const u32 va = pts[static_cast<usize>(a)].v[static_cast<usize>(axis)];
      const u32 vb = pts[static_cast<usize>(b)].v[static_cast<usize>(axis)];
      if (va < vb) return true;
      if (vb < va) return false;
      return a < b;
    };

    std::nth_element(indices_.begin() + static_cast<i64>(begin),
                     indices_.begin() + static_cast<i64>(mid),
                     indices_.begin() + static_cast<i64>(end),
                     comp);

    const u32 handle = indices_[static_cast<usize>(mid)];

    const u32 node_idx = static_cast<u32>(nodes_.size());
    nodes_.push_back(Node{});

    nodes_[static_cast<usize>(node_idx)].handle = handle;
    nodes_[static_cast<usize>(node_idx)].pt = pts[static_cast<usize>(handle)];
    nodes_[static_cast<usize>(node_idx)].parent = parent;

    // Recurse.
    const u32 left = BuildRec(pts, begin, mid, depth + 1, node_idx);
    const u32 right = BuildRec(pts, mid + 1, end, depth + 1, node_idx);

    nodes_[static_cast<usize>(node_idx)].left = left;
    nodes_[static_cast<usize>(node_idx)].right = right;

    // BBox from point union children.
    BoxT bb = BoxFromPoint<K>(nodes_[static_cast<usize>(node_idx)].pt);
    if (left != kNull) UnionInto<K>(&bb, nodes_[static_cast<usize>(left)].bbox);
    if (right != kNull) UnionInto<K>(&bb, nodes_[static_cast<usize>(right)].bbox);
    nodes_[static_cast<usize>(node_idx)].bbox = bb;

    // Handle -> node index.
    node_of_handle_[static_cast<usize>(handle)] = node_idx;

    return node_idx;
  }

  // Collect all active handles in subtree rooted at v (deterministic pre-order).
  void CollectAllActive(u32 v, std::vector<u32>* out) const {
    if (v == kNull || !out) return;
    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(v);
    while (!stack.empty()) {
      const u32 x = stack.back();
      stack.pop_back();
      const Node& n = nodes_[static_cast<usize>(x)];
      if (n.cnt == 0) continue;
      if (n.active) out->push_back(n.handle);
      if (n.right != kNull) stack.push_back(n.right);
      if (n.left != kNull) stack.push_back(n.left);
    }
  }

  // Collect canonical pieces for range sampling.
  u64 CollectPieces(const BoxT& q, std::vector<Piece>* pieces) const {
    if (!pieces) return 0;
    pieces->clear();
    if (root_ == kNull) return 0;

    u64 total = 0;
    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const u32 v = stack.back();
      stack.pop_back();
      const Node& n = nodes_[static_cast<usize>(v)];

      if (n.cnt == 0) continue;
      if (!IntersectsHalfOpen<K>(n.bbox, q)) continue;

      if (ContainsBoxHalfOpen<K>(q, n.bbox)) {
        pieces->push_back(Piece{v, PieceKind::Subtree, n.cnt});
        total += static_cast<u64>(n.cnt);
        continue;
      }

      if (n.active && ContainsPointHalfOpen<K>(q, n.pt)) {
        pieces->push_back(Piece{v, PieceKind::Point, 1});
        total += 1;
      }

      if (n.right != kNull) stack.push_back(n.right);
      if (n.left != kNull) stack.push_back(n.left);
    }

    return total;
  }

  // Uniform sample among active points in a subtree using cnt-guided descent.
  u32 SampleUniformInSubtree(u32 v, Rng* rng) const {
    SJS_DASSERT(rng);
    SJS_DASSERT(v != kNull);

    u32 cur = v;
    while (true) {
      const Node& n = nodes_[static_cast<usize>(cur)];
      const u32 w_self = n.active ? 1U : 0U;
      const u32 w_left = (n.left == kNull) ? 0U : nodes_[static_cast<usize>(n.left)].cnt;
      const u32 w_right = (n.right == kNull) ? 0U : nodes_[static_cast<usize>(n.right)].cnt;
      const u32 tot = w_self + w_left + w_right;
      SJS_DASSERT(tot == n.cnt);
      SJS_DASSERT(tot > 0);

      u64 r = rng->UniformU64(static_cast<u64>(tot));
      if (r < w_self) return n.handle;
      r -= w_self;

      if (r < w_left) {
        cur = n.left;
      } else {
        cur = n.right;
      }
      SJS_DASSERT(cur != kNull);
    }
  }

  std::vector<Node> nodes_;
  std::vector<u32> node_of_handle_;  // handle -> node index
  std::vector<u32> indices_;         // build-time permutation
  u32 root_ = kNull;
};

// --------------------------
// Deterministic join enumerator (sweep + ActiveKDTree)
// --------------------------

template <int Dim, class T>
class KDJoinEnumerator final : public baselines::IJoinEnumerator {
 public:
  static_assert(Dim >= 2, "KDJoinEnumerator requires Dim >= 2");
  using BoxT = Box<Dim, T>;
  static constexpr int K = 2 * (Dim - 1);

  KDJoinEnumerator(const Relation<Dim, T>* rel_r,
                   const Relation<Dim, T>* rel_s,
                   int axis = 0,
                   join::SideTieBreak tb = join::SideTieBreak::RBeforeS)
      : R_(rel_r), S_(rel_s), axis_(axis), side_order_(tb) {
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

    embed_.Reset();
    kd_r_.Clear();
    kd_s_.Clear();

    if (!R_ || !S_) return;

    // Build embedding and KD-trees (owned by this enumerator).
    std::vector<RankPoint<K>> pts_r;
    std::vector<RankPoint<K>> pts_s;
    if (!embed_.Build(*R_, *S_, &pts_r, &pts_s, /*phases=*/nullptr, /*err=*/nullptr)) {
      // If embedding fails, keep enumerator empty.
      return;
    }
    (void)kd_r_.Build(pts_r);
    (void)kd_s_.Build(pts_s);

    // Free temporary point arrays (KD nodes keep copies).
    std::vector<RankPoint<K>>().swap(pts_r);
    std::vector<RankPoint<K>>().swap(pts_s);

    // Build events.
    events_ = join::BuildSweepEvents<Dim, T>(*R_, *S_, axis_, side_order_);
    stats_.num_events = static_cast<u64>(events_.size());

    // Force empty active set.
    kd_r_.ResetToEmpty();
    kd_s_.ResetToEmpty();
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

        // Finish this START: activate current in its own tree.
        if (cur_side_ == join::Side::R) {
          kd_r_.Activate(cur_index_);
          stats_.active_max_r = std::max(stats_.active_max_r, static_cast<u64>(kd_r_.ActiveCount()));
        } else {
          kd_s_.Activate(cur_index_);
          stats_.active_max_s = std::max(stats_.active_max_s, static_cast<u64>(kd_s_.ActiveCount()));
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
          kd_r_.Deactivate(static_cast<u32>(e.index));
        } else {
          kd_s_.Deactivate(static_cast<u32>(e.index));
        }
        continue;
      }

      // START: query opposite active tree and emit candidates.
      cur_side_ = e.side;
      cur_index_ = static_cast<u32>(e.index);
      cur_id_ = e.id;

      const BoxT& q = (cur_side_ == join::Side::R)
                          ? R_->boxes[static_cast<usize>(cur_index_)]
                          : S_->boxes[static_cast<usize>(cur_index_)];

      RankBox<K> range;
      cur_candidates_.clear();

      typename ActiveKDTree<K>::QueryStats qst;
      if (embed_.MakeQueryRange(q, &range)) {
        if (cur_side_ == join::Side::R) {
          kd_s_.Report(range, &cur_candidates_, &qst);
        } else {
          kd_r_.Report(range, &cur_candidates_, &qst);
        }
      }

      // Interpret "candidate_checks" as point-level membership tests (+ accepted outputs).
      stats_.candidate_checks += qst.point_tests + qst.reported;

      cur_pos_ = 0;
      scanning_ = true;
      // Loop: either output first candidate, or if none, fall through to activate.
    }
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  const Relation<Dim, T>* R_ = nullptr;
  const Relation<Dim, T>* S_ = nullptr;

  int axis_ = 0;
  join::SideTieBreak side_order_ = join::SideTieBreak::RBeforeS;

  join::JoinStats stats_;
  std::vector<join::Event> events_;

  // Embedding + trees.
  GlobalRankEmbedding<Dim, T> embed_;
  ActiveKDTree<K> kd_r_;
  ActiveKDTree<K> kd_s_;

  // Current START emission state.
  std::vector<u32> cur_candidates_;
  usize cur_pos_ = 0;
  bool scanning_ = false;

  join::Side cur_side_ = join::Side::R;
  u32 cur_index_ = 0;
  Id cur_id_ = kInvalidId;

  usize ev_pos_ = 0;
};

}  // namespace detail

// --------------------------
// KDTreeSamplingBaseline
// --------------------------

template <int Dim, class T = Scalar>
class KDTreeSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "KDTreeSamplingBaseline requires Dim >= 2");
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;
  static constexpr int K = 2 * (Dim - 1);

  Method method() const noexcept override { return Method::KDTree; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "kd_tree_sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;

    events_.clear();
    start_id_of_event_.clear();
    w_total_.clear();

    embed_.Reset();
    kd_r_.Clear();
    kd_s_.Clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    ds_ = &ds;

    // Events.
    {
      auto _ = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    // START event ids.
    {
      auto _ = phases ? phases->Scoped("build_start_index") : PhaseRecorder::ScopedPhase(nullptr, "");
      start_id_of_event_.assign(events_.size(), kInvalidStartId);
      u32 sid = 0;
      for (usize i = 0; i < events_.size(); ++i) {
        if (events_[i].kind == join::EventKind::Start) {
          start_id_of_event_[i] = sid++;
        }
      }
      w_total_.assign(static_cast<usize>(sid), 0ULL);
    }

    // Build embedding + KD trees.
    {
      std::vector<detail::RankPoint<K>> pts_r;
      std::vector<detail::RankPoint<K>> pts_s;

      if (!embed_.Build(ds.R, ds.S, &pts_r, &pts_s, phases, err)) {
        return false;
      }

      auto _ = phases ? phases->Scoped("build_kdtrees") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!kd_r_.Build(pts_r, err)) return false;
      if (!kd_s_.Build(pts_s, err)) return false;

      // Free temporary point arrays (KD nodes keep copies).
      std::vector<detail::RankPoint<K>>().swap(pts_r);
      std::vector<detail::RankPoint<K>>().swap(pts_s);
    }

    built_ = true;
    weights_valid_ = false;
    W_ = 0;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;
    if (!built_ || !ds_) {
      if (err) *err = "KDTreeSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "KDTreeSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    // Empty relation => empty join.
    if (ds_->R.Size() == 0 || ds_->S.Size() == 0) {
      W_ = 0;
      weights_valid_ = true;
      std::fill(w_total_.begin(), w_total_.end(), 0ULL);
      *out = MakeExactCount(0);
      return true;
    }

    // Reset weights and active sets.
    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    kd_r_.ResetToEmpty();
    kd_s_.ResetToEmpty();

    u64 W = 0;
    for (usize ev_i = 0; ev_i < events_.size(); ++ev_i) {
      const join::Event& e = events_[ev_i];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) kd_r_.Deactivate(static_cast<u32>(e.index));
        else kd_s_.Deactivate(static_cast<u32>(e.index));
        continue;
      }

      // START
      const u32 sid = start_id_of_event_[ev_i];
      SJS_DASSERT(sid != kInvalidStartId);

      const BoxT& q = (e.side == join::Side::R)
                          ? ds_->R.boxes[static_cast<usize>(e.index)]
                          : ds_->S.boxes[static_cast<usize>(e.index)];

      detail::RankBox<K> range;
      u64 w = 0;
      if (embed_.MakeQueryRange(q, &range)) {
        if (e.side == join::Side::R) {
          w = kd_s_.Count(range);
        } else {
          w = kd_r_.Count(range);
        }
      }

      w_total_[static_cast<usize>(sid)] = w;

      // Accumulate W with overflow guard.
      if (w > 0) {
        if (W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "KDTreeSamplingBaseline::Count: |J| overflowed u64";
          return false;
        }
        W += w;
      }

      // Activate current.
      if (e.side == join::Side::R) kd_r_.Activate(static_cast<u32>(e.index));
      else kd_s_.Activate(static_cast<u32>(e.index));
    }

    W_ = W;
    weights_valid_ = true;
    *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "KDTreeSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "KDTreeSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "KDTreeSamplingBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t64 = cfg.run.t;
    if (t64 == 0) return true;
    if (t64 > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "KDTreeSamplingBaseline::Sample: t too large for u32 slots";
      return false;
    }
    const u32 t = static_cast<u32>(t64);

    // Ensure weights.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    if (W_ == 0) {
      // Empty join.
      return true;
    }

    // -----------------
    // Phase 2: alias + slot assignment (filter w_e==0)
    // -----------------
    {
      auto _ = phases ? phases->Scoped("phase2_alias") : PhaseRecorder::ScopedPhase(nullptr, "");

      // Baseline v2.0 (impl checklist #4): filter out START events with w_e==0
      // before building the alias table, so they are never sampled.
      std::vector<u32> nz_sids;
      std::vector<u64> nz_w;
      nz_sids.reserve(w_total_.size());
      nz_w.reserve(w_total_.size());
      for (u32 sid = 0; sid < static_cast<u32>(w_total_.size()); ++sid) {
        const u64 w = w_total_[static_cast<usize>(sid)];
        if (w == 0) continue;
        nz_sids.push_back(sid);
        nz_w.push_back(w);
      }
      if (nz_sids.empty()) {
        if (err) *err = "KDTreeSamplingBaseline::Sample: internal error (W_>0 but no positive w_e)";
        return false;
      }

      sampling::AliasTable alias;
      if (!alias.BuildFromU64(Span<const u64>(nz_w.data(), nz_w.size()), err)) {
        return false;
      }

      struct SlotAssign {
        u32 sid;   // original START id
        u32 slot;  // output position
      };

      std::vector<SlotAssign> asg;
      asg.reserve(static_cast<usize>(t));
      for (u32 j = 0; j < t; ++j) {
        const u32 p = static_cast<u32>(alias.Sample(rng));  // index into nz_sids / nz_w
        const u32 sid = nz_sids[static_cast<usize>(p)];
        asg.push_back(SlotAssign{sid, j});
      }
      std::sort(asg.begin(), asg.end(), [](const SlotAssign& a, const SlotAssign& b) {
        if (a.sid < b.sid) return true;
        if (b.sid < a.sid) return false;
        return a.slot < b.slot;
      });

      // -----------------
      // Phase 3: second sweep + conditional range sampling
      // -----------------
      auto __ = phases ? phases->Scoped("phase3_sweep") : PhaseRecorder::ScopedPhase(nullptr, "");

      out->pairs.resize(static_cast<usize>(t));

      // Start from empty active sets.
      kd_r_.ResetToEmpty();
      kd_s_.ResetToEmpty();

      usize ptr = 0;
      for (usize ev_i = 0; ev_i < events_.size(); ++ev_i) {
        const join::Event& e = events_[ev_i];

        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) kd_r_.Deactivate(static_cast<u32>(e.index));
          else kd_s_.Deactivate(static_cast<u32>(e.index));
          continue;
        }

        const u32 sid = start_id_of_event_[ev_i];
        SJS_DASSERT(sid != kInvalidStartId);

        // Slots for this START.
        usize end = ptr;
        while (end < asg.size() && asg[end].sid == sid) ++end;
        const u32 k = static_cast<u32>(end - ptr);

        if (k > 0) {
          // By construction, probability of choosing sid with w_e==0 is 0,
          // but keep a defensive guard.
          if (w_total_[static_cast<usize>(sid)] == 0) {
            if (err) *err = "KDTreeSamplingBaseline::Sample: sampled a START event with w_e==0 (unexpected)";
            return false;
          }

          const BoxT& q = (e.side == join::Side::R)
                              ? ds_->R.boxes[static_cast<usize>(e.index)]
                              : ds_->S.boxes[static_cast<usize>(e.index)];

          detail::RankBox<K> range;
          if (!embed_.MakeQueryRange(q, &range)) {
            if (err) *err = "KDTreeSamplingBaseline::Sample: empty query range for START that has slots";
            return false;
          }

          std::vector<u32> picked;
          picked.reserve(static_cast<usize>(k));

          if (e.side == join::Side::R) {
            if (!kd_s_.Sample(range, k, rng, &picked, err)) {
              if (err && err->empty()) *err = "KDTreeSamplingBaseline::Sample: kd_s_.Sample failed";
              return false;
            }
            for (u32 j = 0; j < k; ++j) {
              const u32 slot = asg[ptr + j].slot;
              const u32 s_idx = picked[static_cast<usize>(j)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(e.index)),
                                                          ds_->S.GetId(static_cast<usize>(s_idx))};
            }
          } else {
            if (!kd_r_.Sample(range, k, rng, &picked, err)) {
              if (err && err->empty()) *err = "KDTreeSamplingBaseline::Sample: kd_r_.Sample failed";
              return false;
            }
            for (u32 j = 0; j < k; ++j) {
              const u32 slot = asg[ptr + j].slot;
              const u32 r_idx = picked[static_cast<usize>(j)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(r_idx)),
                                                          ds_->S.GetId(static_cast<usize>(e.index))};
            }
          }

          ptr = end;
        }

        // Activate current.
        if (e.side == join::Side::R) kd_r_.Activate(static_cast<u32>(e.index));
        else kd_s_.Activate(static_cast<u32>(e.index));
      }

      if (ptr != asg.size()) {
        if (err) *err = "KDTreeSamplingBaseline::Sample: internal error (ptr != asg.size())";
        return false;
      }

      return true;
    }
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !ds_) {
      if (err) *err = "KDTreeSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    return std::make_unique<detail::KDJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                              join::SideTieBreak::RBeforeS);
  }

 private:
  static constexpr u32 kInvalidStartId = std::numeric_limits<u32>::max();

  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  // Preprocessing.
  std::vector<join::Event> events_;
  std::vector<u32> start_id_of_event_;

  detail::GlobalRankEmbedding<Dim, T> embed_;
  detail::ActiveKDTree<K> kd_r_;
  detail::ActiveKDTree<K> kd_s_;

  // Phase-1 results.
  bool weights_valid_ = false;
  u64 W_ = 0;
  std::vector<u64> w_total_;  // size = number of START events
};

}  // namespace kd_tree
}  // namespace baselines
}  // namespace sjs
