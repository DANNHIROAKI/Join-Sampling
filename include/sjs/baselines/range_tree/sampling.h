#pragma once
// sjs/baselines/range_tree/sampling.h
//
// Plane Sweep + Dynamic 2D Range-Tree baseline (Variant::Sampling).
//
// This implements the algorithm in "RangeTree Baseline.md" for the
// currently-supported experimental setting (Dim=2 rectangles):
//   - Sweep on axis 0 (x).
//   - Embed rectangles on the remaining axes (here only y) into rank-space
//     points p(r) = (rank(Ly(r)), rank(Ry(r))) in K=2 dims.
//   - Maintain two dynamic 2D range trees (one for R active set, one for S).
//   - Phase 1: for each START event e (query q), compute w_e = |K_e| as an
//     orthogonal range count on the opposite tree.
//   - Phase 2: build alias table over {w_e}, draw t START events (with
//     replacement), and assign sample slots to events.
//   - Phase 3: second sweep; for each START with t_e>0, sample t_e i.i.d.
//     uniform points from K_e using the dynamic range tree.
//
// Notes
// -----
// * Geometry uses half-open boxes [lo, hi) (consistent with sjs::Box).
// * Rank embedding uses a global (value, gid) sort to make strict inequalities
//   deterministic and to avoid floating-point "epsilon" issues.
// * This file also provides a deterministic join enumerator based on the same
//   sweep + range-tree machinery (used by EnumSampling & Adaptive variants).
// * The dynamic range tree implemented here is a baseline-grade structure:
//   it supports Activate/Deactivate, RangeCount, RangeReport, RangeSample.
//   Currently, we provide a specialization for 2D point space (K=2), which is
//   exactly what Dim=2 rectangle join embedding needs (K=2*(Dim-1)=2).
//
// Extensibility
// -------------
// The project is designed to grow to higher dimensions, but a full k-d range
// tree (k = 2*(Dim-1)) has very high memory overhead. For now, we keep this
// baseline restricted to Dim==2 (K==2). To extend later:
//   - Add ActiveRangeTree<K> specializations for K>2,
//   - Or swap to a different orthogonal sampler (e.g., kd-tree) for K>2.

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
namespace range_tree {

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

// --------------------------
// Global rank embedding (shared coordinate system for both relations)
// --------------------------

template <int Dim, class T>
struct GlobalRankEmbedding {
  static_assert(Dim >= 2, "GlobalRankEmbedding requires Dim >= 2");
  static constexpr int M = Dim - 1;       // number of non-sweep dimensions
  static constexpr int K = 2 * (Dim - 1); // embedding dimension

  struct Key {
    T value{};
    u32 gid{}; // global id for tie-break (unique across R and S)
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

      // Collect keys.
      for (u32 i = 0; i < n_r_; ++i) {
        const u32 gid = 1U + i;
        const auto& b = R.boxes[static_cast<usize>(i)];
        for (int ax = 1; ax < Dim; ++ax) {
          const int j = ax - 1;
          A_[j].push_back(Key{b.lo.v[ax], gid});
          B_[j].push_back(Key{b.hi.v[ax], gid});
        }
      }
      for (u32 i = 0; i < n_s_; ++i) {
        const u32 gid = 1U + n_r_ + i;
        const auto& b = S.boxes[static_cast<usize>(i)];
        for (int ax = 1; ax < Dim; ++ax) {
          const int j = ax - 1;
          A_[j].push_back(Key{b.lo.v[ax], gid});
          B_[j].push_back(Key{b.hi.v[ax], gid});
        }
      }

      // Sort each axis' rank arrays.
      for (int i = 0; i < M; ++i) {
        std::sort(A_[i].begin(), A_[i].end(), less);
        std::sort(B_[i].begin(), B_[i].end(), less);
        SJS_DASSERT(A_[i].size() == nt);
        SJS_DASSERT(B_[i].size() == nt);

        // Fill coordinates for relation points.
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
// Dynamic 2D Range Tree on rank-space points
// --------------------------
//
// We implement a classic static range-tree skeleton on x-order (dimension 0),
// where each node stores all points in its x-interval as a y-sorted list
// (dimension 1). We then add a Fenwick tree (BIT) per node to support dynamic
// activation/deactivation of points:
//   - Activate(p): update O(log n) nodes on the leaf-to-root path.
//   - Count([xl,xh) x [yl,yh)): decompose x-range into O(log n) nodes, and
//     query each node's BIT for y-range sum via two lower_bounds + BIT sum.
//   - Sample: use the same decomposition and BIT-select (order statistics).
//
// This is sufficient for the RangeTree baseline when Dim==2 rectangles,
// because embedding dimension K=2*(Dim-1)=2.

template <int K>
class ActiveRangeTree;

template <>
class ActiveRangeTree<2> {
 public:
  static constexpr int K = 2;
  using PointT = RankPoint<K>;
  using BoxT = RankBox<K>;

  struct QueryStats {
    u64 x_nodes = 0;     // number of canonical x-nodes touched
    u64 bit_queries = 0; // number of BIT range queries
    u64 reported = 0;    // number of reported points
  };

  ActiveRangeTree() = default;

  void Clear() {
    built_ = false;
    n_ = 0;
    p_ = 0;
    points_.clear();
    x_order_.clear();
    x_pos_.clear();
    node_.clear();
    y_pool_.clear();
    bit_pool_.clear();
    active_.clear();
    active_count_ = 0;
    active_list_.clear();
  }

  bool Built() const noexcept { return built_; }
  u32 Size() const noexcept { return n_; }
  u32 ActiveCount() const noexcept { return active_count_; }

  // Build the static skeleton on all points (points are copied).
  bool Build(Span<const PointT> points, std::string* err = nullptr) {
    Clear();
    const usize n64 = points.size();
    if (n64 > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "ActiveRangeTree<2>::Build: n too large for u32";
      return false;
    }
    n_ = static_cast<u32>(n64);
    points_.assign(points.begin(), points.end());
    if (n_ == 0) {
      built_ = true;
      p_ = 1;
      node_.assign(2, Node{});
      return true;
    }

    // x_order_ sorts point indices by x (dim0), tie by id.
    x_order_.resize(n_);
    for (u32 i = 0; i < n_; ++i) x_order_[static_cast<usize>(i)] = i;
    std::sort(x_order_.begin(), x_order_.end(), [&](u32 a, u32 b) {
      const u32 ax = points_[static_cast<usize>(a)].v[0];
      const u32 bx = points_[static_cast<usize>(b)].v[0];
      if (ax < bx) return true;
      if (bx < ax) return false;
      return a < b;
    });
    x_pos_.assign(n_, 0);
    for (u32 pos = 0; pos < n_; ++pos) {
      x_pos_[static_cast<usize>(x_order_[static_cast<usize>(pos)])] = pos;
    }

    // Build full segment tree over n leaves.
    p_ = 1;
    while (p_ < n_) p_ <<= 1;
    node_.assign(static_cast<usize>(2 * p_), Node{});

    // Estimate pool sizes (rough): y_pool ~ n * (log2(p_)+1).
    // This is a baseline, so we keep it simple and rely on vector growth.

    // Leaves.
    for (u32 i = 0; i < p_; ++i) {
      const u32 nid = p_ + i;
      if (i < n_) {
        const u32 obj = x_order_[static_cast<usize>(i)];
        AppendYList(nid, &obj, 1);
      } else {
        AppendYList(nid, nullptr, 0);
      }
    }

    // Internal nodes: merge children's y lists.
    std::vector<u32> tmp;
    for (u32 nid = p_ - 1; nid >= 1; --nid) {
      const u32 l = nid * 2;
      const u32 r = nid * 2 + 1;
      const Span<const u32> yl = YList(l);
      const Span<const u32> yr = YList(r);
      tmp.clear();
      tmp.reserve(yl.size() + yr.size());
      std::merge(yl.begin(), yl.end(), yr.begin(), yr.end(), std::back_inserter(tmp),
                 [&](u32 a, u32 b) { return LessY(a, b); });
      AppendYList(nid, tmp.data(), static_cast<u32>(tmp.size()));
      if (nid == 1) break;  // u32 underflow guard
    }

    // Active flags start empty.
    active_.assign(n_, 0);
    active_count_ = 0;
    active_list_.clear();

    built_ = true;
    return true;
  }

  bool Build(const std::vector<PointT>& points, std::string* err = nullptr) {
    return Build(Span<const PointT>(points.data(), points.size()), err);
  }

  // Reset the structure to the empty active set.
  // This is cheap if the set is already empty (common after a full sweep).
  void ResetToEmpty() {
    if (!built_) return;
    if (active_count_ == 0) {
      active_list_.clear();
      return;
    }
    // Fallback: deactivate anything still active.
    for (u32 obj : active_list_) {
      if (obj < n_ && active_[static_cast<usize>(obj)]) {
        Deactivate(obj);
      }
    }
    active_list_.clear();
    SJS_DASSERT(active_count_ == 0);
  }

  bool Contains(u32 obj) const noexcept {
    if (obj >= n_) return false;
    return active_[static_cast<usize>(obj)] != 0;
  }

  void Activate(u32 obj) {
    SJS_DASSERT(built_);
    SJS_DASSERT(obj < n_);
    if (active_[static_cast<usize>(obj)]) return;
    active_[static_cast<usize>(obj)] = 1;
    ++active_count_;
    active_list_.push_back(obj);

    const u32 xpos = x_pos_[static_cast<usize>(obj)];
    u32 nid = p_ + xpos;
    while (nid >= 1) {
      const u32 pos_y = FindPosInNode(nid, obj);
      BitAdd(nid, pos_y, +1);
      if (nid == 1) break;
      nid >>= 1;
    }
  }

  void Deactivate(u32 obj) {
    SJS_DASSERT(built_);
    SJS_DASSERT(obj < n_);
    if (!active_[static_cast<usize>(obj)]) return;
    active_[static_cast<usize>(obj)] = 0;
    SJS_DASSERT(active_count_ > 0);
    --active_count_;

    const u32 xpos = x_pos_[static_cast<usize>(obj)];
    u32 nid = p_ + xpos;
    while (nid >= 1) {
      const u32 pos_y = FindPosInNode(nid, obj);
      BitAdd(nid, pos_y, -1);
      if (nid == 1) break;
      nid >>= 1;
    }
  }

  // Count active points in the half-open query box.
  u64 Count(const BoxT& q, QueryStats* st = nullptr) const {
    if (!built_) return 0;
    if (n_ == 0) return 0;
    if (q.IsEmpty()) return 0;

    const u32 L = LowerBoundX(q.lo[0]);
    const u32 R = LowerBoundX(q.hi[0]);
    if (L >= R) return 0;

    std::vector<u32> nodes;
    CollectXNodes(L, R, &nodes);
    if (st) {
      st->x_nodes += static_cast<u64>(nodes.size());
    }

    u64 total = 0;
    for (u32 nid : nodes) {
      const u32 len = node_[static_cast<usize>(nid)].y_len;
      if (len == 0) continue;
      const u32 a = LowerBoundY(nid, q.lo[1]);
      const u32 b = LowerBoundY(nid, q.hi[1]);
      if (a >= b) continue;
      const u32 c = BitRangeSum(nid, a, b);
      total += static_cast<u64>(c);
      if (st) st->bit_queries += 1;
    }
    return total;
  }

  // Report (append) active point indices in the half-open query box.
  // The order is deterministic (x-canonical nodes left-to-right, and within
  // each node, increasing y-order).
  void Report(const BoxT& q, std::vector<u32>* out, QueryStats* st = nullptr) const {
    if (!out) return;
    if (!built_) return;
    if (n_ == 0) return;
    if (q.IsEmpty()) return;

    const u32 L = LowerBoundX(q.lo[0]);
    const u32 R = LowerBoundX(q.hi[0]);
    if (L >= R) return;

    std::vector<u32> nodes;
    CollectXNodes(L, R, &nodes);
    if (st) st->x_nodes += static_cast<u64>(nodes.size());

    for (u32 nid : nodes) {
      const u32 len = node_[static_cast<usize>(nid)].y_len;
      if (len == 0) continue;
      const u32 a = LowerBoundY(nid, q.lo[1]);
      const u32 b = LowerBoundY(nid, q.hi[1]);
      if (a >= b) continue;
      const u32 c = BitRangeSum(nid, a, b);
      if (st) st->bit_queries += 1;
      if (c == 0) continue;

      const u32 base = BitPrefixSum(nid, a);
      const Span<const u32> yl = YList(nid);
      for (u32 t = 0; t < c; ++t) {
        const u32 target = base + t + 1;  // 1-indexed order statistic
        const u32 pos = BitFindByOrder(nid, target);
        if (pos >= b) break;
        out->push_back(yl[static_cast<usize>(pos)]);
      }
      if (st) st->reported += static_cast<u64>(c);
    }
  }

  // Sample k i.i.d. active point indices uniformly from the query box.
  // Returns false only if k>0 and the query set is empty.
  bool Sample(const BoxT& q, u32 k, Rng* rng, std::vector<u32>* out, std::string* err = nullptr) const {
    if (!out) {
      if (err) *err = "ActiveRangeTree<2>::Sample: out is null";
      return false;
    }
    if (k == 0) return true;
    if (!rng) {
      if (err) *err = "ActiveRangeTree<2>::Sample: rng is null";
      return false;
    }
    if (!built_ || n_ == 0 || q.IsEmpty()) {
      if (err) *err = "ActiveRangeTree<2>::Sample: empty query";
      return false;
    }

    const u32 L = LowerBoundX(q.lo[0]);
    const u32 R = LowerBoundX(q.hi[0]);
    if (L >= R) {
      if (err) *err = "ActiveRangeTree<2>::Sample: empty x-range";
      return false;
    }

    // Canonical x decomposition.
    std::vector<u32> nodes;
    CollectXNodes(L, R, &nodes);

    struct Slice {
      u32 nid;
      u32 a;
      u32 b;
      u32 c;  // active count in [a,b)
    };
    std::vector<Slice> slices;
    slices.reserve(nodes.size());
    std::vector<u64> pref;  // prefix sums of counts
    pref.reserve(nodes.size() + 1);
    pref.push_back(0);

    for (u32 nid : nodes) {
      const u32 len = node_[static_cast<usize>(nid)].y_len;
      if (len == 0) continue;
      const u32 a = LowerBoundY(nid, q.lo[1]);
      const u32 b = LowerBoundY(nid, q.hi[1]);
      if (a >= b) continue;
      const u32 c = BitRangeSum(nid, a, b);
      if (c == 0) continue;
      slices.push_back(Slice{nid, a, b, c});
      pref.push_back(pref.back() + static_cast<u64>(c));
    }

    const u64 total = pref.back();
    if (total == 0) {
      if (err) *err = "ActiveRangeTree<2>::Sample: query set is empty";
      return false;
    }

    out->clear();
    out->reserve(static_cast<usize>(k));

    for (u32 i = 0; i < k; ++i) {
      const u64 r = rng->UniformU64(total);  // 0..total-1
      // Find slice index s.t. pref[j] <= r < pref[j+1].
      const auto it = std::upper_bound(pref.begin() + 1, pref.end(), r);
      const usize j = static_cast<usize>(std::distance(pref.begin(), it) - 1);
      SJS_DASSERT(j < slices.size());
      const Slice sl = slices[j];

      const u32 offset = static_cast<u32>(r - pref[j]);
      const u32 base = BitPrefixSum(sl.nid, sl.a);
      const u32 target = base + offset + 1;
      const u32 pos = BitFindByOrder(sl.nid, target);
      // Defensive: due to inconsistencies, pos might fall outside [a,b).
      if (pos < sl.a || pos >= sl.b) {
        if (err) *err = "ActiveRangeTree<2>::Sample: BIT select out of range (inconsistent state)";
        return false;
      }

      const Span<const u32> yl = YList(sl.nid);
      out->push_back(yl[static_cast<usize>(pos)]);
    }
    return true;
  }

 private:
  static constexpr u32 kNull = std::numeric_limits<u32>::max();

  struct Node {
    u32 y_off = 0;
    u32 y_len = 0;
    u32 bit_off = kNull;  // offset in bit_pool_ (0 means valid), kNull means no BIT
  };

  // Comparator: y (dim1) ascending, then id.
  bool LessY(u32 a, u32 b) const noexcept {
    const u32 ay = points_[static_cast<usize>(a)].v[1];
    const u32 by = points_[static_cast<usize>(b)].v[1];
    if (ay < by) return true;
    if (by < ay) return false;
    return a < b;
  }

  Span<const u32> YList(u32 nid) const noexcept {
    const Node& nd = node_[static_cast<usize>(nid)];
    return Span<const u32>(y_pool_.data() + nd.y_off, nd.y_len);
  }

  Span<u32> BitArr(u32 nid) noexcept {
    Node& nd = node_[static_cast<usize>(nid)];
    if (nd.y_len == 0 || nd.bit_off == kNull) return Span<u32>();
    // BIT length is y_len + 1.
    return Span<u32>(bit_pool_.data() + nd.bit_off, static_cast<usize>(nd.y_len) + 1);
  }

  Span<const u32> BitArr(u32 nid) const noexcept {
    const Node& nd = node_[static_cast<usize>(nid)];
    if (nd.y_len == 0 || nd.bit_off == kNull) return Span<const u32>();
    return Span<const u32>(bit_pool_.data() + nd.bit_off, static_cast<usize>(nd.y_len) + 1);
  }

  void AppendYList(u32 nid, const u32* data, u32 len) {
    Node& nd = node_[static_cast<usize>(nid)];
    nd.y_off = static_cast<u32>(y_pool_.size());
    nd.y_len = len;
    if (len > 0 && data) {
      y_pool_.insert(y_pool_.end(), data, data + len);
      nd.bit_off = static_cast<u32>(bit_pool_.size());
      bit_pool_.insert(bit_pool_.end(), static_cast<usize>(len) + 1, 0U);
    } else {
      nd.bit_off = kNull;
    }
  }

  // Return lower_bound on x-order for x >= value.
  u32 LowerBoundX(u32 value) const noexcept {
    auto it = std::lower_bound(x_order_.begin(), x_order_.end(), value,
                               [&](u32 idx, u32 v) { return points_[static_cast<usize>(idx)].v[0] < v; });
    return static_cast<u32>(std::distance(x_order_.begin(), it));
  }

  // Return lower_bound on node's y-list for y >= value.
  u32 LowerBoundY(u32 nid, u32 value) const noexcept {
    const Span<const u32> yl = YList(nid);
    auto it = std::lower_bound(yl.begin(), yl.end(), value,
                               [&](u32 idx, u32 v) { return points_[static_cast<usize>(idx)].v[1] < v; });
    return static_cast<u32>(std::distance(yl.begin(), it));
  }

  // Find the exact position of obj inside nid's y-list (binary search).
  u32 FindPosInNode(u32 nid, u32 obj) const {
    const Span<const u32> yl = YList(nid);
    auto it = std::lower_bound(yl.begin(), yl.end(), obj, [&](u32 a, u32 b) { return LessY(a, b); });
    // The tree is built from all points, so obj must be present.
    SJS_DASSERT(it != yl.end() && *it == obj);
    if (it == yl.end() || *it != obj) {
      // In release builds, fall back to linear scan (should never happen).
      for (usize i = 0; i < yl.size(); ++i) {
        if (yl[i] == obj) return static_cast<u32>(i);
      }
      return 0;
    }
    return static_cast<u32>(std::distance(yl.begin(), it));
  }

  // Collect canonical segment-tree nodes covering [L,R) in x-order index space.
  // Output order is left-to-right for determinism.
  void CollectXNodes(u32 L, u32 R, std::vector<u32>* out) const {
    SJS_DASSERT(out);
    out->clear();
    if (L >= R) return;
    u32 l = L + p_;
    u32 r = R + p_;
    std::vector<u32> left;
    std::vector<u32> right;
    left.reserve(64);
    right.reserve(64);
    while (l < r) {
      if (l & 1U) left.push_back(l++);
      if (r & 1U) right.push_back(--r);
      l >>= 1;
      r >>= 1;
    }
    std::reverse(right.begin(), right.end());
    out->reserve(left.size() + right.size());
    out->insert(out->end(), left.begin(), left.end());
    out->insert(out->end(), right.begin(), right.end());
  }

  // BIT primitives.
  void BitAdd(u32 nid, u32 pos0, int delta) {
    Node& nd = node_[static_cast<usize>(nid)];
    if (nd.y_len == 0 || nd.bit_off == kNull) return;
    SJS_DASSERT(pos0 < nd.y_len);
    const u32 len = nd.y_len;
    u32 i = pos0 + 1;
    u32* bit = bit_pool_.data() + nd.bit_off;
    while (i <= len) {
      if (delta > 0) {
        bit[i] += static_cast<u32>(delta);
      } else {
        bit[i] -= static_cast<u32>(-delta);
      }
      i += (i & (~i + 1));
    }
  }

  u32 BitPrefixSum(u32 nid, u32 pos_exclusive) const {
    const Node& nd = node_[static_cast<usize>(nid)];
    if (nd.y_len == 0 || nd.bit_off == kNull) return 0;
    const u32 len = nd.y_len;
    if (pos_exclusive > len) pos_exclusive = len;
    const u32* bit = bit_pool_.data() + nd.bit_off;
    u32 s = 0;
    u32 i = pos_exclusive;
    while (i > 0) {
      s += bit[i];
      i -= (i & (~i + 1));
    }
    return s;
  }

  u32 BitRangeSum(u32 nid, u32 a, u32 b) const {
    if (a >= b) return 0;
    return BitPrefixSum(nid, b) - BitPrefixSum(nid, a);
  }

  // Find smallest position pos (0-based) such that prefix(pos+1) >= target.
  // target must be in [1, total].
  u32 BitFindByOrder(u32 nid, u32 target) const {
    const Node& nd = node_[static_cast<usize>(nid)];
    const u32 len = nd.y_len;
    SJS_DASSERT(len > 0);
    const u32* bit = bit_pool_.data() + nd.bit_off;

    // Largest power of two >= len.
    u32 bitmask = 1;
    while (bitmask <= len) bitmask <<= 1;

    u32 idx = 0;
    for (u32 step = bitmask >> 1; step != 0; step >>= 1) {
      const u32 next = idx + step;
      if (next <= len && bit[next] < target) {
        idx = next;
        target -= bit[next];
      }
    }
    // idx is the largest with prefix(idx) < original target. Answer is idx+1.
    const u32 ans1 = idx + 1;
    SJS_DASSERT(ans1 >= 1 && ans1 <= len);
    return ans1 - 1;
  }

  bool built_ = false;
  u32 n_ = 0;  // number of points
  u32 p_ = 0;  // segment-tree base (power of two)

  std::vector<PointT> points_;
  std::vector<u32> x_order_;  // point indices sorted by x
  std::vector<u32> x_pos_;    // point index -> x_order position

  std::vector<Node> node_;     // size 2*p_
  std::vector<u32> y_pool_;    // concatenated y-lists
  std::vector<u32> bit_pool_;  // concatenated BIT arrays

  std::vector<u8> active_;
  u32 active_count_ = 0;
  std::vector<u32> active_list_;  // points activated since last ResetToEmpty()
};

// Default template: not implemented (yet).
template <int K>
class ActiveRangeTree {
  static_assert(K == 2, "ActiveRangeTree<K> is currently implemented only for K==2");
};

// --------------------------
// Deterministic enumerator (sweep + ActiveRangeTree)
// --------------------------

template <int Dim, class T>
class RangeJoinEnumerator final : public baselines::IJoinEnumerator {
 public:
  static_assert(Dim == 2, "RangeJoinEnumerator is currently implemented for Dim==2 only");
  using BoxT = Box<Dim, T>;
  static constexpr int K = 2 * (Dim - 1);

  RangeJoinEnumerator(const Relation<Dim, T>* rel_r,
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
    rt_r_.Clear();
    rt_s_.Clear();

    if (!R_ || !S_) return;

    // Build embedding and range trees (owned by this enumerator).
    std::vector<RankPoint<K>> pts_r;
    std::vector<RankPoint<K>> pts_s;
    if (!embed_.Build(*R_, *S_, &pts_r, &pts_s, /*phases=*/nullptr, /*err=*/nullptr)) {
      return;
    }
    (void)rt_r_.Build(pts_r);
    (void)rt_s_.Build(pts_s);
    std::vector<RankPoint<K>>().swap(pts_r);
    std::vector<RankPoint<K>>().swap(pts_s);

    // Events.
    events_ = join::BuildSweepEvents<Dim, T>(*R_, *S_, axis_, side_order_);
    stats_.num_events = static_cast<u64>(events_.size());

    // Force empty active sets.
    rt_r_.ResetToEmpty();
    rt_s_.ResetToEmpty();
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

        // Finish this START: activate current.
        if (cur_side_ == join::Side::R) {
          rt_r_.Activate(cur_index_);
          stats_.active_max_r = std::max(stats_.active_max_r, static_cast<u64>(rt_r_.ActiveCount()));
        } else {
          rt_s_.Activate(cur_index_);
          stats_.active_max_s = std::max(stats_.active_max_s, static_cast<u64>(rt_s_.ActiveCount()));
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
          rt_r_.Deactivate(static_cast<u32>(e.index));
        } else {
          rt_s_.Deactivate(static_cast<u32>(e.index));
        }
        continue;
      }

      // START: query opposite active set.
      cur_side_ = e.side;
      cur_index_ = static_cast<u32>(e.index);
      cur_id_ = e.id;

      const BoxT& q = (cur_side_ == join::Side::R)
                          ? R_->boxes[static_cast<usize>(cur_index_)]
                          : S_->boxes[static_cast<usize>(cur_index_)];

      RankBox<K> range;
      cur_candidates_.clear();

      typename ActiveRangeTree<K>::QueryStats qst;
      if (embed_.MakeQueryRange(q, &range)) {
        if (cur_side_ == join::Side::R) {
          rt_s_.Report(range, &cur_candidates_, &qst);
        } else {
          rt_r_.Report(range, &cur_candidates_, &qst);
        }
      }

      // For range-tree, candidate_checks is best interpreted as the number of
      // produced candidates (since BIT avoids explicit predicate checks).
      stats_.candidate_checks += static_cast<u64>(cur_candidates_.size());

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

  join::JoinStats stats_;
  std::vector<join::Event> events_;

  GlobalRankEmbedding<Dim, T> embed_;
  ActiveRangeTree<K> rt_r_;
  ActiveRangeTree<K> rt_s_;

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
// RangeTreeSamplingBaseline
// --------------------------

template <int Dim, class T = Scalar>
class RangeTreeSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim == 2, "RangeTreeSamplingBaseline is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;
  static constexpr int K = 2 * (Dim - 1);

  Method method() const noexcept override { return Method::RangeTree; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "range_tree_sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;

    events_.clear();
    start_id_of_event_.clear();
    w_total_.clear();

    embed_.Reset();
    rt_r_.Clear();
    rt_s_.Clear();
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

    // Embedding + range trees.
    {
      std::vector<detail::RankPoint<K>> pts_r;
      std::vector<detail::RankPoint<K>> pts_s;

      if (!embed_.Build(ds.R, ds.S, &pts_r, &pts_s, phases, err)) return false;

      auto _ = phases ? phases->Scoped("build_range_trees") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!rt_r_.Build(pts_r, err)) return false;
      if (!rt_s_.Build(pts_s, err)) return false;

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
    (void)rng;  // deterministic
    if (!built_ || !ds_) {
      if (err) *err = "RangeTreeSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RangeTreeSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    if (ds_->R.Size() == 0 || ds_->S.Size() == 0) {
      W_ = 0;
      weights_valid_ = true;
      std::fill(w_total_.begin(), w_total_.end(), 0ULL);
      *out = MakeExactCount(0);
      return true;
    }

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    rt_r_.ResetToEmpty();
    rt_s_.ResetToEmpty();

    u64 W = 0;
    for (usize ev_i = 0; ev_i < events_.size(); ++ev_i) {
      const join::Event& e = events_[ev_i];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) rt_r_.Deactivate(static_cast<u32>(e.index));
        else rt_s_.Deactivate(static_cast<u32>(e.index));
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
        if (e.side == join::Side::R) w = rt_s_.Count(range);
        else w = rt_r_.Count(range);
      }

      w_total_[static_cast<usize>(sid)] = w;
      if (w > 0) {
        if (W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "RangeTreeSamplingBaseline::Count: |J| overflowed u64";
          return false;
        }
        W += w;
      }

      // Activate current.
      if (e.side == join::Side::R) rt_r_.Activate(static_cast<u32>(e.index));
      else rt_s_.Activate(static_cast<u32>(e.index));
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
      if (err) *err = "RangeTreeSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "RangeTreeSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "RangeTreeSamplingBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t64 = cfg.run.t;
    if (t64 == 0) return true;
    if (t64 > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RangeTreeSamplingBaseline::Sample: t too large for u32 slots";
      return false;
    }
    const u32 t = static_cast<u32>(t64);

    // Ensure weights.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }
    if (W_ == 0) {
      return true;  // empty join
    }

    // -----------------
    // Phase 2: alias + slot assignment
    // -----------------
    {
      auto _ = phases ? phases->Scoped("phase2_alias") : PhaseRecorder::ScopedPhase(nullptr, "");

      sampling::AliasTable alias;
      if (!alias.BuildFromU64(Span<const u64>(w_total_.data(), w_total_.size()), err)) return false;

      struct SlotAssign {
        u32 sid;
        u32 slot;
      };

      std::vector<SlotAssign> asg;
      asg.reserve(static_cast<usize>(t));
      for (u32 j = 0; j < t; ++j) {
        const u32 sid = static_cast<u32>(alias.Sample(rng));
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

      rt_r_.ResetToEmpty();
      rt_s_.ResetToEmpty();

      usize ptr = 0;
      for (usize ev_i = 0; ev_i < events_.size(); ++ev_i) {
        const join::Event& e = events_[ev_i];

        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) rt_r_.Deactivate(static_cast<u32>(e.index));
          else rt_s_.Deactivate(static_cast<u32>(e.index));
          continue;
        }

        const u32 sid = start_id_of_event_[ev_i];
        SJS_DASSERT(sid != kInvalidStartId);

        // Slots for this START.
        usize end = ptr;
        while (end < asg.size() && asg[end].sid == sid) ++end;
        const u32 k = static_cast<u32>(end - ptr);

        if (k > 0) {
          if (w_total_[static_cast<usize>(sid)] == 0) {
            if (err) *err = "RangeTreeSamplingBaseline::Sample: sampled a START with w_e==0 (unexpected)";
            return false;
          }

          const BoxT& q = (e.side == join::Side::R)
                              ? ds_->R.boxes[static_cast<usize>(e.index)]
                              : ds_->S.boxes[static_cast<usize>(e.index)];

          detail::RankBox<K> range;
          if (!embed_.MakeQueryRange(q, &range)) {
            if (err) *err = "RangeTreeSamplingBaseline::Sample: empty query range for START that has slots";
            return false;
          }

          std::vector<u32> picked;
          picked.reserve(static_cast<usize>(k));

          if (e.side == join::Side::R) {
            if (!rt_s_.Sample(range, k, rng, &picked, err)) return false;
            for (u32 j = 0; j < k; ++j) {
              const u32 slot = asg[ptr + j].slot;
              const u32 s_idx = picked[static_cast<usize>(j)];
              out->pairs[static_cast<usize>(slot)] = PairId{ds_->R.GetId(static_cast<usize>(e.index)),
                                                          ds_->S.GetId(static_cast<usize>(s_idx))};
            }
          } else {
            if (!rt_r_.Sample(range, k, rng, &picked, err)) return false;
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
        if (e.side == join::Side::R) rt_r_.Activate(static_cast<u32>(e.index));
        else rt_s_.Activate(static_cast<u32>(e.index));
      }

      if (ptr != asg.size()) {
        if (err) *err = "RangeTreeSamplingBaseline::Sample: internal error (ptr != asg.size())";
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
      if (err) *err = "RangeTreeSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    return std::make_unique<detail::RangeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, /*axis=*/0,
                                                                 join::SideTieBreak::RBeforeS);
  }

 private:
  static constexpr u32 kInvalidStartId = std::numeric_limits<u32>::max();

  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  std::vector<join::Event> events_;
  std::vector<u32> start_id_of_event_;

  detail::GlobalRankEmbedding<Dim, T> embed_;
  detail::ActiveRangeTree<K> rt_r_;
  detail::ActiveRangeTree<K> rt_s_;

  bool weights_valid_ = false;
  u64 W_ = 0;
  std::vector<u64> w_total_;
};

}  // namespace range_tree
}  // namespace baselines
}  // namespace sjs
