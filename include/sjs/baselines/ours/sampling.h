#pragma once
// sjs/baselines/ours/sampling.h
//
// Our method (Sampling variant) for uniform sampling from the spatial join.
//
// Implementation notes
// --------------------
// * The project currently focuses on 2D. This header implements an optimized
//   Dim=2 version, but is templated so you can extend to higher dimensions by
//   swapping in higher-dimensional pattern indices later.
// * Geometry uses half-open rectangles: [lo, hi) in each dimension.
// * We use a deterministic plane sweep on axis 0 (x). END events are processed
//   before START events at the same x (half-open), and START tie-break follows
//   join::SideTieBreak (default: R before S).
//
// Algorithm (2D)
// --------------
// Let q be the rectangle of a START event on one side; the opposite side has an
// active set A_x of rectangles whose x-interval contains q.lo.x.
// Since A_x implies overlap on the sweep axis, an (R,S) pair intersects iff its
// y-intervals intersect.
//
// For 1D y-interval intersection under half-open semantics:
//   [a0,a1) intersects [b0,b1) iff a0 < b1 AND b0 < a1.
// Partition by comparing lower endpoints (ties go to pattern A):
//   Pattern A: b0 <= a0 < b1
//   Pattern B: a0 < b0 < a1
// These two cases are disjoint and cover all intersections.
//
// Sampling scheme
// ---------------
// Phase 1: For each START event e with rectangle q, compute
//   w_e^A = |{ r in A_x : y_r contains y_q_lo }|
//   w_e^B = |{ r in A_x : y_r_lo in (y_q_lo, y_q_hi) }|
// and w_e = w_e^A + w_e^B. Let W = sum_e w_e = |J|.
// Phase 2: Build alias table on (w_e). For each sample slot j:
//   pick an event e with prob w_e / W
//   pick pattern g in {A,B} with prob w_e^g / w_e
//   record assignment (e,g,slot=j)
// Phase 3: Re-sweep. When reaching START event e for q, for each assigned slot
// in group g, sample uniformly from the corresponding set and output the pair.
//
// Correctness follows because each join pair is assigned to exactly one START
// event (the later START in sweep order) and exactly one pattern (A/B).

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/join/join_enumerator.h"
#include "sjs/join/sweep_events.h"
#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace ours {

namespace detail {

// --------------------------
// Fenwick tree (BIT) for u64 counts
// --------------------------
class FenwickU64 {
 public:
  FenwickU64() = default;

  void Init(u32 n) {
    n_ = n;
    bit_.assign(static_cast<usize>(n_) + 1, 0ULL);
  }

  void Clear() {
    n_ = 0;
    bit_.clear();
  }

  u32 Size() const noexcept { return n_; }

  // Add delta to index idx (0-based).
  void Add(u32 idx, i64 delta) {
    SJS_DASSERT(idx < n_);
    const i64 d = delta;
    u32 i = idx + 1;
    while (i <= n_) {
      // We store u64; allow negative delta (erase).
      const i64 cur = static_cast<i64>(bit_[i]);
      const i64 next = cur + d;
      SJS_DASSERT(next >= 0);
      bit_[i] = static_cast<u64>(next);
      i += (i & (~i + 1));
    }
  }

  // Sum of [0, idx) (idx is exclusive, 0..n).
  u64 SumPrefix(u32 idx_exclusive) const {
    SJS_DASSERT(idx_exclusive <= n_);
    u64 res = 0;
    u32 i = idx_exclusive;
    while (i > 0) {
      res += bit_[i];
      i -= (i & (~i + 1));
    }
    return res;
  }

  u64 SumRange(u32 l, u32 r) const {
    if (r <= l) return 0ULL;
    return SumPrefix(r) - SumPrefix(l);
  }

  u64 Total() const { return SumPrefix(n_); }

  // Find smallest idx such that SumPrefix(idx+1) >= target, where target is in [1, Total()].
  // Returns idx in [0, n_-1].
  u32 LowerBound(u64 target) const {
    SJS_DASSERT(target >= 1);
    SJS_DASSERT(target <= Total());

    // Largest power of two <= n_.
    u32 bit = 1;
    while ((bit << 1) <= n_) bit <<= 1;

    u32 idx = 0;
    u64 cur = 0;
    while (bit != 0) {
      const u32 next = idx + bit;
      if (next <= n_ && cur + bit_[next] < target) {
        idx = next;
        cur += bit_[next];
      }
      bit >>= 1;
    }
    // idx is the largest Fenwick index with prefix sum < target.
    // Convert to 0-based element index.
    return idx;
  }

 private:
  u32 n_{0};
  std::vector<u64> bit_;
};

// --------------------------
// Range-point structure for Pattern B (2D)
//
// Maintains active points y = lo_y(r) with multiplicities (multiple rectangles
// can share the same y). Supports:
//  - Count in [l, r) over ranks
//  - Sample uniformly from points in [l, r)
//  - Report (for debugging / small-join enumeration)
//
// We store handles per rank in a bucket vector and maintain counts in a BIT.
// Sampling selects a rank by order statistic (BIT lower_bound) then a random
// element from that rank's bucket.
// --------------------------
class RangePointIndex {
 public:
  void Init(u32 num_handles, u32 num_ranks) {
    n_handles_ = num_handles;
    m_ = num_ranks;
    bit_.Init(m_);
    buckets_.assign(static_cast<usize>(m_), {});
    pos_.assign(static_cast<usize>(n_handles_), -1);
  }

  void Clear() {
    n_handles_ = 0;
    m_ = 0;
    bit_.Clear();
    buckets_.clear();
    pos_.clear();
  }

  void Insert(u32 handle, u32 rank) {
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(rank < m_);
    SJS_DASSERT(pos_[handle] == -1);
    auto& b = buckets_[rank];
    pos_[handle] = static_cast<i32>(b.size());
    b.push_back(handle);
    bit_.Add(rank, +1);
  }

  void Erase(u32 handle, u32 rank) {
    SJS_DASSERT(handle < n_handles_);
    SJS_DASSERT(rank < m_);
    const i32 p = pos_[handle];
    SJS_DASSERT(p >= 0);
    auto& b = buckets_[rank];
    const usize up = static_cast<usize>(p);
    SJS_DASSERT(up < b.size());
    const u32 last = b.back();
    b[up] = last;
    pos_[last] = static_cast<i32>(up);
    b.pop_back();
    pos_[handle] = -1;
    bit_.Add(rank, -1);
  }

  u64 CountRange(u32 l, u32 r) const {
    if (r <= l) return 0ULL;
    l = std::min(l, m_);
    r = std::min(r, m_);
    if (r <= l) return 0ULL;
    return bit_.SumRange(l, r);
  }

  // Sample k handles uniformly from points with ranks in [l, r).
  // Returns false if the range is empty.
  bool SampleRange(u32 l, u32 r, u32 k, Rng* rng, std::vector<u32>* out) const {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(out != nullptr);
    out->clear();
    out->reserve(k);

    if (r <= l) return false;
    l = std::min(l, m_);
    r = std::min(r, m_);
    if (r <= l) return false;

    const u64 total = CountRange(l, r);
    if (total == 0) return false;

    const u64 prefix_l = bit_.SumPrefix(l);

    for (u32 i = 0; i < k; ++i) {
      const u64 pick = rng->UniformU64(total);  // [0,total)
      const u64 target = prefix_l + pick + 1;  // 1-based target in [prefix_l+1, prefix_l+total]
      const u32 rank = bit_.LowerBound(target);
      SJS_DASSERT(rank >= l && rank < r);
      const auto& b = buckets_[rank];
      SJS_DASSERT(!b.empty());
      const u32 idx = rng->UniformU32(static_cast<u32>(b.size()));
      out->push_back(b[idx]);
    }

    return true;
  }

  // Report all handles in [l, r) (deterministic order by increasing rank).
  void ReportRange(u32 l, u32 r, std::vector<u32>* out) const {
    SJS_DASSERT(out != nullptr);
    if (r <= l) return;
    l = std::min(l, m_);
    r = std::min(r, m_);
    for (u32 rank = l; rank < r; ++rank) {
      const auto& b = buckets_[rank];
      out->insert(out->end(), b.begin(), b.end());
    }
  }

 private:
  u32 n_handles_{0};
  u32 m_{0};
  FenwickU64 bit_;
  std::vector<std::vector<u32>> buckets_;  // per-rank bucket of handles
  std::vector<i32> pos_;                   // handle -> position in its bucket, -1 if inactive
};

// --------------------------
// Stabbing segment tree for Pattern A (2D)
//
// Maintains active y-intervals [lo, hi) represented on the discrete domain of
// query points (ranks). Supports:
//  - Count intervals containing a query point rank
//  - Sample uniformly from those intervals
//  - Report all intervals containing the query point
//
// Implementation:
//  - Iterative segment tree over ranks [0, P) where P is the next power-of-two.
//  - For each interval range [L, R) on ranks, we add the handle to O(log P)
//    canonical nodes (standard segment tree cover).
//  - For point queries, we traverse the path from the leaf to the root and
//    take the union of buckets along the path. Because canonical cover nodes are
//    disjoint, each interval appears in exactly one bucket on the path.
//  - To support O(1) deletions, each inserted bucket item stores a backref
//    (index into a per-handle placement list), and we update that backref on
//    swap-delete.
//
// Memory/performance trade-off:
//  - We avoid per-handle heap allocations by storing a fixed-size placement pool
//    sized (num_handles * max_refs_per_handle).
// --------------------------
class StabbingSegTree {
 public:
  struct Placement {
    u32 node = 0;  // node index in the implicit tree (1..2P-1)
    u32 pos = 0;   // position inside nodes_[node].items
  };

  struct Item {
    u32 handle = 0;
    u32 backref = 0;  // index in the handle's placement list
  };

  void Init(u32 num_handles, u32 num_points) {
    n_handles_ = num_handles;
    m_ = num_points;

    // Compute P = next power-of-two >= m_. If m_==0, keep P=1.
    p_ = 1;
    while (p_ < m_) p_ <<= 1;

    nodes_.assign(static_cast<usize>(2 * p_), Node{});

    // log2(p_)
    u32 logp = 0;
    for (u32 x = p_; x > 1; x >>= 1) ++logp;

    // Range decomposition can touch at most 2*log2(p_) nodes.
    max_refs_ = static_cast<u32>(2 * logp + 4);

    placement_size_.assign(static_cast<usize>(n_handles_), 0);
    placements_.assign(static_cast<usize>(n_handles_) * static_cast<usize>(max_refs_), Placement{});
  }

  void Clear() {
    n_handles_ = 0;
    m_ = 0;
    p_ = 0;
    max_refs_ = 0;
    nodes_.clear();
    placements_.clear();
    placement_size_.clear();
  }

  // Insert interval [L, R) on ranks (0 <= L <= R <= m_).
  void Insert(u32 handle, u32 L, u32 R) {
    SJS_DASSERT(handle < n_handles_);
    if (R <= L) return;  // empty on the discrete query domain
    L = std::min(L, m_);
    R = std::min(R, m_);
    if (R <= L) return;

    u32 l = L + p_;
    u32 r = R + p_;

    while (l < r) {
      if (l & 1U) {
        AddToNode(handle, l);
        ++l;
      }
      if (r & 1U) {
        --r;
        AddToNode(handle, r);
      }
      l >>= 1;
      r >>= 1;
    }
  }

  // Erase all placements of handle.
  void Erase(u32 handle) {
    SJS_DASSERT(handle < n_handles_);
    const u32 sz = placement_size_[handle];
    // Note: Even if sz==0, erase is valid (idempotent in release builds).
    for (u32 i = 0; i < sz; ++i) {
      const Placement p = placements_[PlacementIndex(handle, i)];
      RemoveFromNode(p.node, p.pos);
    }
    placement_size_[handle] = 0;
  }

  // Count intervals containing point rank q (0 <= q < m_).
  u64 Count(u32 q) const {
    if (m_ == 0) return 0ULL;
    if (q >= m_) return 0ULL;
    u64 total = 0;
    u32 idx = q + p_;
    while (idx > 0) {
      total += static_cast<u64>(nodes_[idx].items.size());
      idx >>= 1;
    }
    return total;
  }

  // Sample k handles uniformly from intervals containing q.
  // Returns false if empty.
  bool Sample(u32 q, u32 k, Rng* rng, std::vector<u32>* out) const {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(out != nullptr);
    out->clear();
    out->reserve(k);

    if (m_ == 0 || q >= m_) return false;

    // Collect non-empty buckets on the root-to-leaf path.
    std::vector<u32> path_nodes;
    std::vector<u64> weights;
    path_nodes.reserve(32);
    weights.reserve(32);

    u32 idx = q + p_;
    u64 total = 0;
    while (idx > 0) {
      const auto sz = static_cast<u64>(nodes_[idx].items.size());
      if (sz > 0) {
        path_nodes.push_back(idx);
        weights.push_back(sz);
        total += sz;
      }
      idx >>= 1;
    }

    if (total == 0) return false;

    sampling::AliasTable alias;
    (void)alias.BuildFromU64(Span<const u64>(weights));

    for (u32 i = 0; i < k; ++i) {
      const usize bi = alias.Sample(rng);
      const u32 node = path_nodes[bi];
      const auto& b = nodes_[node].items;
      SJS_DASSERT(!b.empty());
      const u32 pos = rng->UniformU32(static_cast<u32>(b.size()));
      out->push_back(b[pos].handle);
    }

    return true;
  }

  // Report all handles containing q (deterministic in path order).
  void Report(u32 q, std::vector<u32>* out) const {
    SJS_DASSERT(out != nullptr);
    if (m_ == 0 || q >= m_) return;

    u32 idx = q + p_;
    while (idx > 0) {
      const auto& b = nodes_[idx].items;
      for (const auto& it : b) out->push_back(it.handle);
      idx >>= 1;
    }
  }

 private:
  struct Node {
    std::vector<Item> items;
  };

  usize PlacementIndex(u32 handle, u32 backref) const {
    return static_cast<usize>(handle) * static_cast<usize>(max_refs_) + static_cast<usize>(backref);
  }

  void AddToNode(u32 handle, u32 node) {
    SJS_DASSERT(node < nodes_.size());
    u32& sz = placement_size_[handle];
    SJS_DASSERT(sz < max_refs_);

    auto& bucket = nodes_[node].items;
    const u32 pos = static_cast<u32>(bucket.size());

    placements_[PlacementIndex(handle, sz)] = Placement{node, pos};
    bucket.push_back(Item{handle, sz});
    ++sz;
  }

  void RemoveFromNode(u32 node, u32 pos) {
    SJS_DASSERT(node < nodes_.size());
    auto& bucket = nodes_[node].items;
    SJS_DASSERT(!bucket.empty());
    SJS_DASSERT(pos < bucket.size());

    const usize last_pos = bucket.size() - 1;
    if (static_cast<usize>(pos) != last_pos) {
      const Item swapped = bucket[last_pos];
      bucket[pos] = swapped;
      // Fix swapped handle's placement position.
      placements_[PlacementIndex(swapped.handle, swapped.backref)].pos = pos;
    }
    bucket.pop_back();
  }

  u32 n_handles_{0};
  u32 m_{0};
  u32 p_{0};

  u32 max_refs_{0};

  std::vector<Node> nodes_;               // size 2*p_
  std::vector<Placement> placements_;     // size n_handles_ * max_refs_
  std::vector<u32> placement_size_;       // per-handle placement count
};

// --------------------------
// Active index for one side (2D)
// --------------------------
class ActiveIndex2D {
 public:
  void Init(u32 num_handles, u32 num_ylo_ranks) {
    stab_.Init(num_handles, num_ylo_ranks);
    pts_.Init(num_handles, num_ylo_ranks);
    n_ = num_handles;
    m_ = num_ylo_ranks;
  }

  void Clear() {
    stab_.Clear();
    pts_.Clear();
    n_ = 0;
    m_ = 0;
  }

  void Insert(u32 handle, u32 ylo_rank, u32 yhi_lb_rank) {
    SJS_DASSERT(handle < n_);
    SJS_DASSERT(ylo_rank < m_);
    // yhi_lb_rank is in [0, m_]
    stab_.Insert(handle, ylo_rank, yhi_lb_rank);
    pts_.Insert(handle, ylo_rank);
  }

  void Erase(u32 handle, u32 ylo_rank) {
    SJS_DASSERT(handle < n_);
    SJS_DASSERT(ylo_rank < m_);
    stab_.Erase(handle);
    pts_.Erase(handle, ylo_rank);
  }

  // Pattern A: y_r contains y_q_lo.
  u64 CountA(u32 y_q_lo_rank) const { return stab_.Count(y_q_lo_rank); }

  // Pattern B: y_r_lo in (y_q_lo, y_q_hi).
  u64 CountB(u32 y_q_lo_rank, u32 y_q_hi_lb_rank) const {
    // Exclusive lower bound: since ranks are unique ylo values, > is +1.
    const u32 l = (y_q_lo_rank < m_) ? (y_q_lo_rank + 1) : m_;
    const u32 r = std::min(y_q_hi_lb_rank, m_);
    return pts_.CountRange(l, r);
  }

  bool SampleA(u32 y_q_lo_rank, u32 k, Rng* rng, std::vector<u32>* out) const {
    return stab_.Sample(y_q_lo_rank, k, rng, out);
  }

  bool SampleB(u32 y_q_lo_rank, u32 y_q_hi_lb_rank, u32 k, Rng* rng, std::vector<u32>* out) const {
    const u32 l = (y_q_lo_rank < m_) ? (y_q_lo_rank + 1) : m_;
    const u32 r = std::min(y_q_hi_lb_rank, m_);
    return pts_.SampleRange(l, r, k, rng, out);
  }

  void ReportA(u32 y_q_lo_rank, std::vector<u32>* out) const { stab_.Report(y_q_lo_rank, out); }

  void ReportB(u32 y_q_lo_rank, u32 y_q_hi_lb_rank, std::vector<u32>* out) const {
    const u32 l = (y_q_lo_rank < m_) ? (y_q_lo_rank + 1) : m_;
    const u32 r = std::min(y_q_hi_lb_rank, m_);
    pts_.ReportRange(l, r, out);
  }

 private:
  u32 n_{0};
  u32 m_{0};
  StabbingSegTree stab_;
  RangePointIndex pts_;
};

// Wrapper enumerator around join::PlaneSweepJoinStream so it satisfies baselines::IJoinEnumerator.
template <int Dim, class T>
class PlaneSweepEnumeratorWrapper final : public baselines::IJoinEnumerator {
 public:
  PlaneSweepEnumeratorWrapper(const Relation<Dim, T>& R,
                              const Relation<Dim, T>& S,
                              join::PlaneSweepOptions opt)
      : stream_(R, S, opt) {}

  void Reset() override { stream_.Reset(); }

  bool Next(PairId* out) override { return stream_.Next(out); }

  const join::JoinStats& Stats() const noexcept override { return stream_.Stats(); }

 private:
  join::PlaneSweepJoinStream<Dim, T> stream_;
};

}  // namespace detail

// --------------------------
// OursSamplingBaseline (Dim=2)
// --------------------------
template <int Dim, class T = Scalar>
class OursSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "ours_sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;

    // Keep allocated memory but drop contents.
    events_.clear();
    start_id_of_event_.clear();
    start_event_pos_.clear();

    y_coords_.clear();
    ylo_rank_r_.clear();
    yhi_lb_rank_r_.clear();
    ylo_rank_s_.clear();
    yhi_lb_rank_s_.clear();

    active_r_.Clear();
    active_s_.Clear();

    w_total_.clear();
    w_a_.clear();
    w_b_.clear();
  }

  bool Build(const DatasetT& ds,
             const Config& cfg,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    if constexpr (Dim != 2) {
      if (err) *err = "OursSamplingBaseline: currently only supports Dim=2";
      return false;
    }

    ds_ = &ds;
    built_ = false;
    weights_valid_ = false;
    W_ = 0;

    const usize nR = ds.R.Size();
    const usize nS = ds.S.Size();
    if (nR > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        nS > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "OursSamplingBaseline: relation size exceeds u32 capacity";
      return false;
    }

    // Build sweep events on axis 0.
    {
      auto scoped = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = join::BuildSweepEvents<Dim, T>(ds, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    // Map event position -> start-event id (dense over START events).
    start_id_of_event_.assign(events_.size(), -1);
    start_event_pos_.clear();
    start_event_pos_.reserve(ds.TotalSize());
    for (usize i = 0; i < events_.size(); ++i) {
      if (events_[i].kind == join::EventKind::Start) {
        start_id_of_event_[i] = static_cast<i32>(start_event_pos_.size());
        start_event_pos_.push_back(i);
      }
    }

    // Coordinate compression domain: all y-lower endpoints from both relations.
    {
      auto scoped = phases ? phases->Scoped("build_y_domain") : PhaseRecorder::ScopedPhase(nullptr, "");
      y_coords_.clear();
      y_coords_.reserve(nR + nS);
      for (const auto& b : ds.R.boxes) y_coords_.push_back(b.lo.v[1]);
      for (const auto& b : ds.S.boxes) y_coords_.push_back(b.lo.v[1]);
      std::sort(y_coords_.begin(), y_coords_.end());
      y_coords_.erase(std::unique(y_coords_.begin(), y_coords_.end()), y_coords_.end());
    }

    if (y_coords_.empty()) {
      if (err) *err = "OursSamplingBaseline: empty y-domain (no boxes?)";
      return false;
    }

    const u32 m = static_cast<u32>(y_coords_.size());

    auto lb_rank = [&](T v) -> u32 {
      const auto it = std::lower_bound(y_coords_.begin(), y_coords_.end(), v);
      return static_cast<u32>(std::distance(y_coords_.begin(), it));
    };

    // Precompute ranks per box for fast sweeps.
    {
      auto scoped = phases ? phases->Scoped("build_ranks") : PhaseRecorder::ScopedPhase(nullptr, "");
      ylo_rank_r_.resize(nR);
      yhi_lb_rank_r_.resize(nR);
      for (usize i = 0; i < nR; ++i) {
        const auto& b = ds.R.boxes[i];
        const u32 lo = lb_rank(b.lo.v[1]);
        // Debug: lo should be exact match because we built domain from lo's.
        SJS_DASSERT(lo < m && y_coords_[lo] == b.lo.v[1]);
        ylo_rank_r_[i] = lo;
        yhi_lb_rank_r_[i] = lb_rank(b.hi.v[1]);  // can be m
      }

      ylo_rank_s_.resize(nS);
      yhi_lb_rank_s_.resize(nS);
      for (usize i = 0; i < nS; ++i) {
        const auto& b = ds.S.boxes[i];
        const u32 lo = lb_rank(b.lo.v[1]);
        SJS_DASSERT(lo < m && y_coords_[lo] == b.lo.v[1]);
        ylo_rank_s_[i] = lo;
        yhi_lb_rank_s_[i] = lb_rank(b.hi.v[1]);
      }
    }

    // Init active indices.
    {
      auto scoped = phases ? phases->Scoped("build_active_index") : PhaseRecorder::ScopedPhase(nullptr, "");
      active_r_.Init(static_cast<u32>(nR), m);
      active_s_.Init(static_cast<u32>(nS), m);
    }

    // Allocate weights.
    const usize E = start_event_pos_.size();
    w_total_.assign(E, 0ULL);
    w_a_.assign(E, 0ULL);
    w_b_.assign(E, 0ULL);

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;  // Count is deterministic in our method.

    if (!built_ || ds_ == nullptr) {
      if (err) *err = "OursSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "OursSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    std::fill(w_a_.begin(), w_a_.end(), 0ULL);
    std::fill(w_b_.begin(), w_b_.end(), 0ULL);

    u64 W = 0;

    // #region agent log
    {
      std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
      if (log.is_open()) {
        log << "{\"id\":\"log_count_start\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:746\",\"message\":\"Count phase1 start\",\"data\":{\"events_size\":" << events_.size() << ",\"start_events\":" << start_event_pos_.size() << "},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"A\"}\n";
      }
    }
    // #endregion

    // Sweep state is stored inside active_*; since each box is inserted once and
    // erased once, the structures end empty.
    for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
      const auto& e = events_[ev_pos];

      if (e.kind == join::EventKind::End) {
        if (e.side == join::Side::R) {
          const u32 h = static_cast<u32>(e.index);
          active_r_.Erase(h, ylo_rank_r_[e.index]);
        } else {
          const u32 h = static_cast<u32>(e.index);
          active_s_.Erase(h, ylo_rank_s_[e.index]);
        }
        continue;
      }

      // START event.
      const i32 eid_i32 = start_id_of_event_[ev_pos];
      SJS_DASSERT(eid_i32 >= 0);
      const usize eid = static_cast<usize>(eid_i32);

      // #region agent log
      {
        std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
        if (log.is_open()) {
          log << "{\"id\":\"log_event_mapping\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:767\",\"message\":\"Event index mapping\",\"data\":{\"ev_pos\":" << ev_pos << ",\"eid_i32\":" << eid_i32 << ",\"eid\":" << eid << ",\"side\":\"" << (e.side == join::Side::R ? "R" : "S") << "\"},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"C\"}\n";
        }
      }
      // #endregion

      if (e.side == join::Side::R) {
        const u32 q_ylo = ylo_rank_r_[e.index];
        const u32 q_yhi = yhi_lb_rank_r_[e.index];

        const u64 wa = active_s_.CountA(q_ylo);
        const u64 wb = active_s_.CountB(q_ylo, q_yhi);
        const u64 w = wa + wb;

        // #region agent log
        {
          std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
          if (log.is_open()) {
            log << "{\"id\":\"log_weight_calc\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:775\",\"message\":\"Weight calculation R-side\",\"data\":{\"eid\":" << eid << ",\"wa\":" << wa << ",\"wb\":" << wb << ",\"w\":" << w << ",\"W_before\":" << W << "},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"E\"}\n";
          }
        }
        // #endregion

        w_a_[eid] = wa;
        w_b_[eid] = wb;
        w_total_[eid] = w;
        
        // #region agent log
        {
          std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
          if (log.is_open()) {
            const u64 W_before = W;
            const bool would_overflow = (W > std::numeric_limits<u64>::max() - w);
            log << "{\"id\":\"log_overflow_check\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:780\",\"message\":\"W accumulation overflow check\",\"data\":{\"W_before\":" << W_before << ",\"w\":" << w << ",\"would_overflow\":" << (would_overflow ? "true" : "false") << ",\"max_u64\":" << std::numeric_limits<u64>::max() << "},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"A\"}\n";
          }
        }
        // #endregion
        
        W += w;

        // Insert q into its side.
        active_r_.Insert(static_cast<u32>(e.index), q_ylo, q_yhi);
      } else {
        const u32 q_ylo = ylo_rank_s_[e.index];
        const u32 q_yhi = yhi_lb_rank_s_[e.index];

        const u64 wa = active_r_.CountA(q_ylo);
        const u64 wb = active_r_.CountB(q_ylo, q_yhi);
        const u64 w = wa + wb;

        // #region agent log
        {
          std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
          if (log.is_open()) {
            log << "{\"id\":\"log_weight_calc\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:790\",\"message\":\"Weight calculation S-side\",\"data\":{\"eid\":" << eid << ",\"wa\":" << wa << ",\"wb\":" << wb << ",\"w\":" << w << ",\"W_before\":" << W << "},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"E\"}\n";
          }
        }
        // #endregion

        w_a_[eid] = wa;
        w_b_[eid] = wb;
        w_total_[eid] = w;
        
        // #region agent log
        {
          std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
          if (log.is_open()) {
            const u64 W_before = W;
            const bool would_overflow = (W > std::numeric_limits<u64>::max() - w);
            log << "{\"id\":\"log_overflow_check\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:795\",\"message\":\"W accumulation overflow check\",\"data\":{\"W_before\":" << W_before << ",\"w\":" << w << ",\"would_overflow\":" << (would_overflow ? "true" : "false") << ",\"max_u64\":" << std::numeric_limits<u64>::max() << "},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"A\"}\n";
          }
        }
        // #endregion
        
        W += w;

        active_s_.Insert(static_cast<u32>(e.index), q_ylo, q_yhi);
      }
    }

    // #region agent log
    {
      std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
      if (log.is_open()) {
        log << "{\"id\":\"log_count_end\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:801\",\"message\":\"Count phase1 end\",\"data\":{\"W\":" << W << "},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"A\"}\n";
      }
    }
    // #endregion

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
    if (!built_ || ds_ == nullptr) {
      if (err) *err = "OursSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "OursSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "OursSamplingBaseline::Sample: out is null";
      return false;
    }

    const u32 t = cfg.run.t;

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    if (t == 0) return true;

    // Ensure weights are available.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    if (W_ == 0) {
      // Empty join: nothing to sample.
      return true;
    }

    // Phase 2: sample (event, pattern) assignments.
    struct Assignment {
      u32 eid;
      u8 pat;   // 0=A, 1=B
      u32 slot; // output slot
    };

    std::vector<Assignment> assign;
    {
      auto scoped = phases ? phases->Scoped("phase2_assign") : PhaseRecorder::ScopedPhase(nullptr, "");

      sampling::AliasTable alias;
      if (!alias.BuildFromU64(Span<const u64>(w_total_), err)) {
        if (err && err->empty()) *err = "OursSamplingBaseline::Sample: failed to build alias table";
        return false;
      }

      assign.reserve(t);
      for (u32 j = 0; j < t; ++j) {
        // Robustly avoid any accidental zero-weight bucket (should not happen).
        u32 eid = 0;
        u64 w = 0;
        int tries_count = 0;
        for (int tries = 0; tries < 16; ++tries) {
          eid = static_cast<u32>(alias.Sample(rng));
          w = w_total_[eid];
          tries_count = tries + 1;
          if (w > 0) break;
        }
        
        // #region agent log
        {
          std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
          if (log.is_open()) {
            log << "{\"id\":\"log_alias_sample\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:871\",\"message\":\"Alias table sampling\",\"data\":{\"slot\":" << j << ",\"eid\":" << eid << ",\"w\":" << w << ",\"tries\":" << tries_count << ",\"W_\":" << W_ << "},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"B\"}\n";
          }
        }
        // #endregion
        
        if (w == 0) {
          if (err) *err = "OursSamplingBaseline::Sample: alias produced only zero-weight events (unexpected)";
          return false;
        }

        const u64 wa = w_a_[eid];
        const u64 wb = w_b_[eid];
        SJS_DASSERT(wa + wb == w);

        u8 pat = 0;
        if (wa == 0) {
          pat = 1;
        } else if (wb == 0) {
          pat = 0;
        } else {
          const u64 r = rng->UniformU64(w);
          pat = (r < wa) ? 0 : 1;
        }

        // #region agent log
        {
          std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
          if (log.is_open()) {
            log << "{\"id\":\"log_pattern_choice\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:888\",\"message\":\"Pattern A/B choice\",\"data\":{\"slot\":" << j << ",\"eid\":" << eid << ",\"wa\":" << wa << ",\"wb\":" << wb << ",\"pat\":" << static_cast<int>(pat) << "},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"E\"}\n";
          }
        }
        // #endregion

        assign.push_back(Assignment{eid, pat, j});
      }

      std::sort(assign.begin(), assign.end(), [](const Assignment& a, const Assignment& b) {
        if (a.eid != b.eid) return a.eid < b.eid;
        if (a.pat != b.pat) return a.pat < b.pat;
        return a.slot < b.slot;
      });
    }

    // Prepare output.
    out->pairs.resize(t);

    // Phase 3: second sweep, fulfill assignments per event/pattern.
    {
      auto scoped = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

      usize ptr = 0;
      std::vector<u32> tmp_handles;

      for (usize ev_pos = 0; ev_pos < events_.size(); ++ev_pos) {
        const auto& e = events_[ev_pos];

        if (e.kind == join::EventKind::End) {
          if (e.side == join::Side::R) {
            const u32 h = static_cast<u32>(e.index);
            active_r_.Erase(h, ylo_rank_r_[e.index]);
          } else {
            const u32 h = static_cast<u32>(e.index);
            active_s_.Erase(h, ylo_rank_s_[e.index]);
          }
          continue;
        }

        // START.
        const i32 eid_i32 = start_id_of_event_[ev_pos];
        SJS_DASSERT(eid_i32 >= 0);
        const u32 eid = static_cast<u32>(eid_i32);

        // #region agent log
        {
          std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
          if (log.is_open()) {
            log << "{\"id\":\"log_phase3_event\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:927\",\"message\":\"Phase3 event processing\",\"data\":{\"ev_pos\":" << ev_pos << ",\"eid\":" << eid << ",\"side\":\"" << (e.side == join::Side::R ? "R" : "S") << "\",\"ptr\":" << ptr << ",\"assign_size\":" << assign.size() << "},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"C\"}\n";
          }
        }
        // #endregion

        // Consume all assignments for this event.
        while (ptr < assign.size() && assign[ptr].eid == eid) {
          const u8 pat = assign[ptr].pat;
          const usize begin = ptr;
          while (ptr < assign.size() && assign[ptr].eid == eid && assign[ptr].pat == pat) {
            ++ptr;
          }
          const u32 k = static_cast<u32>(ptr - begin);
          if (k == 0) continue;

          bool ok = false;

          if (e.side == join::Side::R) {
            const u32 q_ylo = ylo_rank_r_[e.index];
            const u32 q_yhi = yhi_lb_rank_r_[e.index];

            if (pat == 0) {
              ok = active_s_.SampleA(q_ylo, k, rng, &tmp_handles);
            } else {
              ok = active_s_.SampleB(q_ylo, q_yhi, k, rng, &tmp_handles);
            }
            
            // #region agent log
            {
              std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
              if (log.is_open()) {
                log << "{\"id\":\"log_sample_result\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:950\",\"message\":\"SampleA/B result R-side\",\"data\":{\"eid\":" << eid << ",\"pat\":" << static_cast<int>(pat) << ",\"k\":" << k << ",\"ok\":" << (ok ? "true" : "false") << "},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"D\"}\n";
              }
            }
            // #endregion
            
            if (!ok) {
              if (err) *err = "OursSamplingBaseline::Sample: unexpected empty candidate set in phase3 (R-start)";
              return false;
            }

            for (u32 i = 0; i < k; ++i) {
              const u32 slot = assign[begin + i].slot;
              const u32 other_h = tmp_handles[i];
              out->pairs[slot] = PairId{ds_->R.GetId(e.index), ds_->S.GetId(other_h)};
            }

          } else {
            const u32 q_ylo = ylo_rank_s_[e.index];
            const u32 q_yhi = yhi_lb_rank_s_[e.index];

            if (pat == 0) {
              ok = active_r_.SampleA(q_ylo, k, rng, &tmp_handles);
            } else {
              ok = active_r_.SampleB(q_ylo, q_yhi, k, rng, &tmp_handles);
            }
            
            // #region agent log
            {
              std::ofstream log("/home/dhy/PhD/Join-Sampling/.cursor/debug.log", std::ios::app);
              if (log.is_open()) {
                log << "{\"id\":\"log_sample_result\",\"timestamp\":" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << ",\"location\":\"sampling.h:970\",\"message\":\"SampleA/B result S-side\",\"data\":{\"eid\":" << eid << ",\"pat\":" << static_cast<int>(pat) << ",\"k\":" << k << ",\"ok\":" << (ok ? "true" : "false") << "},\"sessionId\":\"debug-session\",\"runId\":\"run1\",\"hypothesisId\":\"D\"}\n";
              }
            }
            // #endregion
            
            if (!ok) {
              if (err) *err = "OursSamplingBaseline::Sample: unexpected empty candidate set in phase3 (S-start)";
              return false;
            }

            for (u32 i = 0; i < k; ++i) {
              const u32 slot = assign[begin + i].slot;
              const u32 other_h = tmp_handles[i];
              out->pairs[slot] = PairId{ds_->R.GetId(other_h), ds_->S.GetId(e.index)};
            }
          }
        }

        // Insert q into its side.
        if (e.side == join::Side::R) {
          const u32 q_ylo = ylo_rank_r_[e.index];
          const u32 q_yhi = yhi_lb_rank_r_[e.index];
          active_r_.Insert(static_cast<u32>(e.index), q_ylo, q_yhi);
        } else {
          const u32 q_ylo = ylo_rank_s_[e.index];
          const u32 q_yhi = yhi_lb_rank_s_[e.index];
          active_s_.Insert(static_cast<u32>(e.index), q_ylo, q_yhi);
        }
      }

      if (ptr != assign.size()) {
        if (err) *err = "OursSamplingBaseline::Sample: internal error (not all assignments consumed)";
        return false;
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || ds_ == nullptr) {
      if (err) *err = "OursSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    join::PlaneSweepOptions opt;
    opt.axis = 0;
    opt.side_order = join::SideTieBreak::RBeforeS;
    opt.skip_axis_check = true;

    return std::make_unique<detail::PlaneSweepEnumeratorWrapper<Dim, T>>(ds_->R, ds_->S, opt);
  }

 private:
  const DatasetT* ds_{nullptr};
  bool built_{false};

  std::vector<join::Event> events_;
  std::vector<i32> start_id_of_event_;  // -1 for END events
  std::vector<usize> start_event_pos_;  // start-id -> event position

  // y-domain (unique y-lower values); ranks are indices into this vector.
  std::vector<T> y_coords_;

  // Precomputed ranks per rectangle for fast sweeps.
  std::vector<u32> ylo_rank_r_;
  std::vector<u32> yhi_lb_rank_r_;
  std::vector<u32> ylo_rank_s_;
  std::vector<u32> yhi_lb_rank_s_;

  // Active indices per side.
  detail::ActiveIndex2D active_r_;
  detail::ActiveIndex2D active_s_;

  // Phase-1 weights per START event.
  std::vector<u64> w_total_;
  std::vector<u64> w_a_;
  std::vector<u64> w_b_;

  u64 W_{0};
  bool weights_valid_{false};
};

}  // namespace ours
}  // namespace baselines
}  // namespace sjs
