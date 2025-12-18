#pragma once
// sjs/baselines/tsunami/sampling.h
//
// Tsunami baseline (Variant::Sampling).
//
// This baseline follows the idea in Tsunami’20-based writeup: reduce box
// intersection join to a 2*Dim dimensional orthogonal range filtering problem
// by encoding each box r as a point p(r) = (L_1..L_D, R_1..R_D), and for each
// query box q run a range query equivalent to:
//   for all i: L_i(r) < R_i(q)  AND  R_i(r) > L_i(q)
//
// Sampling algorithm (exact, i.i.d., uniform WITH replacement):
//   Pass 1: compute deg[q] for each query q and W = sum deg[q] = |J|.
//   Draw t global ranks U_j ~ Unif{1..W} and sort.
//   Pass 2: re-scan only those queries whose blocks contain at least one rank;
//           within such a query, stop early once the maximum requested local
//           rank is reached.
//
// INDEX IMPLEMENTATION NOTE:
//   The original Tsunami'20 paper uses Grid Tree + Augmented Grid as the learned
//   multi-dimensional index. For this baseline, we use a StaticKDTree as a
//   practical replacement that provides the same interface (range filtering over
//   high-dimensional points) and deterministic query results. This substitution
//   maintains the algorithmic correctness while making the implementation more
//   accessible and reproducible.
//
// Implementation note (strict inequalities / half-open semantics):
//   The Tsunami'20 paper suggests using integer coordinates with ±1 boundary
//   adjustments (scaling floats to integers if needed). However, we use a
//   floating-point approach with nextafter() because:
//   1. It is exact for floating-point data without scaling errors
//   2. It avoids the need for coordinate scaling/compression preprocessing
//   3. It preserves the original data precision
//   The condition R_i(r) > L_i(q) is encoded as:
//     R_i(r) >= nextafter(L_i(q), +inf)
//   so that a half-open range query [lo, +inf) matches the strict inequality.
//   This is exact for floating-point Scalars.

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/point.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace tsunami {

namespace detail {

// --------------------------
// Helpers: config parsing
// --------------------------

inline bool EqualsIgnoreCase(std::string_view a, std::string_view b) {
  if (a.size() != b.size()) return false;
  for (usize i = 0; i < a.size(); ++i) {
    const char ca = a[i];
    const char cb = b[i];
    const char la = (ca >= 'A' && ca <= 'Z') ? static_cast<char>(ca - 'A' + 'a') : ca;
    const char lb = (cb >= 'A' && cb <= 'Z') ? static_cast<char>(cb - 'A' + 'a') : cb;
    if (la != lb) return false;
  }
  return true;
}

inline bool ParseSide(std::string_view s, join::Side* out) {
  if (!out) return false;
  if (EqualsIgnoreCase(s, "r") || EqualsIgnoreCase(s, "R")) {
    *out = join::Side::R;
    return true;
  }
  if (EqualsIgnoreCase(s, "s") || EqualsIgnoreCase(s, "S")) {
    *out = join::Side::S;
    return true;
  }
  return false;
}

inline join::Side ChooseIndexSideDefault(const Dataset<2>& ds) {
  // Default: index the larger relation.
  return (ds.S.Size() >= ds.R.Size()) ? join::Side::S : join::Side::R;
}

template <int Dim, class T>
join::Side ChooseIndexSide(const Dataset<Dim, T>& ds, const Config& cfg) {
  // Config override: run.extra["tsunami.index_side"] in {"R","S","larger","smaller"}.
  auto it = cfg.run.extra.find("tsunami.index_side");
  if (it != cfg.run.extra.end()) {
    const std::string_view v = it->second;
    join::Side s{};
    if (ParseSide(v, &s)) return s;
    if (EqualsIgnoreCase(v, "larger")) {
      return (ds.S.Size() >= ds.R.Size()) ? join::Side::S : join::Side::R;
    }
    if (EqualsIgnoreCase(v, "smaller")) {
      return (ds.S.Size() <= ds.R.Size()) ? join::Side::S : join::Side::R;
    }
  }

  // Default heuristic.
  return (ds.S.Size() >= ds.R.Size()) ? join::Side::S : join::Side::R;
}

// --------------------------
// Static KD-tree with streaming range iterator.
// --------------------------
// We keep a small KD-tree implementation local to this baseline because we
// need an incremental iterator (Next()) to implement deterministic join
// enumeration without materializing full query results.
//
// NOTE: This serves as a practical replacement for Tsunami's Grid Tree +
// Augmented Grid index, providing the same range filtering interface while
// being simpler to implement and verify. The deterministic traversal order
// ensures consistent results across runs.

template <int K, class T>
class StaticKDTree {
 public:
  static_assert(std::is_floating_point_v<T>, "StaticKDTree requires floating-point T (for +/-inf and nextafter)");
  using PointT = Point<K, T>;
  using BoxT = Box<K, T>;

  struct BuildOptions {
    u32 leaf_size = 16;
  };

  StaticKDTree() = default;

  void Clear() {
    points_.clear();
    nodes_.clear();
    indices_.clear();
    root_ = kNull;
  }

  usize Size() const noexcept { return points_.size(); }
  bool Empty() const noexcept { return points_.empty(); }

  void Build(const std::vector<PointT>& points, BuildOptions opt = BuildOptions{}) {
    Clear();
    points_ = points;
    opt_ = opt;
    if (points_.empty()) return;

    indices_.resize(points_.size());
    for (usize i = 0; i < points_.size(); ++i) indices_[i] = static_cast<u32>(i);

    nodes_.reserve(points_.size() * 2);
    root_ = BuildRec(/*begin=*/0, /*end=*/static_cast<u32>(points_.size()));
  }

  struct RangeIter {
    const StaticKDTree* tree = nullptr;
    BoxT range = BoxT::Empty();
    std::vector<u32> stack;

    // Leaf scan state.
    bool scanning = false;
    u32 leaf_pos = 0;
    u32 leaf_end = 0;

    // Optional counter for diagnostics.
    u64* candidate_checks = nullptr;

    RangeIter() = default;

    RangeIter(const StaticKDTree* t, const BoxT& r, u64* cand = nullptr) { Reset(t, r, cand); }

    void Reset(const StaticKDTree* t, const BoxT& r, u64* cand = nullptr) {
      tree = t;
      range = r;
      candidate_checks = cand;
      stack.clear();
      scanning = false;
      leaf_pos = leaf_end = 0;

      if (!tree || tree->root_ == kNull || range.IsEmpty()) return;
      stack.reserve(64);
      stack.push_back(tree->root_);
    }

    bool Next(usize* out_point_index) {
      if (!tree || tree->root_ == kNull) return false;
      if (!out_point_index) return false;

      while (true) {
        // Continue scanning current leaf.
        if (scanning) {
          while (leaf_pos < leaf_end) {
            const u32 pidx = tree->indices_[leaf_pos++];
            if (candidate_checks) ++(*candidate_checks);
            if (range.Contains(tree->points_[pidx])) {
              *out_point_index = static_cast<usize>(pidx);
              return true;
            }
          }
          scanning = false;
          continue;
        }

        // Pop next node.
        if (stack.empty()) return false;
        const u32 ni = stack.back();
        stack.pop_back();
        const Node& node = tree->nodes_[ni];
        if (!node.bbox.Intersects(range)) {
          continue;
        }
        if (node.IsLeaf()) {
          scanning = true;
          leaf_pos = node.begin;
          leaf_end = node.end;
          continue;
        }
        // Deterministic traversal: push left then right so right subtree is processed first.
        stack.push_back(node.left);
        stack.push_back(node.right);
      }
    }
  };

 private:
  static constexpr u32 kNull = std::numeric_limits<u32>::max();

  struct Node {
    BoxT bbox = BoxT::Empty();
    u32 left = kNull;
    u32 right = kNull;
    u32 begin = 0;
    u32 end = 0;
    u8 axis = 0;

    bool IsLeaf() const noexcept { return left == kNull && right == kNull; }
  };

  std::vector<PointT> points_;
  std::vector<Node> nodes_;
  std::vector<u32> indices_;
  u32 root_ = kNull;
  BuildOptions opt_{};

  BoxT ComputeBounds(u32 begin, u32 end) const {
    BoxT b = BoxT::Empty();
    for (u32 i = begin; i < end; ++i) {
      b.ExpandToIncludePoint(points_[indices_[i]]);
    }
    return b;
  }

  int ChooseSplitAxis(const BoxT& bbox) const noexcept {
    int axis = 0;
    T best = bbox.Width(0);
    for (int d = 1; d < K; ++d) {
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

    const u32 me = static_cast<u32>(nodes_.size());
    nodes_.push_back(node);

    if (count <= opt_.leaf_size) {
      return me;
    }

    const int axis = ChooseSplitAxis(node.bbox);
    nodes_[me].axis = static_cast<u8>(axis);

    const u32 mid = begin + count / 2;
    std::nth_element(indices_.begin() + begin,
                     indices_.begin() + mid,
                     indices_.begin() + end,
                     [&](u32 ia, u32 ib) {
                       const T a = points_[ia].v[axis];
                       const T b = points_[ib].v[axis];
                       if (a < b) return true;
                       if (b < a) return false;
                       return ia < ib;
                     });

    // Handle degenerate split: if all values equal along axis, split by index.
    bool all_equal = true;
    {
      const T pivot = points_[indices_[mid]].v[axis];
      for (u32 i = begin; i < end; ++i) {
        if (points_[indices_[i]].v[axis] != pivot) {
          all_equal = false;
          break;
        }
      }
    }

    u32 left_begin = begin;
    u32 left_end = mid;
    u32 right_begin = mid;
    u32 right_end = end;
    if (all_equal) {
      // Just ensure both sides are non-empty.
      left_end = begin + count / 2;
      right_begin = left_end;
    } else {
      // Tighten so that left and right are non-empty.
      if (left_end == left_begin) {
        left_end = left_begin + 1;
        right_begin = left_end;
      }
      if (right_begin >= right_end) {
        right_begin = right_end - 1;
        left_end = right_begin;
      }
    }

    nodes_[me].left = BuildRec(left_begin, left_end);
    nodes_[me].right = BuildRec(right_begin, right_end);
    return me;
  }
};

// --------------------------
// Box -> point embedding for Tsunami baseline.
// --------------------------

template <int Dim, class T>
struct TsunamiEmbedding {
  static_assert(Dim >= 1, "TsunamiEmbedding: Dim must be >= 1");
  static constexpr int K = 2 * Dim;
  static_assert(K <= kMaxSupportedDim,
                "TsunamiEmbedding uses Point<2*Dim>; increase kMaxSupportedDim or reduce Dim");
  static_assert(std::is_floating_point_v<T>, "TsunamiEmbedding requires floating-point T (Scalar)");

  using BoxT = Box<Dim, T>;
  using PointK = Point<K, T>;
  using RangeBox = Box<K, T>;

  static PointK Encode(const BoxT& b) {
    PointK p;
    for (int i = 0; i < Dim; ++i) {
      p.v[static_cast<usize>(i)] = b.lo.v[static_cast<usize>(i)];
      p.v[static_cast<usize>(Dim + i)] = b.hi.v[static_cast<usize>(i)];
    }
    return p;
  }

  // Build the 2*Dim query range for a query box q.
  // Returns false if the implied range is empty.
  static bool MakeQueryRange(const BoxT& q, RangeBox* out) {
    if (!out) return false;
    if (q.IsEmpty()) {
      *out = RangeBox::Empty();
      return false;
    }

    const T pinf = std::numeric_limits<T>::infinity();
    const T ninf = -std::numeric_limits<T>::infinity();

    PointK lo;
    PointK hi;
    for (int i = 0; i < Dim; ++i) {
      // L_i(r) < R_i(q)  ->  L_i(r) in (-inf, R_i(q))
      lo.v[static_cast<usize>(i)] = ninf;
      hi.v[static_cast<usize>(i)] = q.hi.v[static_cast<usize>(i)];

      // R_i(r) > L_i(q)  ->  R_i(r) in (L_i(q), +inf)
      // Represent strict lower bound using nextafter.
      const T open_lo = std::nextafter(q.lo.v[static_cast<usize>(i)], pinf);
      lo.v[static_cast<usize>(Dim + i)] = open_lo;
      hi.v[static_cast<usize>(Dim + i)] = pinf;
    }

    RangeBox range(lo, hi);
    if (range.IsEmpty()) {
      *out = RangeBox::Empty();
      return false;
    }
    *out = range;
    return true;
  }
};

// --------------------------
// Preprocessing state shared by Tsunami variants.
// --------------------------

template <int Dim, class T>
class TsunamiPreproc {
 public:
  using DatasetT = Dataset<Dim, T>;
  using RelT = Relation<Dim, T>;
  using BoxT = Box<Dim, T>;
  static constexpr int K = 2 * Dim;
  using Embed = TsunamiEmbedding<Dim, T>;
  using PointK = typename Embed::PointK;
  using RangeBox = typename Embed::RangeBox;
  using KD = StaticKDTree<K, T>;

  void Reset() {
    ds_ = nullptr;
    built_ = false;
    query_side_ = join::Side::R;
    index_side_ = join::Side::S;
    rel_query_ = nullptr;
    rel_data_ = nullptr;
    query_order_.clear();
    points_.clear();
    kd_.Clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) {
    (void)err;
    auto scoped = phases ? phases->Scoped("tsunami_build") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    ds_ = &ds;

    index_side_ = ChooseIndexSide(ds, cfg);
    query_side_ = join::OtherSide(index_side_);

    rel_query_ = (query_side_ == join::Side::R) ? &ds.R : &ds.S;
    rel_data_ = (index_side_ == join::Side::R) ? &ds.R : &ds.S;

    // Deterministic query order: by stable object id ascending.
    {
      auto _ = phases ? phases->Scoped("tsunami_build_query_order") : PhaseRecorder::ScopedPhase(nullptr, "");
      query_order_.resize(rel_query_->Size());
      for (usize i = 0; i < rel_query_->Size(); ++i) query_order_[i] = static_cast<u32>(i);
      std::sort(query_order_.begin(), query_order_.end(), [&](u32 a, u32 b) {
        const Id ida = rel_query_->GetId(static_cast<usize>(a));
        const Id idb = rel_query_->GetId(static_cast<usize>(b));
        if (ida < idb) return true;
        if (idb < ida) return false;
        return a < b;
      });
    }

    // Encode data side boxes to points and build KD-tree.
    {
      auto _ = phases ? phases->Scoped("tsunami_build_index") : PhaseRecorder::ScopedPhase(nullptr, "");
      points_.reserve(rel_data_->Size());
      for (const auto& b : rel_data_->boxes) {
        points_.push_back(Embed::Encode(b));
      }
      kd_.Build(points_, typename KD::BuildOptions{});
    }

    built_ = true;
    return true;
  }

  bool built() const noexcept { return built_; }
  const DatasetT* ds() const noexcept { return ds_; }

  join::Side query_side() const noexcept { return query_side_; }
  join::Side index_side() const noexcept { return index_side_; }
  const RelT& QueryRel() const {
    SJS_DASSERT(rel_query_);
    return *rel_query_;
  }
  const RelT& DataRel() const {
    SJS_DASSERT(rel_data_);
    return *rel_data_;
  }
  const std::vector<u32>& QueryOrder() const noexcept { return query_order_; }
  const KD& Index() const noexcept { return kd_; }

  bool MakeRangeForQuery(const BoxT& q, RangeBox* out) const { return Embed::MakeQueryRange(q, out); }

  // Convert (query index into QueryRel, data index into DataRel) into PairId in (R,S) order.
  PairId MakePair(u32 q_idx, u32 data_idx) const {
    SJS_DASSERT(ds_);
    if (query_side_ == join::Side::R) {
      // Query is R, data is S.
      return PairId{ds_->R.GetId(static_cast<usize>(q_idx)), ds_->S.GetId(static_cast<usize>(data_idx))};
    }
    // Query is S, data is R.
    return PairId{ds_->R.GetId(static_cast<usize>(data_idx)), ds_->S.GetId(static_cast<usize>(q_idx))};
  }

 private:
  const DatasetT* ds_ = nullptr;
  bool built_ = false;

  join::Side query_side_{join::Side::R};
  join::Side index_side_{join::Side::S};
  const RelT* rel_query_ = nullptr;
  const RelT* rel_data_ = nullptr;

  std::vector<u32> query_order_;

  std::vector<PointK> points_;
  KD kd_;
};

// --------------------------
// Deterministic join enumerator using TsunamiPreproc.
// --------------------------

template <int Dim, class T>
class TsunamiJoinEnumerator final : public IJoinEnumerator {
 public:
  using BoxT = Box<Dim, T>;
  using Preproc = TsunamiPreproc<Dim, T>;
  static constexpr int K = 2 * Dim;
  using Embed = TsunamiEmbedding<Dim, T>;
  using RangeBox = typename Embed::RangeBox;
  using KD = typename Preproc::KD;

  TsunamiJoinEnumerator(const Preproc* prep)
      : prep_(prep) {
    Reset();
  }

  void Reset() override {
    stats_.Reset();
    q_pos_ = 0;
    cur_q_idx_ = 0;
    has_iter_ = false;
    iter_ = typename KD::RangeIter{};
    cand_ = 0;
  }

  bool Next(PairId* out) override {
    if (!out) return false;
    if (!prep_ || !prep_->built()) return false;

    const auto& rel_q = prep_->QueryRel();
    const auto& order = prep_->QueryOrder();

    while (true) {
      if (has_iter_) {
        usize data_pt = 0;
        if (iter_.Next(&data_pt)) {
          // data_pt is an index into DataRel().boxes
          *out = prep_->MakePair(cur_q_idx_, static_cast<u32>(data_pt));
          ++stats_.output_pairs;
          stats_.candidate_checks = cand_;
          return true;
        }
        has_iter_ = false;
        continue;
      }

      if (q_pos_ >= order.size()) {
        stats_.candidate_checks = cand_;
        return false;
      }

      cur_q_idx_ = order[q_pos_++];
      ++stats_.num_events;  // treat each query as an "event" for bookkeeping

      RangeBox range;
      const BoxT& q = rel_q.boxes[static_cast<usize>(cur_q_idx_)];
      if (!prep_->MakeRangeForQuery(q, &range)) {
        // Empty range => no matches.
        continue;
      }

      iter_.Reset(&prep_->Index(), range, &cand_);
      has_iter_ = true;
    }
  }

  const join::JoinStats& Stats() const noexcept override {
    // Keep candidate_checks in sync with the iterator's running counter.
    stats_.candidate_checks = cand_;
    return stats_;
  }

 private:
  const Preproc* prep_ = nullptr;

  // NOTE: Stats() is const but the interface returns a const reference; we
  // keep candidate_checks in sync by marking stats_ mutable.
  mutable join::JoinStats stats_;
  usize q_pos_ = 0;
  u32 cur_q_idx_ = 0;

  bool has_iter_ = false;
  typename KD::RangeIter iter_;
  u64 cand_ = 0;
};

// --------------------------
// Two-pass rank sampling using precomputed deg[] (Sampling/Adaptive fallback).
// --------------------------

struct RankRec {
  u64 U = 0;     // 1-based global rank in [1..W]
  u64 slot = 0;  // output position
};

inline bool RankRecLess(const RankRec& a, const RankRec& b) noexcept {
  if (a.U < b.U) return true;
  if (b.U < a.U) return false;
  return a.slot < b.slot;
}

struct LocalRec {
  u64 u = 0;     // 1-based local rank in [1..deg]
  u64 slot = 0;  // output position
};

inline bool LocalRecLess(const LocalRec& a, const LocalRec& b) noexcept {
  if (a.u < b.u) return true;
  if (b.u < a.u) return false;
  return a.slot < b.slot;
}

template <int Dim, class T>
bool TwoPassRankSampleUsingDeg(const TsunamiPreproc<Dim, T>& prep,
                               const std::vector<u64>& deg_order,
                               u64 W,
                               u64 t,
                               Rng* rng,
                               std::vector<PairId>* out_pairs,
                               PhaseRecorder* phases,
                               std::string* err) {
  if (!rng || !out_pairs) {
    if (err) *err = "TwoPassRankSampleUsingDeg: rng/out_pairs is null";
    return false;
  }
  out_pairs->clear();
  if (t == 0 || W == 0) return true;

  if (deg_order.size() != prep.QueryOrder().size()) {
    if (err) *err = "TwoPassRankSampleUsingDeg: deg_order size mismatch";
    return false;
  }

  auto scoped = phases ? phases->Scoped("tsunami_pass2_rank_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

  // Draw ranks.
  std::vector<RankRec> ranks;
  ranks.resize(static_cast<usize>(t));
  for (u64 i = 0; i < t; ++i) {
    const u64 u = rng->UniformU64(W) + 1;  // [1..W]
    ranks[static_cast<usize>(i)] = RankRec{u, i};
  }
  std::sort(ranks.begin(), ranks.end(), RankRecLess);

  out_pairs->assign(static_cast<usize>(t), PairId{});

  const auto& rel_q = prep.QueryRel();
  const auto& order = prep.QueryOrder();

  u64 g = 0;      // processed pairs so far
  usize p = 0;    // pointer into ranks

  std::vector<LocalRec> local;
  local.reserve(64);

  for (usize qi = 0; qi < order.size(); ++qi) {
    if (p >= ranks.size()) break;

    const u64 deg = deg_order[qi];
    if (deg == 0) {
      // No pairs in this block.
      continue;
    }

    // Skip this query if next rank is beyond this block.
    if (ranks[p].U > g + deg) {
      g += deg;
      continue;
    }

    // Collect all local ranks within this block.
    local.clear();
    while (p < ranks.size() && ranks[p].U <= g + deg) {
      const u64 u_local = ranks[p].U - g;  // 1..deg
      local.push_back(LocalRec{u_local, ranks[p].slot});
      ++p;
    }
    std::sort(local.begin(), local.end(), LocalRecLess);
    const u64 umax = local.back().u;

    // Execute query and stop early at umax.
    const u32 q_idx = order[qi];
    typename TsunamiEmbedding<Dim, T>::RangeBox range;
    const auto& qbox = rel_q.boxes[static_cast<usize>(q_idx)];
    if (!prep.MakeRangeForQuery(qbox, &range)) {
      if (err) *err = "TwoPassRankSampleUsingDeg: empty range for a query that should have matches";
      return false;
    }

    typename TsunamiPreproc<Dim, T>::KD::RangeIter it(&prep.Index(), range);
    u64 c = 0;
    usize j = 0;
    usize data_pt = 0;
    while (it.Next(&data_pt)) {
      ++c;
      while (j < local.size() && local[j].u == c) {
        (*out_pairs)[static_cast<usize>(local[j].slot)] = prep.MakePair(q_idx, static_cast<u32>(data_pt));
        ++j;
      }
      if (c == umax) break;
    }

    if (c != umax) {
      if (err) *err = "TwoPassRankSampleUsingDeg: second pass ended early (non-deterministic stream or wrong deg?)";
      return false;
    }

    g += deg;
  }

  if (p != ranks.size()) {
    if (err) *err = "TwoPassRankSampleUsingDeg: not all ranks were assigned (deg/W mismatch?)";
    return false;
  }

  return true;
}

}  // namespace detail

// --------------------------
// TsunamiSamplingBaseline
// --------------------------

template <int Dim, class T = Scalar>
class TsunamiSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 1, "TsunamiSamplingBaseline requires Dim >= 1");
  static_assert(2 * Dim <= kMaxSupportedDim,
                "TsunamiSamplingBaseline uses Point<2*Dim>; increase kMaxSupportedDim or reduce Dim");

  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;

  // NOTE: Method::Tsunami is not in core/types.h yet. We return Unknown for now.
  Method method() const noexcept override { return Method::Tsunami; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "tsunami_sampling"; }

  void Reset() override {
    prep_.Reset();
    counted_ = false;
    W_ = 0;
    deg_order_.clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    Reset();
    if (!prep_.Build(ds, cfg, phases, err)) return false;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;
    if (!out) {
      if (err) *err = "TsunamiSamplingBaseline::Count: out is null";
      return false;
    }
    *out = MakeExactCount(0);

    if (!prep_.built()) {
      if (err) *err = "TsunamiSamplingBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    const auto& rel_q = prep_.QueryRel();
    const auto& order = prep_.QueryOrder();

    deg_order_.assign(order.size(), 0ULL);
    u64 W = 0;

    for (usize qi = 0; qi < order.size(); ++qi) {
      const u32 q_idx = order[qi];
      const BoxT& q = rel_q.boxes[static_cast<usize>(q_idx)];

      typename detail::TsunamiEmbedding<Dim, T>::RangeBox range;
      if (!prep_.MakeRangeForQuery(q, &range)) {
        deg_order_[qi] = 0;
        continue;
      }

      typename detail::TsunamiPreproc<Dim, T>::KD::RangeIter it(&prep_.Index(), range);
      u64 deg = 0;
      usize data_pt = 0;
      while (it.Next(&data_pt)) {
        ++deg;
      }

      deg_order_[qi] = deg;

      if (deg > 0) {
        if (W > std::numeric_limits<u64>::max() - deg) {
          if (err) *err = "TsunamiSamplingBaseline::Count: |J| overflowed u64";
          return false;
        }
        W += deg;
      }
    }

    W_ = W;
    counted_ = true;
    *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!prep_.built()) {
      if (err) *err = "TsunamiSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "TsunamiSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "TsunamiSamplingBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    if (!counted_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    if (W_ == 0) {
      // Empty join.
      return true;
    }

    auto scoped = phases ? phases->Scoped("phase2_3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

    std::vector<PairId> pairs;
    if (!detail::TwoPassRankSampleUsingDeg(prep_, deg_order_, W_, t, rng, &pairs, phases, err)) {
      return false;
    }

    out->pairs = std::move(pairs);
    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");

    if (!prep_.built()) {
      if (err) *err = "TsunamiSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    return std::make_unique<detail::TsunamiJoinEnumerator<Dim, T>>(&prep_);
  }

 private:
  detail::TsunamiPreproc<Dim, T> prep_;

  bool counted_ = false;
  u64 W_ = 0;
  std::vector<u64> deg_order_;  // aligned with prep_.QueryOrder()
};

}  // namespace tsunami
}  // namespace baselines
}  // namespace sjs
