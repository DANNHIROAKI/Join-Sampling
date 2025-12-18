#pragma once
// sjs/baselines/tsunami/sampling.h
//
// Tsunami baseline (Variant::Sampling) — strictly aligned with docs/Baseline/Tsunami’20 Baseline.md.
//
// Goal: exact i.i.d. uniform sampling WITH replacement on the spatial intersection join result J.
//
// Baseline idea:
//   1) Reduce box intersection join to a 2*Dim dimensional orthogonal range filtering problem:
//        encode each data-side box r as a point
//          p(r) = (L_1..L_D, R_1..R_D)
//        and for each query-side box q build a 2*D dimensional range Q(q) such that:
//          r intersects q  <=>  p(r) ∈ Q(q)
//   2) Treat Tsunami as a black-box range filtering engine with interface:
//        BuildTsunamiIndex(Points, WorkloadQueries) -> Index
//        Query(Index, Q(q)) -> stream of matched ids
//      In this repository we provide a small deterministic in-memory iterator to stand in for that
//      interface; it can be swapped with a real Tsunami implementation without changing the sampling logic.
//   3) Implement TSUNAMI-2Pass-RankSample (exact, i.i.d., uniform with replacement):
//        Pass 1: compute deg[q] and W=Σdeg[q]=|J|
//        Draw t global ranks U_j ~ Unif{1..W} and sort
//        Pass 2: re-scan only the q-blocks that contain at least one rank; within a q-block stop early at umax.
//
// Important correctness requirements (mirrors the baseline doc):
//   - Deterministic join stream: query order A is fixed (we sort by object id), and for a fixed query q the
//     range engine enumerates matches in a repeatable order.
//   - Strict inequalities for half-open boxes are handled via "rank + secondary key" coordinate compression
//     (Baseline §3.3, Scheme B). This avoids relying on +/-1 or floating epsilons and is provably exact.

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/point.h"

#include <algorithm>
#include <array>
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
  // Default heuristic: index the larger relation.
  return (ds.S.Size() >= ds.R.Size()) ? join::Side::S : join::Side::R;
}

// --------------------------
// Tsunami black-box adapter (streaming range filter).
// --------------------------
// This is an adapter/stub that matches the baseline doc's black-box interface.
// It is intentionally deterministic so that Pass1 and Pass2 see the same join stream order.
//
// If you later integrate a real Tsunami implementation, keep this interface and swap the internals.

template <int K, class T>
class TsunamiIndexStub {
 public:
  static_assert(std::is_floating_point_v<T>,
                "TsunamiIndexStub currently uses Box<K,T> with +/-inf Empty() sentinels; instantiate with floating T");
  using PointT = Point<K, T>;
  using RangeT = Box<K, T>;

  struct BuildOptions {
    u32 leaf_size = 16;
  };

  TsunamiIndexStub() = default;

  void Clear() {
    points_.clear();
    nodes_.clear();
    indices_.clear();
    root_ = kNull;
  }

  usize Size() const noexcept { return points_.size(); }
  bool Empty() const noexcept { return points_.empty(); }

  // BuildTsunamiIndex(Points, WorkloadQueries) -> Index
  // NOTE: The stub ignores WorkloadQueries, but we keep the signature to align with the baseline doc.
  void Build(const std::vector<PointT>& points,
             const std::vector<RangeT>& workload_queries,
             BuildOptions opt = BuildOptions{}) {
    (void)workload_queries;
    Clear();
    points_ = points;
    opt_ = opt;
    if (points_.empty()) return;

    indices_.resize(points_.size());
    for (usize i = 0; i < points_.size(); ++i) indices_[i] = static_cast<u32>(i);

    nodes_.reserve(points_.size() * 2);
    root_ = BuildRec(/*begin=*/0, /*end=*/static_cast<u32>(points_.size()));
  }

  // Query(Index, range) -> stream of matched ids (indices into Points)
  struct QueryIter {
    const TsunamiIndexStub* tree = nullptr;
    RangeT range = RangeT::Empty();
    std::vector<u32> stack;

    // Leaf scan state.
    bool scanning = false;
    u32 leaf_pos = 0;
    u32 leaf_end = 0;

    // Optional counter for diagnostics.
    u64* candidate_checks = nullptr;

    QueryIter() = default;

    QueryIter(const TsunamiIndexStub* t, const RangeT& r, u64* cand = nullptr) { Reset(t, r, cand); }

    void Reset(const TsunamiIndexStub* t, const RangeT& r, u64* cand = nullptr) {
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
    RangeT bbox = RangeT::Empty();
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

  RangeT ComputeBounds(u32 begin, u32 end) const {
    RangeT b = RangeT::Empty();
    for (u32 i = begin; i < end; ++i) {
      b.ExpandToIncludePoint(points_[indices_[i]]);
    }
    return b;
  }

  int ChooseSplitAxis(const RangeT& bbox) const noexcept {
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
// Rank + secondary-key coordinate encoding (Baseline §3.3, Scheme B).
// --------------------------
//
// For each of the K=2*Dim axes, build a strict total order of the data-side attribute values
// using the key (value, record_id). Each data record gets a unique rank per axis, even when
// values repeat. Strict inequalities in the original domain are mapped to half-open rank ranges:
//   X < b  -> rank(X) in [0, lower_bound((b, MIN_ID)))
//   X > a  -> rank(X) in [upper_bound((a, MAX_ID)), N)
// where N is the number of data records.
// This makes half-open box intersection conditions provably exact without epsilons.
//
// NOTE: This transform is only used to make strict inequalities executable and exact.
// The range filtering engine indexes rank-points and answers rank-ranges.

template <class T>
inline bool KeyLess(const std::pair<T, Id>& a, const std::pair<T, Id>& b) noexcept {
  if (a.first < b.first) return true;
  if (b.first < a.first) return false;
  return a.second < b.second;
}

template <int Dim, class T>
class TsunamiRankEncoding {
 public:
  static_assert(Dim >= 1, "TsunamiRankEncoding requires Dim >= 1");
  static constexpr int K = 2 * Dim;
  using PointK = Point<K, T>;
  using RangeBox = Box<K, T>;
  using BoxT = Box<Dim, T>;

  struct Axis {
    // Sorted (value, id) keys for this axis; size == n_data
    std::vector<std::pair<T, Id>> keys;
    // rank_of_data_idx[data_idx] = rank position in keys[]
    std::vector<u32> rank_of_data_idx;

    void Clear() {
      keys.clear();
      rank_of_data_idx.clear();
    }

    usize Size() const noexcept { return keys.size(); }

    // lower_bound((bound, MIN_ID))
    usize LowerBoundByValue(const T& bound) const {
      const Id kMinId = std::numeric_limits<Id>::lowest();
      const std::pair<T, Id> probe{bound, kMinId};
      return static_cast<usize>(std::lower_bound(keys.begin(), keys.end(), probe, KeyLess<T>) - keys.begin());
    }

    // upper_bound((bound, MAX_ID))
    usize UpperBoundByValue(const T& bound) const {
      const Id kMaxId = std::numeric_limits<Id>::max();
      const std::pair<T, Id> probe{bound, kMaxId};
      return static_cast<usize>(std::upper_bound(keys.begin(), keys.end(), probe, KeyLess<T>) - keys.begin());
    }
  };

  void Reset() {
    n_data_ = 0;
    for (auto& a : axes_) a.Clear();
  }

  usize n_data() const noexcept { return n_data_; }
  const std::array<Axis, K>& axes() const noexcept { return axes_; }

  // Build axis ranks from the data-side relation.
  template <class RelT>
  void BuildFromDataRel(const RelT& rel_data) {
    Reset();
    n_data_ = rel_data.Size();
    for (int axis = 0; axis < K; ++axis) {
      BuildAxis(rel_data, axis);
    }
  }

  // Encode a data-side record (by index) into a K-dim rank point.
  PointK EncodeDataPoint(u32 data_idx) const {
    PointK p;
    for (int axis = 0; axis < K; ++axis) {
      const u32 r = axes_[axis].rank_of_data_idx[static_cast<usize>(data_idx)];
      p.v[static_cast<usize>(axis)] = static_cast<T>(r);
    }
    return p;
  }

  // Build the K-dim rank range for a query box q.
  // Returns false if the implied rank-range is empty.
  bool MakeQueryRange(const BoxT& q, RangeBox* out) const {
    if (!out) return false;
    if (q.IsEmpty()) {
      *out = RangeBox::Empty();
      return false;
    }

    // N can be 0.
    const usize N = n_data_;
    PointK lo;
    PointK hi;

    // First Dim axes are L_i(r) values: require L_i(r) < R_i(q).
    // Second Dim axes are R_i(r) values: require R_i(r) > L_i(q).
    for (int i = 0; i < Dim; ++i) {
      const usize axis_L = static_cast<usize>(i);
      const usize axis_R = static_cast<usize>(Dim + i);

      const T bound_hi = q.hi.v[static_cast<usize>(i)];
      const T bound_lo = q.lo.v[static_cast<usize>(i)];

      const usize hi_rank = axes_[axis_L].LowerBoundByValue(bound_hi);    // [0..N]
      const usize lo_rank = axes_[axis_R].UpperBoundByValue(bound_lo);    // [0..N]

      lo.v[axis_L] = static_cast<T>(0);
      hi.v[axis_L] = static_cast<T>(hi_rank);

      lo.v[axis_R] = static_cast<T>(lo_rank);
      hi.v[axis_R] = static_cast<T>(N);
    }

    RangeBox range(lo, hi);
    if (range.IsEmpty()) {
      *out = RangeBox::Empty();
      return false;
    }
    *out = range;
    return true;
  }

 private:
  usize n_data_ = 0;
  std::array<Axis, K> axes_{};

  template <class RelT>
  void BuildAxis(const RelT& rel_data, int axis) {
    Axis& A = axes_[static_cast<usize>(axis)];
    A.Clear();
    A.keys.reserve(n_data_);
    A.rank_of_data_idx.assign(n_data_, 0U);

    // Temporary array with data_idx so we can write rank_of_data_idx.
    struct Tmp {
      T v;
      Id id;
      u32 data_idx;
    };
    std::vector<Tmp> tmp;
    tmp.reserve(n_data_);

    for (usize i = 0; i < n_data_; ++i) {
      const auto& b = rel_data.boxes[i];
      const Id id = rel_data.GetId(i);

      T v{};
      if (axis < Dim) {
        v = b.lo.v[static_cast<usize>(axis)];
      } else {
        v = b.hi.v[static_cast<usize>(axis - Dim)];
      }
      tmp.push_back(Tmp{v, id, static_cast<u32>(i)});
    }

    std::sort(tmp.begin(), tmp.end(), [](const Tmp& a, const Tmp& b) {
      if (a.v < b.v) return true;
      if (b.v < a.v) return false;
      return a.id < b.id;
    });

    A.keys.resize(n_data_);
    for (usize pos = 0; pos < n_data_; ++pos) {
      A.keys[pos] = std::pair<T, Id>{tmp[pos].v, tmp[pos].id};
      A.rank_of_data_idx[static_cast<usize>(tmp[pos].data_idx)] = static_cast<u32>(pos);
    }
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

  using PointK = Point<K, T>;
  using RangeBox = Box<K, T>;
  using IndexT = TsunamiIndexStub<K, T>;
  using RankEnc = TsunamiRankEncoding<Dim, T>;

  void Reset() {
    ds_ = nullptr;
    built_ = false;
    query_side_ = join::Side::R;
    index_side_ = join::Side::S;
    rel_query_ = nullptr;
    rel_data_ = nullptr;
    query_order_.clear();
    rank_points_.clear();
    workload_ranges_.clear();
    index_.Clear();
    rank_enc_.Reset();
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

    // Deterministic query order A: by stable object id ascending.
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

    // Build rank encoding (Scheme B) from the data-side relation.
    {
      auto _ = phases ? phases->Scoped("tsunami_build_rank_encoding") : PhaseRecorder::ScopedPhase(nullptr, "");
      rank_enc_.BuildFromDataRel(*rel_data_);
    }

    // Encode data-side boxes to rank points.
    {
      auto _ = phases ? phases->Scoped("tsunami_build_points") : PhaseRecorder::ScopedPhase(nullptr, "");
      rank_points_.reserve(rel_data_->Size());
      for (usize i = 0; i < rel_data_->Size(); ++i) {
        rank_points_.push_back(rank_enc_.EncodeDataPoint(static_cast<u32>(i)));
      }
    }

    // Build workload query set WorkloadQueries = {Q(q) | q in A}.
    {
      auto _ = phases ? phases->Scoped("tsunami_build_workload_queries") : PhaseRecorder::ScopedPhase(nullptr, "");
      workload_ranges_.resize(query_order_.size());
      const auto& rel_q = *rel_query_;
      for (usize qi = 0; qi < query_order_.size(); ++qi) {
        const u32 q_idx = query_order_[qi];
        const BoxT& q = rel_q.boxes[static_cast<usize>(q_idx)];
        RangeBox range;
        if (!rank_enc_.MakeQueryRange(q, &range)) {
          workload_ranges_[qi] = RangeBox::Empty();
        } else {
          workload_ranges_[qi] = range;
        }
      }
    }

    // Build index: BuildTsunamiIndex(Points, WorkloadQueries).
    {
      auto _ = phases ? phases->Scoped("tsunami_build_index") : PhaseRecorder::ScopedPhase(nullptr, "");
      typename IndexT::BuildOptions opt{};
      // Optional: cfg.run.extra["tsunami.leaf_size"] (stub only).
      auto it = cfg.run.extra.find("tsunami.leaf_size");
      if (it != cfg.run.extra.end()) {
        try {
          const int v = std::stoi(it->second);
          if (v > 0) opt.leaf_size = static_cast<u32>(v);
        } catch (...) {
          // ignore parse failure; keep default
        }
      }
      index_.Build(rank_points_, workload_ranges_, opt);
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
  const std::vector<RangeBox>& WorkloadRanges() const noexcept { return workload_ranges_; }
  const IndexT& Index() const noexcept { return index_; }

  // Range for a query block position qi (aligned with QueryOrder()).
  const RangeBox& RangeAt(usize qi) const {
    SJS_DASSERT(qi < workload_ranges_.size());
    return workload_ranges_[qi];
  }

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

  RankEnc rank_enc_;
  std::vector<PointK> rank_points_;
  std::vector<RangeBox> workload_ranges_;  // aligned with query_order_ (i.e., A[1..n1] order)
  IndexT index_;
};

// --------------------------
// Deterministic join enumerator using TsunamiPreproc.
// --------------------------

template <int Dim, class T>
class TsunamiJoinEnumerator final : public IJoinEnumerator {
 public:
  using Preproc = TsunamiPreproc<Dim, T>;
  using IndexT = typename Preproc::IndexT;
  using RangeBox = typename Preproc::RangeBox;

  explicit TsunamiJoinEnumerator(const Preproc* prep)
      : prep_(prep) {
    Reset();
  }

  void Reset() override {
    stats_.Reset();
    q_pos_ = 0;
    cur_q_idx_ = 0;
    has_iter_ = false;
    iter_ = typename IndexT::QueryIter{};
    cand_ = 0;
  }

  bool Next(PairId* out) override {
    if (!out) return false;
    if (!prep_ || !prep_->built()) return false;

    const auto& order = prep_->QueryOrder();

    while (true) {
      if (has_iter_) {
        usize data_pt = 0;
        if (iter_.Next(&data_pt)) {
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

      const usize qi = q_pos_;
      cur_q_idx_ = order[q_pos_++];
      ++stats_.num_events;  // treat each query as an "event" for bookkeeping

      const RangeBox& range = prep_->RangeAt(qi);
      if (range.IsEmpty()) {
        continue;  // no matches for this query
      }

      iter_.Reset(&prep_->Index(), range, &cand_);
      has_iter_ = true;
    }
  }

  const join::JoinStats& Stats() const noexcept override {
    stats_.candidate_checks = cand_;
    return stats_;
  }

 private:
  const Preproc* prep_ = nullptr;

  mutable join::JoinStats stats_;
  usize q_pos_ = 0;
  u32 cur_q_idx_ = 0;

  bool has_iter_ = false;
  typename IndexT::QueryIter iter_;
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

  // Draw i.i.d. global ranks with replacement.
  std::vector<RankRec> ranks;
  ranks.resize(static_cast<usize>(t));
  for (u64 i = 0; i < t; ++i) {
    const u64 u = rng->UniformU64(W) + 1;  // [1..W]
    ranks[static_cast<usize>(i)] = RankRec{u, i};
  }
  std::sort(ranks.begin(), ranks.end(), RankRecLess);

  out_pairs->assign(static_cast<usize>(t), PairId{});

  const auto& order = prep.QueryOrder();
  const auto& ranges = prep.WorkloadRanges();

  u64 g = 0;      // processed pairs so far
  usize p = 0;    // pointer into ranks

  std::vector<LocalRec> local;
  local.reserve(64);

  for (usize qi = 0; qi < order.size(); ++qi) {
    if (p >= ranks.size()) break;

    const u64 deg = deg_order[qi];
    if (deg == 0) {
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
    const auto& range = ranges[qi];
    if (range.IsEmpty()) {
      if (err) *err = "TwoPassRankSampleUsingDeg: empty range for a query that should have matches";
      return false;
    }

    typename TsunamiPreproc<Dim, T>::IndexT::QueryIter it(&prep.Index(), range);
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
// TsunamiSamplingBaseline  (TSUNAMI-2Pass-RankSample)
// --------------------------

template <int Dim, class T = Scalar>
class TsunamiSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 1, "TsunamiSamplingBaseline requires Dim >= 1");
  static_assert(2 * Dim <= kMaxSupportedDim,
                "TsunamiSamplingBaseline uses Point<2*Dim>; increase kMaxSupportedDim or reduce Dim");

  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;

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
    return prep_.Build(ds, cfg, phases, err);
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

    const auto& order = prep_.QueryOrder();
    const auto& ranges = prep_.WorkloadRanges();

    deg_order_.assign(order.size(), 0ULL);
    u64 W = 0;

    for (usize qi = 0; qi < order.size(); ++qi) {
      const auto& range = ranges[qi];
      if (range.IsEmpty()) {
        deg_order_[qi] = 0;
        continue;
      }

      typename detail::TsunamiPreproc<Dim, T>::IndexT::QueryIter it(&prep_.Index(), range);
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

    auto scoped = phases ? phases->Scoped("phase2_pass2_locate") : PhaseRecorder::ScopedPhase(nullptr, "");

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
