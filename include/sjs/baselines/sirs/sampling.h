#pragma once
// sjs/baselines/sirs/sampling.h
//
// SIRS baseline (Variant::Sampling).
//
// This baseline implements the reduction described in "SIRS’21.md":
//   - Choose one relation as OUTER and the other as INNER.
//   - Embed each INNER rectangle b into a 2*Dim point phi(b).
//   - For each OUTER rectangle a, construct a 2*Dim range Q(a) such that:
//         b intersects a  <=>  phi(b) \in Q(a)
//     under half-open rectangle semantics [lo,hi).
//   - Phase 1 (Count): compute deg(a)=|{b in INNER : b intersects a}| for every
//     a in OUTER using a point range-count index; W = sum_a deg(a) = |J|.
//   - Phase 2/3 (Sample): sample OUTER a with probability deg(a)/W, then sample
//     one b uniformly from INNER\cap Q(a) (using a range sampler). This yields
//     uniform i.i.d. join samples with replacement.
//
// Practical note:
//   The original SIRS paper provides a specialized independent range sampler.
//   For this SIGMOD-grade experimental harness we implement a deterministic
//   KD-tree based range sampler that supports:
//     Count(Q), Enumerate(Q), and Sample(Q,k) uniformly with replacement.
//   This preserves the *algorithmic contract* used by the baseline.

#include "sjs/baselines/baseline_api.h"
#include "sjs/core/assert.h"
#include "sjs/geometry/embedding.h"
#include "sjs/geometry/predicates.h"
#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace sirs {

namespace detail {

// --------------------------
// Small helpers for cfg.run.extra
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

inline bool ParseBool(std::string_view s, bool* out) {
  if (!out) return false;
  if (s.empty()) return false;
  std::string t(s);
  std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (t == "1" || t == "true" || t == "yes" || t == "y" || t == "on") {
    *out = true;
    return true;
  }
  if (t == "0" || t == "false" || t == "no" || t == "n" || t == "off") {
    *out = false;
    return true;
  }
  return false;
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

inline bool ExtraBoolOr(const std::unordered_map<std::string, std::string>& extra,
                        std::string_view key,
                        bool default_v) {
  auto it = extra.find(std::string(key));
  if (it == extra.end()) return default_v;
  bool v = default_v;
  if (!ParseBool(it->second, &v)) return default_v;
  return v;
}

inline std::string ExtraStringOr(const std::unordered_map<std::string, std::string>& extra,
                                 std::string_view key,
                                 std::string default_v) {
  auto it = extra.find(std::string(key));
  if (it == extra.end()) return default_v;
  return it->second;
}

// Parse which relation should be treated as OUTER.
// Supported values (case-insensitive): "R" or "S".
// Returns true if parsed, false otherwise.
inline bool ParseOuterSide(std::string_view s, bool* out_outer_is_r) {
  if (!out_outer_is_r) return false;
  if (s.empty()) return false;
  if (s.size() == 1) {
    const char c = static_cast<char>(std::tolower(static_cast<unsigned char>(s[0])));
    if (c == 'r') {
      *out_outer_is_r = true;
      return true;
    }
    if (c == 's') {
      *out_outer_is_r = false;
      return true;
    }
  }
  std::string t(s);
  std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (t == "r" || t == "outer_r") {
    *out_outer_is_r = true;
    return true;
  }
  if (t == "s" || t == "outer_s") {
    *out_outer_is_r = false;
    return true;
  }
  return false;
}

// --------------------------
// Deterministic KD-tree range sampler
// --------------------------

// A KD-tree over points with:
//   - RangeCount(Q)
//   - RangeSample(Q,k) : i.i.d uniform samples with replacement from points in Q
//   - A streaming Cursor to enumerate points in Q without materializing all.
//
// The implementation is *static*: all points are present; node weights are
// subtree sizes.

template <int K, class T>
class KDRangeSampler {
 public:
  static_assert(K >= 1, "KDRangeSampler requires K>=1");
  using PointT = Point<K, T>;
  using BoxT = Box<K, T>;

  struct Options {
    // Leaf threshold (scan points directly).
    u32 leaf_size = 32;
  };

  KDRangeSampler() = default;

  void Reset() {
    pts_ = nullptr;
    n_ = 0;
    opt_ = Options{};
    indices_.clear();
    nodes_.clear();
    root_ = kNull;
  }

  bool Built() const noexcept { return pts_ != nullptr; }
  u32 Size() const noexcept { return n_; }

  bool Build(Span<const PointT> pts, Options opt, std::string* err) {
    Reset();
    if (opt.leaf_size == 0) opt.leaf_size = 1;
    opt_ = opt;

    if (pts.size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "KDRangeSampler::Build: too many points for u32 indexing";
      return false;
    }
    n_ = static_cast<u32>(pts.size());
    pts_ = pts.data();

    indices_.resize(static_cast<usize>(n_));
    std::iota(indices_.begin(), indices_.end(), 0u);
    nodes_.reserve(static_cast<usize>(n_ == 0 ? 0 : (2 * n_)));

    if (n_ == 0) {
      root_ = kNull;
      return true;
    }

    root_ = BuildRec(/*begin=*/0, /*end=*/n_, /*depth=*/0);
    return true;
  }

  // Exact count of points in Q.
  u64 Count(const BoxT& q) const noexcept {
    if (n_ == 0 || root_ == kNull) return 0;
    if (q.IsEmpty()) return 0;

    u64 total = 0;
    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const u32 v = stack.back();
      stack.pop_back();
      const Node& node = nodes_[static_cast<usize>(v)];

      if (!IntersectsHalfOpen<K, T>(node.bbox, q)) continue;
      if (ContainsBoxHalfOpen<K, T>(q, node.bbox)) {
        total += static_cast<u64>(node.end - node.begin);
        continue;
      }

      if (node.IsLeaf()) {
        for (u32 i = node.begin; i < node.end; ++i) {
          const u32 h = indices_[static_cast<usize>(i)];
          if (ContainsHalfOpen<K, T>(q, pts_[static_cast<usize>(h)])) {
            ++total;
          }
        }
        continue;
      }

      // Internal: test pivot + traverse children.
      const u32 mid = Mid(node);
      const u32 h = indices_[static_cast<usize>(mid)];
      if (ContainsHalfOpen<K, T>(q, pts_[static_cast<usize>(h)])) ++total;

      if (node.right != kNull) stack.push_back(node.right);
      if (node.left != kNull) stack.push_back(node.left);
    }

    return total;
  }

  // Streaming enumerator over points in Q.
  class Cursor {
   public:
    Cursor() = default;

    Cursor(const KDRangeSampler* tree, const BoxT& q)
        : tree_(tree), q_(q) {
      Reset();
    }

    void Reset() {
      stack_.clear();
      seg_pos_ = seg_end_ = 0;
      leaf_pos_ = leaf_end_ = 0;
      done_ = true;

      if (!tree_ || tree_->n_ == 0 || tree_->root_ == kNull || q_.IsEmpty()) return;

      done_ = false;
      stack_.reserve(64);
      stack_.push_back(tree_->root_);
    }

    // Produce next handle in Q. Returns false at end.
    // If candidate_checks != nullptr, we increment it for each point predicate
    // evaluation (ContainsHalfOpen) performed on boundary nodes.
    bool Next(u32* out_handle, u64* candidate_checks = nullptr) {
      if (!out_handle || done_) return false;

      while (true) {
        // 1) Emit from a fully-contained segment.
        if (seg_pos_ < seg_end_) {
          const u32 h = tree_->indices_[static_cast<usize>(seg_pos_++)];
          *out_handle = h;
          return true;
        }

        // 2) Emit from an active leaf scan.
        while (leaf_pos_ < leaf_end_) {
          const u32 h = tree_->indices_[static_cast<usize>(leaf_pos_++)];
          if (candidate_checks) ++(*candidate_checks);
          if (ContainsHalfOpen<K, T>(q_, tree_->pts_[static_cast<usize>(h)])) {
            *out_handle = h;
            return true;
          }
        }

        // 3) Advance the DFS stack.
        if (stack_.empty()) {
          done_ = true;
          return false;
        }

        const u32 v = stack_.back();
        stack_.pop_back();
        const Node& node = tree_->nodes_[static_cast<usize>(v)];

        if (!IntersectsHalfOpen<K, T>(node.bbox, q_)) continue;

        if (ContainsBoxHalfOpen<K, T>(q_, node.bbox)) {
          // Whole subtree is inside: scan the contiguous index segment.
          seg_pos_ = node.begin;
          seg_end_ = node.end;
          continue;
        }

        if (node.IsLeaf()) {
          leaf_pos_ = node.begin;
          leaf_end_ = node.end;
          continue;
        }

        // Internal node: pivot is a single point that is not included in children.
        const u32 mid = tree_->Mid(node);
        const u32 h = tree_->indices_[static_cast<usize>(mid)];
        if (candidate_checks) ++(*candidate_checks);
        if (ContainsHalfOpen<K, T>(q_, tree_->pts_[static_cast<usize>(h)])) {
          *out_handle = h;
          // Keep traversing children on future calls.
          // Push right then left for deterministic left-first DFS.
          if (node.right != kNull) stack_.push_back(node.right);
          if (node.left != kNull) stack_.push_back(node.left);
          return true;
        }

        if (node.right != kNull) stack_.push_back(node.right);
        if (node.left != kNull) stack_.push_back(node.left);
      }
    }

   private:
    const KDRangeSampler* tree_ = nullptr;
    BoxT q_{};

    std::vector<u32> stack_;

    // Current fully-contained segment [seg_pos_, seg_end_)
    u32 seg_pos_ = 0;
    u32 seg_end_ = 0;

    // Current leaf scan [leaf_pos_, leaf_end_)
    u32 leaf_pos_ = 0;
    u32 leaf_end_ = 0;

    bool done_ = true;
  };

  Cursor MakeCursor(const BoxT& q) const { return Cursor(this, q); }

  // i.i.d uniform samples (with replacement) from points in Q.
  // Returns true on success. If Q is empty (0 points), returns true and sets
  // out_handles empty.
  bool Sample(const BoxT& q, u32 k, Rng* rng, std::vector<u32>* out_handles, std::string* err) const {
    if (!rng || !out_handles) {
      if (err) *err = "KDRangeSampler::Sample: null rng/out_handles";
      return false;
    }
    out_handles->clear();
    if (k == 0) return true;
    if (n_ == 0 || root_ == kNull || q.IsEmpty()) return true;

    // Collect disjoint pieces: fully-contained subtrees (segments) + boundary points.
    std::vector<Piece> pieces;
    std::vector<u64> weights;
    u64 total_pts = 0;
    CollectPieces(q, &pieces, &weights, &total_pts);

    if (total_pts == 0) {
      return true;
    }

    sampling::AliasTable alias;
    if (!alias.BuildFromU64(Span<const u64>(weights), err)) return false;

    struct Assign {
      u32 piece;
      u32 slot;
    };

    std::vector<Assign> asg;
    asg.resize(static_cast<usize>(k));
    for (u32 i = 0; i < k; ++i) {
      const usize p = alias.Sample(rng);
      asg[static_cast<usize>(i)] = Assign{static_cast<u32>(p), i};
    }

    std::sort(asg.begin(), asg.end(), [](const Assign& a, const Assign& b) {
      if (a.piece < b.piece) return true;
      if (b.piece < a.piece) return false;
      return a.slot < b.slot;
    });

    out_handles->resize(static_cast<usize>(k));

    usize i = 0;
    while (i < asg.size()) {
      const u32 pidx = asg[i].piece;
      usize j = i + 1;
      while (j < asg.size() && asg[j].piece == pidx) ++j;

      const Piece& pc = pieces[static_cast<usize>(pidx)];

      if (pc.kind == Piece::Kind::Point) {
        for (usize t = i; t < j; ++t) {
          (*out_handles)[static_cast<usize>(asg[t].slot)] = pc.handle;
        }
      } else {
        const u32 len = pc.end - pc.begin;
        SJS_DASSERT(len > 0);
        for (usize t = i; t < j; ++t) {
          const u32 off = static_cast<u32>(rng->UniformU64(static_cast<u64>(len)));
          const u32 h = indices_[static_cast<usize>(pc.begin + off)];
          (*out_handles)[static_cast<usize>(asg[t].slot)] = h;
        }
      }

      i = j;
    }

    return true;
  }

 private:
  static constexpr u32 kNull = std::numeric_limits<u32>::max();

  struct Node {
    u32 left = kNull;
    u32 right = kNull;
    u32 begin = 0;
    u32 end = 0;
    BoxT bbox = BoxT::Empty();

    bool IsLeaf() const noexcept { return left == kNull && right == kNull; }
  };

  struct Piece {
    enum class Kind : u8 { Segment = 0, Point = 1 };
    Kind kind{Kind::Point};
    u32 begin = 0;
    u32 end = 0;
    u32 handle = 0;
    u64 weight = 0;
  };

  u32 Mid(const Node& n) const noexcept {
    // BuildRec uses this exact mid definition for pivot placement.
    const u32 sz = n.end - n.begin;
    return n.begin + (sz / 2);
  }

  BoxT BoxOfPoint(const PointT& p) const noexcept {
    BoxT b = BoxT::Empty();
    b.ExpandToIncludePoint(p);
    return b;
  }

  u32 BuildRec(u32 begin, u32 end, int depth) {
    const u32 node_id = static_cast<u32>(nodes_.size());
    nodes_.push_back(Node{});
    Node& node = nodes_.back();
    node.begin = begin;
    node.end = end;

    const u32 sz = end - begin;
    if (sz <= opt_.leaf_size) {
      // Leaf: compute bbox by scanning points.
      BoxT bb = BoxT::Empty();
      for (u32 i = begin; i < end; ++i) {
        const u32 h = indices_[static_cast<usize>(i)];
        bb.ExpandToIncludePoint(pts_[static_cast<usize>(h)]);
      }
      node.bbox = bb;
      node.left = kNull;
      node.right = kNull;
      return node_id;
    }

    const int axis = depth % K;
    const u32 mid = begin + (sz / 2);

    auto comp = [&](u32 a, u32 b) {
      const T va = pts_[static_cast<usize>(a)].v[static_cast<usize>(axis)];
      const T vb = pts_[static_cast<usize>(b)].v[static_cast<usize>(axis)];
      if (va < vb) return true;
      if (vb < va) return false;
      return a < b;
    };

    std::nth_element(indices_.begin() + static_cast<usize>(begin),
                     indices_.begin() + static_cast<usize>(mid),
                     indices_.begin() + static_cast<usize>(end),
                     comp);

    // Children exclude the pivot at mid.
    if (mid > begin) node.left = BuildRec(begin, mid, depth + 1);
    if (mid + 1 < end) node.right = BuildRec(mid + 1, end, depth + 1);

    // bbox = union(children bbox, pivot point)
    BoxT bb = BoxOfPoint(pts_[static_cast<usize>(indices_[static_cast<usize>(mid)])]);
    if (node.left != kNull) bb.ExpandToIncludeBox(nodes_[static_cast<usize>(node.left)].bbox);
    if (node.right != kNull) bb.ExpandToIncludeBox(nodes_[static_cast<usize>(node.right)].bbox);
    node.bbox = bb;

    return node_id;
  }

  void CollectPieces(const BoxT& q,
                     std::vector<Piece>* out_pieces,
                     std::vector<u64>* out_weights,
                     u64* out_total_pts) const {
    SJS_DASSERT(out_pieces != nullptr);
    SJS_DASSERT(out_weights != nullptr);
    SJS_DASSERT(out_total_pts != nullptr);

    out_pieces->clear();
    out_weights->clear();
    *out_total_pts = 0;

    if (n_ == 0 || root_ == kNull || q.IsEmpty()) return;

    std::vector<u32> stack;
    stack.reserve(64);
    stack.push_back(root_);

    while (!stack.empty()) {
      const u32 v = stack.back();
      stack.pop_back();
      const Node& node = nodes_[static_cast<usize>(v)];

      if (!IntersectsHalfOpen<K, T>(node.bbox, q)) continue;

      if (ContainsBoxHalfOpen<K, T>(q, node.bbox)) {
        const u32 len = node.end - node.begin;
        if (len > 0) {
          out_pieces->push_back(Piece{Piece::Kind::Segment, node.begin, node.end, /*handle=*/0u, /*weight=*/len});
          out_weights->push_back(static_cast<u64>(len));
          *out_total_pts += static_cast<u64>(len);
        }
        continue;
      }

      if (node.IsLeaf()) {
        for (u32 i = node.begin; i < node.end; ++i) {
          const u32 h = indices_[static_cast<usize>(i)];
          if (ContainsHalfOpen<K, T>(q, pts_[static_cast<usize>(h)])) {
            out_pieces->push_back(Piece{Piece::Kind::Point, 0u, 0u, h, 1u});
            out_weights->push_back(1ULL);
            *out_total_pts += 1ULL;
          }
        }
        continue;
      }

      // Internal node: add pivot point if it falls in Q, then explore children.
      const u32 mid = Mid(node);
      const u32 h = indices_[static_cast<usize>(mid)];
      if (ContainsHalfOpen<K, T>(q, pts_[static_cast<usize>(h)])) {
        out_pieces->push_back(Piece{Piece::Kind::Point, 0u, 0u, h, 1u});
        out_weights->push_back(1ULL);
        *out_total_pts += 1ULL;
      }

      if (node.right != kNull) stack.push_back(node.right);
      if (node.left != kNull) stack.push_back(node.left);
    }
  }

  const PointT* pts_ = nullptr;
  u32 n_ = 0;
  Options opt_{};

  std::vector<u32> indices_;  // permutation of [0,n)
  std::vector<Node> nodes_;
  u32 root_ = kNull;
};

// --------------------------
// Shared SIRS state (built index + orientation)
// --------------------------

template <int Dim, class T>
struct SIRSState {
  static_assert(Dim >= 1, "SIRSState requires Dim>=1");

  using DatasetT = Dataset<Dim, T>;
  using RelT = Relation<Dim, T>;
  using BoxT = Box<Dim, T>;

  static constexpr int K = 2 * Dim;
  using EmbPointT = EmbeddedPoint<Dim, T>;
  using EmbBoxT = EmbeddedBox<Dim, T>;
  using KD = KDRangeSampler<K, T>;

  const DatasetT* ds = nullptr;
  const RelT* outer = nullptr;
  const RelT* inner = nullptr;
  bool outer_is_r = true;

  DomainBounds<Dim, T> inner_domain;

  KD kd;
  typename KD::Options kd_opt;

  std::vector<EmbPointT> inner_points;  // phi(b) in original order (handle = inner index)

  // Phase-1 weights (deg per OUTER)
  std::vector<u64> deg;
  sampling::AliasTable outer_alias;
  u64 W = 0;
  bool built = false;
  bool weights_valid = false;

  void Reset() {
    ds = nullptr;
    outer = nullptr;
    inner = nullptr;
    outer_is_r = true;

    inner_domain = DomainBounds<Dim, T>{};

    kd.Reset();
    kd_opt = typename KD::Options{};
    inner_points.clear();

    deg.clear();
    outer_alias.Clear();
    W = 0;
    built = false;
    weights_valid = false;
  }

  // Build the embedding and the point index for the chosen INNER relation.
  bool BuildIndex(const DatasetT& dataset, const Config& cfg, PhaseRecorder* phases, std::string* err) {
    auto scoped = phases ? phases->Scoped("sirs_build_index") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    ds = &dataset;

    // Determine OUTER side.
    // Priority:
    //   1) cfg.run.extra["sirs_outer"] = "R" or "S".
    //   2) otherwise choose the smaller relation as OUTER (reduces number of range-count queries).
    bool outer_r = (dataset.R.Size() <= dataset.S.Size());
    {
      const std::string v = ExtraStringOr(cfg.run.extra, "sirs_outer", "");
      bool parsed = false;
      if (!v.empty()) parsed = ParseOuterSide(v, &outer_r);
      (void)parsed;
    }

    outer_is_r = outer_r;
    outer = outer_is_r ? &dataset.R : &dataset.S;
    inner = outer_is_r ? &dataset.S : &dataset.R;

    // KD-tree leaf size.
    kd_opt.leaf_size = ExtraU32Or(cfg.run.extra, "sirs_leaf_size", kd_opt.leaf_size);
    if (kd_opt.leaf_size == 0) kd_opt.leaf_size = 1;

    // Domain bounds from INNER only (sufficient for the reduction).
    inner_domain = DomainBounds<Dim, T>::FromBoxes(Span<const BoxT>(inner->boxes));

    // Build embedded points for INNER in original order.
    {
      auto _ = phases ? phases->Scoped("sirs_embed_inner") : PhaseRecorder::ScopedPhase(nullptr, "");
      inner_points.resize(inner->Size());
      for (usize i = 0; i < inner->Size(); ++i) {
        inner_points[i] = EmbedLowerUpper<Dim, T>(inner->boxes[i]);
      }
    }

    // Build KD-tree over embedded points.
    {
      auto _ = phases ? phases->Scoped("sirs_build_kd") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!kd.Build(Span<const EmbPointT>(inner_points), kd_opt, err)) return false;
    }

    built = true;
    weights_valid = false;
    W = 0;
    deg.clear();
    outer_alias.Clear();
    return true;
  }

  // Convert an OUTER rectangle a to the embedded query range Q(a).
  EmbBoxT QueryForOuter(const BoxT& a) const noexcept {
    return MakeIntersectQueryRange<Dim, T>(a, inner_domain);
  }

  // Map (outer index, inner handle) to (R,S) ids.
  PairId MakePair(u32 outer_idx, u32 inner_handle) const noexcept {
    const Id out_id = outer->GetId(static_cast<usize>(outer_idx));
    const Id in_id = inner->GetId(static_cast<usize>(inner_handle));
    if (outer_is_r) {
      // OUTER is R, INNER is S.
      return PairId{out_id, in_id};
    }
    // OUTER is S, INNER is R.
    return PairId{in_id, out_id};
  }

  // Phase 1: compute deg for each OUTER and build alias table.
  bool ComputeDegreesAndAlias(const Config& cfg, PhaseRecorder* phases, std::string* err) {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("sirs_phase1_degrees") : PhaseRecorder::ScopedPhase(nullptr, "");

    if (!built || !outer || !inner) {
      if (err) *err = "SIRSState::ComputeDegreesAndAlias: call BuildIndex() first";
      return false;
    }

    const usize n_outer = outer->Size();
    deg.assign(n_outer, 0ULL);
    W = 0;

    if (n_outer == 0 || inner->Size() == 0) {
      outer_alias.Clear();
      weights_valid = true;
      return true;
    }

    for (usize i = 0; i < n_outer; ++i) {
      const BoxT& a = outer->boxes[i];
      const EmbBoxT q = QueryForOuter(a);
      const u64 w = q.IsEmpty() ? 0ULL : kd.Count(q);
      deg[i] = w;

      if (w > 0) {
        if (W > std::numeric_limits<u64>::max() - w) {
          if (err) *err = "SIRSState::ComputeDegreesAndAlias: |J| overflowed u64";
          return false;
        }
        W += w;
      }
    }

    // Alias table over OUTER indices.
    outer_alias.Clear();
    if (!outer_alias.BuildFromU64(Span<const u64>(deg), err)) return false;

    weights_valid = true;
    return true;
  }
};

// --------------------------
// Deterministic join enumerator for SIRS (outer-loop over OUTER, range-scan INNER)
// --------------------------

template <int Dim, class T>
class SIRSJoinEnumerator final : public IJoinEnumerator {
 public:
  using State = SIRSState<Dim, T>;
  using BoxT = typename State::BoxT;
  using EmbBoxT = typename State::EmbBoxT;
  using KD = typename State::KD;

  explicit SIRSJoinEnumerator(const State* st) : st_(st) {
    Reset();
  }

  void Reset() override {
    stats_.Reset();
    outer_i_ = 0;
    cursor_ = typename KD::Cursor();
    cursor_valid_ = false;

    if (st_ && st_->built && st_->outer) {
      stats_.num_events = static_cast<u64>(st_->outer->Size());
    }
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

  bool Next(PairId* out) override {
    if (!out) return false;
    if (!st_ || !st_->built || !st_->outer || !st_->inner) return false;

    const usize n_outer = st_->outer->Size();

    while (true) {
      if (!cursor_valid_) {
        if (outer_i_ >= n_outer) return false;

        // Initialize cursor for current OUTER rectangle.
        const BoxT& a = st_->outer->boxes[outer_i_];
        const EmbBoxT q = st_->QueryForOuter(a);
        cursor_ = st_->kd.MakeCursor(q);
        cursor_valid_ = true;
      }

      u32 inner_handle = 0;
      if (cursor_.Next(&inner_handle, &stats_.candidate_checks)) {
        *out = st_->MakePair(static_cast<u32>(outer_i_), inner_handle);
        ++stats_.output_pairs;
        return true;
      }

      // Current OUTER exhausted.
      cursor_valid_ = false;
      ++outer_i_;
    }
  }

 private:
  const State* st_ = nullptr;

  usize outer_i_ = 0;
  typename KD::Cursor cursor_{};
  bool cursor_valid_ = false;

  join::JoinStats stats_{};
};

}  // namespace detail

// --------------------------
// Baseline: Sampling variant
// --------------------------

template <int Dim, class T = Scalar>
class SIRSSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::SIRS; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "sirs_sampling"; }

  void Reset() override { st_.Reset(); }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    return st_.BuildIndex(ds, cfg, phases, err);
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic
    if (!out) {
      if (err) *err = "SIRSSamplingBaseline::Count: out is null";
      return false;
    }
    if (!st_.built) {
      if (err) *err = "SIRSSamplingBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");

    if (!st_.ComputeDegreesAndAlias(cfg, phases, err)) return false;

    *out = MakeExactCount(st_.W);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!rng) {
      if (err) *err = "SIRSSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "SIRSSamplingBaseline::Sample: out is null";
      return false;
    }
    if (!st_.built) {
      if (err) *err = "SIRSSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!st_.weights_valid) {
      if (err) *err = "SIRSSamplingBaseline::Sample: call Count() first";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    const u64 t = cfg.run.t;
    if (t == 0 || st_.W == 0) return true;

    if (t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "SIRSSamplingBaseline::Sample: t exceeds u32 limit";
      return false;
    }
    const u32 t32 = static_cast<u32>(t);

    auto scoped = phases ? phases->Scoped("sample") : PhaseRecorder::ScopedPhase(nullptr, "");

    // Phase 2: assign each sample slot to an OUTER index by alias.
    struct Assign {
      u32 outer;
      u32 slot;
    };

    std::vector<Assign> asg;
    asg.resize(static_cast<usize>(t32));
    for (u32 i = 0; i < t32; ++i) {
      const usize oi = st_.outer_alias.Sample(rng);
      asg[static_cast<usize>(i)] = Assign{static_cast<u32>(oi), i};
    }

    std::sort(asg.begin(), asg.end(), [](const Assign& a, const Assign& b) {
      if (a.outer < b.outer) return true;
      if (b.outer < a.outer) return false;
      return a.slot < b.slot;
    });

    out->pairs.resize(static_cast<usize>(t32));

    // Phase 3: for each OUTER group, sample INNER within its query range.
    std::vector<u32> inner_handles;
    inner_handles.reserve(256);

    usize i = 0;
    while (i < asg.size()) {
      const u32 outer_idx = asg[i].outer;
      usize j = i + 1;
      while (j < asg.size() && asg[j].outer == outer_idx) ++j;

      const u32 k = static_cast<u32>(j - i);

      // Defensive: alias should never pick a zero-degree OUTER unless W==0.
      if (outer_idx >= st_.deg.size() || st_.deg[static_cast<usize>(outer_idx)] == 0) {
        if (err) *err = "SIRSSamplingBaseline::Sample: selected an outer with deg=0 (alias table bug?)";
        return false;
      }

      const auto& a = st_.outer->boxes[static_cast<usize>(outer_idx)];
      const auto q = st_.QueryForOuter(a);

      inner_handles.clear();
      if (!st_.kd.Sample(q, k, rng, &inner_handles, err)) return false;
      if (inner_handles.size() != k) {
        if (err) *err = "SIRSSamplingBaseline::Sample: kd.Sample produced wrong number of samples";
        return false;
      }

      for (u32 tslot = 0; tslot < k; ++tslot) {
        const u32 slot = asg[i + tslot].slot;
        out->pairs[static_cast<usize>(slot)] = st_.MakePair(outer_idx, inner_handles[static_cast<usize>(tslot)]);
      }

      i = j;
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!st_.built) {
      if (err) *err = "SIRSSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::SIRSJoinEnumerator<Dim, T>>(&st_);
  }

 private:
  detail::SIRSState<Dim, T> st_;
};

}  // namespace sirs
}  // namespace baselines
}  // namespace sjs
