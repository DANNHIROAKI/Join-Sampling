#pragma once
// utils/kdtree.hpp (modified)
//
// A lightweight static KD-tree supporting orthogonal range COUNT and range SEARCH
// in configurable dimensionality D (for this repository typically D = 2 * Dim
// under the standard box-lifting embedding).
//
// Design notes
// ------------
// - Header-only, no global RNG.
// - Points are stored in-place (std::nth_element) so each subtree corresponds
//   to a contiguous index interval [l, r) in the points array.
// - Node bounding boxes are computed bottom-up by merging children (O(n) total),
//   instead of scanning all points at each node (O(n log n)).
// - Range semantics are CLOSED (inclusive): qlo[d] <= x[d] <= qhi[d].
//
// This matches the original RS-over-SRJ code style and allows expressing
// strict inequalities using std::nextafter (see isrjs_kds.hpp modifications).

#include <array>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

// --------------------------
// KDPoint
// --------------------------
// A D-dimensional point with a payload index.
// In our integration, payload will store an index into relation S.

template <int D, class T>
struct KDPoint {
  static_assert(D >= 1, "KDPoint<D,T>: D must be >= 1");
  std::array<T, D> coord{};
  std::uint32_t idx = 0;
};

// --------------------------
// KDTree
// --------------------------

template <int D, class T>
class KDTree {
 public:
  static_assert(D >= 1, "KDTree<D,T>: D must be >= 1");

  using PointT = KDPoint<D, T>;

  static constexpr std::uint32_t kNull = std::numeric_limits<std::uint32_t>::max();

  struct Node {
    // Children in nodes_ (kNull if absent)
    std::uint32_t left = kNull;
    std::uint32_t right = kNull;

    // Contiguous interval of points_ covered by this subtree: [l, r)
    std::uint32_t l = 0;
    std::uint32_t r = 0;

    // Median index in points_ (the actual point stored at this node)
    std::uint32_t mid = 0;

    // Split axis used to choose the median
    std::uint8_t axis = 0;

    // Bounding box of this subtree (inclusive min/max)
    std::array<T, D> lo{};
    std::array<T, D> hi{};

    std::uint32_t Size() const noexcept { return r - l; }
    bool IsLeaf() const noexcept { return left == kNull && right == kNull; }
  };

  struct SearchItem {
    std::uint32_t node = kNull;
    bool fully_contained = false;  // true => subtree fully contained in query
  };

  KDTree() = default;

  void Clear() {
    points_.clear();
    nodes_.clear();
    root_ = kNull;
  }

  bool Empty() const noexcept { return root_ == kNull; }
  std::uint32_t Size() const noexcept { return static_cast<std::uint32_t>(points_.size()); }
  std::uint32_t Root() const noexcept { return root_; }

  const std::vector<PointT>& Points() const noexcept { return points_; }
  const std::vector<Node>& Nodes() const noexcept { return nodes_; }

  const Node& GetNode(std::uint32_t node) const noexcept { return nodes_[node]; }

  // Payload stored at the median point of a node.
  std::uint32_t NodePointPayload(std::uint32_t node) const noexcept {
    return points_[nodes_[node].mid].idx;
  }

  // Build from a point vector (moved in). After Build(), points_ are permuted.
  void Build(std::vector<PointT>&& pts) {
    Clear();
    points_ = std::move(pts);
    if (points_.empty()) return;
    nodes_.reserve(points_.size());
    root_ = BuildRec(/*l=*/0U, /*r=*/static_cast<std::uint32_t>(points_.size()), /*depth=*/0);
  }

  // Range COUNT (closed intervals): qlo[d] <= x[d] <= qhi[d].
  std::uint64_t Count(const std::array<T, D>& qlo, const std::array<T, D>& qhi) const {
    if (root_ == kNull) return 0;
    return CountRec(root_, qlo, qhi);
  }

  // Range SEARCH: returns a decomposition of the query result set into
  //  - fully contained subtrees (fully_contained=true), and
  //  - single points stored at nodes (fully_contained=false).
  //
  // This structure is designed for uniform sampling over the query set:
  // weight(subtree)=subtree_size, weight(point)=1.
  void Search(const std::array<T, D>& qlo,
              const std::array<T, D>& qhi,
              std::vector<SearchItem>* out) const {
    if (!out) return;
    out->clear();
    if (root_ == kNull) return;
    SearchRec(root_, qlo, qhi, out);
  }

 private:
  std::uint32_t root_ = kNull;
  std::vector<PointT> points_;
  std::vector<Node> nodes_;

  static inline bool PointInRange(const PointT& p,
                                  const std::array<T, D>& qlo,
                                  const std::array<T, D>& qhi) noexcept {
    for (int d = 0; d < D; ++d) {
      const std::size_t j = static_cast<std::size_t>(d);
      const T v = p.coord[j];
      if (v < qlo[j] || v > qhi[j]) return false;
    }
    return true;
  }

  static inline bool Disjoint(const Node& nd,
                              const std::array<T, D>& qlo,
                              const std::array<T, D>& qhi) noexcept {
    for (int d = 0; d < D; ++d) {
      const std::size_t j = static_cast<std::size_t>(d);
      if (nd.hi[j] < qlo[j] || nd.lo[j] > qhi[j]) return true;
    }
    return false;
  }

  static inline bool Contained(const Node& nd,
                               const std::array<T, D>& qlo,
                               const std::array<T, D>& qhi) noexcept {
    for (int d = 0; d < D; ++d) {
      const std::size_t j = static_cast<std::size_t>(d);
      if (nd.lo[j] < qlo[j]) return false;
      if (nd.hi[j] > qhi[j]) return false;
    }
    return true;
  }

  std::uint32_t BuildRec(std::uint32_t l, std::uint32_t r, int depth) {
    if (l >= r) return kNull;

    const std::uint32_t n = r - l;
    const std::uint32_t mid = l + (n / 2U);
    const std::uint8_t axis = static_cast<std::uint8_t>(depth % D);

    auto comp = [axis](const PointT& a, const PointT& b) {
      return a.coord[static_cast<std::size_t>(axis)] < b.coord[static_cast<std::size_t>(axis)];
    };

    // Partition so that points_[mid] is the median w.r.t. axis within [l, r).
    std::nth_element(points_.begin() + static_cast<std::ptrdiff_t>(l),
                     points_.begin() + static_cast<std::ptrdiff_t>(mid),
                     points_.begin() + static_cast<std::ptrdiff_t>(r),
                     comp);

    // Allocate node.
    const std::uint32_t node_idx = static_cast<std::uint32_t>(nodes_.size());
    nodes_.push_back(Node{});

    // Recurse.
    const std::uint32_t left = BuildRec(l, mid, depth + 1);
    const std::uint32_t right = BuildRec(mid + 1U, r, depth + 1);

    // Fill node fields (do NOT keep references across recursive calls).
    Node& nd = nodes_[node_idx];
    nd.l = l;
    nd.r = r;
    nd.mid = mid;
    nd.axis = axis;
    nd.left = left;
    nd.right = right;

    // Bottom-up bbox merge: start with the median point.
    for (int d = 0; d < D; ++d) {
      const std::size_t j = static_cast<std::size_t>(d);
      nd.lo[j] = points_[mid].coord[j];
      nd.hi[j] = points_[mid].coord[j];
    }

    auto merge_child = [&](std::uint32_t child) {
      if (child == kNull) return;
      const Node& c = nodes_[child];
      for (int d = 0; d < D; ++d) {
        const std::size_t j = static_cast<std::size_t>(d);
        if (c.lo[j] < nd.lo[j]) nd.lo[j] = c.lo[j];
        if (c.hi[j] > nd.hi[j]) nd.hi[j] = c.hi[j];
      }
    };

    merge_child(left);
    merge_child(right);

    return node_idx;
  }

  std::uint64_t CountRec(std::uint32_t node,
                         const std::array<T, D>& qlo,
                         const std::array<T, D>& qhi) const {
    if (node == kNull) return 0;
    const Node& nd = nodes_[node];

    if (Disjoint(nd, qlo, qhi)) return 0;
    if (Contained(nd, qlo, qhi)) return static_cast<std::uint64_t>(nd.Size());

    std::uint64_t cnt = 0;
    if (PointInRange(points_[nd.mid], qlo, qhi)) ++cnt;

    if (nd.left != kNull) cnt += CountRec(nd.left, qlo, qhi);
    if (nd.right != kNull) cnt += CountRec(nd.right, qlo, qhi);
    return cnt;
  }

  void SearchRec(std::uint32_t node,
                 const std::array<T, D>& qlo,
                 const std::array<T, D>& qhi,
                 std::vector<SearchItem>* out) const {
    if (node == kNull) return;
    const Node& nd = nodes_[node];

    if (Disjoint(nd, qlo, qhi)) return;

    // Leaf => the single point is in range (by not-disjoint).
    if (nd.IsLeaf()) {
      out->push_back(SearchItem{node, /*fully_contained=*/false});
      return;
    }

    // Internal node: include its median point if it is in range.
    if (PointInRange(points_[nd.mid], qlo, qhi)) {
      out->push_back(SearchItem{node, /*fully_contained=*/false});
    }

    // Recurse on children; if a child subtree is fully contained, add it as one block.
    if (nd.left != kNull) {
      const Node& c = nodes_[nd.left];
      if (Contained(c, qlo, qhi)) {
        out->push_back(SearchItem{nd.left, /*fully_contained=*/true});
      } else {
        SearchRec(nd.left, qlo, qhi, out);
      }
    }

    if (nd.right != kNull) {
      const Node& c = nodes_[nd.right];
      if (Contained(c, qlo, qhi)) {
        out->push_back(SearchItem{nd.right, /*fully_contained=*/true});
      } else {
        SearchRec(nd.right, qlo, qhi, out);
      }
    }
  }
};