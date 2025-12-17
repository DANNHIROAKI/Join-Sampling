// src/index/index_factory_2d.cpp
//
// A small, self-contained factory for Dim=2 index objects.
//
// Why a factory?
// -------------
// Many baselines share the same "index-style" building blocks (AABB-tree,
// interval tree, KD-tree on embeddings, etc.). For micro-benchmarks and
// verification utilities, it is often convenient to select an index at runtime
// by a string name.
//
// This TU implements a minimal type-erased wrapper `BoxIndex2D` that:
//  - can Build() over a set of 2D boxes,
//  - can QueryIntersecting() for a 2D query box,
//  - can CountIntersecting().
//
// Notes
// -----
// - Some indices (interval/KD/Range) are primarily candidate generators; for a
//   unified interface we filter candidates with an exact box-box intersection
//   check (half-open semantics).
// - This file lives under src/ because it is not essential to baseline APIs; it
//   is mainly for tooling and future extensions.

#include "sjs/index/aabb_tree.h"
#include "sjs/index/interval_tree.h"
#include "sjs/index/kd_tree.h"
#include "sjs/index/r_tree.h"
#include "sjs/index/range_tree.h"

#include "sjs/core/assert.h"
#include "sjs/core/types.h"
#include "sjs/geometry/predicates.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace index {

namespace {

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

inline bool EqualsIgnoreCase(std::string_view a, std::string_view b) {
  if (a.size() != b.size()) return false;
  for (usize i = 0; i < a.size(); ++i) {
    char ca = a[i];
    char cb = b[i];
    if (ca >= 'A' && ca <= 'Z') ca = static_cast<char>(ca - 'A' + 'a');
    if (cb >= 'A' && cb <= 'Z') cb = static_cast<char>(cb - 'A' + 'a');
    if (ca != cb) return false;
  }
  return true;
}

}  // namespace

// --------------------------
// Type-erased index wrapper
// --------------------------

struct BoxIndex2D {
  using BoxT = Box<2, Scalar>;

  // Type erasure: a raw pointer + function table.
  void* impl = nullptr;
  void (*destroy)(void*) = nullptr;
  void (*build)(void*, Span<const BoxT>) = nullptr;
  void (*query)(const void*, const BoxT&, std::vector<usize>*) = nullptr;
  u64 (*count)(const void*, const BoxT&) = nullptr;
  std::string name;

  BoxIndex2D() = default;
  ~BoxIndex2D() { Reset(); }

  BoxIndex2D(BoxIndex2D&& other) noexcept { *this = std::move(other); }
  BoxIndex2D& operator=(BoxIndex2D&& other) noexcept {
    if (this == &other) return *this;
    Reset();
    impl = other.impl;
    destroy = other.destroy;
    build = other.build;
    query = other.query;
    count = other.count;
    name = std::move(other.name);

    other.impl = nullptr;
    other.destroy = nullptr;
    other.build = nullptr;
    other.query = nullptr;
    other.count = nullptr;
    other.name.clear();
    return *this;
  }

  BoxIndex2D(const BoxIndex2D&) = delete;
  BoxIndex2D& operator=(const BoxIndex2D&) = delete;

  void Reset() {
    if (impl && destroy) destroy(impl);
    impl = nullptr;
    destroy = nullptr;
    build = nullptr;
    query = nullptr;
    count = nullptr;
    name.clear();
  }

  bool Valid() const noexcept {
    return impl != nullptr && destroy != nullptr && build != nullptr && query != nullptr && count != nullptr;
  }

  void Build(Span<const BoxT> boxes) {
    SJS_ASSERT(Valid());
    build(impl, boxes);
  }

  void QueryIntersecting(const BoxT& q, std::vector<usize>* out) const {
    SJS_ASSERT(Valid());
    SJS_ASSERT(out != nullptr);
    query(impl, q, out);
  }

  u64 CountIntersecting(const BoxT& q) const {
    SJS_ASSERT(Valid());
    return count(impl, q);
  }
};

// --------------------------
// Concrete implementations
// --------------------------

class IntervalTreeIndex2DImpl {
 public:
  using BoxT = BoxIndex2D::BoxT;

  void Build(Span<const BoxT> boxes) {
    boxes_ = boxes.data();
    n_ = boxes.size();

    typename IntervalTreeOnAxis<2, Scalar>::Options opt;
    opt.axis = 1;  // y-axis
    // Keep default IntervalTree options (randomized treap). For reproducibility
    // you may want to seed IntervalTree::Options::seed.
    tree_.Init(n_, opt);

    // Insert all boxes.
    for (usize i = 0; i < n_; ++i) {
      // IntervalTree uses handles in [0,capacity). We use i as handle.
      std::string err;
      const bool ok = tree_.InsertBox(i, boxes_[i], &err);
      (void)ok;
      SJS_ASSERT_MSG(ok, err.c_str());
    }
  }

  void Query(const BoxT& q, std::vector<usize>* out) const {
    out->clear();
    if (n_ == 0) return;
    if (q.IsEmpty()) return;

    tree_.QueryBox(q, [&](usize h) {
      // Filter to exact 2D intersection.
      if (h < n_ && IntersectsHalfOpen<2, Scalar>(boxes_[h], q)) {
        out->push_back(h);
      }
      return true;
    });
  }

  u64 Count(const BoxT& q) const {
    if (n_ == 0 || q.IsEmpty()) return 0;
    u64 c = 0;
    tree_.QueryBox(q, [&](usize h) {
      if (h < n_ && IntersectsHalfOpen<2, Scalar>(boxes_[h], q)) ++c;
      return true;
    });
    return c;
  }

 private:
  const BoxT* boxes_{nullptr};
  usize n_{0};
  IntervalTreeOnAxis<2, Scalar> tree_{};
};

class AABBTreeIndex2DImpl {
 public:
  using BoxT = BoxIndex2D::BoxT;

  void Build(Span<const BoxT> boxes) {
    boxes_ = boxes.data();
    n_ = boxes.size();
    tree_.Build(boxes);
  }

  void Query(const BoxT& q, std::vector<usize>* out) const {
    out->clear();
    if (n_ == 0 || q.IsEmpty()) return;
    tree_.QueryIntersections(q, [&](usize idx) {
      out->push_back(idx);
      return true;
    });
  }

  u64 Count(const BoxT& q) const {
    if (n_ == 0 || q.IsEmpty()) return 0;
    return tree_.CountIntersections(q);
  }

 private:
  const BoxT* boxes_{nullptr};
  usize n_{0};
  AABBTree<2, Scalar> tree_{};
};

class RTreeIndex2DImpl {
 public:
  using BoxT = BoxIndex2D::BoxT;

  void Build(Span<const BoxT> boxes) {
    boxes_ = boxes.data();
    n_ = boxes.size();
    tree_.Build(boxes);
  }

  void Query(const BoxT& q, std::vector<usize>* out) const {
    out->clear();
    if (n_ == 0 || q.IsEmpty()) return;
    tree_.QueryIntersections(q, [&](usize idx) {
      out->push_back(idx);
      return true;
    });
  }

  u64 Count(const BoxT& q) const {
    if (n_ == 0 || q.IsEmpty()) return 0;
    return tree_.CountIntersections(q);
  }

 private:
  const BoxT* boxes_{nullptr};
  usize n_{0};
  RTree<2, Scalar> tree_{};
};

// KD/Range tree work on embedded points. For Dim=2 join, a common embedding is:
//   p(r) = (L_y(r), R_y(r))
// A query box q overlaps r in y iff:
//   L_y(r) < R_y(q)  AND  R_y(r) > L_y(q)
// which is a 2D orthogonal range query.

inline Box<2, Scalar> YOverlapQueryToPointRange(const Box<2, Scalar>& q) {
  const Scalar ninf = -std::numeric_limits<Scalar>::infinity();
  const Scalar pinf = std::numeric_limits<Scalar>::infinity();

  // Strict: R_y(r) > L_y(q)  ->  R_y(r) >= nextafter(L_y(q), +inf)
  const Scalar y_lo = std::nextafter(q.lo.v[1], pinf);
  const Scalar y_hi = pinf;

  // Strict: L_y(r) < R_y(q) -> L_y(r) in (-inf, R_y(q)) => [ninf, R_y(q))
  const Scalar x_lo = ninf;
  const Scalar x_hi = q.hi.v[1];

  Box<2, Scalar> range;
  range.lo.v[0] = x_lo;
  range.hi.v[0] = x_hi;
  range.lo.v[1] = y_lo;
  range.hi.v[1] = y_hi;
  return range;
}

class KDTreeIndex2DImpl {
 public:
  using BoxT = BoxIndex2D::BoxT;
  using PointT = Point<2, Scalar>;

  void Build(Span<const BoxT> boxes) {
    boxes_ = boxes.data();
    n_ = boxes.size();

    pts_.clear();
    pts_.resize(n_);
    for (usize i = 0; i < n_; ++i) {
      pts_[i].v[0] = boxes_[i].lo.v[1];
      pts_[i].v[1] = boxes_[i].hi.v[1];
    }
    tree_.Build(Span<const PointT>(pts_.data(), pts_.size()));
  }

  void Query(const BoxT& q, std::vector<usize>* out) const {
    out->clear();
    if (n_ == 0 || q.IsEmpty()) return;

    const Box<2, Scalar> range = YOverlapQueryToPointRange(q);
    if (range.IsEmpty()) return;

    tree_.RangeQuery(range, [&](usize idx) {
      // idx is the point index, which matches box index.
      if (idx < n_ && IntersectsHalfOpen<2, Scalar>(boxes_[idx], q)) {
        out->push_back(idx);
      }
      return true;
    });
  }

  u64 Count(const BoxT& q) const {
    if (n_ == 0 || q.IsEmpty()) return 0;
    const Box<2, Scalar> range = YOverlapQueryToPointRange(q);
    if (range.IsEmpty()) return 0;

    u64 c = 0;
    tree_.RangeQuery(range, [&](usize idx) {
      if (idx < n_ && IntersectsHalfOpen<2, Scalar>(boxes_[idx], q)) ++c;
      return true;
    });
    return c;
  }

 private:
  const BoxT* boxes_{nullptr};
  usize n_{0};
  std::vector<PointT> pts_;
  KDTree<2, Scalar> tree_{};
};

class RangeTreeIndex2DImpl {
 public:
  using BoxT = BoxIndex2D::BoxT;
  using PointT = Point<2, Scalar>;

  void Build(Span<const BoxT> boxes) {
    boxes_ = boxes.data();
    n_ = boxes.size();

    pts_.clear();
    pts_.resize(n_);
    for (usize i = 0; i < n_; ++i) {
      pts_[i].v[0] = boxes_[i].lo.v[1];
      pts_[i].v[1] = boxes_[i].hi.v[1];
    }
    tree_.Build(Span<const PointT>(pts_.data(), pts_.size()));
  }

  void Query(const BoxT& q, std::vector<usize>* out) const {
    out->clear();
    if (n_ == 0 || q.IsEmpty()) return;

    const Box<2, Scalar> range = YOverlapQueryToPointRange(q);
    if (range.IsEmpty()) return;

    tree_.RangeQuery(range, [&](usize idx) {
      if (idx < n_ && IntersectsHalfOpen<2, Scalar>(boxes_[idx], q)) {
        out->push_back(idx);
      }
      return true;
    });
  }

  u64 Count(const BoxT& q) const {
    if (n_ == 0 || q.IsEmpty()) return 0;
    const Box<2, Scalar> range = YOverlapQueryToPointRange(q);
    if (range.IsEmpty()) return 0;

    u64 c = 0;
    tree_.RangeQuery(range, [&](usize idx) {
      if (idx < n_ && IntersectsHalfOpen<2, Scalar>(boxes_[idx], q)) ++c;
      return true;
    });
    return c;
  }

 private:
  const BoxT* boxes_{nullptr};
  usize n_{0};
  std::vector<PointT> pts_;
  RangeTree<2, Scalar> tree_{};
};

// --------------------------
// Factory
// --------------------------

BoxIndex2D MakeIndex2D(std::string_view name, std::string* err) {
  BoxIndex2D out;

  auto make = [&](auto* impl, std::string pretty_name) {
    using ImplT = std::decay_t<decltype(*impl)>;
    out.impl = impl;
    out.destroy = [](void* p) { delete static_cast<ImplT*>(p); };
    out.build = [](void* p, Span<const BoxIndex2D::BoxT> boxes) {
      static_cast<ImplT*>(p)->Build(boxes);
    };
    out.query = [](const void* p, const BoxIndex2D::BoxT& q, std::vector<usize>* out_vec) {
      static_cast<const ImplT*>(p)->Query(q, out_vec);
    };
    out.count = [](const void* p, const BoxIndex2D::BoxT& q) -> u64 {
      return static_cast<const ImplT*>(p)->Count(q);
    };
    out.name = std::move(pretty_name);
  };

  // Normalize a few common aliases.
  if (EqualsIgnoreCase(name, "interval") || EqualsIgnoreCase(name, "interval_tree") || EqualsIgnoreCase(name, "it")) {
    make(new IntervalTreeIndex2DImpl(), "interval_tree");
    return out;
  }
  if (EqualsIgnoreCase(name, "aabb") || EqualsIgnoreCase(name, "aabb_tree") || EqualsIgnoreCase(name, "bvh")) {
    make(new AABBTreeIndex2DImpl(), "aabb_tree");
    return out;
  }
  if (EqualsIgnoreCase(name, "r") || EqualsIgnoreCase(name, "rtree") || EqualsIgnoreCase(name, "r_tree")) {
    make(new RTreeIndex2DImpl(), "r_tree");
    return out;
  }
  if (EqualsIgnoreCase(name, "kd") || EqualsIgnoreCase(name, "kdtree") || EqualsIgnoreCase(name, "kd_tree")) {
    make(new KDTreeIndex2DImpl(), "kd_tree");
    return out;
  }
  if (EqualsIgnoreCase(name, "range") || EqualsIgnoreCase(name, "range_tree") || EqualsIgnoreCase(name, "rangetree")) {
    make(new RangeTreeIndex2DImpl(), "range_tree");
    return out;
  }

  SetErr(err,
         "Unknown index name '" + std::string(name) +
         "'. Supported: interval_tree, aabb_tree, r_tree, kd_tree, range_tree.");
  return out;
}

// For convenience in interactive tooling:
std::string IndexHelp2D() {
  return std::string(
      "Supported index names (Dim=2):\n"
      "  interval_tree  : 1D interval tree on y (filters to exact 2D)\n"
      "  aabb_tree      : static AABB tree (BVH) over 2D boxes\n"
      "  r_tree         : packed R-tree over 2D boxes\n"
      "  kd_tree        : KD-tree over embedding p=(L_y,R_y) (filters to exact 2D)\n"
      "  range_tree     : range tree over embedding p=(L_y,R_y) (filters to exact 2D)\n");
}

}  // namespace index
}  // namespace sjs
