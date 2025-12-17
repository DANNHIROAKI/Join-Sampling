// tests/test_geometry.cpp
//
// Geometry sanity tests:
//  - Point / Box basic operations
//  - Half-open semantics ([lo, hi))
//  - Predicates: contains/intersects
//  - Embedding helpers (join -> range query)

#include "sjs/core/types.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/embedding.h"
#include "sjs/geometry/point.h"
#include "sjs/geometry/predicates.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct TestContext {
  int fails = 0;

  void Check(bool ok, const char* expr, const char* file, int line) {
    if (ok) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line << "  CHECK(" << expr << ")\n";
  }

  template <class A, class B>
  void CheckEq(const A& a, const B& b, const char* ea, const char* eb, const char* file, int line) {
    if (a == b) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line << "  CHECK_EQ(" << ea << ", " << eb
              << ")  got " << a << " vs " << b << "\n";
  }

  void CheckNear(double a, double b, double rel_eps, const char* ea, const char* eb,
                 const char* file, int line) {
    const double diff = std::fabs(a - b);
    const double scale = std::max({1.0, std::fabs(a), std::fabs(b)});
    if (diff <= rel_eps * scale) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line
              << "  CHECK_NEAR(" << ea << ", " << eb << ")  got " << a << " vs " << b
              << "  diff=" << diff << "\n";
  }
};

#define CHECK(ctx, expr) (ctx).Check((expr), #expr, __FILE__, __LINE__)
#define CHECK_EQ(ctx, a, b) (ctx).CheckEq((a), (b), #a, #b, __FILE__, __LINE__)
#define CHECK_NEAR(ctx, a, b, eps) (ctx).CheckNear((a), (b), (eps), #a, #b, __FILE__, __LINE__)

using sjs::Box;
using sjs::Point;
using sjs::Scalar;

void TestHalfOpenContainsIntersects(TestContext& t) {
  using P2 = Point<2, Scalar>;
  using B2 = Box<2, Scalar>;

  const B2 b(P2{0.0, 0.0}, P2{1.0, 1.0});

  // Half-open: lo is included, hi is excluded.
  CHECK(t, sjs::ContainsHalfOpen(b, P2{0.0, 0.0}));
  CHECK(t, sjs::ContainsHalfOpen(b, P2{0.999, 0.999}));
  CHECK(t, !sjs::ContainsHalfOpen(b, P2{1.0, 0.0}));
  CHECK(t, !sjs::ContainsHalfOpen(b, P2{0.0, 1.0}));
  CHECK(t, !sjs::ContainsHalfOpen(b, P2{1.0, 1.0}));

  const B2 c(P2{1.0, 0.0}, P2{2.0, 1.0});
  // Touching at x=1 should NOT intersect under half-open.
  CHECK(t, !sjs::IntersectsHalfOpen(b, c));

  const B2 d(P2{0.999, 0.0}, P2{2.0, 1.0});
  CHECK(t, sjs::IntersectsHalfOpen(b, d));

  // Self intersection.
  CHECK(t, sjs::IntersectsHalfOpen(b, b));

  // Empty boxes should never intersect.
  const B2 e(P2{0.0, 0.0}, P2{0.0, 1.0});
  CHECK(t, e.IsEmpty());
  CHECK(t, !sjs::IntersectsHalfOpen(b, e));
}

void TestIntersectionVolume(TestContext& t) {
  using P2 = Point<2, Scalar>;
  using B2 = Box<2, Scalar>;

  const B2 a(P2{0.0, 0.0}, P2{2.0, 2.0});
  const B2 b(P2{1.0, 1.0}, P2{3.0, 4.0});

  const double vol = sjs::IntersectionVolume(a, b);
  // Intersection is [1,2) x [1,2) => area 1.
  CHECK_NEAR(t, vol, 1.0, 1e-12);

  const B2 c(P2{2.0, 0.0}, P2{3.0, 1.0});
  // Touching boundary => zero volume.
  CHECK_NEAR(t, sjs::IntersectionVolume(a, c), 0.0, 1e-12);
}

void TestBoxExpandAndBounds(TestContext& t) {
  using P2 = Point<2, Scalar>;
  using B2 = Box<2, Scalar>;

  B2 b = B2::Empty();
  CHECK(t, b.IsEmpty());

  b.ExpandToIncludePoint(P2{1.0, 2.0});
  // ExpandToIncludePoint uses NextUp for hi, so it becomes a proper half-open box.
  CHECK(t, !b.IsEmpty());
  CHECK(t, sjs::ContainsHalfOpen(b, P2{1.0, 2.0}));

  // Expanding again should keep containing previous points.
  b.ExpandToIncludePoint(P2{-1.0, -2.0});
  CHECK(t, sjs::ContainsHalfOpen(b, P2{1.0, 2.0}));
  CHECK(t, sjs::ContainsHalfOpen(b, P2{-1.0, -2.0}));
}

void TestEmbeddingIntersectRange(TestContext& t) {
  using P2 = Point<2, Scalar>;
  using B2 = Box<2, Scalar>;

  // Query box q.
  const B2 q(P2{0.0, 0.0}, P2{1.0, 1.0});

  // r1 intersects q.
  const B2 r1(P2{0.5, 0.5}, P2{1.5, 1.5});

  // r2 does NOT intersect q (touch at boundary x=1).
  const B2 r2(P2{1.0, 0.0}, P2{2.0, 1.0});

  // r3 does NOT intersect q (touch at boundary y=1) but DOES overlap in x.
  // This is useful for testing SkipDim0 embeddings, where x overlap is assumed
  // to have been enforced by a sweep.
  const B2 r3(P2{0.5, 1.0}, P2{1.5, 2.0});

  // Build embedding domain bounds.
  sjs::DomainBounds<2, Scalar> dom;
  dom.ExpandToInclude(q);
  dom.ExpandToInclude(r1);
  dom.ExpandToInclude(r2);
  dom.ExpandToInclude(r3);

  const auto p1 = sjs::EmbedLowerUpper<2, Scalar>(r1);
  const auto p2 = sjs::EmbedLowerUpper<2, Scalar>(r2);
  const auto p3 = sjs::EmbedLowerUpper<2, Scalar>(r3);

  const auto range = sjs::MakeIntersectQueryRange<2, Scalar>(q, dom);

  // For this embedding, intersection implies the embedded point lies in range.
  CHECK(t, sjs::ContainsHalfOpen(range, p1));
  CHECK(t, !sjs::ContainsHalfOpen(range, p2));
  CHECK(t, !sjs::ContainsHalfOpen(range, p3));

  // Skip-dim0 embedding (drop x): this embedding only enforces dims 1..Dim-1.
  // It must be used together with an x-sweep (or other pre-filter) that ensures x-overlap.
  // Therefore we use r3 (which overlaps in x but fails in y) to test rejection.
  const auto p1s = sjs::EmbedLowerUpperSkipDim0<2, Scalar>(r1);
  const auto p3s = sjs::EmbedLowerUpperSkipDim0<2, Scalar>(r3);
  const auto ranges = sjs::MakeIntersectQueryRangeSkipDim0<2, Scalar>(q, dom);

  CHECK(t, sjs::ContainsHalfOpen(ranges, p1s));
  CHECK(t, !sjs::ContainsHalfOpen(ranges, p3s));
}

}  // namespace

int main() {
  TestContext t;

  TestHalfOpenContainsIntersects(t);
  TestIntersectionVolume(t);
  TestBoxExpandAndBounds(t);
  TestEmbeddingIntersectRange(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_geometry\n";
    return 0;
  }
  std::cerr << "[FAILED] test_geometry: " << t.fails << " failure(s)\n";
  return 1;
}
