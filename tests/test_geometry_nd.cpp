// tests/test_geometry_nd.cpp
//
// ND geometry sanity tests (Dim >= 2 typical):
//  - Point / Box basic operations
//  - Half-open semantics ([lo, hi))
//  - Predicates: contains/intersects + intersection volume
//  - Embedding helpers (box-intersection -> orthogonal range query)
//
// This is an extension of the 2D version to multiple dimensions.

#include "sjs/core/types.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/embedding.h"
#include "sjs/geometry/point.h"
#include "sjs/geometry/predicates.h"

#include <algorithm>
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

template <int Dim>
Point<Dim, Scalar> PFill(Scalar x) {
  return Point<Dim, Scalar>::Filled(x);
}

template <int Dim>
void TestPointOps(TestContext& t) {
  using P = Point<Dim, Scalar>;

  P a = PFill<Dim>(0.0);
  P b = PFill<Dim>(1.0);
  a[0] = 2.0;

  const P c = a + b;
  CHECK_NEAR(t, c[0], 3.0, 1e-12);
  CHECK_NEAR(t, c[1], 1.0, 1e-12);

  const P d = c - b;
  CHECK(t, d == a);

  const P e = b * 2.0;
  CHECK_NEAR(t, e[0], 2.0, 1e-12);
  CHECK_NEAR(t, e[Dim - 1], 2.0, 1e-12);

  CHECK_NEAR(t, e.NormInf(), 2.0, 1e-12);
}

template <int Dim>
void TestHalfOpenContainsIntersects(TestContext& t) {
  using P = Point<Dim, Scalar>;
  using B = Box<Dim, Scalar>;

  const B b(PFill<Dim>(0.0), PFill<Dim>(1.0));

  // Half-open: lo included, hi excluded.
  CHECK(t, sjs::ContainsHalfOpen(b, PFill<Dim>(0.0)));

  P almost = PFill<Dim>(0.999);
  CHECK(t, sjs::ContainsHalfOpen(b, almost));

  P on_hi = PFill<Dim>(0.0);
  on_hi[0] = 1.0;
  CHECK(t, !sjs::ContainsHalfOpen(b, on_hi));

  // Touching at boundary in axis 0 should NOT intersect.
  P clo = PFill<Dim>(0.0);
  P chi = PFill<Dim>(1.0);
  clo[0] = 1.0;
  chi[0] = 2.0;
  const B c(clo, chi);
  CHECK(t, !sjs::IntersectsHalfOpen(b, c));

  // Small overlap in axis 0 should intersect.
  P dlo = PFill<Dim>(0.0);
  P dhi = PFill<Dim>(1.0);
  dlo[0] = 0.999;
  dhi[0] = 2.0;
  const B d(dlo, dhi);
  CHECK(t, sjs::IntersectsHalfOpen(b, d));

  // Empty box should never intersect.
  P elo = PFill<Dim>(0.0);
  P ehi = PFill<Dim>(1.0);
  ehi[0] = elo[0];  // width 0 in axis 0
  const B e(elo, ehi);
  CHECK(t, e.IsEmpty());
  CHECK(t, !sjs::IntersectsHalfOpen(b, e));
}

template <int Dim>
void TestIntersectionVolume(TestContext& t) {
  using P = Point<Dim, Scalar>;
  using B = Box<Dim, Scalar>;

  // a = [0,2)^Dim, b = [1,3)^Dim -> intersection is [1,2)^Dim => volume 1.
  const B a(PFill<Dim>(0.0), PFill<Dim>(2.0));
  const B b(PFill<Dim>(1.0), PFill<Dim>(3.0));

  const double vol = static_cast<double>(sjs::IntersectionVolume(a, b));
  CHECK_NEAR(t, vol, 1.0, 1e-12);

  // Touching boundary: b2 = [2,3)^Dim -> intersection is empty (half-open), volume 0.
  const B b2(PFill<Dim>(2.0), PFill<Dim>(3.0));
  CHECK_NEAR(t, static_cast<double>(sjs::IntersectionVolume(a, b2)), 0.0, 1e-12);
}

template <int Dim>
void TestBoxExpand(TestContext& t) {
  using P = Point<Dim, Scalar>;
  using B = Box<Dim, Scalar>;

  B b = B::Empty();
  CHECK(t, b.IsEmpty());

  P p = PFill<Dim>(0.0);
  p[0] = 1.0;
  if constexpr (Dim >= 2) p[1] = 2.0;

  b.ExpandToIncludePoint(p);
  CHECK(t, !b.IsEmpty());
  CHECK(t, sjs::ContainsHalfOpen(b, p));

  P q = PFill<Dim>(0.0);
  q[0] = -1.0;
  if constexpr (Dim >= 2) q[1] = -2.0;

  b.ExpandToIncludePoint(q);
  CHECK(t, sjs::ContainsHalfOpen(b, p));
  CHECK(t, sjs::ContainsHalfOpen(b, q));
}

template <int Dim>
void TestEmbeddingIntersectRange(TestContext& t) {
  static_assert(Dim >= 2, "Embedding tests assume Dim >= 2");
  using P = Point<Dim, Scalar>;
  using B = Box<Dim, Scalar>;

  // Query box q = [0,1)^Dim
  const B q(PFill<Dim>(0.0), PFill<Dim>(1.0));

  // r1 intersects q
  P r1lo = PFill<Dim>(0.5);
  P r1hi = PFill<Dim>(1.5);
  const B r1(r1lo, r1hi);

  // r2 touches boundary at dim0 (x=1) => NOT intersect
  P r2lo = PFill<Dim>(0.0);
  P r2hi = PFill<Dim>(1.0);
  r2lo[0] = 1.0;
  r2hi[0] = 2.0;
  const B r2(r2lo, r2hi);

  // r3 overlaps in dim0 but touches boundary at dim1 (y=1) => NOT intersect
  P r3lo = PFill<Dim>(0.0);
  P r3hi = PFill<Dim>(1.0);
  r3lo[0] = 0.5;
  r3hi[0] = 1.5;
  r3lo[1] = 1.0;
  r3hi[1] = 2.0;
  const B r3(r3lo, r3hi);

  // Domain bounds for embedding.
  sjs::DomainBounds<Dim, Scalar> dom;
  dom.ExpandToInclude(q);
  dom.ExpandToInclude(r1);
  dom.ExpandToInclude(r2);
  dom.ExpandToInclude(r3);

  const auto p1 = sjs::EmbedLowerUpper<Dim, Scalar>(r1);
  const auto p2 = sjs::EmbedLowerUpper<Dim, Scalar>(r2);
  const auto p3 = sjs::EmbedLowerUpper<Dim, Scalar>(r3);

  const auto range = sjs::MakeIntersectQueryRange<Dim, Scalar>(q, dom);

  // Intersection => embedded point in range.
  CHECK(t, sjs::ContainsHalfOpen(range, p1));
  CHECK(t, !sjs::ContainsHalfOpen(range, p2));
  CHECK(t, !sjs::ContainsHalfOpen(range, p3));

  // Skip dim0 embedding: only enforces dims 1..Dim-1 (needs x-overlap externally).
  const auto p1s = sjs::EmbedLowerUpperSkipDim0<Dim, Scalar>(r1);
  const auto p3s = sjs::EmbedLowerUpperSkipDim0<Dim, Scalar>(r3);
  const auto ranges = sjs::MakeIntersectQueryRangeSkipDim0<Dim, Scalar>(q, dom);

  CHECK(t, sjs::ContainsHalfOpen(ranges, p1s));
  CHECK(t, !sjs::ContainsHalfOpen(ranges, p3s));
}

}  // namespace

int main() {
  TestContext t;

  // A few representative dimensions.
  TestPointOps<2>(t);
  TestHalfOpenContainsIntersects<2>(t);
  TestIntersectionVolume<2>(t);
  TestBoxExpand<2>(t);
  TestEmbeddingIntersectRange<2>(t);

  TestPointOps<3>(t);
  TestHalfOpenContainsIntersects<3>(t);
  TestIntersectionVolume<3>(t);
  TestBoxExpand<3>(t);
  TestEmbeddingIntersectRange<3>(t);

  TestPointOps<5>(t);
  TestHalfOpenContainsIntersects<5>(t);
  TestIntersectionVolume<5>(t);
  TestBoxExpand<5>(t);
  TestEmbeddingIntersectRange<5>(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_geometry_nd\n";
    return 0;
  }
  std::cerr << "[FAILED] test_geometry_nd: " << t.fails << " failure(s)\n";
  return 1;
}
