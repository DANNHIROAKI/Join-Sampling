// tests/test_rtree_split_regression.cpp
//
// Regression test: DynamicRTree must not keep dangling references across NewNode().
// This test forces many splits (max_children=2) and exercises insert/query/remove.
// The test must not crash and must keep correct counts.

#include "sjs/baselines/r_tree/sampling.h"
#include "sjs/core/types.h"
#include "sjs/geometry/box.h"

#include <iostream>

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
};

#define CHECK(ctx, expr) (ctx).Check((expr), #expr, __FILE__, __LINE__)
#define CHECK_EQ(ctx, a, b) (ctx).CheckEq((a), (b), #a, #b, __FILE__, __LINE__)

using Scalar = sjs::Scalar;
using sjs::Box;
using sjs::u32;
using sjs::u64;

void TestManySplitsNoCrash(TestContext& t) {
  using Tree = sjs::baselines::r_tree::detail::DynamicRTree<1, Scalar>;

  Tree tree;
  Tree::Options opt;
  opt.max_children = 2;   // force frequent splits
  opt.min_children = 1;
  opt.ignore_duplicate_insert = false;

  const u32 N = 2000;
  tree.Init(N, opt);

  // Insert N intervals: [i, i+0.5)
  for (u32 i = 0; i < N; ++i) {
    Box<1, Scalar> b;
    b.lo.v[0] = static_cast<Scalar>(i);
    b.hi.v[0] = static_cast<Scalar>(i) + static_cast<Scalar>(0.5);
    CHECK(t, tree.Insert(i, b));
  }
  CHECK_EQ(t, tree.Size(), N);

  // Query that intersects all inserted intervals
  Box<1, Scalar> q;
  q.lo.v[0] = static_cast<Scalar>(0);
  q.hi.v[0] = static_cast<Scalar>(N);
  const u64 cnt = tree.CountIntersect(q);
  CHECK_EQ(t, cnt, static_cast<u64>(N));

  // Remove all
  for (u32 i = 0; i < N; ++i) {
    CHECK(t, tree.Remove(i));
  }
  CHECK_EQ(t, tree.Size(), 0u);

  // Repeat a few rounds to exercise Clear() + reuse paths
  for (int round = 0; round < 3; ++round) {
    tree.Clear();
    const u32 M = 500;
    for (u32 i = 0; i < M; ++i) {
      Box<1, Scalar> b;
      b.lo.v[0] = static_cast<Scalar>(i);
      b.hi.v[0] = static_cast<Scalar>(i) + static_cast<Scalar>(1);
      CHECK(t, tree.Insert(i, b));
    }
    for (u32 i = 0; i < M; ++i) CHECK(t, tree.Remove(i));
    CHECK_EQ(t, tree.Size(), 0u);
  }
}

}  // namespace

int main() {
  TestContext t;
  TestManySplitsNoCrash(t);
  if (t.fails) return 1;
  std::cerr << "[OK] test_rtree_split_regression\n";
  return 0;
}
