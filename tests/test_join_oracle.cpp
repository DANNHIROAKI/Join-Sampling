// tests/test_join_oracle.cpp
//
// Correctness tests for join oracles and enumerators:
//  - Naive O(|R||S|) join count / enumeration
//  - Plane-sweep join enumerator
//
// Goal: ensure the library's basic join machinery is correct and consistent.

#include "sjs/core/types.h"
#include "sjs/io/dataset.h"
#include "sjs/join/join_enumerator.h"
#include "sjs/join/join_oracle.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_set>
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
  void CheckEq(const A& a, const B& b, const char* ea, const char* eb,
               const char* file, int line) {
    if (a == b) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line << "  CHECK_EQ(" << ea << ", " << eb
              << ")  got " << a << " vs " << b << "\n";
  }
};

#define CHECK(ctx, expr) (ctx).Check((expr), #expr, __FILE__, __LINE__)
#define CHECK_EQ(ctx, a, b) (ctx).CheckEq((a), (b), #a, #b, __FILE__, __LINE__)

using sjs::Box;
using sjs::Dataset;
using sjs::PairId;
using sjs::PairIdHash;
using sjs::Point;
using sjs::Scalar;

static bool PairLess(const PairId& a, const PairId& b) {
  if (a.r != b.r) return a.r < b.r;
  return a.s < b.s;
}

static Dataset<2, Scalar> MakeTinyDataset() {
  using P2 = Point<2, Scalar>;
  using B2 = Box<2, Scalar>;

  Dataset<2, Scalar> ds;
  ds.name = "tiny_oracle";

  // R: three rectangles
  ds.R.name = "R";
  ds.R.Add(B2(P2{0.0, 0.0}, P2{1.0, 1.0}));    // r0
  ds.R.Add(B2(P2{0.5, 0.5}, P2{1.5, 1.5}));    // r1
  ds.R.Add(B2(P2{2.0, 2.0}, P2{3.0, 3.0}));    // r2

  // S: three rectangles
  ds.S.name = "S";
  ds.S.Add(B2(P2{0.2, 0.2}, P2{0.8, 0.8}));     // s0 intersects r0,r1
  ds.S.Add(B2(P2{1.0, 1.0}, P2{2.5, 2.5}));     // s1 intersects r1,r2 (NOT r0 because boundary)
  ds.S.Add(B2(P2{3.0, 0.0}, P2{4.0, 1.0}));     // s2 disjoint

  ds.EnsureIds();
  return ds;
}

static void TestCountsMatch(TestContext& t) {
  const auto ds = MakeTinyDataset();
  std::string err;
  CHECK(t, ds.Validate(true, &err));

  sjs::join::JoinStats st_naive;
  const sjs::u64 c_naive = sjs::join::CountNaive<2, Scalar>(ds.R, ds.S, &st_naive);
  CHECK_EQ(t, c_naive, 4u);

  sjs::join::JoinStats st_sweep;
  sjs::join::PlaneSweepOptions opt;
  opt.axis = 0;
  opt.side_order = sjs::join::SideTieBreak::RBeforeS;
  opt.skip_axis_check = false;

  const sjs::u64 c_sweep = sjs::join::CountPlaneSweep<2, Scalar>(ds.R, ds.S, opt, &st_sweep);
  CHECK_EQ(t, c_sweep, c_naive);

  // Naive must check all pairs.
  CHECK_EQ(t, st_naive.candidate_checks, ds.R.Size() * ds.S.Size());
  // Sweep should not exceed naive.
  CHECK(t, st_sweep.candidate_checks <= st_naive.candidate_checks);
}

static void TestEnumerationMatches(TestContext& t) {
  const auto ds = MakeTinyDataset();

  // Collect oracle pairs.
  auto oracle_pairs = sjs::join::CollectNaivePairs<2, Scalar>(ds.R, ds.S, /*cap=*/0);
  std::sort(oracle_pairs.begin(), oracle_pairs.end(), PairLess);
  CHECK_EQ(t, oracle_pairs.size(), 4u);

  // Collect plane-sweep enumerator pairs.
  sjs::join::PlaneSweepOptions opt;
  opt.axis = 0;
  opt.side_order = sjs::join::SideTieBreak::RBeforeS;
  opt.skip_axis_check = true;

  sjs::join::PlaneSweepJoinStream<2, Scalar> stream(ds.R, ds.S, opt);

  std::vector<PairId> sweep_pairs;
  PairId p;
  while (stream.Next(&p)) {
    sweep_pairs.push_back(p);
  }

  std::sort(sweep_pairs.begin(), sweep_pairs.end(), PairLess);
  CHECK_EQ(t, sweep_pairs.size(), oracle_pairs.size());

  for (size_t i = 0; i < oracle_pairs.size(); ++i) {
    CHECK(t, oracle_pairs[i] == sweep_pairs[i]);
  }

  // Also verify there are no duplicates in the enumerated stream.
  std::unordered_set<PairId, PairIdHash> uniq;
  for (const auto& x : sweep_pairs) {
    CHECK(t, uniq.insert(x).second);
  }
}

static void TestBoundaryTouchIsNotIntersect(TestContext& t) {
  using P2 = Point<2, Scalar>;
  using B2 = Box<2, Scalar>;

  Dataset<2, Scalar> ds;
  ds.name = "touch_only";

  ds.R.Add(B2(P2{0.0, 0.0}, P2{1.0, 1.0}));
  ds.S.Add(B2(P2{1.0, 0.0}, P2{2.0, 1.0}));
  ds.EnsureIds();

  const sjs::u64 c = sjs::join::CountNaive<2, Scalar>(ds.R, ds.S);
  CHECK_EQ(t, c, 0u);

  sjs::join::PlaneSweepOptions opt;
  opt.axis = 0;
  opt.side_order = sjs::join::SideTieBreak::RBeforeS;
  opt.skip_axis_check = true;

  const sjs::u64 c2 = sjs::join::CountPlaneSweep<2, Scalar>(ds.R, ds.S, opt);
  CHECK_EQ(t, c2, 0u);
}

}  // namespace

int main() {
  TestContext t;

  TestCountsMatch(t);
  TestEnumerationMatches(t);
  TestBoundaryTouchIsNotIntersect(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_join_oracle\n";
    return 0;
  }
  std::cerr << "[FAILED] test_join_oracle: " << t.fails << " failure(s)\n";
  return 1;
}
