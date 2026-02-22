// tests/test_rank_space.cpp
//
// Validate RankSpace (α=(a,id) lexicographic rank domain) and
// mapping from open interval (ℓ,r) to rank range [L,R).
//
// This test is aligned with the current RankSpace<Value,Handle> API:
//   - BuildFromValues(a_by_handle)
//   - BuildFromKeys(std::vector<Key>)
//   - RankOf(id), UpperBound(x), LowerBound(x), OpenToRankRange(ℓ,r)

#include "sjs/core/types.h"
#include "sjs/index/rank_space.h"

#include <iostream>
#include <utility>
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
    std::cerr << "[FAIL] " << file << ":" << line << "  CHECK_EQ(" << ea << ", " << eb << ")  got " << a
              << " vs " << b << "\n";
  }
};

#define CHECK(ctx, expr) (ctx).Check((expr), #expr, __FILE__, __LINE__)
#define CHECK_EQ(ctx, a, b) (ctx).CheckEq((a), (b), #a, #b, __FILE__, __LINE__)

using sjs::Id;
using sjs::Scalar;
namespace idx = sjs::index;

static void TestBuildFromValues(TestContext& t) {
  using RS = idx::RankSpace<Scalar, Id>;

  // Values with duplicates.
  // With dense handles [0..n), BuildFromValues uses handle=i.
  // Sorted by α=(a,id) lexicographic should be:
  //   (0,id0),(0,id1),(1,id2),(2,id3),(2,id4)
  const std::vector<Scalar> a_by_handle = {0.0, 0.0, 1.0, 2.0, 2.0};

  RS rs;
  rs.BuildFromValues(a_by_handle);

  CHECK_EQ(t, rs.Size(), 5u);
  CHECK(t, !rs.Empty());

  // Rank order among duplicates must follow id.
  CHECK_EQ(t, rs.RankOf(0u), 0u);
  CHECK_EQ(t, rs.RankOf(1u), 1u);
  CHECK_EQ(t, rs.RankOf(2u), 2u);

  // upper_bound(0) => first a>0 => element with a=1.
  CHECK_EQ(t, rs.UpperBound(0.0), 2u);

  // lower_bound(2) => first a>=2 => first 2.0
  CHECK_EQ(t, rs.LowerBound(2.0), 3u);

  // Out-of-range should return N.
  CHECK_EQ(t, rs.LowerBound(100.0), rs.Size());
  CHECK_EQ(t, rs.UpperBound(100.0), rs.Size());

  // Open interval (0,2): expect only a=1.0 (id2)
  {
    const auto lr = rs.OpenToRankRange(0.0, 2.0);
    const sjs::u32 L = lr.first;
    const sjs::u32 R = lr.second;
    CHECK(t, L < R);
    CHECK_EQ(t, R - L, 1u);
    const sjs::u32 r2 = rs.RankOf(2u);
    CHECK(t, L <= r2 && r2 < R);

    // Boundary values excluded: id0 (a=0) and id3 (a=2) are outside (0,2)
    const sjs::u32 r0 = rs.RankOf(0u);
    const sjs::u32 r3 = rs.RankOf(3u);
    CHECK(t, !(L <= r0 && r0 < R));
    CHECK(t, !(L <= r3 && r3 < R));
  }

  // A few more strict-interval sanity checks.
  // (-1,0) excludes a==0 due to strict < r
  {
    const auto lr = rs.OpenToRankRange(-1.0, 0.0);
    CHECK_EQ(t, lr.first, 0u);
    CHECK_EQ(t, lr.second, 0u);
  }
  // (0,100) includes all a>0 => ids 2,3,4
  {
    const auto lr = rs.OpenToRankRange(0.0, 100.0);
    CHECK_EQ(t, lr.first, 2u);
    CHECK_EQ(t, lr.second, rs.Size());
    CHECK_EQ(t, lr.second - lr.first, 3u);
  }
  // (2,2) empty
  {
    const auto lr = rs.OpenToRankRange(2.0, 2.0);
    CHECK(t, lr.first >= lr.second);
  }
}

static void TestBuildFromKeys(TestContext& t) {
  using RS = idx::RankSpace<Scalar, Id>;
  RS rs;

  // Explicit key list with non-dense ids, to exercise BuildFromKeys.
  // keys are intentionally unsorted.
  std::vector<RS::Key> keys;
  keys.push_back(RS::Key{2.0, 10u});
  keys.push_back(RS::Key{0.0, 11u});
  keys.push_back(RS::Key{0.0, 7u});
  keys.push_back(RS::Key{1.0, 3u});

  rs.BuildFromKeys(std::move(keys));
  CHECK_EQ(t, rs.Size(), 4u);

  // Sorted order should be (0,7),(0,11),(1,3),(2,10).
  CHECK_EQ(t, rs.At(0).a, 0.0);
  CHECK_EQ(t, rs.At(0).id, 7u);
  CHECK_EQ(t, rs.At(1).a, 0.0);
  CHECK_EQ(t, rs.At(1).id, 11u);
  CHECK_EQ(t, rs.At(2).a, 1.0);
  CHECK_EQ(t, rs.At(2).id, 3u);
  CHECK_EQ(t, rs.At(3).a, 2.0);
  CHECK_EQ(t, rs.At(3).id, 10u);

  CHECK_EQ(t, rs.RankOf(7u), 0u);
  CHECK_EQ(t, rs.RankOf(11u), 1u);
  CHECK_EQ(t, rs.RankOf(3u), 2u);
  CHECK_EQ(t, rs.RankOf(10u), 3u);

  // Open interval (0,2) should contain only (1,3).
  {
    const auto lr = rs.OpenToRankRange(0.0, 2.0);
    CHECK_EQ(t, lr.second - lr.first, 1u);
    const sjs::u32 r3 = rs.RankOf(3u);
    CHECK(t, lr.first <= r3 && r3 < lr.second);
  }
}

}  // namespace

int main() {
  TestContext t;

  TestBuildFromValues(t);
  TestBuildFromKeys(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_rank_space\n";
    return 0;
  }
  std::cerr << "[FAILED] test_rank_space: " << t.fails << " failure(s)\n";
  return 1;
}
