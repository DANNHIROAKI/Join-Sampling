// tests/test_segtree_1d.cpp
//
// 1D segment-tree building blocks tests (aligned with current headers):
//  - Dynamic stabbing segment tree (SJS-HighDims §2.2): intervals stabbed by a point x
//      * SegTreeStabbingDyn<Value,Handle>::Init(shared_ptr<const vector<Value>>, OccPool*)
//      * InsertInterval(a,b,h), Count(x), Report(x,out), Sample(x,k,rng,out)
//      * Deletion is via OccPool::EraseAll(handle) (shared occurrence pool)
//  - Dynamic rank-range segment tree (SJS-HighDims §2.3): open interval (ℓ,r) mapped to rank range
//      * SegTreeRankRangeDyn<Value,Handle>::Init(shared_ptr<const RankSpace<Value,Handle>>, OccPool*)
//      * Insert(handle), Count(ℓ,r), Report(ℓ,r,out), Sample(ℓ,r,k,rng,out)
//      * Deletion is via OccPool::EraseAll(handle)

#include "sjs/core/rng.h"
#include "sjs/core/types.h"

#include "sjs/index/occ_pool.h"
#include "sjs/index/rank_space.h"
#include "sjs/index/segtree_rankrange_dyn.h"
#include "sjs/index/segtree_stabbing_dyn.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_set>
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

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

static std::vector<Id> SortedUnique(std::vector<Id> v) {
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return v;
}

struct Interval {
  Id id;
  Scalar lo;
  Scalar hi;
};

static bool ContainsHalfOpen1D(Scalar lo, Scalar hi, Scalar x) {
  return (lo <= x) && (x < hi);
}

static std::vector<Id> NaiveStab(const std::vector<Interval>& active, Scalar x) {
  std::vector<Id> out;
  for (const auto& it : active) {
    if (ContainsHalfOpen1D(it.lo, it.hi, x)) out.push_back(it.id);
  }
  return SortedUnique(std::move(out));
}

static std::vector<Id> NaiveOpenRange(const std::vector<std::pair<Id, Scalar>>& active, Scalar ell, Scalar r) {
  std::vector<Id> out;
  for (const auto& kv : active) {
    const Scalar a = kv.second;
    if (ell < a && a < r) out.push_back(kv.first);
  }
  return SortedUnique(std::move(out));
}

// -----------------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------------

static void TestStabbingSegTreeBasic(TestContext& t) {
  // Coordinate grid for segment tree: leaves correspond to
  // [0,1), [1,2), [2,3), [3,4)
  const std::vector<Scalar> coords = {0.0, 1.0, 2.0, 3.0, 4.0};
  auto coords_sp = std::make_shared<const std::vector<Scalar>>(coords);

  // Use dense handles [0..3] like the real codebase.
  idx::OccPool pool;
  pool.Init(/*num_handles=*/4u, /*reserve_occ=*/128u);

  using ST = idx::SegTreeStabbingDyn<Scalar, Id>;
  ST st;
  st.Init(coords_sp, &pool);
  CHECK(t, st.IsInitialized());

  std::vector<Interval> active;

  auto add = [&](Id id, Scalar lo, Scalar hi) {
    st.InsertInterval(lo, hi, id);
    active.push_back(Interval{id, lo, hi});
  };
  auto del = [&](Id id) {
    pool.EraseAll(id);
    active.erase(std::remove_if(active.begin(), active.end(), [&](const Interval& it) { return it.id == id; }),
                 active.end());
  };

  // Insert a few intervals (endpoints must appear in coords).
  add(0u, 0.0, 2.0);  // A
  add(1u, 1.0, 4.0);  // B
  add(2u, 2.0, 3.0);  // C
  add(3u, 0.0, 1.0);  // D

  const std::vector<Scalar> qs = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.5, 4.0};
  for (Scalar x : qs) {
    const auto ref = NaiveStab(active, x);

    std::vector<Id> raw;
    st.Report(x, &raw);
    const auto got = SortedUnique(raw);

    CHECK(t, got == ref);
    CHECK_EQ(t, st.Count(x), static_cast<sjs::u64>(ref.size()));

    // REPORT should not contain duplicates even before sorting/unique.
    CHECK_EQ(t, raw.size(), got.size());
  }

  // SAMPLE sanity: samples must be drawn from the reported set.
  {
    const Scalar x = 1.5;
    const auto ref = NaiveStab(active, x);
    std::unordered_set<Id> ref_set(ref.begin(), ref.end());

    sjs::Rng rng(123);
    std::vector<Id> samp;
    st.Sample(x, /*k=*/32, &rng, &samp);
    CHECK_EQ(t, samp.size(), 32u);
    for (Id id : samp) {
      CHECK(t, ref_set.find(id) != ref_set.end());
    }
  }

  // Delete one interval and re-check.
  del(1u);  // remove B
  for (Scalar x : qs) {
    const auto ref = NaiveStab(active, x);
    std::vector<Id> raw;
    st.Report(x, &raw);
    const auto got = SortedUnique(raw);
    CHECK(t, got == ref);
    CHECK_EQ(t, st.Count(x), static_cast<sjs::u64>(ref.size()));
    CHECK_EQ(t, raw.size(), got.size());
  }
}

static void TestRankRangeSegTreeBasic(TestContext& t) {
  // Universe for rank space (a(u)=L(u) in §2.3).
  // Values with duplicates:
  //   id0:0, id1:0, id2:1, id3:2, id4:3, id5:3
  const std::vector<Scalar> a_by_handle = {0.0, 0.0, 1.0, 2.0, 3.0, 3.0};

  using RS = idx::RankSpace<Scalar, Id>;
  auto rs = std::make_shared<RS>();
  rs->BuildFromValues(a_by_handle);

  idx::OccPool pool;
  pool.Init(/*num_handles=*/static_cast<sjs::u32>(a_by_handle.size()), /*reserve_occ=*/256u);

  using RT = idx::SegTreeRankRangeDyn<Scalar, Id>;
  RT rt;
  rt.Init(rs, &pool);
  CHECK(t, rt.IsInitialized());

  // Active subset: ids {0,2,4,5}
  std::vector<std::pair<Id, Scalar>> active;
  auto ins = [&](Id id) {
    rt.Insert(id);
    active.emplace_back(id, a_by_handle[static_cast<std::size_t>(id)]);
  };
  auto del = [&](Id id) {
    pool.EraseAll(id);
    active.erase(std::remove_if(active.begin(), active.end(), [&](const auto& kv) { return kv.first == id; }),
                 active.end());
  };

  ins(0u);
  ins(2u);
  ins(4u);
  ins(5u);

  struct Q {
    Scalar ell;
    Scalar r;
  };
  const std::vector<Q> queries = {
      {-1.0, 4.0},  // should include all active (0,2,4,5)
      {0.0, 3.0},   // only a in (0,3): id2
      {2.0, 4.0},   // only a in (2,4): id4,id5
      {2.0, 3.0},   // empty (strict): excludes 3
      {0.0, 0.1},   // empty (strict): excludes 0
  };

  for (const auto& q : queries) {
    const auto ref = NaiveOpenRange(active, q.ell, q.r);

    std::vector<Id> raw;
    rt.Report(q.ell, q.r, &raw);
    const auto got = SortedUnique(raw);

    CHECK(t, got == ref);
    CHECK_EQ(t, rt.Count(q.ell, q.r), static_cast<sjs::u64>(ref.size()));
    CHECK_EQ(t, raw.size(), got.size());
  }

  // SAMPLE sanity on a non-empty query.
  {
    const Scalar ell = -1.0;
    const Scalar r = 4.0;
    const auto ref = NaiveOpenRange(active, ell, r);
    std::unordered_set<Id> ref_set(ref.begin(), ref.end());

    sjs::Rng rng(999);
    std::vector<Id> samp;
    rt.Sample(ell, r, /*k=*/64, &rng, &samp);
    CHECK_EQ(t, samp.size(), 64u);
    for (Id id : samp) {
      CHECK(t, ref_set.find(id) != ref_set.end());
    }
  }

  // Delete an element and re-check.
  del(2u);  // remove id2 (a=1)
  for (const auto& q : queries) {
    const auto ref = NaiveOpenRange(active, q.ell, q.r);

    std::vector<Id> raw;
    rt.Report(q.ell, q.r, &raw);
    const auto got = SortedUnique(raw);

    CHECK(t, got == ref);
    CHECK_EQ(t, rt.Count(q.ell, q.r), static_cast<sjs::u64>(ref.size()));
    CHECK_EQ(t, raw.size(), got.size());
  }
}

}  // namespace

int main() {
  TestContext t;

  TestStabbingSegTreeBasic(t);
  TestRankRangeSegTreeBasic(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_segtree_1d\n";
    return 0;
  }
  std::cerr << "[FAILED] test_segtree_1d: " << t.fails << " failure(s)\n";
  return 1;
}
