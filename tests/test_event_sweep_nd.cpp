// tests/test_event_sweep_nd.cpp
//
// ND sweep-line event tests (Dim >= 2 typical):
//  - START/END events per box on chosen sweep axis
//  - Half-open tie-break: END before START at same coordinate
//  - Deterministic total order tie-break (id first, then side if ids tie)
//
// This generalizes the 2D test to Dim=3 (and can be extended further).

#include "sjs/core/types.h"
#include "sjs/io/dataset.h"
#include "sjs/join/sweep_events.h"

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
using sjs::Point;
using sjs::Scalar;
using sjs::join::EventKind;
using sjs::join::Side;
using sjs::join::SideTieBreak;

template <int Dim>
static Dataset<Dim, Scalar> MakeTinyDatasetForSweepND() {
  using P = Point<Dim, Scalar>;
  using B = Box<Dim, Scalar>;

  auto make_box = [](Scalar x0, Scalar x1) -> B {
    P lo = P::Filled(0.0);
    P hi = P::Filled(1.0);
    lo[0] = x0;
    hi[0] = x1;
    return B(lo, hi);
  };

  Dataset<Dim, Scalar> ds;
  ds.name = "tiny_sweep_nd";
  ds.R.name = "R";
  ds.S.name = "S";

  // Two R boxes starting at x=0.
  ds.R.Add(make_box(0.0, 2.0));  // r0 (id 0)
  ds.R.Add(make_box(0.0, 1.0));  // r1 (id 1) ends at x=1

  // Two S boxes: one starting at 0 and ending at 1 (id 0),
  // and one starting exactly at x=1 to test END-before-START (id 1).
  ds.S.Add(make_box(0.0, 1.0));  // s0 (id 0)
  ds.S.Add(make_box(1.0, 2.0));  // s1 (id 1)

  ds.EnsureIds();
  return ds;
}

template <int Dim>
static void TestBasicCounts(TestContext& t) {
  auto ds = MakeTinyDatasetForSweepND<Dim>();
  std::string err;
  CHECK(t, ds.Validate(true, &err));

  const auto events = sjs::join::BuildSweepEvents<Dim, Scalar>(ds, /*axis=*/0, SideTieBreak::RBeforeS);

  const size_t expected = (ds.R.Size() + ds.S.Size()) * 2;
  CHECK_EQ(t, events.size(), expected);

  size_t r_start = 0, r_end = 0, s_start = 0, s_end = 0;
  for (const auto& e : events) {
    if (e.side == Side::R) {
      if (e.kind == EventKind::Start) ++r_start;
      else ++r_end;
    } else {
      if (e.kind == EventKind::Start) ++s_start;
      else ++s_end;
    }
  }
  CHECK_EQ(t, r_start, ds.R.Size());
  CHECK_EQ(t, r_end, ds.R.Size());
  CHECK_EQ(t, s_start, ds.S.Size());
  CHECK_EQ(t, s_end, ds.S.Size());
}

template <int Dim>
static void TestOrderingRBeforeS(TestContext& t) {
  auto ds = MakeTinyDatasetForSweepND<Dim>();
  auto events = sjs::join::BuildSweepEvents<Dim, Scalar>(ds, /*axis=*/0, SideTieBreak::RBeforeS);
  sjs::join::SortSweepEvents(&events, SideTieBreak::RBeforeS);

  // At x=0: START for r0(id0), r1(id1), s0(id0).
  // id-first: id0 events come before id1; for equal id0, RBeforeS => R then S.
  CHECK_EQ(t, events[0].x, 0.0);
  CHECK_EQ(t, events[0].kind, EventKind::Start);
  CHECK_EQ(t, events[0].side, Side::R);
  CHECK_EQ(t, events[0].id, 0u);

  CHECK_EQ(t, events[1].x, 0.0);
  CHECK_EQ(t, events[1].kind, EventKind::Start);
  CHECK_EQ(t, events[1].side, Side::S);
  CHECK_EQ(t, events[1].id, 0u);

  CHECK_EQ(t, events[2].x, 0.0);
  CHECK_EQ(t, events[2].kind, EventKind::Start);
  CHECK_EQ(t, events[2].side, Side::R);
  CHECK_EQ(t, events[2].id, 1u);

  // At x=1: END for r1(id1), END for s0(id0), START for s1(id1).
  // END must come before START at x=1.
  size_t b = 0;
  while (b < events.size() && events[b].x < 1.0) ++b;
  CHECK(t, b + 2 < events.size());
  CHECK_EQ(t, events[b].x, 1.0);

  // id-first among END at x=1: s0(id0) END comes before r1(id1) END.
  CHECK_EQ(t, events[b].kind, EventKind::End);
  CHECK_EQ(t, events[b].side, Side::S);
  CHECK_EQ(t, events[b].id, 0u);

  CHECK_EQ(t, events[b + 1].x, 1.0);
  CHECK_EQ(t, events[b + 1].kind, EventKind::End);
  CHECK_EQ(t, events[b + 1].side, Side::R);
  CHECK_EQ(t, events[b + 1].id, 1u);

  CHECK_EQ(t, events[b + 2].x, 1.0);
  CHECK_EQ(t, events[b + 2].kind, EventKind::Start);
  CHECK_EQ(t, events[b + 2].side, Side::S);
  CHECK_EQ(t, events[b + 2].id, 1u);
}

template <int Dim>
static void TestOrderingSBeforeR(TestContext& t) {
  auto ds = MakeTinyDatasetForSweepND<Dim>();
  auto events = sjs::join::BuildSweepEvents<Dim, Scalar>(ds, /*axis=*/0, SideTieBreak::SBeforeR);
  sjs::join::SortSweepEvents(&events, SideTieBreak::SBeforeR);

  // At x=0, for equal id0 (r0 and s0), SBeforeR => S START first.
  CHECK_EQ(t, events[0].x, 0.0);
  CHECK_EQ(t, events[0].kind, EventKind::Start);
  CHECK_EQ(t, events[0].side, Side::S);
  CHECK_EQ(t, events[0].id, 0u);

  // At x=1, END-before-START must still hold.
  size_t b = 0;
  while (b < events.size() && events[b].x < 1.0) ++b;
  CHECK(t, b + 2 < events.size());
  CHECK_EQ(t, events[b].x, 1.0);
  CHECK_EQ(t, events[b].kind, EventKind::End);
  CHECK_EQ(t, events[b + 1].x, 1.0);
  CHECK_EQ(t, events[b + 1].kind, EventKind::End);
  CHECK_EQ(t, events[b + 2].x, 1.0);
  CHECK_EQ(t, events[b + 2].kind, EventKind::Start);
}

}  // namespace

int main() {
  TestContext t;

  // Use Dim=3 for ND sweep verification.
  TestBasicCounts<3>(t);
  TestOrderingRBeforeS<3>(t);
  TestOrderingSBeforeR<3>(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_event_sweep_nd\n";
    return 0;
  }
  std::cerr << "[FAILED] test_event_sweep_nd: " << t.fails << " failure(s)\n";
  return 1;
}
