// tests/test_event_sweep.cpp
//
// Validate sweep-line event construction and sorting:
//  - START/END events per rectangle
//  - Half-open tie-break: END before START at the same coordinate
//  - Fixed total-order tie-break for identical x/kind (SJS v3 §1.3.1):
//      id first, then side tie-break only if ids tie
//
// These invariants are crucial because several baselines rely on a deterministic
// event ordering to make the join partitioning correct.

#include "core/types.h"
#include "io/dataset.h"
#include "join/sweep_events.h"

#include <algorithm>
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
using sjs::join::Event;
using sjs::join::EventKind;
using sjs::join::Side;
using sjs::join::SideTieBreak;

static Dataset<2, Scalar> MakeTinyDatasetForSweep() {
  using P2 = Point<2, Scalar>;
  using B2 = Box<2, Scalar>;

  Dataset<2, Scalar> ds;
  ds.name = "tiny_sweep";

  // Two R boxes starting at x=0.
  ds.R.name = "R";
  ds.R.Add(B2(P2{0.0, 0.0}, P2{2.0, 1.0}));  // r0
  ds.R.Add(B2(P2{0.0, 0.0}, P2{1.0, 1.0}));  // r1 ends at x=1

  // One S box also starting at x=0, ending at x=1.
  ds.S.name = "S";
  ds.S.Add(B2(P2{0.0, 0.0}, P2{1.0, 1.0}));  // s0 ends at x=1

  // Another S box starting exactly at x=1 (boundary), to test END-before-START.
  ds.S.Add(B2(P2{1.0, 0.0}, P2{2.0, 1.0}));  // s1 starts at x=1

  ds.EnsureIds();
  return ds;
}

static void TestBasicCounts(TestContext& t) {
  const auto ds = MakeTinyDatasetForSweep();
  std::string err;
  CHECK(t, ds.Validate(true, &err));

  const auto events = sjs::join::BuildSweepEvents<2, Scalar>(ds, /*axis=*/0, SideTieBreak::RBeforeS);

  // Each rectangle -> 2 events.
  const size_t expected = (ds.R.Size() + ds.S.Size()) * 2;
  CHECK_EQ(t, events.size(), expected);

  // Check each object appears exactly once with START and once with END.
  // (We only sanity-check counts; event ids are per-side.)
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

static void TestOrderingRBeforeS(TestContext& t) {
  const auto ds = MakeTinyDatasetForSweep();
  auto events = sjs::join::BuildSweepEvents<2, Scalar>(ds, /*axis=*/0, SideTieBreak::RBeforeS);
  sjs::join::SortSweepEvents(&events, SideTieBreak::RBeforeS);

  // At x=0: START events exist for r0(id0), r1(id1), s0(id0).
  // SJS v3 tie-break: id first; for equal ids, SideTieBreak applies.
  CHECK_EQ(t, events[0].x, 0.0);
  CHECK_EQ(t, events[0].kind, EventKind::Start);
  CHECK_EQ(t, events[0].side, Side::R);
  CHECK_EQ(t, events[0].id, 0u);

  // id0 on S follows id0 on R.
  CHECK_EQ(t, events[1].x, 0.0);
  CHECK_EQ(t, events[1].kind, EventKind::Start);
  CHECK_EQ(t, events[1].side, Side::S);
  CHECK_EQ(t, events[1].id, 0u);

  // id1 (r1) comes after all id0 starts.
  CHECK_EQ(t, events[2].x, 0.0);
  CHECK_EQ(t, events[2].kind, EventKind::Start);
  CHECK_EQ(t, events[2].side, Side::R);
  CHECK_EQ(t, events[2].id, 1u);

  // At x=1: END events for r1 and s0, and START for s1.
  // Half-open invariant: END must come before START at the same x.
  // Also, among END at same x, side tie-break applies (R before S).

  // Find the block of events with x==1.
  size_t b = 0;
  while (b < events.size() && events[b].x < 1.0) ++b;
  CHECK(t, b < events.size());
  CHECK_EQ(t, events[b].x, 1.0);

  // Under id-first order, S END of s0(id0) comes before R END of r1(id1).
  CHECK_EQ(t, events[b].kind, EventKind::End);
  CHECK_EQ(t, events[b].side, Side::S);
  CHECK_EQ(t, events[b].id, 0u);

  CHECK_EQ(t, events[b + 1].x, 1.0);
  CHECK_EQ(t, events[b + 1].kind, EventKind::End);
  CHECK_EQ(t, events[b + 1].side, Side::R);
  CHECK_EQ(t, events[b + 1].id, 1u);

  // Third at x=1 should be S START of s1.
  CHECK_EQ(t, events[b + 2].x, 1.0);
  CHECK_EQ(t, events[b + 2].kind, EventKind::Start);
  CHECK_EQ(t, events[b + 2].side, Side::S);
  CHECK_EQ(t, events[b + 2].id, 1u);
}

static void TestOrderingSBeforeR(TestContext& t) {
  const auto ds = MakeTinyDatasetForSweep();
  auto events = sjs::join::BuildSweepEvents<2, Scalar>(ds, /*axis=*/0, SideTieBreak::SBeforeR);
  sjs::join::SortSweepEvents(&events, SideTieBreak::SBeforeR);

  // At x=0, for equal id0 (r0 and s0), SBeforeR makes S START come before R START.
  CHECK_EQ(t, events[0].x, 0.0);
  CHECK_EQ(t, events[0].kind, EventKind::Start);
  CHECK_EQ(t, events[0].side, Side::S);

  // At x=1, END-before-START still must hold regardless of side ordering.
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

  TestBasicCounts(t);
  TestOrderingRBeforeS(t);
  TestOrderingSBeforeR(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_event_sweep\n";
    return 0;
  }
  std::cerr << "[FAILED] test_event_sweep: " << t.fails << " failure(s)\n";
  return 1;
}
