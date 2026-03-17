
#include "core/rng.h"
#include "join/join_types.h"
#include "join/sweep_events.h"

#include <algorithm>
#include <cstdio>
#include <vector>

using namespace sjs;

static void FillRandomEvents(u32 n, Rng* rng, std::vector<join::Event>* out) {
  out->clear();
  out->reserve(n);

  for (u32 i = 0; i < n; ++i) {
    join::Event e;
    // Keep x finite and with both signs.
    const i64 raw = static_cast<i64>(rng->UniformU64(1ULL << 52));
    const double x = (rng->UniformU32(2) == 0) ? double(raw) : -double(raw);
    e.x = x * 1e-3;

    e.kind = (rng->UniformU32(2) == 0) ? join::EventKind::End : join::EventKind::Start;
    e.side = (rng->UniformU32(2) == 0) ? join::Side::R : join::Side::S;

    // Force some ties.
    e.id = static_cast<Id>(rng->UniformU32(1000));
    e.index = rng->UniformU32(1000);

    out->push_back(e);
  }
}

static bool CheckOrder(const std::vector<join::Event>& a, const std::vector<join::Event>& b) {
  if (a.size() != b.size()) return false;
  for (usize i = 0; i < a.size(); ++i) {
    const auto& x = a[i];
    const auto& y = b[i];
    if (!(x.x == y.x) || x.kind != y.kind || x.side != y.side || x.id != y.id || x.index != y.index) {
      return false;
    }
  }
  return true;
}

int main() {
  Rng rng(12345);

  for (u32 rep = 0; rep < 50; ++rep) {
    const u32 n = 5000 + rng.UniformU32(20000);
    std::vector<join::Event> ev;
    FillRandomEvents(n, &rng, &ev);

    for (auto order : {join::SideTieBreak::RBeforeS, join::SideTieBreak::SBeforeR}) {
      std::vector<join::Event> a = ev;
      std::vector<join::Event> b = ev;

      // Reference: comparator sort.
      join::EventLess less{order};
      std::sort(a.begin(), a.end(), less);

      // Under test: SortSweepEvents (may use radix for large sizes).
      join::SortSweepEvents(&b, order);

      if (!CheckOrder(a, b)) {
        std::fprintf(stderr, "Radix sort mismatch (rep=%u)\n", rep);
        return 1;
      }
    }
  }

  return 0;
}
