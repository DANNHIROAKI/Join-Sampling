
#include "baselines/ours/sampling.h"
#include "core/rng.h"

#include <cstdio>
#include <vector>

using namespace sjs;

// A small, deterministic stress test for ActiveIndex2D (Pattern A/B).
// Verifies COUNT_A/COUNT_B against a naive reference and checks that
// SAMPLE_A/SAMPLE_B return only valid handles.

int main() {
  using namespace sjs::baselines::ours::detail;

  static constexpr u32 kM = 64;
  static constexpr u32 kN = 200;
  static constexpr u32 kOps = 5000;

  // Pre-generate rectangle y ranges.
  struct Rect {
    u32 ylo;
    u32 yhi;
  };
  std::vector<Rect> rects(kN);

  Rng rng(777);
  for (u32 h = 0; h < kN; ++h) {
    const u32 a = rng.UniformU32(kM);
    const u32 b = rng.UniformU32(kM + 1);
    const u32 lo = std::min(a, b);
    const u32 hi = std::max(a + 1U, b);
    rects[h] = Rect{lo, std::min<u32>(hi, kM)};
  }

  ActiveIndex2D idx;
  idx.Init(kN, kM);

  std::vector<u8> active(kN, 0);

  auto naive_count_a = [&](u32 q) -> u64 {
    u64 c = 0;
    for (u32 h = 0; h < kN; ++h) {
      if (!active[h]) continue;
      const auto r = rects[h];
      if (r.ylo <= q && q < r.yhi) ++c;
    }
    return c;
  };

  auto naive_count_b = [&](u32 qlo, u32 qhi) -> u64 {
    u64 c = 0;
    for (u32 h = 0; h < kN; ++h) {
      if (!active[h]) continue;
      const auto r = rects[h];
      // Pattern B: qlo < ylo < qhi.
      if (qlo < r.ylo && r.ylo < qhi) ++c;
    }
    return c;
  };

  std::vector<u32> buf;

  for (u32 op = 0; op < kOps; ++op) {
    const bool do_insert = (rng.UniformU32(2) == 0);

    if (do_insert) {
      // Insert a random inactive handle.
      u32 h = rng.UniformU32(kN);
      for (u32 tries = 0; tries < kN && active[h]; ++tries) h = (h + 1U) % kN;
      if (!active[h]) {
        active[h] = 1;
        idx.Insert(h, rects[h].ylo, rects[h].yhi);
      }
    } else {
      // Erase a random active handle.
      u32 h = rng.UniformU32(kN);
      for (u32 tries = 0; tries < kN && !active[h]; ++tries) h = (h + 1U) % kN;
      if (active[h]) {
        active[h] = 0;
        idx.Erase(h);
      }
    }

    // Validate counts for a few random queries.
    for (u32 qrep = 0; qrep < 8; ++qrep) {
      const u32 q = rng.UniformU32(kM);
      const u64 ca = idx.CountA(q);
      const u64 ca_ref = naive_count_a(q);
      if (ca != ca_ref) {
        std::fprintf(stderr, "COUNT_A mismatch: got=%llu ref=%llu\n",
                     (unsigned long long)ca, (unsigned long long)ca_ref);
        return 1;
      }

      const u32 qlo = rng.UniformU32(kM);
      const u32 qhi = rng.UniformU32(kM + 1);
      const u32 lo = std::min(qlo, qhi);
      const u32 hi = std::max(qlo, qhi);

      const u64 cb = idx.CountB(lo, hi);
      const u64 cb_ref = naive_count_b(lo, hi);
      if (cb != cb_ref) {
        std::fprintf(stderr, "COUNT_B mismatch: got=%llu ref=%llu (lo=%u hi=%u)\n",
                     (unsigned long long)cb, (unsigned long long)cb_ref, lo, hi);
        return 1;
      }

      // Validate sampling returns only valid handles.
      if (ca_ref > 0) {
        const u32 k = 32;
        if (!idx.SampleA(q, k, &rng, &buf) || buf.size() != k) {
          std::fprintf(stderr, "SAMPLE_A failed\n");
          return 1;
        }
        for (u32 x : buf) {
          if (x >= kN || !active[x]) return 1;
          const auto r = rects[x];
          if (!(r.ylo <= q && q < r.yhi)) return 1;
        }
      }

      if (cb_ref > 0) {
        const u32 k = 32;
        if (!idx.SampleB(lo, hi, k, &rng, &buf) || buf.size() != k) {
          std::fprintf(stderr, "SAMPLE_B failed\n");
          return 1;
        }
        for (u32 x : buf) {
          if (x >= kN || !active[x]) return 1;
          const auto r = rects[x];
          if (!(lo < r.ylo && r.ylo < hi)) return 1;
        }
      }
    }
  }

  return 0;
}
