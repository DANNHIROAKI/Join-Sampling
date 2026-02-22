// tests/test_mode_index_nd.cpp
//
// Correctness tests for sjs::index::ModeIndexND (\mathcal{D}_{m,g}) in high dimensions.
//
// We validate, for a fixed Dim and each mode-mask g in {A,B}^{Dim-1}:
//  - COUNT_g(q) equals the naive |S_g(q)|
//  - REPORT_g(q) returns exactly S_g(q) (as ids)
//  - SAMPLE_g(q,k) always returns ids in S_g(q) (and roughly uniform when |S_g|>1)
//  - Dynamic updates (Insert + OccPool::EraseAll) preserve correctness
//
// NOTE: ModeIndexND indexes the projection dimensions 1..Dim-1 and assumes
// the scan dimension is axis 0 (as in SJS-HighDims.md). We therefore construct
// boxes so that scan-axis intersection is always guaranteed for all active
// objects at the query start.

#include "sjs/core/rng.h"
#include "sjs/core/types.h"
#include "sjs/geometry/box.h"

#include "sjs/index/mode_index.h"
#include "sjs/index/occ_pool.h"
#include "sjs/index/rank_space.h"
#include "sjs/index/segtree_common.h"

#include "sjs/sampling/sample_quality.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>
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
  void CheckEq(const A& a,
               const B& b,
               const char* ea,
               const char* eb,
               const char* file,
               int line) {
    if (a == b) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line << "  CHECK_EQ(" << ea << ", " << eb
              << ")  got " << a << " vs " << b << "\n";
  }
};

#define CHECK(ctx, expr) (ctx).Check((expr), #expr, __FILE__, __LINE__)
#define CHECK_EQ(ctx, a, b) (ctx).CheckEq((a), (b), #a, #b, __FILE__, __LINE__)

using sjs::Id;
using sjs::Rng;
using sjs::Scalar;
using sjs::u32;
using sjs::u64;
using sjs::usize;

namespace idx = sjs::index;

// --------------------------
// Naive definition of mode-partition S_g(q) (projected dims only)
// --------------------------

template <int Dim, class T>
bool InMode(u32 mask, const sjs::Box<Dim, T>& s, const sjs::Box<Dim, T>& q) {
  static_assert(Dim >= 2, "ModeIndexND test requires Dim>=2");
  // bits correspond to projected dimensions j=1..Dim-1 (bit j-1)
  for (int j = 1; j < Dim; ++j) {
    const usize k = static_cast<usize>(j);
    const bool isB = ((mask >> static_cast<u32>(j - 1)) & 1u) != 0u;
    if (!isB) {
      // A: L_s <= L_q < R_s
      if (!(s.lo.v[k] <= q.lo.v[k] && q.lo.v[k] < s.hi.v[k])) return false;
    } else {
      // B: L_q < L_s < R_q (open interval)
      if (!(q.lo.v[k] < s.lo.v[k] && s.lo.v[k] < q.hi.v[k])) return false;
    }
  }
  return true;
}

template <int Dim, class T>
std::vector<Id> NaiveModeSet(u32 mask,
                             const std::vector<sjs::Box<Dim, T>>& boxes,
                             const std::vector<Id>& active_ids,
                             const sjs::Box<Dim, T>& q) {
  std::vector<Id> out;
  out.reserve(active_ids.size());
  for (Id id : active_ids) {
    const auto& s = boxes[static_cast<usize>(id)];
    if (InMode<Dim>(mask, s, q)) out.push_back(id);
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// --------------------------
// Build the shared per-axis universes expected by ModeIndexND
// --------------------------

template <int Dim>
std::vector<idx::ModeAxisData<Scalar, Id>> BuildAxes(const std::vector<sjs::Box<Dim, Scalar>>& boxes) {
  static_assert(Dim >= 2, "BuildAxes requires Dim>=2");

  const int m = Dim - 1;
  const u32 N = static_cast<u32>(boxes.size());

  std::vector<idx::ModeAxisData<Scalar, Id>> axes;
  axes.resize(static_cast<usize>(m));

  for (int axis_offset = 0; axis_offset < m; ++axis_offset) {
    const int axis = axis_offset + 1;

    idx::ModeAxisData<Scalar, Id> ax;

    // ---- coords for stabbing (A): unique endpoints ----
    {
      std::vector<Scalar> coords;
      coords.reserve(static_cast<usize>(2) * boxes.size());
      for (const auto& b : boxes) {
        coords.push_back(b.lo[axis]);
        coords.push_back(b.hi[axis]);
      }
      std::sort(coords.begin(), coords.end());
      coords.erase(std::unique(coords.begin(), coords.end()), coords.end());

      ax.coords = std::make_shared<const std::vector<Scalar>>(std::move(coords));
      ax.stab_leaves = (ax.coords->size() >= 2) ? static_cast<u32>(ax.coords->size() - 1) : 0u;
      ax.stab_base = idx::NextPow2(ax.stab_leaves);
    }

    // ---- rank space for range (B): alpha(u)=L(u) with tie-break by id ----
    {
      std::vector<Scalar> lefts;
      lefts.resize(boxes.size());
      for (usize i = 0; i < boxes.size(); ++i) {
        lefts[i] = boxes[i].lo[axis];
      }
      auto rank = std::make_shared<idx::RankSpace<Scalar, Id>>();
      rank->BuildFromValues(lefts);
      ax.rank = std::move(rank);

      ax.rank_leaves = N;
      ax.rank_base = idx::NextPow2(ax.rank_leaves);
    }

    axes[static_cast<usize>(axis_offset)] = std::move(ax);
  }

  return axes;
}

// --------------------------
// Random box generator with guaranteed scan-axis stabbing
// --------------------------

template <int Dim>
sjs::Box<Dim, Scalar> RandomBoxWithSweepStab(Rng& rng, double sweep_x) {
  using BoxT = sjs::Box<Dim, Scalar>;
  BoxT b = BoxT::Empty();

  // dim0: ensure L0 <= sweep_x < R0
  {
    const double lo = rng.UniformDouble(0.0, sweep_x);
    const double hi = rng.UniformDouble(sweep_x + 1e-6, 1.0);
    b.lo.v[0] = lo;
    b.hi.v[0] = std::max(hi, lo + 1e-6);
  }

  // projection dims 1..Dim-1
  for (int i = 1; i < Dim; ++i) {
    const usize k = static_cast<usize>(i);
    const double lo = rng.UniformDouble(0.0, 0.9);
    const double w = rng.UniformDouble(0.05, 0.5);
    const double hi = std::min(1.0, lo + w);
    b.lo.v[k] = lo;
    b.hi.v[k] = std::max(hi, lo + 1e-6);
  }

  return b;
}

template <int Dim>
sjs::Box<Dim, Scalar> RandomQueryBox(Rng& rng, double sweep_x) {
  using BoxT = sjs::Box<Dim, Scalar>;
  BoxT q = BoxT::Empty();

  // query starts at sweep_x on scan axis
  q.lo.v[0] = sweep_x;
  q.hi.v[0] = std::min(1.0, sweep_x + rng.UniformDouble(0.05, 0.5));

  for (int i = 1; i < Dim; ++i) {
    const usize k = static_cast<usize>(i);
    const double lo = rng.UniformDouble(0.0, 0.9);
    const double w = rng.UniformDouble(0.05, 0.5);
    const double hi = std::min(1.0, lo + w);
    q.lo.v[k] = lo;
    q.hi.v[k] = std::max(hi, lo + 1e-6);
  }

  return q;
}

// --------------------------
// The actual test routine per mask
// --------------------------

template <int Dim>
void TestOneMask(TestContext& t, u32 mask) {
  using BoxT = sjs::Box<Dim, Scalar>;
  static_assert(Dim >= 2, "Dim must be >=2");

  Rng rng(12345 + static_cast<u64>(mask) * 17ULL);
  const double sweep_x = 0.5;

  const int N = 48;
  std::vector<BoxT> universe;
  universe.reserve(N);

  for (int i = 0; i < N; ++i) {
    universe.push_back(RandomBoxWithSweepStab<Dim>(rng, sweep_x));
  }

  std::vector<Id> ids;
  ids.reserve(N);
  for (int i = 0; i < N; ++i) ids.push_back(static_cast<Id>(i));

  // Active set: start with all.
  std::vector<Id> active = ids;

  // Build shared per-axis data and an occurrence pool.
  auto axes = BuildAxes<Dim>(universe);

  idx::OccPool pool;
  pool.Init(static_cast<u32>(N));

  // Build ModeIndexND for this mask.
  idx::ModeIndexND<Dim, Scalar, Id> mi;
  mi.Init(&universe, &axes, &pool, /*axis_offset=*/0, mask);

  // Insert all actives.
  for (Id id : active) {
    mi.Insert(id);
  }

  // Queries.
  const int Q = 20;
  std::vector<BoxT> queries;
  queries.reserve(Q);
  for (int i = 0; i < Q; ++i) queries.push_back(RandomQueryBox<Dim>(rng, sweep_x));

  // Validate Count/Report for each query.
  for (const auto& q : queries) {
    const auto naive = NaiveModeSet<Dim>(mask, universe, active, q);

    const u64 c = mi.Count(q);
    CHECK_EQ(t, c, static_cast<u64>(naive.size()));

    std::vector<Id> rep;
    mi.Report(q, &rep);
    std::sort(rep.begin(), rep.end());
    rep.erase(std::unique(rep.begin(), rep.end()), rep.end());

    CHECK_EQ(t, rep.size(), naive.size());
    CHECK(t, rep == naive);
  }

  // Dynamic updates: delete first half, re-check, then re-insert.
  std::vector<Id> removed;
  removed.reserve(static_cast<usize>(N / 2));
  for (int i = 0; i < N / 2; ++i) removed.push_back(static_cast<Id>(i));

  for (Id id : removed) {
    pool.EraseAll(id);
  }

  active.erase(std::remove_if(active.begin(), active.end(),
                              [&](Id x) { return x < static_cast<Id>(N / 2); }),
               active.end());

  for (const auto& q : queries) {
    const auto naive = NaiveModeSet<Dim>(mask, universe, active, q);

    const u64 c = mi.Count(q);
    CHECK_EQ(t, c, static_cast<u64>(naive.size()));

    std::vector<Id> rep;
    mi.Report(q, &rep);
    std::sort(rep.begin(), rep.end());
    rep.erase(std::unique(rep.begin(), rep.end()), rep.end());
    CHECK(t, rep == naive);
  }

  // Re-insert removed, check again.
  for (Id id : removed) {
    mi.Insert(id);
  }
  active = ids;

  for (const auto& q : queries) {
    const auto naive = NaiveModeSet<Dim>(mask, universe, active, q);
    const u64 c = mi.Count(q);
    CHECK_EQ(t, c, static_cast<u64>(naive.size()));
  }

  // Sampling sanity: pick a query with |S_g(q)| >= 2, if exists.
  const u64 K = 2000;
  bool did_sampling = false;

  for (const auto& q : queries) {
    const auto naive = NaiveModeSet<Dim>(mask, universe, active, q);
    if (naive.size() < 2) continue;

    std::unordered_map<Id, usize> pos;
    pos.reserve(naive.size());
    for (usize i = 0; i < naive.size(); ++i) pos[naive[i]] = i;

    std::vector<Id> samples;
    Rng rng_s(999 + static_cast<u64>(mask));
    mi.Sample(q, K, &rng_s, &samples);

    CHECK_EQ(t, samples.size(), static_cast<usize>(K));

    std::vector<u64> counts(naive.size(), 0);
    for (Id x : samples) {
      auto it = pos.find(x);
      CHECK(t, it != pos.end());  // must be in set
      if (it != pos.end()) counts[it->second]++;
    }

    // Loose goodness-of-fit: chi-square p-value should not be extremely tiny.
    // Under uniformity, P(p < 1e-6) = 1e-6, so this is robust.
    const auto chi2 = sjs::sampling::quality::ChiSquareUniform(counts);
    CHECK(t, std::isfinite(chi2.p_value));
    CHECK(t, chi2.p_value > 1e-6);

    did_sampling = true;
    break;
  }

  // For some masks on random data, it is possible no query has >=2 candidates.
  (void)did_sampling;
}

}  // namespace

int main() {
  TestContext t;

  // Choose a dimension with m=Dim-1 >= 2 to exercise recursion.
  constexpr int Dim = 4;
  constexpr int m = Dim - 1;
  static_assert(m == 3, "Expected Dim=4 here");

  const u32 num_masks = 1u << static_cast<u32>(m);
  for (u32 mask = 0; mask < num_masks; ++mask) {
    TestOneMask<Dim>(t, mask);
  }

  if (t.fails == 0) {
    std::cout << "[OK] test_mode_index_nd\n";
    return 0;
  }
  std::cerr << "[FAILED] test_mode_index_nd: " << t.fails << " failure(s)\n";
  return 1;
}
