// tests/test_mode_family_nd.cpp
//
// Correctness tests for sjs::index::ModeFamily (across all modes).
//
// We validate for a fixed Dim:
//  - COUNT(q) equals naive |{s in active : s^\downarrow intersects q^\downarrow }|
//  - REPORT(q) returns exactly that set (ids)
//  - SAMPLE(q,k) returns ids within that set, and looks roughly uniform when size>1
//  - Dynamic updates (Insert/Delete) preserve correctness
//
// The test constructs boxes so scan-axis (dim0) intersection at query-start is guaranteed,
// and checks only projection dims 1..Dim-1 for the partner set definition.

#include "sjs/core/rng.h"
#include "sjs/core/types.h"
#include "sjs/geometry/box.h"
#include "sjs/index/mode_family.h"
#include "sjs/sampling/sample_quality.h"

#include <algorithm>
#include <cmath>
#include <iostream>
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

using sjs::Id;
using sjs::Rng;
using sjs::Scalar;
using sjs::u64;
using sjs::usize;

namespace idx = sjs::index;

template <int Dim, class T>
bool ProjIntersects(const sjs::Box<Dim, T>& a, const sjs::Box<Dim, T>& b) {
  for (int i = 1; i < Dim; ++i) {
    const usize k = static_cast<usize>(i);
    if (!(a.lo.v[k] < b.hi.v[k] && b.lo.v[k] < a.hi.v[k])) return false;
  }
  return true;
}

template <int Dim, class T>
std::vector<Id> NaivePartnerSet(const std::vector<sjs::Box<Dim, T>>& boxes,
                                const std::vector<Id>& active_ids,
                                const sjs::Box<Dim, T>& q) {
  std::vector<Id> out;
  for (Id id : active_ids) {
    const auto& s = boxes[static_cast<usize>(id)];
    if (ProjIntersects<Dim>(s, q)) out.push_back(id);
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// Random box generator with guaranteed scan-axis stabbing (dim0 contains sweep_x).
template <int Dim>
sjs::Box<Dim, Scalar> RandomBoxWithSweepStab(Rng& rng, double sweep_x) {
  using BoxT = sjs::Box<Dim, Scalar>;
  BoxT b = BoxT::Empty();

  // dim0 contains sweep_x
  b.lo.v[0] = rng.UniformDouble(0.0, sweep_x);
  b.hi.v[0] = rng.UniformDouble(sweep_x + 1e-6, 1.0);
  if (!(b.hi.v[0] > b.lo.v[0])) b.hi.v[0] = b.lo.v[0] + 1e-6;

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

}  // namespace

int main() {
  TestContext t;

  constexpr int Dim = 4;
  using BoxT = sjs::Box<Dim, Scalar>;

  Rng rng(20240218);
  const double sweep_x = 0.5;

  const int N = 60;
  std::vector<BoxT> universe;
  std::vector<Id> ids;
  universe.reserve(N);
  ids.reserve(N);

  for (int i = 0; i < N; ++i) {
    universe.push_back(RandomBoxWithSweepStab<Dim>(rng, sweep_x));
    ids.push_back(static_cast<Id>(i));
  }

  std::vector<Id> active = ids;

  idx::ModeFamily<Dim, Scalar> fam;
  fam.Init(&universe);

  for (Id id : active) {
    fam.Insert(id);
  }

  const int Q = 25;
  std::vector<BoxT> queries;
  queries.reserve(Q);
  for (int i = 0; i < Q; ++i) queries.push_back(RandomQueryBox<Dim>(rng, sweep_x));

  // Count/Report checks.
  for (const auto& q : queries) {
    const auto naive = NaivePartnerSet<Dim>(universe, active, q);

    const u64 c = fam.Count(q);
    CHECK_EQ(t, c, static_cast<u64>(naive.size()));

    std::vector<Id> rep;
    fam.Report(q, &rep);
    std::sort(rep.begin(), rep.end());
    rep.erase(std::unique(rep.begin(), rep.end()), rep.end());
    CHECK(t, rep == naive);
  }

  // Dynamic updates: delete a chunk, re-check, then re-insert.
  std::vector<Id> removed;
  for (int i = 0; i < N / 3; ++i) removed.push_back(static_cast<Id>(i));

  for (Id id : removed) {
    fam.Delete(id);
  }
  active.erase(std::remove_if(active.begin(), active.end(),
                              [&](Id x) { return x < static_cast<Id>(N / 3); }),
               active.end());

  for (const auto& q : queries) {
    const auto naive = NaivePartnerSet<Dim>(universe, active, q);
    CHECK_EQ(t, fam.Count(q), static_cast<u64>(naive.size()));

    std::vector<Id> rep;
    fam.Report(q, &rep);
    std::sort(rep.begin(), rep.end());
    rep.erase(std::unique(rep.begin(), rep.end()), rep.end());
    CHECK(t, rep == naive);
  }

  for (Id id : removed) {
    fam.Insert(id);
  }
  active = ids;

  // Sampling sanity: choose a query with at least 2 candidates.
  const u64 K = 3000;
  bool did_sampling = false;
  for (const auto& q : queries) {
    const auto naive = NaivePartnerSet<Dim>(universe, active, q);
    if (naive.size() < 2) continue;

    std::unordered_map<Id, usize> pos;
    for (usize i = 0; i < naive.size(); ++i) pos[naive[i]] = i;

    std::vector<Id> samples;
    Rng rs(777);
    const bool ok = fam.Sample(q, K, &rs, &samples);
    CHECK(t, ok);
    CHECK_EQ(t, samples.size(), static_cast<usize>(K));

    std::vector<u64> counts(naive.size(), 0);
    for (Id x : samples) {
      auto it = pos.find(x);
      CHECK(t, it != pos.end());
      if (it != pos.end()) counts[it->second]++;
    }

    const auto chi2 = sjs::sampling::quality::ChiSquareUniform(counts);
    CHECK(t, std::isfinite(chi2.p_value));
    CHECK(t, chi2.p_value > 1e-6);

    did_sampling = true;
    break;
  }
  CHECK(t, did_sampling);

  if (t.fails == 0) {
    std::cout << "[OK] test_mode_family_nd\n";
    return 0;
  }
  std::cerr << "[FAILED] test_mode_family_nd: " << t.fails << " failure(s)\n";
  return 1;
}
