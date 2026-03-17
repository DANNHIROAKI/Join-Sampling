// third_party/RS-over-SRJ/tests/test_kdtree_4d.cpp
//
// Self-test for the generic KDTree<D,T> used by the RS-over-SRJ baseline.
//
// This test validates:
//  - Count(qlo,qhi) matches brute force
//  - Search(qlo,qhi) decomposes the exact query set without duplicates
//  - Sampling from the Search() decomposition never returns out-of-range points

#include "utils/kdtree.hpp"
#include "utils/weighted_sampling.hpp"

#include "core/rng.h"
#include "core/types.h"

#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

template <int D, class T>
bool InRange(const KDPoint<D, T>& p, const std::array<T, D>& qlo, const std::array<T, D>& qhi) {
  for (int d = 0; d < D; ++d) {
    const auto j = static_cast<std::size_t>(d);
    if (p.coord[j] < qlo[j] || p.coord[j] > qhi[j]) return false;
  }
  return true;
}

}  // namespace

int main() {
  using T = double;
  static constexpr int D = 4;

  sjs::Rng rng(1234567);

  // Build a random 4D point set.
  const sjs::usize n = 512;
  std::vector<KDPoint<D, T>> pts;
  pts.reserve(n);
  for (sjs::usize i = 0; i < n; ++i) {
    KDPoint<D, T> p;
    for (int d = 0; d < D; ++d) {
      p.coord[static_cast<std::size_t>(d)] = rng.UniformDouble(-10.0, 10.0);
    }
    p.idx = static_cast<std::uint32_t>(i);
    pts.push_back(p);
  }

  KDTree<D, T> tree;
  tree.Build(std::move(pts));

  if (tree.Size() != static_cast<std::uint32_t>(n)) {
    std::cerr << "FAIL: tree.Size mismatch\n";
    return 1;
  }

  // Map payload idx -> point coordinates (idx is unique 0..n-1).
  std::vector<KDPoint<D, T>> by_id(n);
  for (const auto& p : tree.Points()) {
    const sjs::usize id = static_cast<sjs::usize>(p.idx);
    if (id >= n) {
      std::cerr << "FAIL: payload idx out of range\n";
      return 1;
    }
    by_id[id] = p;
  }

  std::vector<typename KDTree<D, T>::SearchItem> items;
  std::vector<sjs::u8> mark;

  // Random queries.
  const int Q = 300;
  for (int qi = 0; qi < Q; ++qi) {
    std::array<T, D> qlo;
    std::array<T, D> qhi;
    for (int d = 0; d < D; ++d) {
      const double a = rng.UniformDouble(-10.0, 10.0);
      const double b = rng.UniformDouble(-10.0, 10.0);
      qlo[static_cast<std::size_t>(d)] = (a < b) ? a : b;
      qhi[static_cast<std::size_t>(d)] = (a < b) ? b : a;
    }

    // Brute-force count over the KDTree's permuted storage (still exact).
    sjs::u64 brute = 0;
    for (const auto& p : tree.Points()) {
      if (InRange<D, T>(p, qlo, qhi)) ++brute;
    }

    const sjs::u64 cnt = tree.Count(qlo, qhi);
    if (cnt != brute) {
      std::cerr << "FAIL: Count mismatch at query " << qi << ": got " << cnt << ", expected " << brute << "\n";
      return 1;
    }

    tree.Search(qlo, qhi, &items);

    // Expand the decomposition into a set of point positions [0, n).
    mark.assign(n, 0);
    sjs::u64 covered = 0;

    for (const auto& it : items) {
      const auto& nd = tree.GetNode(it.node);
      if (it.fully_contained) {
        for (std::uint32_t pos = nd.l; pos < nd.r; ++pos) {
          if (mark[static_cast<sjs::usize>(pos)] != 0) {
            std::cerr << "FAIL: duplicate point coverage (fully_contained)\n";
            return 1;
          }
          mark[static_cast<sjs::usize>(pos)] = 1;
          ++covered;
        }
      } else {
        const std::uint32_t pos = nd.mid;
        if (mark[static_cast<sjs::usize>(pos)] != 0) {
          std::cerr << "FAIL: duplicate point coverage (mid point)\n";
          return 1;
        }
        mark[static_cast<sjs::usize>(pos)] = 1;
        ++covered;
      }
    }

    // Verify every covered point is in range.
    for (sjs::usize pos = 0; pos < n; ++pos) {
      if (mark[pos] == 0) continue;
      if (!InRange<D, T>(tree.Points()[pos], qlo, qhi)) {
        std::cerr << "FAIL: Search returned an out-of-range point at pos=" << pos << "\n";
        return 1;
      }
    }

    if (covered != brute) {
      std::cerr << "FAIL: Search coverage mismatch: covered=" << covered << ", brute=" << brute << "\n";
      return 1;
    }

    // Sampling sanity: sample many times and ensure in-range.
    if (brute > 0) {
      std::vector<sjs::u64> weights;
      weights.reserve(items.size());
      for (const auto& it : items) {
        if (it.fully_contained) {
          const auto& nd = tree.GetNode(it.node);
          weights.push_back(static_cast<sjs::u64>(nd.Size()));
        } else {
          weights.push_back(1ULL);
        }
      }

      alias a;
      std::string err;
      if (!a.Build(weights, &err)) {
        std::cerr << "FAIL: alias build failed: " << err << "\n";
        return 1;
      }

      for (int k = 0; k < 2000; ++k) {
        const sjs::usize ii = a.Sample(&rng);
        if (ii >= items.size()) {
          std::cerr << "FAIL: alias sampled out of range\n";
          return 1;
        }

        const auto& it = items[ii];
        std::uint32_t s_idx = 0;
        if (it.fully_contained) {
          const auto& nd = tree.GetNode(it.node);
          const sjs::u64 off = rng.UniformU64(static_cast<sjs::u64>(nd.Size()));
          const std::uint32_t pos = nd.l + static_cast<std::uint32_t>(off);
          s_idx = tree.Points()[static_cast<sjs::usize>(pos)].idx;
        } else {
          s_idx = tree.NodePointPayload(it.node);
        }

        const sjs::usize sid = static_cast<sjs::usize>(s_idx);
        if (sid >= n) {
          std::cerr << "FAIL: sampled payload out of range\n";
          return 1;
        }
        if (!InRange<D, T>(by_id[sid], qlo, qhi)) {
          std::cerr << "FAIL: sampled payload corresponds to an out-of-range point\n";
          return 1;
        }
      }
    }
  }

  std::cout << "PASS\n";
  return 0;
}
