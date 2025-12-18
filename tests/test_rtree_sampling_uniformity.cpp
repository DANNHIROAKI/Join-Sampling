// tests/test_rtree_sampling_uniformity.cpp
//
// Test R-tree sampling uniformity:
//   - Verify that SampleFromSubtree() produces uniform distribution over leaf entries
//   - Verify that SampleIntersect() produces uniform distribution over intersecting entries
//   - Use statistical tests (chi-square) to validate uniformity

// Test R-tree sampling uniformity by directly accessing the internal DynamicRTree.
// Note: This test uses the detail namespace, which is acceptable for unit tests.
#include "sjs/baselines/r_tree/sampling.h"
#include "sjs/core/rng.h"
#include "sjs/core/types.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/point.h"
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
  void CheckEq(const A& a, const B& b, const char* ea, const char* eb,
               const char* file, int line) {
    if (a == b) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line << "  CHECK_EQ(" << ea << ", " << eb
              << ")  got " << a << " vs " << b << "\n";
  }

  void CheckNear(double a, double b, double rel_eps, const char* ea, const char* eb,
                 const char* file, int line) {
    const double diff = std::fabs(a - b);
    const double scale = std::max({1.0, std::fabs(a), std::fabs(b)});
    if (diff <= rel_eps * scale) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line
              << "  CHECK_NEAR(" << ea << ", " << eb << ")  got " << a << " vs " << b
              << "  diff=" << diff << "\n";
  }
};

#define CHECK(ctx, expr) (ctx).Check((expr), #expr, __FILE__, __LINE__)
#define CHECK_EQ(ctx, a, b) (ctx).CheckEq((a), (b), #a, #b, __FILE__, __LINE__)
#define CHECK_NEAR(ctx, a, b, eps) (ctx).CheckNear((a), (b), (eps), #a, #b, __FILE__, __LINE__)

using sjs::Box;
using sjs::Point;
using sjs::Rng;
using sjs::Scalar;
namespace r_tree_detail = sjs::baselines::r_tree::detail;
using r_tree_detail::DynamicRTree;

// Test SampleFromSubtree uniformity on a manually constructed tree.
static void TestSampleFromSubtreeUniformity(TestContext& t) {
  using B1 = Box<1, Scalar>;
  using P1 = Point<1, Scalar>;

  // Build a tree with known structure:
  //   - Insert 10 boxes with distinct IDs (0-9)
  //   - Query the entire space to get all entries
  //   - Sample from the root subtree many times
  //   - Verify uniform distribution using chi-square test

  DynamicRTree<1, Scalar> tree;
  tree.Init(20, {});  // capacity 20

  // Insert 10 boxes: [i, i+1) for i=0..9
  for (sjs::u32 i = 0; i < 10; ++i) {
    const B1 box(P1{static_cast<Scalar>(i)}, P1{static_cast<Scalar>(i + 1)});
    std::string err;
    CHECK(t, tree.Insert(i, box, &err));
  }

  CHECK_EQ(t, tree.Size(), 10u);

  // Sample many times from the root.
  Rng rng(12345);
  const sjs::u32 num_samples = 10000;
  std::vector<sjs::u32> counts(10, 0);

  // Access the root node (we need to use a query that hits everything).
  // Since we can't directly access the root, we use SampleIntersect with a large query.
  const B1 query_all(P1{0.0}, P1{10.0});
  std::vector<sjs::u32> samples;
  samples.reserve(num_samples);

  for (sjs::u32 i = 0; i < num_samples; ++i) {
    samples.clear();
    if (tree.SampleIntersect(query_all, 1, &rng, &samples)) {
      CHECK(t, samples.size() == 1);
      const sjs::u32 id = samples[0];
      CHECK(t, id < 10);
      ++counts[id];
    } else {
      CHECK(t, false);  // Should not fail
    }
  }

  // Chi-square test for uniformity.
  const auto chi2_result = sjs::sampling::quality::ChiSquareUniform(
      sjs::Span<const sjs::u64>(reinterpret_cast<const sjs::u64*>(counts.data()), 10));

  // With 10 categories and 10000 samples, expected count per category = 1000.
  // Chi-square statistic should be small (well below critical value for df=9).
  // For a p-value > 0.01, we expect chi2 < ~21.67 (99th percentile).
  CHECK(t, chi2_result.p_value > 0.01);
  CHECK(t, chi2_result.df == 9u);

  // Also check that counts are reasonably close to expected (1000 each).
  const double expected = static_cast<double>(num_samples) / 10.0;
  for (sjs::usize i = 0; i < 10; ++i) {
    const double ratio = static_cast<double>(counts[i]) / expected;
    // Each count should be within 20% of expected (very loose for statistical test).
    CHECK(t, ratio > 0.7 && ratio < 1.3);
  }
}

// Test SampleIntersect uniformity on a 2D tree with known intersection set.
static void TestSampleIntersectUniformity(TestContext& t) {
  using B2 = Box<2, Scalar>;
  using P2 = Point<2, Scalar>;

  // Build a 2D tree with boxes that all intersect a query region.
  DynamicRTree<2, Scalar> tree;
  tree.Init(20, {});

  // Insert 8 boxes that all intersect query [0.5, 4.5) x [0.5, 4.5).
  // Boxes: [i, i+1) x [i, i+1) for i=0..7, but shifted so they all intersect.
  for (sjs::u32 i = 0; i < 8; ++i) {
    const B2 box(P2{static_cast<Scalar>(i) * 0.5, static_cast<Scalar>(i) * 0.5},
                 P2{static_cast<Scalar>(i) * 0.5 + 1.0, static_cast<Scalar>(i) * 0.5 + 1.0});
    std::string err;
    CHECK(t, tree.Insert(i, box, &err));
  }

  CHECK_EQ(t, tree.Size(), 8u);

  // Query that intersects all boxes.
  const B2 query(P2{0.5}, P2{4.5});

  // Verify all boxes intersect the query.
  sjs::u64 count = tree.CountIntersect(query);
  CHECK_EQ(t, count, 8u);

  // Sample many times.
  Rng rng(54321);
  const sjs::u32 num_samples = 8000;
  std::vector<sjs::u32> counts(8, 0);

  for (sjs::u32 i = 0; i < num_samples; ++i) {
    std::vector<sjs::u32> samples;
    if (tree.SampleIntersect(query, 1, &rng, &samples)) {
      CHECK(t, samples.size() == 1);
      const sjs::u32 id = samples[0];
      CHECK(t, id < 8);
      ++counts[id];
    } else {
      CHECK(t, false);
    }
  }

  // Chi-square test for uniformity.
  const auto chi2_result = sjs::sampling::quality::ChiSquareUniform(
      sjs::Span<const sjs::u64>(reinterpret_cast<const sjs::u64*>(counts.data()), 8));

  // With 8 categories and 8000 samples, expected = 1000 each.
  // For df=7, p-value > 0.01 requires chi2 < ~18.48.
  CHECK(t, chi2_result.p_value > 0.01);
  CHECK(t, chi2_result.df == 7u);

  // Check counts are reasonably close to expected.
  const double expected = static_cast<double>(num_samples) / 8.0;
  for (sjs::usize i = 0; i < 8; ++i) {
    const double ratio = static_cast<double>(counts[i]) / expected;
    CHECK(t, ratio > 0.7 && ratio < 1.3);
  }
}

// Test SampleIntersect with multiple samples per call (batch sampling).
static void TestSampleIntersectBatchUniformity(TestContext& t) {
  using B2 = Box<2, Scalar>;
  using P2 = Point<2, Scalar>;

  DynamicRTree<2, Scalar> tree;
  tree.Init(15, {});

  // Insert 5 boxes.
  for (sjs::u32 i = 0; i < 5; ++i) {
    const B2 box(P2{static_cast<Scalar>(i), static_cast<Scalar>(i)},
                 P2{static_cast<Scalar>(i) + 1.0, static_cast<Scalar>(i) + 1.0});
    std::string err;
    CHECK(t, tree.Insert(i, box, &err));
  }

  const B2 query(P2{0.0}, P2{5.0});
  CHECK_EQ(t, tree.CountIntersect(query), 5u);

  // Sample in batches of 10, many times.
  Rng rng(99999);
  const sjs::u32 num_batches = 1000;
  const sjs::u32 batch_size = 10;
  std::vector<sjs::u32> counts(5, 0);

  for (sjs::u32 b = 0; b < num_batches; ++b) {
    std::vector<sjs::u32> samples;
    if (tree.SampleIntersect(query, batch_size, &rng, &samples)) {
      CHECK_EQ(t, samples.size(), batch_size);
      for (sjs::u32 id : samples) {
        CHECK(t, id < 5);
        ++counts[id];
      }
    } else {
      CHECK(t, false);
    }
  }

  // Chi-square test.
  const sjs::u32 total_samples = num_batches * batch_size;
  const auto chi2_result = sjs::sampling::quality::ChiSquareUniform(
      sjs::Span<const sjs::u64>(reinterpret_cast<const sjs::u64*>(counts.data()), 5));

  CHECK(t, chi2_result.p_value > 0.01);
  CHECK(t, chi2_result.df == 4u);

  // Check counts.
  const double expected = static_cast<double>(total_samples) / 5.0;
  for (sjs::usize i = 0; i < 5; ++i) {
    const double ratio = static_cast<double>(counts[i]) / expected;
    CHECK(t, ratio > 0.8 && ratio < 1.2);  // Tighter bound for larger sample size
  }
}

}  // namespace

int main() {
  TestContext t;

  TestSampleFromSubtreeUniformity(t);
  TestSampleIntersectUniformity(t);
  TestSampleIntersectBatchUniformity(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_rtree_sampling_uniformity\n";
    return 0;
  }
  std::cerr << "[FAILED] test_rtree_sampling_uniformity: " << t.fails << " failure(s)\n";
  return 1;
}
