// tests/test_dataset_validation.cpp
//
// Validation tests for dataset / I/O contracts that are critical to the
// SJS-nD semantics: unique object ids and finite coordinates.

#include "geometry/box.h"
#include "io/binary_io.h"
#include "io/csv_io.h"
#include "io/dataset.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

namespace fs = std::filesystem;

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
using sjs::Relation;
using sjs::Scalar;

Box<2, Scalar> MakeBox(double x0, double y0, double x1, double y1) {
  Point<2, Scalar> lo{};
  Point<2, Scalar> hi{};
  lo[0] = x0; lo[1] = y0;
  hi[0] = x1; hi[1] = y1;
  return Box<2, Scalar>(lo, hi);
}

void TestDuplicateIdsRejected(TestContext& t) {
  Relation<2, Scalar> rel;
  rel.Add(MakeBox(0.0, 0.0, 1.0, 1.0), 7);
  rel.Add(MakeBox(1.0, 1.0, 2.0, 2.0), 7);

  std::string err;
  CHECK(t, !rel.Validate(true, &err));
  CHECK(t, err.find("Duplicate id") != std::string::npos);
}

void TestNonFiniteCoordinatesRejected(TestContext& t) {
  Relation<2, Scalar> rel_nan;
  rel_nan.Add(MakeBox(std::numeric_limits<double>::quiet_NaN(), 0.0, 1.0, 1.0), 0);
  std::string err_nan;
  CHECK(t, !rel_nan.Validate(true, &err_nan));
  CHECK(t, err_nan.find("Non-finite") != std::string::npos);

  Relation<2, Scalar> rel_inf;
  rel_inf.Add(MakeBox(0.0, 0.0, std::numeric_limits<double>::infinity(), 1.0), 1);
  std::string err_inf;
  CHECK(t, !rel_inf.Validate(true, &err_inf));
  CHECK(t, err_inf.find("Non-finite") != std::string::npos);
}

void TestMixedExplicitImplicitIdsStayUnique(TestContext& t) {
  Relation<2, Scalar> rel;
  rel.Add(MakeBox(0.0, 0.0, 1.0, 1.0), 0);
  rel.Add(MakeBox(1.0, 1.0, 2.0, 2.0));
  rel.Add(MakeBox(2.0, 2.0, 3.0, 3.0), 5);
  rel.Add(MakeBox(3.0, 3.0, 4.0, 4.0));

  std::string err;
  CHECK(t, rel.Validate(true, &err));
  CHECK_EQ(t, rel.ids.size(), 4u);
  CHECK(t, rel.ids[0] != rel.ids[1]);
  CHECK(t, rel.ids[1] != rel.ids[2]);
  CHECK(t, rel.ids[2] != rel.ids[3]);
}

void TestCsvRejectsInvalidRows(TestContext& t) {
  const fs::path tmpdir = fs::temp_directory_path();
  const fs::path dup_path = tmpdir / "sjs_test_dup_ids.csv";
  const fs::path nan_path = tmpdir / "sjs_test_nan.csv";

  {
    std::ofstream out(dup_path);
    out << "id,lo0,lo1,hi0,hi1\n";
    out << "3,0,0,1,1\n";
    out << "3,1,1,2,2\n";
  }
  {
    std::ofstream out(nan_path);
    out << "id,lo0,lo1,hi0,hi1\n";
    out << "1,nan,0,1,1\n";
  }

  Relation<2, Scalar> rel;
  std::string err;
  CHECK(t, (!sjs::csv::ReadBoxesSimple<2, Scalar>(dup_path.string(), &rel, ',', true, &err)));
  CHECK(t, err.find("Duplicate id") != std::string::npos);

  err.clear();
  CHECK(t, (!sjs::csv::ReadBoxesSimple<2, Scalar>(nan_path.string(), &rel, ',', true, &err)));
  CHECK(t, err.find("parse") != std::string::npos || err.find("Non-finite") != std::string::npos);

  std::error_code ec;
  fs::remove(dup_path, ec);
  fs::remove(nan_path, ec);
}

void TestBinaryWriterRejectsInvalidRelation(TestContext& t) {
  Relation<2, Scalar> rel;
  rel.Add(MakeBox(0.0, 0.0, 1.0, 1.0), 4);
  rel.Add(MakeBox(1.0, 1.0, 2.0, 2.0), 4);

  const fs::path out_path = fs::temp_directory_path() / "sjs_test_invalid_relation.bin";
  std::string err;
  CHECK(t, (!sjs::binary::WriteRelationBinary<2, Scalar>(out_path.string(), rel, sjs::binary::BinaryWriteOptions{}, &err)));
  CHECK(t, err.find("invalid relation") != std::string::npos || err.find("Duplicate id") != std::string::npos);

  std::error_code ec;
  fs::remove(out_path, ec);
}

void TestDatasetValidationMessagePrefixesRelation(TestContext& t) {
  Dataset<2, Scalar> ds;
  ds.name = "bad_ds";
  ds.R.Add(MakeBox(0.0, 0.0, 1.0, 1.0), 2);
  ds.R.Add(MakeBox(1.0, 1.0, 2.0, 2.0), 2);
  ds.S.Add(MakeBox(0.0, 0.0, 1.0, 1.0), 0);

  std::string err;
  CHECK(t, !ds.Validate(true, &err));
  CHECK(t, err.rfind("R:", 0) == 0);
}

}  // namespace

int main() {
  TestContext t;
  TestDuplicateIdsRejected(t);
  TestNonFiniteCoordinatesRejected(t);
  TestMixedExplicitImplicitIdsStayUnique(t);
  TestCsvRejectsInvalidRows(t);
  TestBinaryWriterRejectsInvalidRelation(t);
  TestDatasetValidationMessagePrefixesRelation(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_dataset_validation\n";
    return 0;
  }
  std::cerr << "[FAILED] test_dataset_validation: " << t.fails << " failure(s)\n";
  return 1;
}
