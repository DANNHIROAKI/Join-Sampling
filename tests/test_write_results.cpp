// tests/test_write_results.cpp
//
// Regression tests for write_results utilities.

#include "sjs/baselines/baseline_api.h"
#include "sjs/core/types.h"
#include "sjs/io/write_results.h"

#include <filesystem>
#include <fstream>
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
};

#define CHECK(ctx, expr) (ctx).Check((expr), #expr, __FILE__, __LINE__)

}  // namespace

int main() {
  TestContext t;

  const fs::path tmp_root = fs::temp_directory_path() / "sjs_write_results_test";
  std::error_code ec;
  fs::remove_all(tmp_root, ec);
  fs::create_directories(tmp_root);

  // Create a file that will be (incorrectly) treated as a directory path.
  const fs::path blocker = tmp_root / "blocker.txt";
  {
    std::ofstream f(blocker.string());
    f << "x";
  }

  const fs::path target = blocker / "run.csv";

  sjs::baselines::RunReport report;
  report.ok = true;
  report.dataset_name = "demo";
  report.baseline_name = "dummy";

  std::string err;
  const bool ok = sjs::io::AppendRunReportCSV(target.string(), report, &err);

  CHECK(t, !ok);
  CHECK(t, !err.empty());
  CHECK(t, err.find("not a directory") != std::string::npos);
  CHECK(t, !fs::exists(target));

  // Summaries should reject datasets with no successful runs.
  const fs::path summary_path = tmp_root / "summary.tsv";
  std::vector<sjs::baselines::RunReport> runs(2);
  runs[0].ok = false;
  runs[1].ok = false;
  std::string summary_err;
  const bool summary_ok = sjs::io::WriteSummaryTSV(
      summary_path.string(),
      sjs::Span<const sjs::baselines::RunReport>(runs.data(), runs.size()),
      &summary_err);

  CHECK(t, !summary_ok);
  CHECK(t, !summary_err.empty());
  CHECK(t, summary_err.find("No successful runs") != std::string::npos);
  CHECK(t, !fs::exists(summary_path));

  fs::remove_all(tmp_root, ec);

  if (t.fails == 0) {
    std::cout << "[OK] test_write_results\n";
    return 0;
  }
  std::cerr << "[FAILED] test_write_results: " << t.fails << " failure(s)\n";
  return 1;
}

