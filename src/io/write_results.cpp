// src/io/write_results.cpp
//
// Experiment output utilities.
//
// This TU provides small, dependency-light writers for:
//  - Per-run results (TSV/CSV)
//  - Per-run JSONL (for programmatic post-processing)
//  - Optional: sampled join pairs (TSV)
//  - Optional: summary rows (mean/std/quantiles) over multiple repeats
//
// We intentionally keep the output schema stable and append-friendly.

#include "sjs/core/assert.h"
#include "sjs/core/logging.h"
#include "sjs/core/stats.h"
#include "sjs/core/types.h"

#include "sjs/baselines/baseline_api.h"
#include "sjs/io/csv_io.h"
#include "sjs/io/write_results.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace io {

namespace {

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

inline bool EnsureDirExists(const std::filesystem::path& dir, std::string* err) {
  std::error_code ec;
  if (dir.empty()) return true;
  if (std::filesystem::exists(dir, ec)) {
    if (std::filesystem::is_directory(dir, ec)) return true;
    SetErr(err, "Path exists but is not a directory: " + dir.string());
    return false;
  }
  if (!std::filesystem::create_directories(dir, ec)) {
    SetErr(err, "Failed to create directory: " + dir.string() + " (" + ec.message() + ")");
    return false;
  }
  return true;
}

inline bool FileNonEmpty(const std::filesystem::path& p) {
  std::error_code ec;
  if (!std::filesystem::exists(p, ec)) return false;
  const auto sz = std::filesystem::file_size(p, ec);
  return !ec && sz > 0;
}

inline std::string SanitizeFilename(std::string_view s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    // Replace path separators and whitespace with '_'.
    if (c == '/' || c == '\\' || c == ':' || c == '*' || c == '?' || c == '"' || c == '<' || c == '>' || c == '|') {
      out.push_back('_');
    } else if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
      out.push_back('_');
    } else {
      out.push_back(c);
    }
  }
  if (out.empty()) out = "out";
  return out;
}

inline std::string NowIso8601Local() {
  using clock = std::chrono::system_clock;
  const auto now = clock::now();
  const auto t = clock::to_time_t(now);

  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &t);
#else
  localtime_r(&t, &tm);
#endif

  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S");
  return oss.str();
}

inline std::string JsonEscape(std::string_view s) {
  std::ostringstream oss;
  for (char c : s) {
    switch (c) {
      case '"': oss << "\\\""; break;
      case '\\': oss << "\\\\"; break;
      case '\b': oss << "\\b"; break;
      case '\f': oss << "\\f"; break;
      case '\n': oss << "\\n"; break;
      case '\r': oss << "\\r"; break;
      case '\t': oss << "\\t"; break;
      default:
        // Control chars -> \u00XX
        if (static_cast<unsigned char>(c) < 0x20) {
          oss << "\\u" << std::hex << std::setw(4) << std::setfill('0')
              << static_cast<int>(static_cast<unsigned char>(c))
              << std::dec;
        } else {
          oss << c;
        }
    }
  }
  return oss.str();
}

inline std::string ToStringOrEmpty(bool ok, const std::string& s) {
  return ok ? s : std::string();
}

inline double GetTopLevelRunTotalMs(const baselines::RunReport& r) {
  // Sum phases that start with "run_". This avoids double counting when
  // baselines record additional fine-grained phases.
  double ms = 0.0;
  for (const auto& kv : r.phases.SnapshotNanosSorted()) {
    const std::string& name = kv.first;
    if (name.rfind("run_", 0) == 0) {
      ms += static_cast<double>(kv.second) / 1e6;
    }
  }
  return ms;
}

inline std::vector<std::string> RunReportHeader() {
  return {
      "ts_local",
      "dataset",
      "method",
      "variant",
      "baseline",
      "seed",
      "t",
      "ok",
      "error",
      "count_value",
      "count_exact",
      "count_stderr",
      "count_ci_low",
      "count_ci_high",
      "count_aux_draws",
      "used_enumeration",
      "enumeration_truncated",
      "enumeration_cap",
      "enum_pairs_pass1",
      "enum_pairs_pass2",
      "adaptive_branch",
      "adaptive_pilot_pairs",
      "run_total_ms",
      "phases_json",
      "note",
  };
}

inline std::vector<std::string> RunReportRow(const baselines::RunReport& r) {
  std::vector<std::string> row;
  row.reserve(32);

  row.push_back(NowIso8601Local());
  row.push_back(r.dataset_name);
  row.push_back(std::string(ToString(r.method)));
  row.push_back(std::string(ToString(r.variant)));
  row.push_back(r.baseline_name);
  row.push_back(std::to_string(r.seed));
  row.push_back(std::to_string(r.t));
  row.push_back(r.ok ? "1" : "0");
  row.push_back(ToStringOrEmpty(!r.ok, r.error));

  // CountResult
  row.push_back(std::to_string(static_cast<double>(r.count.value)));
  row.push_back(r.count.exact ? "1" : "0");

  auto num_or_empty = [](long double x) -> std::string {
    if (!std::isfinite(static_cast<double>(x))) return std::string();
    return std::to_string(static_cast<double>(x));
  };
  row.push_back(num_or_empty(r.count.stderr));
  row.push_back(num_or_empty(r.count.ci_low));
  row.push_back(num_or_empty(r.count.ci_high));
  row.push_back(std::to_string(r.count.aux_draws));

  row.push_back(r.used_enumeration ? "1" : "0");
  row.push_back(r.enumeration_truncated ? "1" : "0");
  row.push_back(std::to_string(r.enumeration_cap));
  row.push_back(std::to_string(r.enumeration_pairs_pass1));
  row.push_back(std::to_string(r.enumeration_pairs_pass2));
  row.push_back(r.adaptive_branch);
  row.push_back(std::to_string(r.adaptive_pilot_pairs));

  row.push_back(std::to_string(GetTopLevelRunTotalMs(r)));

  // Phases: store as a compact JSON string in a single cell.
  row.push_back(r.phases.ToJsonMillis());

  row.push_back(r.note);

  return row;
}

inline bool WriteRow(std::ofstream& out,
                     const std::vector<std::string>& cols,
                     const csv::Dialect& d,
                     std::string* err) {
  if (!out) {
    SetErr(err, "WriteRow: output stream not ready");
    return false;
  }
  for (usize i = 0; i < cols.size(); ++i) {
    if (i) out << d.sep;
    out << csv::EscapeCell(cols[i], d);
  }
  out << "\n";
  if (!out) {
    SetErr(err, "WriteRow: write failed");
    return false;
  }
  return true;
}

inline std::string RunReportToJsonLine(const baselines::RunReport& r) {
  // Minimal JSON schema (one object per line).
  // Strings are escaped; numeric fields are written as numbers.
  // Note: phases_json is already JSON, so we embed it as a JSON object (not a string).
  std::ostringstream oss;
  oss << "{";
  oss << "\"ts_local\":\"" << JsonEscape(NowIso8601Local()) << "\",";
  oss << "\"dataset\":\"" << JsonEscape(r.dataset_name) << "\",";
  oss << "\"method\":\"" << JsonEscape(std::string(ToString(r.method))) << "\",";
  oss << "\"variant\":\"" << JsonEscape(std::string(ToString(r.variant))) << "\",";
  oss << "\"baseline\":\"" << JsonEscape(r.baseline_name) << "\",";
  oss << "\"seed\":" << r.seed << ",";
  oss << "\"t\":" << r.t << ",";
  oss << "\"ok\":" << (r.ok ? "true" : "false") << ",";
  oss << "\"error\":\"" << JsonEscape(r.error) << "\",";

  // Count
  oss << "\"count\":" << r.count.ToJsonLite() << ",";

  // Metadata
  oss << "\"used_enumeration\":" << (r.used_enumeration ? "true" : "false") << ",";
  oss << "\"enumeration_truncated\":" << (r.enumeration_truncated ? "true" : "false") << ",";
  oss << "\"enumeration_cap\":" << r.enumeration_cap << ",";
  oss << "\"enum_pairs_pass1\":" << r.enumeration_pairs_pass1 << ",";
  oss << "\"enum_pairs_pass2\":" << r.enumeration_pairs_pass2 << ",";
  oss << "\"adaptive_branch\":\"" << JsonEscape(r.adaptive_branch) << "\",";
  oss << "\"adaptive_pilot_pairs\":" << r.adaptive_pilot_pairs << ",";

  // Derived metrics
  oss << "\"run_total_ms\":" << GetTopLevelRunTotalMs(r) << ",";

  // Phases: embed as JSON object.
  oss << "\"phases\":" << r.phases.ToJsonMillis() << ",";

  oss << "\"note\":\"" << JsonEscape(r.note) << "\"";
  oss << "}";
  return oss.str();
}

}  // namespace

// --------------------------
// Public writers
// --------------------------

namespace {

bool AppendRunReportDelimited(const std::string& path,
                             const baselines::RunReport& report,
                             char sep,
                             std::string* err) {
  const std::filesystem::path p(path);
  if (!EnsureDirExists(p.parent_path(), err)) return false;

  const bool need_header = !FileNonEmpty(p);

  std::ofstream out(path, std::ios::out | std::ios::app);
  if (!out) {
    SetErr(err, "Cannot open results file for append: " + path);
    return false;
  }

  csv::Dialect d;
  d.sep = sep;
  d.quote = '"';
  d.always_quote = false;
  d.write_header = true;

  if (need_header) {
    if (!WriteRow(out, RunReportHeader(), d, err)) return false;
  }
  if (!WriteRow(out, RunReportRow(report), d, err)) return false;
  return true;
}

}  // namespace

bool AppendRunReportTSV(const std::string& path,
                        const baselines::RunReport& report,
                        std::string* err) {
  return AppendRunReportDelimited(path, report, '\t', err);
}

bool AppendRunReportCSV(const std::string& path,
                        const baselines::RunReport& report,
                        std::string* err) {
  return AppendRunReportDelimited(path, report, ',', err);
}

bool AppendRunReportJSONL(const std::string& path,
                          const baselines::RunReport& report,
                          std::string* err) {
  const std::filesystem::path p(path);
  if (!EnsureDirExists(p.parent_path(), err)) return false;

  std::ofstream out(path, std::ios::out | std::ios::app);
  if (!out) {
    SetErr(err, "Cannot open JSONL for append: " + path);
    return false;
  }

  out << RunReportToJsonLine(report) << "\n";
  if (!out) {
    SetErr(err, "JSONL write failed: " + path);
    return false;
  }
  return true;
}

// Write sampled join pairs (r_id,s_id) as a TSV file under out_dir.
// Returns the written path via out_path if provided.
bool WriteSamplesTSV(const std::string& out_dir,
                     const baselines::RunReport& report,
                     std::string* out_path,
                     std::string* err) {
  if (report.samples.pairs.empty()) {
    if (out_path) out_path->clear();
    return true;
  }

  const std::filesystem::path dir = std::filesystem::path(out_dir) / "samples";
  if (!EnsureDirExists(dir, err)) return false;

  std::ostringstream fname;
  fname << SanitizeFilename(report.dataset_name) << "_"
        << SanitizeFilename(std::string(ToString(report.method))) << "_"
        << SanitizeFilename(std::string(ToString(report.variant))) << "_seed" << report.seed
        << "_t" << report.t
        << ".tsv";

  const std::filesystem::path path = dir / fname.str();

  std::string io_err;
  if (!csv::WritePairs(path.string(), Span<const PairId>(report.samples.pairs),
                       csv::Dialect{'\t'},
                       &io_err)) {
    SetErr(err, "WriteSamplesTSV failed: " + io_err);
    return false;
  }

  if (out_path) *out_path = path.string();
  return true;
}

// A small summary writer for multiple repeats.
// The caller groups runs as desired (e.g., same dataset/method/variant).
bool WriteSummaryTSV(const std::string& path,
                     Span<const baselines::RunReport> runs,
                     std::string* err) {
  const std::filesystem::path p(path);
  if (!EnsureDirExists(p.parent_path(), err)) return false;

  // Compute derived per-run metrics.
  std::vector<double> total_ms;
  total_ms.reserve(runs.size());

  std::vector<double> count_est;
  count_est.reserve(runs.size());

  std::string dataset;
  std::string method;
  std::string variant;
  std::string baseline;

  for (usize i = 0; i < runs.size(); ++i) {
    const auto& r = runs[i];
    if (!r.ok) continue;  // Only summarize successful runs by default.
    total_ms.push_back(GetTopLevelRunTotalMs(r));
    count_est.push_back(static_cast<double>(r.count.value));

    if (dataset.empty()) dataset = r.dataset_name;
    if (method.empty()) method = std::string(ToString(r.method));
    if (variant.empty()) variant = std::string(ToString(r.variant));
    if (baseline.empty()) baseline = r.baseline_name;
  }

  if (total_ms.empty()) {
    SetErr(err, "No successful runs to summarize");
    return false;
  }

  const Summary s_ms = Summarize(total_ms);
  const Summary s_cnt = Summarize(count_est);

  const bool need_header = !FileNonEmpty(p);
  std::ofstream out(path, std::ios::out | std::ios::app);
  if (!out) {
    SetErr(err, "Cannot open summary TSV for append: " + path);
    return false;
  }

  csv::Dialect d;
  d.sep = '\t';

  auto header = std::vector<std::string>{
      "ts_local",
      "dataset",
      "method",
      "variant",
      "baseline",
      "n_ok",
      "total_ms_mean",
      "total_ms_stdev",
      "total_ms_median",
      "total_ms_p90",
      "total_ms_p95",
      "total_ms_p99",
      "count_mean",
      "count_stdev",
      "count_median",
  };

  if (need_header) {
    if (!WriteRow(out, header, d, err)) return false;
  }

  std::vector<std::string> row;
  row.reserve(header.size());
  row.push_back(NowIso8601Local());
  row.push_back(dataset);
  row.push_back(method);
  row.push_back(variant);
  row.push_back(baseline);
  row.push_back(std::to_string(s_ms.n));
  row.push_back(std::to_string(s_ms.mean));
  row.push_back(std::to_string(s_ms.stdev));
  row.push_back(std::to_string(s_ms.median));
  row.push_back(std::to_string(s_ms.p90));
  row.push_back(std::to_string(s_ms.p95));
  row.push_back(std::to_string(s_ms.p99));
  row.push_back(std::to_string(s_cnt.mean));
  row.push_back(std::to_string(s_cnt.stdev));
  row.push_back(std::to_string(s_cnt.median));

  if (!WriteRow(out, row, d, err)) return false;
  return true;
}

}  // namespace io
}  // namespace sjs
