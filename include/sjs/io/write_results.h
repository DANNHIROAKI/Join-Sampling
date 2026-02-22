#pragma once
// sjs/io/write_results.h
//
// Experiment output utilities.
//
// This header intentionally contains declarations only.
// Implementations live in: src/io/write_results.cpp
//
// Writers provided:
//  - AppendRunReportTSV / CSV: append one per-run row (header written on first write)
//  - AppendRunReportJSONL: append one JSON object per line (JSONL)
//  - WriteSamplesTSV: write sampled join pairs to out_dir/samples/
//  - WriteSummaryTSV: append aggregated stats across repeats
//
// HighDims note:
//   baselines::RunReport records the compile-time dimension used by the run.
//   The writers emit a "dim" column/field accordingly.

#include "sjs/baselines/baseline_api.h"  // baselines::RunReport, sjs::Span

#include <string>

namespace sjs::io {

// Append one RunReport as a TSV row. If the file does not exist (or is empty),
// a header row is written first.
bool AppendRunReportTSV(const std::string& path,
                        const baselines::RunReport& report,
                        std::string* err = nullptr);

// Append one RunReport as a CSV row. If the file does not exist (or is empty),
// a header row is written first.
bool AppendRunReportCSV(const std::string& path,
                        const baselines::RunReport& report,
                        std::string* err = nullptr);

// Append one RunReport as a JSON line (JSONL). Each line is one JSON object.
bool AppendRunReportJSONL(const std::string& path,
                          const baselines::RunReport& report,
                          std::string* err = nullptr);

// Write sampled join pairs (r_id,s_id) as a TSV file under out_dir/samples/.
// - If report.samples is empty, returns true and writes nothing.
// - If out_path != nullptr, returns the written path (or empty if nothing written).
bool WriteSamplesTSV(const std::string& out_dir,
                     const baselines::RunReport& report,
                     std::string* out_path = nullptr,
                     std::string* err = nullptr);

// Append a summary row (mean/stdev/quantiles) over multiple runs to a TSV file.
// Caller should group runs (e.g. same dataset/method/variant/dim) before calling.
bool WriteSummaryTSV(const std::string& path,
                     Span<const baselines::RunReport> runs,
                     std::string* err = nullptr);

}  // namespace sjs::io
