#pragma once
// sjs/io/write_results.h
//
// Light-weight writers for experiment outputs (CSV/TSV/JSONL) and sampled pairs.
// These helpers keep the output schema consistent across apps and tests while
// avoiding heavy dependencies.

#include "sjs/baselines/baseline_api.h"

#include <string>

namespace sjs {
namespace io {

bool AppendRunReportTSV(const std::string& path,
                        const baselines::RunReport& report,
                        std::string* err = nullptr);

bool AppendRunReportCSV(const std::string& path,
                        const baselines::RunReport& report,
                        std::string* err = nullptr);

bool AppendRunReportJSONL(const std::string& path,
                          const baselines::RunReport& report,
                          std::string* err = nullptr);

// Write sampled join pairs (if any) to a TSV file under out_dir.
// Returns the written path via out_path when provided.
bool WriteSamplesTSV(const std::string& out_dir,
                     const baselines::RunReport& report,
                     std::string* out_path = nullptr,
                     std::string* err = nullptr);

// Append a summary row (mean/stdev/quantiles) over multiple runs to a TSV file.
bool WriteSummaryTSV(const std::string& path,
                     Span<const baselines::RunReport> runs,
                     std::string* err = nullptr);

}  // namespace io
}  // namespace sjs

