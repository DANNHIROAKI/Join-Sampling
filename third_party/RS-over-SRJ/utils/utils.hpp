#pragma once
// utils/utils.hpp (modified)
//
// The original RS-over-SRJ repository bundled a full demo program with:
//  - global datasets / parameters
//  - CSV I/O
//  - result logging
//
// When vendoring into the main repository as a third_party baseline, those globals
// and I/O routines become problematic (ODR violations, hard-coded paths, etc.).
//
// This header is therefore slimmed down to *only* small, reusable utilities.
// If you still want the original standalone demo, keep a copy of the upstream
// version elsewhere or gate it behind a separate build target.

#include <fstream>
#include <string>
#include <unistd.h>

// Returns resident set size in MB (approx).
inline double process_mem_usage_mb() {
  double resident_set_kb = 0.0;

  // /proc/self/stat: we need vsize (unused) and rss.
  unsigned long vsize = 0;
  long rss = 0;
  {
    std::string ignore;
    std::ifstream ifs("/proc/self/stat", std::ios_base::in);
    if (!ifs.good()) return 0.0;
    ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
        >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
        >> ignore >> ignore >> vsize >> rss;
  }

  const long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
  resident_set_kb = static_cast<double>(rss) * static_cast<double>(page_size_kb);

  return resident_set_kb / 1000.0;  // MB
}