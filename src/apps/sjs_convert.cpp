// src/apps/sjs_convert.cpp
//
// Dataset format conversion utility.
// Supports converting *pairs* of relations (R,S) between:
//   - CSV/TSV (simple box format)
//   - SJS binary (.bin) format
//
// Example:
//   ./sjs_convert --in_fmt=csv --out_fmt=binary \
//                 --in_r=R.csv --in_s=S.csv \
//                 --out_r=R.bin --out_s=S.bin
//
// Example (binary -> csv):
//   ./sjs_convert --in_fmt=binary --out_fmt=csv \
//                 --in_r=R.bin --in_s=S.bin \
//                 --out_r=R.csv --out_s=S.csv --sep=tab

#include "sjs/core/logging.h"
#include "sjs/core/config.h"
#include "sjs/core/types.h"

#include "sjs/io/binary_io.h"
#include "sjs/io/csv_io.h"
#include "sjs/io/dataset.h"

#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>

namespace fs = std::filesystem;

namespace sjs {
namespace apps {

namespace {

inline bool IsHelpRequested(const sjs::ArgMap& args) {
  return args.Has("help") || args.Has("h") || args.Has("-h") || args.Has("--help");
}

inline void PrintUsage() {
  std::cerr
      << "sjs_convert: convert dataset formats (R,S)\n\n"
      << "Required:\n"
      << "  --in_fmt=csv|binary\n"
      << "  --out_fmt=csv|binary\n"
      << "  --in_r=<path> --in_s=<path>\n"
      << "  --out_r=<path> --out_s=<path>\n"
      << "\nOptional:\n"
      << "  --dim=2                (default 2; this build supports only 2)\n"
      << "  --sep=,|tab|\\t          (for CSV read/write; default ',')\n"
      << "\nCSV box format (simple):\n"
      << "  id, lo0, lo1, ..., hi0, hi1, ...\n"
      << "Header optional; quotes not supported.\n\n";
}

inline std::string ToLower(std::string s) {
  for (char& c : s) {
    if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
  }
  return s;
}

inline char ParseSep(const std::optional<std::string_view>& v) {
  if (!v) return ',';
  if (*v == "tab" || *v == "\\t") return '\t';
  if (!v->empty()) return (*v)[0];
  return ',';
}

inline bool ParseI32(const std::optional<std::string_view>& v, i32* out) {
  if (!v || !out) return false;
  try {
    const std::string sv(*v);
    *out = static_cast<sjs::i32>(std::stoi(sv));
    return true;
  } catch (...) {
    return false;
  }
}

inline bool EnsureDir(const fs::path& p) {
  if (p.empty()) return true;
  std::error_code ec;
  fs::create_directories(p, ec);
  return !ec;
}

}  // namespace

}  // namespace apps
}  // namespace sjs

int main(int argc, char** argv) {
  sjs::ArgMap args = sjs::ArgMap::FromArgv(argc, argv);
  if (sjs::apps::IsHelpRequested(args)) {
    sjs::apps::PrintUsage();
    return 0;
  }

  const auto in_fmt_opt = args.Get("in_fmt");
  const auto out_fmt_opt = args.Get("out_fmt");
  const auto in_r = args.Get("in_r");
  const auto in_s = args.Get("in_s");
  const auto out_r = args.Get("out_r");
  const auto out_s = args.Get("out_s");

  if (!in_fmt_opt || !out_fmt_opt || !in_r || !in_s || !out_r || !out_s) {
    SJS_LOG_ERROR("Missing required flags.");
    sjs::apps::PrintUsage();
    return 2;
  }

  sjs::i32 dim = 2;
  sjs::apps::ParseI32(args.Get("dim"), &dim);
  if (dim != 2) {
    SJS_LOG_ERROR("This build supports only --dim=2 for conversion.");
    return 2;
  }

  const std::string in_fmt = sjs::apps::ToLower(std::string(*in_fmt_opt));
  const std::string out_fmt = sjs::apps::ToLower(std::string(*out_fmt_opt));
  const char sep = sjs::apps::ParseSep(args.Get("sep"));

  // Ensure output dirs.
  if (!sjs::apps::EnsureDir(fs::path(std::string(*out_r)).parent_path()) ||
      !sjs::apps::EnsureDir(fs::path(std::string(*out_s)).parent_path())) {
    SJS_LOG_ERROR("Failed to create output directories.");
    return 3;
  }

  // Load relations (Dim=2).
  sjs::Relation<2, sjs::Scalar> R, S;
  std::string err;

  if (in_fmt == "binary") {
    sjs::binary::BinaryReadOptions opt;
    opt.generate_ids_if_missing = true;
    opt.drop_empty = false;
    if (!sjs::binary::ReadRelationBinary<2, sjs::Scalar>(std::string(*in_r), &R, nullptr, opt, &err)) {
      SJS_LOG_ERROR("Read binary R failed:", err);
      return 4;
    }
    if (!sjs::binary::ReadRelationBinary<2, sjs::Scalar>(std::string(*in_s), &S, nullptr, opt, &err)) {
      SJS_LOG_ERROR("Read binary S failed:", err);
      return 4;
    }
  } else if (in_fmt == "csv") {
    if (!sjs::csv::ReadBoxesSimple<2, sjs::Scalar>(std::string(*in_r), &R, sep, /*has_header=*/true, &err)) {
      SJS_LOG_ERROR("Read CSV R failed:", err);
      return 4;
    }
    if (!sjs::csv::ReadBoxesSimple<2, sjs::Scalar>(std::string(*in_s), &S, sep, /*has_header=*/true, &err)) {
      SJS_LOG_ERROR("Read CSV S failed:", err);
      return 4;
    }
  } else {
    SJS_LOG_ERROR("Unknown --in_fmt:", in_fmt);
    return 2;
  }

  // Write.
  if (out_fmt == "binary") {
    sjs::binary::BinaryWriteOptions opt;
    opt.half_open = true;
    opt.write_ids = true;
    opt.scalar = sjs::binary::ScalarEncoding::Float64;
    opt.write_name = true;
    if (!sjs::binary::WriteRelationBinary<2, sjs::Scalar>(std::string(*out_r), R, opt, &err)) {
      SJS_LOG_ERROR("Write binary R failed:", err);
      return 5;
    }
    if (!sjs::binary::WriteRelationBinary<2, sjs::Scalar>(std::string(*out_s), S, opt, &err)) {
      SJS_LOG_ERROR("Write binary S failed:", err);
      return 5;
    }
  } else if (out_fmt == "csv") {
    sjs::csv::Dialect d;
    d.sep = sep;
    d.write_header = true;
    if (!sjs::csv::WriteBoxes<2, sjs::Scalar>(std::string(*out_r), R, d, &err)) {
      SJS_LOG_ERROR("Write CSV R failed:", err);
      return 5;
    }
    if (!sjs::csv::WriteBoxes<2, sjs::Scalar>(std::string(*out_s), S, d, &err)) {
      SJS_LOG_ERROR("Write CSV S failed:", err);
      return 5;
    }
  } else {
    SJS_LOG_ERROR("Unknown --out_fmt:", out_fmt);
    return 2;
  }

  SJS_LOG_INFO("Converted OK. out_r=", std::string(*out_r), " out_s=", std::string(*out_s));
  return 0;
}
