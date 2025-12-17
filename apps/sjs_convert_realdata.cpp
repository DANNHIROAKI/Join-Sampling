// apps/sjs_convert_realdata.cpp
//
// Real-data conversion entry point.
//
// Goal:
//   Convert real-world spatial datasets (OSM/TIGER/GeoJSON/...) into the
//   internal SJS binary format used by the experiment framework.
//
// Status:
//   The real-data loaders are currently stubs (see include/sjs/io/realdata_stub.h).
//   This app wires CLI + export plumbing so that once you implement a loader
//   (e.g., backed by GDAL/OGR or libosmium), the rest of the system works
//   unchanged.
//
// Example (future, after implementing loaders):
//   ./sjs_convert_realdata --src_r=geojson --in_r=roads.geojson \
//                         --src_s=geojson --in_s=buildings.geojson \
//                         --out_r=data/real/roads.bin --out_s=data/real/buildings.bin

#include "sjs/core/logging.h"
#include "sjs/core/config.h"
#include "sjs/core/types.h"

#include "sjs/io/binary_io.h"
#include "sjs/io/dataset.h"
#include "sjs/io/realdata_stub.h"

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
      << "sjs_convert_realdata: convert real datasets to SJS binary (R,S)\n\n"
      << "Required:\n"
      << "  --src_r=<osm_pbf|tiger_shp|geojson|wkt_csv>  --in_r=<path>\n"
      << "  --src_s=<osm_pbf|tiger_shp|geojson|wkt_csv>  --in_s=<path>\n"
      << "  --out_r=<path.bin> --out_s=<path.bin>\n"
      << "\nOptional:\n"
      << "  --dataset=<name>           (default: realdata)\n"
      << "  --limit_r=<u64> --limit_s=<u64>   (0 = no limit)\n"
      << "  --drop_empty=0|1           (default 1)\n"
      << "  --to_mbr=0|1               (default 1)\n"
      << "  --dim=2                    (this build supports only 2)\n"
      << "\nNotes:\n"
      << "  Real-data importers are stubs right now. Implement them in\n"
      << "  include/sjs/io/realdata_stub.h (or replace the stub with GDAL/libosmium).\n";
}

inline std::string ToLower(std::string s) {
  for (char& c : s) {
    if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
  }
  return s;
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

inline bool ParseU64(const std::optional<std::string_view>& v, u64* out) {
  if (!v || !out) return false;
  try {
    std::size_t idx = 0;
    const std::string sv(*v);
    unsigned long long x = std::stoull(sv, &idx, 10);
    if (idx != sv.size()) return false;
    *out = static_cast<u64>(x);
    return true;
  } catch (...) {
    return false;
  }
}

inline bool ParseBool(const std::optional<std::string_view>& v, bool def) {
  if (!v) return def;
  if (*v == "1" || *v == "true" || *v == "yes" || *v == "on") return true;
  if (*v == "0" || *v == "false" || *v == "no" || *v == "off") return false;
  return def;
}

inline bool ParseSource(std::string_view s, sjs::realdata::Source* out) {
  if (!out) return false;
  const std::string t = ToLower(std::string(s));
  if (t == "osm_pbf" || t == "osmpbf" || t == "pbf") {
    *out = sjs::realdata::Source::OSM_PBF;
    return true;
  }
  if (t == "tiger_shp" || t == "tiger" || t == "shp" || t == "shapefile") {
    *out = sjs::realdata::Source::TIGER_SHP;
    return true;
  }
  if (t == "geojson" || t == "json") {
    *out = sjs::realdata::Source::GEOJSON;
    return true;
  }
  if (t == "wkt_csv" || t == "wkt" || t == "csv_wkt") {
    *out = sjs::realdata::Source::WKT_CSV;
    return true;
  }
  *out = sjs::realdata::Source::Unknown;
  return false;
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

  const auto src_r_s = args.Get("src_r");
  const auto src_s_s = args.Get("src_s");
  const auto in_r = args.Get("in_r");
  const auto in_s = args.Get("in_s");
  const auto out_r = args.Get("out_r");
  const auto out_s = args.Get("out_s");

  if (!src_r_s || !src_s_s || !in_r || !in_s || !out_r || !out_s) {
    SJS_LOG_ERROR("Missing required flags.");
    sjs::apps::PrintUsage();
    return 2;
  }

  sjs::i32 dim = 2;
  sjs::apps::ParseI32(args.Get("dim"), &dim);
  if (dim != 2) {
    SJS_LOG_ERROR("This build supports only --dim=2.");
    return 2;
  }

  sjs::realdata::Source src_r = sjs::realdata::Source::Unknown;
  sjs::realdata::Source src_s = sjs::realdata::Source::Unknown;
  if (!sjs::apps::ParseSource(*src_r_s, &src_r) || src_r == sjs::realdata::Source::Unknown) {
    SJS_LOG_ERROR("Unknown --src_r=", std::string(*src_r_s));
    return 2;
  }
  if (!sjs::apps::ParseSource(*src_s_s, &src_s) || src_s == sjs::realdata::Source::Unknown) {
    SJS_LOG_ERROR("Unknown --src_s=", std::string(*src_s_s));
    return 2;
  }

  // Ensure output dirs.
  if (!sjs::apps::EnsureDir(fs::path(std::string(*out_r)).parent_path()) ||
      !sjs::apps::EnsureDir(fs::path(std::string(*out_s)).parent_path())) {
    SJS_LOG_ERROR("Failed to create output directories.");
    return 3;
  }

  // Options.
  sjs::realdata::LoadOptions opt_r;
  sjs::realdata::LoadOptions opt_s;
  opt_r.drop_empty = sjs::apps::ParseBool(args.Get("drop_empty"), /*def=*/true);
  opt_s.drop_empty = opt_r.drop_empty;
  opt_r.to_mbr = sjs::apps::ParseBool(args.Get("to_mbr"), /*def=*/true);
  opt_s.to_mbr = opt_r.to_mbr;

  sjs::apps::ParseU64(args.Get("limit_r"), &opt_r.limit);
  sjs::apps::ParseU64(args.Get("limit_s"), &opt_s.limit);

  // Load.
  sjs::Dataset<2, sjs::Scalar> ds;
  ds.name = args.Get("dataset") ? std::string(*args.Get("dataset")) : "realdata";

  std::string err;
  if (!sjs::realdata::LoadDatasetPair<2, sjs::Scalar>(
          src_r, std::string(*in_r), src_s, std::string(*in_s),
          &ds, opt_r, opt_s, &err)) {
    SJS_LOG_ERROR("Realdata loading failed: ", err);
    SJS_LOG_ERROR("(Expected if you have not implemented the realdata loaders yet.)");
    return 4;
  }

  // Validate.
  {
    std::string v;
    if (!ds.Validate(/*require_proper=*/true, &v)) {
      SJS_LOG_ERROR("Loaded dataset failed validation: ", v);
      return 5;
    }
  }

  // Write binary.
  {
    sjs::binary::BinaryWriteOptions opt;
    opt.half_open = true;
    opt.write_ids = true;
    opt.scalar = sjs::binary::ScalarEncoding::Float64;
    opt.write_name = true;

    if (!sjs::binary::WriteDatasetBinaryPair<2, sjs::Scalar>(
            std::string(*out_r), std::string(*out_s), ds, opt, &err)) {
      SJS_LOG_ERROR("Binary write failed: ", err);
      return 6;
    }
  }

  SJS_LOG_INFO("Converted real data OK. out_r=", std::string(*out_r), " out_s=", std::string(*out_s));
  return 0;
}
