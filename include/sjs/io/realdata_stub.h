#pragma once
// sjs/io/realdata_stub.h
//
// Stub interfaces for loading real-world spatial datasets (OSM / TIGER / etc.).
//
// This project intentionally avoids heavy dependencies in the core codebase.
// When you are ready to run on real datasets, you can implement these functions
// using libraries such as:
//  - GDAL/OGR (for Shapefile/GeoJSON/Geopackage etc.)
//  - libosmium / osmium-tool (for OSM PBF)
//  - custom parsers
//
// The intent is to keep the rest of the system stable: baselines and runners
// should only depend on Relation<Dim> / Dataset<Dim> as defined in dataset.h.
// Real data loaders can live behind these thin functions.
//
// For now, these functions return false with a clear error message.

#include "sjs/core/types.h"
#include "sjs/io/dataset.h"

#include <string>
#include <string_view>

namespace sjs {
namespace realdata {

// Identifies a real dataset source type (so config/CLI can select importer).
enum class Source : u8 {
  OSM_PBF = 0,     // .osm.pbf (needs libosmium)
  TIGER_SHP = 1,   // ESRI Shapefile (needs GDAL/OGR or shapefile parser)
  GEOJSON = 2,     // GeoJSON/JSON (needs JSON parser)
  WKT_CSV = 3,     // CSV with WKT geometry (needs WKT parser)
  Unknown = 255,
};

inline constexpr std::string_view ToString(Source s) noexcept {
  switch (s) {
    case Source::OSM_PBF: return "osm_pbf";
    case Source::TIGER_SHP: return "tiger_shp";
    case Source::GEOJSON: return "geojson";
    case Source::WKT_CSV: return "wkt_csv";
    case Source::Unknown: return "unknown";
  }
  return "unknown";
}

// Options common to many real dataset loaders.
struct LoadOptions {
  // If true, compute MBRs (bounding boxes) from input geometries (polygons/lines).
  // Many file formats store non-rectangular geometries.
  bool to_mbr = true;

  // Filter empty MBRs.
  bool drop_empty = true;

  // Optional spatial clipping window (only keep objects whose MBR intersects clip).
  // If clip is empty, ignore clipping.
  // (Dim is 2 for spatial; for spatio-temporal you can extend.)
  Box<2, Scalar> clip = Box<2, Scalar>::Empty();

  // Limit number of features (0 = no limit).
  u64 limit = 0;
};

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

// Generic entry point: load a relation as MBRs.
// For now, it's a stub; later you can dispatch by Source to concrete loaders.
template <int Dim, class T>
bool LoadRelation(Source /*src*/,
                  const std::string& /*path*/,
                  Relation<Dim, T>* /*out_rel*/,
                  const LoadOptions& /*opt*/ = {},
                  std::string* err = nullptr) {
  SetErr(err, "realdata::LoadRelation is not implemented yet. "
               "Integrate GDAL/OGR or libosmium and implement this function.");
  return false;
}

// Generic entry point: load a dataset pair (R,S) from two files.
template <int Dim, class T>
bool LoadDatasetPair(Source src_r,
                     const std::string& path_r,
                     Source src_s,
                     const std::string& path_s,
                     Dataset<Dim, T>* out_ds,
                     const LoadOptions& opt_r = {},
                     const LoadOptions& opt_s = {},
                     std::string* err = nullptr) {
  if (!out_ds) {
    SetErr(err, "out_ds is null");
    return false;
  }
  Dataset<Dim, T> ds;
  ds.half_open = true;
  if (!LoadRelation<Dim, T>(src_r, path_r, &ds.R, opt_r, err)) return false;
  if (!LoadRelation<Dim, T>(src_s, path_s, &ds.S, opt_s, err)) return false;
  *out_ds = std::move(ds);
  return true;
}

}  // namespace realdata
}  // namespace sjs
