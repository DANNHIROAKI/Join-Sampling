// src/io/realdata_stub.cpp
//
// Real-data import stubs (Dim=2).
//
// The project is currently focused on synthetic datasets. We still keep a
// placeholder API for importing real-world spatial data (OSM, TIGER, etc.) so
// that future extensions can plug in without refactoring the experiment runner.
//
// This file provides non-template convenience wrappers for Dim=2 that forward
// to the header-only stubs in include/sjs/io/realdata_stub.h.

#include "sjs/core/types.h"
#include "sjs/io/dataset.h"
#include "sjs/io/realdata_stub.h"

#include <string>

namespace sjs {
namespace realdata {

bool LoadRelation2D(Source src,
                    const std::string& path,
                    Relation<2, Scalar>* out_rel,
                    const LoadOptions& opt,
                    std::string* err) {
  return LoadRelation<2, Scalar>(src, path, out_rel, opt, err);
}

bool LoadDatasetPair2D(Source src_r,
                       const std::string& path_r,
                       Source src_s,
                       const std::string& path_s,
                       Dataset<2, Scalar>* out_ds,
                       const LoadOptions& opt,
                       std::string* err) {
  // The template API supports separate options for R and S.
  // For now we pass the same options to both sides.
  return LoadDatasetPair<2, Scalar>(src_r, path_r, src_s, path_s, out_ds, opt, opt, err);
}

}  // namespace realdata
}  // namespace sjs
