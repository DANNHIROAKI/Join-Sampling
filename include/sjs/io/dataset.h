#pragma once
// sjs/io/dataset.h
//
// Dataset container types for Spatial Join Sampling experiments.
//
// Design goals:
//  - Minimal dependencies (core + geometry only).
//  - Support both 2D and future higher-D (templated by Dim).
//  - Keep object IDs stable across reordering (optional ids vector).
//  - Provide common metadata (name, domain bounds, etc.) for experiment logging.
//
// Notes:
//  - Geometry uses half-open boxes [lo, hi) in each dimension.
//  - Many algorithms will reorder boxes internally; call EnsureIds() to
//    preserve stable IDs (index becomes unstable after permutation).
//
// This header defines:
//  - Relation<Dim>: a collection of boxes with optional explicit IDs.
//  - Dataset<Dim>: a pair of relations R and S plus metadata.
//  - Lightweight validation helpers.

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/embedding.h"

#include <algorithm>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {

// A relation is a set of axis-aligned boxes (MBRs), plus optional stable ids.
// If ids is empty, id(i) is implicitly i (0-based).
template <int Dim, class T = Scalar>
struct Relation {
  using BoxT = Box<Dim, T>;

  std::string name;               // optional; useful for logs
  std::vector<BoxT> boxes;        // geometry
  std::vector<Id> ids;            // optional stable ids; may be empty

  void Clear() {
    name.clear();
    boxes.clear();
    ids.clear();
  }

  void Reserve(usize n) {
    boxes.reserve(n);
    ids.reserve(n);
  }

  usize Size() const noexcept { return boxes.size(); }
  bool Empty() const noexcept { return boxes.empty(); }

  bool HasExplicitIds() const noexcept { return !ids.empty(); }

  // Returns the stable ID of the i-th stored object.
  // If ids[] is empty, this returns i (narrowed to Id).
  Id GetId(usize i) const noexcept {
    if (ids.empty()) return static_cast<Id>(i);
    return ids[i];
  }

  // Ensure ids exist and are stable (0..n-1 by default).
  // Useful before any algorithm that permutes / partitions the boxes.
  void EnsureIds() {
    if (!ids.empty()) return;
    ids.resize(boxes.size());
    for (usize i = 0; i < boxes.size(); ++i) ids[i] = static_cast<Id>(i);
  }

  // Set ids to sequential 0..n-1 (overwriting any existing ids).
  void ForceSequentialIds() {
    ids.resize(boxes.size());
    for (usize i = 0; i < boxes.size(); ++i) ids[i] = static_cast<Id>(i);
  }

  // Add a box; if you pass an id, ids will become explicit.
  // If some objects have ids and some don't, we normalize by filling missing ones.
  void Add(const BoxT& b, std::optional<Id> id = std::nullopt) {
    boxes.push_back(b);
    if (id.has_value()) {
      if (ids.empty() && boxes.size() > 1) {
        // Previously implicit; materialize old implicit ids.
        ids.resize(boxes.size() - 1);
        for (usize i = 0; i + 1 < boxes.size(); ++i) ids[i] = static_cast<Id>(i);
      }
      ids.push_back(*id);
    } else {
      if (!ids.empty()) {
        // Previously explicit; assign a default id for this new element.
        ids.push_back(static_cast<Id>(boxes.size() - 1));
      }
    }
  }

  // Bounding box of all boxes in the relation. Returns Box::Empty() if no objects.
  BoxT Bounds() const noexcept {
    BoxT b = BoxT::Empty();
    for (const auto& x : boxes) b.ExpandToIncludeBox(x);
    return b;
  }

  // Remove empty boxes (lo >= hi in any dimension).
  // If ids are present, they are removed correspondingly.
  void RemoveEmptyBoxes() {
    if (boxes.empty()) return;
    const bool has_ids = !ids.empty();
    usize w = 0;
    for (usize i = 0; i < boxes.size(); ++i) {
      if (!boxes[i].IsEmpty()) {
        if (w != i) {
          boxes[w] = boxes[i];
          if (has_ids) ids[w] = ids[i];
        }
        ++w;
      }
    }
    boxes.resize(w);
    if (has_ids) ids.resize(w);
  }

  // Basic validation:
  //  - if require_proper: each box must satisfy lo < hi in all dimensions.
  //  - otherwise we only require lo <= hi (valid but possibly empty).
  bool Validate(bool require_proper = true, std::string* err = nullptr) const {
    for (usize i = 0; i < boxes.size(); ++i) {
      const auto& b = boxes[i];
      if (!b.IsValid()) {
        if (err) *err = "Invalid box (lo > hi) at index " + std::to_string(i);
        return false;
      }
      if (require_proper && b.IsEmpty()) {
        if (err) *err = "Empty box (lo >= hi) at index " + std::to_string(i);
        return false;
      }
    }
    if (!ids.empty() && ids.size() != boxes.size()) {
      if (err) *err = "ids.size() != boxes.size()";
      return false;
    }
    return true;
  }
};

// Dataset: a pair of relations R and S plus common metadata.
template <int Dim, class T = Scalar>
struct Dataset {
  using BoxT = Box<Dim, T>;
  using RelT = Relation<Dim, T>;

  std::string name;   // dataset name / tag (e.g., "OSM_buildings_vs_roads")
  bool half_open = true;  // geometry convention (should remain true in this project)

  RelT R;
  RelT S;

  // Cached global domain bounds (union of R and S). Optional; computed on demand.
  mutable bool domain_cached = false;
  mutable BoxT domain_cache = BoxT::Empty();

  void Clear() {
    name.clear();
    half_open = true;
    R.Clear();
    S.Clear();
    domain_cached = false;
    domain_cache = BoxT::Empty();
  }

  usize SizeR() const noexcept { return R.Size(); }
  usize SizeS() const noexcept { return S.Size(); }
  usize TotalSize() const noexcept { return R.Size() + S.Size(); }

  // Compute and cache the union bounds across both relations.
  BoxT Domain() const noexcept {
    if (!domain_cached) {
      BoxT d = BoxT::Empty();
      d.ExpandToIncludeBox(R.Bounds());
      d.ExpandToIncludeBox(S.Bounds());
      domain_cache = d;
      domain_cached = true;
    }
    return domain_cache;
  }

  // Domain bounds for embedding-based baselines (KD/RangeTree/SIRS).
  DomainBounds<Dim, T> DomainForEmbedding() const noexcept {
    DomainBounds<Dim, T> b;
    for (const auto& x : R.boxes) b.ExpandToInclude(x);
    for (const auto& x : S.boxes) b.ExpandToInclude(x);
    if (!b.IsInitialized()) {
      b.min_lo = Point<Dim, T>::Zero();
      b.max_hi = Point<Dim, T>::Zero();
    }
    return b;
  }

  // Ensure stable IDs exist for both relations.
  void EnsureIds() {
    R.EnsureIds();
    S.EnsureIds();
  }

  // Remove empty boxes in both relations.
  void RemoveEmptyBoxes() {
    R.RemoveEmptyBoxes();
    S.RemoveEmptyBoxes();
    domain_cached = false;
  }

  bool Validate(bool require_proper = true, std::string* err = nullptr) const {
    if (!half_open) {
      if (err) *err = "Dataset is not marked as half-open";
      return false;
    }
    std::string tmp;
    if (!R.Validate(require_proper, &tmp)) {
      if (err) *err = "R: " + tmp;
      return false;
    }
    if (!S.Validate(require_proper, &tmp)) {
      if (err) *err = "S: " + tmp;
      return false;
    }
    return true;
  }
};

}  // namespace sjs
