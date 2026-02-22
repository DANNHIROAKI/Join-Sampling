#pragma once
// sjs/join/sweep_events.h
//
// Plane-sweep event generation and ordering.
//
// This module standardizes event creation for half-open boxes.
// Given a sweep axis `axis`:
//   - Create START event at lo[axis]
//   - Create END event at hi[axis]
// Sort order (see EventLess):
//   1) coordinate ascending
//   2) END before START at the same coordinate (half-open)
//   3) object id ascending (fixed total order; SJS-HighDims §1.3.1)
//   4) side tie-break (only if ids tie)
//   5) index (final deterministic tie-break)

#include "sjs/core/assert.h"
#include "sjs/io/dataset.h"
#include "sjs/join/join_types.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace join {

// Append sweep events for one relation into `events`.
// - side: indicates whether these boxes are from R or S.
// - axis: sweep axis.
//
// Note: This function does not validate boxes. For correctness, your dataset should
// use proper boxes (lo < hi in every dimension). You can call Dataset::Validate().
template <int Dim, class T>
inline void AppendRelationEvents(const Relation<Dim, T>& rel,
                                Side side,
                                int axis,
                                std::vector<Event>* events) {
  SJS_ASSERT(events != nullptr);
  SJS_ASSERT(axis >= 0 && axis < Dim);

  events->reserve(events->size() + rel.boxes.size() * 2);

  for (usize i = 0; i < rel.boxes.size(); ++i) {
    const auto& b = rel.boxes[i];
    const Scalar start = static_cast<Scalar>(b.lo[axis]);
    const Scalar end = static_cast<Scalar>(b.hi[axis]);
    const Id id = rel.GetId(i);

    // Skip boxes that are empty on the sweep axis (start >= end).
    // Such boxes cannot intersect any proper box under half-open semantics.
    if (!(start < end)) continue;

    // Create both events; ordering is handled by sorting.
    events->push_back(Event{start, EventKind::Start, side, id, i});
    events->push_back(Event{end, EventKind::End, side, id, i});
  }
}

inline void SortSweepEvents(std::vector<Event>* events,
                           SideTieBreak side_order = SideTieBreak::RBeforeS) {
  SJS_ASSERT(events != nullptr);
  std::sort(events->begin(), events->end(), EventLess{side_order});
}

// Build and sort sweep events for a pair of relations.
template <int Dim, class T>
inline std::vector<Event> BuildSweepEvents(const Relation<Dim, T>& R,
                                           const Relation<Dim, T>& S,
                                           int axis = 0,
                                           SideTieBreak side_order = SideTieBreak::RBeforeS) {
  SJS_ASSERT(axis >= 0 && axis < Dim);
  std::vector<Event> events;
  events.reserve((R.boxes.size() + S.boxes.size()) * 2);
  AppendRelationEvents<Dim, T>(R, Side::R, axis, &events);
  AppendRelationEvents<Dim, T>(S, Side::S, axis, &events);
  SortSweepEvents(&events, side_order);
  return events;
}

// Convenience: build events for a Dataset.
template <int Dim, class T>
inline std::vector<Event> BuildSweepEvents(const Dataset<Dim, T>& ds,
                                           int axis = 0,
                                           SideTieBreak side_order = SideTieBreak::RBeforeS) {
  return BuildSweepEvents<Dim, T>(ds.R, ds.S, axis, side_order);
}

}  // namespace join
}  // namespace sjs
