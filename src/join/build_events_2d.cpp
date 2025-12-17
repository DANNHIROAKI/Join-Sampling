// src/join/build_events_2d.cpp
//
// Dim=2 wrappers for sweep event construction.
//
// Motivation
// ----------
// The library provides header-only templates for event building
//   (see include/sjs/join/sweep_events.h).
// For apps that operate only on 2D datasets, it is convenient to have small
// non-templated wrappers to reduce compile times and avoid template plumbing.
//
// These wrappers are intentionally minimal and simply forward to the template
// implementation.

#include "sjs/join/sweep_events.h"

#include "sjs/core/assert.h"
#include "sjs/core/types.h"
#include "sjs/io/dataset.h"
#include "sjs/join/join_types.h"

#include <vector>

namespace sjs {
namespace join {

std::vector<Event> BuildSweepEvents2D(const Relation<2, Scalar>& R,
                                     const Relation<2, Scalar>& S,
                                     int axis,
                                     SideTieBreak side_order) {
  if (axis < 0 || axis >= 2) {
    // Defensive: keep behavior predictable for tooling.
    axis = 0;
  }
  return BuildSweepEvents<2, Scalar>(R, S, axis, side_order);
}

std::vector<Event> BuildSweepEvents2D(const Dataset<2, Scalar>& ds,
                                     int axis,
                                     SideTieBreak side_order) {
  if (axis < 0 || axis >= 2) axis = 0;
  return BuildSweepEvents<2, Scalar>(ds, axis, side_order);
}

void SortSweepEvents2D(std::vector<Event>* events, SideTieBreak side_order) {
  SJS_ASSERT(events != nullptr);
  SortSweepEvents(events, side_order);
}

}  // namespace join
}  // namespace sjs
