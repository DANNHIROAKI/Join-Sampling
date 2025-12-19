#pragma once
// sjs/join/join_enumerator.h
//
// Generic join enumerator (plane sweep) with deterministic ordering.
//
// This is intended as a reusable component for:
//  - Enum+Sampling variants (two-pass algorithms that need a stable join stream)
//  - Adaptive algorithms that interleave enumeration and sampling
//  - Large correctness checks that are faster than O(|R|*|S|)
//
// Implementation notes
// --------------------
// We use a standard plane sweep on a chosen axis. We generate START/END events
// for both relations, sort them with half-open tie-break (END before START at the same coordinate), and maintain two active sets.
// When processing a START event for one side, we test against all active boxes in the other side and report intersecting pairs.
//
// This implementation is a generic enumerator. It is not meant to beat highly optimized join baselines (PBSM, R-tree join, etc.), but it provides:
//   - deterministic join stream order,
//   - streaming enumeration (no need to store all pairs),
//   - simple stats instrumentation.

#include "sjs/core/assert.h"
#include "sjs/core/types.h"
#include "sjs/geometry/predicates.h"
#include "sjs/io/dataset.h"
#include "sjs/join/join_types.h"
#include "sjs/join/sweep_events.h"

#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace join {

// Options for the plane-sweep enumerator.
struct PlaneSweepOptions {
  // Sweep axis. For 2D datasets, axis=0 is x and axis=1 is y.
  int axis = 0;

  // Deterministic tie-break between R/S when x and kind tie.
  SideTieBreak side_order = SideTieBreak::RBeforeS;

  // If true, intersection checks skip the sweep axis because active-set
  // membership implies overlap on the sweep axis for proper (non-empty) boxes.
  // Set false if you want the most conservative check.
  bool skip_axis_check = true;
};

namespace detail {

// Half-open interval overlap on one axis: [a0,a1) intersects [b0,b1) iff a0 < b1 && b0 < a1.
template <class T>
inline bool Intersects1D(T a0, T a1, T b0, T b1) noexcept {
  return (a0 < b1) && (b0 < a1);
}

// Check intersection between two boxes but optionally skip one axis.
// Returns true iff they intersect under half-open semantics.
template <int Dim, class T>
inline bool IntersectsBox(const Box<Dim, T>& a,
                          const Box<Dim, T>& b,
                          int skip_axis,
                          bool skip_enabled) noexcept {
  // Optional quick empty check.
  if (a.IsEmpty() || b.IsEmpty()) return false;

  for (int d = 0; d < Dim; ++d) {
    if (skip_enabled && d == skip_axis) continue;
    const usize d_idx = static_cast<usize>(d);
    if (!Intersects1D<T>(a.lo.v[d_idx], a.hi.v[d_idx], b.lo.v[d_idx], b.hi.v[d_idx])) return false;
  }
  return true;
}

// O(1) active-set container with swap-delete.
struct ActiveSet {
  std::vector<usize> items;
  std::vector<i64> pos;  // pos[idx] = position in items, -1 if inactive

  void Init(usize n) {
    items.clear();
    items.reserve(n);
    pos.assign(n, -1);
  }

  bool Contains(usize idx) const noexcept {
    return idx < pos.size() && pos[idx] >= 0;
  }

  void Insert(usize idx) {
    SJS_DASSERT(idx < pos.size());
    SJS_DASSERT(pos[idx] == -1);
    pos[idx] = static_cast<i64>(items.size());
    items.push_back(idx);
  }

  void Erase(usize idx) {
    SJS_DASSERT(idx < pos.size());
    const i64 p = pos[idx];
    SJS_DASSERT(p >= 0);
    const usize up = static_cast<usize>(p);
    const usize last = items.back();
    items[up] = last;
    pos[last] = p;
    items.pop_back();
    pos[idx] = -1;
  }

  usize Size() const noexcept { return items.size(); }
};

}  // namespace detail

// A deterministic stream of join pairs produced by plane sweep.
//
// Stream order:
//  - Event order: BuildSweepEvents + EventLess (axis + half-open + side tie-break).
//  - Within a START event: iterate active opposite set in the internal order of
//    the active container (swap-delete). This order is deterministic for a given
//    dataset and side_order.
//
// The stream does NOT store all pairs; it yields one pair at a time.
//
// NOTE: The stream processes all pairs for a START event *before* inserting that
// box into its active set. This ensures that when two boxes start at the same
// coordinate, the pair is reported by the "later" START event (as defined by side_order).

template <int Dim, class T = Scalar>
class PlaneSweepJoinStream final : public IJoinStream {
 public:
  PlaneSweepJoinStream(const Relation<Dim, T>& R,
                       const Relation<Dim, T>& S,
                       PlaneSweepOptions opt = {})
      : R_(&R), S_(&S), opt_(opt) {
    SJS_ASSERT(opt_.axis >= 0 && opt_.axis < Dim);
    Reset();
  }

  void Reset() override {
    events_ = BuildSweepEvents<Dim, T>(*R_, *S_, opt_.axis, opt_.side_order);
    stats_.Reset();
    stats_.num_events = static_cast<u64>(events_.size());

    active_r_.Init(R_->boxes.size());
    active_s_.Init(S_->boxes.size());

    ev_i_ = 0;

    scanning_ = false;
    scan_other_pos_ = 0;
    scan_start_side_ = Side::R;
    scan_start_idx_ = 0;
  }

  const JoinStats& Stats() const noexcept { return stats_; }

  bool Next(PairId* out) override {
    if (!out) return false;

    while (true) {
      // If we are in the middle of scanning candidates for a START event,
      // continue scanning until we find the next intersecting pair.
      if (scanning_) {
        const auto& start_box = (scan_start_side_ == Side::R)
                                    ? R_->boxes[scan_start_idx_]
                                    : S_->boxes[scan_start_idx_];
        const Id start_id = (scan_start_side_ == Side::R)
                                ? R_->GetId(scan_start_idx_)
                                : S_->GetId(scan_start_idx_);

        detail::ActiveSet& other_as = (scan_start_side_ == Side::R) ? active_s_ : active_r_;
        const auto& other_boxes = (scan_start_side_ == Side::R) ? S_->boxes : R_->boxes;
        const auto& other_rel = (scan_start_side_ == Side::R) ? *S_ : *R_;

        while (scan_other_pos_ < other_as.items.size()) {
          const usize other_idx = other_as.items[scan_other_pos_++];
          stats_.candidate_checks++;
          const auto& other_box = other_boxes[other_idx];

          const bool ok = detail::IntersectsBox<Dim, T>(start_box, other_box, opt_.axis, opt_.skip_axis_check);
          if (!ok) continue;

          // Emit pair in (R,S) order.
          if (scan_start_side_ == Side::R) {
            *out = PairId{start_id, other_rel.GetId(other_idx)};  // start is R, other is S
          } else {
            *out = PairId{other_rel.GetId(other_idx), start_id};  // start is S, other is R
          }
          stats_.output_pairs++;
          return true;
        }

        // Done scanning this START event: now insert the start box into its active set.
        if (scan_start_side_ == Side::R) {
          active_r_.Insert(scan_start_idx_);
          if (static_cast<u64>(active_r_.Size()) > stats_.active_max_r) stats_.active_max_r = active_r_.Size();
        } else {
          active_s_.Insert(scan_start_idx_);
          if (static_cast<u64>(active_s_.Size()) > stats_.active_max_s) stats_.active_max_s = active_s_.Size();
        }
        scanning_ = false;
        // Continue to next event.
        continue;
      }

      // No active scanning: advance events.
      if (ev_i_ >= events_.size()) {
        return false;
      }

      const Event e = events_[ev_i_++];
      if (e.kind == EventKind::End) {
        // Remove from active set.
        if (e.side == Side::R) {
          if (active_r_.Contains(e.index)) active_r_.Erase(e.index);
        } else {
          if (active_s_.Contains(e.index)) active_s_.Erase(e.index);
        }
        continue;
      }

      // START: begin scanning candidates from opposite active set.
      scanning_ = true;
      scan_other_pos_ = 0;
      scan_start_side_ = e.side;
      scan_start_idx_ = e.index;
      // We do NOT insert e.index yet; it will be inserted after scanning completes.
      continue;
    }
  }

 private:
  const Relation<Dim, T>* R_;
  const Relation<Dim, T>* S_;
  PlaneSweepOptions opt_;

  std::vector<Event> events_;
  usize ev_i_{0};

  detail::ActiveSet active_r_;
  detail::ActiveSet active_s_;

  // Streaming state for the current START event scanning.
  bool scanning_{false};
  usize scan_other_pos_{0};
  Side scan_start_side_{Side::R};
  usize scan_start_idx_{0};

  JoinStats stats_;
};

// Convenience: enumerate join pairs via plane sweep and call a callback.
// Callback signature:
//   bool emit(PairId p)
// Return false from emit to stop early.
//
// Returns true if fully enumerated, false if stopped early.
template <int Dim, class T, class EmitFn>
inline bool EnumeratePlaneSweep(const Relation<Dim, T>& R,
                                const Relation<Dim, T>& S,
                                EmitFn&& emit,
                                const PlaneSweepOptions& opt = {},
                                JoinStats* stats = nullptr) {
  PlaneSweepJoinStream<Dim, T> stream(R, S, opt);
  PairId p;

  while (stream.Next(&p)) {
    if (!emit(p)) {
      if (stats) *stats = stream.Stats();
      return false;
    }
  }

  if (stats) *stats = stream.Stats();
  return true;
}

// Convenience: count join pairs via plane sweep.
template <int Dim, class T>
inline u64 CountPlaneSweep(const Relation<Dim, T>& R,
                           const Relation<Dim, T>& S,
                           const PlaneSweepOptions& opt = {},
                           JoinStats* stats = nullptr) {
  PlaneSweepJoinStream<Dim, T> stream(R, S, opt);
  u64 cnt = 0;
  PairId p;
  while (stream.Next(&p)) {
    ++cnt;
  }
  if (stats) *stats = stream.Stats();
  return cnt;
}

}  // namespace join
}  // namespace sjs
