#pragma once
// join/join_types.h
//
// Public definitions shared by join algorithms / baselines.
//
// Conventions
// -----------
// Geometry uses half-open boxes: [lo, hi) on each axis.
// For a plane sweep on axis a:
//   - START event at lo[a]
//   - END   event at hi[a]
//   - END events are processed before START events at the same coordinate
// This ordering matches half-open semantics: objects ending at x are not
// considered active for objects starting at x.

#include "core/types.h"

#include <cstdint>
#include <ostream>
#include <string>
#include <string_view>

namespace sjs {
namespace join {

// Which relation a box/event belongs to.
enum class Side : u8 {
  R = 0,
  S = 1,
};

inline constexpr Side OtherSide(Side s) noexcept {
  return (s == Side::R) ? Side::S : Side::R;
}

inline constexpr std::string_view ToString(Side s) noexcept {
  switch (s) {
    case Side::R: return "R";
    case Side::S: return "S";
  }
  return "?";
}

// Event kind for plane sweep.
// NOTE: numeric order is chosen so that sorting by kind ascending yields
// END before START (to respect half-open boxes).
enum class EventKind : u8 {
  End = 0,
  Start = 1,
};

inline constexpr std::string_view ToString(EventKind k) noexcept {
  switch (k) {
    case EventKind::End: return "end";
    case EventKind::Start: return "start";
  }
  return "?";
}

// Tie-break preference when multiple events share the same coordinate and kind.
// This affects determinism and (if you later decompose join pairs per event)
// affects which START event is considered "later" when lo coordinates tie.
enum class SideTieBreak : u8 {
  RBeforeS = 0,
  SBeforeR = 1,
};

inline constexpr bool IsFirst(SideTieBreak tb, Side a, Side b) noexcept {
  if (a == b) return false;
  if (tb == SideTieBreak::RBeforeS) {
    return a == Side::R && b == Side::S;
  }
  return a == Side::S && b == Side::R;
}

inline constexpr std::string_view ToString(SideTieBreak tb) noexcept {
  switch (tb) {
    case SideTieBreak::RBeforeS: return "R_before_S";
    case SideTieBreak::SBeforeR: return "S_before_R";
  }
  return "?";
}

// A sweep event associated with a specific box in R or S.
// - x: coordinate on the sweep axis.
// - kind: Start or End.
// - side: R or S.
// - id: stable object ID (from Relation::GetId()).
// - index: index into Relation::boxes.
struct Event {
  Scalar x{};
  EventKind kind{EventKind::Start};
  Side side{Side::R};
  Id id{kInvalidId};
  u32 index{0};
};

// Comparator implementing deterministic sweep ordering.
//
// SJS v3 alignment note
// --------------------
// To match the event-order requirement for half-open semantics (END before START)
// *and* the "fixed total order" tie-break described in SJS v3 (§1.3.1), we
// break ties by object id before breaking ties by side.
//
// Order:
//   1) x ascending
//   2) kind ascending (End before Start)
//   3) id ascending (global tie-break)
//   4) side according to SideTieBreak (only if ids tie)
//   5) index (final deterministic tie-break)
struct EventLess {
  SideTieBreak side_order{SideTieBreak::RBeforeS};

  bool operator()(const Event& a, const Event& b) const noexcept {
    if (a.x < b.x) return true;
    if (b.x < a.x) return false;

    if (static_cast<u8>(a.kind) < static_cast<u8>(b.kind)) return true;
    if (static_cast<u8>(b.kind) < static_cast<u8>(a.kind)) return false;

    if (a.id < b.id) return true;
    if (b.id < a.id) return false;

    if (a.side != b.side) {
      return IsFirst(side_order, a.side, b.side);
    }

    return a.index < b.index;
  }
};

// Lightweight stats for join enumeration (plane sweep or naive oracle).
struct JoinStats {
  u64 num_events = 0;
  u64 candidate_checks = 0;  // number of candidate intersection predicate evaluations
  u64 output_pairs = 0;      // number of reported intersecting (R,S) pairs

  u64 active_max_r = 0;
  u64 active_max_s = 0;

  void Reset() {
    num_events = 0;
    candidate_checks = 0;
    output_pairs = 0;
    active_max_r = 0;
    active_max_s = 0;
  }

  std::string ToJsonLite() const {
    return std::string("{") +
           "\"num_events\":" + std::to_string(num_events) + "," +
           "\"candidate_checks\":" + std::to_string(candidate_checks) + "," +
           "\"output_pairs\":" + std::to_string(output_pairs) + "," +
           "\"active_max_r\":" + std::to_string(active_max_r) + "," +
           "\"active_max_s\":" + std::to_string(active_max_s) +
           "}";
  }
};

// Optional abstraction: a deterministic stream of join pairs.
// This is useful for two-pass algorithms (e.g., enumerate + rank sampling).
class IJoinStream {
 public:
  virtual ~IJoinStream() = default;

  // Reset to stream start.
  virtual void Reset() = 0;

  // Produce next intersecting pair. Returns false at end.
  virtual bool Next(PairId* out) = 0;
};

inline std::ostream& operator<<(std::ostream& os, Side s) {
  os << ToString(s);
  return os;
}

inline std::ostream& operator<<(std::ostream& os, EventKind k) {
  os << ToString(k);
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const Event& e) {
  os << "Event{x=" << e.x << ", kind=" << ToString(e.kind) << ", side=" << ToString(e.side)
     << ", id=" << e.id << ", idx=" << e.index << "}";
  return os;
}

}  // namespace join
}  // namespace sjs
