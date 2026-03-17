#pragma once
// join/sweep_events.h
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
//   3) object id ascending (fixed total order; SJS v3 §1.3.1)
//   4) side tie-break (only if ids tie)
//   5) index (final deterministic tie-break)

#include "core/assert.h"
#include "io/dataset.h"
#include "join/join_types.h"

#include <algorithm>
#include <cstring>
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
    events->push_back(Event{start, EventKind::Start, side, id, static_cast<u32>(i)});
    events->push_back(Event{end, EventKind::End, side, id, static_cast<u32>(i)});
  }
}


// Sort events in-place (deterministic total order).
//
// Stage3 (Opt): Radix sort for large event arrays
// ----------------------------------------------
// The comparator-based std::sort becomes a major build-time hotspot for very
// large datasets (tens of millions of events). We provide a stable multi-key
// radix sort that exactly matches EventLess ordering:
//   (x, kind, id, side, index)
// where kind is End(0) before Start(1), and side order depends on SideTieBreak.
//
// Notes:
// - We normalize -0.0 to +0.0 to match the comparator's "equal" treatment.
// - Inputs are assumed finite (no NaN). If NaN is observed, we fall back to
//   std::sort to preserve strict-weak-ordering requirements.

namespace detail {

inline u64 OrderedKeyFromDouble(Scalar x) noexcept {
  // IEEE-754 double -> ordered u64 key.
  // For ascending order:
  //   - negatives are bitwise inverted
  //   - non-negatives flip sign bit
  u64 bits = 0;
  static_assert(sizeof(bits) == sizeof(x), "u64/double size mismatch");
  std::memcpy(&bits, &x, sizeof(bits));

  // Normalize -0.0 to +0.0 so it ties with +0.0 under the key.
  if (bits == 0x8000000000000000ULL) bits = 0ULL;

  const u64 sign = bits >> 63;
  if (sign) return ~bits;
  return bits ^ 0x8000000000000000ULL;
}

inline u8 SideKey(Side s, SideTieBreak order) noexcept {
  const u8 v = static_cast<u8>(s);  // R=0, S=1
  if (order == SideTieBreak::RBeforeS) return v;
  return static_cast<u8>(1U - v);
}

// Stable counting sort on 16-bit digit.
template <class KeyFn>
inline void CountingPass16(const std::vector<Event>& src,
                           std::vector<Event>* dst,
                           KeyFn key_fn,
                           int shift) {
  static constexpr usize kBuckets = 1ULL << 16;
  static thread_local std::vector<usize> counts;
  if (counts.size() != kBuckets) counts.assign(kBuckets, 0);
  std::fill(counts.begin(), counts.end(), 0);

  for (const auto& e : src) {
    const u64 key = static_cast<u64>(key_fn(e));
    const usize b = static_cast<usize>((key >> shift) & 0xFFFFULL);
    ++counts[b];
  }

  usize sum = 0;
  for (usize i = 0; i < kBuckets; ++i) {
    const usize c = counts[i];
    counts[i] = sum;
    sum += c;
  }

  dst->resize(src.size());
  for (const auto& e : src) {
    const u64 key = static_cast<u64>(key_fn(e));
    const usize b = static_cast<usize>((key >> shift) & 0xFFFFULL);
    (*dst)[counts[b]++] = e;
  }
}

// Stable counting sort on 2 buckets (0/1 keys).
template <class KeyFn>
inline void CountingPass2(const std::vector<Event>& src,
                          std::vector<Event>* dst,
                          KeyFn key_fn) {
  usize c0 = 0;
  usize c1 = 0;
  for (const auto& e : src) {
    const u8 k = static_cast<u8>(key_fn(e));
    if (k == 0) ++c0;
    else ++c1;
  }

  dst->resize(src.size());
  usize p0 = 0;
  usize p1 = c0;
  for (const auto& e : src) {
    const u8 k = static_cast<u8>(key_fn(e));
    if (k == 0) (*dst)[p0++] = e;
    else (*dst)[p1++] = e;
  }
}

inline bool HasNaN(const std::vector<Event>& ev) {
  for (const auto& e : ev) {
    // NaN != NaN.
    if (!(e.x == e.x)) return true;
  }
  return false;
}

inline void SortSweepEventsRadix(std::vector<Event>* events, SideTieBreak side_order) {
  SJS_ASSERT(events != nullptr);
  if (events->size() <= 1) return;

  // Defensive: if NaNs exist, comparator order is not strict; fall back.
  if (HasNaN(*events)) {
    EventLess less{side_order};
    std::sort(events->begin(), events->end(), less);
    return;
  }

  std::vector<Event> tmp;
  tmp.reserve(events->size());

  std::vector<Event>* src = events;
  std::vector<Event>* dst = &tmp;

  auto swap_buffers = [&]() {
    std::swap(src, dst);
  };

  // Lowest-priority key first (LSD): index, side, id, kind, x.
  // index (u32): 2 passes of 16-bit.
  CountingPass16(*src, dst, [](const Event& e) -> u64 { return static_cast<u64>(e.index); }, 0);
  swap_buffers();
  CountingPass16(*src, dst, [](const Event& e) -> u64 { return static_cast<u64>(e.index); }, 16);
  swap_buffers();

  // side (0/1)
  CountingPass2(*src, dst, [side_order](const Event& e) -> u8 { return SideKey(e.side, side_order); });
  swap_buffers();

  // id (u32): 2 passes
  CountingPass16(*src, dst, [](const Event& e) -> u64 { return static_cast<u64>(e.id); }, 0);
  swap_buffers();
  CountingPass16(*src, dst, [](const Event& e) -> u64 { return static_cast<u64>(e.id); }, 16);
  swap_buffers();

  // kind (0/1)
  CountingPass2(*src, dst, [](const Event& e) -> u8 { return static_cast<u8>(e.kind); });
  swap_buffers();

  // x (u64): 4 passes of 16-bit
  CountingPass16(*src, dst, [](const Event& e) -> u64 { return OrderedKeyFromDouble(e.x); }, 0);
  swap_buffers();
  CountingPass16(*src, dst, [](const Event& e) -> u64 { return OrderedKeyFromDouble(e.x); }, 16);
  swap_buffers();
  CountingPass16(*src, dst, [](const Event& e) -> u64 { return OrderedKeyFromDouble(e.x); }, 32);
  swap_buffers();
  CountingPass16(*src, dst, [](const Event& e) -> u64 { return OrderedKeyFromDouble(e.x); }, 48);
  swap_buffers();

  // Result might be in tmp buffer.
  if (src != events) {
    events->swap(*src);
  }
}

}  // namespace detail

inline void SortSweepEvents(std::vector<Event>* events,
                           SideTieBreak side_order = SideTieBreak::RBeforeS) {
  if (!events) return;
  if (events->size() <= 1) return;

  // Heuristic: radix wins only for sufficiently large arrays.
  // (For small sizes, std::sort is usually faster due to lower overhead.)
  static constexpr usize kRadixThreshold = 1ULL << 18;  // ~262k

  if (events->size() >= kRadixThreshold) {
    detail::SortSweepEventsRadix(events, side_order);
    return;
  }

  EventLess less{side_order};
  std::sort(events->begin(), events->end(), less);
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
