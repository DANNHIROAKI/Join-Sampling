#pragma once
// baselines/ours/detail/slot_plan_nd.h
//
// Generic event/mode slot planner for the ND SJS sampling path.

#include "baselines/ours/detail/mode_mask.h"
#include "sampling/alias_table.h"

#include <algorithm>
#include <array>
#include <limits>
#include <string>
#include <vector>

namespace sjs {
namespace baselines {
namespace ours {
namespace detail {

template <int Dim>
struct SlotPlanND {
  static constexpr u32 kModeCount = kModeCountV<Dim>;

  u32 num_events = 0;
  std::array<std::vector<u32>, kModeCount> slots_by_mode;
  std::array<std::vector<u32>, kModeCount> offsets_by_mode;

  void Clear() {
    num_events = 0;
    for (auto& v : slots_by_mode) v.clear();
    for (auto& v : offsets_by_mode) v.clear();
  }
};

template <int Dim>
inline bool BuildSlotPlanND(u32 t,
                            Rng* rng,
                            const std::vector<ModeWeights<Dim>>& w_modes,
                            SlotPlanND<Dim>* out,
                            std::string* err) {
  static constexpr u32 kModeCount = kModeCountV<Dim>;
  if (!out || !rng) {
    if (err) *err = "BuildSlotPlanND: null out/rng";
    return false;
  }
  out->Clear();
  out->num_events = static_cast<u32>(w_modes.size());
  if (t == 0) {
    for (auto& off : out->offsets_by_mode) off.assign(static_cast<usize>(out->num_events) + 1U, 0U);
    return true;
  }
  if (w_modes.empty()) {
    if (err) *err = "BuildSlotPlanND: empty event set";
    return false;
  }

  std::vector<u64> event_weights(w_modes.size(), 0ULL);
  for (usize e = 0; e < w_modes.size(); ++e) {
    u64 total = 0ULL;
    for (u32 m = 0; m < kModeCount; ++m) total += w_modes[e][m];
    event_weights[e] = total;
  }

  sampling::AliasTable event_alias;
  if (!event_alias.BuildFromU64(Span<const u64>(event_weights), err)) {
    if (err && err->empty()) *err = "BuildSlotPlanND: failed to build event alias";
    return false;
  }

  std::vector<u32> event_of_slot(static_cast<usize>(t), 0U);
  std::vector<u32> mode_of_slot(static_cast<usize>(t), 0U);
  std::array<std::vector<u32>, kModeCount> counts_by_mode;
  for (auto& c : counts_by_mode) c.assign(w_modes.size(), 0U);

  for (u32 slot = 0; slot < t; ++slot) {
    const u32 e = static_cast<u32>(event_alias.Sample(rng));
    event_of_slot[static_cast<usize>(slot)] = e;

    const u64 total = event_weights[static_cast<usize>(e)];
    if (total == 0ULL) {
      if (err) *err = "BuildSlotPlanND: sampled zero-weight event";
      return false;
    }

    const u64 x = rng->UniformU64(total);
    u64 cum = 0ULL;
    u32 chosen = 0;
    for (; chosen < kModeCount; ++chosen) {
      cum += w_modes[static_cast<usize>(e)][chosen];
      if (x < cum) break;
    }
    if (chosen >= kModeCount) chosen = kModeCount - 1;
    mode_of_slot[static_cast<usize>(slot)] = chosen;
    ++counts_by_mode[chosen][static_cast<usize>(e)];
  }

  for (u32 m = 0; m < kModeCount; ++m) {
    auto& offsets = out->offsets_by_mode[m];
    offsets.assign(w_modes.size() + 1U, 0U);
    for (usize e = 0; e < w_modes.size(); ++e) {
      offsets[e + 1] = offsets[e] + counts_by_mode[m][e];
    }
    out->slots_by_mode[m].assign(static_cast<usize>(offsets.back()), 0U);
  }

  auto write_ptr = out->offsets_by_mode;
  for (u32 slot = 0; slot < t; ++slot) {
    const u32 e = event_of_slot[static_cast<usize>(slot)];
    const u32 m = mode_of_slot[static_cast<usize>(slot)];
    const u32 pos = write_ptr[m][static_cast<usize>(e)]++;
    out->slots_by_mode[m][static_cast<usize>(pos)] = slot;
  }

  return true;
}

}  // namespace detail
}  // namespace ours
}  // namespace baselines
}  // namespace sjs
