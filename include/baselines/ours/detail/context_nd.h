#pragma once
// baselines/ours/detail/context_nd.h
//
// Shared ND preprocessing context for the 3D/4D/5D SJS extension.
//
// It keeps the common sweep-event cache and one mode family per relation side.
// Each mode family contains the fixed-mode local recursive structures
// \mathcal D_{Dim-1,g} described in the research outline.

#include "baselines/ours/detail/local_mode_index.h"
#include "core/assert.h"
#include "core/types.h"
#include "core/timer.h"
#include "io/dataset.h"
#include "join/sweep_events.h"

#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace ours {
namespace detail {

class ActiveHandleSetND {
 public:
  void Init(usize n) {
    items_.clear();
    items_.reserve(n);
    pos_.assign(n, -1);
  }

  void Clear() {
    items_.clear();
    pos_.clear();
  }

  void Reset() {
    items_.clear();
    std::fill(pos_.begin(), pos_.end(), -1);
  }

  void Insert(u32 handle) {
    SJS_DASSERT(static_cast<usize>(handle) < pos_.size());
    SJS_DASSERT(pos_[handle] < 0);
    pos_[handle] = static_cast<i64>(items_.size());
    items_.push_back(handle);
  }

  void Erase(u32 handle) {
    SJS_DASSERT(static_cast<usize>(handle) < pos_.size());
    const i64 p = pos_[handle];
    SJS_DASSERT(p >= 0);
    const usize up = static_cast<usize>(p);
    const u32 last = items_.back();
    items_[up] = last;
    pos_[last] = p;
    items_.pop_back();
    pos_[handle] = -1;
  }

  bool Contains(u32 handle) const noexcept {
    return static_cast<usize>(handle) < pos_.size() && pos_[handle] >= 0;
  }

  const std::vector<u32>& Items() const noexcept { return items_; }
  usize Size() const noexcept { return items_.size(); }

 private:
  std::vector<u32> items_;
  std::vector<i64> pos_;
};

template <int Dim, class T>
class OursNDContext {
 public:
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;
  using ModeFamilyT = ModeFamilyND<Dim, T>;

  void Reset() {
    ds_ = nullptr;
    built_ = false;
    start_id_of_event_.clear();
    ev_kind_side_.clear();
    ev_handle_.clear();
    num_start_events_ = 0;
    active_r_.Clear();
    active_s_.Clear();
    modes_r_.Clear();
    modes_s_.Clear();
  }

  bool Build(const DatasetT& ds, PhaseRecorder* phases, std::string* err) {
    static_assert(Dim >= 3, "OursNDContext is intended for Dim >= 3");
    if (err) err->clear();
    Reset();
    ds_ = &ds;

    std::vector<join::Event> events;
    {
      auto scoped = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events = join::BuildSweepEvents<Dim, T>(ds.R, ds.S, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    start_id_of_event_.assign(events.size(), -1);
    ev_kind_side_.resize(events.size());
    ev_handle_.resize(events.size());

    for (usize pos = 0; pos < events.size(); ++pos) {
      const auto& ev = events[pos];
      ev_kind_side_[pos] = static_cast<u8>((static_cast<u8>(ev.kind) << 1) | static_cast<u8>(ev.side));
      ev_handle_[pos] = static_cast<u32>(ev.index);
      if (ev.kind == join::EventKind::Start) {
        start_id_of_event_[pos] = static_cast<i32>(num_start_events_++);
      }
    }

    {
      auto scoped = phases ? phases->Scoped("build_mode_families") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!modes_r_.Build(ds.R, err)) return false;
      if (!modes_s_.Build(ds.S, err)) return false;
    }

    active_r_.Init(ds.R.Size());
    active_s_.Init(ds.S.Size());
    built_ = true;
    return true;
  }

  void ResetActive() {
    // The full sweep should normally end with empty active sets. This reset is
    // defensive and also supports abort/restart scenarios.
    const auto r_items = active_r_.Items();
    for (u32 h : r_items) modes_r_.EraseAll(h);
    const auto s_items = active_s_.Items();
    for (u32 h : s_items) modes_s_.EraseAll(h);
    active_r_.Reset();
    active_s_.Reset();
  }

  bool built() const noexcept { return built_; }
  const DatasetT* dataset() const noexcept { return ds_; }

  const std::vector<i32>& start_id_of_event() const noexcept { return start_id_of_event_; }
  const std::vector<u8>& ev_kind_side() const noexcept { return ev_kind_side_; }
  const std::vector<u32>& ev_handle() const noexcept { return ev_handle_; }

  usize num_events() const noexcept { return ev_kind_side_.size(); }
  usize num_start_events() const noexcept { return static_cast<usize>(num_start_events_); }

  ActiveHandleSetND& active_r() noexcept { return active_r_; }
  ActiveHandleSetND& active_s() noexcept { return active_s_; }
  const ActiveHandleSetND& active_r() const noexcept { return active_r_; }
  const ActiveHandleSetND& active_s() const noexcept { return active_s_; }

  const BoxT& BoxR(u32 handle) const noexcept { return ds_->R.boxes[static_cast<usize>(handle)]; }
  const BoxT& BoxS(u32 handle) const noexcept { return ds_->S.boxes[static_cast<usize>(handle)]; }
  Id IdR(u32 handle) const noexcept { return ds_->R.GetId(static_cast<usize>(handle)); }
  Id IdS(u32 handle) const noexcept { return ds_->S.GetId(static_cast<usize>(handle)); }

  ModeFamilyT& modes_r() noexcept { return modes_r_; }
  ModeFamilyT& modes_s() noexcept { return modes_s_; }
  const ModeFamilyT& modes_r() const noexcept { return modes_r_; }
  const ModeFamilyT& modes_s() const noexcept { return modes_s_; }

  void InsertR(u32 handle) {
    modes_r_.InsertAll(handle);
    active_r_.Insert(handle);
  }
  void InsertS(u32 handle) {
    modes_s_.InsertAll(handle);
    active_s_.Insert(handle);
  }
  void EraseR(u32 handle) {
    modes_r_.EraseAll(handle);
    active_r_.Erase(handle);
  }
  void EraseS(u32 handle) {
    modes_s_.EraseAll(handle);
    active_s_.Erase(handle);
  }

 private:
  const DatasetT* ds_{nullptr};
  bool built_{false};

  std::vector<i32> start_id_of_event_;
  std::vector<u8> ev_kind_side_;
  std::vector<u32> ev_handle_;
  u32 num_start_events_{0};

  ActiveHandleSetND active_r_;
  ActiveHandleSetND active_s_;
  ModeFamilyT modes_r_;
  ModeFamilyT modes_s_;
};

}  // namespace detail
}  // namespace ours
}  // namespace baselines
}  // namespace sjs
