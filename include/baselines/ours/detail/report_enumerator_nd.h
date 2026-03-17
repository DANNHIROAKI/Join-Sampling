#pragma once
// baselines/ours/detail/report_enumerator_nd.h
//
// Deterministic join enumerator for Dim >= 3 based on the fixed-mode local
// recursive structures. The constructor materializes the join output once using
// the same event-block decomposition as the SJS sampling path.

#include "baselines/baseline_api.h"
#include "baselines/ours/detail/context_nd.h"
#include "join/join_types.h"

#include <array>
#include <vector>

namespace sjs {
namespace baselines {
namespace ours {
namespace detail {

template <int Dim, class T>
class OursReportJoinEnumeratorND final : public baselines::IJoinEnumerator {
 public:
  static_assert(Dim >= 3, "OursReportJoinEnumeratorND is intended for Dim >= 3");
  static constexpr u32 kModeCount = kModeCountV<Dim>;
  using Ctx = OursNDContext<Dim, T>;

  explicit OursReportJoinEnumeratorND(Ctx* ctx) : ctx_(ctx) {
    SJS_DASSERT(ctx_ != nullptr);
    Materialize();
    Reset();
  }

  void Reset() override {
    pos_ = 0;
    stats_.Reset();
    stats_.num_events = static_cast<u64>(ctx_ ? ctx_->num_events() : 0U);
  }

  bool Next(PairId* out) override {
    if (pos_ >= pairs_.size()) return false;
    if (out) *out = pairs_[pos_];
    ++pos_;
    stats_.output_pairs = static_cast<u64>(pos_);
    return true;
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  void Materialize() {
    pairs_.clear();
    if (!ctx_ || !ctx_->built() || ctx_->dataset() == nullptr) return;

    ctx_->ResetActive();

    std::array<std::vector<u32>, kModeCount> by_mode;
    const auto& ks = ctx_->ev_kind_side();
    const auto& handles = ctx_->ev_handle();

    for (usize pos = 0; pos < ks.size(); ++pos) {
      const u8 kind_side = ks[pos];
      const bool is_start = ((kind_side >> 1U) != 0U);
      const bool is_r = ((kind_side & 1U) == 0U);
      const u32 handle = handles[pos];

      if (!is_start) {
        if (is_r) ctx_->EraseR(handle);
        else ctx_->EraseS(handle);
        continue;
      }

      const auto& q = is_r ? ctx_->BoxR(handle) : ctx_->BoxS(handle);
      const auto& other = is_r ? ctx_->modes_s() : ctx_->modes_r();
      const Id q_id = is_r ? ctx_->IdR(handle) : ctx_->IdS(handle);

      for (auto& v : by_mode) v.clear();
      for (u32 mode = 0; mode < kModeCount; ++mode) {
        other.ReportMode(mode, q, &by_mode[static_cast<usize>(mode)]);
      }

      for (u32 mode = 0; mode < kModeCount; ++mode) {
        for (u32 oh : by_mode[static_cast<usize>(mode)]) {
          if (is_r) {
            pairs_.push_back(PairId{q_id, ctx_->IdS(oh)});
          } else {
            pairs_.push_back(PairId{ctx_->IdR(oh), q_id});
          }
        }
      }

      if (is_r) ctx_->InsertR(handle);
      else ctx_->InsertS(handle);
    }

    ctx_->ResetActive();
  }

  Ctx* ctx_{nullptr};
  std::vector<PairId> pairs_;
  usize pos_{0};
  join::JoinStats stats_{};
};

}  // namespace detail
}  // namespace ours
}  // namespace baselines
}  // namespace sjs
