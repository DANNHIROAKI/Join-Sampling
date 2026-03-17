#pragma once
// baselines/detail/vector_join_enumerator.h
//
// Helper enumerator over a pre-materialized vector of join pairs.
//
// Motivation
// ----------
// In SJS v3 Framework I (Enumerate+Sampling), the algorithm materializes all
// join pairs into memory once, and then draws i.i.d. uniform samples with
// replacement by indexing that materialized array.  When a baseline keeps such
// a cached vector internally, providing an enumerator over the cached vector
// avoids re-running another sweep just to enumerate again.

#include "baselines/baseline_api.h"

namespace sjs {
namespace baselines {
namespace detail {

class VectorJoinEnumerator final : public IJoinEnumerator {
 public:
  explicit VectorJoinEnumerator(const std::vector<PairId>* pairs) : pairs_(pairs) {
    stats_.Reset();
  }

  void Reset() override {
    pos_ = 0;
    stats_.Reset();
  }

  bool Next(PairId* out) override {
    if (!pairs_) return false;
    if (pos_ >= pairs_->size()) return false;
    if (out) *out = (*pairs_)[pos_];
    ++pos_;
    // Keep a lightweight "so far" counter.
    stats_.output_pairs = static_cast<u64>(pos_);
    return true;
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  const std::vector<PairId>* pairs_ = nullptr;
  usize pos_ = 0;
  join::JoinStats stats_{};
};

}  // namespace detail
}  // namespace baselines
}  // namespace sjs
