#pragma once
// sjs/join/join_oracle.h
//
// Small-scale correctness oracle: brute-force join enumeration.
//
// Intended uses:
//  - Verify baseline outputs on small datasets (exact join size/pairs).
//  - Statistical tests for sampling (chi-square, KS) by providing ground truth.
//
// The oracle uses a simple O(|R|*|S|) double loop and the project-standard
// half-open intersection predicate.
//
// To avoid high memory use, you can either:
//  - Count only
//  - Enumerate via callback
//  - Use NaiveJoinStream (streaming, no storage)

#include "sjs/core/assert.h"
#include "sjs/core/types.h"
#include "sjs/geometry/predicates.h"
#include "sjs/io/dataset.h"
#include "sjs/join/join_types.h"

#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace join {

// Count join pairs with a brute-force oracle. Returns |J|.
//
// If stats != nullptr, fills candidate_checks and output_pairs.
template <int Dim, class T>
inline u64 CountNaive(const Relation<Dim, T>& R,
                      const Relation<Dim, T>& S,
                      JoinStats* stats = nullptr) {
  if (stats) stats->Reset();

  u64 count = 0;
  const usize nr = R.boxes.size();
  const usize ns = S.boxes.size();

  for (usize i = 0; i < nr; ++i) {
    const auto& rb = R.boxes[i];
    for (usize j = 0; j < ns; ++j) {
      if (stats) stats->candidate_checks++;
      const auto& sb = S.boxes[j];
      if (IntersectsHalfOpen<Dim, T>(rb, sb)) {
        ++count;
      }
    }
  }

  if (stats) {
    stats->output_pairs = count;
    stats->num_events = 0;
    stats->active_max_r = 0;
    stats->active_max_s = 0;
  }
  return count;
}

// Enumerate join pairs with a brute-force oracle and send them to `emit`.
//
// Callback signature:
//   bool emit(PairId p)
// Return false to stop early.
//
// Returns true if finished full enumeration, false if stopped early.
template <int Dim, class T, class EmitFn>
inline bool EnumerateNaive(const Relation<Dim, T>& R,
                           const Relation<Dim, T>& S,
                           EmitFn&& emit,
                           JoinStats* stats = nullptr) {
  if (stats) stats->Reset();

  const usize nr = R.boxes.size();
  const usize ns = S.boxes.size();

  for (usize i = 0; i < nr; ++i) {
    const auto& rb = R.boxes[i];
    const Id rid = R.GetId(i);
    for (usize j = 0; j < ns; ++j) {
      if (stats) stats->candidate_checks++;
      const auto& sb = S.boxes[j];
      if (IntersectsHalfOpen<Dim, T>(rb, sb)) {
        const PairId p{rid, S.GetId(j)};
        if (stats) stats->output_pairs++;
        if (!emit(p)) {
          // stopped early
          return false;
        }
      }
    }
  }

  return true;
}

// Collect all join pairs (oracle) into a vector.
// Optional cap: if cap>0, stop after collecting `cap` pairs.
//
// WARNING: For large datasets or large joins this can be huge; prefer callback/stream.
template <int Dim, class T>
inline std::vector<PairId> CollectNaivePairs(const Relation<Dim, T>& R,
                                            const Relation<Dim, T>& S,
                                            u64 cap = 0,
                                            JoinStats* stats = nullptr) {
  std::vector<PairId> out;
  if (stats) stats->Reset();

  const usize nr = R.boxes.size();
  const usize ns = S.boxes.size();

  for (usize i = 0; i < nr; ++i) {
    const auto& rb = R.boxes[i];
    const Id rid = R.GetId(i);
    for (usize j = 0; j < ns; ++j) {
      if (stats) stats->candidate_checks++;
      const auto& sb = S.boxes[j];
      if (IntersectsHalfOpen<Dim, T>(rb, sb)) {
        out.push_back(PairId{rid, S.GetId(j)});
        if (stats) stats->output_pairs++;
        if (cap > 0 && out.size() >= static_cast<usize>(cap)) {
          return out;
        }
      }
    }
  }
  return out;
}

// A streaming brute-force join enumerator (lexicographic order in (i,j)).
//
// This is useful when you want a deterministic join stream on small datasets
// without storing all pairs.
template <int Dim, class T = Scalar>
class NaiveJoinStream final : public IJoinStream {
 public:
  NaiveJoinStream(const Relation<Dim, T>& R, const Relation<Dim, T>& S)
      : R_(&R), S_(&S) {
    Reset();
  }

  void Reset() override {
    i_ = 0;
    j_ = 0;
    done_ = false;
  }

  bool Next(PairId* out) override {
    if (!out || done_) return false;

    const usize nr = R_->boxes.size();
    const usize ns = S_->boxes.size();

    while (i_ < nr) {
      const auto& rb = R_->boxes[i_];
      const Id rid = R_->GetId(i_);
      while (j_ < ns) {
        const auto& sb = S_->boxes[j_];
        const Id sid = S_->GetId(j_);
        ++j_;
        if (IntersectsHalfOpen<Dim, T>(rb, sb)) {
          *out = PairId{rid, sid};
          return true;
        }
      }
      ++i_;
      j_ = 0;
    }

    done_ = true;
    return false;
  }

 private:
  const Relation<Dim, T>* R_;
  const Relation<Dim, T>* S_;
  usize i_{0};
  usize j_{0};
  bool done_{false};
};

}  // namespace join
}  // namespace sjs
