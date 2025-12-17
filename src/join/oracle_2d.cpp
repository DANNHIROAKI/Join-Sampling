// src/join/oracle_2d.cpp
//
// Dim=2 correctness oracle wrappers.
//
// The oracle logic itself is defined as header-only templates in
//   include/sjs/join/join_oracle.h
//
// This compilation unit provides thin, non-templated wrappers for common
// tooling code (verification, sampling-quality checks) that only needs Dim=2.

#include "sjs/join/join_oracle.h"

#include "sjs/core/types.h"
#include "sjs/io/dataset.h"
#include "sjs/join/join_types.h"

#include <functional>
#include <vector>

namespace sjs {
namespace join {

u64 CountNaive2D(const Relation<2, Scalar>& R,
                const Relation<2, Scalar>& S,
                JoinStats* stats) {
  return CountNaive<2, Scalar>(R, S, stats);
}

u64 CountNaive2D(const Dataset<2, Scalar>& ds, JoinStats* stats) {
  return CountNaive<2, Scalar>(ds.R, ds.S, stats);
}

bool EnumerateNaive2D(const Relation<2, Scalar>& R,
                      const Relation<2, Scalar>& S,
                      const std::function<bool(PairId)>& emit,
                      JoinStats* stats) {
  return EnumerateNaive<2, Scalar>(R, S, emit, stats);
}

bool EnumerateNaive2D(const Dataset<2, Scalar>& ds,
                      const std::function<bool(PairId)>& emit,
                      JoinStats* stats) {
  return EnumerateNaive<2, Scalar>(ds.R, ds.S, emit, stats);
}

std::vector<PairId> CollectNaivePairs2D(const Relation<2, Scalar>& R,
                                       const Relation<2, Scalar>& S,
                                       u64 cap,
                                       JoinStats* stats) {
  return CollectNaivePairs<2, Scalar>(R, S, cap, stats);
}

std::vector<PairId> CollectNaivePairs2D(const Dataset<2, Scalar>& ds,
                                       u64 cap,
                                       JoinStats* stats) {
  return CollectNaivePairs<2, Scalar>(ds.R, ds.S, cap, stats);
}

}  // namespace join
}  // namespace sjs
