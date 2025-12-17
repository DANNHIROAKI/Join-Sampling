#pragma once
// sjs/sampling/rank_sampling.h
//
// Two-pass rank sampling from an enumeration stream.
//
// Motivation (for Enum+Sampling baselines):
//  - You have a deterministic stream of join pairs (or any items) produced by an enumerator.
//  - You want i.i.d. uniform samples *with replacement* from the stream items,
//    without storing the full stream.
//
// Algorithm (two-pass):
//  1) Pass 1: count N = stream length.
//  2) Draw t i.i.d. integers r_j ~ Unif{0,...,N-1}.
//     Store pairs (r_j, j) and sort by r_j.
//  3) Pass 2: re-scan the stream; when reaching position r_j, output that item
//     to sample slot j. Duplicates are naturally handled.
//
// Complexity:
//  - Time: O(N + t log t)
//  - Memory: O(t)
//
// Requirements:
//  - Stream provides Reset() and Next(T* out) -> bool.
//  - Enumeration order must be deterministic across passes (our join streams are).
//
// This module is generic; it works for join::IJoinStream (PairId) or any similar stream.

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/core/rng.h"

#include <algorithm>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace sampling {

namespace rank_detail {

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

// Count items in stream by full scan.
template <class Stream, class T>
inline bool CountStream(Stream* stream, u64* out_n, std::string* err) {
  if (!stream || !out_n) {
    SetErr(err, "CountStream: null stream/out_n");
    return false;
  }
  stream->Reset();
  u64 n = 0;
  T tmp{};
  while (stream->Next(&tmp)) {
    ++n;
  }
  *out_n = n;
  return true;
}

struct RankRec {
  u64 rank;   // [0,N)
  u64 pos;    // original sample position
};

inline bool RankLess(const RankRec& a, const RankRec& b) noexcept {
  if (a.rank < b.rank) return true;
  if (b.rank < a.rank) return false;
  return a.pos < b.pos;
}

}  // namespace rank_detail

// Rank sampling output diagnostics.
struct RankSamplingInfo {
  u64 universe_size = 0;   // N
  u64 requested = 0;       // t
  u64 produced = 0;        // should equal t on success (unless N=0)
};

// Sample t items uniformly with replacement from a stream.
// Returns false only on errors (e.g., null pointers). If the stream is empty (N=0),
// returns true and produces 0 samples (info->produced = 0).
template <class Stream, class T>
inline bool RankSampleWithReplacement(Stream* stream,
                                      u64 t,
                                      Rng* rng,
                                      std::vector<T>* out_samples,
                                      RankSamplingInfo* info = nullptr,
                                      std::string* err = nullptr) {
  if (!stream || !rng || !out_samples) {
    rank_detail::SetErr(err, "RankSampleWithReplacement: null stream/rng/out_samples");
    return false;
  }

  out_samples->clear();
  out_samples->resize(static_cast<usize>(t));

  RankSamplingInfo local_info;
  local_info.requested = t;

  // Pass 1: count N.
  u64 N = 0;
  if (!rank_detail::CountStream<Stream, T>(stream, &N, err)) return false;
  local_info.universe_size = N;

  if (N == 0 || t == 0) {
    out_samples->clear();
    local_info.produced = 0;
    if (info) *info = local_info;
    return true;
  }

  // Draw t ranks.
  std::vector<rank_detail::RankRec> ranks;
  ranks.resize(static_cast<usize>(t));
  for (u64 i = 0; i < t; ++i) {
    ranks[static_cast<usize>(i)] = rank_detail::RankRec{
        /*rank=*/rng->UniformU64(N),
        /*pos=*/static_cast<u64>(i)};
  }
  std::sort(ranks.begin(), ranks.end(), rank_detail::RankLess);

  // Pass 2: pick items at the ranks.
  stream->Reset();

  u64 cur = 0;         // current stream index
  usize k = 0;         // pointer into ranks[]
  T item{};

  while (stream->Next(&item)) {
    // Assign this item to any samples that ask for rank==cur.
    while (k < ranks.size() && ranks[k].rank == cur) {
      (*out_samples)[static_cast<usize>(ranks[k].pos)] = item;
      ++k;
    }
    ++cur;
    if (k == ranks.size()) break;  // done
  }

  if (k != ranks.size()) {
    // This indicates stream length changed between passes (non-deterministic stream),
    // or a bug in stream counting.
    rank_detail::SetErr(err, "RankSampleWithReplacement: stream ended early in second pass (non-deterministic order?)");
    return false;
  }

  local_info.produced = t;
  if (info) *info = local_info;
  return true;
}

// Convenience: return only N (stream length).
template <class Stream, class T>
inline bool CountStreamItems(Stream* stream, u64* out_n, std::string* err = nullptr) {
  return rank_detail::CountStream<Stream, T>(stream, out_n, err);
}

}  // namespace sampling
}  // namespace sjs
