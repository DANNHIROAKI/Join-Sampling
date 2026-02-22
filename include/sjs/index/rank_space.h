#pragma once
// sjs/index/rank_space.h
//
// RankSpace: global rank domain for strict open-interval queries.
//
// In §2.3 (range-type primitive), we need to answer queries of the form
//   l < L(u) < r
// where L(u) may have duplicates across different objects. To implement strict
// inequalities robustly, we define a unique key alpha(u) = (L(u), id(u)) and
// sort lexicographically. Then an open interval (l,r) maps to a half-open rank
// interval [L,R) via
//   L = upper_bound((l, +inf))   => first key with coordinate > l
//   R = lower_bound((r, 0))      => first key with coordinate >= r
//
// This header provides that mapping and per-object rank lookup.

#include "sjs/core/assert.h"
#include "sjs/core/types.h"

#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace sjs {
namespace index {

template <class ValueT = Scalar, class HandleT = Id>
class RankSpace {
 public:
  static_assert(std::is_arithmetic_v<ValueT>, "RankSpace expects an arithmetic value type");
  static_assert(std::is_integral_v<HandleT> || std::is_enum_v<HandleT>,
                "RankSpace expects an integral-like handle type");

  using Value = ValueT;
  using Handle = HandleT;

  struct Key {
    Value a{};
    Handle id{};
  };

  RankSpace() = default;

  void Clear() {
    keys_.clear();
    rank_of_.clear();
  }

  // Build from per-handle values, where handle is the index in [0,n).
  // This is the common case in this project: handle ids are dense.
  void BuildFromValues(const std::vector<Value>& a_by_handle) {
    const usize n = a_by_handle.size();
    keys_.resize(n);
    rank_of_.assign(n, 0);

    for (usize i = 0; i < n; ++i) {
      keys_[i] = Key{a_by_handle[i], static_cast<Handle>(i)};
    }

    std::sort(keys_.begin(), keys_.end(), [](const Key& x, const Key& y) {
      if (x.a < y.a) return true;
      if (x.a > y.a) return false;
      return static_cast<u64>(x.id) < static_cast<u64>(y.id);
    });

    for (usize r = 0; r < n; ++r) {
      const Handle id = keys_[r].id;
      const usize idx = static_cast<usize>(static_cast<u64>(id));
      SJS_DASSERT(idx < n);
      rank_of_[idx] = static_cast<u32>(r);
    }
  }

  // Build from an explicit (a,id) list.
  // NOTE: If ids are not dense, rank_of_ will be sized to max_id+1.
  void BuildFromKeys(std::vector<Key> keys) {
    keys_ = std::move(keys);
    std::sort(keys_.begin(), keys_.end(), [](const Key& x, const Key& y) {
      if (x.a < y.a) return true;
      if (x.a > y.a) return false;
      return static_cast<u64>(x.id) < static_cast<u64>(y.id);
    });

    Handle max_id = 0;
    for (const auto& k : keys_) max_id = std::max(max_id, k.id);
    const usize n = static_cast<usize>(static_cast<u64>(max_id)) + 1;

    rank_of_.assign(n, 0);
    for (usize r = 0; r < keys_.size(); ++r) {
      const usize idx = static_cast<usize>(static_cast<u64>(keys_[r].id));
      if (idx >= rank_of_.size()) rank_of_.resize(idx + 1, 0);
      rank_of_[idx] = static_cast<u32>(r);
    }
  }

  u32 Size() const noexcept { return static_cast<u32>(keys_.size()); }
  bool Empty() const noexcept { return keys_.empty(); }

  // 0-based rank in [0, N).
  u32 RankOf(Handle id) const {
    const usize idx = static_cast<usize>(static_cast<u64>(id));
    SJS_DASSERT(idx < rank_of_.size());
    return rank_of_[idx];
  }

  // Access by rank.
  const Key& At(u32 rank) const {
    SJS_DASSERT(rank < keys_.size());
    return keys_[rank];
  }

  // First rank with a > ell. Returns N if all a <= ell.
  u32 UpperBound(Value ell) const {
    // Use a probe that is >= every (ell, id) with a==ell by choosing max id.
    const Key probe{ell, static_cast<Handle>(std::numeric_limits<u32>::max())};
    const auto it = std::upper_bound(keys_.begin(), keys_.end(), probe, [](const Key& x, const Key& y) {
      if (x.a < y.a) return true;
      if (x.a > y.a) return false;
      return static_cast<u64>(x.id) < static_cast<u64>(y.id);
    });
    return static_cast<u32>(std::distance(keys_.begin(), it));
  }

  // First rank with a >= r. Returns N if all a < r.
  u32 LowerBound(Value r) const {
    // Use a probe that is <= every (r, id) with a==r by choosing id=0.
    const Key probe{r, static_cast<Handle>(0)};
    const auto it = std::lower_bound(keys_.begin(), keys_.end(), probe, [](const Key& x, const Key& y) {
      if (x.a < y.a) return true;
      if (x.a > y.a) return false;
      return static_cast<u64>(x.id) < static_cast<u64>(y.id);
    });
    return static_cast<u32>(std::distance(keys_.begin(), it));
  }

  // Map open interval (ell, r) to rank interval [L, R).
  // When L >= R, the interval is empty.
  std::pair<u32, u32> OpenToRankRange(Value ell, Value r) const {
    const u32 L = UpperBound(ell);
    const u32 R = LowerBound(r);
    return {L, R};
  }

 private:
  std::vector<Key> keys_;      // sorted by (a, id)
  std::vector<u32> rank_of_;   // indexed by handle id (dense in this project)
};

}  // namespace index
}  // namespace sjs
