#pragma once
// sjs/index/mode_family.h
//
// §2.5 Mode family: across-mode aggregation for event-block primitives.
//
// For a fixed projected dimension count m=Dim-1, modes are g in {A,B}^m.
// For each mode g, we maintain a ModeIndexND<Dim> instance D_{m,g}.
// Modes partition the intersection answer set (Theorem 2.2), thus:
//
// COUNT(q)   = sum_g COUNT_g(q)
// REPORT(q)  = concat_g REPORT_g(q)   (no duplicates)
// SAMPLE(q,k): alias over modes by w^{(g)}(q)=COUNT_g(q), then grouped mode sampling,
//             yielding i.i.d. uniform samples over K(q) (Proposition 2.11).

#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/core/types.h"

#include "sjs/geometry/box.h"

#include "sjs/index/mode_index.h"
#include "sjs/index/rank_space.h"
#include "sjs/index/segtree_common.h"

#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace index {

template <int Dim, class ValueT = Scalar, class HandleT = Id>
class ModeFamily {
 public:
  static_assert(Dim >= 2, "ModeFamily requires Dim >= 2.");

  using Value = ValueT;
  using Handle = HandleT;
  using BoxT = Box<Dim, Value>;
  using AxisData = ModeAxisData<Value, Handle>;
  using RankSpaceT = RankSpace<Value, Handle>;
  using IndexT = ModeIndexND<Dim, Value, Handle>;

  ModeFamily() = default;

  ModeFamily(const ModeFamily&) = delete;
  ModeFamily& operator=(const ModeFamily&) = delete;
  ModeFamily(ModeFamily&&) noexcept = default;
  ModeFamily& operator=(ModeFamily&&) noexcept = default;

  // Initialize with the static universe of boxes that will ever be inserted.
  // Handles are assumed to be dense indices [0, boxes.size()).
  void Init(const std::vector<BoxT>* boxes) {
    boxes_ = boxes;
    SJS_ASSERT(boxes_ != nullptr);

    const u32 N = static_cast<u32>(boxes_->size());

    // Shared occurrence pool: deletion across all modes and all recursion nodes.
    pool_.Init(N);

    // Build per-axis universes for projected axes 1..Dim-1 (m = Dim-1).
    const int m = Dim - 1;
    axes_.clear();
    axes_.resize(static_cast<usize>(m));

    for (int axis_offset = 0; axis_offset < m; ++axis_offset) {
      const int axis = axis_offset + 1;

      AxisData ax;

      // ---- coords for stabbing (A) ----
      {
        std::vector<Value> coords;
        coords.reserve(static_cast<usize>(2) * boxes_->size());
        for (usize i = 0; i < boxes_->size(); ++i) {
          const auto& b = (*boxes_)[i];
          coords.push_back(b.lo[axis]);
          coords.push_back(b.hi[axis]);
        }
        std::sort(coords.begin(), coords.end());
        coords.erase(std::unique(coords.begin(), coords.end()), coords.end());

        ax.coords = std::make_shared<const std::vector<Value>>(std::move(coords));
        ax.stab_leaves = (ax.coords->size() >= 2) ? static_cast<u32>(ax.coords->size() - 1) : 0;
        ax.stab_base = NextPow2(ax.stab_leaves);
      }

      // ---- rank space for range (B) ----
      {
        std::vector<Value> lefts;
        lefts.resize(boxes_->size());
        for (usize i = 0; i < boxes_->size(); ++i) {
          lefts[i] = (*boxes_)[i].lo[axis];
        }
        auto rank = std::make_shared<RankSpaceT>();
        // NOTE: RankSpace exposes BuildFromValues (dense handles) in this repo.
        rank->BuildFromValues(lefts);
        ax.rank = std::move(rank);

        ax.rank_leaves = N;
        ax.rank_base = NextPow2(ax.rank_leaves);
      }

      axes_[static_cast<usize>(axis_offset)] = std::move(ax);
    }

    // Build mode indices: 2^m modes.
    // NOTE: if m is too large, this is infeasible; keep Dim small in practice.
    const u64 num_modes = (m >= 63) ? 0ULL : (1ULL << static_cast<u64>(m));
    SJS_ASSERT_MSG(num_modes > 0, "ModeFamily: Dim too large for mode enumeration");
    SJS_ASSERT_MSG(num_modes <= 1'048'576ULL, "ModeFamily: too many modes (Dim too large in practice)");

    modes_.clear();
    modes_.resize(static_cast<usize>(num_modes));

    for (u64 mask = 0; mask < num_modes; ++mask) {
      modes_[static_cast<usize>(mask)].Init(
          boxes_, &axes_, &pool_, /*axis_offset=*/0, static_cast<u32>(mask));
    }
  }

  bool IsInitialized() const noexcept { return boxes_ != nullptr; }

  // Reset the active set while keeping the static universes (boxes_/axes_) intact.
  //
  // This is used by sweep-style algorithms that repeatedly scan events and want a clean
  // empty active set between passes.
  //
  // Implementation note:
  //  - We must clear *all* buckets/children that currently store Occ ids before calling
  //    OccPool::ResetActive(), otherwise those buckets would keep dangling Occ ids.
  //  - ModeIndexND does not expose a deep clear API here, so we re-init each mode index
  //    in-place (destroying all buckets/children), then reset the pool.
  void ResetActive() {
    if (!IsInitialized()) return;

    const u64 num_modes = static_cast<u64>(modes_.size());
    for (u64 mask = 0; mask < num_modes; ++mask) {
      modes_[static_cast<usize>(mask)] = IndexT{};
      modes_[static_cast<usize>(mask)].Init(
          boxes_, &axes_, &pool_, /*axis_offset=*/0, static_cast<u32>(mask));
    }
    pool_.ResetActive();
  }

  // Insert handle into all modes (active set insertion).
  void Insert(Handle h) {
    if (!IsInitialized()) return;
    for (auto& idx : modes_) idx.Insert(h);
  }

  // Compatibility overloads (some test harnesses expect a box argument even when the
  // index has access to a static box universe by handle).
  void Insert(Handle h, const BoxT&) { Insert(h); }
  void Insert(const BoxT&, Handle h) { Insert(h); }

  // Delete handle from all buckets in all modes (active set deletion).
  // This is O(#occurrences(handle)) and does not require traversing indices.
  void Delete(Handle h) { pool_.EraseAll(static_cast<u32>(h)); }

  // Compatibility alias (baselines often use Erase()).
  void Erase(Handle h) { Delete(h); }

  // COUNT(q): total partners = sum_g COUNT_g(q)
  u64 Count(const BoxT& q) const {
    if (!IsInitialized()) return 0;
    u64 sum = 0;
    for (const auto& idx : modes_) sum += idx.Count(q);
    return sum;
  }

  // REPORT(q): enumerate all partners without duplicates (mode partition).
  void Report(const BoxT& q, std::vector<Handle>* out) const {
    SJS_ASSERT(out != nullptr);
    out->clear();
    if (!IsInitialized()) return;
    for (const auto& idx : modes_) idx.Report(q, out);
  }

  // Convenience return-by-value report (useful for tests / prototyping).
  std::vector<Handle> Report(const BoxT& q) const {
    std::vector<Handle> out;
    Report(q, &out);
    return out;
  }

  // SAMPLE(q,k): i.i.d uniform samples from K(q) (partners), with replacement.
  //
  // Returns true iff it produced exactly k samples (out->size()==k).
  // Returns false if sampling is impossible (e.g., |K(q)|==0 while k>0) or uninitialized.
  //
  // Implements Proposition 2.11: alias over modes + grouped mode sampling.
  bool Sample(const BoxT& q, u64 k, Rng* rng, std::vector<Handle>* out) const {
    SJS_ASSERT(rng != nullptr);
    SJS_ASSERT(out != nullptr);
    out->clear();

    if (k == 0) return true;
    if (!IsInitialized()) return false;

    const usize M = modes_.size();
    if (M == 0) return false;

    std::vector<u64> w(M, 0);
    u64 W = 0;
    for (usize i = 0; i < M; ++i) {
      const u64 c = modes_[i].Count(q);
      w[i] = c;
      W += c;
    }
    if (W == 0) return false;

    sampling::AliasTable alias;
    (void)alias.BuildFromU64(Span<const u64>(w.data(), w.size()));

    out->assign(static_cast<usize>(k), Handle{});

    // Assign each sample position to a mode.
    std::vector<std::vector<usize>> pos(M);
    for (u64 t = 0; t < k; ++t) {
      const usize mi = alias.Sample(rng);
      pos[mi].push_back(static_cast<usize>(t));
    }

    // Grouped mode sampling.
    for (usize mi = 0; mi < M; ++mi) {
      const usize ki = pos[mi].size();
      if (ki == 0) continue;

      std::vector<Handle> tmp;
      modes_[mi].Sample(q, static_cast<u64>(ki), rng, &tmp);
      if (tmp.size() != ki) {
        // This indicates an inconsistency between Count() and Sample() or a logic error
        // upstream (e.g., requesting samples from an empty partner set).
        out->clear();
        return false;
      }

      for (usize j = 0; j < ki; ++j) {
        (*out)[pos[mi][j]] = tmp[j];
      }
    }

    return out->size() == static_cast<usize>(k);
  }

  // Convenience single-sample API (optional).
  bool Sample(const BoxT& q, Rng* rng, Handle* out) const {
    SJS_ASSERT(rng != nullptr);
    SJS_ASSERT(out != nullptr);

    std::vector<Handle> tmp;
    if (!Sample(q, /*k=*/1, rng, &tmp) || tmp.size() != 1) return false;
    *out = tmp[0];
    return true;
  }

  // Expose pool pointer if needed by sweep (optional).
  OccPool* OccurrencePool() noexcept { return &pool_; }
  const OccPool* OccurrencePool() const noexcept { return &pool_; }

  const std::vector<AxisData>& Axes() const noexcept { return axes_; }

 private:
  const std::vector<BoxT>* boxes_{nullptr};

  std::vector<AxisData> axes_;   // size m=Dim-1
  OccPool pool_;                 // shared deletion pool

  std::vector<IndexT> modes_;    // size 2^m
};

}  // namespace index
}  // namespace sjs
