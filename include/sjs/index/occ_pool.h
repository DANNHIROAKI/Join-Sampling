#pragma once
// sjs/index/occ_pool.h
//
// A shared "occurrence" pool for dynamic segment-tree style indices.
//
// In SJS-HighDims (§2), each active object may be inserted into O(log N)
// buckets per level and recursively into many sub-structures. To support
// O(1) deletion from each bucket (swap-with-last) *without* storing per-bucket
// backref arrays for every object, we keep a global pool of occurrences:
//
//  - Each time an object handle is inserted into a DenseBucket, we allocate an
//    Occ entry that records (bucket pointer, position in bucket) and link it
//    into a per-handle singly-linked list.
//  - To delete an object from the whole index family, we traverse that list and
//    ask each bucket to erase the occurrence in O(1).
//
// This design matches the "dense array + position handle" idea in the paper,
// but makes the handles allocation-friendly and shareable across all modes.

#include "sjs/core/assert.h"
#include "sjs/core/types.h"

#include <limits>
#include <vector>

namespace sjs {
namespace index {

class OccPool {
 public:
  using OccId = u32;
  using EraseFn = void (*)(void* bucket, OccId occ);

  static constexpr OccId kNull = std::numeric_limits<OccId>::max();

  OccPool() = default;

  // Initialize the pool for handles in [0, num_handles).
  // reserve_occ is optional and helps avoid reallocations for heavy workloads.
  void Init(u32 num_handles, u32 reserve_occ = 0) {
    num_handles_ = num_handles;
    head_.assign(static_cast<usize>(num_handles_), kNull);
    occs_.clear();
    occs_.shrink_to_fit();
    if (reserve_occ > 0) occs_.reserve(static_cast<usize>(reserve_occ));
    free_head_ = kNull;
  }

  void Clear() {
    num_handles_ = 0;
    head_.clear();
    occs_.clear();
    free_head_ = kNull;
  }

  // Reset active state while keeping the handle universe size.
  // IMPORTANT: Caller must ensure all buckets have been cleared separately
  // (otherwise buckets would keep dangling occ ids).
  void ResetActive() {
    std::fill(head_.begin(), head_.end(), kNull);
    occs_.clear();
    free_head_ = kNull;
  }

  u32 num_handles() const noexcept { return num_handles_; }
  usize num_occs() const noexcept { return occs_.size(); }

  // Allocate a new occurrence for `handle` and link it into the handle list.
  // `bucket` must remain valid until this occurrence is erased.
  OccId Alloc(u32 handle, void* bucket, EraseFn erase_fn, u32 pos_in_bucket) {
    SJS_DASSERT(handle < num_handles_);
    SJS_DASSERT(bucket != nullptr);
    SJS_DASSERT(erase_fn != nullptr);

    OccId id = kNull;
    if (free_head_ != kNull) {
      id = free_head_;
      free_head_ = occs_[id].next;
    } else {
      id = static_cast<OccId>(occs_.size());
      occs_.push_back(Occ{});
    }

    Occ& o = occs_[id];
    o.next = head_[handle];
    o.pos = pos_in_bucket;
    o.bucket = bucket;
    o.erase_fn = erase_fn;
    o.handle = handle;

    head_[handle] = id;
    return id;
  }

  // Erase all occurrences currently linked to `handle`.
  // This is the preferred deletion path for dynamic indices.
  void EraseAll(u32 handle) {
    SJS_DASSERT(handle < num_handles_);

    OccId cur = head_[handle];
    head_[handle] = kNull;

    while (cur != kNull) {
      const OccId next = occs_[cur].next;  // save before erasing
      Occ& o = occs_[cur];

      SJS_DASSERT(o.erase_fn != nullptr);
      SJS_DASSERT(o.bucket != nullptr);

      // Ask the bucket to remove the element at this occurrence.
      o.erase_fn(o.bucket, cur);

      Free(cur);
      cur = next;
    }
  }

  // Occurrence -> bucket position.
  u32 Pos(OccId occ) const {
    SJS_DASSERT(occ < occs_.size());
    return occs_[occ].pos;
  }

  void SetPos(OccId occ, u32 new_pos) {
    SJS_DASSERT(occ < occs_.size());
    occs_[occ].pos = new_pos;
  }

  // Occurrence -> owning handle (debug / checks).
  u32 HandleOf(OccId occ) const {
    SJS_DASSERT(occ < occs_.size());
    return occs_[occ].handle;
  }

 private:
  struct Occ {
    OccId next{kNull};
    u32 pos{0};
    void* bucket{nullptr};
    EraseFn erase_fn{nullptr};
    u32 handle{0};
  };

  void Free(OccId occ) {
    SJS_DASSERT(occ < occs_.size());
    Occ& o = occs_[occ];
    o.bucket = nullptr;
    o.erase_fn = nullptr;
    o.handle = 0;
    o.pos = 0;
    o.next = free_head_;
    free_head_ = occ;
  }

  u32 num_handles_{0};
  std::vector<OccId> head_;   // per-handle head of occurrence list

  std::vector<Occ> occs_;     // occurrence storage
  OccId free_head_{kNull};    // free-list head (reuses Occ::next)
};

}  // namespace index
}  // namespace sjs
