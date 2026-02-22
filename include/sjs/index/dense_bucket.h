#pragma once
// sjs/index/dense_bucket.h
//
// DenseBucket: a "dense array + position handle" bucket with O(1) erase.
//
// This is the concrete bucket used by dynamic segment tree indices in
// SJS-HighDims (§2). Elements are stored in a dense vector; deletion is
// implemented via swap-with-last and uses OccPool to update the moved
// element's position.

#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/core/types.h"
#include "sjs/index/occ_pool.h"

#include <type_traits>
#include <vector>

namespace sjs {
namespace index {

template <class HandleT = Id>
class DenseBucket {
 public:
  static_assert(std::is_integral_v<HandleT> || std::is_enum_v<HandleT>,
                "DenseBucket expects an integral-like handle type");

  using Handle = HandleT;
  using OccId = OccPool::OccId;

  struct Item {
    Handle handle{};
    OccId occ{OccPool::kNull};
  };

  DenseBucket() = default;
  explicit DenseBucket(OccPool* pool) : pool_(pool) {}

  void SetPool(OccPool* pool) { pool_ = pool; }
  OccPool* pool() const noexcept { return pool_; }

  void Clear() { items_.clear(); }
  void Reserve(usize n) { items_.reserve(n); }

  usize Size() const noexcept { return items_.size(); }
  bool Empty() const noexcept { return items_.empty(); }

  const std::vector<Item>& Items() const noexcept { return items_; }

  // Insert a handle and return its occurrence id (mainly for debugging).
  OccId Insert(Handle h) {
    SJS_DASSERT(pool_ != nullptr);
    const u32 pos = static_cast<u32>(items_.size());
    const u32 hid = static_cast<u32>(h);
    const OccId occ = pool_->Alloc(hid, static_cast<void*>(this), &DenseBucket::EraseThunk, pos);
    items_.push_back(Item{h, occ});
    return occ;
  }

  // Uniformly sample one element from this bucket.
  // Returns false if the bucket is empty.
  bool SampleOne(Rng* rng, Handle* out) const {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(out != nullptr);
    if (items_.empty()) return false;
    const u32 pos = rng->UniformU32(static_cast<u32>(items_.size()));
    *out = items_[pos].handle;
    return true;
  }

  // Append all handles to `out`.
  void Report(std::vector<Handle>* out) const {
    SJS_DASSERT(out != nullptr);
    out->reserve(out->size() + items_.size());
    for (const auto& it : items_) out->push_back(it.handle);
  }

  // --- Called by OccPool via type-erased EraseFn ---

  void EraseOcc(OccId occ) {
    SJS_DASSERT(pool_ != nullptr);
    const u32 pos = pool_->Pos(occ);
    SJS_DASSERT(pos < items_.size());
    SJS_DASSERT(items_[pos].occ == occ);

    const usize last_pos = items_.size() - 1;
    if (static_cast<usize>(pos) != last_pos) {
      const Item moved = items_[last_pos];
      items_[pos] = moved;
      pool_->SetPos(moved.occ, pos);
    }
    items_.pop_back();
  }

 private:
  static void EraseThunk(void* bucket, OccId occ) {
    auto* self = static_cast<DenseBucket*>(bucket);
    self->EraseOcc(occ);
  }

  OccPool* pool_{nullptr};
  std::vector<Item> items_;
};

}  // namespace index
}  // namespace sjs
