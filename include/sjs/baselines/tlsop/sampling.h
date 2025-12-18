#pragma once
// sjs/baselines/tlsop/sampling.h
//
// TLSOP baseline (Two-layer SOP) from Tsitsigkos et al. (2023), summarized in
// the uploaded note "Tsitsigkos’23.md".
//
// What this header provides
// -------------------------
//  1) A deterministic join enumerator implementing the *unique* TLSOP join stream:
//       - A regular grid of tiles over the dataset domain.
//       - SOP multi-assignment: each rectangle is assigned to all intersecting tiles.
//       - Two-layer classification per tile into A/B/C/D based on (x_l,y_l) vs tile origin.
//       - Only 9 of the 16 class-to-class mini-joins are evaluated (paper's key result),
//         yielding no duplicates *without* reference-point elimination.
//       - Each mini-join is enumerated by a deterministic forward-scan plane-sweep.
//
//  2) Variant::Sampling baseline for TLSOP:
//       - Count(): enumerates the join stream once to compute |J| exactly.
//       - Sample(): draws t i.i.d. ranks in [0,|J|) and re-enumerates the same
//         deterministic stream to pick the requested pairs (uniform with replacement).
//
// Notes
// -----
//  - Geometry uses half-open boxes [lo,hi) in each dimension.
//  - This baseline is currently implemented for Dim==2 (as in Tsitsigkos'23).
//    The index layout is written to be extendable, but the classification and
//    mini-join set are 2D-specific.

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/core/config.h"
#include "sjs/core/rng.h"
#include "sjs/core/timer.h"
#include "sjs/geometry/box.h"
#include "sjs/geometry/predicates.h"
#include "sjs/io/dataset.h"
#include "sjs/join/join_types.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace tlsop {

namespace detail {

// --------------------------
// Extra-map parsing helpers
// --------------------------

inline std::optional<std::string_view> GetExtra(const std::unordered_map<std::string, std::string>& extra,
                                                std::string_view key) {
  auto it = extra.find(std::string(key));
  if (it == extra.end()) return std::nullopt;
  return std::string_view(it->second);
}

inline bool ParseI32(std::string_view s, i32* out) {
  if (!out || s.empty()) return false;
  try {
    std::size_t idx = 0;
    long v = std::stol(std::string(s), &idx, 10);
    if (idx != s.size()) return false;
    *out = static_cast<i32>(v);
    return true;
  } catch (...) {
    return false;
  }
}

inline bool ParseU64(std::string_view s, u64* out) {
  if (!out || s.empty()) return false;
  try {
    std::size_t idx = 0;
    unsigned long long v = std::stoull(std::string(s), &idx, 10);
    if (idx != s.size()) return false;
    *out = static_cast<u64>(v);
    return true;
  } catch (...) {
    return false;
  }
}

inline bool ParseBool(std::string_view s, bool* out) {
  if (!out || s.empty()) return false;
  auto eq = sjs::detail::EqualsIgnoreCase;
  if (eq(s, "1") || eq(s, "true") || eq(s, "yes") || eq(s, "y") || eq(s, "on")) {
    *out = true;
    return true;
  }
  if (eq(s, "0") || eq(s, "false") || eq(s, "no") || eq(s, "n") || eq(s, "off")) {
    *out = false;
    return true;
  }
  return false;
}

inline i32 GetI32Or(const std::unordered_map<std::string, std::string>& extra,
                    std::string_view key,
                    i32 def) {
  if (auto v = GetExtra(extra, key)) {
    i32 x = 0;
    if (ParseI32(*v, &x)) return x;
  }
  return def;
}

inline u64 GetU64Or(const std::unordered_map<std::string, std::string>& extra,
                    std::string_view key,
                    u64 def) {
  if (auto v = GetExtra(extra, key)) {
    u64 x = 0;
    if (ParseU64(*v, &x)) return x;
  }
  return def;
}

inline bool GetBoolOr(const std::unordered_map<std::string, std::string>& extra,
                      std::string_view key,
                      bool def) {
  if (auto v = GetExtra(extra, key)) {
    bool x = false;
    if (ParseBool(*v, &x)) return x;
  }
  return def;
}

// --------------------------
// TLSOP parameter bundle
// --------------------------

struct TLSOPParams {
  // Grid resolution.
  // If 0, we use a conservative default (64).
  u64 nx = 64;
  u64 ny = 64;

  // Safety cap to prevent accidental huge allocations.
  u64 max_tiles = 1'000'000ULL;

  // Sweep axis inside mini-joins (0 = x).
  int sweep_axis = 0;

  // Sort lists by (x_l, id, idx) once during Build and reuse.
  bool reuse_sorted_lists = true;
};

inline TLSOPParams ParseTLSOPParams(const Config& cfg) {
  TLSOPParams p;
  const auto& extra = cfg.run.extra;

  // Naming convention: tlsop_* keys (avoid collisions with other baselines).
  p.nx = GetU64Or(extra, "tlsop_nx", p.nx);
  p.ny = GetU64Or(extra, "tlsop_ny", p.ny);
  p.max_tiles = GetU64Or(extra, "tlsop_max_tiles", p.max_tiles);
  p.reuse_sorted_lists = GetBoolOr(extra, "tlsop_reuse_sort", p.reuse_sorted_lists);

  // Optional: allow tlsop_sweep_axis (currently only 0 supported in 2D).
  p.sweep_axis = GetI32Or(extra, "tlsop_sweep_axis", p.sweep_axis);
  return p;
}

// --------------------------
// Tile + class lists
// --------------------------

enum class ABClass : u8 { A = 0, B = 1, C = 2, D = 3 };

inline constexpr const char* ToString(ABClass c) noexcept {
  switch (c) {
    case ABClass::A: return "A";
    case ABClass::B: return "B";
    case ABClass::C: return "C";
    case ABClass::D: return "D";
  }
  return "?";
}

// A tile stores 8 lists: R^A..R^D and S^A..S^D.
// Indices are into ds.R.boxes / ds.S.boxes.
// We store u32 indices to save memory.
template <int Dim, class T>
struct Tile {
  static_assert(Dim == 2, "TLSOP Tile is currently specialized for Dim==2");
  using BoxT = Box<Dim, T>;

  BoxT region{};  // tile region [lo,hi)

  std::array<std::vector<u32>, 4> r_cls;  // A,B,C,D
  std::array<std::vector<u32>, 4> s_cls;  // A,B,C,D

  void Clear() {
    for (auto& v : r_cls) v.clear();
    for (auto& v : s_cls) v.clear();
  }
};

// Classify a rectangle into A/B/C/D for a given tile origin (x0,y0).
// Rules (Tsitsigkos'23):
//  A: x_l >= x0 && y_l >= y0
//  B: x_l >= x0 && y_l <  y0
//  C: x_l <  x0 && y_l >= y0
//  D: x_l <  x0 && y_l <  y0

template <class T>
inline ABClass Classify2D(const Box<2, T>& b, T tile_x0, T tile_y0) noexcept {
  const bool x_ge = (b.lo.v[0] >= tile_x0);
  const bool y_ge = (b.lo.v[1] >= tile_y0);
  if (x_ge) {
    return y_ge ? ABClass::A : ABClass::B;
  }
  return y_ge ? ABClass::C : ABClass::D;
}

// --------------------------
// Index: build tiles + SOP multi-assignment + A/B/C/D classification
// --------------------------

template <int Dim, class T>
class TLSOPIndex {
 public:
  static_assert(Dim == 2, "TLSOPIndex is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;

  void Reset() {
    ds_ = nullptr;
    built_ = false;
    params_ = TLSOPParams{};
    domain_ = BoxT::Empty();
    nx_ = ny_ = 0;
    x_edges_.clear();
    y_edges_.clear();
    tiles_.clear();
  }

  bool Built() const noexcept { return built_; }
  const DatasetT* DatasetPtr() const noexcept { return ds_; }
  const TLSOPParams& Params() const noexcept { return params_; }
  u64 Nx() const noexcept { return nx_; }
  u64 Ny() const noexcept { return ny_; }
  const BoxT& Domain() const noexcept { return domain_; }
  const std::vector<Tile<Dim, T>>& Tiles() const noexcept { return tiles_; }

  // Build the TLSOP tile index.
  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) {
    auto scoped = phases ? phases->Scoped("tlsop_build") : PhaseRecorder::ScopedPhase(nullptr, "");

    Reset();
    ds_ = &ds;
    params_ = ParseTLSOPParams(cfg);

    // Validate dims/params.
    if (params_.nx == 0) params_.nx = 64;
    if (params_.ny == 0) params_.ny = 64;

    if (params_.sweep_axis != 0) {
      // For 2D TLSOP, the paper's presentation and our implementation assume x-sweep.
      // Keeping this restriction avoids subtle boundary issues.
      if (err) *err = "TLSOPIndex: currently only tlsop_sweep_axis=0 (x) is supported";
      return false;
    }

    if (params_.nx > std::numeric_limits<u32>::max() || params_.ny > std::numeric_limits<u32>::max()) {
      if (err) *err = "TLSOPIndex: nx/ny too large";
      return false;
    }

    nx_ = params_.nx;
    ny_ = params_.ny;

    const u64 num_tiles = nx_ * ny_;
    if (nx_ == 0 || ny_ == 0 || num_tiles == 0) {
      built_ = true;
      return true;
    }
    if (num_tiles > params_.max_tiles) {
      if (err) {
        *err = "TLSOPIndex: requested nx*ny=" + std::to_string(num_tiles) +
               " exceeds tlsop_max_tiles=" + std::to_string(params_.max_tiles);
      }
      return false;
    }

    // Compute domain.
    domain_ = ds.Domain();
    if (domain_.IsEmpty()) {
      // Degenerate dataset (likely empty); build empty tiles with empty regions.
      tiles_.assign(static_cast<usize>(num_tiles), Tile<Dim, T>{});
      built_ = true;
      return true;
    }

    const T x_min = domain_.lo.v[0];
    const T x_max = domain_.hi.v[0];
    const T y_min = domain_.lo.v[1];
    const T y_max = domain_.hi.v[1];

    const T dx = x_max - x_min;
    const T dy = y_max - y_min;

    if (!(dx > T(0)) || !(dy > T(0))) {
      if (err) *err = "TLSOPIndex: domain has non-positive extent";
      return false;
    }

    // Precompute tile edges.
    {
      auto _ = phases ? phases->Scoped("tlsop_build_edges") : PhaseRecorder::ScopedPhase(nullptr, "");
      x_edges_.resize(static_cast<usize>(nx_ + 1));
      y_edges_.resize(static_cast<usize>(ny_ + 1));
      const T wx = dx / static_cast<T>(nx_);
      const T wy = dy / static_cast<T>(ny_);

      for (u64 i = 0; i <= nx_; ++i) {
        if (i == nx_) {
          x_edges_[static_cast<usize>(i)] = x_max;
        } else {
          x_edges_[static_cast<usize>(i)] = x_min + static_cast<T>(i) * wx;
        }
      }
      for (u64 j = 0; j <= ny_; ++j) {
        if (j == ny_) {
          y_edges_[static_cast<usize>(j)] = y_max;
        } else {
          y_edges_[static_cast<usize>(j)] = y_min + static_cast<T>(j) * wy;
        }
      }
    }

    // Allocate tiles and set their regions.
    {
      auto _ = phases ? phases->Scoped("tlsop_alloc_tiles") : PhaseRecorder::ScopedPhase(nullptr, "");
      tiles_.assign(static_cast<usize>(num_tiles), Tile<Dim, T>{});
      for (u64 j = 0; j < ny_; ++j) {
        for (u64 i = 0; i < nx_; ++i) {
          const usize tid = static_cast<usize>(i + nx_ * j);
          Tile<Dim, T>& Tt = tiles_[tid];
          Point<Dim, T> lo;
          Point<Dim, T> hi;
          lo.v[0] = x_edges_[static_cast<usize>(i)];
          hi.v[0] = x_edges_[static_cast<usize>(i + 1)];
          lo.v[1] = y_edges_[static_cast<usize>(j)];
          hi.v[1] = y_edges_[static_cast<usize>(j + 1)];
          Tt.region = BoxT(lo, hi);
        }
      }
    }

    // Two-pass multi-assignment: pass1 counts capacities; pass2 fills vectors.
    std::vector<std::array<u64, 8>> caps;  // [RA,RB,RC,RD, SA,SB,SC,SD]
    caps.assign(static_cast<usize>(num_tiles), std::array<u64, 8>{});

    // Helper: compute tile index range for one axis under half-open semantics.
    auto tile_range_1d = [&](T lo, T hi, const std::vector<T>& edges, u64 n, u64* out_a, u64* out_b) {
      if (n == 0) {
        *out_a = *out_b = 0;
        return;
      }
      const T minv = edges.front();
      const T maxv = edges.back();
      const T d = maxv - minv;
      const T w = d / static_cast<T>(n);
      const T invw = static_cast<T>(1) / w;

      // Lower index: floor((lo - min)/w)
      long long a = static_cast<long long>(std::floor((lo - minv) * invw));

      // Upper index: floor((prev(hi) - min)/w)
      T hi_prev = hi;
      if constexpr (std::numeric_limits<T>::has_infinity) {
        hi_prev = std::nextafter(hi, -std::numeric_limits<T>::infinity());
      } else {
        hi_prev = std::nextafter(hi, std::numeric_limits<T>::lowest());
      }
      long long b = static_cast<long long>(std::floor((hi_prev - minv) * invw));

      if (a < 0) a = 0;
      if (b < 0) b = 0;
      if (a >= static_cast<long long>(n)) a = static_cast<long long>(n) - 1;
      if (b >= static_cast<long long>(n)) b = static_cast<long long>(n) - 1;

      if (b < a) {
        // No intersection; produce empty range.
        *out_a = 1;
        *out_b = 0;
        return;
      }

      *out_a = static_cast<u64>(a);
      *out_b = static_cast<u64>(b);
    };

    auto count_side = [&](const Relation<Dim, T>& rel, bool is_r_side) {
      const int base = is_r_side ? 0 : 4;
      for (usize idx = 0; idx < rel.boxes.size(); ++idx) {
        if (idx > static_cast<usize>(std::numeric_limits<u32>::max())) {
          if (err) *err = "TLSOPIndex: relation too large for u32 indices";
          throw std::runtime_error("TLSOPIndex: relation too large");
        }
        const auto& b = rel.boxes[idx];
        if (b.IsEmpty()) continue;

        u64 ix0 = 0, ix1 = 0, iy0 = 0, iy1 = 0;
        tile_range_1d(b.lo.v[0], b.hi.v[0], x_edges_, nx_, &ix0, &ix1);
        tile_range_1d(b.lo.v[1], b.hi.v[1], y_edges_, ny_, &iy0, &iy1);
        if (ix0 > ix1 || iy0 > iy1) continue;

        for (u64 j = iy0; j <= iy1; ++j) {
          const T ty0 = y_edges_[static_cast<usize>(j)];
          for (u64 i = ix0; i <= ix1; ++i) {
            const T tx0 = x_edges_[static_cast<usize>(i)];
            const usize tid = static_cast<usize>(i + nx_ * j);
            const ABClass c = Classify2D<T>(reinterpret_cast<const Box<2, T>&>(b), tx0, ty0);
            const int ci = static_cast<int>(c);
            caps[tid][static_cast<usize>(base + ci)] += 1;
          }
        }
      }
    };

    // Pass 1: counts.
    {
      auto _ = phases ? phases->Scoped("tlsop_pass1_count_caps") : PhaseRecorder::ScopedPhase(nullptr, "");
      try {
        count_side(ds.R, /*is_r_side=*/true);
        count_side(ds.S, /*is_r_side=*/false);
      } catch (const std::exception& e) {
        if (err && err->empty()) *err = e.what();
        return false;
      }
    }

    // Reserve capacities.
    {
      auto _ = phases ? phases->Scoped("tlsop_pass1_reserve") : PhaseRecorder::ScopedPhase(nullptr, "");
      for (usize tid = 0; tid < tiles_.size(); ++tid) {
        for (int c = 0; c < 4; ++c) {
          tiles_[tid].r_cls[static_cast<usize>(c)].reserve(static_cast<usize>(caps[tid][static_cast<usize>(c)]));
          tiles_[tid].s_cls[static_cast<usize>(c)].reserve(static_cast<usize>(caps[tid][static_cast<usize>(4 + c)]));
        }
      }
    }

    auto fill_side = [&](const Relation<Dim, T>& rel, bool is_r_side) {
      for (usize idx = 0; idx < rel.boxes.size(); ++idx) {
        const auto& b = rel.boxes[idx];
        if (b.IsEmpty()) continue;

        u64 ix0 = 0, ix1 = 0, iy0 = 0, iy1 = 0;
        tile_range_1d(b.lo.v[0], b.hi.v[0], x_edges_, nx_, &ix0, &ix1);
        tile_range_1d(b.lo.v[1], b.hi.v[1], y_edges_, ny_, &iy0, &iy1);
        if (ix0 > ix1 || iy0 > iy1) continue;

        const u32 uidx = static_cast<u32>(idx);

        for (u64 j = iy0; j <= iy1; ++j) {
          const T ty0 = y_edges_[static_cast<usize>(j)];
          for (u64 i = ix0; i <= ix1; ++i) {
            const T tx0 = x_edges_[static_cast<usize>(i)];
            const usize tid = static_cast<usize>(i + nx_ * j);
            const ABClass c = Classify2D<T>(reinterpret_cast<const Box<2, T>&>(b), tx0, ty0);
            const int ci = static_cast<int>(c);
            if (is_r_side) {
              tiles_[tid].r_cls[static_cast<usize>(ci)].push_back(uidx);
            } else {
              tiles_[tid].s_cls[static_cast<usize>(ci)].push_back(uidx);
            }
          }
        }
      }
    };

    // Pass 2: fill.
    {
      auto _ = phases ? phases->Scoped("tlsop_pass2_fill") : PhaseRecorder::ScopedPhase(nullptr, "");
      fill_side(ds.R, /*is_r_side=*/true);
      fill_side(ds.S, /*is_r_side=*/false);
    }

    // Sort lists (optional but recommended for deterministic plane-sweep mini-join).
    if (params_.reuse_sorted_lists) {
      auto _ = phases ? phases->Scoped("tlsop_sort_lists") : PhaseRecorder::ScopedPhase(nullptr, "");

      auto sort_rel_list = [&](std::vector<u32>& idxs, const Relation<Dim, T>& rel) {
        std::sort(idxs.begin(), idxs.end(), [&](u32 a, u32 b) {
          const auto& ba = rel.boxes[static_cast<usize>(a)];
          const auto& bb = rel.boxes[static_cast<usize>(b)];
          const T ax = ba.lo.v[0];
          const T bx = bb.lo.v[0];
          if (ax < bx) return true;
          if (bx < ax) return false;
          const Id ida = rel.GetId(static_cast<usize>(a));
          const Id idb = rel.GetId(static_cast<usize>(b));
          if (ida < idb) return true;
          if (idb < ida) return false;
          return a < b;
        });
      };

      for (auto& tile : tiles_) {
        for (int c = 0; c < 4; ++c) {
          sort_rel_list(tile.r_cls[static_cast<usize>(c)], ds.R);
          sort_rel_list(tile.s_cls[static_cast<usize>(c)], ds.S);
        }
      }
    }

    built_ = true;
    return true;
  }

 private:
  const DatasetT* ds_ = nullptr;
  bool built_ = false;
  TLSOPParams params_{};

  BoxT domain_ = BoxT::Empty();
  u64 nx_ = 0;
  u64 ny_ = 0;
  std::vector<T> x_edges_;
  std::vector<T> y_edges_;
  std::vector<Tile<Dim, T>> tiles_;
};

// --------------------------
// Mini-join forward-scan (deterministic)
// --------------------------

// Check intersection on all axes except sweep_axis (assumes overlap on sweep axis
// implied by the forward-scan loop).
template <int Dim, class T>
inline bool IntersectsOtherAxes(const Box<Dim, T>& a, const Box<Dim, T>& b, int sweep_axis) {
  if (a.IsEmpty() || b.IsEmpty()) return false;
  for (int d = 0; d < Dim; ++d) {
    if (d == sweep_axis) continue;
    const T alo = a.lo.v[static_cast<usize>(d)];
    const T ahi = a.hi.v[static_cast<usize>(d)];
    const T blo = b.lo.v[static_cast<usize>(d)];
    const T bhi = b.hi.v[static_cast<usize>(d)];
    if (!(alo < bhi && blo < ahi)) return false;
  }
  return true;
}

// A stateful forward-scan mini-join stream over two lists (R_idx, S_idx).
// Lists must be sorted by x_l (lo on sweep axis) with deterministic tie-break.

template <int Dim, class T>
class MiniJoinStream {
 public:
  using DatasetT = Dataset<Dim, T>;

  MiniJoinStream() = default;

  void Init(const DatasetT* ds,
            const std::vector<u32>* r_idx,
            const std::vector<u32>* s_idx,
            int sweep_axis,
            join::JoinStats* stats) {
    ds_ = ds;
    r_idx_ = r_idx;
    s_idx_ = s_idx;
    sweep_axis_ = sweep_axis;
    stats_ = stats;
    Reset();
  }

  void Reset() {
    i_ = 0;
    j_ = 0;
    scanning_ = false;
    fixed_side_r_ = true;
    fixed_idx_ = 0;
    scan_pos_ = 0;
  }

  bool Next(PairId* out) {
    if (!out) return false;
    if (!ds_ || !r_idx_ || !s_idx_) return false;

    const auto& ds = *ds_;
    const auto& R = ds.R;
    const auto& S = ds.S;
    const auto& r_list = *r_idx_;
    const auto& s_list = *s_idx_;

    if (r_list.empty() || s_list.empty()) return false;

    while (true) {
      if (scanning_) {
        if (fixed_side_r_) {
          if (fixed_idx_ >= r_list.size()) {
            scanning_ = false;
            continue;
          }
          const u32 r_i = r_list[fixed_idx_];
          const auto& rb = R.boxes[static_cast<usize>(r_i)];
          const T r_hi = rb.hi.v[static_cast<usize>(sweep_axis_)];

          if (scan_pos_ >= s_list.size()) {
            scanning_ = false;
            ++i_;
            continue;
          }

          const u32 s_j = s_list[scan_pos_];
          const auto& sb = S.boxes[static_cast<usize>(s_j)];
          const T s_lo = sb.lo.v[static_cast<usize>(sweep_axis_)];

          // Stop: s_lo >= r_hi => no further S can overlap this R on sweep axis.
          if (!(s_lo < r_hi)) {
            scanning_ = false;
            ++i_;
            continue;
          }

          ++scan_pos_;
          if (stats_) stats_->candidate_checks++;

          if (!IntersectsOtherAxes<Dim, T>(rb, sb, sweep_axis_)) {
            continue;
          }

          *out = PairId{R.GetId(static_cast<usize>(r_i)), S.GetId(static_cast<usize>(s_j))};
          if (stats_) stats_->output_pairs++;
          return true;
        } else {
          // fixed side is S
          if (fixed_idx_ >= s_list.size()) {
            scanning_ = false;
            continue;
          }
          const u32 s_j = s_list[fixed_idx_];
          const auto& sb = S.boxes[static_cast<usize>(s_j)];
          const T s_hi = sb.hi.v[static_cast<usize>(sweep_axis_)];

          if (scan_pos_ >= r_list.size()) {
            scanning_ = false;
            ++j_;
            continue;
          }

          const u32 r_i = r_list[scan_pos_];
          const auto& rb = R.boxes[static_cast<usize>(r_i)];
          const T r_lo = rb.lo.v[static_cast<usize>(sweep_axis_)];

          if (!(r_lo < s_hi)) {
            scanning_ = false;
            ++j_;
            continue;
          }

          ++scan_pos_;
          if (stats_) stats_->candidate_checks++;

          if (!IntersectsOtherAxes<Dim, T>(rb, sb, sweep_axis_)) {
            continue;
          }

          *out = PairId{R.GetId(static_cast<usize>(r_i)), S.GetId(static_cast<usize>(s_j))};
          if (stats_) stats_->output_pairs++;
          return true;
        }
      }

      // Not scanning: choose next step.
      if (i_ >= r_list.size() || j_ >= s_list.size()) {
        return false;
      }

      const u32 r_i = r_list[i_];
      const u32 s_j = s_list[j_];
      const auto& rb = R.boxes[static_cast<usize>(r_i)];
      const auto& sb = S.boxes[static_cast<usize>(s_j)];
      const T r_lo = rb.lo.v[static_cast<usize>(sweep_axis_)];
      const T s_lo = sb.lo.v[static_cast<usize>(sweep_axis_)];

      if (r_lo < s_lo) {
        scanning_ = true;
        fixed_side_r_ = true;
        fixed_idx_ = i_;
        scan_pos_ = j_;
      } else {
        scanning_ = true;
        fixed_side_r_ = false;
        fixed_idx_ = j_;
        scan_pos_ = i_;
      }
      // loop continues
    }
  }

 private:
  const DatasetT* ds_ = nullptr;
  const std::vector<u32>* r_idx_ = nullptr;
  const std::vector<u32>* s_idx_ = nullptr;
  int sweep_axis_ = 0;
  join::JoinStats* stats_ = nullptr;  // accumulates candidate/output stats

  usize i_ = 0;
  usize j_ = 0;
  bool scanning_ = false;
  bool fixed_side_r_ = true;
  usize fixed_idx_ = 0;
  usize scan_pos_ = 0;
};

// --------------------------
// TLSOP join enumerator (streaming, deterministic)
// --------------------------

struct MiniJoinSpec {
  ABClass r;
  ABClass s;
};

inline constexpr MiniJoinSpec kMiniJoins9[9] = {
    // 1) RA ⋈ SA
    {ABClass::A, ABClass::A},
    // 2) RA ⋈ SB
    {ABClass::A, ABClass::B},
    // 3) RA ⋈ SC
    {ABClass::A, ABClass::C},
    // 4) RA ⋈ SD
    {ABClass::A, ABClass::D},
    // 5) RB ⋈ SA
    {ABClass::B, ABClass::A},
    // 6) RB ⋈ SC
    {ABClass::B, ABClass::C},
    // 7) RC ⋈ SA
    {ABClass::C, ABClass::A},
    // 8) RC ⋈ SB
    {ABClass::C, ABClass::B},
    // 9) RD ⋈ SA
    {ABClass::D, ABClass::A},
};

// A deterministic stream of unique join pairs according to TLSOP.

template <int Dim, class T>
class TLSOPJoinEnumerator final : public baselines::IJoinEnumerator {
 public:
  static_assert(Dim == 2, "TLSOPJoinEnumerator is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;

  explicit TLSOPJoinEnumerator(const TLSOPIndex<Dim, T>* index)
      : index_(index) {
    SJS_ASSERT(index_ != nullptr);
    Reset();
  }

  void Reset() override {
    stats_.Reset();
    tile_pos_ = 0;
    mj_pos_ = 0;
    mj_active_ = false;
    mini_.Reset();

    // A deterministic proxy for "events": total assigned objects across tiles.
    u64 n_events = 0;
    if (index_ && index_->Built()) {
      for (const auto& tile : index_->Tiles()) {
        for (int c = 0; c < 4; ++c) {
          n_events += static_cast<u64>(tile.r_cls[static_cast<usize>(c)].size());
          n_events += static_cast<u64>(tile.s_cls[static_cast<usize>(c)].size());
        }
      }
    }
    stats_.num_events = n_events;
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

  bool Next(PairId* out) override {
    if (!out) return false;
    if (!index_ || !index_->Built() || !index_->DatasetPtr()) return false;

    const auto* ds = index_->DatasetPtr();
    const auto& tiles = index_->Tiles();
    const int sweep_axis = index_->Params().sweep_axis;

    while (tile_pos_ < tiles.size()) {
      const auto& tile = tiles[tile_pos_];

      // If no active mini-join, initialize the current mj_pos_.
      if (!mj_active_) {
        if (mj_pos_ >= 9) {
          // move to next tile
          ++tile_pos_;
          mj_pos_ = 0;
          continue;
        }

        const MiniJoinSpec spec = kMiniJoins9[mj_pos_];
        const auto& r_list = tile.r_cls[static_cast<usize>(static_cast<int>(spec.r))];
        const auto& s_list = tile.s_cls[static_cast<usize>(static_cast<int>(spec.s))];

        // Skip empty mini-joins.
        if (r_list.empty() || s_list.empty()) {
          ++mj_pos_;
          continue;
        }

        mini_.Init(ds, &r_list, &s_list, sweep_axis, &stats_);
        mj_active_ = true;
      }

      // Try to emit next pair from current mini-join.
      if (mini_.Next(out)) {
        return true;
      }

      // mini-join exhausted
      mj_active_ = false;
      ++mj_pos_;
    }

    return false;
  }

 private:
  const TLSOPIndex<Dim, T>* index_ = nullptr;

  // Iteration state
  usize tile_pos_ = 0;
  usize mj_pos_ = 0;  // 0..8
  bool mj_active_ = false;
  MiniJoinStream<Dim, T> mini_;

  join::JoinStats stats_{};
};

}  // namespace detail

// --------------------------
// TLSOP Sampling baseline
// --------------------------

template <int Dim, class T = Scalar>
class TLSOPSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim == 2, "TLSOPSamplingBaseline is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::TLSOP; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "tlsop_sampling"; }

  void Reset() override {
    index_.Reset();
    built_ = false;
    have_count_ = false;
    W_ = 0;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    built_ = false;
    have_count_ = false;
    W_ = 0;
    if (!index_.Build(ds, cfg, phases, err)) return false;
    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;  // deterministic
    if (!built_ || !index_.Built()) {
      if (err) *err = "TLSOPSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "TLSOPSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("tlsop_count_pass1") : PhaseRecorder::ScopedPhase(nullptr, "");

    detail::TLSOPJoinEnumerator<Dim, T> stream(&index_);
    stream.Reset();

    u64 W = 0;
    PairId tmp;
    while (stream.Next(&tmp)) {
      ++W;
    }

    W_ = W;
    have_count_ = true;
    *out = MakeExactCount(W);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !index_.Built()) {
      if (err) *err = "TLSOPSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "TLSOPSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "TLSOPSamplingBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->weighted = false;
    out->with_replacement = true;

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    // Ensure we know W.
    if (!have_count_) {
      CountResult cr;
      if (!Count(cfg, /*rng=*/nullptr, &cr, phases, err)) return false;
    }
    const u64 W = W_;
    if (W == 0) return true;

    // Draw ranks in [0,W) with replacement and sort them.
    struct RankReq {
      u64 rank;
      u64 slot;
    };
    std::vector<RankReq> req;
    req.resize(static_cast<usize>(t));
    {
      auto scoped2 = phases ? phases->Scoped("tlsop_sample_draw_ranks") : PhaseRecorder::ScopedPhase(nullptr, "");
      for (u64 i = 0; i < t; ++i) {
        req[static_cast<usize>(i)] = RankReq{rng->UniformU64(W), i};
      }
      std::sort(req.begin(), req.end(), [](const RankReq& a, const RankReq& b) {
        if (a.rank < b.rank) return true;
        if (b.rank < a.rank) return false;
        return a.slot < b.slot;
      });
    }

    // Pass 2: scan join stream and pick pairs at requested ranks.
    out->pairs.assign(static_cast<usize>(t), PairId{});
    {
      auto scoped2 = phases ? phases->Scoped("tlsop_sample_pass2_select") : PhaseRecorder::ScopedPhase(nullptr, "");
      detail::TLSOPJoinEnumerator<Dim, T> stream(&index_);
      stream.Reset();
      PairId p;
      u64 idx = 0;
      usize j = 0;
      while (stream.Next(&p)) {
        while (j < req.size() && req[j].rank == idx) {
          out->pairs[static_cast<usize>(req[j].slot)] = p;
          ++j;
        }
        if (j >= req.size()) break;
        ++idx;
      }
      if (j != req.size()) {
        if (err) *err = "TLSOPSamplingBaseline::Sample: stream ended early in pass2 (non-deterministic enumerate?)";
        return false;
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("tlsop_enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !index_.Built()) {
      if (err) *err = "TLSOPSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::TLSOPJoinEnumerator<Dim, T>>(&index_);
  }

 private:
  detail::TLSOPIndex<Dim, T> index_;
  bool built_ = false;

  bool have_count_ = false;
  u64 W_ = 0;
};

}  // namespace tlsop
}  // namespace baselines
}  // namespace sjs
