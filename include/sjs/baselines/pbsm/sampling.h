#pragma once
// sjs/baselines/pbsm/sampling.h
//
// PBSM (Partition-Based Spatial-Merge) baseline (Variant::Sampling).
//
// Reference
// ---------
// This implementation is an engineering-oriented version of the PBSM family of
// algorithms described in Tsitsigkos et al. (2019) and summarized in
// "Tsitsigkos’19.md" (uploaded by the user).
//
// We implement the *unique* join enumerator:
//   - Multi-assignment of rectangles to partitions (1D stripes by default,
//     optional 2D grid for Dim==2).
//   - Per-partition forward-scan plane sweep (Algorithm-1 style).
//   - Duplicate elimination using the reference-point test (Eq.(1)).
//
// Then the Sampling variant implements *exact* i.i.d. uniform sampling WITH
// replacement from J = R ⋈ S using a two-pass rank-selection scheme:
//   - Count(): enumerates the unique join stream once and returns |J| exactly.
//   - Sample(): draws t i.i.d. ranks in [0,|J|) and re-enumerates the same
//     deterministic stream, selecting the pairs at those ranks.
//
// Notes
// -----
// - Geometry uses half-open boxes [lo, hi) in each dimension.
// - Correctness of the two-pass *rank* sampling (Baseline §4.4 / Theorem 2) only
//   requires Pass2 to be a one-to-one enumeration of J; the enumeration order
//   does NOT need to match Pass1 (Pass1 is only used to obtain W = |J|).
//   For reproducibility and debugging, we still enforce determinism by:
//     * fixed partition traversal order,
//     * sorting per-partition lists by (lower endpoint, id) as a total order,
//     * a deterministic forward-scan state machine.
// - The current project focus is 2D, but the code is written to be extensible
//   to Dim>2 (stripes partitioning is naturally extensible).

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/core/config.h"
#include "sjs/core/rng.h"
#include "sjs/core/timer.h"
#include "sjs/geometry/box.h"
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
namespace pbsm {

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

inline std::string GetStringOr(const std::unordered_map<std::string, std::string>& extra,
                               std::string_view key,
                               std::string def) {
  if (auto v = GetExtra(extra, key)) return std::string(*v);
  return def;
}

// --------------------------
// PBSM parameter bundle
// --------------------------

enum class PartitionScheme : u8 {
  Stripes1D = 0,  // partition along one axis into K stripes
  Grid2D = 1,     // partition 2D domain into KxK tiles (Dim must be 2)
};

enum class SweepAxisPolicy : u8 {
  Orthogonal = 0,  // sweep axis = other axis (for stripes)
  Fixed = 1,       // sweep axis fixed to a given axis
  AutoHist = 2,    // choose sweep axis per partition using histogram cost model
};

struct PBSMParams {
  PartitionScheme scheme = PartitionScheme::Stripes1D;

  // K = number of stripes for Stripes1D; for Grid2D, K is tiles per dimension.
  // If 0, we auto-pick a reasonable K using the "tile ~ 10x object" heuristic.
  u64 K = 0;

  // For Stripes1D: which axis to partition.
  int part_axis = 0;

  // Sweep axis selection.
  SweepAxisPolicy sweep_policy = SweepAxisPolicy::Orthogonal;
  int fixed_sweep_axis = 1;

  // Histogram model knobs for AutoHist.
  int hist_bins = 32;       // number of bins
  int hist_phi = 100;       // subsampling rate: process 1 out of phi items (>=1)
  bool hist_use_domain = true;  // use global domain range for hist instead of partition range

  // Internal safety/perf caps.
  u64 max_partitions = 1'000'000ULL;  // avoid insane memory usage by default
};

inline PartitionScheme ParseScheme(std::string_view s) {
  auto eq = sjs::detail::EqualsIgnoreCase;
  if (eq(s, "grid") || eq(s, "2d") || eq(s, "grid2d") || eq(s, "tile")) {
    return PartitionScheme::Grid2D;
  }
  // default
  return PartitionScheme::Stripes1D;
}

inline SweepAxisPolicy ParseSweepPolicy(std::string_view s) {
  auto eq = sjs::detail::EqualsIgnoreCase;
  if (eq(s, "auto") || eq(s, "hist") || eq(s, "autohist")) return SweepAxisPolicy::AutoHist;
  if (eq(s, "fixed") || eq(s, "x") || eq(s, "y") || eq(s, "0") || eq(s, "1")) return SweepAxisPolicy::Fixed;
  return SweepAxisPolicy::Orthogonal;
}

template <int Dim, class T>
inline PBSMParams ReadPBSMParams(const Dataset<Dim, T>& ds, const Config& cfg) {
  (void)ds;
  PBSMParams p;
  const auto& extra = cfg.run.extra;

  // Partition scheme
  p.scheme = ParseScheme(GetStringOr(extra, "pbsm_scheme", "stripes"));

  // Partitions
  p.K = GetU64Or(extra, "pbsm_k", 0);
  p.part_axis = GetI32Or(extra, "pbsm_part_axis", 0);
  p.max_partitions = GetU64Or(extra, "pbsm_max_partitions", p.max_partitions);

  // Sweep axis policy
  const std::string sweep_s = GetStringOr(extra, "pbsm_sweep", "orthogonal");
  p.sweep_policy = ParseSweepPolicy(sweep_s);
  if (p.sweep_policy == SweepAxisPolicy::Fixed) {
    // Accept pbsm_sweep as "x"/"y" too.
    if (sjs::detail::EqualsIgnoreCase(sweep_s, "x") || sjs::detail::EqualsIgnoreCase(sweep_s, "0")) {
      p.fixed_sweep_axis = 0;
    } else if (sjs::detail::EqualsIgnoreCase(sweep_s, "y") || sjs::detail::EqualsIgnoreCase(sweep_s, "1")) {
      p.fixed_sweep_axis = 1;
    } else {
      p.fixed_sweep_axis = GetI32Or(extra, "pbsm_sweep_axis", 1);
    }
  } else {
    p.fixed_sweep_axis = GetI32Or(extra, "pbsm_sweep_axis", 1);
  }

  p.hist_bins = GetI32Or(extra, "pbsm_hist_bins", p.hist_bins);
  p.hist_phi = GetI32Or(extra, "pbsm_hist_phi", p.hist_phi);
  p.hist_use_domain = GetBoolOr(extra, "pbsm_hist_use_domain", p.hist_use_domain);

  // Sanity clamps
  if (p.hist_bins < 1) p.hist_bins = 1;
  if (p.hist_bins > 4096) p.hist_bins = 4096;
  if (p.hist_phi < 1) p.hist_phi = 1;

  // Clamp axes
  if (p.part_axis < 0) p.part_axis = 0;
  if (p.part_axis >= Dim) p.part_axis = Dim - 1;
  if (p.fixed_sweep_axis < 0) p.fixed_sweep_axis = 0;
  if (p.fixed_sweep_axis >= Dim) p.fixed_sweep_axis = Dim - 1;
  if (p.fixed_sweep_axis == p.part_axis && p.sweep_policy == SweepAxisPolicy::Orthogonal) {
    // For orthogonal policy, try to pick a different axis if possible.
    if (Dim >= 2) p.fixed_sweep_axis = (p.part_axis == 0) ? 1 : 0;
  }

  return p;
}

// --------------------------
// Partition representation
// --------------------------
// We store partitions as a bounding box region (used for duplicate tests) and
// two sorted index lists (R and S) belonging to that partition.

template <int Dim, class T>
struct Partition {
  // Partition region bounds (half-open).
  Box<Dim, T> region;

  // Sweep axis chosen for this partition.
  int sweep_axis = 0;

  // Sorted indices into ds.R.boxes / ds.S.boxes.
  std::vector<u32> r_idx;
  std::vector<u32> s_idx;
};

// --------------------------
// Utility: auto-pick K using "tile width ~ 10x avg object width" heuristic.
// --------------------------

template <int Dim, class T>
inline u64 AutoPickKStripes(const Dataset<Dim, T>& ds, int part_axis, u64 max_partitions) {
  const auto dom = ds.Domain();
  const T lo = dom.lo.v[static_cast<usize>(part_axis)];
  const T hi = dom.hi.v[static_cast<usize>(part_axis)];
  const T range = hi - lo;
  if (!(range > T(0))) return 1;

  long double sum_w = 0.0L;
  u64 cnt = 0;
  auto add_rel = [&](const Relation<Dim, T>& rel) {
    for (const auto& b : rel.boxes) {
      const T w = b.hi.v[static_cast<usize>(part_axis)] - b.lo.v[static_cast<usize>(part_axis)];
      if (w > T(0)) {
        sum_w += static_cast<long double>(w);
        ++cnt;
      }
    }
  };
  add_rel(ds.R);
  add_rel(ds.S);
  if (cnt == 0) return 1;

  const long double avg_w = sum_w / static_cast<long double>(cnt);
  if (!(avg_w > 0.0L)) return 1;

  const long double target = 10.0L * avg_w;
  if (!(target > 0.0L)) return 1;

  long double K = std::ceil(static_cast<long double>(range) / target);
  if (!(K >= 1.0L)) K = 1.0L;
  if (K > static_cast<long double>(max_partitions)) K = static_cast<long double>(max_partitions);
  return static_cast<u64>(K);
}

template <int Dim, class T>
inline u64 AutoPickKGrid2D(const Dataset<Dim, T>& ds, u64 max_partitions) {
  static_assert(Dim >= 2, "AutoPickKGrid2D requires Dim>=2");
  const auto dom = ds.Domain();
  const T rx = dom.hi.v[0] - dom.lo.v[0];
  const T ry = dom.hi.v[1] - dom.lo.v[1];
  if (!(rx > T(0)) || !(ry > T(0))) return 1;

  // Estimate typical "linear" size per axis.
  long double sum_wx = 0.0L, sum_wy = 0.0L;
  u64 cnt = 0;
  auto add_rel = [&](const Relation<Dim, T>& rel) {
    for (const auto& b : rel.boxes) {
      const T wx = b.hi.v[0] - b.lo.v[0];
      const T wy = b.hi.v[1] - b.lo.v[1];
      if (wx > T(0) && wy > T(0)) {
        sum_wx += static_cast<long double>(wx);
        sum_wy += static_cast<long double>(wy);
        ++cnt;
      }
    }
  };
  add_rel(ds.R);
  add_rel(ds.S);
  if (cnt == 0) return 1;
  const long double ax = sum_wx / static_cast<long double>(cnt);
  const long double ay = sum_wy / static_cast<long double>(cnt);
  if (!(ax > 0.0L) || !(ay > 0.0L)) return 1;

  // Rule of thumb: tile width per axis ~ 10x avg object width along that axis.
  const long double tx = 10.0L * ax;
  const long double ty = 10.0L * ay;
  if (!(tx > 0.0L) || !(ty > 0.0L)) return 1;

  const long double Kx = std::ceil(static_cast<long double>(rx) / tx);
  const long double Ky = std::ceil(static_cast<long double>(ry) / ty);
  long double K = std::max(Kx, Ky);
  if (!(K >= 1.0L)) K = 1.0L;

  // Cap K so that total tiles K^2 does not exceed max_partitions.
  const long double K_cap = std::floor(std::sqrt(static_cast<long double>(max_partitions)));
  if (K > K_cap) K = K_cap;
  if (!(K >= 1.0L)) K = 1.0L;
  return static_cast<u64>(K);
}

// --------------------------
// Partition builder
// --------------------------

inline u32 ClampU32(i64 x, u32 lo, u32 hi) {
  if (x < static_cast<i64>(lo)) return lo;
  if (x > static_cast<i64>(hi)) return hi;
  return static_cast<u32>(x);
}

// Compute [start,end] stripe indices overlapped by [lo,hi) in a 1D partition.
//
// - Domain is [dom_lo, dom_hi) (dom_hi is a max-hi; rectangles satisfy lo < hi <= dom_hi).
// - There are K stripes with width w = (dom_hi - dom_lo)/K, except last stripe hi == dom_hi.
//
// We use nextafter(hi, -inf) so that if hi is exactly on a stripe boundary,
// we do NOT include the stripe starting at hi (half-open).
template <class T>
inline std::pair<u32, u32> IntervalToStripeRange(T lo,
                                                 T hi,
                                                 T dom_lo,
                                                 T dom_hi,
                                                 u32 K) {
  if (K == 0) return {0U, 0U};
  if (!(dom_hi > dom_lo)) return {0U, 0U};
  const T len = dom_hi - dom_lo;
  const T w = len / static_cast<T>(K);
  if (!(w > T(0))) return {0U, 0U};

  // Clamp to domain.
  if (lo < dom_lo) lo = dom_lo;
  if (hi > dom_hi) hi = dom_hi;

  const long double invw = 1.0L / static_cast<long double>(w);
  const long double a = (static_cast<long double>(lo) - static_cast<long double>(dom_lo)) * invw;
  const long double b = (static_cast<long double>(std::nextafter(hi, std::numeric_limits<T>::lowest())) -
                         static_cast<long double>(dom_lo)) * invw;
  i64 i0 = static_cast<i64>(std::floor(a));
  i64 i1 = static_cast<i64>(std::floor(b));
  const u32 start = ClampU32(i0, 0U, K - 1);
  const u32 end = ClampU32(i1, 0U, K - 1);
  return {start, end};
}

// Histogram cost model: estimate candidate comparisons for a given axis.
// Implements a simple version of the I_T^axis score in Tsitsigkos'19.
//
// For each relation list (indices in this partition), we build a binned
// coverage histogram along the given axis:
//   hist[b] = #rectangles whose projection overlaps bin b.
// Then score = sum_b histR[b] * histS[b]. Lower is better.
//
// To keep overhead low, we optionally subsample by taking every phi-th item.
// The score is still comparable between axes for choosing the min.
template <int Dim, class T>
inline u64 HistScoreAxis(const Relation<Dim, T>& R,
                         const Relation<Dim, T>& S,
                         const std::vector<u32>& r_idx,
                         const std::vector<u32>& s_idx,
                         int axis,
                         const Box<Dim, T>& range_box,
                         int bins,
                         int phi) {
  SJS_DASSERT(axis >= 0 && axis < Dim);
  if (bins <= 0) return 0;
  if (phi <= 0) phi = 1;

  const T lo = range_box.lo.v[static_cast<usize>(axis)];
  const T hi = range_box.hi.v[static_cast<usize>(axis)];
  const T len = hi - lo;
  if (!(len > T(0))) return 0;

  std::vector<u32> hR(static_cast<usize>(bins), 0U);
  std::vector<u32> hS(static_cast<usize>(bins), 0U);

  auto add_interval = [&](std::vector<u32>& h, T a, T b) {
    // clamp
    if (a < lo) a = lo;
    if (b > hi) b = hi;
    if (!(a < b)) return;

    const long double inv = static_cast<long double>(bins) / static_cast<long double>(len);
    const long double fa = (static_cast<long double>(a) - static_cast<long double>(lo)) * inv;
    const long double fb = (static_cast<long double>(std::nextafter(b, std::numeric_limits<T>::lowest())) -
                            static_cast<long double>(lo)) * inv;

    i64 i0 = static_cast<i64>(std::floor(fa));
    i64 i1 = static_cast<i64>(std::floor(fb));
    if (i0 < 0) i0 = 0;
    if (i1 < 0) i1 = 0;
    if (i0 >= bins) i0 = bins - 1;
    if (i1 >= bins) i1 = bins - 1;
    for (i64 i = i0; i <= i1; ++i) {
      h[static_cast<usize>(i)]++;
    }
  };

  // R
  for (usize i = 0; i < r_idx.size(); i += static_cast<usize>(phi)) {
    const u32 idx = r_idx[i];
    const auto& b = R.boxes[static_cast<usize>(idx)];
    add_interval(hR, b.lo.v[static_cast<usize>(axis)], b.hi.v[static_cast<usize>(axis)]);
  }
  // S
  for (usize i = 0; i < s_idx.size(); i += static_cast<usize>(phi)) {
    const u32 idx = s_idx[i];
    const auto& b = S.boxes[static_cast<usize>(idx)];
    add_interval(hS, b.lo.v[static_cast<usize>(axis)], b.hi.v[static_cast<usize>(axis)]);
  }

  // score
  unsigned long long score = 0ULL;
  for (int i = 0; i < bins; ++i) {
    score += static_cast<unsigned long long>(hR[static_cast<usize>(i)]) *
             static_cast<unsigned long long>(hS[static_cast<usize>(i)]);
  }
  return static_cast<u64>(score);
}

// Choose sweep axis for this partition.
template <int Dim, class T>
inline int ChooseSweepAxis(const Relation<Dim, T>& R,
                           const Relation<Dim, T>& S,
                           const Partition<Dim, T>& P,
                           const PBSMParams& params,
                           int part_axis,
                           const Box<Dim, T>& domain_box) {
  if (params.sweep_policy == SweepAxisPolicy::Fixed) {
    return params.fixed_sweep_axis;
  }
  if (params.sweep_policy == SweepAxisPolicy::Orthogonal) {
    if (Dim >= 2) return (part_axis == 0) ? 1 : 0;
    return 0;
  }

  // AutoHist: evaluate all axes (but for 2D we only compare 0 vs 1).
  // For Dim>2, we can still choose among all axes.
  const Box<Dim, T>& range_box = params.hist_use_domain ? domain_box : P.region;

  u64 best_score = std::numeric_limits<u64>::max();
  int best_axis = 0;
  for (int ax = 0; ax < Dim; ++ax) {
    const u64 sc = HistScoreAxis<Dim, T>(R, S, P.r_idx, P.s_idx, ax, range_box,
                                        params.hist_bins, params.hist_phi);
    if (sc < best_score) {
      best_score = sc;
      best_axis = ax;
    }
  }
  return best_axis;
}

// Build partitions and per-partition sorted lists.
template <int Dim, class T>
class PBSMIndex {
 public:
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;

  void Reset() {
    ds_ = nullptr;
    built_ = false;
    params_ = PBSMParams{};
    domain_ = BoxT::Empty();
    partitions_.clear();
  }

  bool Built() const noexcept { return built_; }
  const DatasetT* DatasetPtr() const noexcept { return ds_; }
  const PBSMParams& Params() const noexcept { return params_; }
  const BoxT& Domain() const noexcept { return domain_; }
  const std::vector<Partition<Dim, T>>& Partitions() const noexcept { return partitions_; }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) {
    auto scoped = phases ? phases->Scoped("pbsm_build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = false;
    partitions_.clear();
    domain_ = ds.Domain();
    params_ = ReadPBSMParams<Dim, T>(ds, cfg);

    if (ds.R.Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        ds.S.Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "PBSMIndex::Build: relation size exceeds u32 index limit";
      return false;
    }

    // Decide K if needed.
    if (params_.K == 0) {
      if (params_.scheme == PartitionScheme::Grid2D && Dim >= 2) {
        params_.K = AutoPickKGrid2D<Dim, T>(ds, params_.max_partitions);
      } else {
        params_.K = AutoPickKStripes<Dim, T>(ds, params_.part_axis, params_.max_partitions);
      }
    }
    if (params_.K == 0) params_.K = 1;

    // Cap K to avoid explosive memory usage.
    if (params_.scheme == PartitionScheme::Grid2D && Dim >= 2) {
      // total tiles = K^2
      const long double kk = static_cast<long double>(params_.K);
      const long double total = kk * kk;
      if (total > static_cast<long double>(params_.max_partitions)) {
        const u64 Kcap = static_cast<u64>(std::floor(std::sqrt(static_cast<long double>(params_.max_partitions))));
        params_.K = (Kcap == 0) ? 1 : Kcap;
      }
    } else {
      if (params_.K > params_.max_partitions) params_.K = params_.max_partitions;
    }
    if (params_.K == 0) params_.K = 1;

    // Build partitions.
    if (params_.scheme == PartitionScheme::Grid2D) {
      if constexpr (Dim == 2) {
        return BuildGrid2D(phases, err);
      } else {
        // Fallback to stripes for Dim != 2.
        params_.scheme = PartitionScheme::Stripes1D;
        return BuildStripes1D(phases, err);
      }
    }
    return BuildStripes1D(phases, err);
  }

 private:
  bool BuildStripes1D(PhaseRecorder* phases, std::string* err) {
    auto scoped = phases ? phases->Scoped("pbsm_partitions_stripes") : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 K64 = params_.K;
    if (K64 == 0) {
      if (err) *err = "PBSMIndex::BuildStripes1D: K=0";
      return false;
    }
    if (K64 > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "PBSMIndex::BuildStripes1D: K too large";
      return false;
    }
    const u32 K = static_cast<u32>(K64);
    const int a = params_.part_axis;
    if (a < 0 || a >= Dim) {
      if (err) *err = "PBSMIndex::BuildStripes1D: invalid part_axis";
      return false;
    }

    const T dom_lo = domain_.lo.v[static_cast<usize>(a)];
    const T dom_hi = domain_.hi.v[static_cast<usize>(a)];
    const T len = dom_hi - dom_lo;
    if (!(len > T(0))) {
      // Degenerate domain: still make one partition.
      partitions_.resize(1);
      partitions_[0].region = domain_;
      partitions_[0].sweep_axis = (Dim >= 2) ? ((a == 0) ? 1 : 0) : 0;
      FillPartitionListsStripes(/*K=*/1, /*a=*/a, /*dom_lo=*/dom_lo, /*dom_hi=*/dom_hi);
      SortAndChooseAxesStripes(a);
      built_ = true;
      return true;
    }

    const T w = len / static_cast<T>(K);
    if (!(w > T(0))) {
      // Too many stripes for tiny range; fallback.
      partitions_.resize(1);
      partitions_[0].region = domain_;
      partitions_[0].sweep_axis = (Dim >= 2) ? ((a == 0) ? 1 : 0) : 0;
      FillPartitionListsStripes(/*K=*/1, /*a=*/a, /*dom_lo=*/dom_lo, /*dom_hi=*/dom_hi);
      SortAndChooseAxesStripes(a);
      built_ = true;
      return true;
    }

    partitions_.resize(static_cast<usize>(K));
    // Initialize partition regions.
    for (u32 i = 0; i < K; ++i) {
      Partition<Dim, T>& P = partitions_[static_cast<usize>(i)];
      P.region = domain_;  // start with full domain
      const T lo = dom_lo + static_cast<T>(i) * w;
      const T hi = (i + 1 == K) ? dom_hi : (dom_lo + static_cast<T>(i + 1) * w);
      P.region.lo.v[static_cast<usize>(a)] = lo;
      P.region.hi.v[static_cast<usize>(a)] = hi;
      // sweep axis filled later.
      P.sweep_axis = 0;
      P.r_idx.clear();
      P.s_idx.clear();
    }

    FillPartitionListsStripes(K, a, dom_lo, dom_hi);
    SortAndChooseAxesStripes(a);
    built_ = true;
    return true;
  }

  void FillPartitionListsStripes(u32 K, int a, T dom_lo, T dom_hi) {
    SJS_DASSERT(K >= 1);

    auto add_rel = [&](const Relation<Dim, T>& rel, bool is_R) {
      for (usize i = 0; i < rel.boxes.size(); ++i) {
        const auto& b = rel.boxes[i];
        const T lo = b.lo.v[static_cast<usize>(a)];
        const T hi = b.hi.v[static_cast<usize>(a)];
        if (!(lo < hi)) continue;
        const auto range = IntervalToStripeRange<T>(lo, hi, dom_lo, dom_hi, K);
        for (u32 p = range.first; p <= range.second; ++p) {
          auto& part = partitions_[static_cast<usize>(p)];
          if (is_R) part.r_idx.push_back(static_cast<u32>(i));
          else part.s_idx.push_back(static_cast<u32>(i));
        }
      }
    };
    add_rel(ds_->R, /*is_R=*/true);
    add_rel(ds_->S, /*is_R=*/false);
  }

  void SortAndChooseAxesStripes(int part_axis) {
    // Choose sweep axis (per partition if AutoHist).
    for (auto& P : partitions_) {
      P.sweep_axis = ChooseSweepAxis<Dim, T>(ds_->R, ds_->S, P, params_, part_axis, domain_);

      auto sort_rel = [&](const Relation<Dim, T>& rel, std::vector<u32>& idxs) {
        const int ax = P.sweep_axis;
        auto cmp = [&](u32 ia, u32 ib) {
          const auto& a = rel.boxes[static_cast<usize>(ia)];
          const auto& b = rel.boxes[static_cast<usize>(ib)];
          const T la = a.lo.v[static_cast<usize>(ax)];
          const T lb = b.lo.v[static_cast<usize>(ax)];
          if (la < lb) return true;
          if (lb < la) return false;
          const Id ida = rel.GetId(static_cast<usize>(ia));
          const Id idb = rel.GetId(static_cast<usize>(ib));
          if (ida < idb) return true;
          if (idb < ida) return false;
          return ia < ib;
        };
        std::sort(idxs.begin(), idxs.end(), cmp);
      };

      sort_rel(ds_->R, P.r_idx);
      sort_rel(ds_->S, P.s_idx);
    }
  }

  bool BuildGrid2D(PhaseRecorder* phases, std::string* err) {
    static_assert(Dim == 2, "BuildGrid2D only for Dim==2");
    auto scoped = phases ? phases->Scoped("pbsm_partitions_grid2d") : PhaseRecorder::ScopedPhase(nullptr, "");

    const u64 K64 = params_.K;
    if (K64 == 0) {
      if (err) *err = "PBSMIndex::BuildGrid2D: K=0";
      return false;
    }
    if (K64 > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "PBSMIndex::BuildGrid2D: K too large";
      return false;
    }
    const u32 K = static_cast<u32>(K64);
    const T x0 = domain_.lo.v[0];
    const T x1 = domain_.hi.v[0];
    const T y0 = domain_.lo.v[1];
    const T y1 = domain_.hi.v[1];
    const T rx = x1 - x0;
    const T ry = y1 - y0;
    if (!(rx > T(0)) || !(ry > T(0))) {
      // Degenerate: fall back to a single tile.
      partitions_.resize(1);
      partitions_[0].region = domain_;
      partitions_[0].sweep_axis = 0;
      FillPartitionListsGrid2D(/*K=*/1, x0, x1, y0, y1);
      SortAndChooseAxesGrid2D();
      built_ = true;
      return true;
    }
    const T wx = rx / static_cast<T>(K);
    const T wy = ry / static_cast<T>(K);
    if (!(wx > T(0)) || !(wy > T(0))) {
      partitions_.resize(1);
      partitions_[0].region = domain_;
      partitions_[0].sweep_axis = 0;
      FillPartitionListsGrid2D(/*K=*/1, x0, x1, y0, y1);
      SortAndChooseAxesGrid2D();
      built_ = true;
      return true;
    }

    partitions_.resize(static_cast<usize>(K) * static_cast<usize>(K));

    // Init tile regions in row-major (y-major) order.
    for (u32 iy = 0; iy < K; ++iy) {
      for (u32 ix = 0; ix < K; ++ix) {
        const usize id = static_cast<usize>(iy) * static_cast<usize>(K) + static_cast<usize>(ix);
        auto& P = partitions_[id];
        P.region = domain_;
        const T xl = x0 + static_cast<T>(ix) * wx;
        const T xh = (ix + 1 == K) ? x1 : (x0 + static_cast<T>(ix + 1) * wx);
        const T yl = y0 + static_cast<T>(iy) * wy;
        const T yh = (iy + 1 == K) ? y1 : (y0 + static_cast<T>(iy + 1) * wy);
        P.region.lo.v[0] = xl;
        P.region.hi.v[0] = xh;
        P.region.lo.v[1] = yl;
        P.region.hi.v[1] = yh;
        P.sweep_axis = 0;
        P.r_idx.clear();
        P.s_idx.clear();
      }
    }

    FillPartitionListsGrid2D(K, x0, x1, y0, y1);
    SortAndChooseAxesGrid2D();
    built_ = true;
    return true;
  }

  void FillPartitionListsGrid2D(u32 K, T x0, T x1, T y0, T y1) {
    SJS_DASSERT(K >= 1);
    auto add_rel = [&](const Relation<Dim, T>& rel, bool is_R) {
      for (usize i = 0; i < rel.boxes.size(); ++i) {
        const auto& b = rel.boxes[i];
        if (b.IsEmpty()) continue;
        const auto xr = IntervalToStripeRange<T>(b.lo.v[0], b.hi.v[0], x0, x1, K);
        const auto yr = IntervalToStripeRange<T>(b.lo.v[1], b.hi.v[1], y0, y1, K);
        for (u32 iy = yr.first; iy <= yr.second; ++iy) {
          for (u32 ix = xr.first; ix <= xr.second; ++ix) {
            const usize pid = static_cast<usize>(iy) * static_cast<usize>(K) + static_cast<usize>(ix);
            auto& P = partitions_[pid];
            if (is_R) P.r_idx.push_back(static_cast<u32>(i));
            else P.s_idx.push_back(static_cast<u32>(i));
          }
        }
      }
    };
    add_rel(ds_->R, /*is_R=*/true);
    add_rel(ds_->S, /*is_R=*/false);
  }

  void SortAndChooseAxesGrid2D() {
    // Choose sweep axis per partition if AutoHist, else use fixed/orthogonal.
    // For grid, "orthogonal" is treated as fixed_sweep_axis (default 0) unless AutoHist.
    for (auto& P : partitions_) {
      int part_axis_dummy = 0;
      if (params_.sweep_policy == SweepAxisPolicy::Orthogonal) {
        // For grid, orthogonal has no obvious meaning; keep fixed_sweep_axis.
        P.sweep_axis = params_.fixed_sweep_axis;
      } else {
        P.sweep_axis = ChooseSweepAxis<Dim, T>(ds_->R, ds_->S, P, params_, part_axis_dummy, domain_);
      }

      auto sort_rel = [&](const Relation<Dim, T>& rel, std::vector<u32>& idxs) {
        const int ax = P.sweep_axis;
        auto cmp = [&](u32 ia, u32 ib) {
          const auto& a = rel.boxes[static_cast<usize>(ia)];
          const auto& b = rel.boxes[static_cast<usize>(ib)];
          const T la = a.lo.v[static_cast<usize>(ax)];
          const T lb = b.lo.v[static_cast<usize>(ax)];
          if (la < lb) return true;
          if (lb < la) return false;
          const Id ida = rel.GetId(static_cast<usize>(ia));
          const Id idb = rel.GetId(static_cast<usize>(ib));
          if (ida < idb) return true;
          if (idb < ida) return false;
          return ia < ib;
        };
        std::sort(idxs.begin(), idxs.end(), cmp);
      };

      sort_rel(ds_->R, P.r_idx);
      sort_rel(ds_->S, P.s_idx);
    }
  }

 private:
  const DatasetT* ds_ = nullptr;
  bool built_ = false;
  PBSMParams params_{};
  BoxT domain_ = BoxT::Empty();
  std::vector<Partition<Dim, T>> partitions_;
};

// --------------------------
// Duplicate elimination (reference point)
// --------------------------
// Reference point (Baseline §3.5 / Lemma 3):
//   p(r,s) = ( max(x_l(r), x_l(s)),  max(y_l(r), y_l(s)) ).
//
// 2D tiles (Eq.(1) in Tsitsigkos'19 / Baseline):
//   report/count in tile T iff  p_x >= T.x_l  &&  p_y >= T.y_l.
// 1D stripes (Tsitsigkos'19 simplification / Baseline):
//   report/count in stripe T iff  p_axis >= T.axis_l  (only check the partition axis).
//
// Note: we intentionally do NOT check the upper bounds (p_x < T.x_u, p_y < T.y_u) here;
// if a pair is discovered inside a partition, multi-assignment implies both rectangles overlap
// that partition, which in turn implies the upper-bound conditions automatically.

template <int Dim, class T>
inline bool DuplicateTestStripes(const Box<Dim, T>& r,
                                 const Box<Dim, T>& s,
                                 const Box<Dim, T>& stripe_region,
                                 int part_axis) {
  const usize ax = static_cast<usize>(part_axis);
  const T rp = (r.lo.v[ax] > s.lo.v[ax]) ? r.lo.v[ax] : s.lo.v[ax];
  const T lo = stripe_region.lo.v[ax];
  // Baseline: only lower-bound test on the partition axis.
  return (rp >= lo);
}

template <class T>
inline bool DuplicateTestGrid2D(const Box<2, T>& r, const Box<2, T>& s, const Box<2, T>& tile_region) {
  const T rpx = (r.lo.v[0] > s.lo.v[0]) ? r.lo.v[0] : s.lo.v[0];
  const T rpy = (r.lo.v[1] > s.lo.v[1]) ? r.lo.v[1] : s.lo.v[1];
  // Baseline Eq.(1): only lower-bound test on x/y.
  return (rpx >= tile_region.lo.v[0]) && (rpy >= tile_region.lo.v[1]);
}

// --------------------------
// Forward-scan join enumerator per partition (streaming)
// --------------------------

// Helper: intersection check assuming overlap on sweep axis.
template <int Dim, class T>
inline bool IntersectsOtherAxes(const Box<Dim, T>& a, const Box<Dim, T>& b, int sweep_axis) {
  // If either box is empty, treat as no intersection.
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

// Deterministic join stream over all partitions.
//
// Stream order:
//   - partitions in stored order (stripes increasing index, grid row-major)
//   - within a partition: forward-scan as in Tsitsigkos’19 Algorithm 1 style
//     (ties resolved by the fixed "if lower(R) < lower(S) else" rule)
template <int Dim, class T>
class PBSMJoinEnumerator final : public baselines::IJoinEnumerator {
 public:
  using DatasetT = Dataset<Dim, T>;

  PBSMJoinEnumerator(const PBSMIndex<Dim, T>* index)
      : index_(index) {
    SJS_ASSERT(index_ != nullptr);
    Reset();
  }

  void Reset() override {
    stats_.Reset();
    part_pos_ = 0;
    i_ = 0;
    j_ = 0;
    scanning_ = false;
    fixed_side_r_ = true;
    fixed_idx_ = 0;
    scan_pos_ = 0;

    // Precompute deterministic total "events" proxy: number of assigned objects.
    // Not a real sweep event count, but provides a stable instrumentation signal.
    u64 n_events = 0;
    if (index_) {
      for (const auto& P : index_->Partitions()) {
        n_events += static_cast<u64>(P.r_idx.size()) + static_cast<u64>(P.s_idx.size());
      }
    }
    stats_.num_events = n_events;
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

  bool Next(PairId* out) override {
    if (!out) return false;
    if (!index_ || !index_->Built() || !index_->DatasetPtr()) return false;

    const auto& ds = *index_->DatasetPtr();
    const auto& parts = index_->Partitions();
    const auto& params = index_->Params();

    while (part_pos_ < parts.size()) {
      const auto& P = parts[part_pos_];
      const int sweep_axis = P.sweep_axis;
      const int other_axis = (Dim >= 2) ? ((sweep_axis == 0) ? 1 : 0) : 0;
      (void)other_axis;

      // Skip empty partitions quickly.
      if (P.r_idx.empty() || P.s_idx.empty()) {
        AdvancePartition();
        continue;
      }

      // Main loop within partition.
      while (true) {
        if (scanning_) {
          if (fixed_side_r_) {
            // Fixed R object; scanning S from scan_pos_.
            if (fixed_idx_ >= P.r_idx.size()) {
              // Should not happen, but recover.
              scanning_ = false;
              continue;
            }
            const u32 r_i = P.r_idx[fixed_idx_];
            const auto& rb = ds.R.boxes[static_cast<usize>(r_i)];
            const T r_hi = rb.hi.v[static_cast<usize>(sweep_axis)];

            if (scan_pos_ >= P.s_idx.size()) {
              // Done scanning this R; advance i.
              scanning_ = false;
              ++i_;
              continue;
            }
            const u32 s_j = P.s_idx[scan_pos_];
            const auto& sb = ds.S.boxes[static_cast<usize>(s_j)];
            const T s_lo = sb.lo.v[static_cast<usize>(sweep_axis)];

            // Stop condition: s_lo >= r_hi => no further S can overlap this R on sweep axis.
            if (!(s_lo < r_hi)) {
              scanning_ = false;
              ++i_;
              continue;
            }

            // Candidate pair (r_i, s_j)
            ++scan_pos_;
            stats_.candidate_checks++;
            if (!IntersectsOtherAxes<Dim, T>(rb, sb, sweep_axis)) {
              continue;
            }

            // Duplicate elimination.
            bool ok_dup = true;
            if (params.scheme == PartitionScheme::Stripes1D) {
              ok_dup = DuplicateTestStripes<Dim, T>(rb, sb, P.region, params.part_axis);
            } else {
              if constexpr (Dim == 2) {
                ok_dup = DuplicateTestGrid2D<T>(rb, sb, P.region);
              }
            }
            if (!ok_dup) continue;

            *out = PairId{ds.R.GetId(static_cast<usize>(r_i)), ds.S.GetId(static_cast<usize>(s_j))};
            stats_.output_pairs++;
            return true;
          } else {
            // Fixed S object; scanning R from scan_pos_.
            if (fixed_idx_ >= P.s_idx.size()) {
              scanning_ = false;
              continue;
            }
            const u32 s_j = P.s_idx[fixed_idx_];
            const auto& sb = ds.S.boxes[static_cast<usize>(s_j)];
            const T s_hi = sb.hi.v[static_cast<usize>(sweep_axis)];

            if (scan_pos_ >= P.r_idx.size()) {
              scanning_ = false;
              ++j_;
              continue;
            }
            const u32 r_i = P.r_idx[scan_pos_];
            const auto& rb = ds.R.boxes[static_cast<usize>(r_i)];
            const T r_lo = rb.lo.v[static_cast<usize>(sweep_axis)];
            if (!(r_lo < s_hi)) {
              scanning_ = false;
              ++j_;
              continue;
            }

            ++scan_pos_;
            stats_.candidate_checks++;
            if (!IntersectsOtherAxes<Dim, T>(rb, sb, sweep_axis)) {
              continue;
            }

            bool ok_dup = true;
            if (params.scheme == PartitionScheme::Stripes1D) {
              ok_dup = DuplicateTestStripes<Dim, T>(rb, sb, P.region, params.part_axis);
            } else {
              if constexpr (Dim == 2) {
                ok_dup = DuplicateTestGrid2D<T>(rb, sb, P.region);
              }
            }
            if (!ok_dup) continue;

            *out = PairId{ds.R.GetId(static_cast<usize>(r_i)), ds.S.GetId(static_cast<usize>(s_j))};
            stats_.output_pairs++;
            return true;
          }
        }

        // Not currently scanning: choose next forward-scan step.
        if (i_ >= P.r_idx.size() || j_ >= P.s_idx.size()) {
          // Partition finished.
          AdvancePartition();
          break;
        }

        const u32 r_i = P.r_idx[i_];
        const u32 s_j = P.s_idx[j_];
        const auto& rb = ds.R.boxes[static_cast<usize>(r_i)];
        const auto& sb = ds.S.boxes[static_cast<usize>(s_j)];
        const T r_lo = rb.lo.v[static_cast<usize>(sweep_axis)];
        const T s_lo = sb.lo.v[static_cast<usize>(sweep_axis)];

        // Tsitsigkos-style merge decision: if lower(R) < lower(S) then scan forward in S,
        // else scan forward in R. This is deterministic.
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
        // Loop continues; the scanning branch will produce (or skip) candidates.
      }
    }

    return false;
  }

 private:
  void AdvancePartition() {
    ++part_pos_;
    i_ = 0;
    j_ = 0;
    scanning_ = false;
    fixed_side_r_ = true;
    fixed_idx_ = 0;
    scan_pos_ = 0;
  }

  const PBSMIndex<Dim, T>* index_ = nullptr;

  // Global iteration state
  usize part_pos_ = 0;

  // Per-partition state
  usize i_ = 0;  // pointer into R list
  usize j_ = 0;  // pointer into S list

  // Current scanning state
  bool scanning_ = false;
  bool fixed_side_r_ = true;  // true => fixed is R[i_]; false => fixed is S[j_]
  usize fixed_idx_ = 0;       // fixed index within the fixed side list (i_ or j_ snapshot)
  usize scan_pos_ = 0;        // scanning position in the opposite list

  join::JoinStats stats_{};
};

}  // namespace detail

// --------------------------
// PBSM Sampling baseline
// --------------------------

template <int Dim, class T = Scalar>
class PBSMSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "PBSMSamplingBaseline requires Dim>=2");
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::PBSM; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "pbsm_sampling"; }

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
      if (err) *err = "PBSMSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "PBSMSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("pbsm_count_pass1") : PhaseRecorder::ScopedPhase(nullptr, "");

    detail::PBSMJoinEnumerator<Dim, T> stream(&index_);
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
      if (err) *err = "PBSMSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "PBSMSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "PBSMSamplingBaseline::Sample: out is null";
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
      auto scoped2 = phases ? phases->Scoped("pbsm_sample_draw_ranks") : PhaseRecorder::ScopedPhase(nullptr, "");
      for (u64 i = 0; i < t; ++i) {
        req[static_cast<usize>(i)] = RankReq{rng->UniformU64(W), i};
      }
      std::sort(req.begin(), req.end(), [](const RankReq& a, const RankReq& b) {
        if (a.rank < b.rank) return true;
        if (b.rank < a.rank) return false;
        return a.slot < b.slot;
      });
    }

    // Pass2: scan join stream and pick pairs at requested ranks.
    out->pairs.assign(static_cast<usize>(t), PairId{});
    {
      auto scoped2 = phases ? phases->Scoped("pbsm_sample_pass2_select") : PhaseRecorder::ScopedPhase(nullptr, "");
      detail::PBSMJoinEnumerator<Dim, T> stream(&index_);
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
        if (err) *err = "PBSMSamplingBaseline::Sample: stream ended early in pass2 (cardinality mismatch: expected W pairs from EnumerateUniquePairs)";
        return false;
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("pbsm_enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !index_.Built()) {
      if (err) *err = "PBSMSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::PBSMJoinEnumerator<Dim, T>>(&index_);
  }

 private:
  detail::PBSMIndex<Dim, T> index_;
  bool built_ = false;

  // cached count from last Count() call
  bool have_count_ = false;
  u64 W_ = 0;
};

}  // namespace pbsm
}  // namespace baselines
}  // namespace sjs
