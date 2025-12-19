#pragma once
// sjs/baselines/rejection/sampling.h
//
// AGR-BoxJoin baseline (Amagata'25-style) — Variant::Sampling (AGR-S)
//
// Goal
// ----
// Output t i.i.d. uniform (with replacement) samples from the true box-intersection
// join result:
//   J = {(r,s) | r in A, s in B, Intersect(r,s)=true}.
//
// Design summary (matches docs/Baseline/Amagata’25 Baseline.md v2)
// ---------------------------------------------------------------
// 1) Half-open boxes: [lo,hi) and strict-intersection predicate.
// 2) Fix representative point for indexed side B:
//      rep(s) := lo(s)   (lower-left / minimum corner)   [fixed; no alternatives here]
// 3) Sparse hash grid on rep(s):
//      key_i(x) = floor(x_i / g_i).
// 4) Upper bound construction:
//      M_i := max_{s in B} (hi_i(s) - lo_i(s))
//      E(r) := Π_i [lo_i(r) - M_i, hi_i(r))
//      C(r) := { k in KeysNonEmpty | cell(k) ∩ E(r) != ∅ }   (cell-intersection definition)
//      μ(r) := Σ_{k in C(r)} |CellMap[k]|
// 5) Two-level alias sampling over the candidate superset U_prop:
//      r ~ AliasR  with weight μ(r)
//      k ~ AliasCell[r] with weight |CellMap[k]|
//      s ~ Uniform(CellMap[k])
//      accept iff Intersect(r,s)
//
// Correctness:
//  - Proposal probability is constant over U_prop: 1 / MuSum, where MuSum=Σ_r μ(r).
//  - Conditioning on accept yields uniform over J.
//  - Independent proposals => i.i.d.
//
// Engineering boundary (robust termination):
//  - If MuSum==0, then |J|==0 (necessary condition) => return empty.
//  - If MuSum>0 but |J|==0 (possible in practice), naive AGR-S would loop forever.
//    We prevent non-termination while preserving correctness by doing a
//    deterministic emptiness check *only if we fail to accept for a while*.
//    If the check proves that |J|==0, we return empty; otherwise we continue
//    with standard AGR-S sampling.
//
//    Knob (run.extra):
//      rej_max_factor : number of rejected proposals allowed before the
//                       deterministic empty-join check (default: 10000).
//
//  - AGR-AS (Adaptive) remains the recommended entry point when you want
//    predictable performance on very sparse joins.
//
// Candidate-cell construction (doc §3.8):
//  - Impl A (Naive): enumerate key hyper-rectangle and hash lookup.
//  - Impl B (Slab, recommended): iterate only non-empty keys by bucketing on k1.
//
// This header supports both; default is Slab.
//
// NOTE: This baseline is intentionally "clean": we do NOT support experimental
// representative points (e.g., center) in this baseline, because the coverage
// proof (Lemma 1) is tied to rep(s)=lo(s).

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/core/config.h"
#include "sjs/core/rng.h"
#include "sjs/geometry/box.h"
#include "sjs/io/dataset.h"
#include "sjs/join/join_types.h"
#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <array>
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
namespace rejection {

namespace detail {

// ---------------------------------
// Extra-map parsing helpers (local)
// ---------------------------------

inline char ToLowerAscii(char c) noexcept {
  if (c >= 'A' && c <= 'Z') return static_cast<char>(c - 'A' + 'a');
  return c;
}

inline bool EqualsIgnoreCase(std::string_view a, std::string_view b) noexcept {
  if (a.size() != b.size()) return false;
  for (usize i = 0; i < a.size(); ++i) {
    if (ToLowerAscii(a[i]) != ToLowerAscii(b[i])) return false;
  }
  return true;
}

inline std::optional<std::string_view> GetExtra(const std::unordered_map<std::string, std::string>& extra,
                                                std::string_view key) {
  auto it = extra.find(std::string(key));
  if (it == extra.end()) return std::nullopt;
  return std::string_view(it->second);
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

inline bool ParseDouble(std::string_view s, double* out) {
  if (!out || s.empty()) return false;
  try {
    std::size_t idx = 0;
    double v = std::stod(std::string(s), &idx);
    if (idx != s.size()) return false;
    *out = v;
    return true;
  } catch (...) {
    return false;
  }
}

inline bool ParseBool(std::string_view s, bool* out) {
  if (!out || s.empty()) return false;
  if (EqualsIgnoreCase(s, "1") || EqualsIgnoreCase(s, "true") || EqualsIgnoreCase(s, "yes") ||
      EqualsIgnoreCase(s, "y") || EqualsIgnoreCase(s, "on")) {
    *out = true;
    return true;
  }
  if (EqualsIgnoreCase(s, "0") || EqualsIgnoreCase(s, "false") || EqualsIgnoreCase(s, "no") ||
      EqualsIgnoreCase(s, "n") || EqualsIgnoreCase(s, "off")) {
    *out = false;
    return true;
  }
  return false;
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

inline double GetDoubleOr(const std::unordered_map<std::string, std::string>& extra,
                          std::string_view key,
                          double def) {
  if (auto v = GetExtra(extra, key)) {
    double x = 0;
    if (ParseDouble(*v, &x)) return x;
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

// ---------------------------------
// Cell key & hashing (sparse grid)
// ---------------------------------

template <int Dim>
struct CellKey {
  // Suppress sign-conversion warning: Dim is compile-time constant and always positive
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsign-conversion"
  std::array<i64, Dim> k{};
  #pragma GCC diagnostic pop
  bool operator==(const CellKey& o) const noexcept { return k == o.k; }
};

template <int Dim>
struct CellKeyHash {
  // splitmix64 finalizer
  static inline u64 Mix64(u64 x) noexcept {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x ^= (x >> 31);
    return x;
  }

  usize operator()(const CellKey<Dim>& key) const noexcept {
    u64 h = 0x243f6a8885a308d3ULL;
    for (int d = 0; d < Dim; ++d) {
      const u64 x = static_cast<u64>(key.k[d]);
      h ^= Mix64(x + static_cast<u64>(d) * 0x9e3779b97f4a7c15ULL);
      h = Mix64(h);
    }
    return static_cast<usize>(h);
  }
};

template <int Dim>
inline bool KeyLexLess(const CellKey<Dim>& a, const CellKey<Dim>& b) noexcept {
  for (int d = 0; d < Dim; ++d) {
    if (a.k[d] < b.k[d]) return true;
    if (b.k[d] < a.k[d]) return false;
  }
  return false;
}

// ---------------------------------
// Numeric helpers for grid indexing
// ---------------------------------

inline i64 ClampToI64(long double v) noexcept {
  const long double lo = static_cast<long double>(std::numeric_limits<i64>::min());
  const long double hi = static_cast<long double>(std::numeric_limits<i64>::max());
  if (v < lo) return std::numeric_limits<i64>::min();
  if (v > hi) return std::numeric_limits<i64>::max();
  return static_cast<i64>(v);
}

// floor(x/g) with mathematical floor (correct for negative x).
inline i64 FloorDiv(long double x, long double g) noexcept { return ClampToI64(std::floor(x / g)); }

// For half-open [*, x): k_max = ceil(x/g) - 1.
inline i64 CeilDivMinus1(long double x, long double g) noexcept { return ClampToI64(std::ceil(x / g) - 1.0L); }

// ---------------------------------
// Parameters (strictly baseline-aligned)
// ---------------------------------

enum class CandidateCellsImpl : u8 {
  // Doc §3.8 Impl B (recommended): iterate non-empty keys via slab index (bucket by k1).
  Slab = 0,
  // Doc §3.8 Impl A: enumerate hyper-rectangle keys and hash lookup.
  Naive = 1,
};

template <int Dim>
struct RejectionParams {
  // Grid cell size per dimension. Must be > 0.
  // Suppress sign-conversion warning: Dim is compile-time constant and always positive
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsign-conversion"
  std::array<double, Dim> g{};
  #pragma GCC diagnostic pop

  // If true: swap sides (index on R and query S). Default false (index on S).
  bool swap_sides = false;

  // Candidate cell enumeration implementation (doc §3.8). Default Slab.
  CandidateCellsImpl cand_impl = CandidateCellsImpl::Slab;

  // Pilot draws for estimating |J| (doc §5.4). 0 disables pilot.
  u64 count_draws = 50'000ULL;
};

// Parse params from cfg (and reject unsupported experimental knobs).
template <int Dim, class T>
inline bool ReadRejectionParams(const Dataset<Dim, T>& ds,
                                const Config& cfg,
                                RejectionParams<Dim>* out,
                                std::string* err) {
  if (!out) {
    if (err) *err = "ReadRejectionParams: out is null";
    return false;
  }
  RejectionParams<Dim> p;
  const auto& extra = cfg.run.extra;

  // Strict baseline: reject experimental representative-point knobs if enabled.
  // (rep(s) must be lo(s) for the coverage proof in the baseline doc.)
  if (auto v = GetExtra(extra, "rej_rep_center")) {
    bool b = false;
    if (ParseBool(*v, &b) && b) {
      if (err) *err = "AGR-BoxJoin baseline: rej_rep_center=true is not supported (baseline fixes rep(s)=lo(s))";
      return false;
    }
  }

  // Side selection
  // Accepted keys:
  //   - rej_swap = true/false
  //   - rej_index_side = "R" or "S" (which relation is indexed B)
  p.swap_sides = GetBoolOr(extra, "rej_swap", false);
  const std::string idx_side = GetStringOr(extra, "rej_index_side", GetStringOr(extra, "rej_index", "S"));
  if (EqualsIgnoreCase(idx_side, "r")) p.swap_sides = true;
  if (EqualsIgnoreCase(idx_side, "s")) p.swap_sides = false;

  // Candidate-cell impl
  // Accepted keys:
  //   - rej_cand_impl = "slab" / "naive" / "A" / "B"
  {
    const std::string impl = GetStringOr(extra, "rej_cand_impl", "slab");
    if (EqualsIgnoreCase(impl, "slab") || EqualsIgnoreCase(impl, "b")) {
      p.cand_impl = CandidateCellsImpl::Slab;
    } else if (EqualsIgnoreCase(impl, "naive") || EqualsIgnoreCase(impl, "a")) {
      p.cand_impl = CandidateCellsImpl::Naive;
    } else {
      if (err) *err = "ReadRejectionParams: unknown rej_cand_impl (expected slab/naive/A/B)";
      return false;
    }
  }

  // Pilot draws
  p.count_draws = GetU64Or(extra, "rej_count_draws", p.count_draws);

  // Grid size defaults: derive from dataset domain.
  const auto dom = ds.Domain();
  for (int d = 0; d < Dim; ++d) {
    const long double w = static_cast<long double>(dom.hi.v[d]) - static_cast<long double>(dom.lo.v[d]);
    // Heuristic default: 64 bins per dimension.
    double gd = 1.0;
    if (w > 0.0L) gd = static_cast<double>(w / 64.0L);
    if (!(gd > 0.0) || !std::isfinite(gd)) gd = 1.0;
    p.g[d] = gd;
  }

  // Allow overriding all dims by rej_g.
  if (auto v = GetExtra(extra, "rej_g")) {
    double g_all = 0;
    if (ParseDouble(*v, &g_all) && (g_all > 0.0) && std::isfinite(g_all)) {
      for (int d = 0; d < Dim; ++d) p.g[d] = g_all;
    }
  }
  // Per-dim overrides: rej_g0, rej_g1, ...
  for (int d = 0; d < Dim; ++d) {
    const std::string key = "rej_g" + std::to_string(d);
    double gd = GetDoubleOr(extra, key, p.g[d]);
    if (gd > 0.0 && std::isfinite(gd)) p.g[d] = gd;
  }

  // Final clamp: ensure > 0.
  for (int d = 0; d < Dim; ++d) {
    if (!(p.g[d] > 0.0) || !std::isfinite(p.g[d])) p.g[d] = 1.0;
  }

  *out = p;
  return true;
}

// ---------------------------------
// RejectionState: preprocessing + sampling + enumerator support
// ---------------------------------

template <int Dim, class T>
struct RejectionState {
  using DatasetT = Dataset<Dim, T>;
  using RelT = Relation<Dim, T>;
  using BoxT = Box<Dim, T>;

  using KeyT = CellKey<Dim>;
  using CellMapT = std::unordered_map<KeyT, std::vector<u32>, CellKeyHash<Dim>>;

  // Slab index entry for doc §3.8 Impl B
  struct SlabEntry {
    KeyT key{};
    const std::vector<u32>* cell = nullptr;  // points into cell_map
  };
  using SlabIndexT = std::unordered_map<i64, std::vector<SlabEntry>>;

  struct ActiveA {
    u32 a_index = 0;  // index into A->boxes
    u64 mu = 0;       // μ(a)
    std::vector<const std::vector<u32>*> cells;  // pointers into cell_map (candidate cells C(a))
    sampling::AliasTable alias_cells;            // weights = cell sizes
  };

  // Built state
  bool built = false;
  const DatasetT* ds = nullptr;
  const RelT* A = nullptr;  // query side
  const RelT* B = nullptr;  // indexed side
  bool swapped = false;     // swapped == true means A=ds.S, B=ds.R

  RejectionParams<Dim> params{};

  // Suppress sign-conversion warning: Dim is compile-time constant and always positive
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsign-conversion"
  std::array<T, Dim> M{};  // per-dim max length on B
  #pragma GCC diagnostic pop
  CellMapT cell_map;

  // Optional acceleration for candidate-cell enumeration (doc §3.8 Impl B).
  SlabIndexT slab_index;

  std::vector<ActiveA> active_a;   // A' = {a | μ(a)>0}
  sampling::AliasTable alias_a;    // AliasR over active_a with weight μ(a)
  u64 mu_sum = 0;                  // MuSum = Σ_{a in A} μ(a) (actually over active_a)

  // Diagnostics from the last pilot count.
  u64 last_pilot_draws = 0;
  u64 last_pilot_accepts = 0;

  void Reset() {
    built = false;
    ds = nullptr;
    A = nullptr;
    B = nullptr;
    swapped = false;
    params = RejectionParams<Dim>{};
    M = {};
    cell_map.clear();
    slab_index.clear();
    active_a.clear();
    alias_a.Clear();
    mu_sum = 0;
    last_pilot_draws = 0;
    last_pilot_accepts = 0;
  }

  // Map (a_index,b_index) (in A/B coordinates) to output PairId (R,S order).
  PairId MakeOutputPair(u32 a_index, u32 b_index) const {
    SJS_DASSERT(ds != nullptr);
    if (!swapped) {
      // A=R, B=S
      return PairId{ds->R.GetId(static_cast<usize>(a_index)), ds->S.GetId(static_cast<usize>(b_index))};
    }
    // A=S, B=R
    return PairId{ds->R.GetId(static_cast<usize>(b_index)), ds->S.GetId(static_cast<usize>(a_index))};
  }

  const BoxT& ABox(u32 a_index) const { return A->boxes[static_cast<usize>(a_index)]; }
  const BoxT& BBox(u32 b_index) const { return B->boxes[static_cast<usize>(b_index)]; }

  // Compute the cell key for an indexed box b in B using rep(b)=lo(b).
  KeyT KeyOfB(const BoxT& b) const {
    KeyT key;
    for (int d = 0; d < Dim; ++d) {
      const long double x = static_cast<long double>(b.lo.v[d]);               // rep(b) := lo(b)
      const long double g = static_cast<long double>(params.g[d]);             // g_i > 0
      key.k[d] = FloorDiv(x, g);
    }
    return key;
  }

  // Candidate cell index range for E(a)=Π_i [lo_i(a)-M_i, hi_i(a)) (half-open).
  // Suppress sign-conversion warning: Dim is compile-time constant and always positive
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsign-conversion"
  void CandidateRange(const BoxT& a,
                      std::array<i64, Dim>* kmin,
                      std::array<i64, Dim>* kmax) const {
    SJS_DASSERT(kmin && kmax);
    for (int d = 0; d < Dim; ++d) {
      const long double g = static_cast<long double>(params.g[d]);
      const long double lo = static_cast<long double>(a.lo.v[d]) - static_cast<long double>(M[d]);
      const long double hi = static_cast<long double>(a.hi.v[d]);
      (*kmin)[d] = FloorDiv(lo, g);
      (*kmax)[d] = CeilDivMinus1(hi, g);
    }
  }
  #pragma GCC diagnostic pop

  // Enumerate all non-empty candidate cells for a (deterministic order).
  // For each found cell, calls emit(ptr_to_vector).
  template <class EmitFn>
  void ForEachCandidateCell(const BoxT& a, EmitFn&& emit) const {
    // Suppress sign-conversion warning: Dim is compile-time constant and always positive
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wsign-conversion"
    std::array<i64, Dim> mn{};
    std::array<i64, Dim> mx{};
    #pragma GCC diagnostic pop
    CandidateRange(a, &mn, &mx);

    // Empty range quick check.
    for (int d = 0; d < Dim; ++d) {
      if (mn[d] > mx[d]) return;
    }

    // Impl B (Slab): iterate only non-empty keys by k1 buckets.
    if (params.cand_impl == CandidateCellsImpl::Slab && !slab_index.empty()) {
      for (i64 k1 = mn[0]; k1 <= mx[0]; ++k1) {
        auto it = slab_index.find(k1);
        if (it == slab_index.end()) continue;
        const auto& bucket = it->second;
        for (const SlabEntry& e : bucket) {
          bool ok = true;
          for (int d = 1; d < Dim; ++d) {
            const i64 kd = e.key.k[d];
            if (kd < mn[d] || kd > mx[d]) {
              ok = false;
              break;
            }
          }
          if (!ok) continue;
          if (e.cell && !e.cell->empty()) emit(e.cell);
        }
      }
      return;
    }

    // Impl A (Naive): enumerate key hyper-rectangle and lookup cell_map.
    KeyT key;
    std::array<i64, Dim> cur{};

    const auto rec = [&](auto&& self, int dim) -> void {
      if (dim == Dim) {
        key.k = cur;
        auto it = cell_map.find(key);
        if (it != cell_map.end() && !it->second.empty()) emit(&it->second);
        return;
      }
      for (i64 v = mn[dim]; v <= mx[dim]; ++v) {
        cur[dim] = v;
        self(self, dim + 1);
      }
    };

    rec(rec, 0);
  }

  // Build preprocessing structures: M_i, CellMap, (optional) SlabIndex, μ/alias tables.
  bool BuildIndex(const DatasetT& ds_in, const Config& cfg, PhaseRecorder* phases, std::string* err) {
    Reset();

    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    ds = &ds_in;

    if (!ReadRejectionParams<Dim, T>(ds_in, cfg, &params, err)) return false;

    // Choose A/B.
    swapped = params.swap_sides;
    if (!swapped) {
      A = &ds->R;
      B = &ds->S;
    } else {
      A = &ds->S;
      B = &ds->R;
    }

    // Capacity checks: we store indices as u32.
    if (A->Size() > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        B->Size() > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RejectionState::BuildIndex: relation size exceeds u32 capacity";
      return false;
    }

    // (1) Compute M_i on B.
    {
      auto sc = phases ? phases->Scoped("build_M") : PhaseRecorder::ScopedPhase(nullptr, "");
      for (int d = 0; d < Dim; ++d) M[d] = static_cast<T>(0);
      for (usize i = 0; i < B->boxes.size(); ++i) {
        const BoxT& b = B->boxes[i];
        if (b.IsEmpty()) continue;
        for (int d = 0; d < Dim; ++d) {
          const T len = static_cast<T>(b.hi.v[d] - b.lo.v[d]);
          if (len > M[d]) M[d] = len;
        }
      }
    }

    // (2) Build CellMap on B using rep(b)=lo(b).
    {
      auto sc = phases ? phases->Scoped("build_cellmap") : PhaseRecorder::ScopedPhase(nullptr, "");
      cell_map.clear();
      cell_map.reserve(std::max<usize>(B->Size(), 1) * 2);

      for (usize i = 0; i < B->boxes.size(); ++i) {
        const BoxT& b = B->boxes[i];
        if (b.IsEmpty()) continue;
        const KeyT key = KeyOfB(b);
        cell_map[key].push_back(static_cast<u32>(i));
      }
    }

    // (3) Optional SlabIndex (doc §3.8 Impl B).
    slab_index.clear();
    if (params.cand_impl == CandidateCellsImpl::Slab) {
      auto sc = phases ? phases->Scoped("build_slab_index") : PhaseRecorder::ScopedPhase(nullptr, "");
      slab_index.reserve(std::max<usize>(cell_map.size(), 1));
      for (const auto& kv : cell_map) {
        const KeyT& key = kv.first;
        const auto& vec = kv.second;
        if (vec.empty()) continue;
        slab_index[key.k[0]].push_back(SlabEntry{key, &vec});
      }
      // Sort each bucket for deterministic enumeration order.
      for (auto& kv : slab_index) {
        auto& bucket = kv.second;
        std::sort(bucket.begin(), bucket.end(), [](const SlabEntry& a, const SlabEntry& b) {
          return KeyLexLess<Dim>(a.key, b.key);
        });
      }
    }

    // (4) For each a in A: build candidate cells C(a), μ(a), AliasCell[a].
    {
      auto sc = phases ? phases->Scoped("build_mu_alias") : PhaseRecorder::ScopedPhase(nullptr, "");
      active_a.clear();
      std::vector<u64> w_a;
      w_a.reserve(A->Size());

      mu_sum = 0;

      for (u32 ai = 0; ai < static_cast<u32>(A->Size()); ++ai) {
        const BoxT& a = ABox(ai);
        if (a.IsEmpty()) continue;

        ActiveA rec;
        rec.a_index = ai;
        rec.mu = 0;

        std::vector<u64> cell_w;

        ForEachCandidateCell(a, [&](const std::vector<u32>* cell) {
          if (!cell || cell->empty()) return;
          const u64 w = static_cast<u64>(cell->size());
          rec.cells.push_back(cell);
          cell_w.push_back(w);
          rec.mu += w;
        });

        if (rec.mu == 0) continue;

        // Build alias on C(a) with weights |CellMap[k]|.
        std::string aerr;
        if (!rec.alias_cells.BuildFromU64(Span<const u64>(cell_w), &aerr)) {
          if (err) *err = "RejectionState::BuildIndex: AliasCell build failed: " + aerr;
          return false;
        }

        mu_sum += rec.mu;
        w_a.push_back(rec.mu);
        active_a.push_back(std::move(rec));
      }

      if (active_a.empty()) {
        // By doc §4.0 necessary condition: MuSum==0 => join empty.
        alias_a.Clear();
        built = true;
        return true;
      }

      // (5) Global AliasR over A' with weights μ(a).
      std::string aerr;
      if (!alias_a.BuildFromU64(Span<const u64>(w_a), &aerr)) {
        if (err) *err = "RejectionState::BuildIndex: AliasR build failed: " + aerr;
        return false;
      }
    }

    built = true;
    return true;
  }

  // Draw one proposal (a_index, b_index) from U_prop.
  // Returns false if MuSum==0 (no proposals possible).
  bool DrawProposal(Rng* rng, u32* out_a_index, u32* out_b_index) const {
    SJS_DASSERT(rng != nullptr);
    SJS_DASSERT(out_a_index != nullptr);
    SJS_DASSERT(out_b_index != nullptr);
    if (mu_sum == 0 || active_a.empty()) return false;

    const usize ai_pos = alias_a.Sample(rng);
    SJS_DASSERT(ai_pos < active_a.size());
    const ActiveA& ar = active_a[ai_pos];

    const usize ci = ar.alias_cells.Sample(rng);
    SJS_DASSERT(ci < ar.cells.size());
    const std::vector<u32>* cell = ar.cells[ci];
    SJS_DASSERT(cell != nullptr);
    const u32 n = static_cast<u32>(cell->size());
    SJS_DASSERT(n > 0);

    const u32 bj = (*cell)[static_cast<usize>(rng->UniformU32(n))];

    *out_a_index = ar.a_index;
    *out_b_index = bj;
    return true;
  }

  // Pilot estimate of |J| (doc §5.4): |J| = MuSum * p_acc, p_acc ≈ accepts/draws.
  // IMPORTANT: caller should pass an RNG copy to keep the main RNG stream unchanged.
  bool EstimateCountByPilot(Rng* rng,
                            CountResult* out,
                            PhaseRecorder* phases,
                            std::string* err) {
    if (!built || !ds || !A || !B) {
      if (err) *err = "RejectionState::EstimateCountByPilot: call BuildIndex() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "RejectionState::EstimateCountByPilot: null rng/out";
      return false;
    }

    auto scoped = phases ? phases->Scoped("count_pilot") : PhaseRecorder::ScopedPhase(nullptr, "");

    if (mu_sum == 0) {
      *out = MakeExactCount(0);
      last_pilot_draws = 0;
      last_pilot_accepts = 0;
      return true;
    }

    const u64 m = params.count_draws;
    if (m == 0) {
      // Pilot disabled: return "unknown" estimate.
      *out = MakeEstimateCount(std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              0);
      last_pilot_draws = 0;
      last_pilot_accepts = 0;
      return true;
    }

    u64 acc = 0;
    u64 draws = 0;
    for (u64 i = 0; i < m; ++i) {
      u32 ai = 0, bi = 0;
      if (!DrawProposal(rng, &ai, &bi)) break;
      ++draws;
      if (ABox(ai).Intersects(BBox(bi))) ++acc;
    }

    last_pilot_draws = draws;
    last_pilot_accepts = acc;

    if (draws == 0) {
      *out = MakeEstimateCount(0.0L, 0.0L, 0.0L, 0.0L, 0);
      return true;
    }

    const long double p_hat = static_cast<long double>(acc) / static_cast<long double>(draws);
    const long double mu = static_cast<long double>(mu_sum);
    const long double est = mu * p_hat;

    // Binomial stderr for p_hat.
    const long double var = (p_hat * (1.0L - p_hat)) / static_cast<long double>(draws);
    const long double se_p = (var > 0.0L) ? std::sqrt(var) : 0.0L;
    const long double se_est = mu * se_p;

    // Normal approximation CI.
    const long double z = 1.9599639845400542L;  // ~N(0,1) 97.5% quantile
    long double lo = p_hat - z * se_p;
    long double hi = p_hat + z * se_p;
    if (lo < 0.0L) lo = 0.0L;
    if (hi > 1.0L) hi = 1.0L;

    *out = MakeEstimateCount(est, se_est, mu * lo, mu * hi, draws);
    return true;
  }
};

// ---------------------------------
// Deterministic join enumerator (grid candidate scan)
// ---------------------------------

template <int Dim, class T>
class RejectionJoinEnumerator final : public IJoinEnumerator {
 public:
  using StateT = RejectionState<Dim, T>;
  using BoxT = Box<Dim, T>;

  explicit RejectionJoinEnumerator(const StateT* st) : st_(st) { Reset(); }

  void Reset() override {
    stats_ = join::JoinStats{};
    a_pos_ = 0;
    cell_i_ = 0;
    b_pos_ = 0;
    cells_.clear();
    cur_cell_ = nullptr;

    if (!st_ || !st_->built || !st_->A || !st_->B) return;
    PrepareNextA();
  }

  bool Next(PairId* out) override {
    if (!out) return false;
    if (!st_ || !st_->built || !st_->A || !st_->B) return false;

    while (true) {
      if (a_pos_ >= st_->A->Size()) return false;  // exhausted A

      if (cur_cell_ == nullptr) {
        if (!AdvanceCell()) {
          ++a_pos_;
          PrepareNextA();
          continue;
        }
      }

      SJS_DASSERT(cur_cell_ != nullptr);
      while (b_pos_ < cur_cell_->size()) {
        const u32 bi = (*cur_cell_)[b_pos_++];
        stats_.candidate_checks++;

        const BoxT& a = st_->A->boxes[a_pos_];
        const BoxT& b = st_->B->boxes[static_cast<usize>(bi)];
        if (!a.Intersects(b)) continue;

        *out = st_->MakeOutputPair(static_cast<u32>(a_pos_), bi);
        stats_.output_pairs++;
        return true;
      }

      cur_cell_ = nullptr;  // done with this cell
    }
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  void PrepareNextA() {
    cells_.clear();
    cur_cell_ = nullptr;
    cell_i_ = 0;
    b_pos_ = 0;

    while (st_ && st_->A && a_pos_ < st_->A->Size()) {
      const BoxT& a = st_->A->boxes[a_pos_];
      if (!a.IsEmpty()) {
        st_->ForEachCandidateCell(a, [&](const std::vector<u32>* cell) {
          if (cell && !cell->empty()) cells_.push_back(cell);
        });
      }
      if (!cells_.empty()) break;
      ++a_pos_;
    }
  }

  bool AdvanceCell() {
    while (cell_i_ < cells_.size()) {
      cur_cell_ = cells_[cell_i_++];
      b_pos_ = 0;
      if (cur_cell_ && !cur_cell_->empty()) return true;
    }
    cur_cell_ = nullptr;
    return false;
  }

  const StateT* st_{nullptr};
  join::JoinStats stats_;

  usize a_pos_{0};
  std::vector<const std::vector<u32>*> cells_;
  usize cell_i_{0};
  const std::vector<u32>* cur_cell_{nullptr};
  usize b_pos_{0};
};

}  // namespace detail

// ---------------------------------
// RejectionSamplingBaseline (AGR-S)
// ---------------------------------

template <int Dim, class T = Scalar>
class RejectionSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "RejectionSamplingBaseline expects Dim>=2");
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Rejection; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "rejection_sampling"; }

  void Reset() override { st_.Reset(); }

  bool Build(const DatasetT& ds,
             const Config& cfg,
             PhaseRecorder* phases,
             std::string* err) override {
    return st_.BuildIndex(ds, cfg, phases, err);
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;  // cfg already used in Build; pilot uses stored params (count_draws)
    if (!st_.built) {
      if (err) *err = "RejectionSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RejectionSamplingBaseline::Count: out is null";
      return false;
    }

    // Doc §5.4: pilot must be RNG-isolated so Sample() remains reproducible.
    if (!rng) {
      *out = MakeEstimateCount(std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              0);
      return true;
    }
    Rng tmp = *rng;
    return st_.EstimateCountByPilot(&tmp, out, phases, err);
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!st_.built || !st_.ds || !st_.A || !st_.B) {
      if (err) *err = "RejectionSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "RejectionSamplingBaseline::Sample: null rng/out";
      return false;
    }

    auto scoped = phases ? phases->Scoped("sample") : PhaseRecorder::ScopedPhase(nullptr, "");

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    // Doc §4.0: MuSum==0 is a necessary condition for empty join.
    if (st_.mu_sum == 0) return true;

    out->pairs.reserve(static_cast<usize>(t));

    // AGR-S: repeat proposals until t accepts.
    //
    // Robustness: if MuSum>0 but the true join is empty, naive rejection
    // sampling would not terminate. We therefore *prove emptiness* (or confirm
    // non-emptiness) once if we observe too many rejected proposals before the
    // first accept.
    const u64 empty_check_after = detail::GetU64Or(cfg.run.extra, "rej_max_factor", 10000);

    auto prove_nonempty_or_empty = [&]() -> bool {
      // Returns true iff it finds at least one true intersecting pair.
      // Deterministic scan over the candidate superset U_prop.
      for (const auto& ar : st_.active_a) {
        const auto& a = st_.ABox(ar.a_index);
        for (const std::vector<u32>* cell : ar.cells) {
          if (!cell || cell->empty()) continue;
          for (u32 bi : *cell) {
            if (a.Intersects(st_.BBox(bi))) return true;
          }
        }
      }
      return false;
    };

    u64 rejected_before_first_accept = 0;
    bool checked_empty = false;

    while (out->pairs.size() < static_cast<usize>(t)) {
      u32 ai = 0, bi = 0;
      if (!st_.DrawProposal(rng, &ai, &bi)) return true;  // should only happen if MuSum==0

      if (!st_.ABox(ai).Intersects(st_.BBox(bi))) {
        // Still searching for the first accept.
        if (out->pairs.empty()) {
          ++rejected_before_first_accept;
          if (!checked_empty && (empty_check_after == 0 || rejected_before_first_accept >= empty_check_after)) {
            checked_empty = true;
            auto sc = phases ? phases->Scoped("sample_empty_check")
                              : PhaseRecorder::ScopedPhase(nullptr, "");
            if (!prove_nonempty_or_empty()) {
              // Certified empty join: return empty sample set.
              return true;
            }
            // Join is non-empty; continue with standard rejection sampling.
          }
        }
        continue;
      }

      // Accept.
      out->pairs.push_back(st_.MakeOutputPair(ai, bi));
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!st_.built) {
      if (err) *err = "RejectionSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RejectionJoinEnumerator<Dim, T>>(&st_);
  }

 private:
  detail::RejectionState<Dim, T> st_;
};

}  // namespace rejection
}  // namespace baselines
}  // namespace sjs
