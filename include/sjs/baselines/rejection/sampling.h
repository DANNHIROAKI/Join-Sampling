#pragma once
// sjs/baselines/rejection/sampling.h
//
// AGR-BoxJoin baseline (Amagata'25 style): grid upper bound \mu + alias + rejection.
// Variant::Sampling (AGR-S).
//
// High-level idea
// ---------------
// Build():
//   1) Pick a query side A and an indexed side B (default: A=R, B=S).
//   2) Compute per-dimension max side length on B: M_i = max_b (hi_i - lo_i).
//   3) Build a sparse hash grid on rep(b) (default rep(b)=b.lo).
//   4) For each a in A, compute candidate cells C(a) that intersect
//      E(a) = \prod_i [lo_i(a) - M_i, hi_i(a)).
//      Define \mu(a) = \sum_{cell in C(a)} |CellMap[cell]|.
//      Build AliasCell[a] on those cell sizes.
//   5) Build global AliasA over all a with \mu(a) > 0.
//
// Sample():
//   Repeat until t accepted samples:
//     a ~ AliasA (weight \mu(a))
//     cell ~ AliasCell[a] (weight cell size)
//     b ~ Uniform(CellMap[cell])
//     if Intersects(a,b): accept, output the join pair; else reject.
//
// Correctness sketch:
//   For any (a,b) such that b lies in a candidate cell of a, the proposal
//   probability is 1 / \sum_a \mu(a) (a constant). Conditioning on accept
//   (which is exactly the true intersection predicate) yields uniform sampling
//   over the true join J. Independent RNG draws per proposal imply i.i.d.
//
// Notes
// -----
// * Geometry uses half-open boxes [lo,hi) as in the project.
// * This is a baseline: performance depends heavily on how tight \mu is.
// * Empty-join corner case: if |J|=0 but \sum \mu > 0, pure rejection would
//   not terminate. We therefore include a configurable maximum-proposals guard
//   (and the Adaptive variant should be used for robust empty-join handling).

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

// --------------------------
// Extra-map parsing helpers
// --------------------------

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

inline bool ParseI64(std::string_view s, i64* out) {
  if (!out || s.empty()) return false;
  try {
    std::size_t idx = 0;
    long long v = std::stoll(std::string(s), &idx, 10);
    if (idx != s.size()) return false;
    *out = static_cast<i64>(v);
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

// --------------------------
// Cell-key & hashing (sparse hash grid)
// --------------------------

template <int Dim>
struct CellKey {
  std::array<i64, Dim> k{};

  bool operator==(const CellKey& o) const noexcept { return k == o.k; }
};

template <int Dim>
struct CellKeyHash {
  // splitmix64 finalizer
  static inline u64 Mix64(u64 x) noexcept {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x = x ^ (x >> 31);
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

// --------------------------
// Numeric helpers for grid indexing
// --------------------------

inline i64 ClampToI64(long double v) noexcept {
  const long double lo = static_cast<long double>(std::numeric_limits<i64>::min());
  const long double hi = static_cast<long double>(std::numeric_limits<i64>::max());
  if (v < lo) return std::numeric_limits<i64>::min();
  if (v > hi) return std::numeric_limits<i64>::max();
  return static_cast<i64>(v);
}

inline i64 FloorDiv(long double x, long double g) noexcept {
  // g must be > 0.
  return ClampToI64(std::floor(x / g));
}

inline i64 CeilDivMinus1(long double x, long double g) noexcept {
  // For half-open [*, x): k_max = ceil(x/g) - 1.
  return ClampToI64(std::ceil(x / g) - 1.0L);
}

// Representative point for an indexed box (default: lo corner).
// If center=true, uses the box center.
// NOTE: The correctness argument (coverage via E(a)) in Amagata'25-style proof
// assumes rep(b)=lo(b). Center is provided for future experimentation and
// may require adapting the upper-bound construction.
template <int Dim, class T>
inline Point<Dim, T> RepPoint(const Box<Dim, T>& b, bool center) noexcept {
  if (!center) return b.lo;
  Point<Dim, T> p;
  for (int d = 0; d < Dim; ++d) {
    p.v[d] = static_cast<T>((static_cast<long double>(b.lo.v[d]) + static_cast<long double>(b.hi.v[d])) * 0.5L);
  }
  return p;
}

// --------------------------
// Parameter bundle
// --------------------------

template <int Dim>
struct RejectionParams {
  // Grid cell size per dimension. Must be > 0.
  std::array<double, Dim> g{};

  // If true: index on R and query S (swap sides). Default false (index on S).
  bool swap_sides = false;

  // Representative point choice.
  bool rep_center = false;

  // Count() pilot draws for estimating |J| (Sampling variant).
  u64 count_draws = 50'000ULL;

  // Maximum proposals allowed in Sample() to avoid infinite loops.
  // If 0, we use max_factor * t.
  u64 max_proposals = 0ULL;

  // Default multiplier for max proposals when max_proposals==0.
  u64 max_factor = 10'000ULL;
};

template <int Dim, class T>
inline RejectionParams<Dim> ReadRejectionParams(const Dataset<Dim, T>& ds, const Config& cfg) {
  RejectionParams<Dim> p;
  const auto& extra = cfg.run.extra;

  // Side selection
  // Accepted keys:
  //   - rej_swap = true/false
  //   - rej_index_side = "R" or "S" (which relation is indexed B)
  p.swap_sides = GetBoolOr(extra, "rej_swap", false);
  const std::string idx_side = GetStringOr(extra, "rej_index_side", GetStringOr(extra, "rej_index", "S"));
  if (sjs::detail::EqualsIgnoreCase(idx_side, "r")) p.swap_sides = true;
  if (sjs::detail::EqualsIgnoreCase(idx_side, "s")) p.swap_sides = false;

  // Representative point
  p.rep_center = GetBoolOr(extra, "rej_rep_center", false);

  // Count pilot
  p.count_draws = GetU64Or(extra, "rej_count_draws", p.count_draws);

  // Sample guard
  p.max_proposals = GetU64Or(extra, "rej_max_proposals", p.max_proposals);
  p.max_factor = GetU64Or(extra, "rej_max_factor", p.max_factor);
  if (p.max_factor == 0) p.max_factor = 1;

  // Grid size defaults: derive from dataset domain.
  const auto dom = ds.Domain();
  for (int d = 0; d < Dim; ++d) {
    const long double w = static_cast<long double>(dom.hi.v[d]) - static_cast<long double>(dom.lo.v[d]);
    // Heuristic default: 64 bins per dimension.
    double gd = 1.0;
    if (w > 0.0L) gd = static_cast<double>(w / 64.0L);
    if (!(gd > 0.0)) gd = 1.0;
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

  return p;
}

// --------------------------
// RejectionState: shared preprocessing + sampling + enumerator support
// --------------------------

template <int Dim, class T>
struct RejectionState {
  using DatasetT = Dataset<Dim, T>;
  using RelT = Relation<Dim, T>;
  using BoxT = Box<Dim, T>;

  // Indexed side B is built into CellMap.
  using KeyT = CellKey<Dim>;
  using CellMapT = std::unordered_map<KeyT, std::vector<u32>, CellKeyHash<Dim>>;

  struct ActiveA {
    u32 a_index = 0;  // index into A->boxes
    u64 mu = 0;       // \mu(a)
    std::vector<const std::vector<u32>*> cells;  // candidate cell lists (pointers into cell_map)
    sampling::AliasTable alias_cells;            // weights = cell sizes
  };

  // Built state
  bool built = false;
  const DatasetT* ds = nullptr;
  const RelT* A = nullptr;  // query side
  const RelT* B = nullptr;  // indexed side
  bool swapped = false;     // swapped == true means A=ds.S, B=ds.R

  RejectionParams<Dim> params{};

  std::array<T, Dim> M{};  // per-dim max length on B
  CellMapT cell_map;

  std::vector<ActiveA> active_a;
  sampling::AliasTable alias_a;
  u64 mu_sum = 0;  // sum over active A of \mu(a)

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

  // Compute the cell key for an indexed box b in B.
  KeyT KeyOfB(const BoxT& b) const {
    KeyT key;
    const auto p = RepPoint<Dim, T>(b, params.rep_center);
    for (int d = 0; d < Dim; ++d) {
      const long double x = static_cast<long double>(p.v[d]);
      const long double g = static_cast<long double>(params.g[d]);
      key.k[d] = FloorDiv(x, g);
    }
    return key;
  }

  // Candidate cell range for an a-box (half-open interval logic).
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

  // Enumerate all non-empty candidate cells for a (in deterministic lexicographic order).
  // For each found cell, calls emit(ptr_to_vector).
  template <class EmitFn>
  void ForEachCandidateCell(const BoxT& a, EmitFn&& emit) const {
    std::array<i64, Dim> mn{};
    std::array<i64, Dim> mx{};
    CandidateRange(a, &mn, &mx);

    // Empty range quick check.
    for (int d = 0; d < Dim; ++d) {
      if (mn[d] > mx[d]) return;
    }

    KeyT key;
    std::array<i64, Dim> cur{};

    // Recursive enumeration over dimensions.
    const auto rec = [&](auto&& self, int dim) -> void {
      if (dim == Dim) {
        key.k = cur;
        auto it = cell_map.find(key);
        if (it != cell_map.end() && !it->second.empty()) {
          emit(&it->second);
        }
        return;
      }
      for (i64 v = mn[dim]; v <= mx[dim]; ++v) {
        cur[dim] = v;
        self(self, dim + 1);
      }
    };

    rec(rec, 0);
  }

  // Build preprocessing structures.
  bool BuildIndex(const DatasetT& ds_in, const Config& cfg, PhaseRecorder* phases, std::string* err) {
    Reset();

    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds = &ds_in;
    params = ReadRejectionParams<Dim, T>(ds_in, cfg);

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

    // Compute M_i on B.
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

    // Build CellMap on B.
    {
      auto sc = phases ? phases->Scoped("build_cellmap") : PhaseRecorder::ScopedPhase(nullptr, "");
      cell_map.clear();
      cell_map.reserve(std::max<usize>(B->Size(), 1) * 2);

      for (usize i = 0; i < B->boxes.size(); ++i) {
        const BoxT& b = B->boxes[i];
        if (b.IsEmpty()) continue;
        const KeyT key = KeyOfB(b);
        auto& vec = cell_map[key];
        vec.push_back(static_cast<u32>(i));
      }
    }

    // Build per-a candidate cells + AliasCell[a] for a with \mu(a) > 0.
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

        // Collect candidate cells and weights.
        std::vector<u64> cell_w;

        ForEachCandidateCell(a, [&](const std::vector<u32>* cell) {
          SJS_DASSERT(cell != nullptr);
          const u64 w = static_cast<u64>(cell->size());
          if (w == 0) return;
          rec.cells.push_back(cell);
          cell_w.push_back(w);
          rec.mu += w;
        });

        if (rec.mu == 0) continue;

        // Build alias for this a.
        std::string aerr;
        if (!rec.alias_cells.BuildFromU64(Span<const u64>(cell_w), &aerr)) {
          if (err) *err = "RejectionState::BuildIndex: AliasCell build failed: " + aerr;
          return false;
        }

        mu_sum += rec.mu;
        w_a.push_back(rec.mu);
        active_a.push_back(std::move(rec));
      }

      // If mu_sum==0 => join must be empty.
      if (active_a.empty()) {
        alias_a.Clear();
        built = true;
        return true;
      }

      // Global alias over active A.
      std::string aerr;
      if (!alias_a.BuildFromU64(Span<const u64>(w_a), &aerr)) {
        if (err) *err = "RejectionState::BuildIndex: AliasA build failed: " + aerr;
        return false;
      }
    }

    built = true;
    return true;
  }

  // Draw one proposal (a_index, b_index) from the candidate superset.
  // Returns false if mu_sum==0 (no proposals possible).
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
    if (n == 0) return false;  // should not happen

    const u32 bj = (*cell)[static_cast<usize>(rng->UniformU32(n))];

    *out_a_index = ar.a_index;
    *out_b_index = bj;
    return true;
  }

  // Estimate |J| using a fixed number of proposals (pilot) and intersection tests.
  bool EstimateCountByPilot(const Config& cfg,
                            Rng* rng,
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

    const u64 m = ReadRejectionParams<Dim, T>(*ds, cfg).count_draws;
    if (m == 0) {
      // No pilot requested.
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

    // Normal approximation CI (good when draws is not too small).
    const long double z = 1.9599639845400542L;  // ~N(0,1) 97.5% quantile
    long double lo = p_hat - z * se_p;
    long double hi = p_hat + z * se_p;
    if (lo < 0.0L) lo = 0.0L;
    if (hi > 1.0L) hi = 1.0L;

    *out = MakeEstimateCount(est, se_est, mu * lo, mu * hi, draws);
    return true;
  }
};

// --------------------------
// Deterministic join enumerator (grid candidate scan)
// --------------------------

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
      // Exhausted all A.
      if (a_pos_ >= st_->A->Size()) return false;

      // If no current cell, advance.
      if (cur_cell_ == nullptr) {
        if (!AdvanceCell()) {
          // Move to next A.
          ++a_pos_;
          PrepareNextA();
          continue;
        }
      }

      // Scan within current cell.
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

      // Done with this cell; move to next cell.
      cur_cell_ = nullptr;
    }
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  void PrepareNextA() {
    cells_.clear();
    cur_cell_ = nullptr;
    cell_i_ = 0;
    b_pos_ = 0;

    // Find first A with non-empty candidate cells.
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

  void PrepareNextAStub() { (void)0; }

  const StateT* st_{nullptr};
  join::JoinStats stats_;

  usize a_pos_{0};
  std::vector<const std::vector<u32>*> cells_;
  usize cell_i_{0};
  const std::vector<u32>* cur_cell_{nullptr};
  usize b_pos_{0};
};

}  // namespace detail

// --------------------------
// RejectionSamplingBaseline (AGR-S)
// --------------------------

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
    // Pilot-based estimate by default.
    return st_.EstimateCountByPilot(cfg, rng, out, phases, err);
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

    if (st_.mu_sum == 0) {
      // Necessary condition for non-empty join.
      return true;
    }

    // Proposal guard.
    const detail::RejectionParams<Dim> p = detail::ReadRejectionParams<Dim, T>(*st_.ds, cfg);
    u64 max_props = p.max_proposals;
    if (max_props == 0) {
      // Default: max_factor * t, clamped.
      long double m = static_cast<long double>(p.max_factor) * static_cast<long double>(t);
      const long double cap = static_cast<long double>(std::numeric_limits<u64>::max());
      if (m > cap) m = cap;
      max_props = static_cast<u64>(m);
      // Ensure a small absolute floor so tiny t doesn't set a tiny budget.
      if (max_props < 1000ULL) max_props = 1000ULL;
    }

    out->pairs.reserve(static_cast<usize>(t));

    u64 props = 0;
    while (out->pairs.size() < static_cast<usize>(t)) {
      if (max_props > 0 && props >= max_props) {
        if (err) {
          *err = "RejectionSamplingBaseline::Sample: reached max_proposals without collecting enough accepts; "
                 "acceptance rate may be ~0 (empty join?) or extremely small";
        }
        return false;
      }
      ++props;

      u32 ai = 0, bi = 0;
      if (!st_.DrawProposal(rng, &ai, &bi)) {
        // Should only happen if mu_sum==0.
        return true;
      }

      if (!st_.ABox(ai).Intersects(st_.BBox(bi))) continue;
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
