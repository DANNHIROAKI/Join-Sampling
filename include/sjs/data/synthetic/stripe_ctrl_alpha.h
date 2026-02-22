#pragma once
// sjs/data/synthetic/stripe_ctrl_alpha.h
//
// Stripe-controlled alpha generator (SCIRG-style).
//
// Construction idea (dimension-agnostic):
//  - Pick one "control axis" (default axis=1; requires Dim>=2).
//  - On all other axes, every box covers a shared "core" interval -> guaranteed overlap.
//  - On control axis, S is made of non-overlapping strips separated by gaps.
//  - Each R box spans exactly d_i consecutive strips (or is placed inside a gap if d_i=0).
//  - Then |J| = Σ_i d_i exactly.
//
// IMPORTANT: alpha here follows SJS-HighDims convention used in your docs:
//   alpha := |J| / (|R| + |S|)
// not |J|/(|R||S|).
//
// Parameters (spec.params):
//  - "control_axis" (int, default 1): axis for strip layout.
//  - "core_lo" / "core_hi" (double in [0,1], default 0.45/0.55): core interval fractions for non-control axes.
//  - Gap/strip geometry:
//      Either provide absolute "g" (or "gap") as gap length,
//      OR provide "gap_factor" in (0,1) (default 0.1) and compute:
//          g = gap_factor * L / (n_s + 1)
//  - "delta_factor" (double in (0,0.5), default 0.25):
//      delta = min(g, h) * delta_factor
//      used as safety margin so half-open semantics never accidentally counts touching edges.
//  - "shuffle_strips" (bool, default true): shuffle S order.
//  - "shuffle_r" (bool, default false): shuffle R order.
//  - "swap_sides" (bool, default false): swap which side is "R" vs "S" in output.
//  - k override:
//      "k_target" or "k" (u64). If absent, k = round(alpha * (n_r+n_s)).

#include "sjs/data/synthetic/generator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace synthetic {

namespace detail {

// Sample a random bounded composition (d_0..d_{n_r-1}) such that
//   0 <= d_i <= n_s and Σ d_i = k.
// We sequentially sample d_i uniformly from the feasible range [low,high] given remaining mass.
inline bool RandomBoundedComposition(u64 k,
                                    u64 n_r,
                                    u64 n_s,
                                    Rng* rng,
                                    std::vector<u32>* out_d,
                                    std::string* err) {
  if (!rng || !out_d) {
    SetErr(err, "RandomBoundedComposition: null rng/out");
    return false;
  }
  out_d->assign(static_cast<usize>(n_r), 0u);

  u64 remaining = k;
  for (u64 i = 0; i < n_r; ++i) {
    const u64 left = n_r - i - 1;

#if defined(__SIZEOF_INT128__)
    const __uint128_t max_future = static_cast<__uint128_t>(left) * static_cast<__uint128_t>(n_s);
    u64 low = 0;
    if (static_cast<__uint128_t>(remaining) > max_future) {
      const __uint128_t diff = static_cast<__uint128_t>(remaining) - max_future;
      if (diff > static_cast<__uint128_t>(std::numeric_limits<u64>::max())) {
        SetErr(err, "RandomBoundedComposition: overflow computing low");
        return false;
      }
      low = static_cast<u64>(diff);
    }
#else
    // Fallback: use long double.
    long double max_future = static_cast<long double>(left) * static_cast<long double>(n_s);
    u64 low = 0;
    if (static_cast<long double>(remaining) > max_future) {
      low = static_cast<u64>(static_cast<long double>(remaining) - max_future);
    }
#endif

    const u64 high = (remaining < n_s) ? remaining : n_s;
    if (low > high) {
      SetErr(err, "RandomBoundedComposition: infeasible (low>high).");
      return false;
    }

    const u64 span = high - low + 1;
    const u64 di = low + rng->UniformU64(span);
    (*out_d)[static_cast<usize>(i)] = static_cast<u32>(di);
    remaining -= di;
  }

  if (remaining != 0) {
    SetErr(err, "RandomBoundedComposition: remaining != 0 at end (bug).");
    return false;
  }
  return true;
}

inline u64 RoundToU64(double x) {
  if (!(x >= 0.0) || !std::isfinite(x)) return 0;
  const long double y = static_cast<long double>(x);
  const long double r = std::llround(y);
  if (r <= 0) return 0;
  if (r >= static_cast<long double>(std::numeric_limits<u64>::max())) return std::numeric_limits<u64>::max();
  return static_cast<u64>(r);
}

}  // namespace detail

template <int Dim, class T = Scalar>
class StripeCtrlAlphaGenerator final : public ISyntheticGenerator<Dim, T> {
 public:
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;

  std::string_view Name() const noexcept override { return "stripe_ctrl_alpha"; }

  bool Generate(const DatasetSpec& spec, DatasetT* out_ds, Report* report, std::string* err) override {
    if (!out_ds) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: out_ds is null");
      return false;
    }
    if (spec.n_r == 0 || spec.n_s == 0) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: n_r and n_s must be > 0");
      return false;
    }
    if (!detail::CheckIdFits(spec.n_r, err, "StripeCtrlAlphaGenerator(R)")) return false;
    if (!detail::CheckIdFits(spec.n_s, err, "StripeCtrlAlphaGenerator(S)")) return false;

    if (Dim < 2) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator requires Dim >= 2");
      return false;
    }

    const double dom_lo = spec.domain_lo;
    const double dom_hi = spec.domain_hi;
    if (!(dom_hi > dom_lo)) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: domain_hi must be > domain_lo");
      return false;
    }
    const double L = dom_hi - dom_lo;

    // Parameters.
    const i32 control_axis = detail::GetI32(spec.params, "control_axis", 1);
    if (control_axis < 0 || control_axis >= Dim) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: invalid control_axis");
      return false;
    }

    const double core_lo_frac = detail::GetDouble(spec.params, "core_lo", 0.45);
    const double core_hi_frac = detail::GetDouble(spec.params, "core_hi", 0.55);
    if (!(core_lo_frac >= 0.0 && core_hi_frac <= 1.0 && core_lo_frac < core_hi_frac)) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: need 0<=core_lo<core_hi<=1");
      return false;
    }
    const double core_lo = dom_lo + core_lo_frac * L;
    const double core_hi = dom_lo + core_hi_frac * L;

    // Gap size g: either absolute g/gap, or computed by gap_factor.
    double g = 0.0;
    if (detail::FindParam(spec.params, "g")) {
      g = detail::GetDouble(spec.params, "g", 0.0);
    } else if (detail::FindParam(spec.params, "gap")) {
      g = detail::GetDouble(spec.params, "gap", 0.0);
    } else {
      const double gap_factor = detail::GetDouble(spec.params, "gap_factor", 0.1);
      if (!(gap_factor > 0.0 && gap_factor < 1.0)) {
        detail::SetErr(err, "StripeCtrlAlphaGenerator: gap_factor must be in (0,1)");
        return false;
      }
      g = (gap_factor * L) / static_cast<double>(spec.n_s + 1);
    }

    if (!(g > 0.0)) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: gap g must be > 0");
      return false;
    }
    if (!(static_cast<double>(spec.n_s + 1) * g < L)) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: (n_s+1)*g must be < domain length");
      return false;
    }
    const double h = (L - static_cast<double>(spec.n_s + 1) * g) / static_cast<double>(spec.n_s);
    if (!(h > 0.0)) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: strip height h must be > 0 (check g)");
      return false;
    }

    const double delta_factor = detail::GetDouble(spec.params, "delta_factor", 0.25);
    if (!(delta_factor > 0.0 && delta_factor < 0.5)) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: delta_factor must be in (0,0.5)");
      return false;
    }
    const double delta = std::min(g, h) * delta_factor;

    if (!(delta > 0.0 && delta < 0.5 * g && delta < 0.5 * h)) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: need 0<delta<min(g/2,h/2)");
      return false;
    }

    const bool shuffle_strips = detail::GetBool(spec.params, "shuffle_strips", true);
    const bool shuffle_r = detail::GetBool(spec.params, "shuffle_r", false);
    const bool swap_sides = detail::GetBool(spec.params, "swap_sides", false);

    // Determine target k.
    u64 k_target = 0;
    bool k_overridden = false;
    if (auto v = detail::FindParam(spec.params, "k_target")) {
      if (!detail::TryParseU64(*v, &k_target)) {
        detail::SetErr(err, "StripeCtrlAlphaGenerator: failed to parse k_target");
        return false;
      }
      k_overridden = true;
    } else if (auto v2 = detail::FindParam(spec.params, "k")) {
      if (!detail::TryParseU64(*v2, &k_target)) {
        detail::SetErr(err, "StripeCtrlAlphaGenerator: failed to parse k");
        return false;
      }
      k_overridden = true;
    } else {
      k_target = detail::RoundToU64(spec.alpha * static_cast<double>(spec.n_r + spec.n_s));
    }

    // Feasibility: k <= n_r*n_s
#if defined(__SIZEOF_INT128__)
    const __uint128_t max_k = static_cast<__uint128_t>(spec.n_r) * static_cast<__uint128_t>(spec.n_s);
    if (static_cast<__uint128_t>(k_target) > max_k) {
      std::ostringstream oss;
      oss << "StripeCtrlAlphaGenerator: infeasible k_target=" << k_target << " > n_r*n_s";
      detail::SetErr(err, oss.str());
      return false;
    }
#else
    // Conservative fallback check using long double.
    const long double max_k = static_cast<long double>(spec.n_r) * static_cast<long double>(spec.n_s);
    if (static_cast<long double>(k_target) > max_k) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: infeasible k_target > n_r*n_s");
      return false;
    }
#endif

    // Sample degrees d_i with sum k_target.
    Rng rng(spec.seed);
    std::vector<u32> degrees;
    if (!detail::RandomBoundedComposition(k_target, spec.n_r, spec.n_s, &rng, &degrees, err)) return false;

    // Strip positions on control axis.
    std::vector<double> strip_lo(spec.n_s);
    std::vector<double> strip_hi(spec.n_s);
    for (u64 j = 0; j < spec.n_s; ++j) {
      const double yb = dom_lo + g + static_cast<double>(j) * (h + g);
      const double yt = yb + h;
      strip_lo[static_cast<usize>(j)] = yb;
      strip_hi[static_cast<usize>(j)] = yt;
    }

    auto sample_core_interval = [&](int axis, double* out_lo, double* out_hi) {
      SJS_ASSERT(out_lo && out_hi);
      SJS_ASSERT(axis != control_axis);
      const double lo = rng.UniformDouble(dom_lo, core_lo);   // <= core_lo
      const double hi = rng.UniformDouble(core_hi, dom_hi);   // >= core_hi
      *out_lo = lo;
      *out_hi = hi;
    };

    // Build S: disjoint strips on control axis, and "core-covering" intervals on other axes.
    Relation<Dim, T> S;
    S.name = "S";
    S.boxes.resize(static_cast<usize>(spec.n_s));
    S.ids.resize(static_cast<usize>(spec.n_s));

    for (u64 j = 0; j < spec.n_s; ++j) {
      BoxT b;
      for (int axis = 0; axis < Dim; ++axis) {
        const usize a = static_cast<usize>(axis);
        if (axis == control_axis) {
          b.lo.v[a] = static_cast<T>(strip_lo[static_cast<usize>(j)]);
          b.hi.v[a] = static_cast<T>(strip_hi[static_cast<usize>(j)]);
        } else {
          double lo, hi;
          sample_core_interval(axis, &lo, &hi);
          b.lo.v[a] = static_cast<T>(lo);
          b.hi.v[a] = static_cast<T>(hi);
        }
      }
      S.boxes[static_cast<usize>(j)] = b;
      S.ids[static_cast<usize>(j)] = static_cast<Id>(j);
    }

    if (shuffle_strips) {
      detail::ShuffleRelation(&S, &rng);
    }

    // Build R with degrees.
    Relation<Dim, T> R;
    R.name = "R";
    R.boxes.resize(static_cast<usize>(spec.n_r));
    R.ids.resize(static_cast<usize>(spec.n_r));

    for (u64 i = 0; i < spec.n_r; ++i) {
      const u64 di = static_cast<u64>(degrees[static_cast<usize>(i)]);
      BoxT b;

      // Non-control axes: cover the core.
      for (int axis = 0; axis < Dim; ++axis) {
        if (axis == control_axis) continue;
        double lo, hi;
        sample_core_interval(axis, &lo, &hi);
        const usize a = static_cast<usize>(axis);
        b.lo.v[a] = static_cast<T>(lo);
        b.hi.v[a] = static_cast<T>(hi);
      }

      // Control axis interval decides how many strips it intersects.
      const usize ca = static_cast<usize>(control_axis);

      if (di == 0) {
        // Place inside a random gap u in [0..n_s].
        const u64 u = rng.UniformU64(spec.n_s + 1);

        double gap_lo = 0.0, gap_hi = 0.0;
        if (u == 0) {
          gap_lo = dom_lo;
          gap_hi = dom_lo + g;
        } else if (u == spec.n_s) {
          gap_lo = dom_hi - g;
          gap_hi = dom_hi;
        } else {
          gap_lo = strip_hi[static_cast<usize>(u - 1)];
          gap_hi = strip_lo[static_cast<usize>(u)];
        }

        // Put a tiny interval [y0, y0+delta) strictly inside the gap (with margins).
        const double lo_y = gap_lo + delta;
        const double hi_y = gap_hi - 2.0 * delta;
        if (!(hi_y > lo_y)) {
          detail::SetErr(err, "StripeCtrlAlphaGenerator: gap too small for delta (reduce delta_factor or increase g)");
          return false;
        }
        const double y0 = rng.UniformDouble(lo_y, hi_y);
        const double y1 = y0 + delta;
        b.lo.v[ca] = static_cast<T>(y0);
        b.hi.v[ca] = static_cast<T>(y1);
      } else {
        if (di > spec.n_s) {
          detail::SetErr(err, "StripeCtrlAlphaGenerator: degree di > n_s (bug)");
          return false;
        }
        const u64 max_start = spec.n_s - di;
        const u64 s = rng.UniformU64(max_start + 1);  // start strip index
        const u64 e = s + di - 1;

        const double y0 = strip_lo[static_cast<usize>(s)] + delta;
        const double y1 = strip_hi[static_cast<usize>(e)] - delta;
        if (!(y1 > y0)) {
          detail::SetErr(err, "StripeCtrlAlphaGenerator: invalid y interval (delta too large?)");
          return false;
        }
        b.lo.v[ca] = static_cast<T>(y0);
        b.hi.v[ca] = static_cast<T>(y1);
      }

      R.boxes[static_cast<usize>(i)] = b;
      R.ids[static_cast<usize>(i)] = static_cast<Id>(i);
    }

    if (shuffle_r) {
      detail::ShuffleRelation(&R, &rng);
    }

    DatasetT ds;
    ds.name = spec.name;
    ds.half_open = true;

    if (!swap_sides) {
      ds.R = std::move(R);
      ds.S = std::move(S);
    } else {
      ds.R = std::move(S);
      ds.S = std::move(R);
      std::swap(ds.R.name, ds.S.name);
    }

    // Validate
    {
      std::string tmp;
      if (!ds.Validate(/*require_proper=*/true, &tmp)) {
        detail::SetErr(err, "StripeCtrlAlphaGenerator produced invalid dataset: " + tmp);
        return false;
      }
    }

    if (report) {
      report->generator = std::string(Name());
      report->dataset_name = spec.name;
      report->n_r = spec.n_r;
      report->n_s = spec.n_s;

      report->has_exact_k = true;
      report->k_target = k_target;
      report->k_achieved = k_target;
      report->alpha_target = spec.alpha;
      report->alpha_achieved = static_cast<double>(k_target) / static_cast<double>(spec.n_r + spec.n_s);

      std::ostringstream notes;
      notes << "control_axis=" << control_axis
            << ", core=[" << core_lo_frac << "," << core_hi_frac << "]"
            << ", g=" << g << ", h=" << h << ", delta=" << delta
            << (shuffle_strips ? ", shuffle_strips=1" : ", shuffle_strips=0")
            << (shuffle_r ? ", shuffle_r=1" : ", shuffle_r=0")
            << (k_overridden ? ", k_overridden=1" : ", k_overridden=0");
      report->notes = notes.str();
    }

    *out_ds = std::move(ds);
    return true;
  }
};

}  // namespace synthetic
}  // namespace sjs
