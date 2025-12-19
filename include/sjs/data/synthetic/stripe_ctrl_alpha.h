#pragma once
// sjs/data/synthetic/stripe_ctrl_alpha.h
//
// Stripe-controlled alpha generator.
//
// This implements the construction described in 数据生成.md (SCIRG):
//  - All boxes share overlap in all axes except one "control axis".
//  - On the control axis, S forms non-overlapping strips separated by gaps.
//  - Each R box spans exactly d_i consecutive strips (or lies in a gap when d_i=0).
//  - By sampling a feasible degree sequence {d_i} with sum k, we can control
//    the exact number of intersecting pairs |J| = k.
//
// Default choices match the markdown defaults on domain [0,1):
//  - core interval: [0.45, 0.55] (fraction of domain length)
//  - gap_factor: 0.1  (so total gap length is 10% of domain, split into n_s+1 gaps)
//  - delta_factor: 0.25 (safety margin relative to min(gap, strip_height))
//
// Key parameters (spec.params):
//  - "control_axis" (int, default 1): which axis to use as strips (y in 2D).
//  - "core_lo" / "core_hi" (double in [0,1], default 0.45/0.55):
//      core interval fractions for non-control axes.
//  - "gap_factor" (double in (0,1), default 0.1):
//      total gaps occupy gap_factor * domain_length.
//  - "delta_factor" (double in (0,0.5), default 0.25):
//      delta = min(gap, h) * delta_factor.
//  - "shuffle_strips" (bool, default true): randomly permute S order.
//  - "shuffle_r" (bool, default false): randomly permute R order.
//  - "k_target" (u64, optional): override k target directly; otherwise
//      k = round(alpha * (n_r + n_s)) as in the markdown.
//
// IMPORTANT: Here alpha is defined as k / (n_r + n_s) (per the uploaded markdown),
// not k / (n_r*n_s). This is intentional to match that generator spec.

#include "sjs/data/synthetic/generator.h"
#include "sjs/core/logging.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace synthetic {

namespace detail {

// Sample a random bounded composition (d_0..d_{n_r-1}) s.t.
//   0<=d_i<=n_s and sum d_i = k.
// The distribution matches the markdown: sequentially sample d_i uniformly
// from feasible [low, high] given remaining.
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

    // low = max(0, remaining - left*n_s)
    __uint128_t max_future = static_cast<__uint128_t>(left) * static_cast<__uint128_t>(n_s);
    u64 low = 0;
    if (static_cast<__uint128_t>(remaining) > max_future) {
      const __uint128_t diff = static_cast<__uint128_t>(remaining) - max_future;
      if (diff > static_cast<__uint128_t>(std::numeric_limits<u64>::max())) {
        SetErr(err, "RandomBoundedComposition: overflow computing low");
        return false;
      }
      low = static_cast<u64>(diff);
    }

    const u64 high = (remaining < n_s) ? remaining : n_s;
    if (low > high) {
      SetErr(err, "RandomBoundedComposition: infeasible bounds (low>high).");
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
  if (!(x >= 0.0)) return 0;
  const long double y = static_cast<long double>(x);
  const long double r = std::llround(y);
  if (r <= 0) return 0;
  if (r >= static_cast<long double>(std::numeric_limits<u64>::max())) {
    return std::numeric_limits<u64>::max();
  }
  return static_cast<u64>(r);
}

}  // namespace detail

template <int Dim, class T = Scalar>
class StripeCtrlAlphaGenerator final : public ISyntheticGenerator<Dim, T> {
 public:
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;
  using PointT = Point<Dim, T>;

  std::string_view Name() const noexcept override { return "stripe_ctrl_alpha"; }

  bool Generate(const DatasetSpec& spec,
                DatasetT* out_ds,
                Report* report,
                std::string* err) override {
    if (!out_ds) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: out_ds is null");
      return false;
    }
    if (spec.n_r == 0 || spec.n_s == 0) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: n_r and n_s must be > 0");
      return false;
    }
    const double dom_lo = spec.domain_lo;
    const double dom_hi = spec.domain_hi;
    if (!(dom_hi > dom_lo)) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: domain_hi must be > domain_lo");
      return false;
    }
    const double L = dom_hi - dom_lo;

    // Parameters
    const i32 control_axis = detail::GetI32(spec.params, "control_axis", 1);
    if (control_axis < 0 || control_axis >= Dim) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: invalid control_axis");
      return false;
    }
    if (Dim < 2) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator requires Dim >= 2");
      return false;
    }

    const double core_lo_frac = detail::GetDouble(spec.params, "core_lo", 0.45);
    const double core_hi_frac = detail::GetDouble(spec.params, "core_hi", 0.55);
    if (!(core_lo_frac >= 0.0 && core_hi_frac <= 1.0 && core_lo_frac < core_hi_frac)) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: core_lo/core_hi must satisfy 0<=lo<hi<=1");
      return false;
    }
    const double core_lo = dom_lo + core_lo_frac * L;
    const double core_hi = dom_lo + core_hi_frac * L;

    // Gap/strip geometry on control axis.
    // Either accept an absolute gap size "g", or compute from gap_factor/(n_s+1).
    double g = 0.0;
    if (auto v = detail::FindParam(spec.params, "g")) {
      g = detail::GetDouble(spec.params, "g", 0.0);
    } else if (auto v2 = detail::FindParam(spec.params, "gap")) {
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
      detail::SetErr(err, "StripeCtrlAlphaGenerator: gap size g must be > 0");
      return false;
    }
    if (!(static_cast<double>(spec.n_s + 1) * g < L)) {
      detail::SetErr(err, "StripeCtrlAlphaGenerator: (n_s+1)*g must be < domain_length");
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
      detail::SetErr(err, "StripeCtrlAlphaGenerator: delta must satisfy 0<delta<min(g/2,h/2)");
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
      // As in the markdown: k = round(alpha * (n_r + n_s)).
      k_target = detail::RoundToU64(spec.alpha * static_cast<double>(spec.n_r + spec.n_s));
    }

    // Feasibility: 0 <= k <= n_r*n_s.
    const __uint128_t max_k128 =
        static_cast<__uint128_t>(spec.n_r) * static_cast<__uint128_t>(spec.n_s);
    if (static_cast<__uint128_t>(k_target) > max_k128) {
      std::ostringstream oss;
      oss << "StripeCtrlAlphaGenerator: infeasible k_target=" << k_target
          << " > n_r*n_s=" << static_cast<unsigned long long>(max_k128)
          << ". (alpha too large?)";
      detail::SetErr(err, oss.str());
      return false;
    }

    // Sample degrees d_i (size n_r) summing to k_target.
    Rng rng(spec.seed);
    std::vector<u32> degrees;
    if (!detail::RandomBoundedComposition(k_target, spec.n_r, spec.n_s, &rng, &degrees, err)) {
      return false;
    }

    // Prepare y positions of strips (on control axis).
    // For strip j in [0..n_s-1]:
    //   yb = dom_lo + g + j*(h+g)
    //   yt = yb + h
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
      if (axis == control_axis) {
        SJS_UNREACHABLE();
      }
      const double lo = rng.UniformDouble(dom_lo, core_lo);
      const double hi = rng.UniformDouble(core_hi, dom_hi);
      *out_lo = lo;
      *out_hi = hi;
    };

    // Build S (strips)
    Relation<Dim, T> S;
    S.name = "S";
    S.boxes.resize(static_cast<usize>(spec.n_s));
    S.ids.resize(static_cast<usize>(spec.n_s));

    for (u64 j = 0; j < spec.n_s; ++j) {
      BoxT b;
      for (int axis = 0; axis < Dim; ++axis) {
        const usize axis_idx = static_cast<usize>(axis);
        if (axis == control_axis) {
          b.lo.v[axis_idx] = static_cast<T>(strip_lo[static_cast<usize>(j)]);
          b.hi.v[axis_idx] = static_cast<T>(strip_hi[static_cast<usize>(j)]);
        } else {
          double lo, hi;
          sample_core_interval(axis, &lo, &hi);
          b.lo.v[axis_idx] = static_cast<T>(lo);
          b.hi.v[axis_idx] = static_cast<T>(hi);
        }
      }
      S.boxes[static_cast<usize>(j)] = b;
      S.ids[static_cast<usize>(j)] = static_cast<Id>(j);
    }

    if (shuffle_strips) {
      // Shuffle S boxes and ids together.
      std::vector<usize> perm(spec.n_s);
      for (u64 i = 0; i < spec.n_s; ++i) perm[static_cast<usize>(i)] = static_cast<usize>(i);
      detail::ShuffleInPlace(&perm, &rng);

      Relation<Dim, T> S2;
      S2.name = S.name;
      S2.boxes.resize(S.boxes.size());
      S2.ids.resize(S.ids.size());
      for (usize i = 0; i < perm.size(); ++i) {
        S2.boxes[i] = S.boxes[perm[i]];
        S2.ids[i] = S.ids[perm[i]];
      }
      S = std::move(S2);
    }

    // Build R
    Relation<Dim, T> R;
    R.name = "R";
    R.boxes.resize(static_cast<usize>(spec.n_r));
    R.ids.resize(static_cast<usize>(spec.n_r));

    for (u64 i = 0; i < spec.n_r; ++i) {
      const u64 di = static_cast<u64>(degrees[static_cast<usize>(i)]);
      BoxT b;

      // Non-control axes: sample intervals that cover the core.
      for (int axis = 0; axis < Dim; ++axis) {
        if (axis == control_axis) continue;
        double lo, hi;
        sample_core_interval(axis, &lo, &hi);
        const usize axis_idx = static_cast<usize>(axis);
        b.lo.v[axis_idx] = static_cast<T>(lo);
        b.hi.v[axis_idx] = static_cast<T>(hi);
      }

      // Control axis interval.
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
          // gap between strip (u-1) and strip u
          gap_lo = strip_hi[static_cast<usize>(u - 1)];
          gap_hi = strip_lo[static_cast<usize>(u)];
        }

        // Choose a tiny interval [y, y+delta) inside the gap with safety margins.
        const double lo_y = gap_lo + delta;
        const double hi_y = gap_hi - 2.0 * delta;
        if (!(hi_y > lo_y)) {
          detail::SetErr(err, "StripeCtrlAlphaGenerator: gap too small for chosen delta");
          return false;
        }
        const double y0 = rng.UniformDouble(lo_y, hi_y);
        const double y1 = y0 + delta;
        const usize control_axis_idx = static_cast<usize>(control_axis);
        b.lo.v[control_axis_idx] = static_cast<T>(y0);
        b.hi.v[control_axis_idx] = static_cast<T>(y1);
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
        const usize control_axis_idx = static_cast<usize>(control_axis);
        b.lo.v[control_axis_idx] = static_cast<T>(y0);
        b.hi.v[control_axis_idx] = static_cast<T>(y1);
      }

      R.boxes[static_cast<usize>(i)] = b;
      R.ids[static_cast<usize>(i)] = static_cast<Id>(i);
    }

    if (shuffle_r) {
      // Shuffle R boxes and ids together.
      std::vector<usize> perm(spec.n_r);
      for (u64 i = 0; i < spec.n_r; ++i) perm[static_cast<usize>(i)] = static_cast<usize>(i);
      detail::ShuffleInPlace(&perm, &rng);

      Relation<Dim, T> R2;
      R2.name = R.name;
      R2.boxes.resize(R.boxes.size());
      R2.ids.resize(R.ids.size());
      for (usize i = 0; i < perm.size(); ++i) {
        R2.boxes[i] = R.boxes[perm[i]];
        R2.ids[i] = R.ids[perm[i]];
      }
      R = std::move(R2);
    }

    // Assemble dataset.
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

    // Ensure proper boxes (debug).
    {
      std::string tmp;
      if (!ds.Validate(/*require_proper=*/true, &tmp)) {
        detail::SetErr(err, "StripeCtrlAlphaGenerator produced invalid dataset: " + tmp);
        return false;
      }
    }

    // Fill report.
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
