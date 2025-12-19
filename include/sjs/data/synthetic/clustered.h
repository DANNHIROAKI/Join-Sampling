#pragma once
// sjs/data/synthetic/clustered.h
//
// Clustered rectangles (hotspots).
//
// Intuition:
//  - Choose K cluster centers uniformly in the domain.
//  - Each rectangle samples a cluster id uniformly and jitters its center around
//    that cluster by a Gaussian with stddev sigma.
//  - Rectangle sizes are uniform in [w_min, w_max] fraction of domain length.
//
// This generator is useful to stress algorithms under skew and local density.
//
// Parameters (spec.params):
//  - "clusters" (int, default 10): number of clusters.
//  - "sigma" (double, default 0.05): stddev as fraction of domain length.
//  - "w_min" / "w_max" (double, default 0.003..0.02): width fraction range.
//    Side-specific overrides: r_w_min/r_w_max and s_w_min/s_w_max.
//  - "share_clusters" (bool, default true): if true, R and S share the same cluster centers.
//  - "shuffle_r" / "shuffle_s" (bool, default false).
//
// Notes:
//  - Uses a deterministic Box-Muller transform driven by sjs::Rng.

#include "sjs/data/synthetic/generator.h"

#include <cmath>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace synthetic {
namespace detail {

class Normal01 {
 public:
  explicit Normal01(Rng* rng) : rng_(rng) { SJS_ASSERT(rng_ != nullptr); }

  // Standard normal N(0,1).
  double Next() {
    if (has_spare_) {
      has_spare_ = false;
      return spare_;
    }
    // Box-Muller (polar form)
    while (true) {
      const double u = 2.0 * rng_->NextDouble() - 1.0;
      const double v = 2.0 * rng_->NextDouble() - 1.0;
      const double s = u * u + v * v;
      if (s <= 0.0 || s >= 1.0) continue;
      const double mul = std::sqrt(-2.0 * std::log(s) / s);
      spare_ = v * mul;
      has_spare_ = true;
      return u * mul;
    }
  }

 private:
  Rng* rng_;
  bool has_spare_{false};
  double spare_{0.0};
};

template <int Dim, class T>
inline void ClampPointToDomain(Point<Dim, T>* p, double dom_lo, double dom_hi) {
  for (int i = 0; i < Dim; ++i) {
    const usize idx = static_cast<usize>(i);
    if (p->v[idx] < static_cast<T>(dom_lo)) p->v[idx] = static_cast<T>(dom_lo);
    if (p->v[idx] > static_cast<T>(dom_hi)) p->v[idx] = static_cast<T>(dom_hi);
  }
}

// Make a box centered at c with width w along each axis, clipped to domain while preserving width when possible.
template <int Dim, class T>
inline Box<Dim, T> BoxFromCenterWidth(const Point<Dim, T>& c, const Point<Dim, T>& w,
                                      double dom_lo, double dom_hi) {
  Box<Dim, T> b;
  for (int axis = 0; axis < Dim; ++axis) {
    const usize axis_idx = static_cast<usize>(axis);
    const double width = static_cast<double>(w.v[axis_idx]);
    const double half = width * 0.5;
    double lo = static_cast<double>(c.v[axis_idx]) - half;
    double hi = static_cast<double>(c.v[axis_idx]) + half;

    // If width exceeds domain, clamp to full domain.
    if (width >= (dom_hi - dom_lo)) {
      lo = dom_lo;
      hi = dom_hi;
    } else {
      // Shift to fit.
      if (lo < dom_lo) {
        const double shift = dom_lo - lo;
        lo += shift;
        hi += shift;
      }
      if (hi > dom_hi) {
        const double shift = hi - dom_hi;
        lo -= shift;
        hi -= shift;
      }
      // Final clamp for numerical safety.
      if (lo < dom_lo) lo = dom_lo;
      if (hi > dom_hi) hi = dom_hi;
    }

    // Ensure proper interval; if degenerates due to numeric issues, nudge.
    if (!(hi > lo)) {
      hi = std::nextafter(lo, dom_hi);
      if (!(hi > lo)) hi = lo + 1e-12;
    }

    b.lo.v[axis_idx] = static_cast<T>(lo);
    b.hi.v[axis_idx] = static_cast<T>(hi);
  }
  return b;
}

}  // namespace detail

template <int Dim, class T = Scalar>
class ClusteredGenerator final : public ISyntheticGenerator<Dim, T> {
 public:
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;
  using PointT = Point<Dim, T>;

  std::string_view Name() const noexcept override { return "clustered"; }

  bool Generate(const DatasetSpec& spec,
                DatasetT* out_ds,
                Report* report,
                std::string* err) override {
    if (!out_ds) {
      detail::SetErr(err, "ClusteredGenerator: out_ds is null");
      return false;
    }
    if (spec.n_r == 0 || spec.n_s == 0) {
      detail::SetErr(err, "ClusteredGenerator: n_r and n_s must be > 0");
      return false;
    }
    const double dom_lo = spec.domain_lo;
    const double dom_hi = spec.domain_hi;
    if (!(dom_hi > dom_lo)) {
      detail::SetErr(err, "ClusteredGenerator: domain_hi must be > domain_lo");
      return false;
    }
    const double L = dom_hi - dom_lo;

    const i32 K = detail::GetI32(spec.params, "clusters", 10);
    if (K <= 0) {
      detail::SetErr(err, "ClusteredGenerator: clusters must be > 0");
      return false;
    }
    const double sigma_frac = detail::GetDouble(spec.params, "sigma", 0.05);
    if (!(sigma_frac >= 0.0 && sigma_frac < 1.0)) {
      detail::SetErr(err, "ClusteredGenerator: sigma must be in [0,1)");
      return false;
    }
    const double sigma = sigma_frac * L;

    const bool share_clusters = detail::GetBool(spec.params, "share_clusters", true);
    const bool shuffle_r = detail::GetBool(spec.params, "shuffle_r", false);
    const bool shuffle_s = detail::GetBool(spec.params, "shuffle_s", false);

    // Cluster centers.
    std::vector<PointT> centers;
    centers.resize(static_cast<usize>(K));
    {
      Rng rng_centers(spec.seed ^ 0x13579BDFULL);
      for (int k = 0; k < K; ++k) {
        PointT c;
        for (int axis = 0; axis < Dim; ++axis) {
          const usize axis_idx = static_cast<usize>(axis);
          c.v[axis_idx] = static_cast<T>(rng_centers.UniformDouble(dom_lo, dom_hi));
        }
        centers[static_cast<usize>(k)] = c;
      }
    }

    auto gen_relation = [&](u64 n,
                            std::string_view prefix,
                            Relation<Dim, T>* out_rel,
                            u64 seed_salt,
                            const std::vector<PointT>& use_centers) -> bool {
      const double w_min_frac = detail::GetDoubleSide(spec.params, prefix, "w_min", 0.003);
      const double w_max_frac = detail::GetDoubleSide(spec.params, prefix, "w_max", 0.02);
      if (!(w_min_frac > 0.0 && w_max_frac >= w_min_frac && w_max_frac < 1.0)) {
        detail::SetErr(err, "ClusteredGenerator: invalid width fraction range");
        return false;
      }

      out_rel->boxes.resize(static_cast<usize>(n));
      out_rel->ids.resize(static_cast<usize>(n));

      Rng rng(spec.seed ^ seed_salt);
      detail::Normal01 normal(&rng);

      for (u64 i = 0; i < n; ++i) {
        const u64 cid = rng.UniformU64(static_cast<u64>(K));
        PointT center = use_centers[static_cast<usize>(cid)];

        // Jitter
        for (int axis = 0; axis < Dim; ++axis) {
          const usize axis_idx = static_cast<usize>(axis);
          const double jitter = sigma * normal.Next();
          double v = static_cast<double>(center.v[axis_idx]) + jitter;
          if (v < dom_lo) v = dom_lo;
          if (v > dom_hi) v = dom_hi;
          center.v[axis_idx] = static_cast<T>(v);
        }

        // Width per axis
        PointT w;
        for (int axis = 0; axis < Dim; ++axis) {
          const usize axis_idx = static_cast<usize>(axis);
          const double wf = rng.UniformDouble(w_min_frac, w_max_frac);
          w.v[axis_idx] = static_cast<T>(wf * L);
        }

        BoxT b = detail::BoxFromCenterWidth<Dim, T>(center, w, dom_lo, dom_hi);
        out_rel->boxes[static_cast<usize>(i)] = b;
        out_rel->ids[static_cast<usize>(i)] = static_cast<Id>(i);
      }

      return true;
    };

    Relation<Dim, T> R, S;
    R.name = "R";
    S.name = "S";

    if (!gen_relation(spec.n_r, "r_", &R, 0x2468ACE0ULL, centers)) return false;

    if (share_clusters) {
      if (!gen_relation(spec.n_s, "s_", &S, 0x0ECA8642ULL, centers)) return false;
    } else {
      // Independent S centers
      std::vector<PointT> centers_s;
      centers_s.resize(static_cast<usize>(K));
      Rng rng_centers_s(spec.seed ^ 0xCAFEBABEU);
      for (int k = 0; k < K; ++k) {
        PointT c;
        for (int axis = 0; axis < Dim; ++axis) {
          const usize axis_idx = static_cast<usize>(axis);
          c.v[axis_idx] = static_cast<T>(rng_centers_s.UniformDouble(dom_lo, dom_hi));
        }
        centers_s[static_cast<usize>(k)] = c;
      }
      if (!gen_relation(spec.n_s, "s_", &S, 0x0ECA8642ULL, centers_s)) return false;
    }

    if (shuffle_r) {
      Rng rng_shuffle(spec.seed + 17);
      std::vector<usize> perm(R.Size());
      for (usize i = 0; i < perm.size(); ++i) perm[i] = i;
      detail::ShuffleInPlace(&perm, &rng_shuffle);
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

    if (shuffle_s) {
      Rng rng_shuffle(spec.seed + 31);
      std::vector<usize> perm(S.Size());
      for (usize i = 0; i < perm.size(); ++i) perm[i] = i;
      detail::ShuffleInPlace(&perm, &rng_shuffle);
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

    DatasetT ds;
    ds.name = spec.name;
    ds.half_open = true;
    ds.R = std::move(R);
    ds.S = std::move(S);

    // Validate
    {
      std::string tmp;
      if (!ds.Validate(/*require_proper=*/true, &tmp)) {
        detail::SetErr(err, "ClusteredGenerator produced invalid dataset: " + tmp);
        return false;
      }
    }

    if (report) {
      report->generator = std::string(Name());
      report->dataset_name = spec.name;
      report->n_r = spec.n_r;
      report->n_s = spec.n_s;
      report->has_exact_k = false;
      report->alpha_target = spec.alpha;
      report->alpha_achieved = std::numeric_limits<double>::quiet_NaN();

      std::ostringstream notes;
      notes << "clusters=" << K << ", sigma=" << sigma_frac
            << ", w_min_r=" << detail::GetDoubleSide(spec.params, "r_", "w_min", 0.003)
            << ", w_max_r=" << detail::GetDoubleSide(spec.params, "r_", "w_max", 0.02)
            << ", share_clusters=" << (share_clusters ? 1 : 0);
      report->notes = notes.str();
    }

    *out_ds = std::move(ds);
    return true;
  }
};

}  // namespace synthetic
}  // namespace sjs
