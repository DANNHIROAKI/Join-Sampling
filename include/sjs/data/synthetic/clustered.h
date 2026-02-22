#pragma once
// sjs/data/synthetic/clustered.h
//
// Clustered boxes (hotspots).
//
// Intuition:
//  - Choose K cluster centers uniformly in the domain.
//  - Each box chooses a cluster id uniformly and jitters its center around it by Gaussian N(0, sigma^2).
//  - Box sizes are uniform in [w_min, w_max] (fractions of domain length).
//
// Parameters (spec.params):
//  - "clusters" (int, default 10)
//  - "sigma" (double, default 0.05): as fraction of domain length.
//  - "w_min" / "w_max" (default 0.003..0.02), side-specific overrides r_*, s_*.
//  - "share_clusters" (bool, default true): if true, R and S share the same cluster centers.
//  - "shuffle_r" / "shuffle_s" (bool, default false).
//
// Uses a deterministic Box-Muller transform driven by sjs::Rng.

#include "sjs/data/synthetic/generator.h"

#include <algorithm>
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

  double Next() {
    if (has_spare_) {
      has_spare_ = false;
      return spare_;
    }
    // Polar Box-Muller.
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
inline Box<Dim, T> BoxFromCenterWidth(const Point<Dim, T>& c,
                                      const Point<Dim, T>& w,
                                      double dom_lo,
                                      double dom_hi) {
  Box<Dim, T> b;
  const double L = dom_hi - dom_lo;

  for (int axis = 0; axis < Dim; ++axis) {
    const usize a = static_cast<usize>(axis);
    const double width = static_cast<double>(w.v[a]);
    const double half = 0.5 * width;

    double lo = static_cast<double>(c.v[a]) - half;
    double hi = static_cast<double>(c.v[a]) + half;

    if (width >= L) {
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
      // Final clamp.
      if (lo < dom_lo) lo = dom_lo;
      if (hi > dom_hi) hi = dom_hi;
    }

    if (!(hi > lo)) {
      hi = std::nextafter(lo, dom_hi);
      if (!(hi > lo)) hi = lo + 1e-12;
    }

    b.lo.v[a] = static_cast<T>(lo);
    b.hi.v[a] = static_cast<T>(hi);
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

  bool Generate(const DatasetSpec& spec, DatasetT* out_ds, Report* report, std::string* err) override {
    if (!out_ds) {
      detail::SetErr(err, "ClusteredGenerator: out_ds is null");
      return false;
    }
    if (spec.n_r == 0 || spec.n_s == 0) {
      detail::SetErr(err, "ClusteredGenerator: n_r and n_s must be > 0");
      return false;
    }
    if (!detail::CheckIdFits(spec.n_r, err, "ClusteredGenerator(R)")) return false;
    if (!detail::CheckIdFits(spec.n_s, err, "ClusteredGenerator(S)")) return false;

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

    // Cluster centers for R (and optionally S).
    std::vector<PointT> centers;
    centers.resize(static_cast<usize>(K));
    {
      Rng rng_centers(DeriveSeed(spec.seed, 0x13579BDFULL));
      for (int k = 0; k < K; ++k) {
        PointT c;
        for (int axis = 0; axis < Dim; ++axis) {
          const usize a = static_cast<usize>(axis);
          c.v[a] = static_cast<T>(rng_centers.UniformDouble(dom_lo, dom_hi));
        }
        centers[static_cast<usize>(k)] = c;
      }
    }

    auto gen_relation = [&](u64 n,
                            std::string_view prefix,
                            Relation<Dim, T>* out_rel,
                            u64 salt,
                            const std::vector<PointT>& use_centers) -> bool {
      const double w_min_frac = detail::GetDoubleSide(spec.params, prefix, "w_min", 0.003);
      const double w_max_frac = detail::GetDoubleSide(spec.params, prefix, "w_max", 0.02);
      if (!(w_min_frac > 0.0 && w_max_frac >= w_min_frac && w_max_frac < 1.0)) {
        detail::SetErr(err, "ClusteredGenerator: invalid width fraction range (need 0<w_min<=w_max<1)");
        return false;
      }

      out_rel->boxes.resize(static_cast<usize>(n));
      out_rel->ids.resize(static_cast<usize>(n));

      Rng rng(DeriveSeed(spec.seed, salt));
      detail::Normal01 normal(&rng);

      for (u64 i = 0; i < n; ++i) {
        const u64 cid = rng.UniformU64(static_cast<u64>(K));
        PointT center = use_centers[static_cast<usize>(cid)];

        // Jitter center by Gaussian.
        for (int axis = 0; axis < Dim; ++axis) {
          const usize a = static_cast<usize>(axis);
          double v = static_cast<double>(center.v[a]) + sigma * normal.Next();
          if (v < dom_lo) v = dom_lo;
          if (v > dom_hi) v = dom_hi;
          center.v[a] = static_cast<T>(v);
        }

        // Width per axis.
        PointT w;
        for (int axis = 0; axis < Dim; ++axis) {
          const usize a = static_cast<usize>(axis);
          w.v[a] = static_cast<T>(rng.UniformDouble(w_min_frac, w_max_frac) * L);
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
      std::vector<PointT> centers_s;
      centers_s.resize(static_cast<usize>(K));
      Rng rng_centers_s(DeriveSeed(spec.seed, 0xCAFEBABEU));
      for (int k = 0; k < K; ++k) {
        PointT c;
        for (int axis = 0; axis < Dim; ++axis) {
          const usize a = static_cast<usize>(axis);
          c.v[a] = static_cast<T>(rng_centers_s.UniformDouble(dom_lo, dom_hi));
        }
        centers_s[static_cast<usize>(k)] = c;
      }
      if (!gen_relation(spec.n_s, "s_", &S, 0x0ECA8642ULL, centers_s)) return false;
    }

    if (shuffle_r) {
      Rng rng_shuf(DeriveSeed(spec.seed, 0x33333333ULL));
      detail::ShuffleRelation(&R, &rng_shuf);
    }
    if (shuffle_s) {
      Rng rng_shuf(DeriveSeed(spec.seed, 0x44444444ULL));
      detail::ShuffleRelation(&S, &rng_shuf);
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
      notes << "clusters=" << K
            << ", sigma=" << sigma_frac
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
