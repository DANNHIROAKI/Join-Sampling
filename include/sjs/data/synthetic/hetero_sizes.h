#pragma once
// sjs/data/synthetic/hetero_sizes.h
//
// Heterogeneous rectangle sizes (mixture of small and large).
//
// Motivation:
//  - Stress "active set" / pruning behavior: a few very large boxes can
//    create a large active set and many candidate checks.
//  - More realistic than purely uniform sizes.
//
// Parameters (spec.params):
//  - "p_large" (double, default 0.1): probability of sampling a "large" box.
//    Side-specific: r_p_large, s_p_large.
//  - Small size range (fractions of domain length):
//      "w_small_min" / "w_small_max" (default 0.002..0.01)
//      Side-specific: r_w_small_min, ...
//  - Large size range (fractions of domain length):
//      "w_large_min" / "w_large_max" (default 0.1..0.3)
//      Side-specific: r_w_large_min, ...
//  - "anisotropic" (bool, default true):
//      If true, when sampling a large box, pick ONE axis to be large and
//      keep other axes small (long-thin rectangles). In 2D, that yields
//      vertical/horizontal slabs.
//  - "shuffle_r" / "shuffle_s" (bool, default false).

#include "sjs/data/synthetic/generator.h"

#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace synthetic {

template <int Dim, class T = Scalar>
class HeteroSizesGenerator final : public ISyntheticGenerator<Dim, T> {
 public:
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;

  std::string_view Name() const noexcept override { return "hetero_sizes"; }

  bool Generate(const DatasetSpec& spec,
                DatasetT* out_ds,
                Report* report,
                std::string* err) override {
    if (!out_ds) {
      detail::SetErr(err, "HeteroSizesGenerator: out_ds is null");
      return false;
    }
    if (spec.n_r == 0 || spec.n_s == 0) {
      detail::SetErr(err, "HeteroSizesGenerator: n_r and n_s must be > 0");
      return false;
    }
    const double dom_lo = spec.domain_lo;
    const double dom_hi = spec.domain_hi;
    if (!(dom_hi > dom_lo)) {
      detail::SetErr(err, "HeteroSizesGenerator: domain_hi must be > domain_lo");
      return false;
    }
    const double L = dom_hi - dom_lo;

    const bool anisotropic = detail::GetBool(spec.params, "anisotropic", true);
    const bool shuffle_r = detail::GetBool(spec.params, "shuffle_r", false);
    const bool shuffle_s = detail::GetBool(spec.params, "shuffle_s", false);

    auto gen_relation = [&](u64 n, std::string_view prefix, Relation<Dim, T>* out_rel, u64 seed_salt) -> bool {
      const double p_large = detail::GetDoubleSide(spec.params, prefix, "p_large", 0.1);
      if (!(p_large >= 0.0 && p_large <= 1.0)) {
        detail::SetErr(err, "HeteroSizesGenerator: p_large must be in [0,1]");
        return false;
      }

      const double ws_min = detail::GetDoubleSide(spec.params, prefix, "w_small_min", 0.002);
      const double ws_max = detail::GetDoubleSide(spec.params, prefix, "w_small_max", 0.01);
      const double wl_min = detail::GetDoubleSide(spec.params, prefix, "w_large_min", 0.1);
      const double wl_max = detail::GetDoubleSide(spec.params, prefix, "w_large_max", 0.3);

      auto check_range = [&](double a, double b, const char* name) -> bool {
        if (!(a > 0.0 && b >= a && b < 1.0)) {
          detail::SetErr(err, std::string("HeteroSizesGenerator: invalid range for ") + name);
          return false;
        }
        return true;
      };
      if (!check_range(ws_min, ws_max, "w_small_*")) return false;
      if (!check_range(wl_min, wl_max, "w_large_*")) return false;

      out_rel->boxes.resize(static_cast<usize>(n));
      out_rel->ids.resize(static_cast<usize>(n));

      Rng rng(spec.seed ^ seed_salt);

      for (u64 i = 0; i < n; ++i) {
        const bool is_large = rng.Bernoulli(p_large);
        int large_axis = -1;
        if (is_large && anisotropic) {
          large_axis = static_cast<int>(rng.UniformU64(static_cast<u64>(Dim)));
        }

        BoxT b;
        for (int axis = 0; axis < Dim; ++axis) {
          double w_frac = 0.0;
          if (!is_large) {
            w_frac = rng.UniformDouble(ws_min, ws_max);
          } else {
            if (!anisotropic || axis == large_axis) {
              w_frac = rng.UniformDouble(wl_min, wl_max);
            } else {
              // Keep thin along other axes.
              w_frac = rng.UniformDouble(ws_min, ws_max);
            }
          }
          const double w = w_frac * L;

          const double max_lo = dom_hi - w;
          double lo = rng.UniformDouble(dom_lo, max_lo);
          double hi = lo + w;

          // Safety
          if (!(hi > lo)) {
            hi = std::nextafter(lo, dom_hi);
            if (!(hi > lo)) hi = lo + 1e-12;
          }

          const usize axis_idx = static_cast<usize>(axis);
          b.lo.v[axis_idx] = static_cast<T>(lo);
          b.hi.v[axis_idx] = static_cast<T>(hi);
        }

        out_rel->boxes[static_cast<usize>(i)] = b;
        out_rel->ids[static_cast<usize>(i)] = static_cast<Id>(i);
      }
      return true;
    };

    Relation<Dim, T> R, S;
    R.name = "R";
    S.name = "S";

    if (!gen_relation(spec.n_r, "r_", &R, 0x11223344ULL)) return false;
    if (!gen_relation(spec.n_s, "s_", &S, 0x55667788ULL)) return false;

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
        detail::SetErr(err, "HeteroSizesGenerator produced invalid dataset: " + tmp);
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
      notes << "p_large_r=" << detail::GetDoubleSide(spec.params, "r_", "p_large", 0.1)
            << ", w_small=[" << detail::GetDoubleSide(spec.params, "r_", "w_small_min", 0.002)
            << "," << detail::GetDoubleSide(spec.params, "r_", "w_small_max", 0.01) << "]"
            << ", w_large=[" << detail::GetDoubleSide(spec.params, "r_", "w_large_min", 0.1)
            << "," << detail::GetDoubleSide(spec.params, "r_", "w_large_max", 0.3) << "]"
            << ", anisotropic=" << (anisotropic ? 1 : 0);
      report->notes = notes.str();
    }

    *out_ds = std::move(ds);
    return true;
  }
};

}  // namespace synthetic
}  // namespace sjs
