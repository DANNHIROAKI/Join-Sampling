#pragma once
// sjs/data/synthetic/uniform.h
//
// Uniform random boxes within a bounded domain.
//
// Parameters (spec.params):
//  - "w_min" / "w_max": width fraction range in each axis (default 0.005..0.02).
//  - Side-specific overrides:
//      "r_w_min", "r_w_max", "s_w_min", "s_w_max".
//  - "same_size_all_dims" (bool, default false): if true, sample one width and reuse for all dims.
//  - "shuffle_r" / "shuffle_s" (bool, default false): permute order in each relation.

#include "sjs/data/synthetic/generator.h"

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

namespace sjs {
namespace synthetic {

template <int Dim, class T = Scalar>
class UniformGenerator final : public ISyntheticGenerator<Dim, T> {
 public:
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;

  std::string_view Name() const noexcept override { return "uniform"; }

  bool Generate(const DatasetSpec& spec, DatasetT* out_ds, Report* report, std::string* err) override {
    if (!out_ds) {
      detail::SetErr(err, "UniformGenerator: out_ds is null");
      return false;
    }
    if (spec.n_r == 0 || spec.n_s == 0) {
      detail::SetErr(err, "UniformGenerator: n_r and n_s must be > 0");
      return false;
    }
    if (!detail::CheckIdFits(spec.n_r, err, "UniformGenerator(R)")) return false;
    if (!detail::CheckIdFits(spec.n_s, err, "UniformGenerator(S)")) return false;

    const double dom_lo = spec.domain_lo;
    const double dom_hi = spec.domain_hi;
    if (!(dom_hi > dom_lo)) {
      detail::SetErr(err, "UniformGenerator: domain_hi must be > domain_lo");
      return false;
    }
    const double L = dom_hi - dom_lo;

    const bool same_size = detail::GetBool(spec.params, "same_size_all_dims", false);
    const bool shuffle_r = detail::GetBool(spec.params, "shuffle_r", false);
    const bool shuffle_s = detail::GetBool(spec.params, "shuffle_s", false);

    auto gen_relation = [&](u64 n, std::string_view prefix, Relation<Dim, T>* out_rel, u64 salt) -> bool {
      const double w_min_frac = detail::GetDoubleSide(spec.params, prefix, "w_min", 0.005);
      const double w_max_frac = detail::GetDoubleSide(spec.params, prefix, "w_max", 0.02);
      if (!(w_min_frac > 0.0 && w_max_frac >= w_min_frac && w_max_frac < 1.0)) {
        detail::SetErr(err, "UniformGenerator: invalid width fraction range (need 0<w_min<=w_max<1)");
        return false;
      }

      out_rel->boxes.resize(static_cast<usize>(n));
      out_rel->ids.resize(static_cast<usize>(n));

      Rng rng(DeriveSeed(spec.seed, salt));

      for (u64 i = 0; i < n; ++i) {
        BoxT b;
        double shared_w = 0.0;
        if (same_size) {
          const double wf = rng.UniformDouble(w_min_frac, w_max_frac);
          shared_w = wf * L;
        }

        for (int axis = 0; axis < Dim; ++axis) {
          const double w = same_size ? shared_w : (rng.UniformDouble(w_min_frac, w_max_frac) * L);
          const double max_lo = dom_hi - w;
          double lo = rng.UniformDouble(dom_lo, max_lo);
          double hi = lo + w;

          // Numerical safety: enforce hi>lo.
          if (!(hi > lo)) {
            hi = std::nextafter(lo, dom_hi);
            if (!(hi > lo)) hi = lo + 1e-12;
          }

          const usize a = static_cast<usize>(axis);
          b.lo.v[a] = static_cast<T>(lo);
          b.hi.v[a] = static_cast<T>(hi);
        }

        out_rel->boxes[static_cast<usize>(i)] = b;
        out_rel->ids[static_cast<usize>(i)] = static_cast<Id>(i);
      }
      return true;
    };

    Relation<Dim, T> R, S;
    R.name = "R";
    S.name = "S";

    if (!gen_relation(spec.n_r, "r_", &R, 0xA1B2C3D4ULL)) return false;
    if (!gen_relation(spec.n_s, "s_", &S, 0xC3D4E5F6ULL)) return false;

    if (shuffle_r) {
      Rng rng_shuf(DeriveSeed(spec.seed, 0x11111111ULL));
      detail::ShuffleRelation(&R, &rng_shuf);
    }
    if (shuffle_s) {
      Rng rng_shuf(DeriveSeed(spec.seed, 0x22222222ULL));
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
        detail::SetErr(err, "UniformGenerator produced invalid dataset: " + tmp);
        return false;
      }
    }

    if (report) {
      report->generator = std::string(Name());
      report->dataset_name = spec.name;
      report->n_r = spec.n_r;
      report->n_s = spec.n_s;
      report->has_exact_k = false;
      report->k_target = 0;
      report->k_achieved = 0;
      report->alpha_target = spec.alpha;
      report->alpha_achieved = std::numeric_limits<double>::quiet_NaN();

      std::ostringstream notes;
      notes << "w_min_r=" << detail::GetDoubleSide(spec.params, "r_", "w_min", 0.005)
            << ", w_max_r=" << detail::GetDoubleSide(spec.params, "r_", "w_max", 0.02)
            << ", w_min_s=" << detail::GetDoubleSide(spec.params, "s_", "w_min", 0.005)
            << ", w_max_s=" << detail::GetDoubleSide(spec.params, "s_", "w_max", 0.02)
            << ", same_size_all_dims=" << (same_size ? 1 : 0);
      report->notes = notes.str();
    }

    *out_ds = std::move(ds);
    return true;
  }
};

}  // namespace synthetic
}  // namespace sjs
