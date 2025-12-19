#pragma once
// sjs/data/synthetic/uniform.h
//
// Uniform random rectangles within a bounded domain.
//
// This generator is NOT meant to control join size; it provides a simple,
// commonly-used baseline distribution for stress tests.
//
// Parameters (spec.params):
//  - "w_min" / "w_max": width fraction range in each axis (default 0.005..0.02).
//  - Side-specific overrides:
//      "r_w_min", "r_w_max", "s_w_min", "s_w_max".
//  - "same_size_all_dims" (bool, default false): if true, sample one width and reuse for all dims.
//  - "shuffle_r" / "shuffle_s" (bool, default false): permute order in each relation.

#include "sjs/data/synthetic/generator.h"

#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace sjs {
namespace synthetic {

template <int Dim, class T = Scalar>
class UniformGenerator final : public ISyntheticGenerator<Dim, T> {
 public:
  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;

  std::string_view Name() const noexcept override { return "uniform"; }

  bool Generate(const DatasetSpec& spec,
                DatasetT* out_ds,
                Report* report,
                std::string* err) override {
    if (!out_ds) {
      detail::SetErr(err, "UniformGenerator: out_ds is null");
      return false;
    }
    if (spec.n_r == 0 || spec.n_s == 0) {
      detail::SetErr(err, "UniformGenerator: n_r and n_s must be > 0");
      return false;
    }
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

    auto gen_relation = [&](u64 n, std::string_view prefix, Relation<Dim, T>* out_rel, u64 seed_salt) -> bool {
      const double w_min_frac = detail::GetDoubleSide(spec.params, prefix, "w_min", 0.005);
      const double w_max_frac = detail::GetDoubleSide(spec.params, prefix, "w_max", 0.02);
      if (!(w_min_frac > 0.0 && w_max_frac >= w_min_frac && w_max_frac < 1.0)) {
        detail::SetErr(err, "UniformGenerator: invalid width fraction range");
        return false;
      }

      out_rel->boxes.resize(static_cast<usize>(n));
      out_rel->ids.resize(static_cast<usize>(n));

      Rng rng(spec.seed ^ seed_salt);

      for (u64 i = 0; i < n; ++i) {
        BoxT b;
        double shared_w = 0.0;
        if (same_size) {
          const double wf = rng.UniformDouble(w_min_frac, w_max_frac);
          shared_w = wf * L;
        }
        for (int axis = 0; axis < Dim; ++axis) {
          const double w = same_size ? shared_w : rng.UniformDouble(w_min_frac, w_max_frac) * L;
          const double max_lo = dom_hi - w;
          const double lo = rng.UniformDouble(dom_lo, max_lo);
          const double hi = lo + w;
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

    if (!gen_relation(spec.n_r, "r_", &R, 0xA1B2C3D4ULL)) return false;
    if (!gen_relation(spec.n_s, "s_", &S, 0xC3D4E5F6ULL)) return false;

    if (shuffle_r) {
      Rng rng_shuffle(spec.seed + 17);
      // Shuffle boxes and ids together via a permutation.
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
