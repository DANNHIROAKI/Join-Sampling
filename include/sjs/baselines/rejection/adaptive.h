#pragma once
// sjs/baselines/rejection/adaptive.h
//
// AGR-BoxJoin baseline (Variant::Adaptive).
//
// This follows the adaptive strategy described in "Amagata’25.md" (AGR-AS):
//   - Build() performs the same preprocessing as the Sampling variant (AGR-S).
//   - Count() enumerates join pairs until either:
//       * the stream ends  (|J| <= J*): keep all pairs in memory (AllPairs)
//         and return the exact join size;
//       * more than J* pairs are seen (|J| > J*): discard AllPairs and switch
//         to rejection sampling (AGR-S) for Sample(). Count() returns a pilot
//         estimate of |J| in this case.
//   - Sample() does:
//       * if small join: i.i.d. uniform sampling with replacement from AllPairs
//       * else: rejection sampling using the prebuilt alias tables
//
// This class is still useful even if you also use the generic adaptive runner
// in sjs/baselines/runners/adaptive_runner.h.

#include "sjs/baselines/rejection/sampling.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace rejection {

template <int Dim, class T = Scalar>
class RejectionAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "RejectionAdaptiveBaseline expects Dim>=2");
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Rejection; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "rejection_adaptive"; }

  void Reset() override {
    st_.Reset();
    mode_ = Mode::Unknown;
    all_pairs_.clear();
    all_pairs_.shrink_to_fit();
    exact_n_ = 0;
  }

  bool Build(const DatasetT& ds,
             const Config& cfg,
             PhaseRecorder* phases,
             std::string* err) override {
    Reset();
    // Build the shared preprocessing (cell map + alias).
    if (!st_.BuildIndex(ds, cfg, phases, err)) return false;
    mode_ = Mode::Unknown;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    if (!st_.built) {
      if (err) *err = "RejectionAdaptiveBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RejectionAdaptiveBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("adaptive_phase1") : PhaseRecorder::ScopedPhase(nullptr, "");

    // Effective J* bound.
    u64 J_star = cfg.run.j_star;
    if (cfg.run.enum_cap > 0) J_star = std::min(J_star, cfg.run.enum_cap);

    // Special-case: J_star==0 still performs a *minimal* enumeration to
    // distinguish empty join from non-empty join, preventing non-termination.
    const u64 limit = (J_star == 0) ? 1ULL : (J_star + 1ULL);

    all_pairs_.clear();
    exact_n_ = 0;

    // Enumerate up to limit pairs.
    auto stream = std::make_unique<detail::RejectionJoinEnumerator<Dim, T>>(&st_);
    PairId p;
    u64 seen = 0;
    bool hit_limit = false;

    while (stream->Next(&p)) {
      ++seen;
      if (seen <= J_star) {
        // Store up to J* pairs.
        all_pairs_.push_back(p);
      } else {
        // We have seen more than J*.
        hit_limit = true;
        break;
      }
      if (seen >= limit) {
        // For J_star==0, limit==1: stop after finding one pair.
        if (J_star == 0) {
          hit_limit = true;  // treat as non-empty and switch to sampling
        }
        break;
      }
    }

    if (!hit_limit) {
      // Fully enumerated (|J| <= J*).
      mode_ = Mode::Enumerate;
      exact_n_ = seen;
      // In this branch, all_pairs_ holds all pairs (since seen<=J*).
      *out = MakeExactCount(static_cast<u64>(exact_n_));
      return true;
    }

    // Switch to sampling branch.
    mode_ = Mode::Sampling;
    all_pairs_.clear();

    // Provide a pilot estimate of |J| (without consuming the caller's RNG
    // sequence, so Sample() remains reproducible given a seed).
    if (!rng) {
      // If rng is missing, fall back to an unknown estimate.
      *out = MakeEstimateCount(std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              std::numeric_limits<long double>::quiet_NaN(),
                              0);
      return true;
    }

    Rng tmp = *rng;
    return st_.EstimateCountByPilot(cfg, &tmp, out, phases, err);
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!st_.built) {
      if (err) *err = "RejectionAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "RejectionAdaptiveBaseline::Sample: null rng/out";
      return false;
    }

    auto scoped = phases ? phases->Scoped("adaptive_sample") : PhaseRecorder::ScopedPhase(nullptr, "");

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    if (mode_ == Mode::Unknown) {
      if (err) *err = "RejectionAdaptiveBaseline::Sample: Count() must be called before Sample()";
      return false;
    }

    if (mode_ == Mode::Enumerate) {
      // Sample directly from AllPairs.
      if (all_pairs_.empty()) return true;
      out->pairs.resize(static_cast<usize>(t));
      const u64 N = static_cast<u64>(all_pairs_.size());
      for (u64 i = 0; i < t; ++i) {
        out->pairs[static_cast<usize>(i)] = all_pairs_[static_cast<usize>(rng->UniformU64(N))];
      }
      return true;
    }

    // mode_ == Sampling: run rejection sampling (AGR-S).
    // Since we only switch here after observing >J* pairs (or at least 1 pair
    // when J*=0), the join is guaranteed non-empty.

    if (st_.mu_sum == 0) {
      // Defensive: if mu_sum==0, the join must be empty.
      return true;
    }

    // Guard against pathological acceptance rates.
    const detail::RejectionParams<Dim> p = detail::ReadRejectionParams<Dim, T>(*st_.ds, cfg);
    u64 max_props = p.max_proposals;
    if (max_props == 0) {
      long double m = static_cast<long double>(p.max_factor) * static_cast<long double>(t);
      const long double cap = static_cast<long double>(std::numeric_limits<u64>::max());
      if (m > cap) m = cap;
      max_props = static_cast<u64>(m);
      if (max_props < 1000ULL) max_props = 1000ULL;
    }

    out->pairs.reserve(static_cast<usize>(t));

    u64 props = 0;
    while (out->pairs.size() < static_cast<usize>(t)) {
      if (max_props > 0 && props >= max_props) {
        if (err) {
          *err = "RejectionAdaptiveBaseline::Sample: reached max_proposals; acceptance rate extremely small";
        }
        return false;
      }
      ++props;

      u32 ai = 0, bi = 0;
      if (!st_.DrawProposal(rng, &ai, &bi)) {
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
      if (err) *err = "RejectionAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RejectionJoinEnumerator<Dim, T>>(&st_);
  }

 private:
  enum class Mode : u8 { Unknown = 0, Enumerate = 1, Sampling = 2 };

  detail::RejectionState<Dim, T> st_;
  Mode mode_{Mode::Unknown};

  std::vector<PairId> all_pairs_;  // used when |J| <= J*
  u64 exact_n_{0};
};

}  // namespace rejection
}  // namespace baselines
}  // namespace sjs
