#pragma once
// baselines/ours/enum_sampling.h
//
// "Our Method" — Enumerate+Sampling variant (SJS v3 Framework I).
//
// Framework I:
//   1) Enumerate *all* join pairs once and materialize them into memory.
//   2) Draw t i.i.d. uniform samples with replacement by indexing that array.
//
// Memory: Theta(|J|).

#include "baselines/ours/sampling.h"  // reuse 2D preprocessing + enumerator

#include "baselines/detail/vector_join_enumerator.h"

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace ours {


template <int Dim, class T = Scalar>
class OursEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "ours_enum+sampling"; }

  void Reset() override {
    ctx_.Reset();
    nd_ctx_.Reset();
    built_ = false;
    pairs_cached_ = false;
    pairs_.clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    Reset();

    if constexpr (Dim == 2) {
      if (!ctx_.Build(ds, phases, err)) return false;
    } else {
      if (!nd_ctx_.Build(ds, phases, err)) return false;
    }

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;

    if constexpr (Dim == 2) {
      if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
        if (err) *err = "OursEnumSamplingBaseline::Count: call Build() first";
        return false;
      }
    } else {
      if (!built_ || !nd_ctx_.built() || nd_ctx_.dataset() == nullptr) {
        if (err) *err = "OursEnumSamplingBaseline::Count: call Build() first";
        return false;
      }
    }

    if (!pairs_cached_) {
      auto scoped = phases ? phases->Scoped("phase1_enumerate_materialize")
                           : PhaseRecorder::ScopedPhase(nullptr, "");

      pairs_.clear();
      const u64 cap = cfg.run.enum_cap;
      if (cap > 0) {
        const u64 reserve_n64 =
            (cap > static_cast<u64>(std::numeric_limits<usize>::max()))
                ? static_cast<u64>(std::numeric_limits<usize>::max())
                : cap;
        try {
          pairs_.reserve(static_cast<usize>(reserve_n64));
        } catch (...) {
        }
      }

      PairId p;
      if constexpr (Dim == 2) {
        detail::OursReportJoinEnumerator2D<Dim, T> it(&ctx_);
        while (it.Next(&p)) {
          pairs_.push_back(p);
          if (cap > 0 && static_cast<u64>(pairs_.size()) > cap) {
            pairs_.clear();
            pairs_cached_ = false;
            if (err) *err = "OursEnumSamplingBaseline: join size exceeds enum_cap; refusing to materialize";
            return false;
          }
        }
      } else {
        detail::OursReportJoinEnumeratorND<Dim, T> it(&nd_ctx_);
        while (it.Next(&p)) {
          pairs_.push_back(p);
          if (cap > 0 && static_cast<u64>(pairs_.size()) > cap) {
            pairs_.clear();
            pairs_cached_ = false;
            if (err) *err = "OursEnumSamplingBaseline: join size exceeds enum_cap; refusing to materialize";
            return false;
          }
        }
      }

      pairs_cached_ = true;
    }

    if (out) *out = MakeExactCount(static_cast<u64>(pairs_.size()));
    return true;
  }

  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    if constexpr (Dim == 2) {
      if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
        if (err) *err = "OursEnumSamplingBaseline::Sample: call Build() first";
        return false;
      }
    } else {
      if (!built_ || !nd_ctx_.built() || nd_ctx_.dataset() == nullptr) {
        if (err) *err = "OursEnumSamplingBaseline::Sample: call Build() first";
        return false;
      }
    }
    if (!rng || !out) {
      if (err) *err = "OursEnumSamplingBaseline::Sample: null rng/out";
      return false;
    }

    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "OursEnumSamplingBaseline::Sample: run.t too large (must fit in u32)";
      return false;
    }

    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    if (t == 0) return true;

    if (!pairs_cached_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    const u64 W = static_cast<u64>(pairs_.size());
    if (W == 0) return true;

    {
      auto scoped = phases ? phases->Scoped("phase2_resample") : PhaseRecorder::ScopedPhase(nullptr, "");
      out->pairs.resize(t);
      for (u32 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(W);
        out->pairs[i] = pairs_[static_cast<usize>(idx)];
      }
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;

    if constexpr (Dim == 2) {
      if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
        if (err) *err = "OursEnumSamplingBaseline::Enumerate: call Build() first";
        return nullptr;
      }
    } else {
      if (!built_ || !nd_ctx_.built() || nd_ctx_.dataset() == nullptr) {
        if (err) *err = "OursEnumSamplingBaseline::Enumerate: call Build() first";
        return nullptr;
      }
    }

    if (pairs_cached_) {
      return std::make_unique<baselines::detail::VectorJoinEnumerator>(&pairs_);
    }
    if constexpr (Dim == 2) {
      return std::make_unique<detail::OursReportJoinEnumerator2D<Dim, T>>(&ctx_);
    } else {
      return std::make_unique<detail::OursReportJoinEnumeratorND<Dim, T>>(&nd_ctx_);
    }
  }

 private:
  bool built_{false};

  bool pairs_cached_{false};
  std::vector<PairId> pairs_;

  detail::Ours2DContext<Dim, T> ctx_;
  detail::OursNDContext<Dim, T> nd_ctx_;
};


extern template class OursEnumSamplingBaseline<2, sjs::Scalar>;
extern template class OursEnumSamplingBaseline<3, sjs::Scalar>;
extern template class OursEnumSamplingBaseline<4, sjs::Scalar>;
extern template class OursEnumSamplingBaseline<5, sjs::Scalar>;

}  // namespace ours
}  // namespace baselines
}  // namespace sjs
