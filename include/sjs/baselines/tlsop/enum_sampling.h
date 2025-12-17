#pragma once
// sjs/baselines/tlsop/enum_sampling.h
//
// TLSOP baseline (Two-layer SOP) (Variant::EnumSampling).
//
// This variant exposes the deterministic TLSOP join enumerator (JoinStream)
// defined in tlsop/sampling.h and performs uniform sampling WITH replacement
// via generic two-pass rank sampling over that enumerator.
//
// In the experiment harness, EnumSampling variants are typically executed via
// baselines::runners::RunEnumSamplingOnce(), which directly calls Enumerate()
// and performs rank sampling itself. We still provide Count() and Sample()
// here for completeness.

#include "sjs/baselines/tlsop/sampling.h"  // TLSOPIndex + TLSOPJoinEnumerator

#include "sjs/sampling/rank_sampling.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace tlsop {

template <int Dim, class T = Scalar>
class TLSOPEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim == 2, "TLSOPEnumSamplingBaseline is currently implemented for Dim==2 only");
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::TLSOP; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "tlsop_enum+sampling"; }

  void Reset() override {
    index_.Reset();
    built_ = false;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");
    built_ = false;
    if (!index_.Build(ds, cfg, phases, err)) return false;
    built_ = true;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)cfg;
    (void)rng;
    if (!built_ || !index_.Built()) {
      if (err) *err = "TLSOPEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "TLSOPEnumSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");

    detail::TLSOPJoinEnumerator<Dim, T> stream(&index_);
    u64 W = 0;
    PairId tmp;
    while (stream.Next(&tmp)) ++W;
    *out = MakeExactCount(W);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !index_.Built()) {
      if (err) *err = "TLSOPEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "TLSOPEnumSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "TLSOPEnumSamplingBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    auto scoped = phases ? phases->Scoped("rank_sampling") : PhaseRecorder::ScopedPhase(nullptr, "");

    detail::TLSOPJoinEnumerator<Dim, T> stream(&index_);
    std::vector<PairId> samples;
    sampling::RankSamplingInfo info;

    if (!sampling::RankSampleWithReplacement<detail::TLSOPJoinEnumerator<Dim, T>, PairId>(
            &stream, t, rng, &samples, &info, err)) {
      return false;
    }

    out->pairs = std::move(samples);
    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !index_.Built()) {
      if (err) *err = "TLSOPEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::TLSOPJoinEnumerator<Dim, T>>(&index_);
  }

 private:
  detail::TLSOPIndex<Dim, T> index_;
  bool built_ = false;
};

}  // namespace tlsop
}  // namespace baselines
}  // namespace sjs
