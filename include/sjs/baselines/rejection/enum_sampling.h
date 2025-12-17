#pragma once
// sjs/baselines/rejection/enum_sampling.h
//
// AGR-BoxJoin baseline (Variant::EnumSampling).
//
// This variant exposes a deterministic join enumerator (grid candidate scan)
// and performs uniform sampling WITH replacement using two-pass rank sampling
// over that stream.
//
// See:
//   - sjs/baselines/rejection/sampling.h (shared preprocessing + enumerator)
//   - sjs/sampling/rank_sampling.h (generic two-pass rank sampling)

#include "sjs/baselines/rejection/sampling.h"  // provides RejectionState + enumerator
#include "sjs/sampling/rank_sampling.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace rejection {

template <int Dim, class T = Scalar>
class RejectionEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "RejectionEnumSamplingBaseline expects Dim>=2");
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Rejection; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "rejection_enum+sampling"; }

  void Reset() override { st_.Reset(); }

  bool Build(const DatasetT& ds,
             const Config& cfg,
             PhaseRecorder* phases,
             std::string* err) override {
    return st_.BuildIndex(ds, cfg, phases, err);
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;  // deterministic
    if (!st_.built) {
      if (err) *err = "RejectionEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RejectionEnumSamplingBaseline::Count: out is null";
      return false;
    }

    auto scoped = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");

    // Exact count by full enumeration.
    auto stream = Enumerate(cfg, nullptr, err);
    if (!stream) return false;

    u64 n = 0;
    PairId p;
    while (stream->Next(&p)) {
      ++n;
    }

    *out = MakeExactCount(n);
    return true;
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!st_.built) {
      if (err) *err = "RejectionEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "RejectionEnumSamplingBaseline::Sample: null rng/out";
      return false;
    }

    auto scoped = phases ? phases->Scoped("sample") : PhaseRecorder::ScopedPhase(nullptr, "");

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    std::string enum_err;
    auto stream = Enumerate(cfg, nullptr, &enum_err);
    if (!stream) {
      if (err) *err = "RejectionEnumSamplingBaseline::Sample: Enumerate failed: " + enum_err;
      return false;
    }

    sampling::RankSamplingInfo info;
    std::vector<PairId> samples;
    if (!sampling::RankSampleWithReplacement<IJoinEnumerator, PairId>(stream.get(), t, rng, &samples, &info, err)) {
      return false;
    }

    out->pairs = std::move(samples);
    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!st_.built) {
      if (err) *err = "RejectionEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::RejectionJoinEnumerator<Dim, T>>(&st_);
  }

 private:
  detail::RejectionState<Dim, T> st_;
};

}  // namespace rejection
}  // namespace baselines
}  // namespace sjs
