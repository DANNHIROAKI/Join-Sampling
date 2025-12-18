#pragma once
// sjs/baselines/tlsop/enum_sampling.h
//
// TLSOP baseline (Two-layer SOP) (Variant::EnumSampling).
//
// This header implements **TS23-TLSOP-EnumerateThenSample** as described in
// docs/Baseline/Tsitsigkos’23 Baseline.md §4.5:
//   1) BuildIndex (two-layer grid + multi-assignment + A/B/C/D)
//   2) One pass JOIN_STREAM(Index) to materialize AllPairs = J (no duplicates)
//   3) i.i.d. uniform sampling with replacement from AllPairs by array indexing
//
// IMPORTANT: This variant intentionally uses \Theta(|J|) memory and may OOM for
// very large joins. For memory-safe sampling, use Variant::Sampling (2-pass)
// or Variant::Adaptive.

#include "sjs/baselines/tlsop/sampling.h"  // TLSOPIndex + TLSOPJoinEnumerator

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
  // Kept for backward-compatibility with the harness naming.
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
    (void)rng;  // deterministic
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

    // 1) One pass materialization: AllPairs = JOIN_STREAM(Index)
    std::vector<PairId> all_pairs;
    {
      auto scoped = phases ? phases->Scoped("tlsop_enum_materialize") : PhaseRecorder::ScopedPhase(nullptr, "");
      detail::TLSOPJoinEnumerator<Dim, T> stream(&index_);
      PairId p;
      while (stream.Next(&p)) {
        all_pairs.push_back(p);
      }
    }

    const u64 N = static_cast<u64>(all_pairs.size());
    if (N == 0) return true;

    // 2) i.i.d. uniform sampling with replacement from the array.
    {
      auto scoped = phases ? phases->Scoped("tlsop_enum_sample_array") : PhaseRecorder::ScopedPhase(nullptr, "");
      out->pairs.resize(static_cast<usize>(t));
      for (u64 i = 0; i < t; ++i) {
        const u64 idx = rng->UniformU64(N);  // [0, N)
        out->pairs[static_cast<usize>(i)] = all_pairs[static_cast<usize>(idx)];
      }
    }

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
