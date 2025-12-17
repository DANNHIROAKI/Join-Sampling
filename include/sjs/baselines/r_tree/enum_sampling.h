#pragma once
// sjs/baselines/r_tree/enum_sampling.h
//
// Plane Sweep + Dynamic R-Tree baseline (Variant::EnumSampling).
//
// This variant exposes a deterministic enumerator (sweep + R-tree) and
// implements Sample() by two-pass rank sampling over that enumerator.
//
// In the experimental harness, EnumSampling variants are typically executed via
// baselines::runners::RunEnumSamplingOnce(), which directly calls Enumerate()
// and performs rank sampling itself. Still, we provide Sample() here for
// completeness.

#include "sjs/baselines/r_tree/sampling.h"
#include "sjs/sampling/rank_sampling.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace r_tree {

template <int Dim, class T>
class RTreeEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "RTreeEnumSamplingBaseline requires Dim >= 2");
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::RTree; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "rtree_enum+sampling"; }

  void Reset() override {
    ds_ = nullptr;
    built_ = false;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    auto scoped = phases ? phases->Scoped("build") : PhaseRecorder::ScopedPhase(nullptr, "");

    ds_ = &ds;
    built_ = true;
    // Empty relations are fine.
    (void)err;
    return true;
  }

  bool Count(const Config& cfg,
             Rng* rng,
             CountResult* out,
             PhaseRecorder* phases,
             std::string* err) override {
    (void)rng;
    if (!built_ || !ds_) {
      if (err) *err = "RTreeEnumSamplingBaseline::Count: call Build() first";
      return false;
    }
    if (!out) {
      if (err) *err = "RTreeEnumSamplingBaseline::Count: out is null";
      return false;
    }

    // Efficient exact count using the same Phase-1 scheme as the Sampling variant.
    auto scoped = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");
    RTreeSamplingBaseline<Dim, T> counter;
    if (!counter.Build(*ds_, /*cfg=*/cfg, phases, err)) return false;
    return counter.Count(cfg, /*rng=*/nullptr, out, phases, err);
  }

  bool Sample(const Config& cfg,
              Rng* rng,
              SampleSet* out,
              PhaseRecorder* phases,
              std::string* err) override {
    if (!built_ || !ds_) {
      if (err) *err = "RTreeEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng) {
      if (err) *err = "RTreeEnumSamplingBaseline::Sample: rng is null";
      return false;
    }
    if (!out) {
      if (err) *err = "RTreeEnumSamplingBaseline::Sample: out is null";
      return false;
    }

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    out->weights.clear();

    const u64 t = cfg.run.t;
    if (t == 0) return true;

    auto scoped = phases ? phases->Scoped("rank_sampling") : PhaseRecorder::ScopedPhase(nullptr, "");

    auto stream = Enumerate(cfg, phases, err);
    if (!stream) return false;

    std::vector<PairId> samples;
    sampling::RankSamplingInfo info;
    if (!sampling::RankSampleWithReplacement<IJoinEnumerator, PairId>(stream.get(), t, rng, &samples, &info, err)) {
      return false;
    }

    out->pairs = std::move(samples);
    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg,
                                             PhaseRecorder* phases,
                                             std::string* err) override {
    auto scoped = phases ? phases->Scoped("enumerate_prepare") : PhaseRecorder::ScopedPhase(nullptr, "");
    if (!built_ || !ds_) {
      if (err) *err = "RTreeEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    // Parse the same R-tree options as the sampling variant so that
    // enumeration is consistent with the chosen configuration.
    typename detail::DynamicRTree<Dim - 1, T>::Options opt;
    {
      u32 m = 32;
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_max_children", m);
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_M", m);
      m = detail::ExtraU32Or(cfg.run.extra, "rtree_m", m);
      if (m < 2) m = 2;
      opt.max_children = m;

      u32 minc = 0;
      minc = detail::ExtraU32Or(cfg.run.extra, "rtree_min_children", minc);
      if (minc > 0) opt.min_children = minc;
      opt.ignore_duplicate_insert = true;
    }

    return std::make_unique<detail::RTreeJoinEnumerator<Dim, T>>(&ds_->R, &ds_->S, opt, /*axis=*/0,
                                                                 join::SideTieBreak::RBeforeS);
  }

 private:
  const DatasetT* ds_ = nullptr;
  bool built_ = false;
};

}  // namespace r_tree
}  // namespace baselines
}  // namespace sjs
