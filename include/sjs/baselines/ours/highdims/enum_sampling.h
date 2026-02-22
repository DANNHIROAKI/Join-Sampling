#pragma once
// sjs/baselines/ours/highdims/enum_sampling.h
//
// Our method (HighDims) — Enumerate+Sampling variant (Framework I).
//
// Framework I:
//   1) Enumerate and materialize all join pairs once (REPORT-based, event-block decomposition).
//   2) Draw t i.i.d. uniform samples with replacement by indexing that vector.
//
// Memory: Theta(|J|). Intended for small joins / oracle-like runs.

#include "sjs/baselines/ours/highdims/sampling.h"

#include "sjs/core/assert.h"
#include "sjs/core/rng.h"

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace ours {
namespace highdims {

namespace detail {

// Simple enumerator over a cached vector of pairs.
class VectorJoinEnumerator final : public baselines::IJoinEnumerator {
 public:
  explicit VectorJoinEnumerator(const std::vector<PairId>* pairs) : pairs_(pairs) {
    SJS_DASSERT(pairs_ != nullptr);
    Reset();
  }

  void Reset() override {
    i_ = 0;
    stats_.Reset();
    if (pairs_) {
      stats_.output_pairs = static_cast<u64>(pairs_->size());
      stats_.candidate_checks = stats_.output_pairs;
    }
  }

  bool Next(PairId* out) override {
    if (!out || !pairs_) return false;
    if (i_ >= pairs_->size()) return false;
    *out = (*pairs_)[i_++];
    return true;
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  const std::vector<PairId>* pairs_{nullptr};
  usize i_{0};
  join::JoinStats stats_{};
};

}  // namespace detail

template <int Dim, class T = Scalar>
class OursHighDimsEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "ours_highdims_enum+sampling"; }

  void Reset() override {
    built_ = false;
    ctx_.Reset();
    pairs_cached_ = false;
    pairs_.clear();
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    Reset();
    if (!ctx_.Build(ds, phases, err)) return false;
    built_ = true;
    return true;
  }

  bool Count(const Config& cfg, Rng* rng, CountResult* out, PhaseRecorder* phases, std::string* err) override {
    (void)rng;

    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursHighDimsEnumSamplingBaseline::Count: call Build() first";
      return false;
    }

    if (!pairs_cached_) {
      auto scoped = phases ? phases->Scoped("phase1_enumerate_materialize")
                           : PhaseRecorder::ScopedPhase(nullptr, "");

      pairs_.clear();
      const u64 cap = cfg.run.enum_cap;

      detail::OursReportJoinEnumeratorND<Dim, T> it(&ctx_);
      PairId p;
      while (it.Next(&p)) {
        pairs_.push_back(p);
        if (cap > 0 && static_cast<u64>(pairs_.size()) > cap) {
          pairs_.clear();
          pairs_cached_ = false;
          if (err) *err = "OursHighDimsEnumSamplingBaseline: join size exceeds enum_cap; refusing to materialize";
          return false;
        }
      }

      pairs_cached_ = true;
    }

    if (out) *out = MakeExactCount(static_cast<u64>(pairs_.size()));
    return true;
  }

  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursHighDimsEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "OursHighDimsEnumSamplingBaseline::Sample: null rng/out";
      return false;
    }

    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "OursHighDimsEnumSamplingBaseline::Sample: run.t too large (must fit in u32)";
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

    auto scoped = phases ? phases->Scoped("phase2_resample") : PhaseRecorder::ScopedPhase(nullptr, "");
    out->pairs.resize(static_cast<usize>(t));

    for (u32 i = 0; i < t; ++i) {
      const u64 idx = rng->UniformU64(W);
      out->pairs[static_cast<usize>(i)] = pairs_[static_cast<usize>(idx)];
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)phases;

    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursHighDimsEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    if (pairs_cached_) {
      return std::make_unique<detail::VectorJoinEnumerator>(&pairs_);
    }
    return std::make_unique<detail::OursReportJoinEnumeratorND<Dim, T>>(&ctx_);
  }

 private:
  bool built_{false};

  // Reuse the same context and report-based decomposition.
  detail::OursHighDimsContext<Dim, T> ctx_;

  bool pairs_cached_{false};
  std::vector<PairId> pairs_;
};

}  // namespace highdims
}  // namespace ours
}  // namespace baselines
}  // namespace sjs
