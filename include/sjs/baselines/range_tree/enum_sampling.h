#pragma once
// sjs/baselines/range_tree/enum_sampling.h
//
// Baseline: "range_tree" — EnumSampling variant (Framework I).
//
// Materialize the full join J using the Chapter-4 RangeTree REPORT primitive
// (sweep + dynamic 2(d-1)-dim range tree), then sample t i.i.d. uniform pairs
// with replacement by uniform indexing.
//
// Intended for small datasets / oracle verification.

#include "sjs/baselines/range_tree/sampling.h"

#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sjs {
namespace baselines {
namespace range_tree {

namespace detail {

class VectorJoinEnumerator final : public baselines::IJoinEnumerator {
 public:
  explicit VectorJoinEnumerator(const std::vector<PairId>* pairs) : pairs_(pairs) {
    SJS_DASSERT(pairs_ != nullptr);
    Reset();
  }

  void Reset() override {
    i_ = 0;
    stats_ = join::JoinStats{};
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
class RangeTreeEnumSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  using DatasetT = Dataset<Dim, T>;

  Method method() const noexcept override { return Method::RangeTree; }
  Variant variant() const noexcept override { return Variant::EnumSampling; }
  std::string_view Name() const noexcept override { return "range_tree_enum+sampling"; }

  void Reset() override {
    built_ = false;
    ctx_.Reset();
    cached_ = false;
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
      if (err) *err = "RangeTreeEnumSamplingBaseline::Count: call Build() first";
      return false;
    }

    if (!cached_) {
      auto _ = phases ? phases->Scoped("materialize_join") : PhaseRecorder::ScopedPhase(nullptr, "");
      pairs_.clear();
      const u64 cap = cfg.run.enum_cap;

      detail::RangeTreeReportJoinEnumeratorND<Dim, T> it(&ctx_);
      PairId p;
      while (it.Next(&p)) {
        pairs_.push_back(p);
        if (cap > 0 && static_cast<u64>(pairs_.size()) > cap) {
          pairs_.clear();
          cached_ = false;
          if (err) *err = "RangeTreeEnumSamplingBaseline: join size exceeds enum_cap; refusing to materialize";
          return false;
        }
      }
      cached_ = true;
    }

    if (out) *out = MakeExactCount(static_cast<u64>(pairs_.size()));
    return true;
  }

  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    (void)phases;

    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: null rng/out";
      return false;
    }
    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "RangeTreeEnumSamplingBaseline::Sample: run.t too large (must fit in u32)";
      return false;
    }
    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    if (t == 0) return true;

    if (!cached_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }

    const u64 W = static_cast<u64>(pairs_.size());
    if (W == 0) return true;

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
      if (err) *err = "RangeTreeEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    if (cached_) return std::make_unique<detail::VectorJoinEnumerator>(&pairs_);
    return std::make_unique<detail::RangeTreeReportJoinEnumeratorND<Dim, T>>(&ctx_);
  }

 private:
  bool built_{false};
  detail::RangeTreeContext<Dim, T> ctx_;

  bool cached_{false};
  std::vector<PairId> pairs_;
};

}  // namespace range_tree
}  // namespace baselines
}  // namespace sjs
