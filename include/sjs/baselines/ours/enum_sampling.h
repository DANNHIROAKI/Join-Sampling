#pragma once
// sjs/baselines/ours/enum_sampling.h
//
// "Our Method" — Enumerate+Sampling variant.
//
// This variant is the baseline described in "Our Method.md":
//   1) One sweep enumerates *all* join pairs and materializes them into memory.
//   2) Draw t i.i.d. uniform samples with replacement by indexing that array.
//
// Memory: Theta(|J|).

#include "sjs/baselines/ours/sampling.h"  // reuse 2D preprocessing + indices

#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
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
    built_ = false;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    Reset();

    if (!ctx_.Build(ds, phases, err)) return false;

    built_ = true;
    return true;
  }

  // Exact join cardinality using the same sweep+index logic (count only).
  bool Count(const Config& cfg, Rng* rng, CountResult* out, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)rng;

    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursEnumSamplingBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("count") : PhaseRecorder::ScopedPhase(nullptr, "");

    ctx_.ResetActive();

    u64 W = 0;

    auto& ar = ctx_.active_r();
    auto& as = ctx_.active_s();

    const auto& events = ctx_.events();
    const auto& sid_of_pos = ctx_.start_id_of_event();

    const auto& ylo_r = ctx_.ylo_rank_r();
    const auto& yhi_r = ctx_.yhi_lb_rank_r();
    const auto& ylo_s = ctx_.ylo_rank_s();
    const auto& yhi_s = ctx_.yhi_lb_rank_s();

    for (usize pos = 0; pos < events.size(); ++pos) {
      const auto& ev = events[pos];
      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        if (ev.side == join::Side::R) ar.Erase(handle);
        else as.Erase(handle);
        continue;
      }

      const i32 sid_i32 = sid_of_pos[pos];
      (void)sid_i32;  // only for debug sanity
      SJS_DASSERT(sid_i32 >= 0);

      const bool q_is_r = (ev.side == join::Side::R);
      const u32 q_ylo = q_is_r ? ylo_r[ev.index] : ylo_s[ev.index];
      const u32 q_yhi = q_is_r ? yhi_r[ev.index] : yhi_s[ev.index];

      const detail::ActiveIndex2D& other = q_is_r ? as : ar;
      W += other.CountA(q_ylo) + other.CountB(q_ylo, q_yhi);

      if (q_is_r) ar.Insert(handle, q_ylo, q_yhi);
      else as.Insert(handle, q_ylo, q_yhi);
    }

    ctx_.ResetActive();

    if (out) *out = MakeExactCount(W);
    return true;
  }

  // Materialize all pairs then sample by uniform indexing.
  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursEnumSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "OursEnumSamplingBaseline::Sample: null rng/out";
      return false;
    }

    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;

    if (t == 0) return true;

    std::vector<PairId> all_pairs;
    {
      auto scoped = phases ? phases->Scoped("enumerate_all") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!EnumerateAllPairs(&all_pairs, err)) return false;
    }

    const u64 W = static_cast<u64>(all_pairs.size());
    if (W == 0) return true;

    out->pairs.resize(t);
    for (u32 i = 0; i < t; ++i) {
      const u64 idx = rng->UniformU64(W);
      out->pairs[i] = all_pairs[static_cast<usize>(idx)];
    }
    return true;
  }

  // Provide a deterministic enumerator; for consistency with the rest of the
  // codebase we return the generic plane-sweep stream wrapper.
  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursEnumSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }

    join::PlaneSweepOptions opt;
    opt.axis = 0;
    opt.side_order = join::SideTieBreak::RBeforeS;
    opt.skip_axis_check = true;

    const auto* ds = ctx_.dataset();
    return std::make_unique<detail::PlaneSweepEnumeratorWrapper<Dim, T>>(ds->R, ds->S, opt);
  }

 private:
  bool EnumerateAllPairs(std::vector<PairId>* out_pairs, std::string* err) {
    SJS_DASSERT(out_pairs != nullptr);
    out_pairs->clear();

    const auto* ds = ctx_.dataset();
    if (!ds) {
      if (err) *err = "OursEnumSamplingBaseline::EnumerateAllPairs: missing dataset";
      return false;
    }

    ctx_.ResetActive();

    auto& ar = ctx_.active_r();
    auto& as = ctx_.active_s();

    const auto& events = ctx_.events();

    const auto& ylo_r = ctx_.ylo_rank_r();
    const auto& yhi_r = ctx_.yhi_lb_rank_r();
    const auto& ylo_s = ctx_.ylo_rank_s();
    const auto& yhi_s = ctx_.yhi_lb_rank_s();

    std::vector<u32> tmp;
    tmp.reserve(1024);

    for (usize pos = 0; pos < events.size(); ++pos) {
      const auto& ev = events[pos];
      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        if (ev.side == join::Side::R) ar.Erase(handle);
        else as.Erase(handle);
        continue;
      }

      // START
      const bool q_is_r = (ev.side == join::Side::R);
      const u32 q_ylo = q_is_r ? ylo_r[ev.index] : ylo_s[ev.index];
      const u32 q_yhi = q_is_r ? yhi_r[ev.index] : yhi_s[ev.index];

      const detail::ActiveIndex2D& other = q_is_r ? as : ar;

      // Pattern A
      tmp.clear();
      other.ReportA(q_ylo, &tmp);
      for (u32 oh : tmp) {
        if (q_is_r) out_pairs->push_back(PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)});
        else out_pairs->push_back(PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)});
      }

      // Pattern B
      tmp.clear();
      other.ReportB(q_ylo, q_yhi, &tmp);
      for (u32 oh : tmp) {
        if (q_is_r) out_pairs->push_back(PairId{ds->R.GetId(ev.index), ds->S.GetId(oh)});
        else out_pairs->push_back(PairId{ds->R.GetId(oh), ds->S.GetId(ev.index)});
      }

      // Insert q
      if (q_is_r) ar.Insert(handle, q_ylo, q_yhi);
      else as.Insert(handle, q_ylo, q_yhi);
    }

    ctx_.ResetActive();
    return true;
  }

  bool built_{false};
  detail::Ours2DContext<Dim, T> ctx_;
};

}  // namespace ours
}  // namespace baselines
}  // namespace sjs
