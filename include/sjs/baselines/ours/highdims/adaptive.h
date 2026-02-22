#pragma once
// sjs/baselines/ours/highdims/adaptive.h
//
// Our method (HighDims) — Adaptive variant (Framework III).
//
// Implements the budgeted caching idea in SJS-HighDims Framework III:
//
// Pass 1 (sweep):
//   - compute exact w_i = COUNT(K_i)
//   - under budget B (cfg.run.budget; cfg.run.extra["budget"] override; legacy alias: cfg.run.j_star):
//       * full cache C_i for small blocks (w_i <= w_small)
//       * prefetch sample cache S_i by a global min-heap of valued slots
//         (Poisson tail score; performance-only)
//
// Plan + Fill:
//   - assign t output positions to events i with Pr(i)=w_i/W
//   - fill from full cache C_i, then from prefetched samples S_i
//   - record residual slots that still need sampling
//
// Pass 2 (optional):
//   - second sweep to sample residual slots only
//
// Correctness: distribution is identical to Framework II; caching only affects performance.
//
// Knobs:
//   - budget B: cfg.run.budget (or cfg.run.extra["budget"] override; legacy alias: cfg.run.j_star)
//   - small-block threshold w_small: cfg.run.w_small (or cfg.run.extra["w_small"] override)
//   - prefetch on/off: cfg.run.prefetch (or cfg.run.extra["prefetch"] override)

#include "sjs/baselines/ours/highdims/sampling.h"

#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/sampling/alias_table.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace baselines {
namespace ours {
namespace highdims {

namespace detail {

// --------------------------
// parse helpers
// --------------------------
inline bool ParseU64(std::string_view s, u64* out) {
  if (!out) return false;
  if (s.empty()) return false;
  try {
    std::size_t idx = 0;
    const unsigned long long v = std::stoull(std::string(s), &idx, 10);
    if (idx != s.size()) return false;
    *out = static_cast<u64>(v);
    return true;
  } catch (...) {
    return false;
  }
}

inline bool ParseBool(std::string_view s, bool* out) {
  if (!out) return false;
  if (s.empty()) return false;

  auto lower = [](unsigned char c) { return static_cast<char>(std::tolower(c)); };
  std::string t;
  t.reserve(s.size());
  for (char c : s) t.push_back(lower(static_cast<unsigned char>(c)));

  if (t == "1" || t == "true" || t == "yes" || t == "y" || t == "on") {
    *out = true;
    return true;
  }
  if (t == "0" || t == "false" || t == "no" || t == "n" || t == "off") {
    *out = false;
    return true;
  }
  return false;
}

inline u64 ExtraU64Or(const Config& cfg, std::string_view key, u64 def) {
  auto it = cfg.run.extra.find(std::string(key));
  if (it == cfg.run.extra.end()) return def;
  u64 v = def;
  if (!ParseU64(it->second, &v)) return def;
  return v;
}

inline bool ExtraBoolOr(const Config& cfg, std::string_view key, bool def) {
  auto it = cfg.run.extra.find(std::string(key));
  if (it == cfg.run.extra.end()) return def;
  bool v = def;
  if (!ParseBool(it->second, &v)) return def;
  return v;
}

// --------------------------
// prefetch scoring (Poisson tail; performance-only)
// --------------------------
inline long double EstimateTotalW(u32 total_events, u32 i, long double w_sofar) noexcept {
  if (total_events == 0 || i == 0) return w_sofar;
  const long double M = static_cast<long double>(total_events);
  const long double ii = static_cast<long double>(i);
  return (ii > 0.0L) ? (w_sofar * (M / ii)) : w_sofar;
}

inline double PoissonSurvival(double mu, u32 r) noexcept {
  if (r <= 0) return 1.0;
  if (!(mu > 0.0) || !std::isfinite(mu)) return 0.0;

  // Small mu: exact prefix sum.
  if (mu <= 50.0 && r <= 200U) {
    const double p0 = std::exp(-mu);
    double p = p0;
    double cdf = p;
    for (u32 k = 1; k < r; ++k) {
      p *= mu / static_cast<double>(k);
      cdf += p;
    }
    double sf = 1.0 - cdf;
    if (sf < 0.0) sf = 0.0;
    if (sf > 1.0) sf = 1.0;
    return sf;
  }

  // Normal approximation.
  const double sigma = std::sqrt(mu);
  if (!(sigma > 0.0) || !std::isfinite(sigma)) return 0.0;
  const double z = (static_cast<double>(r) - 0.5 - mu) / sigma;
  const double sf = 0.5 * std::erfc(z * 0.70710678118654752440);  // 1/sqrt(2)
  if (sf < 0.0) return 0.0;
  if (sf > 1.0) return 1.0;
  return sf;
}

inline double SlotScorePoisson(u64 w_i,
                               u32 total_events,
                               u32 i,
                               long double w_sofar,
                               u64 t,
                               u32 r) noexcept {
  if (r == 0) return 1.0;
  if (w_i == 0 || t == 0) return 0.0;

  const long double W_hat = EstimateTotalW(total_events, i, w_sofar);
  if (!(W_hat > 0.0L)) return 0.0;

  const long double p_hat = static_cast<long double>(w_i) / W_hat;
  const long double mu_hat = static_cast<long double>(t) * p_hat;
  return PoissonSurvival(static_cast<double>(mu_hat), r);
}

struct PrefetchSlot {
  double score = 0.0;
  u32 sid = 0;
};

struct PrefetchSlotGreater {
  bool operator()(const PrefetchSlot& a, const PrefetchSlot& b) const noexcept {
    // min-heap by score
    if (a.score != b.score) return a.score > b.score;
    return a.sid > b.sid;
  }
};

class PrefetchHeap {
 public:
  void Clear() { while (!h_.empty()) h_.pop(); }
  usize Size() const noexcept { return h_.size(); }
  bool Empty() const noexcept { return h_.empty(); }
  double MinScore() const noexcept {
    if (h_.empty()) return -std::numeric_limits<double>::infinity();
    return h_.top().score;
  }
  void Push(PrefetchSlot s) { h_.push(s); }
  PrefetchSlot PopMin() {
    SJS_DASSERT(!h_.empty());
    PrefetchSlot t = h_.top();
    h_.pop();
    return t;
  }

 private:
  std::priority_queue<PrefetchSlot, std::vector<PrefetchSlot>, PrefetchSlotGreater> h_;
};

}  // namespace detail

template <int Dim, class T = Scalar>
class OursHighDimsAdaptiveBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "OursHighDimsAdaptiveBaseline requires Dim >= 2");

  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;
  using BoxT = Box<Dim, T>;
  using MF = index::ModeFamily<Dim, T>;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::Adaptive; }
  std::string_view Name() const noexcept override { return "ours_highdims_adaptive"; }

  void Reset() override {
    built_ = false;
    ctx_.Reset();

    start_side_.clear();
    start_index_.clear();

    w_total_.clear();
    W_ = 0;
    weights_valid_ = false;

    cached_.clear();
    cache_off_.clear();
    cache_len_.clear();
    cache_partners_.clear();

    prefetch_keep_.clear();
    prefetch_partners_.clear();

    cache_valid_ = false;
    budget_B_ = 0;
    budget_used_ = 0;
    w_small_ = 0;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    Reset();
    if (!ctx_.Build(ds, phases, err)) return false;

    const usize E = ctx_.num_start_events();
    w_total_.assign(E, 0ULL);

    // Mapping sid -> (side,index)
    start_side_.assign(E, join::Side::R);
    start_index_.assign(E, 0U);
    {
      const auto& events = ctx_.events();
      const auto& sid_of_pos = ctx_.start_id_of_event();
      for (usize pos = 0; pos < events.size(); ++pos) {
        if (events[pos].kind != join::EventKind::Start) continue;
        const i32 sid_i32 = sid_of_pos[pos];
        if (sid_i32 < 0) continue;
        const usize sid = static_cast<usize>(sid_i32);
        if (sid >= E) continue;
        start_side_[sid] = events[pos].side;
        start_index_[sid] = static_cast<u32>(events[pos].index);
      }
    }

    cached_.assign(E, 0U);
    cache_off_.assign(E, 0ULL);
    cache_len_.assign(E, 0ULL);
    cache_partners_.clear();

    prefetch_keep_.assign(E, 0U);
    prefetch_partners_.assign(E, std::vector<u32>{});

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg, Rng* rng, CountResult* out, PhaseRecorder* phases, std::string* err) override {
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursHighDimsAdaptiveBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count_and_cache") : PhaseRecorder::ScopedPhase(nullptr, "");

    const usize E = w_total_.size();
    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    std::fill(cached_.begin(), cached_.end(), 0U);
    std::fill(cache_off_.begin(), cache_off_.end(), 0ULL);
    std::fill(cache_len_.begin(), cache_len_.end(), 0ULL);
    cache_partners_.clear();

    std::fill(prefetch_keep_.begin(), prefetch_keep_.end(), 0U);
    for (auto& v : prefetch_partners_) v.clear();

    // Budget & threshold
    budget_B_ = detail::ExtraU64Or(cfg, "budget", cfg.run.budget);
    // Backward compatibility: if a legacy config sets j_star but not budget.
    if (budget_B_ == cfg.run.budget && cfg.run.j_star != cfg.run.budget) {
      budget_B_ = cfg.run.j_star;
    }

    w_small_ = detail::ExtraU64Or(cfg, "w_small", cfg.run.w_small);

    const bool prefetch_on = detail::ExtraBoolOr(cfg, "prefetch", cfg.run.prefetch);
    const bool enable_prefetch = prefetch_on && (rng != nullptr) && (budget_B_ > 0) && (cfg.run.t > 0);
    const u64 base_prefetch_seed = enable_prefetch ? rng->NextU64() : 0ULL;

    detail::PrefetchHeap heap;
    u64 mem_full = 0;

    ctx_.ResetActive();
    u64 W = 0;
    long double W_sofar = 0.0L;

    auto& ar = ctx_.active_r();
    auto& as = ctx_.active_s();

    const auto& events = ctx_.events();
    const auto& sid_of_pos = ctx_.start_id_of_event();

    std::vector<u32> tmp;
    tmp.reserve(1024);

    for (usize pos = 0; pos < events.size(); ++pos) {
      const auto& ev = events[pos];
      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        ctx_.Delete(ev.side, handle);
        continue;
      }

      // START
      const i32 sid_i32 = sid_of_pos[pos];
      SJS_DASSERT(sid_i32 >= 0);
      const u32 sid = static_cast<u32>(sid_i32);
      const usize sid_u = static_cast<usize>(sid);
      if (sid_u >= E) {
        if (err) *err = "OursHighDimsAdaptiveBaseline::Count: sid out of range";
        ctx_.ResetActive();
        return false;
      }

      const bool q_is_r = (ev.side == join::Side::R);
      const BoxT& q = ctx_.BoxOf(ev);
      const MF& other = q_is_r ? as : ar;

      const u64 w = other.Count(q);
      w_total_[sid_u] = w;

      // Total join size W
      if (w > 0) {
        if (W > (std::numeric_limits<u64>::max() - w)) {
          if (err) *err = "OursHighDimsAdaptiveBaseline::Count: |J| overflowed u64";
          ctx_.ResetActive();
          return false;
        }
        W += w;
      }
      W_sofar += static_cast<long double>(w);

      // 1) Full cache for small blocks
      const bool can_full_cache =
          (w_small_ > 0 && w > 0 && w <= w_small_ && (mem_full + w) <= budget_B_);

      if (can_full_cache) {
        cached_[sid_u] = 1U;
        cache_off_[sid_u] = static_cast<u64>(cache_partners_.size());

        tmp.clear();
        other.Report(q, &tmp);
        cache_partners_.insert(cache_partners_.end(), tmp.begin(), tmp.end());

        const u64 off = cache_off_[sid_u];
        const u64 len = static_cast<u64>(cache_partners_.size()) - off;
        cache_len_[sid_u] = len;

        if (len != w) {
          if (err) *err = "OursHighDimsAdaptiveBaseline::Count: full cache REPORT size mismatch (bug)";
          ctx_.ResetActive();
          return false;
        }

        mem_full += w;

        // Shrink heap capacity if full-cache consumed budget.
        const u64 B_samp = (budget_B_ > mem_full) ? (budget_B_ - mem_full) : 0ULL;
        while (heap.Size() > static_cast<usize>(B_samp)) {
          const auto popped = heap.PopMin();
          const usize j = static_cast<usize>(popped.sid);
          if (j < prefetch_keep_.size() && prefetch_keep_[j] > 0) {
            --prefetch_keep_[j];
            if (!prefetch_partners_[j].empty()) prefetch_partners_[j].pop_back();
          }
        }
      } else if (enable_prefetch && w > 0) {
        // 2) Prefetch sample cache under remaining budget
        const u64 B_samp = (budget_B_ > mem_full) ? (budget_B_ - mem_full) : 0ULL;
        if (B_samp > 0) {
          // i is 1-based START index (sid is 0-based in event order).
          const u32 i = sid + 1U;

          while (true) {
            const u32 r = prefetch_keep_[sid_u] + 1U;
            const double score = detail::SlotScorePoisson(
                w, static_cast<u32>(E), i, W_sofar, cfg.run.t, r);

            if (heap.Size() < static_cast<usize>(B_samp)) {
              heap.Push(detail::PrefetchSlot{score, sid});
              ++prefetch_keep_[sid_u];
              continue;
            }

            if (score > heap.MinScore()) {
              heap.Push(detail::PrefetchSlot{score, sid});
              ++prefetch_keep_[sid_u];

              const auto popped = heap.PopMin();
              const usize j = static_cast<usize>(popped.sid);
              if (j < prefetch_keep_.size() && prefetch_keep_[j] > 0) {
                --prefetch_keep_[j];
                if (!prefetch_partners_[j].empty()) prefetch_partners_[j].pop_back();
              }
              continue;
            }

            // Scores are non-increasing in r -> stop.
            break;
          }

          const u32 s_keep = prefetch_keep_[sid_u];
          if (s_keep > 0) {
            // Generate a prefix of i.i.d samples for this event (with replacement).
            const u64 ev_seed = DeriveSeed(base_prefetch_seed, 0xA11C3EULL, static_cast<u64>(sid));
            Rng rng_samp(DeriveSeed(ev_seed, 1));

            auto& dst = prefetch_partners_[sid_u];
            dst.clear();
            dst.reserve(static_cast<usize>(s_keep));

            std::vector<u32> samp;
            samp.clear();
            other.Sample(q, static_cast<u64>(s_keep), &rng_samp, &samp);
            if (samp.size() != static_cast<usize>(s_keep)) {
              if (err) *err = "OursHighDimsAdaptiveBaseline::Count: prefetch SAMPLE failed";
              ctx_.ResetActive();
              return false;
            }
            dst.assign(samp.begin(), samp.end());
          }
        }
      }

      // Insert q after query/caching decision.
      ctx_.Insert(ev.side, handle);
    }

    ctx_.ResetActive();

    W_ = W;
    weights_valid_ = true;
    cache_valid_ = true;

    budget_used_ = mem_full + static_cast<u64>(heap.Size());
    if (budget_used_ > budget_B_) budget_used_ = budget_B_;  // defensive

    if (out) *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursHighDimsAdaptiveBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "OursHighDimsAdaptiveBaseline::Sample: null rng/out";
      return false;
    }
    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "OursHighDimsAdaptiveBaseline::Sample: run.t too large (must fit in u32)";
      return false;
    }
    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    if (t == 0) return true;

    // Ensure pass-1 state.
    if (!weights_valid_ || !cache_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/rng, &tmp, phases, err)) return false;
    }
    if (W_ == 0) return true;

    const usize E = w_total_.size();

    // -----------------
    // Phase 2: assign slots to events using alias over w_total_.
    // -----------------
    struct SlotAssign {
      u32 sid;
      u32 slot;
    };
    auto by_sid_slot = [](const SlotAssign& a, const SlotAssign& b) {
      if (a.sid < b.sid) return true;
      if (b.sid < a.sid) return false;
      return a.slot < b.slot;
    };

    sampling::AliasTable alias;
    {
      auto _ = phases ? phases->Scoped("phase2_alias") : PhaseRecorder::ScopedPhase(nullptr, "");
      if (!alias.BuildFromU64(Span<const u64>(w_total_), err)) {
        if (err && err->empty()) *err = "OursHighDimsAdaptiveBaseline::Sample: failed to build alias";
        return false;
      }
    }

    // Independent seeds for (assignment, cache sampling, residual sweep).
    const u64 seed_assign = rng->NextU64();
    const u64 seed_cache  = rng->NextU64();
    const u64 seed_sweep  = rng->NextU64();

    Rng rng_assign(seed_assign);

    std::vector<SlotAssign> asg;
    asg.reserve(static_cast<usize>(t));
    for (u32 j = 0; j < t; ++j) {
      const u32 sid = static_cast<u32>(alias.Sample(&rng_assign));
      asg.push_back(SlotAssign{sid, j});
    }
    std::sort(asg.begin(), asg.end(), by_sid_slot);

    out->pairs.assign(static_cast<usize>(t), PairId{});

    const auto* ds = ctx_.dataset();
    SJS_DASSERT(ds != nullptr);

    // Residual assignments that still need sampling in pass 2.
    std::vector<SlotAssign> residual;
    residual.reserve(asg.size());

    // Fill from full cache / prefetched samples; record residual.
    {
      auto _ = phases ? phases->Scoped("phase2_fill") : PhaseRecorder::ScopedPhase(nullptr, "");

      usize ptr = 0;
      while (ptr < asg.size()) {
        const u32 sid = asg[ptr].sid;
        usize end = ptr + 1;
        while (end < asg.size() && asg[end].sid == sid) ++end;

        const usize sid_u = static_cast<usize>(sid);
        const u32 k = static_cast<u32>(end - ptr);

        if (sid_u >= E) {
          if (err) *err = "OursHighDimsAdaptiveBaseline::Sample: sid out of range";
          return false;
        }

        const join::Side side = start_side_[sid_u];
        const u32 q_idx = start_index_[sid_u];
        const Id q_id = (side == join::Side::R) ? ds->R.GetId(static_cast<usize>(q_idx))
                                                : ds->S.GetId(static_cast<usize>(q_idx));

        if (cached_[sid_u]) {
          const u64 len = cache_len_[sid_u];
          const u64 off = cache_off_[sid_u];
          if (len == 0) {
            if (err) *err = "OursHighDimsAdaptiveBaseline::Sample: cached sid has empty cache";
            return false;
          }

          Rng crng(DeriveSeed(seed_cache, 0xCA11CEULL, static_cast<u64>(sid)));
          for (u32 j = 0; j < k; ++j) {
            const u32 slot = asg[ptr + j].slot;
            const u64 pick = crng.UniformU64(len);
            const u32 oh = cache_partners_[static_cast<usize>(off + pick)];

            if (side == join::Side::R) {
              out->pairs[static_cast<usize>(slot)] = PairId{q_id, ds->S.GetId(static_cast<usize>(oh))};
            } else {
              out->pairs[static_cast<usize>(slot)] = PairId{ds->R.GetId(static_cast<usize>(oh)), q_id};
            }
          }

          ptr = end;
          continue;
        }

        // Not fully cached: consume prefetched samples if any.
        const auto& pref = prefetch_partners_[sid_u];
        const u32 s = static_cast<u32>(pref.size());
        const u32 a = (s < k) ? s : k;

        for (u32 j = 0; j < a; ++j) {
          const u32 slot = asg[ptr + j].slot;
          const u32 oh = pref[static_cast<usize>(j)];
          if (side == join::Side::R) {
            out->pairs[static_cast<usize>(slot)] = PairId{q_id, ds->S.GetId(static_cast<usize>(oh))};
          } else {
            out->pairs[static_cast<usize>(slot)] = PairId{ds->R.GetId(static_cast<usize>(oh)), q_id};
          }
        }

        // Residual slots for pass 2.
        for (u32 j = a; j < k; ++j) {
          residual.push_back(SlotAssign{sid, asg[ptr + j].slot});
        }

        ptr = end;
      }
    }

    if (residual.empty()) {
      return true;  // one-pass completion
    }

    // -----------------
    // Build residual per-event slot plan for pass 2.
    // -----------------
    std::sort(residual.begin(), residual.end(), by_sid_slot);

    detail::SlotPlanND plan;
    plan.offset.assign(E + 1, 0U);

    for (const auto& x : residual) {
      const usize sid_u = static_cast<usize>(x.sid);
      if (sid_u < E) plan.offset[sid_u + 1]++;
    }
    for (usize i = 0; i < E; ++i) {
      plan.offset[i + 1] += plan.offset[i];
    }

    plan.slots.resize(residual.size());
    {
      std::vector<u32> cursor = plan.offset;
      for (const auto& x : residual) {
        const usize sid_u = static_cast<usize>(x.sid);
        const u32 pos = cursor[sid_u]++;
        plan.slots[static_cast<usize>(pos)] = x.slot;
      }
    }

    // -----------------
    // Pass 2: second sweep, generate residual samples only.
    // -----------------
    {
      auto _ = phases ? phases->Scoped("phase3_fill_residual") : PhaseRecorder::ScopedPhase(nullptr, "");
      ctx_.ResetActive();

      auto& ar = ctx_.active_r();
      auto& as = ctx_.active_s();
      const auto& events = ctx_.events();
      const auto& sid_of_pos = ctx_.start_id_of_event();

      std::vector<u32> sampled;
      sampled.reserve(1024);

      for (usize pos = 0; pos < events.size(); ++pos) {
        const auto& ev = events[pos];
        const u32 handle = static_cast<u32>(ev.index);

        if (ev.kind == join::EventKind::End) {
          ctx_.Delete(ev.side, handle);
          continue;
        }

        const i32 sid_i32 = sid_of_pos[pos];
        SJS_DASSERT(sid_i32 >= 0);
        const u32 sid = static_cast<u32>(sid_i32);
        const usize sid_u = static_cast<usize>(sid);

        const u32 begin = plan.offset[sid_u];
        const u32 end = plan.offset[sid_u + 1];
        const u32 k = end - begin;

        const bool q_is_r = (ev.side == join::Side::R);
        const BoxT& q = ctx_.BoxOf(ev);
        const MF& other = q_is_r ? as : ar;

        if (k > 0) {
          const u64 ev_seed = DeriveSeed(seed_sweep, 0x5EEDULL, static_cast<u64>(sid));
          Rng ev_rng(DeriveSeed(ev_seed, 1));

          sampled.clear();
          other.Sample(q, static_cast<u64>(k), &ev_rng, &sampled);
          if (sampled.size() != static_cast<usize>(k)) {
            if (err) *err = "OursHighDimsAdaptiveBaseline::Sample: residual SAMPLE failed";
            ctx_.ResetActive();
            return false;
          }

          for (u32 i = 0; i < k; ++i) {
            const u32 slot = plan.slots[static_cast<usize>(begin + i)];
            const u32 oh = sampled[static_cast<usize>(i)];

            if (q_is_r) {
              out->pairs[static_cast<usize>(slot)] =
                  PairId{ds->R.GetId(ev.index), ds->S.GetId(static_cast<usize>(oh))};
            } else {
              out->pairs[static_cast<usize>(slot)] =
                  PairId{ds->R.GetId(static_cast<usize>(oh)), ds->S.GetId(ev.index)};
            }
          }
        }

        // Insert q after query.
        ctx_.Insert(ev.side, handle);
      }

      ctx_.ResetActive();
    }

    return true;
  }

  std::unique_ptr<IJoinEnumerator> Enumerate(const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)phases;
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursHighDimsAdaptiveBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::OursReportJoinEnumeratorND<Dim, T>>(&ctx_);
  }

 private:
  bool built_{false};
  detail::OursHighDimsContext<Dim, T> ctx_;

  // sid -> START event identity
  std::vector<join::Side> start_side_;
  std::vector<u32> start_index_;

  // weights
  std::vector<u64> w_total_;
  u64 W_{0};
  bool weights_valid_{false};

  // full cache
  std::vector<u8> cached_;
  std::vector<u64> cache_off_;
  std::vector<u64> cache_len_;
  std::vector<u32> cache_partners_;

  // prefetch sample cache
  std::vector<u32> prefetch_keep_;
  std::vector<std::vector<u32>> prefetch_partners_;
  bool cache_valid_{false};

  // knobs
  u64 budget_B_{0};
  u64 budget_used_{0};
  u64 w_small_{0};
};

}  // namespace highdims
}  // namespace ours
}  // namespace baselines
}  // namespace sjs
