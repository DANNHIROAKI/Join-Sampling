#pragma once
// sjs/baselines/ours/highdims/sampling.h
//
// Our method (HighDims) — Sampling variant (SJS-HighDims Framework II).
//
// Framework II (2-pass):
//   Pass 1: sweep on axis=0, for each START event e compute w_e = COUNT(K_e)
//   Pass 2: assign t output slots to events with Pr(e)=w_e/W (alias), then
//           sweep again and for each START(e) generate k_e samples via SAMPLE(K_e,k_e).
//
// Correctness: outputs t i.i.d. uniform samples from the join J with replacement.
//
// Notes:
//  - This implementation relies on index::ModeFamily<Dim,T> which provides
//    event-block primitives COUNT/REPORT/SAMPLE for the projected (Dim-1)-D query
//    and supports dynamic active set (Insert/Delete).
//  - Enumerate() returns a REPORT-based deterministic enumerator that follows
//    the same event-block decomposition.

#include "sjs/baselines/baseline_api.h"

#include "sjs/core/assert.h"
#include "sjs/core/rng.h"
#include "sjs/io/dataset.h"
#include "sjs/join/sweep_events.h"
#include "sjs/index/mode_family.h"
#include "sjs/sampling/alias_table.h"

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
namespace ours {
namespace highdims {

namespace detail {

// -----------------------------------------------------------------------------
// HighDims shared context: events + two ModeFamily indices (R-active, S-active).
// -----------------------------------------------------------------------------
template <int Dim, class T = Scalar>
class OursHighDimsContext {
 public:
  static_assert(Dim >= 2, "HighDims SJS requires Dim >= 2");

  using DatasetT = Dataset<Dim, T>;
  using BoxT = Box<Dim, T>;
  using MF = index::ModeFamily<Dim, T>;

  // Position sentinel for the active-handle vectors.
  static constexpr u32 kNoPos = std::numeric_limits<u32>::max();

  void Reset() {
    ds_ = nullptr;
    built_ = false;

    events_.clear();
    start_id_of_event_.clear();
    start_event_pos_.clear();

    mf_r_ = MF{};
    mf_s_ = MF{};

    active_r_.clear();
    active_s_.clear();
    pos_r_.clear();
    pos_s_.clear();
  }

  bool Build(const DatasetT& ds, PhaseRecorder* phases, std::string* err) {
    Reset();
    ds_ = &ds;

    // B-mode relies on proper (non-empty) half-open boxes. We validate early to
    // avoid subtle rank-space corner cases when lo >= hi.
    {
      std::string tmp;
      if (!ds.Validate(/*require_proper=*/true, &tmp)) {
        if (err) *err = "OursHighDimsContext: invalid dataset: " + tmp;
        return false;
      }
    }

    // Basic size guards (we use u32 handles).
    const usize nR = ds.R.Size();
    const usize nS = ds.S.Size();
    if (nR > static_cast<usize>(std::numeric_limits<u32>::max()) ||
        nS > static_cast<usize>(std::numeric_limits<u32>::max())) {
      if (err) *err = "OursHighDimsContext: relation size exceeds u32 handle range";
      return false;
    }
    if ((nR + nS) > static_cast<usize>(std::numeric_limits<u32>::max())) {
      // Alias table + slot planning use u32 event ids.
      if (err) *err = "OursHighDimsContext: number of START events exceeds u32 range";
      return false;
    }

    // 1) Build sweep events on axis=0.
    {
      auto scoped = phases ? phases->Scoped("build_events") : PhaseRecorder::ScopedPhase(nullptr, "");
      events_ = join::BuildSweepEvents<Dim, T>(ds, /*axis=*/0, join::SideTieBreak::RBeforeS);
    }

    // 2) Build START id mapping (dense 0..E-1).
    start_id_of_event_.assign(events_.size(), -1);
    start_event_pos_.clear();
    start_event_pos_.reserve(nR + nS);
    for (usize pos = 0; pos < events_.size(); ++pos) {
      if (events_[pos].kind == join::EventKind::Start) {
        start_id_of_event_[pos] = static_cast<i32>(start_event_pos_.size());
        start_event_pos_.push_back(pos);
      }
    }

    // 3) Build ModeFamily universes (static) for R and S.
    {
      auto scoped = phases ? phases->Scoped("build_mode_family") : PhaseRecorder::ScopedPhase(nullptr, "");
      // Universe for mf_r_ is ds.R.boxes; for mf_s_ is ds.S.boxes.
      mf_r_.Init(&ds.R.boxes);
      mf_s_.Init(&ds.S.boxes);

      // Active-set tracking (kept in sync via Insert*/Delete* wrappers below).
      pos_r_.assign(nR, kNoPos);
      pos_s_.assign(nS, kNoPos);
      active_r_.clear();
      active_s_.clear();
      active_r_.reserve(nR);
      active_s_.reserve(nS);
    }

    built_ = true;
    return true;
  }

  bool built() const noexcept { return built_; }
  const DatasetT* dataset() const noexcept { return ds_; }

  const std::vector<join::Event>& events() const noexcept { return events_; }
  const std::vector<i32>& start_id_of_event() const noexcept { return start_id_of_event_; }
  usize num_start_events() const noexcept { return start_event_pos_.size(); }

  MF& active_r() noexcept { return mf_r_; }
  MF& active_s() noexcept { return mf_s_; }
  const MF& active_r() const noexcept { return mf_r_; }
  const MF& active_s() const noexcept { return mf_s_; }

  u64 active_size_r() const noexcept { return static_cast<u64>(active_r_.size()); }
  u64 active_size_s() const noexcept { return static_cast<u64>(active_s_.size()); }

  // Insert handle into the active set (and all modes). Handles are dense indices.
  void InsertR(u32 h) {
    SJS_DASSERT(h < pos_r_.size());
    if (h >= pos_r_.size()) return;
    SJS_DASSERT(pos_r_[h] == kNoPos);
    if (pos_r_[h] != kNoPos) return;  // defensive

    mf_r_.Insert(static_cast<Id>(h));
    pos_r_[h] = static_cast<u32>(active_r_.size());
    active_r_.push_back(h);
  }

  void InsertS(u32 h) {
    SJS_DASSERT(h < pos_s_.size());
    if (h >= pos_s_.size()) return;
    SJS_DASSERT(pos_s_[h] == kNoPos);
    if (pos_s_[h] != kNoPos) return;

    mf_s_.Insert(static_cast<Id>(h));
    pos_s_[h] = static_cast<u32>(active_s_.size());
    active_s_.push_back(h);
  }

  // Delete handle from the active set (O(#occurrences(handle))). Safe to call on
  // non-active handles (becomes a no-op in the underlying OccPool).
  void DeleteR(u32 h) {
    SJS_DASSERT(h < pos_r_.size());
    if (h >= pos_r_.size()) return;

    mf_r_.Delete(static_cast<Id>(h));

    const u32 p = pos_r_[h];
    if (p == kNoPos) return;
    SJS_DASSERT(p < active_r_.size());

    const u32 last = active_r_.back();
    active_r_[static_cast<usize>(p)] = last;
    pos_r_[last] = p;
    active_r_.pop_back();
    pos_r_[h] = kNoPos;
  }

  void DeleteS(u32 h) {
    SJS_DASSERT(h < pos_s_.size());
    if (h >= pos_s_.size()) return;

    mf_s_.Delete(static_cast<Id>(h));

    const u32 p = pos_s_[h];
    if (p == kNoPos) return;
    SJS_DASSERT(p < active_s_.size());

    const u32 last = active_s_.back();
    active_s_[static_cast<usize>(p)] = last;
    pos_s_[last] = p;
    active_s_.pop_back();
    pos_s_[h] = kNoPos;
  }

  void Insert(join::Side side, u32 h) {
    if (side == join::Side::R) InsertR(h);
    else InsertS(h);
  }

  void Delete(join::Side side, u32 h) {
    if (side == join::Side::R) DeleteR(h);
    else DeleteS(h);
  }

  // Reset dynamic active state.
  // This performs correct deletion (via OccPool-linked occurrences), rather than
  // "clearing" internal arrays which would leave dangling occurrence ids in buckets.
  void ResetActive() {
    while (!active_r_.empty()) {
      const u32 h = active_r_.back();
      active_r_.pop_back();
      pos_r_[h] = kNoPos;
      mf_r_.Delete(static_cast<Id>(h));
    }
    while (!active_s_.empty()) {
      const u32 h = active_s_.back();
      active_s_.pop_back();
      pos_s_[h] = kNoPos;
      mf_s_.Delete(static_cast<Id>(h));
    }
  }

  // Get the box referenced by an event.
  const BoxT& BoxOf(const join::Event& ev) const {
    SJS_DASSERT(ds_ != nullptr);
    if (ev.side == join::Side::R) {
      return ds_->R.boxes[ev.index];
    }
    return ds_->S.boxes[ev.index];
  }

 private:
  const DatasetT* ds_{nullptr};
  bool built_{false};

  std::vector<join::Event> events_;
  std::vector<i32> start_id_of_event_;  // size = events_.size(); -1 for non-START
  std::vector<usize> start_event_pos_;  // positions of START events

  MF mf_r_;
  MF mf_s_;

  // Active-handle tracking (to implement ResetActive safely without full rebuild).
  std::vector<u32> active_r_;
  std::vector<u32> active_s_;
  std::vector<u32> pos_r_;  // handle -> position in active_r_ (or kNoPos)
  std::vector<u32> pos_s_;  // handle -> position in active_s_ (or kNoPos)
};

// -----------------------------------------------------------------------------
// Slot plan for Framework II (per-event only, no 2D pattern split).
// -----------------------------------------------------------------------------
struct SlotPlanND {
  std::vector<u32> offset;  // size E+1
  std::vector<u32> slots;   // size t

  void Clear() {
    offset.clear();
    slots.clear();
  }
};

inline bool BuildSlotPlanND(u32 t,
                           Rng* rng,
                           const std::vector<u64>& w_total,
                           SlotPlanND* plan,
                           std::string* err) {
  if (!rng || !plan) {
    if (err) *err = "BuildSlotPlanND: null rng/plan";
    return false;
  }
  plan->Clear();

  const usize E = w_total.size();
  if (E == 0) {
    if (err) *err = "BuildSlotPlanND: empty event weight vector";
    return false;
  }
  if (E > static_cast<usize>(std::numeric_limits<u32>::max())) {
    if (err) *err = "BuildSlotPlanND: E exceeds u32";
    return false;
  }

  sampling::AliasTable alias;
  if (!alias.BuildFromU64(Span<const u64>(w_total), err)) {
    if (err && err->empty()) *err = "BuildSlotPlanND: alias build failed";
    return false;
  }

  // First pass: draw event id per slot and count per event.
  std::vector<u32> eid_of_slot;
  eid_of_slot.resize(t);

  std::vector<u32> cnt(E, 0U);
  for (u32 j = 0; j < t; ++j) {
    const u32 eid = static_cast<u32>(alias.Sample(rng));
    eid_of_slot[j] = eid;
    ++cnt[static_cast<usize>(eid)];
  }

  // Prefix sums -> offsets.
  plan->offset.resize(E + 1);
  plan->offset[0] = 0;
  for (usize e = 0; e < E; ++e) {
    plan->offset[e + 1] = plan->offset[e] + cnt[e];
  }

  plan->slots.resize(t);
  std::vector<u32> cur = plan->offset;

  // Stable-fill slots.
  for (u32 j = 0; j < t; ++j) {
    const u32 eid = eid_of_slot[j];
    const u32 pos = cur[static_cast<usize>(eid)]++;
    plan->slots[static_cast<usize>(pos)] = j;
  }

  return true;
}

// -----------------------------------------------------------------------------
// REPORT-based deterministic join enumerator (HighDims).
// -----------------------------------------------------------------------------
template <int Dim, class T>
class OursReportJoinEnumeratorND final : public baselines::IJoinEnumerator {
 public:
  using Ctx = OursHighDimsContext<Dim, T>;
  using BoxT = Box<Dim, T>;
  using MF = index::ModeFamily<Dim, T>;

  explicit OursReportJoinEnumeratorND(Ctx* ctx) : ctx_(ctx) {
    SJS_DASSERT(ctx_ != nullptr);
    Reset();
  }

  void Reset() override {
    SJS_DASSERT(ctx_ != nullptr);
    ctx_->ResetActive();

    pos_ = 0;
    stage_ = Stage::Scan;
    tmp_.clear();
    tmp_i_ = 0;

    stats_ = join::JoinStats{};
  }

  bool Next(PairId* out) override {
    if (!out) return false;
    SJS_DASSERT(ctx_ != nullptr);

    const auto* ds = ctx_->dataset();
    SJS_DASSERT(ds != nullptr);

    auto& ar = ctx_->active_r();
    auto& as = ctx_->active_s();

    const auto& events = ctx_->events();

    while (true) {
      // Emit stage: output buffered partners for the current START.
      if (stage_ == Stage::Emit) {
        if (tmp_i_ < tmp_.size()) {
          const u32 oh = tmp_[tmp_i_++];
          if (q_is_r_) {
            *out = PairId{ds->R.GetId(static_cast<usize>(q_idx_)),
                          ds->S.GetId(static_cast<usize>(oh))};
          } else {
            *out = PairId{ds->R.GetId(static_cast<usize>(oh)),
                          ds->S.GetId(static_cast<usize>(q_idx_))};
          }
          ++stats_.output_pairs;
          return true;
        }

        // Done emitting this START: resume scan.
        stage_ = Stage::Scan;
        continue;
      }

      // Scan stage.
      if (pos_ >= events.size()) return false;

      const auto& ev = events[pos_++];
      ++stats_.num_events;

      const u32 handle = static_cast<u32>(ev.index);

      if (ev.kind == join::EventKind::End) {
        ctx_->Delete(ev.side, handle);
        continue;
      }

      // START: query before inserting.
      q_is_r_ = (ev.side == join::Side::R);
      q_idx_ = handle;

      const BoxT& q = ctx_->BoxOf(ev);
      const MF& other = q_is_r_ ? as : ar;

      tmp_.clear();
      other.Report(q, &tmp_);
      tmp_i_ = 0;

      stats_.candidate_checks += static_cast<u64>(tmp_.size());

      // Insert q after reporting.
      ctx_->Insert(ev.side, handle);

      // Track peak active set sizes (after insertion).
      stats_.active_max_r = std::max(stats_.active_max_r, ctx_->active_size_r());
      stats_.active_max_s = std::max(stats_.active_max_s, ctx_->active_size_s());

      // Emit buffered partners.
      stage_ = Stage::Emit;
      // Loop; if tmp_ empty we'll continue scanning.
    }
  }

  const join::JoinStats& Stats() const noexcept override { return stats_; }

 private:
  enum class Stage : u8 { Scan = 0, Emit = 1 };

  Ctx* ctx_{nullptr};

  usize pos_{0};
  Stage stage_{Stage::Scan};

  // Current START.
  bool q_is_r_{true};
  u32 q_idx_{0};

  // Buffer.
  std::vector<u32> tmp_;
  usize tmp_i_{0};

  join::JoinStats stats_;
};

}  // namespace detail

// -----------------------------------------------------------------------------
// OursHighDimsSamplingBaseline (Framework II).
// -----------------------------------------------------------------------------
template <int Dim, class T = Scalar>
class OursHighDimsSamplingBaseline final : public IBaseline<Dim, T> {
 public:
  static_assert(Dim >= 2, "OursHighDimsSamplingBaseline requires Dim >= 2");

  using Base = IBaseline<Dim, T>;
  using DatasetT = typename Base::DatasetT;
  using BoxT = Box<Dim, T>;
  using MF = index::ModeFamily<Dim, T>;

  Method method() const noexcept override { return Method::Ours; }
  Variant variant() const noexcept override { return Variant::Sampling; }
  std::string_view Name() const noexcept override { return "ours_highdims_sampling"; }

  void Reset() override {
    built_ = false;
    ctx_.Reset();
    w_total_.clear();
    W_ = 0;
    weights_valid_ = false;
  }

  bool Build(const DatasetT& ds, const Config& cfg, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    Reset();
    if (!ctx_.Build(ds, phases, err)) return false;

    const usize E = ctx_.num_start_events();
    w_total_.assign(E, 0ULL);

    built_ = true;
    return true;
  }

  bool Count(const Config& cfg, Rng* rng, CountResult* out, PhaseRecorder* phases, std::string* err) override {
    (void)cfg;
    (void)rng;

    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursHighDimsSamplingBaseline::Count: call Build() first";
      return false;
    }

    auto scoped = phases ? phases->Scoped("phase1_count") : PhaseRecorder::ScopedPhase(nullptr, "");

    std::fill(w_total_.begin(), w_total_.end(), 0ULL);
    ctx_.ResetActive();

    u64 W = 0;

    auto& ar = ctx_.active_r();
    auto& as = ctx_.active_s();

    const auto& events = ctx_.events();
    const auto& sid_of_pos = ctx_.start_id_of_event();

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

      const bool q_is_r = (ev.side == join::Side::R);
      const BoxT& q = ctx_.BoxOf(ev);
      const MF& other = q_is_r ? as : ar;

      const u64 w = other.Count(q);
      w_total_[static_cast<usize>(sid)] = w;

      // Accumulate W with overflow guard.
      if (w > 0) {
        if (W > (std::numeric_limits<u64>::max() - w)) {
          if (err) *err = "OursHighDimsSamplingBaseline::Count: |J| overflowed u64";
          ctx_.ResetActive();
          return false;
        }
        W += w;
      }

      // Insert q into its side after querying.
      ctx_.Insert(ev.side, handle);
    }

    ctx_.ResetActive();

    W_ = W;
    weights_valid_ = true;

    if (out) *out = MakeExactCount(W_);
    return true;
  }

  bool Sample(const Config& cfg, Rng* rng, SampleSet* out, PhaseRecorder* phases, std::string* err) override {
    if (!built_ || !ctx_.built() || ctx_.dataset() == nullptr) {
      if (err) *err = "OursHighDimsSamplingBaseline::Sample: call Build() first";
      return false;
    }
    if (!rng || !out) {
      if (err) *err = "OursHighDimsSamplingBaseline::Sample: null rng/out";
      return false;
    }

    if (cfg.run.t > static_cast<u64>(std::numeric_limits<u32>::max())) {
      if (err) *err = "OursHighDimsSamplingBaseline::Sample: run.t too large (must fit in u32)";
      return false;
    }
    const u32 t = static_cast<u32>(cfg.run.t);

    out->Clear();
    out->with_replacement = true;
    out->weighted = false;
    if (t == 0) return true;

    // ---------------------------------------------------------------------
    // Randomness protocol (SJS-HighDims §3.1.4):
    //   seed(i, phase_id, ctr) = Hash(master_seed, i, phase_id, ctr)
    //
    // We split randomness into two independent stages:
    //   (1) assign: pick which START event owns each output slot
    //   (2) sweep : for each START event e_i, sample within K_e using a
    //               per-event RNG derived from i (start-event id)
    //
    // This improves reproducibility: changes in slot assignment code do not
    // perturb within-event sampling randomness, and vice versa.
    // ---------------------------------------------------------------------

    // Establish a deterministic per-call master seed WITHOUT coupling it to
    // downstream RNG draw counts.
    const u64 master_seed = rng->NextU64();

    // Fixed phase identifiers for this baseline variant.
    static constexpr u64 kPhaseAssign = 1;  // slot->event assignment
    static constexpr u64 kPhaseSweep  = 2;  // within-event partner sampling
    static constexpr u64 kCtr0 = 0;

    // Ensure weights exist.
    if (!weights_valid_) {
      CountResult tmp;
      if (!Count(cfg, /*rng=*/nullptr, &tmp, phases, err)) return false;
    }
    if (W_ == 0) {
      // Empty join.
      return true;
    }

    // Phase2: slot plan (assignment stage).
    detail::SlotPlanND plan;
    {
      auto scoped = phases ? phases->Scoped("phase2_plan") : PhaseRecorder::ScopedPhase(nullptr, "");
      // Independent randomness stream for assignment.
      Rng rng_assign(DeriveSeed(master_seed, /*i=*/0, kPhaseAssign, kCtr0));
      if (!detail::BuildSlotPlanND(t, &rng_assign, w_total_, &plan, err)) return false;
    }

    out->pairs.resize(static_cast<usize>(t));

    // Phase3: second sweep, generate samples per event in batches.
    {
      auto scoped = phases ? phases->Scoped("phase3_sample") : PhaseRecorder::ScopedPhase(nullptr, "");
      ctx_.ResetActive();

      auto& ar = ctx_.active_r();
      auto& as = ctx_.active_s();

      const auto& events = ctx_.events();
      const auto& sid_of_pos = ctx_.start_id_of_event();
      const auto* ds = ctx_.dataset();
      SJS_DASSERT(ds != nullptr);

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
          sampled.clear();
          // Per-event RNG derived from the START-event id.
          // This makes each event's sampling independent from the others.
          Rng rng_event(DeriveSeed(master_seed, static_cast<u64>(sid), kPhaseSweep, kCtr0));
          other.Sample(q, static_cast<u64>(k), &rng_event, &sampled);
          if (sampled.size() != static_cast<usize>(k)) {
            if (err) *err = "OursHighDimsSamplingBaseline::Sample: SAMPLE failed (weights/index inconsistent)";
            ctx_.ResetActive();
            return false;
          }

          for (u32 i = 0; i < k; ++i) {
            const u32 slot = plan.slots[static_cast<usize>(begin + i)];
            const u32 oh = sampled[static_cast<usize>(i)];
            if (q_is_r) {
              out->pairs[static_cast<usize>(slot)] = PairId{ds->R.GetId(ev.index), ds->S.GetId(static_cast<usize>(oh))};
            } else {
              out->pairs[static_cast<usize>(slot)] = PairId{ds->R.GetId(static_cast<usize>(oh)), ds->S.GetId(ev.index)};
            }
          }
        }

        // Insert q.
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
      if (err) *err = "OursHighDimsSamplingBaseline::Enumerate: call Build() first";
      return nullptr;
    }
    return std::make_unique<detail::OursReportJoinEnumeratorND<Dim, T>>(&ctx_);
  }

 private:
    bool built_{false};
    detail::OursHighDimsContext<Dim, T> ctx_;

    std::vector<u64> w_total_;
    u64 W_{0};
    bool weights_valid_{false};
};

}  // namespace highdims
}  // namespace ours
}  // namespace baselines
}  // namespace sjs
