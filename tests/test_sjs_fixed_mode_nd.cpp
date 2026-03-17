// tests/test_sjs_fixed_mode_nd.cpp
//
// Direct validation of the fixed-mode local recursive structures used by the
// high-dimensional SJS path. For each START event we compare per-mode counts
// against a naive active-set classifier on the current sweep state.

#include "baselines/ours/detail/context_nd.h"
#include "baselines/ours/detail/mode_mask.h"

#include "core/timer.h"
#include "geometry/box.h"
#include "io/dataset.h"

#include <array>
#include <iostream>
#include <string>

namespace {

struct TestContext {
  int fails = 0;
  void Check(bool ok, const char* expr, const char* file, int line) {
    if (ok) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line << "  CHECK(" << expr << ")\n";
  }
  template <class A, class B>
  void CheckEq(const A& a, const B& b, const char* ea, const char* eb,
               const char* file, int line) {
    if (a == b) return;
    ++fails;
    std::cerr << "[FAIL] " << file << ":" << line << "  CHECK_EQ(" << ea << ", " << eb
              << ") got " << a << " vs " << b << "\n";
  }
};

#define CHECK(ctx, expr) (ctx).Check((expr), #expr, __FILE__, __LINE__)
#define CHECK_EQ(ctx, a, b) (ctx).CheckEq((a), (b), #a, #b, __FILE__, __LINE__)

using sjs::Box;
using sjs::Dataset;
using sjs::Point;
using sjs::Scalar;
using sjs::u32;
using sjs::u64;

template <int Dim>
Box<Dim, Scalar> MakeBox(std::array<double, Dim> lo_a, std::array<double, Dim> hi_a) {
  Point<Dim, Scalar> lo{}, hi{};
  for (int i = 0; i < Dim; ++i) {
    lo[i] = lo_a[static_cast<size_t>(i)];
    hi[i] = hi_a[static_cast<size_t>(i)];
  }
  return Box<Dim, Scalar>(lo, hi);
}

template <int Dim>
Dataset<Dim, Scalar> MakeDataset() {
  Dataset<Dim, Scalar> ds;
  ds.name = std::string("mode_check_") + std::to_string(Dim);
  ds.R.name = "R";
  ds.S.name = "S";

  std::array<double, Dim> lo{}, hi{};

  for (int i = 0; i < Dim; ++i) { lo[i] = 0.0; hi[i] = (i == 0 ? 1.1 : 0.9); }
  ds.R.Add(MakeBox<Dim>(lo, hi));

  for (int i = 0; i < Dim; ++i) {
    lo[i] = (i == 0 ? 0.4 : (0.1 * ((i % 3) + 1)));
    hi[i] = lo[i] + 1.0;
  }
  ds.R.Add(MakeBox<Dim>(lo, hi));

  for (int i = 0; i < Dim; ++i) { lo[i] = 2.0; hi[i] = 2.6; }
  ds.R.Add(MakeBox<Dim>(lo, hi));

  for (int i = 0; i < Dim; ++i) {
    lo[i] = (i == 0 ? 0.2 : 0.2 + 0.1 * i);
    hi[i] = lo[i] + 0.7;
  }
  ds.S.Add(MakeBox<Dim>(lo, hi));

  for (int i = 0; i < Dim; ++i) {
    lo[i] = (i == 0 ? 0.7 : (i % 2 == 0 ? 0.8 : 0.1));
    hi[i] = lo[i] + 0.6;
  }
  ds.S.Add(MakeBox<Dim>(lo, hi));

  for (int i = 0; i < Dim; ++i) { lo[i] = 3.0; hi[i] = 3.5; }
  ds.S.Add(MakeBox<Dim>(lo, hi));

  ds.EnsureIds();
  return ds;
}

template <int Dim>
void RunDim(TestContext& t) {
  constexpr u32 kModeCount = sjs::baselines::ours::detail::kModeCountV<Dim>;
  const auto ds = MakeDataset<Dim>();
  std::string err;
  CHECK(t, ds.Validate(true, &err));

  sjs::PhaseRecorder phases;
  sjs::baselines::ours::detail::OursNDContext<Dim, Scalar> ctx;
  CHECK(t, ctx.Build(ds, &phases, &err));

  const auto& ks = ctx.ev_kind_side();
  const auto& handles = ctx.ev_handle();

  for (size_t pos = 0; pos < ks.size(); ++pos) {
    const bool is_start = ((ks[pos] >> 1U) != 0U);
    const bool is_r = ((ks[pos] & 1U) == 0U);
    const u32 handle = handles[pos];

    if (!is_start) {
      if (is_r) ctx.EraseR(handle);
      else ctx.EraseS(handle);
      continue;
    }

    const auto& q = is_r ? ctx.BoxR(handle) : ctx.BoxS(handle);
    const auto& other_rel = is_r ? ds.S : ds.R;
    const auto& other_active = is_r ? ctx.active_s().Items() : ctx.active_r().Items();
    const auto& other_modes = is_r ? ctx.modes_s() : ctx.modes_r();

    std::array<u64, kModeCount> naive{};
    naive.fill(0ULL);
    for (u32 oh : other_active) {
      u32 mode = 0;
      if (!sjs::baselines::ours::detail::ComputeModeMaskSkipAxis0<Dim, Scalar>(
              q, other_rel.boxes[static_cast<size_t>(oh)], &mode)) {
        continue;
      }
      ++naive[static_cast<size_t>(mode)];
    }

    for (u32 mode = 0; mode < kModeCount; ++mode) {
      CHECK_EQ(t, other_modes.CountMode(mode, q), naive[static_cast<size_t>(mode)]);
    }

    if (is_r) ctx.InsertR(handle);
    else ctx.InsertS(handle);
  }

  ctx.ResetActive();
}

}  // namespace

int main() {
  TestContext t;
  RunDim<3>(t);
  RunDim<4>(t);
  RunDim<5>(t);

  if (t.fails == 0) {
    std::cout << "[OK] test_sjs_fixed_mode_nd\n";
    return 0;
  }
  std::cerr << "[FAILED] test_sjs_fixed_mode_nd: " << t.fails << " failure(s)\n";
  return 1;
}
