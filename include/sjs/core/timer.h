#pragma once
// sjs/core/timer.h
//
// Timing utilities for experiments.
// Provides:
//  - Stopwatch: a simple wall-clock timer.
//  - PhaseRecorder: accumulate phase durations by name with RAII ScopedPhase.
//
// Notes:
//  - Uses std::chrono::steady_clock (monotonic).
//  - PhaseRecorder uses std::map with transparent comparator (std::less<>) so
//    lookups with std::string_view do not allocate in C++17.

#include "sjs/core/types.h"

#include <chrono>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {

class Stopwatch {
 public:
  using Clock = std::chrono::steady_clock;

  Stopwatch() : start_(Clock::now()) {}

  void Reset() { start_ = Clock::now(); }

  u64 ElapsedNanos() const {
    return static_cast<u64>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start_).count());
  }

  double ElapsedMillis() const { return static_cast<double>(ElapsedNanos()) / 1e6; }
  double ElapsedSeconds() const { return static_cast<double>(ElapsedNanos()) / 1e9; }

 private:
  Clock::time_point start_;
};

class PhaseRecorder {
 public:
  // std::less<> is transparent; std::map supports heterogeneous lookup in C++17.
  using Map = std::map<std::string, u64, std::less<>>;

  PhaseRecorder() = default;

  void Clear() { nanos_.clear(); }

  // Accumulate a duration (ns) to a phase name.
  void Add(std::string_view name, u64 nanos) {
    auto it = nanos_.find(name);
    if (it == nanos_.end()) {
      nanos_.emplace(std::string(name), nanos);
    } else {
      it->second += nanos;
    }
  }

  // Get phase time in ns; returns 0 if not found.
  u64 GetNanos(std::string_view name) const {
    auto it = nanos_.find(name);
    return (it == nanos_.end()) ? 0ULL : it->second;
  }

  double GetMillis(std::string_view name) const { return static_cast<double>(GetNanos(name)) / 1e6; }

  // Snapshot as (name, nanos) vector (already sorted by name).
  std::vector<std::pair<std::string, u64>> SnapshotNanosSorted() const {
    std::vector<std::pair<std::string, u64>> out;
    out.reserve(nanos_.size());
    for (const auto& kv : nanos_) out.emplace_back(kv.first, kv.second);
    return out;
  }

  // Convenience: produce a tiny JSON object string (no external deps).
  // Example: {"build_ms":12.3,"phase1_ms":4.5}
  std::string ToJsonMillis() const {
    std::ostringstream oss;
    oss << "{";
    usize i = 0;
    for (const auto& kv : nanos_) {
      const double ms = static_cast<double>(kv.second) / 1e6;
      oss << "\"" << kv.first << "_ms\":" << std::fixed << std::setprecision(6) << ms;
      if (++i < nanos_.size()) oss << ",";
    }
    oss << "}";
    return oss.str();
  }

  // RAII helper: measures duration and adds it to the recorder on destruction.
  class ScopedPhase {
   public:
    ScopedPhase(PhaseRecorder* rec, std::string_view name)
        : rec_(rec), name_(name), start_(Stopwatch::Clock::now()) {}

    ScopedPhase(const ScopedPhase&) = delete;
    ScopedPhase& operator=(const ScopedPhase&) = delete;

    ScopedPhase(ScopedPhase&& other) noexcept
        : rec_(other.rec_), name_(other.name_), start_(other.start_) {
      other.rec_ = nullptr;
    }
    ScopedPhase& operator=(ScopedPhase&& other) noexcept {
      if (this == &other) return *this;
      Finish();
      rec_ = other.rec_;
      name_ = other.name_;
      start_ = other.start_;
      other.rec_ = nullptr;
      return *this;
    }

    ~ScopedPhase() { Finish(); }

   private:
    void Finish() {
      if (!rec_) return;
      const auto end = Stopwatch::Clock::now();
      const u64 ns = static_cast<u64>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_).count());
      rec_->Add(name_, ns);
      rec_ = nullptr;
    }

    PhaseRecorder* rec_;
    std::string_view name_;
    Stopwatch::Clock::time_point start_;
  };

  ScopedPhase Scoped(std::string_view name) { return ScopedPhase(this, name); }

 private:
  Map nanos_;
};

}  // namespace sjs
