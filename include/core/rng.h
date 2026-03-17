#pragma once
// core/rng.h
//
// Fast reproducible RNG wrapper.
// We use splitmix64 for seeding and xoshiro256** for generation.
// This is fast, statistically solid for Monte Carlo sampling, and fully deterministic.

#include "core/types.h"

#include <array>
#include <cstdint>

namespace sjs {
namespace detail {

inline constexpr u64 Rotl(const u64 x, int k) noexcept {
  return (x << k) | (x >> (64 - k));
}

// splitmix64: simple mixing generator used to expand a single seed into
// multiple 64-bit values.
class SplitMix64 {
 public:
  explicit SplitMix64(u64 seed) noexcept : state_(seed) {}

  u64 Next() noexcept {
    u64 z = (state_ += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }

 private:
  u64 state_;
};

// xoshiro256**: https://prng.di.unimi.it/
class Xoshiro256StarStar {
 public:
  Xoshiro256StarStar() noexcept : s_{0, 0, 0, 0} {}

  void Seed(u64 seed) noexcept {
    SplitMix64 sm(seed);
    for (usize i = 0; i < 4; ++i) {
      s_[i] = sm.Next();
    }
    // all-zero state is not allowed; splitmix won't produce it in practice,
    // but keep a guard anyway.
    if ((s_[0] | s_[1] | s_[2] | s_[3]) == 0) {
      s_[0] = 0x9e3779b97f4a7c15ULL;
    }
  }

  u64 NextU64() noexcept {
    const u64 result = Rotl(s_[1] * 5ULL, 7) * 9ULL;

    const u64 t = s_[1] << 17;

    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];

    s_[2] ^= t;
    s_[3] = Rotl(s_[3], 45);

    return result;
  }

 private:
  std::array<u64, 4> s_;
};

}  // namespace detail

class Rng {
 public:
  explicit Rng(u64 seed = 0) noexcept { Seed(seed); }

  void Seed(u64 seed) noexcept { gen_.Seed(seed); }

  // Raw draws
  u64 NextU64() noexcept { return gen_.NextU64(); }
  u32 NextU32() noexcept { return static_cast<u32>(gen_.NextU64() >> 32); }

  // [0,1) with 53 bits of precision (like std::uniform_real_distribution<double>)
  double NextDouble() noexcept {
    const u64 x = NextU64();
    // Take top 53 bits and scale.
    return static_cast<double>(x >> 11) * (1.0 / 9007199254740992.0);  // 2^53
  }

  // Unbiased integer in [0, bound). bound must be > 0.
  u64 UniformU64(u64 bound) noexcept {
    // Defensive: avoid UB for bound==0.
    if (bound == 0) return 0;

#if defined(__SIZEOF_INT128__)
    // Lemire's method: https://arxiv.org/abs/1805.10941
    // Compute high 64 bits of (x * bound). Retry if low part is in the "bias" region.
    const u64 threshold = static_cast<u64>(-bound) % bound;
    while (true) {
      const u64 x = NextU64();
      const __uint128_t m = static_cast<__uint128_t>(x) * static_cast<__uint128_t>(bound);
      const u64 l = static_cast<u64>(m);
      if (l >= threshold) {
        return static_cast<u64>(m >> 64);
      }
    }
#else
    // Portable fallback: rejection sampling using modulo.
    const u64 limit = (std::numeric_limits<u64>::max() / bound) * bound;
    while (true) {
      const u64 x = NextU64();
      if (x < limit) return x % bound;
    }
#endif
  }

  u32 UniformU32(u32 bound) noexcept { return static_cast<u32>(UniformU64(bound)); }

  // Uniform in [lo, hi). If hi <= lo, returns lo.
  double UniformDouble(double lo, double hi) noexcept {
    if (!(hi > lo)) return lo;
    return lo + (hi - lo) * NextDouble();
  }

  // Bernoulli(p): true with probability p (clamped to [0,1]).
  bool Bernoulli(double p) noexcept {
    if (p <= 0.0) return false;
    if (p >= 1.0) return true;
    return NextDouble() < p;
  }

 private:
  detail::Xoshiro256StarStar gen_;
};

// -----------------------------------------------------------------------------
// Seed hashing / derivation utilities
// -----------------------------------------------------------------------------
// SJS v3 suggests deriving independent RNG streams for different phases/calls
// from a master seed. These helpers implement a lightweight, deterministic
// 64-bit mixer based on splitmix64.

// Hash a single 64-bit value to a new 64-bit value.
inline u64 HashSeed(u64 x) noexcept {
  // splitmix64 finalizer
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

// Derive a sub-seed from (master, salt).
inline u64 DeriveSeed(u64 master, u64 salt) noexcept {
  return HashSeed(master ^ (salt + 0xD1B54A32D192ED03ULL));
}

// Derive a sub-seed from (master, a, b).
inline u64 DeriveSeed(u64 master, u64 a, u64 b) noexcept {
  return DeriveSeed(DeriveSeed(master, a), b);
}

// Derive a sub-seed from (master, a, b, c).
inline u64 DeriveSeed(u64 master, u64 a, u64 b, u64 c) noexcept {
  return DeriveSeed(DeriveSeed(master, a, b), c);
}

}  // namespace sjs
