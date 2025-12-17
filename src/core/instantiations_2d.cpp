// src/core/instantiations_2d.cpp
//
// Optional compile unit for core headers.
//
// This project is template-heavy (Dim-parametrized geometry/index/join/sampling).
// Many pieces are header-only. Having a small "core TU" can help:
//   - catch header-only compilation issues early,
//   - centralize any future explicit instantiations if needed.
//
// Currently, most core utilities are non-templated, so this file mainly serves
// as an anchor compilation unit.

#include "sjs/core/assert.h"
#include "sjs/core/config.h"
#include "sjs/core/logging.h"
#include "sjs/core/rng.h"
#include "sjs/core/stats.h"
#include "sjs/core/timer.h"
#include "sjs/core/types.h"

namespace sjs {
namespace {

// A tiny sanity check to ensure core typedefs remain consistent.
static_assert(sizeof(Id) == 4, "Id is expected to be 32-bit for memory efficiency");
static_assert(kDefaultDim == 2, "Default Dim expected to be 2 in current project phase");

}  // namespace
}  // namespace sjs
