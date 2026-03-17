#pragma once
// baselines/ours/detail/local_mode_count_index.h
//
// The high-dimensional SJS path now counts directly on the same fixed-mode local
// recursive structures used for REPORT/SAMPLE. We keep this header as a stable
// include point to preserve the original file layout.

#include "baselines/ours/detail/local_mode_index.h"

namespace sjs {
namespace baselines {
namespace ours {
namespace detail {

// Intentionally empty: Dim >= 3 phase-1 counting reuses FixedModeLocalIndexND
// through ModeFamilyND::CountMode().

}  // namespace detail
}  // namespace ours
}  // namespace baselines
}  // namespace sjs
