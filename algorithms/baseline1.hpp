#pragma once

#include <vector>
#include <random>
#include <cstddef>

#include "geometry.hpp"

namespace rect_sampler {
namespace baseline1 {

// 暴力枚举所有跨集合相交对 J
inline void enumerate_all_pairs(
    const RectList& Rc,
    const RectList& Rbar,
    PairList& all_pairs_out
) {
    all_pairs_out.clear();
    const std::size_t n1 = Rc.size();
    const std::size_t n2 = Rbar.size();

    if (n1 > 0 && n2 > 0) {
        all_pairs_out.reserve(std::min<std::size_t>(n1 * 4, n1 * n2));
    }

    for (std::size_t i = 0; i < n1; ++i) {
        for (std::size_t j = 0; j < n2; ++j) {
            if (rectangles_intersect(Rc[i], Rbar[j])) {
                RectPair p;
                p.idx_c   = static_cast<int>(i);
                p.idx_bar = static_cast<int>(j);
                all_pairs_out.push_back(p);
            }
        }
    }
}

// 在 all_pairs 上做 i.i.d. 均匀采样 t 次（有放回）
template <class URNG>
inline void sample_from_all_pairs(
    const PairList& all_pairs,
    std::size_t t,
    URNG& rng,
    PairList& samples_out
) {
    samples_out.clear();
    if (all_pairs.empty() || t == 0) return;

    samples_out.reserve(t);
    std::uniform_int_distribution<std::size_t> dist(0, all_pairs.size() - 1);

    for (std::size_t k = 0; k < t; ++k) {
        std::size_t idx = dist(rng);
        samples_out.push_back(all_pairs[idx]);
    }
}

} // namespace baseline1
} // namespace rect_sampler
