#pragma once

#include <vector>
#include <random>
#include <cstddef>

#include "geometry.hpp"

// =============================
// Baseline1：暴力枚举 + 均匀采样
// =============================
//
// 提供两个核心函数：
//   1) baseline1_enumerate(...)  显式枚举所有跨集合相交对 J
//   2) baseline_sample_uniform(...)  在 J 上 i.i.d. 均匀采样 t 次
//
// 说明：
//  - all_pairs_out / samples_out 在函数内部会被 clear()。
//  - 交给 main.cpp 自己计时、统计内存峰值。

/// 暴力枚举所有相交对，输出为 RectPair(idx_c, idx_bar)
inline void baseline1_enumerate(
    const std::vector<Rect>& Rc,
    const std::vector<Rect>& Rbar,
    std::vector<RectPair>& all_pairs_out
) {
    all_pairs_out.clear();
    const std::size_t n1 = Rc.size();
    const std::size_t n2 = Rbar.size();

    // 粗略预留一点空间，避免频繁扩容（完全可以视情况调整或删掉）
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

/// 在 all_pairs 上做 i.i.d. 均匀采样 t 次（有放回）
/// URNG 可以是 std::mt19937_64 等
template <class URNG>
inline void baseline_sample_uniform(
    const std::vector<RectPair>& all_pairs,
    std::size_t t,
    URNG& rng,
    std::vector<RectPair>& samples_out
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
