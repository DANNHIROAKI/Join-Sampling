#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <cstddef>
#include <cstdint>

#include "geometry.hpp"
#include "structures.hpp"
#include "utils.hpp"

// ========================================
// Ours：O(n log n + t) i.i.d. 均匀采样
// ========================================
//
// 接口：
//   ours_sample(Rc, Rbar, t, rng, samples_out);
//
// 其中：
//   - Rc / Rbar: 输入矩形集合
//   - t: 需要采样的样本数
//   - rng: 随机数引擎（例如 std::mt19937_64）
//   - samples_out: 输出长度 t 的 RectPair 序列（函数内部会 clear）

namespace ours_detail {

    enum EventType : int {
        END_EVENT   = 0,
        START_EVENT = 1
    };

    struct Event {
        double x{0.0};
        EventType type{END_EVENT};
        bool is_c{false};  // true: Rc, false: Rbar
        int  idx{-1};      // 在对应集合中的索引
        int  start_id{-1}; // 仅 START 事件使用：在 START 集合中的编号 [0, num_starts)
    };

    inline bool event_less(const Event& a, const Event& b) {
        if (a.x < b.x) return true;
        if (a.x > b.x) return false;
        // x 相同：END 在 START 前
        if (a.type != b.type) return a.type < b.type;
        // 同为 START：固定 R_c 在 R_{\bar c} 前
        if (a.type == START_EVENT && b.type == START_EVENT) {
            if (a.is_c != b.is_c) return a.is_c > b.is_c; // true 在前
        }
        return false;
    }

} // namespace ours_detail

/// Ours：在不显式枚举 J 的前提下，从 J 上 i.i.d. 均匀采样 t 次
template <class URNG>
inline void ours_sample(
    const std::vector<Rect>& Rc,
    const std::vector<Rect>& Rbar,
    std::size_t t,
    URNG& rng,
    std::vector<RectPair>& samples_out
) {
    using namespace ours_detail;

    samples_out.clear();
    if (t == 0) return;

    const int n1 = static_cast<int>(Rc.size());
    const int n2 = static_cast<int>(Rbar.size());

    // ========== 0. 构建 1D 结构的全集（interval + point） ==========

    // R_c
    std::vector<DynamicIntervalStabbing::Interval> intervals_c;
    intervals_c.reserve(n1);
    std::vector<double> points_c;
    points_c.reserve(n1);
    for (int i = 0; i < n1; ++i) {
        DynamicIntervalStabbing::Interval I;
        I.lo = Rc[i].y_min;
        I.hi = Rc[i].y_max;
        I.id = i;
        intervals_c.push_back(I);
        points_c.push_back(Rc[i].y_min);
    }

    // R_{\bar c}
    std::vector<DynamicIntervalStabbing::Interval> intervals_bar;
    intervals_bar.reserve(n2);
    std::vector<double> points_bar;
    points_bar.reserve(n2);
    for (int i = 0; i < n2; ++i) {
        DynamicIntervalStabbing::Interval I;
        I.lo = Rbar[i].y_min;
        I.hi = Rbar[i].y_max;
        I.id = i;
        intervals_bar.push_back(I);
        points_bar.push_back(Rbar[i].y_min);
    }

    DynamicIntervalStabbing stab_c(intervals_c);
    DynamicIntervalStabbing stab_bar(intervals_bar);
    DynamicRangeTree      range_c(points_c);
    DynamicRangeTree      range_bar(points_bar);

    // ========== 1. 构建事件数组并排序 ==========

    std::vector<Event> events;
    events.reserve((n1 + n2) * 2);

    for (int i = 0; i < n1; ++i) {
        Event s, e;
        s.x    = Rc[i].x_min;
        s.type = START_EVENT;
        s.is_c = true;
        s.idx  = i;

        e.x    = Rc[i].x_max;
        e.type = END_EVENT;
        e.is_c = true;
        e.idx  = i;

        events.push_back(s);
        events.push_back(e);
    }

    for (int i = 0; i < n2; ++i) {
        Event s, e;
        s.x    = Rbar[i].x_min;
        s.type = START_EVENT;
        s.is_c = false;
        s.idx  = i;

        e.x    = Rbar[i].x_max;
        e.type = END_EVENT;
        e.is_c = false;
        e.idx  = i;

        events.push_back(s);
        events.push_back(e);
    }

    std::sort(events.begin(), events.end(), event_less);

    // 给所有 START 事件分配 start_id
    int num_starts = 0;
    for (auto& ev : events) {
        if (ev.type == START_EVENT) {
            ev.start_id = num_starts++;
        }
    }
    if (num_starts == 0) {
        // 没有任何 START 事件 → 没有相交对
        samples_out.resize(0);
        return;
    }

    // ========== 2. 第一次扫描：只计数 k1, k2, w = k1 + k2 = |J_e| = |K_e^{(1)}| + |K_e^{(2)}| ==========

    std::vector<std::uint64_t> k1(num_starts, 0);
    std::vector<std::uint64_t> k2(num_starts, 0);
    std::vector<std::uint64_t> w (num_starts, 0);

    for (const auto& ev : events) {
        if (ev.type == END_EVENT) {
            if (ev.is_c) {
                stab_c.deactivate(ev.idx);
                range_c.deactivate(ev.idx);
            } else {
                stab_bar.deactivate(ev.idx);
                range_bar.deactivate(ev.idx);
            }
        } else {
            // START_EVENT
            const int e_id = ev.start_id;
            if (ev.is_c) {
                const Rect& r = Rc[ev.idx];
                double D = r.y_min;
                double U = r.y_max;

                // K_e^{(1)}：D(r_bar) <= D < U(r_bar)
                std::size_t cnt1 = stab_bar.count(D);

                // K_e^{(2)}：D < D(r_bar) < U
                std::size_t cnt2 = range_bar.count(D, U);

                k1[e_id] = static_cast<std::uint64_t>(cnt1);
                k2[e_id] = static_cast<std::uint64_t>(cnt2);
                w [e_id] = k1[e_id] + k2[e_id];

                // 激活 r ∈ R_c
                stab_c.activate(ev.idx);
                range_c.activate(ev.idx);
            } else {
                const Rect& r = Rbar[ev.idx];
                double D = r.y_min;
                double U = r.y_max;

                std::size_t cnt1 = stab_c.count(D);
                std::size_t cnt2 = range_c.count(D, U);

                k1[e_id] = static_cast<std::uint64_t>(cnt1);
                k2[e_id] = static_cast<std::uint64_t>(cnt2);
                w [e_id] = k1[e_id] + k2[e_id];

                stab_bar.activate(ev.idx);
                range_bar.activate(ev.idx);
            }
        }
    }

    // 计算全局 |J| = sum_e w_e
    std::uint64_t total_pairs = 0;
    for (int e = 0; e < num_starts; ++e) {
        total_pairs += w[e];
    }
    if (total_pairs == 0) {
        // 没有任何相交对，按约定直接返回空
        samples_out.resize(0);
        return;
    }

    // ========== 3. 第二阶段：事件级别 alias + slot 规划 ==========

    // 只在 w_e > 0 的事件上建 alias
    std::vector<int>    nonzero_events;
    std::vector<double> event_weights;
    nonzero_events.reserve(num_starts);
    event_weights.reserve(num_starts);

    for (int e = 0; e < num_starts; ++e) {
        if (w[e] > 0) {
            nonzero_events.push_back(e);
            event_weights.push_back(static_cast<double>(w[e]));
        }
    }

    AliasTable event_alias;
    event_alias.build(event_weights);

    // 为每个事件 e 和类型 g=1,2 维护 slots 集合
    std::vector<std::vector<std::size_t>> slots_k1(num_starts);
    std::vector<std::vector<std::size_t>> slots_k2(num_starts);

    std::uniform_real_distribution<double> dist01(0.0, 1.0);

    for (std::size_t j = 0; j < t; ++j) {
        // step1: 在事件集合（非零权重）上按 w_e 抽一个 e
        int col = event_alias.sample(rng);
        if (col < 0) {
            // 理论上不会发生
            continue;
        }
        int e_id = nonzero_events[col];

        std::uint64_t k1e = k1[e_id];
        std::uint64_t k2e = k2[e_id];
        std::uint64_t we  = w [e_id];

        if (we == 0) continue; // 防守式

        int group = 1;
        if (k1e == 0 && k2e == 0) {
            // 不应发生（we>0），防守
            continue;
        } else if (k1e == 0) {
            group = 2;
        } else if (k2e == 0) {
            group = 1;
        } else {
            double p = static_cast<double>(k1e) / static_cast<double>(we);
            double u = dist01(rng);
            group = (u < p ? 1 : 2);
        }

        if (group == 1) {
            slots_k1[e_id].push_back(j);
        } else {
            slots_k2[e_id].push_back(j);
        }
    }

    // ========== 4. 第二次构建 1D 结构 & 第三阶段第二次扫描 ==========

    // 重新构建 1D 结构，以清空活跃集
    stab_c.build(intervals_c);
    stab_bar.build(intervals_bar);
    range_c.build(points_c);
    range_bar.build(points_bar);

    samples_out.clear();
    samples_out.resize(t);
    // 先填充一个非法值，便于 debug 时检测
    for (std::size_t j = 0; j < t; ++j) {
        samples_out[j].idx_c   = -1;
        samples_out[j].idx_bar = -1;
    }

    std::vector<int> local_ids;   // 复用的局部采样数组

    for (const auto& ev : events) {
        if (ev.type == END_EVENT) {
            if (ev.is_c) {
                stab_c.deactivate(ev.idx);
                range_c.deactivate(ev.idx);
            } else {
                stab_bar.deactivate(ev.idx);
                range_bar.deactivate(ev.idx);
            }
        } else {
            // START_EVENT
            const int e_id = ev.start_id;
            const std::size_t need1 = slots_k1[e_id].size();
            const std::size_t need2 = slots_k2[e_id].size();

            if (ev.is_c) {
                const Rect& r = Rc[ev.idx];
                double D = r.y_min;
                double U = r.y_max;

                // ---- 类型 1：K_e^{(1)}，使用 stabbing 结构在 R_{\bar c} 上采样 need1 次 ----
                if (need1 > 0) {
                    local_ids.clear();
                    stab_bar.sample(D, need1, rng, local_ids);
                    // 写回对应 slots
                    const auto& slots = slots_k1[e_id];
                    const std::size_t m = std::min(need1, local_ids.size());
                    for (std::size_t i = 0; i < m; ++i) {
                        std::size_t j = slots[i];
                        RectPair p;
                        p.idx_c   = ev.idx;
                        p.idx_bar = local_ids[i];
                        samples_out[j] = p;
                    }
                }

                // ---- 类型 2：K_e^{(2)}，使用 range 结构在 R_{\bar c} 上采样 need2 次 ----
                if (need2 > 0) {
                    local_ids.clear();
                    range_bar.sample(D, U, need2, rng, local_ids);
                    const auto& slots = slots_k2[e_id];
                    const std::size_t m = std::min(need2, local_ids.size());
                    for (std::size_t i = 0; i < m; ++i) {
                        std::size_t j = slots[i];
                        RectPair p;
                        p.idx_c   = ev.idx;
                        p.idx_bar = local_ids[i];
                        samples_out[j] = p;
                    }
                }

                // 最后把 r 加入 R_c 的活跃集
                stab_c.activate(ev.idx);
                range_c.activate(ev.idx);
            } else {
                // r ∈ R_{\bar c}，对称处理
                const Rect& r = Rbar[ev.idx];
                double D = r.y_min;
                double U = r.y_max;

                if (need1 > 0) {
                    local_ids.clear();
                    stab_c.sample(D, need1, rng, local_ids);
                    const auto& slots = slots_k1[e_id];
                    const std::size_t m = std::min(need1, local_ids.size());
                    for (std::size_t i = 0; i < m; ++i) {
                        std::size_t j = slots[i];
                        RectPair p;
                        p.idx_c   = local_ids[i];
                        p.idx_bar = ev.idx;
                        samples_out[j] = p;
                    }
                }

                if (need2 > 0) {
                    local_ids.clear();
                    range_c.sample(D, U, need2, rng, local_ids);
                    const auto& slots = slots_k2[e_id];
                    const std::size_t m = std::min(need2, local_ids.size());
                    for (std::size_t i = 0; i < m; ++i) {
                        std::size_t j = slots[i];
                        RectPair p;
                        p.idx_c   = local_ids[i];
                        p.idx_bar = ev.idx;
                        samples_out[j] = p;
                    }
                }

                stab_bar.activate(ev.idx);
                range_bar.activate(ev.idx);
            }
        }
    }

    // （可选）这里可以检查是否所有 samples_out[j] 都被填上了合法值
    // for (std::size_t j = 0; j < t; ++j) {
    //     assert(samples_out[j].idx_c >= 0 && samples_out[j].idx_bar >= 0);
    // }
}
