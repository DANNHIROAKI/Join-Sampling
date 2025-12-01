#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <cstddef>

#include "geometry.hpp"
#include "structures.hpp"
#include "utils.hpp"
#include "baseline1.hpp"   // 复用 baseline_sample_uniform

// ==============================================
// Baseline2：Plane Sweep + 1D 结构显式枚举 J
// ==============================================
//
// 结构与主算法同构，只是这里是“枚举所有相交对”。
// 接口：
//   baseline2_enumerate(Rc, Rbar, all_pairs_out);
//
// 采样直接复用 baseline_sample_uniform(...)。

namespace baseline2_detail {

    enum EventType : int {
        END_EVENT   = 0,
        START_EVENT = 1
    };

    struct Event {
        double x{0.0};
        EventType type{END_EVENT};
        bool is_c{false};   // true: Rc, false: Rbar
        int  idx{-1};       // 在对应集合中的索引
    };

    inline bool event_less(const Event& a, const Event& b) {
        if (a.x < b.x) return true;
        if (a.x > b.x) return false;
        // x 相同：END 在 START 前
        if (a.type != b.type) return a.type < b.type;
        // 若同为 START，则固定 R_c 在 R_{\bar c} 前（也可以反过来，只要固定）
        if (a.type == START_EVENT && b.type == START_EVENT) {
            if (a.is_c != b.is_c) return a.is_c > b.is_c; // true 在前
        }
        return false;
    }

} // namespace baseline2_detail

/// 使用 plane sweep + 1D 结构显式枚举所有跨集合相交对 J
inline void baseline2_enumerate(
    const std::vector<Rect>& Rc,
    const std::vector<Rect>& Rbar,
    std::vector<RectPair>& all_pairs_out
) {
    using namespace baseline2_detail;

    all_pairs_out.clear();

    const int n1 = static_cast<int>(Rc.size());
    const int n2 = static_cast<int>(Rbar.size());

    // ---------- 构建 1D 结构的全集 ----------

    // R_c: 区间 [D, U) = [y_min, y_max)，点 D = y_min
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

    // R_{\bar c}: 同理
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

    // ---------- 构建事件数组并排序 ----------

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

    // ---------- 扫描并显式枚举所有 pair ----------

    std::vector<int> tmp_ids;   // 复用的临时数组

    for (const auto& ev : events) {
        if (ev.type == END_EVENT) {
            // END：从自身集合的一维结构中移除
            if (ev.is_c) {
                // r ∈ R_c
                stab_c.deactivate(ev.idx);
                range_c.deactivate(ev.idx);
            } else {
                // r ∈ R_{\bar c}
                stab_bar.deactivate(ev.idx);
                range_bar.deactivate(ev.idx);
            }
        } else {
            // START_EVENT
            if (ev.is_c) {
                // r ∈ R_c，另一类是 R_{\bar c}
                const Rect& r = Rc[ev.idx];
                double D = r.y_min;
                double U = r.y_max;

                // K_e^{(1)}: D(r_bar) <= D < U(r_bar)
                tmp_ids.clear();
                stab_bar.enumerate(D, tmp_ids);
                for (int idx_bar : tmp_ids) {
                    RectPair p;
                    p.idx_c   = ev.idx;
                    p.idx_bar = idx_bar;
                    all_pairs_out.push_back(p);
                }

                // K_e^{(2)}: D < D(r_bar) < U
                tmp_ids.clear();
                range_bar.enumerate(D, U, tmp_ids); // (D, U) 开区间
                for (int idx_bar : tmp_ids) {
                    RectPair p;
                    p.idx_c   = ev.idx;
                    p.idx_bar = idx_bar;
                    all_pairs_out.push_back(p);
                }

                // 把 r 加入 R_c 的活跃集
                stab_c.activate(ev.idx);
                range_c.activate(ev.idx);
            } else {
                // r ∈ R_{\bar c}，对称处理
                const Rect& r = Rbar[ev.idx];
                double D = r.y_min;
                double U = r.y_max;

                // K_e^{(1)}: D(r_c) <= D < U(r_c)
                tmp_ids.clear();
                stab_c.enumerate(D, tmp_ids);
                for (int idx_c : tmp_ids) {
                    RectPair p;
                    p.idx_c   = idx_c;
                    p.idx_bar = ev.idx;
                    all_pairs_out.push_back(p);
                }

                // K_e^{(2)}: D < D(r_c) < U
                tmp_ids.clear();
                range_c.enumerate(D, U, tmp_ids);
                for (int idx_c : tmp_ids) {
                    RectPair p;
                    p.idx_c   = idx_c;
                    p.idx_bar = ev.idx;
                    all_pairs_out.push_back(p);
                }

                // 把 r 加入 R_{\bar c} 的活跃集
                stab_bar.activate(ev.idx);
                range_bar.activate(ev.idx);
            }
        }
    }
}

