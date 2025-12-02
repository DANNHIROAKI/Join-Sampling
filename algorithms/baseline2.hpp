#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <cstddef>

#include "geometry.hpp"
#include "structures.hpp"
#include "baseline1.hpp"

namespace rect_sampler {
namespace baseline2 {

namespace detail {

enum EventType : int {
    END_EVENT   = 0,
    START_EVENT = 1
};

struct Event {
    double    x{0.0};
    EventType type{END_EVENT};
    bool      is_c{false};
    int       idx{-1};   // 矩形在对应集合中的索引
};

inline bool event_less(const Event& a, const Event& b) {
    if (a.x < b.x) return true;
    if (a.x > b.x) return false;
    // x 相同：END 在 START 前
    if (a.type != b.type) return a.type < b.type;
    // 同为 START：固定 R_c 在 R_{\bar c} 前
    if (a.type == START_EVENT && b.type == START_EVENT) {
        if (a.is_c != b.is_c) return a.is_c > b.is_c;
    }
    return false;
}

} // namespace detail

// plane sweep + 1D 结构显式枚举所有相交对
inline void enumerate_all_pairs(
    const RectList& Rc,
    const RectList& Rbar,
    PairList& all_pairs_out
) {
    using namespace detail;

    all_pairs_out.clear();

    const std::size_t n1 = Rc.size();
    const std::size_t n2 = Rbar.size();

    if (n1 > 0 && n2 > 0) {
        const std::size_t max_reservable = all_pairs_out.max_size();
        // Safely compute n1 * n2 without overflow; fall back to max_size if it would wrap.
        if (n1 > max_reservable / n2) {
            all_pairs_out.reserve(max_reservable);
        } else {
            all_pairs_out.reserve(std::min<std::size_t>(n1 * n2, max_reservable));
        }
    }

    // ---------- 构建 1D 结构的全集 ----------

    // R_c
    std::vector<double> endpoints_c;
    endpoints_c.reserve(2 * n1);
    std::vector<double> points_c;
    points_c.reserve(n1);

    for (std::size_t i = 0; i < n1; ++i) {
        endpoints_c.push_back(Rc[i].y_min);
        endpoints_c.push_back(Rc[i].y_max);
        points_c.push_back(Rc[i].y_min);
    }

    // R_{\bar c}
    std::vector<double> endpoints_bar;
    endpoints_bar.reserve(2 * n2);
    std::vector<double> points_bar;
    points_bar.reserve(n2);

    for (std::size_t i = 0; i < n2; ++i) {
        endpoints_bar.push_back(Rbar[i].y_min);
        endpoints_bar.push_back(Rbar[i].y_max);
        points_bar.push_back(Rbar[i].y_min);
    }

    DynamicStabbing1D stab_c(endpoints_c);
    DynamicStabbing1D stab_bar(endpoints_bar);
    RangeTree1D      range_c(points_c);
    RangeTree1D      range_bar(points_bar);

    // ---------- 构建事件数组并排序 ----------

    std::vector<Event> events;
    events.reserve((n1 + n2) * 2);

    for (std::size_t i = 0; i < n1; ++i) {
        Event s,e;
        s.x    = Rc[i].x_min;
        s.type = START_EVENT;
        s.is_c = true;
        s.idx  = static_cast<int>(i);

        e.x    = Rc[i].x_max;
        e.type = END_EVENT;
        e.is_c = true;
        e.idx  = static_cast<int>(i);

        events.push_back(s);
        events.push_back(e);
    }

    for (std::size_t i = 0; i < n2; ++i) {
        Event s,e;
        s.x    = Rbar[i].x_min;
        s.type = START_EVENT;
        s.is_c = false;
        s.idx  = static_cast<int>(i);

        e.x    = Rbar[i].x_max;
        e.type = END_EVENT;
        e.is_c = false;
        e.idx  = static_cast<int>(i);

        events.push_back(s);
        events.push_back(e);
    }

    std::sort(events.begin(), events.end(), event_less);

    // ---------- 扫描并显式枚举所有 pair ----------

    std::vector<int> tmp_ids;

    for (const auto& ev : events) {
        if (ev.type == END_EVENT) {
            if (ev.is_c) {
                const Rect& r = Rc[ev.idx];
                stab_c.deactivate(ev.idx, r.y_min, r.y_max);
                range_c.deactivate(ev.idx);
            } else {
                const Rect& r = Rbar[ev.idx];
                stab_bar.deactivate(ev.idx, r.y_min, r.y_max);
                range_bar.deactivate(ev.idx);
            }
        } else {
            if (ev.is_c) {
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
                range_bar.enumerate(D, U, tmp_ids); // (D,U)
                for (int idx_bar : tmp_ids) {
                    RectPair p;
                    p.idx_c   = ev.idx;
                    p.idx_bar = idx_bar;
                    all_pairs_out.push_back(p);
                }

                // 激活 r ∈ R_c
                stab_c.activate(ev.idx, r.y_min, r.y_max);
                range_c.activate(ev.idx);
            } else {
                const Rect& r = Rbar[ev.idx];
                double D = r.y_min;
                double U = r.y_max;

                tmp_ids.clear();
                stab_c.enumerate(D, tmp_ids);
                for (int idx_c : tmp_ids) {
                    RectPair p;
                    p.idx_c   = idx_c;
                    p.idx_bar = ev.idx;
                    all_pairs_out.push_back(p);
                }

                tmp_ids.clear();
                range_c.enumerate(D, U, tmp_ids);
                for (int idx_c : tmp_ids) {
                    RectPair p;
                    p.idx_c   = idx_c;
                    p.idx_bar = ev.idx;
                    all_pairs_out.push_back(p);
                }

                stab_bar.activate(ev.idx, r.y_min, r.y_max);
                range_bar.activate(ev.idx);
            }
        }
    }
}

// 采样直接复用 baseline1 的 sample_from_all_pairs
template <class URNG>
inline void sample_from_all_pairs(
    const PairList& all_pairs,
    std::size_t t,
    URNG& rng,
    PairList& samples_out
) {
    baseline1::sample_from_all_pairs(all_pairs, t, rng, samples_out);
}

} // namespace baseline2
} // namespace rect_sampler
