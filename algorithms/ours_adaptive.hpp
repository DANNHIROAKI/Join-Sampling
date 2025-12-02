#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <utility>
#include <vector>

#include "geometry.hpp"
#include "structures.hpp"
#include "baseline1.hpp"

namespace rect_sampler {
namespace ours_adaptive {

namespace detail {

enum EventType : int {
    END_EVENT   = 0,
    START_EVENT = 1
};

struct Event {
    double x{0.0};
    EventType type{END_EVENT};
    bool is_c{false};
    int rect_idx{-1};
    std::size_t start_id{std::numeric_limits<std::size_t>::max()};
};

inline bool event_less(const Event& a, const Event& b) {
    if (a.x < b.x) return true;
    if (a.x > b.x) return false;
    if (a.type != b.type) return a.type < b.type;  // END before START
    if (a.type == START_EVENT && b.type == START_EVENT) {
        if (a.is_c != b.is_c) return a.is_c > b.is_c;  // Rc START before Rbar
    }
    return false;
}

struct AliasTable {
    std::vector<double> prob;
    std::vector<std::size_t> alias;

    void clear() {
        prob.clear();
        alias.clear();
    }

    void build(const std::vector<double>& weights) {
        clear();
        const std::size_t n = weights.size();
        if (n == 0) return;

        prob.resize(n);
        alias.resize(n);

        std::vector<double> scaled(n);
        double sum = 0.0;
        for (double w : weights) sum += w;
        if (sum <= 0.0) {
            for (std::size_t i = 0; i < n; ++i) {
                prob[i]  = 1.0;
                alias[i] = i;
            }
            return;
        }

        for (std::size_t i = 0; i < n; ++i) {
            scaled[i] = weights[i] * static_cast<double>(n) / sum;
        }

        std::vector<std::size_t> small;
        std::vector<std::size_t> large;
        small.reserve(n);
        large.reserve(n);

        for (std::size_t i = 0; i < n; ++i) {
            if (scaled[i] < 1.0) small.push_back(i);
            else large.push_back(i);
        }

        while (!small.empty() && !large.empty()) {
            std::size_t s = small.back();
            small.pop_back();
            std::size_t l = large.back();
            prob[s]  = scaled[s];
            alias[s] = l;
            scaled[l] = (scaled[l] + scaled[s]) - 1.0;
            if (scaled[l] < 1.0) {
                small.push_back(l);
                large.pop_back();
            }
        }

        while (!large.empty()) {
            std::size_t l = large.back();
            large.pop_back();
            prob[l]  = 1.0;
            alias[l] = l;
        }
        while (!small.empty()) {
            std::size_t s = small.back();
            small.pop_back();
            prob[s]  = 1.0;
            alias[s] = s;
        }
    }

    template <class URNG>
    std::size_t sample(URNG& rng) const {
        const std::size_t n = prob.size();
        if (n == 0) return std::numeric_limits<std::size_t>::max();
        std::uniform_int_distribution<std::size_t> col_dist(0, n - 1);
        std::uniform_real_distribution<double> u01(0.0, 1.0);
        std::size_t col = col_dist(rng);
        double u = u01(rng);
        return (u < prob[col]) ? col : alias[col];
    }
};

inline std::size_t safe_threshold(
    std::size_t n,
    std::size_t t,
    std::size_t max_pairs) {
    double nn = static_cast<double>(n);
    double logn = (n > 1) ? std::log(nn) : 0.0;
    double est = nn * logn + static_cast<double>(t);
    if (est < 1.0) est = 1.0;
    if (est > static_cast<double>(std::numeric_limits<std::size_t>::max())) {
        return max_pairs;
    }
    std::size_t th = static_cast<std::size_t>(est);
    if (th > max_pairs) th = max_pairs;
    return th;
}

inline std::size_t safe_pair_capacity(std::size_t n1, std::size_t n2, std::size_t max_pairs) {
    if (n1 == 0 || n2 == 0) return 0;
    if (n1 > max_pairs / n2) return max_pairs;
    return std::min<std::size_t>(n1 * n2, max_pairs);
}

} // namespace detail

// 自适应版本：当跨集合 pair 数较小时显式枚举并在数组上采样；
// 当 pair 数超过阈值时退化为 Ours 的双扫描计数采样。
template <class URNG>
inline void sample_pairs_adaptive(
    const RectList& Rc,
    const RectList& Rbar,
    std::size_t t,
    URNG& rng,
    PairList& samples_out,
    std::size_t pair_threshold = 0
) {
    using namespace detail;

    samples_out.clear();
    const std::size_t max_samples = samples_out.max_size();
    if (t == 0 || max_samples == 0) return;
    if (t > max_samples) t = max_samples;

    const std::size_t n1 = Rc.size();
    const std::size_t n2 = Rbar.size();
    if (n1 == 0 || n2 == 0) return;

    const std::size_t n = n1 + n2;
    const std::size_t max_pairs = samples_out.max_size();
    const std::size_t threshold = (pair_threshold == 0)
        ? safe_threshold(n, t, max_pairs)
        : std::min(pair_threshold, max_pairs);

    // ---------- 1) 准备 y 方向全集 ----------
    std::vector<double> endpoints_c;
    std::vector<double> endpoints_bar;
    std::vector<double> points_c;
    std::vector<double> points_bar;

    const std::size_t cap_end_c = safe_pair_capacity(2, n1, endpoints_c.max_size());
    endpoints_c.reserve(cap_end_c);
    const std::size_t cap_pts_c = std::min<std::size_t>(n1, points_c.max_size());
    points_c.reserve(cap_pts_c);
    for (std::size_t i = 0; i < n1; ++i) {
        endpoints_c.push_back(Rc[i].y_min);
        endpoints_c.push_back(Rc[i].y_max);
        points_c.push_back(Rc[i].y_min);
    }

    const std::size_t cap_end_bar = safe_pair_capacity(2, n2, endpoints_bar.max_size());
    endpoints_bar.reserve(cap_end_bar);
    const std::size_t cap_pts_bar = std::min<std::size_t>(n2, points_bar.max_size());
    points_bar.reserve(cap_pts_bar);
    for (std::size_t i = 0; i < n2; ++i) {
        endpoints_bar.push_back(Rbar[i].y_min);
        endpoints_bar.push_back(Rbar[i].y_max);
        points_bar.push_back(Rbar[i].y_min);
    }

    // ---------- 2) 构建事件 ----------
    std::vector<Event> events;
    const std::size_t total_rects = (n1 > std::numeric_limits<std::size_t>::max() - n2)
        ? std::numeric_limits<std::size_t>::max()
        : n1 + n2;
    const std::size_t event_reserve = safe_pair_capacity(2, total_rects, events.max_size());
    events.reserve(event_reserve);

    for (std::size_t i = 0; i < n1; ++i) {
        Event s, e;
        s.x = Rc[i].x_min;
        s.type = START_EVENT;
        s.is_c = true;
        s.rect_idx = static_cast<int>(i);

        e.x = Rc[i].x_max;
        e.type = END_EVENT;
        e.is_c = true;
        e.rect_idx = static_cast<int>(i);

        events.push_back(s);
        events.push_back(e);
    }

    for (std::size_t i = 0; i < n2; ++i) {
        Event s, e;
        s.x = Rbar[i].x_min;
        s.type = START_EVENT;
        s.is_c = false;
        s.rect_idx = static_cast<int>(i);

        e.x = Rbar[i].x_max;
        e.type = END_EVENT;
        e.is_c = false;
        e.rect_idx = static_cast<int>(i);

        events.push_back(s);
        events.push_back(e);
    }

    std::sort(events.begin(), events.end(), event_less);

    // ---------- 3) 首次扫描：自适应枚举/计数 ----------
    std::size_t num_starts = 0;
    for (const auto& ev : events) {
        if (ev.type == START_EVENT) ++num_starts;
    }
    if (num_starts == 0) return;

    std::vector<std::size_t> k1(num_starts, 0), k2(num_starts, 0), w(num_starts, 0);

    DynamicStabbing1D stab_c(endpoints_c);
    DynamicStabbing1D stab_bar(endpoints_bar);
    RangeTree1D range_c(points_c);
    RangeTree1D range_bar(points_bar);

    PairList enumerated_pairs;
    bool count_only = false;
    std::size_t total_w = 0;
    std::size_t M = 0;
    std::size_t next_start_id = 0;
    std::vector<int> temp_ids;

    if (!count_only) {
        const std::size_t reserve_pairs = std::min(safe_pair_capacity(n1, n2, enumerated_pairs.max_size()), threshold);
        enumerated_pairs.reserve(reserve_pairs);
    }

    for (auto& ev : events) {
        if (ev.type == END_EVENT) {
            if (ev.is_c) {
                const Rect& r = Rc[ev.rect_idx];
                stab_c.deactivate(ev.rect_idx, r.y_min, r.y_max);
                range_c.deactivate(ev.rect_idx);
            } else {
                const Rect& r = Rbar[ev.rect_idx];
                stab_bar.deactivate(ev.rect_idx, r.y_min, r.y_max);
                range_bar.deactivate(ev.rect_idx);
            }
            continue;
        }

        std::size_t sid = next_start_id++;
        ev.start_id = sid;

        const bool is_c = ev.is_c;
        const Rect& r = is_c ? Rc[ev.rect_idx] : Rbar[ev.rect_idx];
        double D = r.y_min;
        double U = r.y_max;

        if (!count_only) {
            // 枚举模式
            temp_ids.clear();
            if (is_c) {
                stab_bar.enumerate(D, temp_ids);
                std::size_t cnt1 = temp_ids.size();
                k1[sid] = cnt1;
                for (int idx_bar : temp_ids) {
                    if (M < threshold) {
                        enumerated_pairs.push_back({ev.rect_idx, idx_bar});
                    }
                    ++M;
                }

                temp_ids.clear();
                range_bar.enumerate(D, U, temp_ids);
                std::size_t cnt2 = temp_ids.size();
                k2[sid] = cnt2;
                for (int idx_bar : temp_ids) {
                    if (M < threshold) {
                        enumerated_pairs.push_back({ev.rect_idx, idx_bar});
                    }
                    ++M;
                }
            } else {
                stab_c.enumerate(D, temp_ids);
                std::size_t cnt1 = temp_ids.size();
                k1[sid] = cnt1;
                for (int idx_c : temp_ids) {
                    if (M < threshold) {
                        enumerated_pairs.push_back({idx_c, ev.rect_idx});
                    }
                    ++M;
                }

                temp_ids.clear();
                range_c.enumerate(D, U, temp_ids);
                std::size_t cnt2 = temp_ids.size();
                k2[sid] = cnt2;
                for (int idx_c : temp_ids) {
                    if (M < threshold) {
                        enumerated_pairs.push_back({idx_c, ev.rect_idx});
                    }
                    ++M;
                }
            }

            w[sid] = k1[sid] + k2[sid];
            total_w += w[sid];

            if (M > threshold) {
                count_only = true;
                PairList().swap(enumerated_pairs);  // 释放内存
            }
        } else {
            std::size_t cnt1 = is_c ? stab_bar.count(D) : stab_c.count(D);
            std::size_t cnt2 = is_c ? range_bar.count(D, U) : range_c.count(D, U);
            k1[sid] = cnt1;
            k2[sid] = cnt2;
            w[sid] = cnt1 + cnt2;
            total_w += w[sid];
        }

        // 激活当前矩形
        if (is_c) {
            stab_c.activate(ev.rect_idx, r.y_min, r.y_max);
            range_c.activate(ev.rect_idx);
        } else {
            stab_bar.activate(ev.rect_idx, r.y_min, r.y_max);
            range_bar.activate(ev.rect_idx);
        }
    }

    if (!count_only) {
        // 情况 A：枚举完毕，使用数组采样
        baseline1::sample_from_all_pairs(enumerated_pairs, t, rng, samples_out);
        return;
    }

    if (total_w == 0) return;  // 没有相交对

    // ---------- 4) Phase2：事件 alias + slot 规划 ----------
    std::vector<std::size_t> active_sids;
    std::vector<double> weights;
    active_sids.reserve(num_starts);
    weights.reserve(num_starts);
    for (std::size_t sid = 0; sid < num_starts; ++sid) {
        if (w[sid] > 0) {
            active_sids.push_back(sid);
            weights.push_back(static_cast<double>(w[sid]));
        }
    }
    if (active_sids.empty()) return;

    AliasTable alias;
    alias.build(weights);

    std::vector<std::vector<std::size_t>> slots1(num_starts);
    std::vector<std::vector<std::size_t>> slots2(num_starts);
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    for (std::size_t j = 0; j < t; ++j) {
        std::size_t col = alias.sample(rng);
        if (col == std::numeric_limits<std::size_t>::max()) continue;
        std::size_t sid = active_sids[col];
        double p1 = (w[sid] == 0)
            ? 0.0
            : static_cast<double>(k1[sid]) / static_cast<double>(w[sid]);
        double u = u01(rng);
        if (u < p1) {
            slots1[sid].push_back(j);
        } else {
            slots2[sid].push_back(j);
        }
    }

    samples_out.resize(t);
    for (std::size_t j = 0; j < t; ++j) {
        samples_out[j].idx_c = -1;
        samples_out[j].idx_bar = -1;
    }

    // ---------- 5) Phase3：第二次扫描 + 局部采样 ----------
    DynamicStabbing1D stab_c2(endpoints_c);
    DynamicStabbing1D stab_bar2(endpoints_bar);
    RangeTree1D range_c2(points_c);
    RangeTree1D range_bar2(points_bar);

    for (const auto& ev : events) {
        if (ev.type == END_EVENT) {
            if (ev.is_c) {
                const Rect& r = Rc[ev.rect_idx];
                stab_c2.deactivate(ev.rect_idx, r.y_min, r.y_max);
                range_c2.deactivate(ev.rect_idx);
            } else {
                const Rect& r = Rbar[ev.rect_idx];
                stab_bar2.deactivate(ev.rect_idx, r.y_min, r.y_max);
                range_bar2.deactivate(ev.rect_idx);
            }
            continue;
        }

        const std::size_t sid = ev.start_id;
        if (sid == std::numeric_limits<std::size_t>::max()) continue;

        const bool is_c = ev.is_c;
        const Rect& r = is_c ? Rc[ev.rect_idx] : Rbar[ev.rect_idx];
        double D = r.y_min;
        double U = r.y_max;

        if (is_c) {
            auto& S1 = slots1[sid];
            if (!S1.empty()) {
                std::vector<int> ids;
                stab_bar2.sample(D, S1.size(), rng, ids);
                for (std::size_t i = 0; i < S1.size() && i < ids.size(); ++i) {
                    samples_out[S1[i]].idx_c = ev.rect_idx;
                    samples_out[S1[i]].idx_bar = ids[i];
                }
                S1.clear();
            }

            auto& S2 = slots2[sid];
            if (!S2.empty()) {
                std::vector<int> ids;
                range_bar2.sample(D, U, S2.size(), rng, ids);
                for (std::size_t i = 0; i < S2.size() && i < ids.size(); ++i) {
                    samples_out[S2[i]].idx_c = ev.rect_idx;
                    samples_out[S2[i]].idx_bar = ids[i];
                }
                S2.clear();
            }

            stab_c2.activate(ev.rect_idx, r.y_min, r.y_max);
            range_c2.activate(ev.rect_idx);
        } else {
            auto& S1 = slots1[sid];
            if (!S1.empty()) {
                std::vector<int> ids;
                stab_c2.sample(D, S1.size(), rng, ids);
                for (std::size_t i = 0; i < S1.size() && i < ids.size(); ++i) {
                    samples_out[S1[i]].idx_c = ids[i];
                    samples_out[S1[i]].idx_bar = ev.rect_idx;
                }
                S1.clear();
            }

            auto& S2 = slots2[sid];
            if (!S2.empty()) {
                std::vector<int> ids;
                range_c2.sample(D, U, S2.size(), rng, ids);
                for (std::size_t i = 0; i < S2.size() && i < ids.size(); ++i) {
                    samples_out[S2[i]].idx_c = ids[i];
                    samples_out[S2[i]].idx_bar = ev.rect_idx;
                }
                S2.clear();
            }

            stab_bar2.activate(ev.rect_idx, r.y_min, r.y_max);
            range_bar2.activate(ev.rect_idx);
        }
    }
}

} // namespace ours_adaptive
} // namespace rect_sampler
