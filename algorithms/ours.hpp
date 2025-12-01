#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <cstddef>

#include "geometry.hpp"
#include "structures.hpp"

namespace rect_sampler {
namespace ours {

namespace detail {

enum EventType : int {
    END_EVENT   = 0,
    START_EVENT = 1
};

struct Event {
    double    x{0.0};
    EventType type{END_EVENT};
    bool      is_c{false};   // true: from Rc, false: from Rbar
    int       rect_idx{-1};  // index in Rc or Rbar
    int       start_id{-1};  // index of START event in [0, num_starts)
};

inline bool event_less(const Event& a, const Event& b) {
    if (a.x < b.x) return true;
    if (a.x > b.x) return false;
    // same x: END before START
    if (a.type != b.type) return a.type < b.type;
    // both START: Rc before Rbar (is_c = true first)
    if (a.type == START_EVENT && b.type == START_EVENT) {
        if (a.is_c != b.is_c) return a.is_c > b.is_c;
    }
    return false;
}

// 简单 alias 表，用于在 START 事件上按权重 w_e 采样
struct AliasTable {
    std::vector<double> prob;   // 列的“保留”概率
    std::vector<int>    alias;  // 列的备选索引

    void clear() {
        prob.clear();
        alias.clear();
    }

    // weights[i] >= 0，且至少有一个 > 0
    void build(const std::vector<double>& weights) {
        clear();
        const int n = static_cast<int>(weights.size());
        if (n == 0) return;

        prob.resize(n);
        alias.resize(n);

        std::vector<double> scaled(n);
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += weights[i];
        }
        if (sum <= 0.0) {
            // 退化情况：全部权重为 0，直接设成均匀（理论上不会走到这）
            for (int i = 0; i < n; ++i) {
                prob[i]  = 1.0;
                alias[i] = i;
            }
            return;
        }

        for (int i = 0; i < n; ++i) {
            scaled[i] = weights[i] * n / sum;
        }

        std::vector<int> small, large;
        small.reserve(n);
        large.reserve(n);
        for (int i = 0; i < n; ++i) {
            if (scaled[i] < 1.0) small.push_back(i);
            else                 large.push_back(i);
        }

        while (!small.empty() && !large.empty()) {
            int s = small.back(); small.pop_back();
            int l = large.back();
            prob[s]  = scaled[s];
            alias[s] = l;
            scaled[l] = (scaled[l] + scaled[s]) - 1.0;
            if (scaled[l] < 1.0) {
                small.push_back(l);
                large.pop_back();
            }
        }
        while (!large.empty()) {
            int l = large.back(); large.pop_back();
            prob[l]  = 1.0;
            alias[l] = l;
        }
        while (!small.empty()) {
            int s = small.back(); small.pop_back();
            prob[s]  = 1.0;
            alias[s] = s;
        }
    }

    template <class URNG>
    int sample(URNG& rng) const {
        const int n = static_cast<int>(prob.size());
        if (n == 0) return -1;
        std::uniform_int_distribution<int> col_dist(0, n - 1);
        std::uniform_real_distribution<double> u01(0.0, 1.0);
        int col = col_dist(rng);
        double u = u01(rng);
        return (u < prob[col]) ? col : alias[col];
    }
};

} // namespace detail

// ======================== 主接口 ========================
//
// 从全体跨集合相交对 J 中 i.i.d. 均匀采样 t 个样本。
// 采样结果存到 samples_out（RectPair 里是索引，而不是几何数据）。
//
template <class URNG>
inline void sample_pairs(
    const RectList& Rc,
    const RectList& Rbar,
    std::size_t t,
    URNG& rng,
    PairList& samples_out
) {
    using namespace detail;

    samples_out.clear();
    if (t == 0) return;

    const int n1 = static_cast<int>(Rc.size());
    const int n2 = static_cast<int>(Rbar.size());
    if (n1 == 0 || n2 == 0) {
        // 没有任何跨集合相交对
        return;
    }

    // ---------- 1. 准备 y 方向的全集（端点 & 点） ----------
    std::vector<double> endpoints_c;
    std::vector<double> endpoints_bar;
    std::vector<double> points_c;
    std::vector<double> points_bar;

    endpoints_c.reserve(2 * n1);
    points_c.reserve(n1);
    for (int i = 0; i < n1; ++i) {
        endpoints_c.push_back(Rc[i].y_min);
        endpoints_c.push_back(Rc[i].y_max);
        points_c.push_back(Rc[i].y_min);
    }

    endpoints_bar.reserve(2 * n2);
    points_bar.reserve(n2);
    for (int i = 0; i < n2; ++i) {
        endpoints_bar.push_back(Rbar[i].y_min);
        endpoints_bar.push_back(Rbar[i].y_max);
        points_bar.push_back(Rbar[i].y_min);
    }

    // ---------- 2. 构建事件数组并排序 ----------
    std::vector<Event> events;
    events.reserve((n1 + n2) * 2);

    // Rc 的 START / END
    for (int i = 0; i < n1; ++i) {
        Event s,e;
        s.x        = Rc[i].x_min;
        s.type     = START_EVENT;
        s.is_c     = true;
        s.rect_idx = i;
        s.start_id = -1;

        e.x        = Rc[i].x_max;
        e.type     = END_EVENT;
        e.is_c     = true;
        e.rect_idx = i;
        e.start_id = -1;

        events.push_back(s);
        events.push_back(e);
    }

    // Rbar 的 START / END
    for (int i = 0; i < n2; ++i) {
        Event s,e;
        s.x        = Rbar[i].x_min;
        s.type     = START_EVENT;
        s.is_c     = false;
        s.rect_idx = i;
        s.start_id = -1;

        e.x        = Rbar[i].x_max;
        e.type     = END_EVENT;
        e.is_c     = false;
        e.rect_idx = i;
        e.start_id = -1;

        events.push_back(s);
        events.push_back(e);
    }

    std::sort(events.begin(), events.end(), event_less);

    // ---------- 3. 第一遍扫描：只做计数，得到每个 START 事件的 w_e, k_e^{(1)}, k_e^{(2)} ----------

    // 统计 START 事件总数
    int num_starts = 0;
    for (const auto& ev : events) {
        if (ev.type == START_EVENT) ++num_starts;
    }
    if (num_starts == 0) {
        // 没有矩形变成 active，就一定没有相交对
        return;
    }

    std::vector<std::size_t> k1(num_starts, 0), k2(num_starts, 0), w(num_starts, 0);

    DynamicStabbing1D stab_c(endpoints_c);
    DynamicStabbing1D stab_bar(endpoints_bar);
    RangeTree1D      range_c(points_c);
    RangeTree1D      range_bar(points_bar);

    int next_start_id = 0;
    std::size_t total_w = 0;

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
        } else { // START_EVENT
            int sid = next_start_id++;
            ev.start_id = sid;

            if (ev.is_c) {
                const Rect& r = Rc[ev.rect_idx];
                double D = r.y_min;
                double U = r.y_max;

                std::size_t cnt1 = stab_bar.count(D);      // K_e^{(1)}
                std::size_t cnt2 = range_bar.count(D, U);  // K_e^{(2)}，(D,U)

                k1[sid] = cnt1;
                k2[sid] = cnt2;
                w[sid]  = cnt1 + cnt2;
                total_w += w[sid];

                // 把 r 插入 R_c 的结构
                stab_c.activate(ev.rect_idx, r.y_min, r.y_max);
                range_c.activate(ev.rect_idx);
            } else {
                const Rect& r = Rbar[ev.rect_idx];
                double D = r.y_min;
                double U = r.y_max;

                std::size_t cnt1 = stab_c.count(D);
                std::size_t cnt2 = range_c.count(D, U);

                k1[sid] = cnt1;
                k2[sid] = cnt2;
                w[sid]  = cnt1 + cnt2;
                total_w += w[sid];

                stab_bar.activate(ev.rect_idx, r.y_min, r.y_max);
                range_bar.activate(ev.rect_idx);
            }
        }
    }

    if (total_w == 0) {
        // 没有相交对，直接返回空
        return;
    }

    // ---------- 4. 第二阶段：采样规划（事件级 alias + 类型分配） ----------

    // 只在 w_e > 0 的事件上做 alias
    std::vector<int>    active_sids;
    std::vector<double> weights;
    active_sids.reserve(num_starts);
    weights.reserve(num_starts);

    for (int sid = 0; sid < num_starts; ++sid) {
        if (w[sid] > 0) {
            active_sids.push_back(sid);
            weights.push_back(static_cast<double>(w[sid]));
        }
    }
    if (active_sids.empty()) {
        // 理论上不会发生，因为 total_w > 0
        return;
    }

    detail::AliasTable alias;
    alias.build(weights);

    // 每个 START 事件 / 类型，对应要填的样本槽位
    std::vector<std::vector<std::size_t>> slots1(num_starts);
    std::vector<std::vector<std::size_t>> slots2(num_starts);

    std::uniform_real_distribution<double> u01(0.0, 1.0);

    for (std::size_t j = 0; j < t; ++j) {
        int col = alias.sample(rng);
        if (col < 0) continue; // 防御式，正常不会触发

        int sid = active_sids[col];

        double p1 = (w[sid] == 0)
            ? 0.0
            : static_cast<double>(k1[sid]) / static_cast<double>(w[sid]);

        double u = u01(rng);
        if (u < p1) {
            // 类型 1：K_e^{(1)}
            slots1[sid].push_back(j);
        } else {
            // 类型 2：K_e^{(2)}
            slots2[sid].push_back(j);
        }
    }

    // 预分配输出数组
    samples_out.resize(t);
    for (std::size_t j = 0; j < t; ++j) {
        samples_out[j].idx_c   = -1;
        samples_out[j].idx_bar = -1;
    }

    // ---------- 5. 第三阶段：第二次扫描 + 局部采样 + 回填 ----------

    DynamicStabbing1D stab_c2(endpoints_c);
    DynamicStabbing1D stab_bar2(endpoints_bar);
    RangeTree1D      range_c2(points_c);
    RangeTree1D      range_bar2(points_bar);

    std::vector<int> local_ids;

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
        } else { // START_EVENT
            int sid = ev.start_id;
            if (sid < 0) continue; // 理论上不会

            if (ev.is_c) {
                const Rect& r = Rc[ev.rect_idx];
                double D = r.y_min;
                double U = r.y_max;

                // 类型 1：从 K_e^{(1)} 中采样（stab_bar2）
                auto& S1 = slots1[sid];
                if (!S1.empty()) {
                    local_ids.clear();
                    stab_bar2.sample(D, S1.size(), rng, local_ids);
                    for (std::size_t i = 0; i < S1.size() && i < local_ids.size(); ++i) {
                        std::size_t slot = S1[i];
                        int idx_bar = local_ids[i];
                        samples_out[slot].idx_c   = ev.rect_idx;
                        samples_out[slot].idx_bar = idx_bar;
                    }
                    S1.clear();
                }

                // 类型 2：从 K_e^{(2)} 中采样（range_bar2）
                auto& S2 = slots2[sid];
                if (!S2.empty()) {
                    local_ids.clear();
                    range_bar2.sample(D, U, S2.size(), rng, local_ids);
                    for (std::size_t i = 0; i < S2.size() && i < local_ids.size(); ++i) {
                        std::size_t slot = S2[i];
                        int idx_bar = local_ids[i];
                        samples_out[slot].idx_c   = ev.rect_idx;
                        samples_out[slot].idx_bar = idx_bar;
                    }
                    S2.clear();
                }

                // 激活 r ∈ R_c
                stab_c2.activate(ev.rect_idx, r.y_min, r.y_max);
                range_c2.activate(ev.rect_idx);

            } else { // r ∈ Rbar
                const Rect& r = Rbar[ev.rect_idx];
                double D = r.y_min;
                double U = r.y_max;

                // 类型 1：K_e^{(1)}，stabbing 在 R_c 上
                auto& S1 = slots1[sid];
                if (!S1.empty()) {
                    local_ids.clear();
                    stab_c2.sample(D, S1.size(), rng, local_ids);
                    for (std::size_t i = 0; i < S1.size() && i < local_ids.size(); ++i) {
                        std::size_t slot = S1[i];
                        int idx_c = local_ids[i];
                        samples_out[slot].idx_c   = idx_c;
                        samples_out[slot].idx_bar = ev.rect_idx;
                    }
                    S1.clear();
                }

                // 类型 2：K_e^{(2)}，range 在 R_c 上
                auto& S2 = slots2[sid];
                if (!S2.empty()) {
                    local_ids.clear();
                    range_c2.sample(D, U, S2.size(), rng, local_ids);
                    for (std::size_t i = 0; i < S2.size() && i < local_ids.size(); ++i) {
                        std::size_t slot = S2[i];
                        int idx_c = local_ids[i];
                        samples_out[slot].idx_c   = idx_c;
                        samples_out[slot].idx_bar = ev.rect_idx;
                    }
                    S2.clear();
                }

                // 激活 r ∈ R_{\bar c}
                stab_bar2.activate(ev.rect_idx, r.y_min, r.y_max);
                range_bar2.activate(ev.rect_idx);
            }
        }
    }

    // 理论上到这里，samples_out[0..t-1] 全部都应该是合法 pair（如果 total_w>0）
}

} // namespace ours
} // namespace rect_sampler
