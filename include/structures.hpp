#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <cstdint>
#include <stdexcept>
#include <limits>

#include "geometry.hpp"

namespace rect_sampler {

// ======================== 动态 stabbing 结构 ========================
// 用于：维护区间 [D(r), U(r)) 的集合，支持：
//  - activate(id, y_low, y_high)
//  - deactivate(id, y_low, y_high)
//  - count(q)
//  - enumerate(q, out_ids)
//  - sample(q, t, rng, out_ids)
// 其中 id 通常对应 “矩形在数组中的索引”。

class DynamicStabbing1D {
public:
    DynamicStabbing1D() = default;

    // endpoints: 全部可能的 interval 端点 (D 或 U)，构造时会去重排序
    explicit DynamicStabbing1D(const std::vector<double>& endpoints) {
        init(endpoints);
    }

    void init(const std::vector<double>& endpoints) {
        coords_ = endpoints;
        std::sort(coords_.begin(), coords_.end());
        coords_.erase(std::unique(coords_.begin(), coords_.end()), coords_.end());

        if (coords_.size() < 2) {
            seg_count_ = 0;
            tree_.clear();
            return;
        }

        seg_count_ = static_cast<int>(coords_.size()) - 1;
        tree_.assign(4 * seg_count_ + 4, Node{});
    }

    bool ready() const {
        return seg_count_ > 0;
    }

    // 插入 / 删除一个区间 [y_low, y_high)
    // 要求 y_low, y_high 都在 endpoints 中（与构造时一致）。
    void activate(int interval_id, double y_low, double y_high) {
        if (!ready()) return;
        int L = coord_to_seg_index(y_low);
        int R = coord_to_seg_index(y_high);
        if (L < 0 || R < 0 || L >= R) return; // 空区间
        activate_internal(1, 0, seg_count_, L, R, interval_id);
    }

    void deactivate(int interval_id, double y_low, double y_high) {
        if (!ready()) return;
        int L = coord_to_seg_index(y_low);
        int R = coord_to_seg_index(y_high);
        if (L < 0 || R < 0 || L >= R) return;
        deactivate_internal(1, 0, seg_count_, L, R, interval_id);
    }

    // 计数：q 被多少个活跃区间刺中
    int count(double q) const {
        if (!ready()) return 0;
        if (q < coords_.front() || q >= coords_.back()) return 0;

        int leaf_seg = seg_index_for_point(q);
        if (leaf_seg < 0) return 0;

        std::vector<int> path;
        path.reserve(32);
        collect_path(1, 0, seg_count_, leaf_seg, path);

        int total = 0;
        for (int node_idx : path) {
            total += static_cast<int>(tree_[node_idx].active.size());
        }
        return total;
    }

    // 枚举：列出所有包含 q 的活跃 interval_id
    void enumerate(double q, std::vector<int>& out_ids) const {
        out_ids.clear();
        if (!ready()) return;
        if (q < coords_.front() || q >= coords_.back()) return;

        int leaf_seg = seg_index_for_point(q);
        if (leaf_seg < 0) return;

        std::vector<int> path;
        path.reserve(32);
        collect_path(1, 0, seg_count_, leaf_seg, path);

        for (int node_idx : path) {
            const Node& nd = tree_[node_idx];
            out_ids.insert(out_ids.end(), nd.active.begin(), nd.active.end());
        }
    }

    // 采样：从所有包含 q 的活跃 interval_id 中，独立均匀采样 sample_count 次（有放回）
    template <typename URNG>
    void sample(double q, int sample_count, URNG& rng,
                std::vector<int>& out_ids) const {
        out_ids.clear();
        if (!ready() || sample_count <= 0) return;
        if (q < coords_.front() || q >= coords_.back()) return;

        int leaf_seg = seg_index_for_point(q);
        if (leaf_seg < 0) return;

        // 收集路径上的节点
        std::vector<int> path;
        path.reserve(32);
        collect_path(1, 0, seg_count_, leaf_seg, path);

        // 只保留 non-empty 节点
        std::vector<int> non_empty_nodes;
        non_empty_nodes.reserve(path.size());
        for (int node_idx : path) {
            if (!tree_[node_idx].active.empty()) {
                non_empty_nodes.push_back(node_idx);
            }
        }
        if (non_empty_nodes.empty()) return;

        const int m = static_cast<int>(non_empty_nodes.size());

        // 构建 alias 表：每个节点权重 = active.size()
        std::vector<double> prob(m);
        std::vector<int> alias(m);
        build_alias_for_nodes(non_empty_nodes, prob, alias);

        std::uniform_int_distribution<int> column_dist(0, m - 1);
        std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

        out_ids.reserve(sample_count);
        for (int s = 0; s < sample_count; ++s) {
            int col = column_dist(rng);
            double p = prob_dist(rng);
            int bucket_idx = (p < prob[col]) ? col : alias[col];
            int node_idx = non_empty_nodes[bucket_idx];

            const Node& nd = tree_[node_idx];
            if (nd.active.empty()) continue; // 理论上不会

            std::uniform_int_distribution<int> item_dist(
                0, static_cast<int>(nd.active.size()) - 1
            );
            int pos = item_dist(rng);
            out_ids.push_back(nd.active[pos]);
        }
    }

    const std::vector<double>& coords() const { return coords_; }

private:
    struct Node {
        std::vector<int> active;                   // 当前活跃 interval_id
        std::unordered_map<int, int> pos;         // interval_id -> active 中的位置
    };

    std::vector<double> coords_;  // 端点坐标（升序去重）
    int seg_count_{0};            // 叶子 segment 数 = coords_.size() - 1
    std::vector<Node> tree_;      // 1-based segment tree

    // 将某个端点映射到 segment index：我们将 y_low, y_high 对应到 coords_ 中的下标
    // 并使用区间 [L, R) 的 segment [L, R)。
    int coord_to_seg_index(double y) const {
        auto it = std::lower_bound(coords_.begin(), coords_.end(), y);
        if (it == coords_.end()) return -1;
        int idx = static_cast<int>(it - coords_.begin());
        // segment index 的有效范围是 [0, seg_count_)
        if (idx < 0 || idx >= static_cast<int>(coords_.size())) return -1;
        return idx;
    }

    // 给定点 q，找到它所在的 segment index k，使得 q ∈ [coords_[k], coords_[k+1))
    int seg_index_for_point(double q) const {
        if (q < coords_.front() || q >= coords_.back()) return -1;
        auto it = std::upper_bound(coords_.begin(), coords_.end(), q);
        int idx = static_cast<int>(it - coords_.begin()) - 1;
        if (idx < 0 || idx >= seg_count_) return -1;
        return idx;
    }

    void activate_internal(int node, int nl, int nr,
                           int L, int R, int interval_id) {
        if (nr <= L || R <= nl) return; // 无交集
        if (L <= nl && nr <= R) {
            Node& nd = tree_[node];
            auto it = nd.pos.find(interval_id);
            if (it == nd.pos.end()) {
                int pos = static_cast<int>(nd.active.size());
                nd.active.push_back(interval_id);
                nd.pos.emplace(interval_id, pos);
            }
            return;
        }
        int mid = (nl + nr) / 2;
        activate_internal(node * 2, nl, mid, L, R, interval_id);
        activate_internal(node * 2 + 1, mid, nr, L, R, interval_id);
    }

    void deactivate_internal(int node, int nl, int nr,
                             int L, int R, int interval_id) {
        if (nr <= L || R <= nl) return;
        if (L <= nl && nr <= R) {
            Node& nd = tree_[node];
            auto it = nd.pos.find(interval_id);
            if (it != nd.pos.end()) {
                int idx = it->second;
                int last_id = nd.active.back();
                nd.active[idx] = last_id;
                nd.pos[last_id] = idx;
                nd.active.pop_back();
                nd.pos.erase(it);
            }
            return;
        }
        int mid = (nl + nr) / 2;
        deactivate_internal(node * 2, nl, mid, L, R, interval_id);
        deactivate_internal(node * 2 + 1, mid, nr, L, R, interval_id);
    }

    // 收集 root -> 叶子的路径
    void collect_path(int node, int nl, int nr, int leaf_seg,
                      std::vector<int>& path) const {
        if (leaf_seg < nl || leaf_seg >= nr) return;
        path.push_back(node);
        if (nl + 1 == nr) return; // 到叶子
        int mid = (nl + nr) / 2;
        if (leaf_seg < mid) {
            collect_path(node * 2, nl, mid, leaf_seg, path);
        } else {
            collect_path(node * 2 + 1, mid, nr, leaf_seg, path);
        }
    }

    void build_alias_for_nodes(const std::vector<int>& nodes,
                               std::vector<double>& prob,
                               std::vector<int>& alias) const {
        const int m = static_cast<int>(nodes.size());
        std::vector<double> w(m);
        double sum = 0.0;
        for (int i = 0; i < m; ++i) {
            int node_idx = nodes[i];
            double wi = static_cast<double>(tree_[node_idx].active.size());
            w[i] = wi;
            sum += wi;
        }
        if (sum <= 0.0) {
            // 不太可能出现（因为 nodes 已经过滤掉了 empty），但防一手
            for (int i = 0; i < m; ++i) {
                prob[i] = 1.0;
                alias[i] = i;
            }
            return;
        }

        std::vector<double> scaled(m);
        for (int i = 0; i < m; ++i) {
            scaled[i] = (w[i] * m) / sum;
        }

        std::vector<int> small, large;
        small.reserve(m);
        large.reserve(m);

        for (int i = 0; i < m; ++i) {
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

        // 剩余的都设为 1
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
};

// ======================== 动态 range 结构 ========================
// 维护点集合 { (value_i, owner_i) }，这里 owner_i 通常是矩形索引。
// 构造时给定所有点的 value，owner 默认 = 0..N-1（按输入顺序）。
// 支持：
//  - activate(owner_id)
//  - deactivate(owner_id)
//  - count(ell, r)        // 开区间 (ell, r) 或按你需要调整成 [ell, r)
//  - enumerate(ell, r, out_owner_ids)
//  - sample(ell, r, t, rng, out_owner_ids)

class RangeTree1D {
public:
    RangeTree1D() = default;

    // values[i] = 点 i 的坐标（例如 D(r_i)），owner_id = i
    explicit RangeTree1D(const std::vector<double>& values) {
        init(values);
    }

    void init(const std::vector<double>& values) {
        const int N = static_cast<int>(values.size());
        N_ = N;
        if (N_ <= 0) {
            xs_.clear();
            index_of_sorted_.clear();
            rank_of_owner_.clear();
            tree_.clear();
            return;
        }

        // 构造 (value, owner) 对
        std::vector<std::pair<double,int>> tmp;
        tmp.reserve(N_);
        for (int i = 0; i < N_; ++i) {
            tmp.emplace_back(values[i], i);
        }
        std::sort(tmp.begin(), tmp.end(),
                  [](const auto& a, const auto& b) {
                      if (a.first < b.first) return true;
                      if (a.first > b.first) return false;
                      return a.second < b.second;
                  });

        xs_.resize(N_);
        index_of_sorted_.resize(N_);
        rank_of_owner_.assign(N_, -1);

        for (int k = 0; k < N_; ++k) {
            xs_[k] = tmp[k].first;
            index_of_sorted_[k] = tmp[k].second;
            rank_of_owner_[tmp[k].second] = k;
        }

        tree_.assign(4 * N_ + 4, Node{});
    }

    bool ready() const { return N_ > 0; }

    // 激活 / 失活 一个 ownerId
    void activate(int owner_id) {
        if (!ready()) return;
        int idx = rank_of_owner_.at(owner_id);
        if (idx < 0 || idx >= N_) return;
        activate_internal(1, 0, N_, idx, owner_id);
    }

    void deactivate(int owner_id) {
        if (!ready()) return;
        int idx = rank_of_owner_.at(owner_id);
        if (idx < 0 || idx >= N_) return;
        deactivate_internal(1, 0, N_, idx, owner_id);
    }

    // 计数：值在 (ell, r) 内的活跃点数（按你需求可以换成 [ell, r)）
    int count(double ell, double r) const {
        if (!ready()) return 0;
        if (ell >= r) return 0;

        int L = static_cast<int>(std::upper_bound(xs_.begin(), xs_.end(), ell) - xs_.begin());
        int R = static_cast<int>(std::lower_bound(xs_.begin(), xs_.end(), r) - xs_.begin());
        if (L >= R) return 0;

        return count_internal(1, 0, N_, L, R);
    }

    // 枚举：值在 (ell, r) 内的所有活跃 owner_id
    void enumerate(double ell, double r, std::vector<int>& out_owner_ids) const {
        out_owner_ids.clear();
        if (!ready()) return;
        if (ell >= r) return;

        int L = static_cast<int>(std::upper_bound(xs_.begin(), xs_.end(), ell) - xs_.begin());
        int R = static_cast<int>(std::lower_bound(xs_.begin(), xs_.end(), r) - xs_.begin());
        if (L >= R) return;

        enumerate_internal(1, 0, N_, L, R, out_owner_ids);
    }

    // 采样：从值在 (ell, r) 内的活跃点集合中，独立均匀采样 sample_count 次（有放回）
    template <typename URNG>
    void sample(double ell, double r, int sample_count,
                URNG& rng, std::vector<int>& out_owner_ids) const {
        out_owner_ids.clear();
        if (!ready() || sample_count <= 0) return;
        if (ell >= r) return;

        int L = static_cast<int>(std::upper_bound(xs_.begin(), xs_.end(), ell) - xs_.begin());
        int R = static_cast<int>(std::lower_bound(xs_.begin(), xs_.end(), r) - xs_.begin());
        if (L >= R) return;

        // 先收集 canonical cover 节点
        std::vector<int> nodes;
        nodes.reserve(32);
        collect_nodes(1, 0, N_, L, R, nodes);

        // 过滤掉空节点
        std::vector<int> non_empty_nodes;
        non_empty_nodes.reserve(nodes.size());
        for (int node_idx : nodes) {
            if (!tree_[node_idx].active.empty()) {
                non_empty_nodes.push_back(node_idx);
            }
        }
        if (non_empty_nodes.empty()) return;

        const int m = static_cast<int>(non_empty_nodes.size());
        std::vector<double> prob(m);
        std::vector<int> alias(m);
        build_alias_for_nodes(non_empty_nodes, prob, alias);

        std::uniform_int_distribution<int> column_dist(0, m - 1);
        std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

        out_owner_ids.reserve(sample_count);
        for (int s = 0; s < sample_count; ++s) {
            int col = column_dist(rng);
            double p = prob_dist(rng);
            int bucket_idx = (p < prob[col]) ? col : alias[col];
            int node_idx = non_empty_nodes[bucket_idx];

            const Node& nd = tree_[node_idx];
            if (nd.active.empty()) continue; // 理论上不会

            std::uniform_int_distribution<int> item_dist(
                0, static_cast<int>(nd.active.size()) - 1
            );
            int pos = item_dist(rng);
            out_owner_ids.push_back(nd.active[pos]);
        }
    }

    const std::vector<double>& xs() const { return xs_; }

private:
    struct Node {
        std::vector<int> active;                   // 当前活跃 owner_id
        std::unordered_map<int, int> pos;         // owner_id -> active 中的位置
    };

    int N_{0};
    std::vector<double> xs_;           // 排序后的值
    std::vector<int> index_of_sorted_; // sorted index -> owner_id（如需可用）
    std::vector<int> rank_of_owner_;   // owner_id -> sorted index
    std::vector<Node> tree_;

    void activate_internal(int node, int nl, int nr,
                           int idx, int owner_id) {
        if (idx < nl || idx >= nr) return;
        Node& nd = tree_[node];
        auto it = nd.pos.find(owner_id);
        if (it == nd.pos.end()) {
            int pos = static_cast<int>(nd.active.size());
            nd.active.push_back(owner_id);
            nd.pos.emplace(owner_id, pos);
        }
        if (nl + 1 == nr) return;
        int mid = (nl + nr) / 2;
        activate_internal(node * 2, nl, mid, idx, owner_id);
        activate_internal(node * 2 + 1, mid, nr, idx, owner_id);
    }

    void deactivate_internal(int node, int nl, int nr,
                             int idx, int owner_id) {
        if (idx < nl || idx >= nr) return;
        Node& nd = tree_[node];
        auto it = nd.pos.find(owner_id);
        if (it != nd.pos.end()) {
            int pos_idx = it->second;
            int last_id = nd.active.back();
            nd.active[pos_idx] = last_id;
            nd.pos[last_id] = pos_idx;
            nd.active.pop_back();
            nd.pos.erase(it);
        }
        if (nl + 1 == nr) return;
        int mid = (nl + nr) / 2;
        deactivate_internal(node * 2, nl, mid, idx, owner_id);
        deactivate_internal(node * 2 + 1, mid, nr, idx, owner_id);
    }

    int count_internal(int node, int nl, int nr,
                       int L, int R) const {
        if (nr <= L || R <= nl) return 0;
        if (L <= nl && nr <= R) {
            return static_cast<int>(tree_[node].active.size());
        }
        int mid = (nl + nr) / 2;
        return count_internal(node * 2, nl, mid, L, R)
             + count_internal(node * 2 + 1, mid, nr, L, R);
    }

    void enumerate_internal(int node, int nl, int nr,
                            int L, int R,
                            std::vector<int>& out_owner_ids) const {
        if (nr <= L || R <= nl) return;
        if (L <= nl && nr <= R) {
            const Node& nd = tree_[node];
            out_owner_ids.insert(out_owner_ids.end(),
                                 nd.active.begin(), nd.active.end());
            return;
        }
        int mid = (nl + nr) / 2;
        enumerate_internal(node * 2, nl, mid, L, R, out_owner_ids);
        enumerate_internal(node * 2 + 1, mid, nr, L, R, out_owner_ids);
    }

    void collect_nodes(int node, int nl, int nr,
                       int L, int R,
                       std::vector<int>& out_nodes) const {
        if (nr <= L || R <= nl) return;
        if (L <= nl && nr <= R) {
            out_nodes.push_back(node);
            return;
        }
        int mid = (nl + nr) / 2;
        collect_nodes(node * 2, nl, mid, L, R, out_nodes);
        collect_nodes(node * 2 + 1, mid, nr, L, R, out_nodes);
    }

    void build_alias_for_nodes(const std::vector<int>& nodes,
                               std::vector<double>& prob,
                               std::vector<int>& alias) const {
        const int m = static_cast<int>(nodes.size());
        std::vector<double> w(m);
        double sum = 0.0;
        for (int i = 0; i < m; ++i) {
            int node_idx = nodes[i];
            double wi = static_cast<double>(tree_[node_idx].active.size());
            w[i] = wi;
            sum += wi;
        }
        if (sum <= 0.0) {
            for (int i = 0; i < m; ++i) {
                prob[i]  = 1.0;
                alias[i] = i;
            }
            return;
        }

        std::vector<double> scaled(m);
        for (int i = 0; i < m; ++i) {
            scaled[i] = (w[i] * m) / sum;
        }

        std::vector<int> small, large;
        small.reserve(m);
        large.reserve(m);
        for (int i = 0; i < m; ++i) {
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
};

} // namespace rect_sampler
