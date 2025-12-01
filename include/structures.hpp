#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <cstddef>

#include "geometry.hpp"

namespace rect_sampler {

// ======================== 动态 stabbing 结构 ========================

class DynamicStabbing1D {
public:
    DynamicStabbing1D() = default;

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

    bool ready() const { return seg_count_ > 0; }

    void activate(int interval_id, double y_low, double y_high) {
        if (!ready()) return;
        int L = coord_to_seg_index(y_low);
        int R = coord_to_seg_index(y_high);
        if (L < 0 || R < 0 || L >= R) return;
        activate_internal(1, 0, seg_count_, L, R, interval_id);
    }

    void deactivate(int interval_id, double y_low, double y_high) {
        if (!ready()) return;
        int L = coord_to_seg_index(y_low);
        int R = coord_to_seg_index(y_high);
        if (L < 0 || R < 0 || L >= R) return;
        deactivate_internal(1, 0, seg_count_, L, R, interval_id);
    }

    std::size_t count(double q) const {
        if (!ready()) return 0;
        if (q < coords_.front() || q >= coords_.back()) return 0;
        int s = seg_index_for_point(q);
        if (s < 0) return 0;

        std::vector<int> path;
        collect_path(1, 0, seg_count_, s, path);

        std::size_t total = 0;
        for (int node : path) {
            total += tree_[node].active.size();
        }
        return total;
    }

    void enumerate(double q, std::vector<int>& out_ids) const {
        out_ids.clear();
        if (!ready()) return;
        if (q < coords_.front() || q >= coords_.back()) return;
        int s = seg_index_for_point(q);
        if (s < 0) return;

        std::vector<int> path;
        collect_path(1, 0, seg_count_, s, path);

        for (int node : path) {
            const Node& nd = tree_[node];
            out_ids.insert(out_ids.end(), nd.active.begin(), nd.active.end());
        }
    }

    template <class URNG>
    void sample(double q, std::size_t sample_count, URNG& rng,
                std::vector<int>& out_ids) const {
        out_ids.clear();
        if (!ready() || sample_count == 0) return;
        if (q < coords_.front() || q >= coords_.back()) return;
        int s = seg_index_for_point(q);
        if (s < 0) return;

        std::vector<int> path;
        collect_path(1, 0, seg_count_, s, path);

        std::vector<int> nodes;
        nodes.reserve(path.size());
        for (int node : path) {
            if (!tree_[node].active.empty()) nodes.push_back(node);
        }
        if (nodes.empty()) return;

        int m = static_cast<int>(nodes.size());
        std::vector<double> prob(m);
        std::vector<int> alias(m);
        build_alias_for_nodes(nodes, prob, alias);

        std::uniform_int_distribution<int> col_dist(0, m - 1);
        std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

        out_ids.reserve(sample_count);
        for (std::size_t k = 0; k < sample_count; ++k) {
            int col = col_dist(rng);
            double p = prob_dist(rng);
            int bucket = (p < prob[col]) ? col : alias[col];
            int node_idx = nodes[bucket];
            const Node& nd = tree_[node_idx];
            std::uniform_int_distribution<int> item_dist(
                0, static_cast<int>(nd.active.size()) - 1);
            int pos = item_dist(rng);
            out_ids.push_back(nd.active[pos]);
        }
    }

private:
    struct Node {
        std::vector<int> active;
        std::unordered_map<int,int> pos;
    };

    std::vector<double> coords_;
    int seg_count_{0};
    std::vector<Node> tree_;

    int coord_to_seg_index(double y) const {
        auto it = std::lower_bound(coords_.begin(), coords_.end(), y);
        if (it == coords_.end()) return -1;
        int idx = static_cast<int>(it - coords_.begin());
        if (idx < 0 || idx >= static_cast<int>(coords_.size())) return -1;
        return idx;
    }

    int seg_index_for_point(double q) const {
        auto it = std::upper_bound(coords_.begin(), coords_.end(), q);
        int idx = static_cast<int>(it - coords_.begin()) - 1;
        if (idx < 0 || idx >= seg_count_) return -1;
        return idx;
    }

    void activate_internal(int node, int nl, int nr,
                           int L, int R, int id) {
        if (nr <= L || R <= nl) return;
        if (L <= nl && nr <= R) {
            Node& nd = tree_[node];
            if (nd.pos.find(id) == nd.pos.end()) {
                int pos = static_cast<int>(nd.active.size());
                nd.active.push_back(id);
                nd.pos.emplace(id, pos);
            }
            return;
        }
        int mid = (nl + nr) / 2;
        activate_internal(node*2, nl, mid, L, R, id);
        activate_internal(node*2+1, mid, nr, L, R, id);
    }

    void deactivate_internal(int node, int nl, int nr,
                             int L, int R, int id) {
        if (nr <= L || R <= nl) return;
        if (L <= nl && nr <= R) {
            Node& nd = tree_[node];
            auto it = nd.pos.find(id);
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
        deactivate_internal(node*2, nl, mid, L, R, id);
        deactivate_internal(node*2+1, mid, nr, L, R, id);
    }

    void collect_path(int node, int nl, int nr,
                      int leaf, std::vector<int>& path) const {
        if (leaf < nl || leaf >= nr) return;
        path.push_back(node);
        if (nl + 1 == nr) return;
        int mid = (nl + nr) / 2;
        if (leaf < mid) collect_path(node*2, nl, mid, leaf, path);
        else            collect_path(node*2+1, mid, nr, leaf, path);
    }

    void build_alias_for_nodes(const std::vector<int>& nodes,
                               std::vector<double>& prob,
                               std::vector<int>& alias) const {
        int m = static_cast<int>(nodes.size());
        prob.assign(m, 0.0);
        alias.assign(m, 0);
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
                prob[i] = 1.0;
                alias[i] = i;
            }
            return;
        }
        std::vector<double> scaled(m);
        for (int i = 0; i < m; ++i) {
            scaled[i] = w[i] * m / sum;
        }
        std::vector<int> small, large;
        small.reserve(m); large.reserve(m);
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

// ======================== 动态 range 结构 ========================

class RangeTree1D {
public:
    RangeTree1D() = default;

    explicit RangeTree1D(const std::vector<double>& values) {
        init(values);
    }

    void init(const std::vector<double>& values) {
        N_ = static_cast<int>(values.size());
        if (N_ <= 0) {
            xs_.clear();
            rank_of_owner_.clear();
            tree_.clear();
            return;
        }

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
        rank_of_owner_.assign(N_, -1);
        for (int k = 0; k < N_; ++k) {
            xs_[k] = tmp[k].first;
            int owner = tmp[k].second;
            rank_of_owner_[owner] = k;
        }

        tree_.assign(4 * N_ + 4, Node{});
    }

    bool ready() const { return N_ > 0; }

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

    std::size_t count(double ell, double r) const {
        if (!ready()) return 0;
        if (ell >= r) return 0;
        int L = static_cast<int>(std::upper_bound(xs_.begin(), xs_.end(), ell) - xs_.begin());
        int R = static_cast<int>(std::lower_bound(xs_.begin(), xs_.end(), r) - xs_.begin());
        if (L >= R) return 0;
        return count_internal(1, 0, N_, L, R);
    }

    void enumerate(double ell, double r, std::vector<int>& out_owner_ids) const {
        out_owner_ids.clear();
        if (!ready()) return;
        if (ell >= r) return;
        int L = static_cast<int>(std::upper_bound(xs_.begin(), xs_.end(), ell) - xs_.begin());
        int R = static_cast<int>(std::lower_bound(xs_.begin(), xs_.end(), r) - xs_.begin());
        if (L >= R) return;
        enumerate_internal(1, 0, N_, L, R, out_owner_ids);
    }

    template <class URNG>
    void sample(double ell, double r, std::size_t sample_count,
                URNG& rng, std::vector<int>& out_owner_ids) const {
        out_owner_ids.clear();
        if (!ready() || sample_count == 0) return;
        if (ell >= r) return;
        int L = static_cast<int>(std::upper_bound(xs_.begin(), xs_.end(), ell) - xs_.begin());
        int R = static_cast<int>(std::lower_bound(xs_.begin(), xs_.end(), r) - xs_.begin());
        if (L >= R) return;

        std::vector<int> nodes;
        collect_nodes(1, 0, N_, L, R, nodes);

        std::vector<int> non_empty;
        non_empty.reserve(nodes.size());
        for (int node : nodes) {
            if (!tree_[node].active.empty()) non_empty.push_back(node);
        }
        if (non_empty.empty()) return;

        int m = static_cast<int>(non_empty.size());
        std::vector<double> prob(m);
        std::vector<int> alias(m);
        build_alias_for_nodes(non_empty, prob, alias);

        std::uniform_int_distribution<int> col_dist(0, m - 1);
        std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

        out_owner_ids.reserve(sample_count);
        for (std::size_t k = 0; k < sample_count; ++k) {
            int col = col_dist(rng);
            double p = prob_dist(rng);
            int bucket = (p < prob[col]) ? col : alias[col];
            int node_idx = non_empty[bucket];
            const Node& nd = tree_[node_idx];
            std::uniform_int_distribution<int> item_dist(
                0, static_cast<int>(nd.active.size()) - 1);
            int pos = item_dist(rng);
            out_owner_ids.push_back(nd.active[pos]);
        }
    }

private:
    struct Node {
        std::vector<int> active;
        std::unordered_map<int,int> pos;
    };

    int N_{0};
    std::vector<double> xs_;
    std::vector<int> rank_of_owner_;
    std::vector<Node> tree_;

    void activate_internal(int node, int nl, int nr,
                           int idx, int owner_id) {
        if (idx < nl || idx >= nr) return;
        Node& nd = tree_[node];
        if (nd.pos.find(owner_id) == nd.pos.end()) {
            int pos = static_cast<int>(nd.active.size());
            nd.active.push_back(owner_id);
            nd.pos.emplace(owner_id, pos);
        }
        if (nl + 1 == nr) return;
        int mid = (nl + nr) / 2;
        activate_internal(node*2, nl, mid, idx, owner_id);
        activate_internal(node*2+1, mid, nr, idx, owner_id);
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
        deactivate_internal(node*2, nl, mid, idx, owner_id);
        deactivate_internal(node*2+1, mid, nr, idx, owner_id);
    }

    std::size_t count_internal(int node, int nl, int nr,
                               int L, int R) const {
        if (nr <= L || R <= nl) return 0;
        if (L <= nl && nr <= R) {
            return tree_[node].active.size();
        }
        int mid = (nl + nr) / 2;
        return count_internal(node*2, nl, mid, L, R)
             + count_internal(node*2+1, mid, nr, L, R);
    }

    void enumerate_internal(int node, int nl, int nr,
                            int L, int R,
                            std::vector<int>& out) const {
        if (nr <= L || R <= nl) return;
        if (L <= nl && nr <= R) {
            const Node& nd = tree_[node];
            out.insert(out.end(), nd.active.begin(), nd.active.end());
            return;
        }
        int mid = (nl + nr) / 2;
        enumerate_internal(node*2, nl, mid, L, R, out);
        enumerate_internal(node*2+1, mid, nr, L, R, out);
    }

    void collect_nodes(int node, int nl, int nr,
                       int L, int R,
                       std::vector<int>& out) const {
        if (nr <= L || R <= nl) return;
        if (L <= nl && nr <= R) {
            out.push_back(node);
            return;
        }
        int mid = (nl + nr) / 2;
        collect_nodes(node*2, nl, mid, L, R, out);
        collect_nodes(node*2+1, mid, nr, L, R, out);
    }

    void build_alias_for_nodes(const std::vector<int>& nodes,
                               std::vector<double>& prob,
                               std::vector<int>& alias) const {
        int m = static_cast<int>(nodes.size());
        prob.assign(m, 0.0);
        alias.assign(m, 0);
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
                prob[i] = 1.0;
                alias[i] = i;
            }
            return;
        }
        std::vector<double> scaled(m);
        for (int i = 0; i < m; ++i) {
            scaled[i] = w[i] * m / sum;
        }
        std::vector<int> small, large;
        small.reserve(m); large.reserve(m);
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
