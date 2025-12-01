#pragma once

#include <vector>
#include <algorithm>

namespace rect_sampler {

struct Rect {
    int    id{ -1 };
    double x_min{0.0};
    double y_min{0.0};
    double x_max{0.0};
    double y_max{0.0};
    bool   is_c{false};   // true => R_c, false => R_{\bar c}
};

// 一个跨集合相交 pair，用索引表示（分别是 Rc、Rbar 中的下标）
struct RectPair {
    int idx_c{ -1 };
    int idx_bar{ -1 };
};

using RectList = std::vector<Rect>;
using PairList = std::vector<RectPair>;

// 半开矩形相交判定：[x_min,x_max) × [y_min,y_max)
inline bool rectangles_intersect(const Rect& a, const Rect& b) {
    if (std::max(a.x_min, b.x_min) >= std::min(a.x_max, b.x_max)) return false;
    if (std::max(a.y_min, b.y_min) >= std::min(a.y_max, b.y_max)) return false;
    return true;
}

} // namespace rect_sampler
