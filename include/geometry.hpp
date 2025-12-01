#pragma once

#include <cmath>

namespace rect_sampler {

// 简单二维点结构（目前主要是占位，后面如果需要可扩展）
struct Point {
    double x{0.0};
    double y{0.0};
};

// 半开矩形 [x_min, x_max) × [y_min, y_max)
struct Rect {
    int    id{ -1 };      // 来自 CSV 的 id（或你自己赋的索引）
    double x_min{0.0};
    double y_min{0.0};
    double x_max{0.0};
    double y_max{0.0};

    bool   is_c{false};   // true 表示属于 R_c，false 表示属于 R_{\bar c}
};

// 一些便捷函数
inline double width(const Rect& r)  { return r.x_max - r.x_min; }
inline double height(const Rect& r) { return r.y_max - r.y_min; }
inline double area(const Rect& r)   { return width(r) * height(r); }

// 半开矩形相交测试：
// x 维：max(L_a, L_b) < min(R_a, R_b)
// y 维：max(D_a, D_b) < min(U_a, U_b)
inline bool intersects_half_open(const Rect& a, const Rect& b) {
    if (std::max(a.x_min, b.x_min) >= std::min(a.x_max, b.x_max)) return false;
    if (std::max(a.y_min, b.y_min) >= std::min(a.y_max, b.y_max)) return false;
    return true;
}

} // namespace rect_sampler
