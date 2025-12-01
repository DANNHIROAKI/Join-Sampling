#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>

#include "geometry.hpp"

namespace rect_sampler {

// ======================== CSV 读取 ========================

// 从 CSV 读取矩形数据：
// 格式：id,x_min,y_min,x_max,y_max,width,height,area
// 只使用前 5 列，其余列直接忽略
inline std::vector<Rect> load_rectangles_csv(const std::string& path,
                                             bool is_c_flag) {
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("Failed to open CSV file: " + path);
    }

    std::vector<Rect> rects;
    std::string line;

    // 跳过表头
    if (!std::getline(ifs, line)) {
        return rects;
    }

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string field;

        Rect r;

        // id
        if (!std::getline(ss, field, ',')) continue;
        r.id = std::stoi(field);

        // x_min
        if (!std::getline(ss, field, ',')) continue;
        r.x_min = std::stod(field);

        // y_min
        if (!std::getline(ss, field, ',')) continue;
        r.y_min = std::stod(field);

        // x_max
        if (!std::getline(ss, field, ',')) continue;
        r.x_max = std::stod(field);

        // y_max
        if (!std::getline(ss, field, ',')) continue;
        r.y_max = std::stod(field);

        r.is_c = is_c_flag;

        rects.push_back(r);
    }

    return rects;
}

// ======================== 简单计时器 ========================

class Timer {
public:
    using clock = std::chrono::high_resolution_clock;

    Timer() : start_(clock::now()) {}

    void reset() {
        start_ = clock::now();
    }

    // 返回毫秒
    double elapsed_ms() const {
        auto now = clock::now();
        std::chrono::duration<double, std::milli> diff = now - start_;
        return diff.count();
    }

    // 返回秒
    double elapsed_sec() const {
        auto now = clock::now();
        std::chrono::duration<double> diff = now - start_;
        return diff.count();
    }

private:
    clock::time_point start_;
};

// RAII 计时 + 自动打印/记录（可选用）
class ScopedTimer {
public:
    ScopedTimer(const std::string& name, bool auto_print = true)
        : name_(name), auto_print_(auto_print) {}

    ~ScopedTimer() {
        if (auto_print_) {
            std::cerr << "[TIMER] " << name_ << " took "
                      << timer_.elapsed_ms() << " ms\n";
        }
    }

    double elapsed_ms() const { return timer_.elapsed_ms(); }

private:
    std::string name_;
    bool auto_print_;
    Timer timer_;
};

// ======================== 内存监控（当前 RSS） ========================

inline std::size_t getCurrentRSS() {
#if defined(__linux__)
    // Linux: 读取 /proc/self/statm
    long rss = 0L;
    FILE* fp = std::fopen("/proc/self/statm", "r");
    if (!fp) {
        return 0;
    }
    long dummy = 0L;
    if (std::fscanf(fp, "%ld %ld", &dummy, &rss) != 2) {
        std::fclose(fp);
        return 0;
    }
    std::fclose(fp);
    long page_size_kb = sysconf(_SC_PAGESIZE); // in bytes
    return static_cast<std::size_t>(rss) * static_cast<std::size_t>(page_size_kb);
#else
    // 其他平台这里留空实现（返回 0），需要的话可以自行扩展
    return 0;
#endif
}

// 在 main 中可以自己维护一个全局 max_rss = max(max_rss, getCurrentRSS());

// ======================== 简单日志工具 ========================

inline void append_line_to_file(const std::string& path,
                                const std::string& line) {
    std::ofstream ofs(path, std::ios::app);
    if (!ofs) {
        throw std::runtime_error("Failed to open log file: " + path);
    }
    ofs << line << '\n';
}

} // namespace rect_sampler
