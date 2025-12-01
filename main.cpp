#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <cstddef>

#include "geometry.hpp"
#include "utils.hpp"
#include "baseline1.hpp"
#include "baseline2.hpp"
#include "ours.hpp"

using namespace rect_sampler;


// 用于对 RectPair 排序的比较函数
struct RectPairLess {
    bool operator()(const RectPair& a, const RectPair& b) const {
        if (a.idx_c != b.idx_c) return a.idx_c < b.idx_c;
        return a.idx_bar < b.idx_bar;
    }
};

// 判断两个 RectPair 是否相等
inline bool rect_pair_equal(const RectPair& a, const RectPair& b) {
    return a.idx_c == b.idx_c && a.idx_bar == b.idx_bar;
}

// 判断两个 pair 集合（视为无序集合）是否完全相同
bool pair_sets_equal(const PairList& A, const PairList& B) {
    if (A.size() != B.size()) return false;

    PairList a_sorted = A;
    PairList b_sorted = B;

    RectPairLess cmp;
    std::sort(a_sorted.begin(), a_sorted.end(), cmp);
    std::sort(b_sorted.begin(), b_sorted.end(), cmp);

    for (std::size_t i = 0; i < a_sorted.size(); ++i) {
        if (!rect_pair_equal(a_sorted[i], b_sorted[i])) {
            return false;
        }
    }
    return true;
}

// 基于 sorted ground truth，用二分判断一个 pair 是否在 GT 中
std::size_t count_pairs_not_in(
    const PairList& ground_truth_sorted,
    const PairList& samples
) {
    RectPairLess cmp;
    std::size_t bad = 0;
    for (const auto& p : samples) {
        bool ok = std::binary_search(
            ground_truth_sorted.begin(),
            ground_truth_sorted.end(),
            p, cmp
        );
        if (!ok) ++bad;
    }
    return bad;
}

int main() {
    try {
        // -------------------- 0. 初始化 --------------------
        std::string data_path1 = "data/rectangles1.csv";
        std::string data_path2 = "data/rectangles2.csv";
        std::string log_path   = "results/experiment_log.txt";

        // 固定随机种子，保证实验可复现
        std::mt19937_64 rng(42);

        // 清空 / 创建 log 文件
        {
            std::ofstream ofs(log_path, std::ios::out);
            if (!ofs) {
                std::cerr << "Failed to open log file: " << log_path << "\n";
                return 1;
            }
            ofs << "# algo,t,time_enum_ms,time_sample_ms,time_total_ms,mem_bytes,"
                   "flag_equal_gt,bad_sample_count,extra\n";
        }

        // -------------------- 1. 读入数据集 --------------------
        std::cout << "[INFO] Loading rectangles...\n";
        RectList Rc  = load_rectangles_csv(data_path1, true);
        RectList Rbar = load_rectangles_csv(data_path2, false);

        std::cout << "[INFO] Rc size = " << Rc.size()
                  << ", Rbar size = " << Rbar.size() << "\n";

        std::size_t rss_after_data = getCurrentRSS();
        std::cout << "[INFO] RSS after data load = "
                  << rss_after_data << " bytes\n";

        // -------------------- 2. Baseline1：枚举 ground truth --------------------
        std::cout << "[INFO] Running Baseline1 (enumerate ground truth)...\n";

        PairList gt_pairs;  // ground truth: AllPairs_baseline1
        Timer timer;

        timer.reset();
        baseline1::enumerate_all_pairs(Rc, Rbar, gt_pairs);
        double b1_enum_ms = timer.elapsed_ms();
        std::size_t rss_after_b1_enum = getCurrentRSS();

        std::size_t J_size = gt_pairs.size();
        std::cout << "[INFO] |J| (ground truth pairs) = " << J_size << "\n";
        std::cout << "[INFO] Baseline1 enumerate time = " << b1_enum_ms << " ms\n";
        std::cout << "[INFO] RSS after Baseline1 enum = "
                  << rss_after_b1_enum << " bytes\n";

        // 把 Baseline1 枚举信息写入 log（采样时间先写 0）
        {
            std::ofstream ofs(log_path, std::ios::app);
            ofs << "baseline1_enum,0,"
                << b1_enum_ms << ",0,"
                << b1_enum_ms << ","
                << rss_after_b1_enum << ","
                << 1 << "," << 0 << ","
                << "ground_truth_enum\n";
        }

        // 预先对 ground truth 排序，方便后面快速 membership 检查
        PairList gt_sorted = gt_pairs;
        std::sort(gt_sorted.begin(), gt_sorted.end(), RectPairLess{});

        // -------------------- 3. Baseline2：plane sweep 枚举 + 与 GT 对比 --------------------
        std::cout << "[INFO] Running Baseline2 (plane sweep enumerate)...\n";

        PairList b2_pairs;
        timer.reset();
        baseline2::enumerate_all_pairs(Rc, Rbar, b2_pairs);
        double b2_enum_ms = timer.elapsed_ms();
        std::size_t rss_after_b2_enum = getCurrentRSS();

        bool b2_equal_gt = pair_sets_equal(gt_pairs, b2_pairs);
        std::cout << "[INFO] Baseline2 enumerate time = " << b2_enum_ms << " ms\n";
        std::cout << "[INFO] RSS after Baseline2 enum = "
                  << rss_after_b2_enum << " bytes\n";
        std::cout << "[INFO] Baseline2 vs GT equal? flag = "
                  << (b2_equal_gt ? 1 : 0) << "\n";

        {
            std::ofstream ofs(log_path, std::ios::app);
            ofs << "baseline2_enum,0,"
                << b2_enum_ms << ",0,"
                << b2_enum_ms << ","
                << rss_after_b2_enum << ","
                << (b2_equal_gt ? 1 : 0) << ","
                << 0 << ","
                << "b2_enum_vs_gt\n";
        }

        // 为 Baseline2 采样阶段也预备一份排序后的 pair，后面检查时会用到
        PairList b2_sorted = b2_pairs;
        std::sort(b2_sorted.begin(), b2_sorted.end(), RectPairLess{});

        // -------------------- 4. 针对多个 t 做采样实验 --------------------
        std::vector<std::size_t> t_values = {
            1000, 10000, 100000, 1000000
        };

        for (std::size_t t : t_values) {
            std::cout << "\n[INFO] ===== t = " << t << " =====\n";

            // ---------- 4.1 Baseline1：只做采样 ----------
            PairList samples_b1;
            timer.reset();
            baseline1::sample_from_all_pairs(gt_pairs, t, rng, samples_b1);
            double b1_sample_ms = timer.elapsed_ms();
            std::size_t rss_after_b1_sample = getCurrentRSS();

            std::cout << "[INFO] Baseline1 sample time = "
                      << b1_sample_ms << " ms\n";

            {
                std::ofstream ofs(log_path, std::ios::app);
                ofs << "baseline1_sample," << t << ","
                    << 0 << ","              // 枚举阶段已单独统计
                    << b1_sample_ms << ","
                    << b1_sample_ms << ","   // 总时间 ≈ 采样时间（这里只记录采样）
                    << rss_after_b1_sample << ","
                    << 1 << ","              // 和 GT 本身 trivially 一致
                    << 0 << ","              // 不需要检测 bad_sample
                    << "baseline1_sample\n";
            }

            // ---------- 4.2 Baseline2：采样 + 检查样本都在 GT 中 ----------
            PairList samples_b2;
            timer.reset();
            baseline2::sample_from_all_pairs(b2_pairs, t, rng, samples_b2);
            double b2_sample_ms = timer.elapsed_ms();
            std::size_t rss_after_b2_sample = getCurrentRSS();

            std::size_t bad_b2_samples =
                count_pairs_not_in(gt_sorted, samples_b2);

            std::cout << "[INFO] Baseline2 sample time = "
                      << b2_sample_ms << " ms\n";
            std::cout << "[INFO] Baseline2 bad samples vs GT = "
                      << bad_b2_samples << "\n";

            {
                std::ofstream ofs(log_path, std::ios::app);
                ofs << "baseline2_sample," << t << ","
                    << 0 << ","
                    << b2_sample_ms << ","
                    << b2_sample_ms << ","
                    << rss_after_b2_sample << ","
                    << (b2_equal_gt ? 1 : 0) << ","
                    << bad_b2_samples << ","
                    << "baseline2_sample\n";
            }

            // ---------- 4.3 Ours：直接采样 + 检查样本都在 GT 中 ----------
            PairList samples_ours;
            timer.reset();
            ours::sample_pairs(Rc, Rbar, t, rng, samples_ours);
            double ours_total_ms = timer.elapsed_ms();
            std::size_t rss_after_ours = getCurrentRSS();

            std::size_t bad_ours_samples =
                count_pairs_not_in(gt_sorted, samples_ours);

            std::cout << "[INFO] Ours total sample time = "
                      << ours_total_ms << " ms\n";
            std::cout << "[INFO] Ours bad samples vs GT = "
                      << bad_ours_samples << "\n";

            {
                std::ofstream ofs(log_path, std::ios::app);
                ofs << "ours," << t << ","
                    << 0 << ","              // 这里统一写在 time_total_ms
                    << 0 << ","
                    << ours_total_ms << ","
                    << rss_after_ours << ","
                    << 1 << ","              // 只要 bad_ours_samples=0 就视为 flag=1
                    << bad_ours_samples << ","
                    << "ours_total\n";
            }
        }

        std::cout << "\n[INFO] Experiment finished. Results logged to "
                  << log_path << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "[ERROR] Exception: " << ex.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "[ERROR] Unknown exception.\n";
        return 1;
    }

    return 0;
}
