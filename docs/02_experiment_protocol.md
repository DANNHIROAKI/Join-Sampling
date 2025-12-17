# 02_experiment_protocol.md

本文件给出一个“可直接写进 SIGMOD 论文”的实验协议（protocol），目标是：

- **可复现**：固定随机种子、记录配置、输出足够多信息
- **可比较**：不同 baseline 在同一数据/同一 t 下，跑法一致
- **可解释**：包含 phase breakdown、是否走枚举、误差条等

---

## 1) 统一实验目标

### 1.1 任务

给定两集合 \(R, S\) 的空间 join 结果：

\[
J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq\emptyset\},
\]

实验目标是输出 \(t\) 个样本 \(Z_1,\dots,Z_t\in J\)，满足：

- **均匀**：\(\Pr(Z_i = p)=1/|J|\)
- **i.i.d. 有放回**：样本之间相互独立

### 1.2 三种版本（variant）

- `sampling`：不显式枚举所有 join 对；适合大 join
- `enum_sampling`：枚举出 `Pairs` 再数组采样；适合小 join
- `adaptive`：阈值触发，兼顾两者

在所有实验图表中，应始终把 **method 与 variant 分开报告**，例如：

- `ours/sampling`
- `r_tree/adaptive`
- `range_tree/enum_sampling`

---

## 2) 评价指标（推荐）

### 2.1 性能

- **Wall-clock time (ms)**：整次 run 的耗时（含 build / count / sample）
- **Phase breakdown**：
  - `build_ms`：建立索引/预处理
  - `count_ms`：计数阶段
  - `enum_ms`：枚举阶段（若发生）
  - `sample_ms`：采样阶段
  - `total_ms`
- **Peak memory（可选）**：若你有测量手段（例如 Linux `/usr/bin/time -v`），建议加到论文里

### 2.2 统计正确性（小规模必须做，大规模建议抽检）

- **Oracle correctness（小规模）**：
  - `count_exact`：是否能给出精确 \(|J|\)
  - `oracle_count`：由全枚举 oracle 得到的 \(|J|\)
  - 误差：`abs_err`, `rel_err`

- **采样分布一致性（小规模）**：
  - χ² / KS：对 join pair 的出现频率做检验（详见 `include/sjs/sampling/sample_quality.h`）
  - 自相关：检查输出序列是否存在明显相关性（抽样 i.i.d. 的必要条件）

> 建议在论文里展示至少 1~2 组小规模数据的分布检验图/表。

---

## 3) 复现性协议（强制）

### 3.1 种子

- `--gen_seed`：用于合成数据生成器（决定数据）
- `--seed`：用于 baseline 的采样/随机过程（决定抽样）

必须把这两个 seed 都记录进 `raw.csv`。

### 3.2 repeats / seeds

两种重复策略（二选一）：

- `--repeats=k`：在固定 dataset 下，跑 k 次（每次用 `seed + rep` 或内部策略派生）
- 或在 sweep JSON 中提供 `seed: [s1, s2, ...]`：明确每次 run 的 seed

论文里建议使用 **显式 seed 列表**（更可复现）。

### 3.3 环境记录

建议在论文实验部分记录：

- CPU 型号 / 核数 / 内存
- OS / 编译器版本
- 编译选项（至少 Release / -O3）
- 是否 pin 到单线程（当前默认 `--threads=1`）

---

## 4) 公平比较（baseline 跑法一致）

框架通过 `include/sjs/baselines/runners/*_runner.h` 强制：

- `sampling`：按“纯采样”流程跑（通常是 Phase1+2+3 或 baseline 等价实现）
- `enum_sampling`：统一做 enumerate -> rank/alias -> sample
- `adaptive`：统一阈值逻辑（先 count，再决定是否保留枚举）

因此 **不要在某个 baseline 里偷换流程**（比如把一些工作挪到 build 阶段不计时），否则需要在论文中明确说明并补齐 phase 计时。

---

## 5) 推荐实验矩阵（SIGMOD 风格）

你可以用 `sjs_sweep` 通过 `config/sweeps/*.json` 跑下面这些实验：

### 5.1 规模扩展（Scale）

- 固定数据分布与密度 \(\alpha\)
- 扫 `n_r = n_s = {1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6}`
- 固定 `t`（例如 1e5 或 5e5）

输出：时间 vs n（log-log / semi-log）。

### 5.2 稠密度扩展（Alpha / density）

- 固定 `n_r, n_s`
- 扫 `alpha = {1e-8, 3e-8, 1e-7, ..., 1e-3}`（按对数间隔）
- 观察各方法在小 join 与大 join 下的行为（尤其 adaptive 的分支比例）

### 5.3 抽样数扩展（t）

- 固定 `n, alpha`
- 扫 `t = {1e3, 1e4, 1e5, 1e6}`
- 观察 `O(n log n + t)` 还是 `O(n log n + |J| + t)` 的差异

### 5.4 分布鲁棒性（Distribution）

用不同合成生成器对比：

- stripe（可控 α）
- uniform
- clustered
- hetero_sizes（active-set 压测）

---

## 6) 输出文件与列字段

`sjs_run` 与 `sjs_sweep` 都会写 CSV（TSV）：

- raw：每次 run/repeat 一行
- summary：按 key 聚合（均值、标准差、分位数、误差条）

常见列（示意）：

- dataset / method / variant / seed / rep
- n_r / n_s / t
- wall_ms
- count_value / count_exact / count_stderr / ci_low / ci_high
- used_enumeration / enum_truncated / enum_cap
- adaptive_branch（enumerate 或 sampling）
- phases_json（JSON-lite，可被后处理脚本解析）

---

## 7) 论文写法建议（可直接复用）

- 正确性：用 oracle + χ²/KS 表/图证明 i.i.d. 均匀
- 性能：用 scale/alpha/t 三个维度展示趋势
- 可解释性：给 1 张 phase breakdown 的堆叠图（build/count/sample）
- 工程细节：说明 half-open 语义与事件 tie-break，避免被审稿人抓 “边界相交” 的漏洞
