# EXP-2 阅读文档：Runtime vs t

## 1. 实验目的与要回答的问题（RQ2）

**EXP-2** 用来回答 RQ2：

> 在固定数据集（固定 $R,S$）与固定密度参数 $\alpha$ 下，随着采样数量 $t$ 从 1k 增加到 1M，各方法的端到端 wall-clock 时间如何变化？是否呈现“采样型随 $t$ 增长；枚举型对 $t$ 的敏感度取决于其额外步骤（如 rank-sampling 排序）”的典型分化？

EXP-2 的核心价值是：**把“生成更多样本”这一维的代价分离出来**，让读者看到方法随 $t$ 的扩展性，并为后续 RQ3/RQ4/RQ6 的趋势解释提供 baseline。

------

## 2. 任务定义与语义约束（必须写清）

- 输入：两集合 $R,S$，元素为二维轴对齐矩形盒，采用 **half-open** 语义：
  $$
  r=[L_x(r),R_x(r))\times[L_y(r),R_y(r)),\quad
  s=[L_x(s),R_x(s))\times[L_y(s),R_y(s)).
  $$

- Join 输出：
  $$
  J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq\emptyset\}.
  $$

- 目标：从 $J$ 中抽取 **$t$ 个 i.i.d.、均匀、有放回**样本。

> 解释时请强调：**边界贴合不算相交**（half-open），所有方法都必须严格遵守同一判定，否则不能放在一起比较。

------

## 3. EXP-2 的实验设计（“固定什么、扫什么”）

### 3.1 固定项

- **固定数据集**：同一组 $R,S$ 在 sweep 全程不变。
- **固定密度参数 $\alpha$**：EXP-2 不扫 $\alpha$。
- **固定线程数**：建议单线程（公平、可解释）。
- **固定是否写样本**：建议关闭 `write_samples`（避免 I/O 噪声）。

### 3.2 扫描项

- 只 sweep **$t$**：典型序列为
  $$
  t\in \{10^3,3\cdot10^3,10^4,3\cdot10^4,10^5,3\cdot10^5,10^6\}.
  $$

### 3.3 方法与三种运行模式（论文推荐固定写法）

每个方法都以三种模式评估：

1. **Sampling**：不物化全部 $J$，直接输出 $t$ 个样本。

2. **Enum+Sampling（枚举基线）**：先枚举 $J$，再从枚举流中做均匀有放回抽样。

3. **Adaptive**：用阈值/探测策略在 “枚举” 与 “采样” 之间切换，目标接近
   $$
   \min\{T(\text{enumerate}),T(\text{sampling})\}.
   $$

------

## 4. 计时口径与统计口径（读图前必须明确）

### 4.1 计时覆盖范围

EXP-2 的主指标是 **端到端 wall-clock 时间**，覆盖方法内部典型阶段（随实现不同会合并/拆分）：

- Build：索引/预处理
- Count：计数或权重计算
- Enumerate（如适用）：枚举 join 流
- Sample：生成 $t$ 个样本（含二次扫描/定位）

**不计入**：数据生成/数据加载、写样本到磁盘等非算法主体开销。

### 4.2 统计方式

- 每个配置点重复运行 $k$ 次（常用 3 或 5；如果要用 p95 更推荐 ≥10）。
- 报告：
  - 主线：`wall_median_ms`
  - 误差条：`p95` 或 `stdev`

------

## 5. 配置文件与运行入口（你应该看哪里）

### 5.1 Sweep 配置文件

EXP-2 推荐直接使用仓库提供的 sweep 配置（或复制后修改）：

- `config/sweeps/sweep_t.json`（固定数据与 $\alpha$，只扫 $t$）

### 5.2 一键运行脚本（推荐）

仓库的 EXP-2 runner 脚本会执行：Release 编译 → sweep → 生成图。
 对应文件：`run/.sh` 

常用方式：

```
bash run/.sh
```

或指定输出目录：

```
bash run/.sh --out_dir results/exp2/runtime_vs_t/my_run
```

> 脚本会把环境信息写入 `logs/env.txt`，便于复现（CPU/编译器/commit 等）。

------

## 6. EXP-2 的输出文件结构（你会拿哪些东西写论文）

一次 EXP-2 运行完成后，输出目录下至少包含：

- `sweep_raw.csv`：**每次 repeat 一行**（含阶段耗时、随机种子、是否成功等）
- `sweep_summary.csv`：**聚合后的统计表**（按 method×variant×t 聚合 median/p95/stdev/ok_rate）

如果启用了绘图脚本（默认启用），还会生成图与说明文件（见下一节）。

------

## 7. 绘图与读图：推荐的“三张图”组合

仓库的绘图脚本 `run/.py`  会自动从 raw+summary 生成以下结果（每个 variant 一张）：

### 图 A：`runtime vs t`（主图）

- 文件名：`exp2_runtime_vs_t_<variant>.png/.pdf`
- x 轴：$t$（log）
- y 轴：wall runtime（log）
- 误差条：默认 **upper-only p95-median**（避免 log-y 下界 <0 的问题）

**读法**：这是“端到端用户感知延迟”曲线，用于论文主结论展示。

------

### 图 B：`Δruntime vs t`（强烈建议论文也给）

- 文件名：`exp2_delta_runtime_vs_t_<variant>.png/.pdf`

- 定义：对每条曲线减去基线点（默认 $t_0=1000$）：
  $$
  \Delta T(t)=T(t)-T(t_0)
  $$

- 作用：**把常数开销（Build/Count/Enumerate）抵消掉**，让“随 $t$ 的增长项”更显著。

**读法**：

- 若 Sampling 的 $\Delta T(t)$ 近线性，说明“生成样本的增量成本”接近 $O(t)$。
- 若 Enum+Sampling 的 $\Delta T(t)$ 明显增长，通常来自 rank-sampling 的排序/定位（常见 $O(t\log t)$）。

------

### 图 C：`sample-phase vs t`（解释型 profiling 图）

- 文件名：`exp2_sample_phase_vs_t_<variant>.png/.pdf`

- 从 `sweep_raw.csv` 的 `phases_json` 中提取 `run_sample_ms`（或兼容键名），对每点取 median/p95。

- 额外输出：`exp2_ns_per_sample.csv`
   用端点斜率估计每样本纳秒级成本：
  $$
  \text{ns/sample}\approx \frac{\Delta (\text{sample\_ms})}{\Delta t}\times 10^6
  $$

**读法**：这是“只看 Sample 阶段”随 $t$ 的 scaling，最适合解释“为什么 wall 主图看起来不怎么变/为什么某条曲线斜率更小”。

------

## 8. 理论预期与“看到什么算正常”

为了写论文时“结论-解释闭环”，EXP-2 推荐用下面这套固定叙事：

### 8.1 Sampling 模式：通常 ~线性增长

Sampling 需要真实产出 $t$ 个样本，因此**随 $t$ 增长是必然的**：

- 预期：`Δruntime`、`sample-phase` 近似线性；
- 解释：Phase-3（或 Sample 阶段）每个样本通常 $O(1)$ 摊还，整体 $O(t)$。

### 8.2 Enum+Sampling 模式：可能出现 $t\log t$ 项

在你们实现框架里，Enum+Sampling 通常包含：

- 生成 $t$ 个 ranks 并排序；
- 再扫描 join 流定位到这些 ranks；

因此：

- 当 $t$ 很大时，可能明显看到 **$t\log t$** 的增长；
- 当 $|J|$ 很大时，枚举扫描成本会主导，曲线可能对 $t$ “不那么敏感”。

> 论文写法建议：不要笼统写“枚举型对 t 不敏感”，而写“枚举型总成本由 $|J|$ 枚举 + $t\log t$ rank sampling 共同决定，在不同区间主导项不同”。

### 8.3 Adaptive：是否“切换”取决于 $J_\star$ 与 $|J|$

- 如果 $|J| < J_\star$：Adaptive 可能全程走枚举分支 → 曲线看起来几乎不随 $t$ 变（或仅很弱）。
- 如果 $|J| > J_\star$：Adaptive 会切到 sampling → 曲线更像 Sampling。

**读图时应同时报告** raw 中的分支统计（例如 `used_enumeration` 或 pilot 规模），否则读者会问“为什么 adaptive 不切换/为什么和某个模式一模一样”。

------

## 9. 公平性与失败处理（读 summary 前先看这些列）

### 9.1 过滤规则（绘图默认做）

- 过滤 `ok_rate < 1.0` 的点（避免失败/跳过点污染）
- 但论文里应在 caption/正文说明：哪些方法/组合被跳过或失败（例如某些 rejection 组合可能被禁用）

### 9.2 `enum_cap` 与“截断点”

如果启用了枚举上限（防止爆内存），那么：

- 图上应把截断点显式标注为失败或缺失；
- 不能静默过滤而不说明（否则不公平）。

------

## 10. 推荐写入论文的 EXP-2 小节结构（模板）

你可以直接按以下结构写实验小节（每段对应 1–2 句话即可）：

1. **Setup**：固定数据集（给出 $n_r,n_s,\alpha$）、单线程、重复 $k$ 次、报告 median+p95（或 stdev）。
2. **Result (Sampling)**：runtime 随 $t$ 增长，Δruntime/sample-phase 近线性；给出 ns/sample 量级。
3. **Result (Enum+Sampling)**：在大 t 区间出现明显增长（解释为 rank sort 的 $t\log t$）；或在高 $|J|$ 区间枚举主导。
4. **Result (Adaptive)**：根据 $|J|$ 与 $J_\star$ 的关系展示分支选择；讨论其是否接近 $\min\{\cdot\}$。
5. **Takeaway**：一句话总结“谁的斜率更小/谁的固定开销更大/哪类方法对 t 更稳”。

------

## 11. 常见坑与排障清单（跑完 EXP-2 必查）

1. `count_mean` 在所有 t 下是否一致？（固定数据集应一致）
2. `ok_rate` 是否为 1？若不是，失败原因是否在 raw 里可追踪？
3. 主图看起来不随 t 变：是否因为常数开销过大？请看 Δruntime 与 sample-phase 图。
4. 使用 p95 误差条但 repeats 太少：建议改用 stdev 或增加 repeats。
5. log-y 下误差条出现负下界：必须用 upper-only（绘图脚本已处理）。

------

## 12. 本实验与整体研究问题的关系（你在论文里怎么“接上”）

- EXP-2 直接回答 **RQ2（Runtime vs t）**；
- EXP-2 的样本阶段斜率（ns/sample）和常数项分解，将在：
  - EXP-3（扫 $\alpha$）解释“稠密导致 Build/Count/Enumerate 退化”
  - EXP-6（adaptive heatmap）解释“为什么 adaptive 在某些区间更接近 min”
  - EXP-7（phase profiling）成为“解释性证据”

------

### 附：与你们脚本输出的文件名对照表（便于读者定位）

- 一键跑：`run/.sh` 
- 绘图与解读：`run/.py` 
- 原始记录：`sweep_raw.csv`
- 聚合统计：`sweep_summary.csv`
- 主图：`exp2_runtime_vs_t_<variant>.(png|pdf)`
- 增量图：`exp2_delta_runtime_vs_t_<variant>.(png|pdf)`
- sample-phase 图：`exp2_sample_phase_vs_t_<variant>.(png|pdf)`
- 每样本成本：`exp2_ns_per_sample.csv`