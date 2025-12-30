# EXP-7 阅读文档：阶段分解与解释性 Profiling（Phase Breakdown）

本阅读文档用于帮助读者**理解、复现并正确解读** EXP-7 的输出（表格 + 图），并将其作为论文实验章中的“解释性证据”来支撑主结论。文档遵循你给出的整份实验总纲的叙事结构（任务定义 → RQ → 计时口径 → 输出解读 → 失败处理 → 论文呈现）。

------

## 1. 实验定位与目的

### 1.1 EXP-7 在整套实验中的位置

EXP-7 是一个**解释性/剖析型实验（explanatory profiling）**，用于增强可信度与可解释性，主要服务于：

- **RQ2（Runtime vs t）**：解释“为什么某方法随 t 的增长更线性/更平坦”
- **RQ3（Sensitivity to α）**：解释“为什么密度升高后某方法退化/拐点出现”
- **RQ6（Adaptivity effectiveness）**：解释“Adaptive 为什么在某些点选择枚举、在另一些点选择采样，以及代价结构如何变化”

> 重要：EXP-7 不替代主线 sweep（EXP-2/3/4/6），它回答的是“**时间花在哪**”与“**谁在什么阶段赢**”，而不是只报总耗时。

### 1.2 EXP-7 的核心输出是什么

在 2–3 个代表性密度点（α）上，对每个 **method × variant** 的端到端运行时间做顶层阶段拆分：

- **Build**：索引/预处理构建（一次性成本）
- **Count**：计数或权重计算（决定采样分布的前置）
- **Enumerate**：枚举连接流（只在 enum/adaptive pilot 等情形出现）
- **Sample**：产生 t 个 i.i.d. 均匀（有放回）样本（含必要二次扫描/定位）

并用堆叠柱状图/分解表回答：“这一点位上瓶颈到底在哪里？”

------

## 2. 输入工作负载与实验设置（与总纲对齐）

### 2.1 任务定义（与整篇一致）

- 输入：二维 half-open 轴对齐矩形集合 $R,S$
- 连接结果：$J=\{(r,s)\mid r\cap s\neq\emptyset\}$，严格 half-open（边界贴合不算相交）
- 目标：从 $J$ 中抽取 **t 个 i.i.d. 均匀、有放回样本**

EXP-7 不改变任务语义，仅改变**测量方式**（增加阶段计时、解释性输出）。

### 2.2 数据分布：可控密度条带（SYN-CTRL）

EXP-7 默认使用可控密度构造（stripe_ctrl_alpha）来确保 **|J| 可控**，从而让阶段解释更“干净”。

- 密度定义使用总纲约定：

$$
\alpha=\frac{|J|}{n_r+n_s}
$$

（强调：不是 $|R||S|$ 归一）

- 典型默认：
  - `NR=NS=100000`
  - `T=100000`
  - `REPEATS=3`
  - `THREADS=1`（公平性：单线程）
- 代表性密度点（默认）：
  - `ALPHAS="0.1 3 30"`（稀疏 / 中等 / 稠密）

> 推荐增强（更能解释 adaptive 拐点）：增加阈值附近点，例如 `ALPHAS="0.1 3 5 30"`（若 $J_\star=1e6$ 且 $n_r+n_s=2e5$，拐点约在 α≈5）。

------

## 3. 方法集合与运行模式（method × variant 矩阵）

### 3.1 方法集合（METHODS）

EXP-7 支持横跨多类基线（索引采样、分区 join、范围查询采样、拒绝采样、映射索引等），默认：

```
METHODS="ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection tsunami"
```

### 3.2 统一三种运行模式（Variants）

EXP-7 固定沿用论文推荐的统一三模式：

1. **sampling**：不显式物化 $J$，直接采样 t 次
2. **enum_sampling**：先枚举全 $J$，再从列表中有放回均匀抽样
3. **adaptive**：先 pilot 枚举探测；若 |J| 小走枚举快；若 |J| 大走采样稳（由阈值 `J_STAR` 控制）

EXP-7 的默认策略是：

- 对所有 α：跑 `sampling + adaptive`
- 仅在稀疏 α：可选跑 `enum_sampling`（避免稠密点枚举爆炸）

### 3.3 重要限制：rejection + sampling 在 sweep 中会被禁用

在 `sjs_sweep` 的大规模 sweep harness 中，**Rejection + Sampling 组合被安全禁用**（会写入 `ok=0` 且 note=SKIPPED）。因此：

- EXP-7 的 sampling+adaptive sweep 里**不会出现** `rejection/sampling` 的有效结果；
- `rejection/adaptive` 仍会跑（建议用它代表 rejection 类方法）。

> 写论文时必须在图注/实验设置里显式说明，否则读者可能以为漏跑或失败。

------

## 4. 计时口径与阶段定义（最关键的“读懂条件”）

### 4.1 两种时间：wall_ms vs phases_json

EXP-7 同时保存两种时间口径：

- `wall_ms`：单次 run 的端到端 wall time（覆盖 runner 的整个执行）
- `phases_json`：细粒度阶段计时（runner 内用 PhaseRecorder 记录）

EXP-7 的阶段图/阶段占比 **默认使用 phases_json 中所有 `run_\*_ms` 的求和**作为 100% 基准（便于“拆解解释”）。

> 注意：`wall_ms` 通常会略大于 `sum(run_*_ms)`，因为有少量开销不在任何 `run_*` scope 内（如 Reset、样本校验 Validate、以及少量框架开销）。
>  建议论文表述：**phase breakdown 用 run_\*；wall_ms 作为 sanity check 与主线性能指标。**

### 4.2 各 variant 的阶段归类（EXP-7 使用的规则）

EXP-7 在后处理脚本里将 runner 的 phase key 映射到四大类：

#### Variant = sampling

- Build = `run_build_ms`
- Count = `run_count_ms`
- Sample = `run_sample_ms`
- Enumerate = 0

#### Variant = enum_sampling

- Build = `run_build_ms`
- Enumerate = `run_enum_prepare_ms + run_enum_pass1_count_ms`
- Sample = `run_draw_ranks_ms + run_enum_pass2_select_ms`
- Count = 0

> 说明：第二遍扫描（pass2_select）严格说是“为抽样服务的选择过程”，本实验将其算在 Sample 类。若你更想强调“双遍枚举”，也可把 pass2_select 归到 Enumerate，但要在图注写清口径。

#### Variant = adaptive

- Build = `run_build_ms`
- Enumerate（pilot） = `run_pilot_enum_prepare_ms + run_pilot_enum_scan_ms`
- Count（fallback） = `run_fallback_count_ms`
- Sample = `run_fallback_sample_ms + run_small_join_sample_from_list_ms`

> Adaptive 有两种分支：
>
> - **enumerate_all**：pilot 直接完成（fallback_count/sample ≈ 0）
> - **fallback_sampling**：pilot 触发阈值后回退（fallback_count/sample 显著 > 0）

------

## 5. 输出文件结构与字段说明（怎么读 raw / summary / breakdown）

EXP-7 脚本会输出一个主目录（默认 `results/sweeps/exp7_profile/`），主要包含：

### 5.1 Sweep 结果目录

#### A) `sampling_adaptive/`

- `sweep_config.json`：本次 sweep 的完整配置（可复现）
- `sweep_raw.csv`：**每次重复一行**（最重要的证据源）
- `sweep_summary.csv`：按 (alpha, method, variant, t, N) 聚合后的 summary（median/mean/p95 等）

#### B) `enum_sparse/`（可选）

- 仅在 `INCLUDE_ENUM_SPARSE=1` 时生成，用于稀疏点的 enum_sampling 对照

### 5.2 EXP-7 后处理输出（解释型产物）

- `exp7_merged_sweep_raw.csv`：把多个 raw 合并（方便统一分析）
- `exp7_breakdown_median.csv`：阶段拆解的核心表（每个配置点一个 median 记录）
- `exp7_breakdown_median.md`：便于直接贴论文/README 的 Markdown 表
- `figs/exp7_phase_breakdown_alpha_*.png`：每个 alpha 一张堆叠柱状图

### 5.3 sweep_raw.csv 关键字段（建议读者优先检查）

- `ok`：该次 run 是否成功（必须为 1 才纳入 breakdown）
- `note`：失败原因（如 SKIPPED、OOM、enum_truncated 等）
- `wall_ms`：端到端时间
- `phases_json`：阶段计时 JSON（EXP-7 分解的来源）
- `count_value` / `count_exact`：计数值及是否 exact（用于 sanity 与解释）
- `used_enumeration`：是否走过枚举流程
- `enum_truncated` / `enum_cap`：枚举是否截断与容量上限（失败处理依据）
- `adaptive_branch`、`adaptive_pilot_pairs`：adaptive 分支信息（解释拐点必看）

------

## 6. 如何阅读 EXP-7 的表与图（建议按这个顺序）

### Step 1：先确认语义一致性与无静默失败

对 `sweep_raw.csv` / merged raw：

1. 检查 `ok==1` 的比例是否接近 100%（除预期 SKIPPED 组合外）
2. 检查是否存在 `enum_truncated==1`（若有，需要在图中标注或剔除并说明）
3. 对 controlled stripe 数据，检查 `count_exact` 是否大多为 1（多数方法应能 exact count；某些估计型/拒绝采样型可能是 approx，需要解释口径）

### Step 2：读 adaptive 的“分支行为”

在同一 alpha 下，查看 `adaptive_branch` 的分布：

- 稀疏点应多为 `enumerate_all`
- 稠密点应多为 `fallback_sampling`
- pilot pairs 通常接近阈值 `J_STAR`（用于说明触发机制）

这一步是回答 RQ6 的“机制证据”。

### Step 3：读堆叠柱图（或 breakdown 表）解释“谁在什么阶段赢”

在固定 alpha 下比较不同方法：

- Build% 高 → 构建重（适合高 t 场景摊薄；或适合重复多 query）
- Count% 高 → 计数/权重计算重（说明瓶颈在统计结构，而不是抽样本身）
- Enumerate% 高 → 枚举/扫描重（典型出现在小 join 或 enum_sampling）
- Sample% 高 → 采样定位/二次扫描主导（随 t 增长更敏感）

建议论文叙述模板（可直接复用）：

> “在 α=… 时，方法 X 的总时间主要由 Count/Sample 主导，这说明其瓶颈在 …；方法 Y 的时间主要消耗在 Build，说明 …；因此在 … 场景下 X/Y 更占优。”

### Step 4：将 EXP-7 与 EXP-2/3/6 的趋势做“闭环解释”

EXP-7 本质是“解释器”：

- EXP-2 看到随 t 的斜率差异 → 用 EXP-7 的 Sample%/Count% 解释
- EXP-3 看到随 α 的退化 → 用 EXP-7 的 Enumerate%（爆炸）或 Count%（退化）解释
- EXP-6 看到 adaptive ratio 接近 1 → 用 EXP-7 的分支与阶段结构解释其接近原因

------

## 7. 失败处理与公平性声明（必须写进论文/README）

### 7.1 枚举内存风险与截断

- 通过 `enum_cap` 统一安全上限（超过即 `enum_truncated=1`）
- 图表中应显式标注失败点（不可静默忽略）

### 7.2 禁用组合（必须声明）

- `rejection + sampling` 在 sweep harness 中被禁用（SKIPPED）
   → EXP-7 只展示 rejection 的 adaptive（或 enum_sampling sparse）结果。

------

## 8. 复现指南（面向读者/审稿人）

### 8.1 一键运行

在仓库根目录：

```
bash run/run_exp7.sh
```

### 8.2 常用覆盖参数（环境变量）

- 代表点选择：

```
ALPHAS="0.1 3 5 30" bash run/run_exp7.sh
```

- 改样本数与重复次数：

```
T=200000 REPEATS=5 bash run/run_exp7.sh
```

- 只跑少数方法（快速验证）：

```
METHODS="ours r_tree pbsm" bash run/run_exp7.sh
```

- 关闭 enum_sampling sparse（更省时间）：

```
INCLUDE_ENUM_SPARSE=0 bash run/run_exp7.sh
```

### 8.3 输出在哪里

默认输出：

- `results/sweeps/exp7_profile/sampling_adaptive/sweep_raw.csv`
- `results/sweeps/exp7_profile/exp7_breakdown_median.md`
- `results/sweeps/exp7_profile/figs/exp7_phase_breakdown_alpha_*.png`

------

## 9. 论文呈现建议（建议放在实验章末尾增强可信度）

### 9.1 推荐主文呈现形式

- **一张表**（`exp7_breakdown_median.md` 的精简版）
  - columns：method、variant、Build%、Count%、Enum%、Sample%、total_ms、wall_ms
- **两到三张堆叠柱图**（按 α 分图）

### 9.2 图注/表注必须包含的口径句

建议直接写成固定句式（避免被 reviewer 挑口径）：

1. “阶段占比以 phases_json 中所有 `run_*` 计时求和为 100%，wall_ms 另行报告作 sanity check。”
2. “数据生成与文件读写不计入阶段时间；仅统计方法内部 Build/Count/Enumerate/Sample。”
3. “Rejection baseline 在 sweep harness 中不运行 sampling 变体（SKIPPED），因此仅展示 adaptive 分支。”

------

## 10. 一页式读者检查清单（建议放在附录/README）

读 EXP-7 时，建议顺序：

1. **raw.csv**：ok 是否全为 1（除 SKIPPED）？有没有 enum_truncated？
2. **adaptive_branch**：稀疏是否 enumerate_all，稠密是否 fallback_sampling？pilot_pairs 是否接近 J_STAR？
3. **count_value/count_exact**：controlled 数据是否自洽（|J| 与 α·(n_r+n_s) 对齐）？
4. **breakdown 表/图**：谁的瓶颈在 Build/Count/Enum/Sample？是否能解释 EXP-2/3/6 的趋势？
5. **wall_ms vs run_total**：差异是否稳定（说明有固定 overhead）？是否需要在图里加 “Other”？