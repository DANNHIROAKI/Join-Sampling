# EXP-5 阅读文档：分布鲁棒性（RQ5）

本阅读文档用于帮助你（以及审稿人/复现者）快速理解 **EXP-5（RQ5: Robustness to distribution shift）** 的实验目的、配置、输出文件含义、以及如何从结果中抽取“可写进论文”的结论。

> 对齐总纲：任务是对空间连接结果集 \(J\) 做 **对 \(J\) 的 i.i.d. 均匀、有放回抽样**（half-open 语义下的二维 AABB 相交）。

---

## 1. EXP-5 的定位与要回答的问题

### 1.1 对应研究问题

- **RQ5：Robustness to distribution shift**
  - 当数据分布从“可控密度构造（SYN-CTRL）”切换到更一般的合成分布（均匀/聚类偏斜/大小混合）后：
    1) RQ2（runtime vs \(t\)）的趋势是否仍成立？
    2) 哪些 baseline 对偏斜或大框混合更敏感？

### 1.2 EXP-5 的核心叙事（写论文时建议固定这 3 句）

1. **目的**：检验不同方法的性能趋势是否对“分布变化”稳健。
2. **设计**：固定规模 \(n_r=n_s\)，切换分布（Uniform / Clustered / Hetero-sizes），扫描样本量 \(t\)（1k→1M），比较 sampling 与 adaptive 两种管线。
3. **结论期待**：
   - 采样类方法应主要随 \(t\) 增长；
   - 拒绝采样类应在“大框混合/高重叠”场景出现明显退化（acceptance 下降）；
   - 自适应应尽量逼近 \(\min\{T(\text{sampling}),T(\text{enum+sampling})\}\) 的下界（但 EXP-5 默认不跑 enum+sampling）。

---

## 2. 实验设计（EXP-5 具体做了什么）

### 2.1 数据分布（SYN-U / SYN-C / SYN-H）

EXP-5 使用 3 类合成生成器（均为 2D）：

- **(SYN-U) Uniform：均匀随机矩形**
  - 通过宽度范围控制重叠概率。
  - 默认参数（来自 `run/run_exp5.sh`）：
    - `w_min=0.005, w_max=0.02`
    - `same_size_all_dims=false`
    - `shuffle_r=false, shuffle_s=false`

- **(SYN-C) Clustered：聚类/热点偏斜**
  - 通过热点簇制造强偏斜与局部高密度。
  - 默认参数：
    - `clusters=10, sigma=0.05, share_clusters=true`
    - `w_min=0.003, w_max=0.02`
    - `shuffle_r=false, shuffle_s=false`

- **(SYN-H) Hetero-sizes：大小混合（少量大框 + 大量小框）**
  - 通过大框比例与大框尺寸制造“活跃集/剪枝退化”的极端场景。
  - 默认参数：
    - `anisotropic=true`
    - `p_large=0.1`
    - `w_small_min=0.002, w_small_max=0.01`
    - `w_large_min=0.1, w_large_max=0.3`
    - `shuffle_r=false, shuffle_s=false`

> 重要说明（必须写进论文/读者须知）：
>
> - EXP-5 的 `alpha=1.0` 仅作为配置字段存在。
> - **对于 Uniform/Clustered/Hetero 这三类 generator，alpha 通常不意味着“真实 \(|J|\)”被控制在某个目标值**。
> - 因此在解释 EXP-5 时，应同时报告每个分布下的 **实际 \(|J|\)**（可用 `count_mean`/exact count 给出），避免把“密度变化”误当成纯“形态变化”。

### 2.2 Sweep 维度与方法集合

- **规模**：默认 `NR=100000, NS=100000`（可通过环境变量覆盖）。
- **样本量 sweep**：
  - `t ∈ {1000, 3000, 10000, 30000, 100000, 300000, 1000000}`
- **方法列表（默认）**：
  - `ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection`
- **运行模式（EXP-5 只跑两种）**：
  1) `sampling`
  2) `adaptive`

> EXP-5 **不跑** `enum_sampling`，并将 `enum_cap=0`：这是为了避免枚举在 \(|J|\) 较大时的内存风险。

### 2.3 随机性与重复次数

- `REPEATS=3`：每个点位（method×variant×t）重复 3 次
- 运行种子：`RUN_SEED`（默认 1），重复时按 rep 偏移
- 生成种子：`GEN_SEED`（默认 1）

> 对“鲁棒性”更强的写法（可选增强）：每个分布使用多个 `GEN_SEED`（例如 1/2/3），跨 seed 聚合 median/p95。

### 2.4 自适应阈值（adaptive）

- `J_STAR=1000000`：adaptive 模式下用于选择“枚举分支/采样分支”的阈值。

---

## 3. 如何运行与复现

### 3.1 一键运行

```bash
bash run/run_exp5.sh
```

### 3.2 常用覆盖参数（环境变量）

- 改规模：
  - `EXP5_NR=200000 EXP5_NS=200000`
- 改重复次数：
  - `EXP5_REPEATS=5`
- 改线程数：
  - `EXP5_THREADS=1`（主文建议固定 1）
- 改方法集合：
  - `EXP5_METHODS="ours r_tree tlsop"`
- 改 t 列表：
  - `EXP5_T_LIST="1000 10000 100000 1000000"`

### 3.3 输出目录结构

运行结束后，结果会写到：

```
results/sweeps/exp5/<TAG>/
  uniform_t/
    sweep_raw.csv
    sweep_summary.csv
    runtime_t.png
    sweep.log
    sweep_config.json
  clustered_t/
    ...
  hetero_t/
    ...
  configs/
    exp5_t_uniform.json
    exp5_t_clustered.json
    exp5_t_hetero.json
```

- `configs/` 目录保存了本次 sweep 使用的完整配置（强烈建议随论文 artifact 一起打包）。

---

## 4. 输出文件解读指南

### 4.1 `sweep_summary.csv`（写论文主要看这个）

每行对应一个（dataset×method×variant×t）配置点的聚合统计。

关键列含义：

- `method`, `variant`, `t`：方法、模式、样本量
- `ok_rate`：重复运行中成功完成的比例（1.0 为满成功）
- `wall_median_ms`, `wall_p95_ms`：端到端 wall time 的 median / p95（单位 ms）
- `count_mean`, `exact_frac`：连接规模 \(|J|\) 的统计（若方法提供 count；`exact_frac=1` 通常表示 exact count）
- `note`：失败/跳过原因（例如某些组合被框架显式禁用）

**论文中建议至少报告：**

- 每个分布一张表/一段文字：`NR, NS, count_mean(|J|)`
- 每个分布一张 `runtime vs t` 曲线（log-x，建议 log-y）

### 4.2 `sweep_raw.csv`（需要排查异常/做更细分析时看）

每行对应一次具体运行（含 rep）。常用列：

- `ok`：该次运行是否成功
- `wall_ms`：该次运行端到端 wall time
- `count_value`, `count_exact`：该次运行的 count
- `adaptive_branch`：adaptive 选择的分支（如有）
- `enum_truncated`, `enum_cap`：枚举截断信息（EXP-5 通常不涉及）
- `note`：失败/跳过原因

> 注意：sweep 输出通常以“压缩可读”为优先，某些 JSON 字段可能被截断；如需完整阶段分解，建议使用 EXP-7 的 `sjs_run`/profiling 流程。

### 4.3 `runtime_t.png`

- 横轴：\(t\)（log）
- 纵轴：`wall_median_ms`（log）
- 每条曲线：`method-variant`

**读图要点：**

- 曲线整体随 \(t\) 增长是否近似线性（log-log 上斜率接近 1）
- 对比不同分布下同一方法曲线的“整体平移/形状变化”是否明显
- 是否存在“某分布突然崩溃/爆炸”的方法（通常对应偏斜/大框导致剪枝失效或 rejection acceptance 下降）

---

## 5. 如何从 EXP-5 写出一段“有说服力”的论文结论

### 5.1 建议的结果组织方式（主文 1 张图 + 1 张表即可）

- **表（dataset stats）**：三行（Uniform/Clustered/Hetero），列出：
  - `n_r`, `n_s`
  - 实际 \(|J|\)（用 `count_mean`，优先选 exact 方法的 count）
- **图（runtime vs t）**：
  - 方案 A：三张子图（每个分布一张）
  - 方案 B：主文放最能区分差异的一张（通常是 Hetero），其余放附录

### 5.2 解释口径（避免审稿人挑刺的关键句）

- **必须说明**：EXP-5 的三类分布 **并不保证相同 \(|J|\)**，因此性能差异来源包括：
  1) 分布形态变化（偏斜、大框混合）
  2) 由形态变化诱发的连接密度变化（\(|J|\)、平均度、active set 大小）

- **对 rejection 类方法的解释模板**：
  - “在大小混合分布下，大框显著扩大候选区域，使得 rejection proposal 的上界变松/acceptance 降低，导致 runtime 爆炸；这符合 rejection 类方法对分布 shift 的已知敏感性。”

- **对索引/分区类 baseline 的解释模板**：
  - “偏斜会导致空间局部拥挤或分区负载不均，削弱剪枝与缓存局部性；大框会提高节点 overlap，削弱树结构 pruning。”

### 5.3 一个很实用的“鲁棒性指标”（可选写进附录）

对同一方法（固定 variant，固定 t），展示跨分布的相对波动：

\[
\text{robustness\_factor} = \frac{\max_{dist} T(dist)}{\min_{dist} T(dist)}.
\]

- 越接近 1 越稳健。

---

## 6. EXP-5 的常见坑与处理建议

### 6.1 `rejection + sampling` 被跳过

- 这是框架的“安全禁用”行为（会在 `sweep.log` 和 `note` 里写 `SKIPPED`）。
- 读图/写论文时要明确：rejection 只比较 adaptive 变体。

### 6.2 `ok_rate` 的过滤标准

- 作图与汇总时，建议只纳入 `ok_rate == 1.0` 的点位。
- 若纳入部分成功点位，必须保证 wall time 的统计只对成功 run 进行聚合，否则可能出现“失败很快→median 被拉低”的假象。

### 6.3 PBSM/TLSOP 等参数的代表性

- EXP-5 默认沿用了某些配置（例如 PBSM 的 scheme），其“典型最优参数”可能与条带构造（SYN-CTRL）不同。
- 若论文希望对 PBSM/TLSOP 给出更公平的对比，建议在附录做一个小的参数敏感性实验（或至少说明采用的是统一默认参数）。

### 6.4 生成种子只有一个（GEN_SEED=1）

- 这会让“robustness”的统计意义偏弱（更像 case study）。
- 若要更硬：每个分布至少 3 个 seed，跨 seed 聚合。

---

## 7. 推荐的“最小增强版 EXP-5”（如果你有预算补实验）

1. 每个分布取 3 个 `GEN_SEED`：`1,2,3`，每个 seed 跑 `REPEATS=1`
2. 输出除 runtime 外，再加一张表：`|J|`（exact count）与 `avg_degree = |J|/(n_r+n_s)`
3. 图中用不同线型或分面区分 seed（或先在 summary 聚合后再画）

这样 EXP-5 的结论会从“看起来合理”提升到“统计意义上更难被挑刺”。

---

## 8. 一段可直接放进论文的 EXP-5 描述模板（可复制）

> **EXP-5 (RQ5: Distribution Robustness).** We evaluate distribution shift robustness by replacing the controllable-density generator with three general synthetic distributions: uniform rectangles (SYN-U), clustered hotspots (SYN-C), and heterogeneous sizes with a fraction of large rectangles (SYN-H). For each distribution, we fix \(n_r=n_s=100k\) and sweep the sample size \(t\in[10^3,10^6]\). We report the median and p95 wall-clock time over 3 runs under a single-thread setting. Since these generators do not enforce a fixed join size \(|J|\), we additionally report the achieved \(|J|\) (exact when available) to contextualize performance differences.