# 03_baselines.md

本文件说明：

1) 本项目中 baseline 的统一接口（保证实验公平）  
2) 每个 baseline 强制实现的三种版本（Sampling / Enum+Sampling / Adaptive）  
3) 已实现 baseline 的方法名与可调参数（`Config.run.extra`）

---

## 1) Baseline 的统一接口

公共接口见 `include/sjs/baselines/baseline_api.h`。

每个 baseline 都要实现：

- `Build(ds, cfg, phases)`：构建索引/预处理
- `Count(cfg, rng, out, phases)`：计数 \(|J|\)（精确或估计）
- `Sample(cfg, rng, out, phases)`：抽样 \(t\) 个 join pair（期望 i.i.d. 均匀）
- `Enumerate(cfg, phases)`：返回确定性枚举器（用于 Enum+Sampling 与小规模验证）
- `Reset()`：repeat 间重置内部状态

输出统一结构：

- `CountResult`：`value/exact/stderr/ci` 等
- `SampleSet`：`pairs` + 可选 `weights`

---

## 2) 三种版本（variant）与 Runner

为了“跑法一致”，项目把三种 variant 的 orchestration 固化在 `include/sjs/baselines/runners/`：

- `sampling_runner.h`：强制走 Sampling 流程（Phase1/2/3 或 baseline 等价流程）
- `enum_sampling_runner.h`：强制先枚举 join，再 rank/alias 采样
- `adaptive_runner.h`：阈值切换逻辑（先 count，再判断是否保留枚举）

baseline 只需要保证自己的 `Count/Sample/Enumerate` 语义正确；流程与计时由 runner 统一处理。

---

## 3) method 名称与 CLI

通过 `--method=<name>` 选择方法，常见值：

- `ours`
- `interval_tree`
- `aabb`
- `kd_tree`
- `r_tree`
- `range_tree`
- `pbsm`
- `tlsop`
- `sirs`
- `rejection`
- `tsunami`（若已注册到工厂/CLI）

通过 `--variant=<sampling|enum_sampling|adaptive>` 选择版本。

> method/variant 的解析与工厂注册在 `src/baselines/baseline_factory_2d.*` 与 `src/baselines/baseline_names.cpp`。

---

## 4) 各 baseline 的实现口径（“你在论文里怎么写”）

### 4.1 ours

**核心思想（Dim=2）**：

- 在 x 维 sweep：把 join pair 分解为“唯一 START 事件归属”
- 在 y 维把相交分成 A/B 两类（由 lower endpoint 的相对关系定义）
- 两遍 sweep：
  - Phase1：对每个 START 事件统计 `w_e^A, w_e^B`
  - Phase2：对事件按 `w_e` 建 alias，给每个样本位置分配 (e, pattern)
  - Phase3：重扫并在对应集合上做局部均匀采样

工程实现使用：
- Fenwick/BIT + buckets 支持 count / sample（A/B 两种模式）

### 4.2 interval_tree

典型 baseline：在 sweep 过程中维护 y 维 interval tree（stabbing 查询），用于快速获取与新加入矩形相交的 active-set 元素。

- Sampling：用计数 + 局部采样组合
- Enum+Sampling：Report 全部相交对
- Adaptive：阈值切换

### 4.3 aabb

AABB-tree（BVH）类 baseline：

- Build：构建两棵 AABB tree
- Enumerate：树-树遍历输出所有相交对
- Sampling：以节点对的覆盖面积/候选数为权重做抽样（实现上在 baseline 内部做等价流程）

### 4.4 kd_tree / range_tree / sirs

这类 baseline 往往依赖 **embedding（join→range）**：

- 将矩形映射为高维点（或 lower/upper 的组合）
- 将相交判定转化为 stabbing / orthant range query
- 用 KD-tree / RangeTree / SIRS 风格的结构做计数与采样

> 目前工程里保留 `include/sjs/geometry/embedding.h` 来统一该映射接口，便于高维扩展。

### 4.5 r_tree

经典 R-tree baseline：

- Build：建 R-tree
- Enumerate：R-tree join（输出全部相交对）
- Sampling：在 R-tree traversal 中用候选对数估计做抽样规划

### 4.6 pbsm / tlsop / tsunami

这三类 baseline 偏“系统/并行 join”方向，常见做法是 **partition + local join**：

- PBSM：按网格/条带分区，把矩形分配到分区，再对分区内做局部 join
- TLSOP：使用排序/扫描的 pipeline 优化与分区复用
- Tsunami：更激进的分区与流水化策略（项目里以 baseline 形式复现）

在本框架中它们仍需支持三种 variant，只是内部实现与纯索引型 baseline 不同。

### 4.7 rejection（重要性/拒绝采样类）

该 baseline 常见于“用易采样的 proposal 分布”近似均匀抽 join pair：

- Count：给出 \(|J|\) 的估计或上界
- Sample：从 proposal 中抽候选，再用 rejection 校正（需确保最终分布正确）

---

## 5) baseline-specific 参数（Config.run.extra）

`Config.run.extra` 是一个 `map<string,string>`，用于在不破坏 CLI/JSON 的前提下扩展参数。

`config/default.json` 给了常用键的示例（你可以直接在 CLI 传 `--key=value` 覆盖）。

### 5.1 PBSM

- `pbsm_scheme`：`stripes` / `grid`（示例）
- `pbsm_k`：分区数（0 表示自动）
- `pbsm_part_axis`：分区轴（0/1）
- `pbsm_sweep`：`orthogonal` 等

### 5.2 TLSOP

- `tlsop_nx`, `tlsop_ny`：网格大小
- `tlsop_reuse_sort`：是否复用排序结果

### 5.3 SIRS

- `sirs_outer`：外层结构选择（留空用默认）
- `sirs_leaf_size`：叶子阈值

### 5.4 Rejection

- `rej_index`：以 R 或 S 为 proposal
- `rej_rep_center`：是否用中心点代表
- `rej_count_draws`：计数/校正抽样次数
- `rej_max_factor`：上界松弛系数

---

## 6) 新增 baseline 的工程步骤

1. 在 `include/sjs/baselines/<name>/` 写 3 个头文件：
   - `sampling.h`
   - `enum_sampling.h`
   - `adaptive.h`

2. 在 `src/baselines/<name>/` 写 3 个 `.cpp` 做 **Dim=2 显式实例化**。

3. 注册到工厂（Dim=2）：
   - `include/sjs/core/types.h`：`enum class Method` 增加新方法
   - `src/baselines/baseline_factory_2d.cpp`：switch 增加 case
   - `src/baselines/baseline_names.cpp`：加入 registry（help 列表会自动包含）

4. 增加 tests（至少 smoke test），并用 `sjs_verify` 做小规模 oracle 对比。
