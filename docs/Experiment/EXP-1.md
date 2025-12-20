# 实验 1 阅读文档：Correctness & Sample Quality（对应 RQ1）

> 本文档用于“读懂并复现解读”**EXP‑1**：它回答 RQ1——所有方法在统一语义下，是否真的输出 **对 Join 结果集 $J$ 的 i.i.d. 均匀（有放回）样本**，并能在统计意义上通过质量检验。

------

## 1. 实验定位与目标

### 1.1 实验在整篇研究中的位置

实验总纲将整章实验按 RQ1–RQ6 组织，其中 **EXP‑1**专门验证 **正确性与采样质量**，是后续性能实验（runtime、scalability、adaptivity）可信度的“地基”。

### 1.2 EXP‑1 的核心目标

在**小规模**数据（可计算 oracle 真值）下，对每个方法、每种运行模式、多个随机种子：

1. **Correctness**：每个输出样本 $(r,s)$ 必须属于真实 join 结果集 $J$。
2. **Uniformity**：样本在 $J$ 上应近似均匀（统计意义上不拒绝均匀假设）。
3. **Independence**：样本间应近似独立（自相关接近 0）。

最终产物是：**采样质量指标表**（missing、χ²/KS、自相关）。

------

## 2. 任务定义与统一语义（必须读懂）

### 2.1 数据与 Join 定义

- 两个关系 $R, S$，元素是 **2D 轴对齐矩形盒**，几何语义是 **half-open**：
  $$
  r=[L_x(r),R_x(r))\times[L_y(r),R_y(r)),\quad
  s=[L_x(s),R_x(s))\times[L_y(s),R_y(s)).
  $$

- 空间连接结果：
  $$
  J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq\emptyset\}
  $$

- 注意：**边界贴合不算相交**（half-open 的关键）。

### 2.2 采样任务（EXP‑1 验证的对象）

在给定工作负载下，从 $J$ 中抽取 **$t$ 个 i.i.d. 均匀、有放回样本**。也即每个样本都是 join pair 的标识，且每次抽样彼此独立、同分布、概率为 $1/|J|$。

------

## 3. 方法集合与运行模式（读表时的维度）

### 3.1 方法集合（methods）

总纲建议覆盖多类代表性 baseline，包括索引类、分区 join 类、拒绝采样类、学习型索引类等。
 仓库当前 build（Dim=2）支持的方法列表可在 README 中看到（如 ours、aabb、interval_tree、kd_tree、r_tree、range_tree、pbsm、tlsop、sirs、rejection、tsunami）。README

### 3.2 三种统一运行模式（variants）

实验总纲要求每种方法都以三种模式评估，以保证语义与公平性一致：

1. **Sampling**：不物化全部 $J$，直接输出 $t$ 个均匀样本
2. **Enum+Sampling**：先枚举全部 $J$，再在枚举结果上均匀有放回抽样
3. **Adaptive**：当 $J$ 小时走枚举更快，否则走采样更稳（阈值/探测决定分支）

> 阅读 EXP‑1 表格时，你应当始终以 “method × variant × seed” 为最小单元看结果。

------

## 4. EXP‑1 的数据设置（为何必须小规模）

### 4.1 为什么必须可算 oracle

EXP‑1 的 correctness 与均匀性检验都依赖 “真实 universe $J$”：

- **missing-in-universe**：输出样本中不属于 $J$ 的比例（必须为 0）
- **χ²**：需要知道每个 join pair 的真实类别集合（或至少能构造可信分桶）

因此 EXP‑1 必须选择小规模（例如 $n_r,n_s$ 在几千量级），确保 oracle 计算可行。

### 4.2 推荐使用的合成生成器

总纲推荐 EXP‑1 使用合成数据，最“paper friendly”的是**可控密度条带构造**（SYN‑CTRL）：能精确控制 $|J|$ 从而精确控制密度参数 $\alpha$。

并且总纲强调密度参数定义为：
$$
\alpha=\frac{|J|}{n_r+n_s}
$$
不是 $|J|/(|R||S|)$。这点在读表/写论文时必须写清楚。

------

## 5. EXP‑1 的执行流程（读懂它产出的每个数字来自哪里）

> 仓库提供 `sjs_verify`，用于**小规模数据的 oracle checks + sampling-quality tests**，是 EXP‑1 的推荐执行入口。README

EXP‑1 的逻辑流程可理解为：

1. **生成/加载数据**（合成或二进制输入；数据生成/IO 不计入主要性能，但 EXP‑1 主要关心质量）
2. **Oracle 枚举**得到真实 $J$ 或足够完整的 universe 表示
3. **运行方法**（指定 method + variant + seed），输出 $t$ 个样本
4. **质量评估**：对样本与 oracle universe 做 correctness 检查与统计检验
5. **重复**：多个 seed / repeats，汇总 median / p95（或误差条）

------

## 6. 输出指标解释（如何读 EXP‑1 质量表）

总纲明确 EXP‑1 至少包含：Oracle $|J|$、missing、χ²/KS p-value、自相关。
 下面给出“指标含义 + 推荐解读方式”。

### 6.1 Oracle $|J|$

- 含义：小规模下用朴素或可靠枚举得到的真实 join size。
- 用途：
  - 为 correctness/χ² 提供 universe
  - 帮你判断 $t/|J|$ 是否足够大（χ² 需要每类期望计数不太小）

### 6.2 missing-in-universe（正确性硬指标）

- 定义：输出样本中“不属于 $J$”的比例。
- 读表规则：
  - **必须为 0**。只要 >0，RQ1 直接不通过（语义不一致）。

### 6.3 χ² p-value（均匀性主证据）

- 含义：在 “样本类别分布 = 对 $J$ 的均匀分布” 这一零假设下，χ² 检验给出的 p-value。
- 推荐解读：
  - p-value 并不是越大越好；它表达“是否拒绝均匀假设”。
  - **主文建议报告 median（跨 seed）**，并报告失败率（多少 seed < 0.05），同时说明多重检验（多个方法/seed）下的偶然拒绝。

> 实务注意：若 $|J|$ 过大导致每类期望计数很小，χ² 会不稳定；此时应采用分桶策略并写清楚（总纲也提示需要分桶与多重检验控制）。

### 6.4 KS p-value（分布 sanity check）

- 总纲建议：将 pair 标识哈希到 $[0,1)$ 后做 KS 检验并报告 p-value。
- 推荐解读方式（写进论文最稳）：
  - 把 KS 定位为 **sanity check**，与 χ² 一起提供“近似均匀”的证据；
  - 如果 $|J|$ 很小、且 hash 映射导致分布支撑离散，KS 对连续均匀的零假设可能出现“系统性拒绝”。此时应解释并配合 χ²/分桶策略呈现（仍以 χ² 为主证据）。

### 6.5 autocorrelation（独立性 sanity check）

- 总纲建议用自相关检验样本间独立性：应接近 0。
- 推荐读法：
  - 报告 `lag=1`（以及可选 lag=2/5）自相关的 median 与 max(|·|)
  - 只要数值很小（接近 0），可作为 “i.i.d.” 的 sanity 支撑。

------

## 7. 通过/失败判定（EXP‑1 结论如何写）

### 7.1 通过标准（建议写法）

- **Correctness**：所有方法、所有模式、所有 seed，missing=0。
- **Uniformity**：χ² 在统计意义上不拒绝均匀假设（p-value 大多 ≥0.05，少量拒绝可归因于多重检验/随机波动，并在附录给出完整 seed 结果）。
- **Independence**：自相关接近 0（sanity check）。

### 7.2 失败处理（总纲要求显式标注）

虽然 EXP‑1 一般不会触发“枚举截断”，但总纲强调：任何失败点都必须显式记录，不可静默忽略。
 因此，如果 `oracle_collect_limit` 或 `enum_cap` 导致 universe 不完整，应在日志/表格中显式写明 “quality skipped / universe truncated”。

------

## 8. 推荐的主文呈现格式（你们论文里该怎么放）

总纲建议 EXP‑1 主文给一张表（方法 × 指标），附录给更全重复结果。
 推荐表格列（最稳健、读者一眼懂）：

- `missing-in-universe`（必须 0）
- `chi2_p`（median across seeds）
- `ks_p`（median across seeds，或附录）
- `autocorr_lag1`（median / max）
- （可选）`oracle |J|`（给读者尺度感）

并在表注中写清楚：

- half-open 语义；
- 采样目标为 i.i.d uniform with replacement；
- seed 重复次数、汇总口径用 median（可加 p95/误差条）。

------

## 9. 一段可直接放进论文的“EXP‑1 描述模板”

> 我们在可计算 oracle 的小规模数据上验证所有方法输出语义一致性与采样质量（RQ1）。任务定义为二维 half‑open 轴对齐矩形的空间连接，目标是从真实 join 结果集 $J$ 中抽取 $t$ 个 i.i.d. 均匀、有放回样本。
>  对每种方法及其三种运行模式（Sampling、Enum+Sampling、Adaptive），我们在多个随机种子下重复运行并对输出样本执行 correctness 检验（missing‑in‑universe 必须为 0）以及均匀性与独立性 sanity checks（χ²/KS p-value、自相关）。我们报告跨 seed 的 median，并在附录提供完整重复结果与统计检验细节。