# EXP-6 阅读文档

## 自适应策略在 $(\alpha,t)$ 网格下的“接近最优”验证（RQ6）

### 对应研究问题（RQ6）

**RQ6：Adaptivity effectiveness**
 自适应策略（Adaptive）能否在不同的 $(\alpha,t)$ 场景下自动选择接近最优的分支（枚举 or 采样），并使得端到端时间接近：
$$
\min\{T(\text{sampling}),\ T(\text{enum+sampling})\}.
$$

------

## 1. 背景与任务定义（与全文一致）

给定两集合 $R,S$，元素是二维轴对齐矩形（**half-open** 语义）：
$$
r=[L_x(r),R_x(r))\times[L_y(r),R_y(r)),\quad
s=[L_x(s),R_x(s))\times[L_y(s),R_y(s)).
$$
空间连接结果为：
$$
J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq\emptyset\},
$$
相交严格遵循 half-open：**边界贴合不算相交**。

EXP-6 的比较目标始终是：返回 **$t$ 个来自 $J$ 的 i.i.d. 均匀、有放回样本**。

------

## 2. EXP-6 的核心目的：把“自适应是否接近最优”画成一张相图

EXP-6 不追求覆盖所有数据分布（那是 EXP-5），也不追求最细的阶段解释（那是 EXP-7）。它专注于一个问题：

> 当输出规模（用 $\alpha$ 控制 $|J|$）和样本量 $t$ 变化时，Adaptive 是否能在每个网格点选到**接近**“sampling 或 enum+sampling 中更快的那个”的分支？

因此 EXP-6 的设计是一个二维网格 sweep：

- 横轴：$t$（要抽多少样本）
- 纵轴：$\alpha$（用可控构造精确控制 $|J|$）

------

## 3. 数据集与工作负载：可控密度条带构造（SYN-CTRL）

### 3.1 为什么 EXP-6 使用 stripe（stripe_ctrl_alpha）

EXP-6 的关键是让“$\alpha$”真正稳定地对应 $|J|$，从而能把自适应拐点画得很清楚。
 因此 EXP-6 使用合成生成器 **stripe（stripe_ctrl_alpha 的别名）**，保证：
$$
k=\mathrm{round}(\alpha\cdot(n_r+n_s)),\qquad |J|=k.
$$
这点非常“paper-friendly”：你能用 $\alpha$ 精确扫过稀疏到稠密的 join 压力，而不受“随机波动导致 $|J|$ 乱跳”的影响。

> **注意（必须写在论文里）**：本文密度定义是

$$
\alpha = \frac{|J|}{n_r+n_s}
$$

不是 $|J|/(|R||S|)$。它代表“平均每个对象的匹配数量级”，更适合做 join 压力对齐。

### 3.2 EXP-6 默认数据规模与生成参数（来自脚本默认值）

- 维度：2D
- 规模：默认 $n_r=n_s=100000$（总对象数 $n_r+n_s=200000$）
- 合成种子：`EXP6_GEN_SEED=1`（固定数据，所有方法/变体共享同一份数据，保证公平）
- stripe 关键参数（默认）：
  - `control_axis=1`（控制轴，一般是 y 轴）
  - `core_lo=0.45, core_hi=0.55`（非控制轴上强制重叠的核心区间）
  - `gap_factor=0.10, delta_factor=0.25`（条带间隔与边界缓冲比例）
  - shuffle：`shuffle_strips=true, shuffle_r=false`（随机打散条带位置但保持结构）

------

## 4. 方法与三种运行模式（EXP-6 固定比较口径）

EXP-6 对每个 `method` 固定比较三种模式（variant）：

1. **sampling**
   - 不物化全部 $J$，直接返回 $t$ 个 i.i.d. 均匀样本
2. **enum_sampling**（枚举基线）
   - 先枚举全部 $J$，再从枚举结果中做均匀有放回抽样
3. **adaptive**（自适应）
   - 当 $J$ 小时走枚举（快且不随 $t$ 增长）
   - 当 $J$ 大时回退 sampling（稳且不会因 $|J|$ 爆炸）

### 4.1 自适应策略的阈值参数 $J_\star$

EXP-6 中 Adaptive 的关键控制参数是：

- `j_star = J_\star`（默认 `1,000,000`）

直观理解：

- 若 join 很小（$|J|\le J_\star$），adaptive 倾向于“枚举分支”
- 若 join 很大（$|J|>J_\star$），adaptive 会“pilot 探测后回退到采样”

在 stripe 构造下 $|J|\approx\alpha(n_r+n_s)$，因此自适应的理论拐点大约在：
$$
\alpha_\star \approx \frac{J_\star}{n_r+n_s}.
$$
例如默认 $n_r+n_s=200000$，$J_\star=10^6$，则 $\alpha_\star\approx 5$。

> **阅读提示**：如果你改变 $n_r,n_s$，要想保持相同的 $\alpha_\star$，就需要按比例调整 $J_\star$：
>  $J_\star \leftarrow \alpha_\star(n_r+n_s)$。

------

## 5. Sweep 设计：$(\alpha,t)$ 网格

EXP-6 默认扫：

- $\alpha$：`0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30`
- $t$：`1k, 3k, 10k, 30k, 100k, 300k, 1M`
- repeats：默认 3（建议论文最终用 5 或 7 提升稳定性）
- 线程：默认 1（公平性）

> **建议（用于论文最终图更“像相变”）**：在 $\alpha_\star$ 附近加密网格，例如加入 $\alpha=4,5,6,8$，更容易展示“自适应拐点”。

------

## 6. 计时口径（EXP-6 只看端到端，但保留阶段信息）

EXP-6 的主比较指标是 **端到端 wall time**（返回 $t$ 个样本的时间），并在 summary 中汇总成 median/p95。

**不计入**：

- 合成数据生成、磁盘读写（脚本“from zero”构建会发生，但算法比较以 sweep 输出的运行计时为准）

**计时覆盖**：

- Build / Count / Enumerate（如适用）/ Sample
   并记录分解（raw.csv）用于解释性能差异（但 EXP-6 的两张主图只用 wall 总时间）。

------

## 7. 输出文件与可复现性

运行脚本后，在 `EXP6_OUT_DIR`（默认 `results/sweeps/exp6_alpha_t/`）下生成：

- `manifest.txt`：完整环境与参数快照（时间戳、git sha、alphas、t、j_star、threads…）
- `exp6_alpha_t.json`：这次 sweep 的自包含配置
- `logs/*.log`：sjs_sweep 运行日志
- `sweep_raw.csv`：每个重复运行的详细记录（参数 + 阶段耗时 + 分支信息）
- `sweep_summary.csv`：按（alpha,t,method,variant）聚合后的 median/p95/ok_rate 等
- `figs/exp6_phase_<method>.png`：phase diagram（网格赢家图）
- `figs/exp6_ratio_<method>.png`：ratio heatmap

------

## 8. 两张核心图该怎么读

### 8.1 Phase diagram（赢家相图）

每个格点 $(\alpha,t)$ 标出三者谁最快：

- 0：sampling 赢
- 1：enum_sampling 赢
- 2：adaptive 赢
- NA：该点失败或缺失（应显式标注）

**你应该从这张图读到的“叙事”**：

- 低 $\alpha$：join 小 → 枚举相关（enum/adaptive）更可能占优
- 高 $\alpha$：join 大 → 采样更稳 → sampling 更可能占优
- 自适应应在两个区域“自动切换”，并在大部分区域赢或接近最优

### 8.2 Ratio heatmap（接近最优程度）

定义：
$$
\text{ratio}=\frac{T(\text{adaptive})}{\min(T(\text{sampling}),T(\text{enum\_sampling}))}.
$$

- ratio 越接近 1 越好
- ratio > 1：自适应比 oracle(min) 慢（存在“pilot 成本/决策不完美”）
- ratio < 1：自适应甚至比两条基线更快（可能来自工程路径复用/缓存/分支避免某些额外开销）

**你应该从这张图读到的“结论句模板”**：

- “Adaptive 在大多数 $(\alpha,t)$ 点位能贴近 oracle(min)（ratio≈1），并在拐点附近出现可解释的额外开销（ratio>1 的区域）。”

> **强烈建议**：在附录额外画一张 “adaptive_branch 热图”（枚举分支 vs 回退采样分支），能非常直观解释 ratio>1 的区域：这通常是 pilot 探测带来的固定成本。

------

## 9. 失败与截断的处理规则（必须明确，否则会被 reviewer 挑刺）

### 9.1 enum_cap 的语义风险

仓库支持 `enum_cap`（枚举安全上限）：

- 一旦开启（enum_cap>0），`enum_sampling` 可能会变成“枚举前缀 + 采样”，这会 **低估 enum 的真实时间**，从而污染 oracle(min)。

**EXP-6 主图的推荐规约：**

- **主实验必须使用 `enum_cap=0`**
- 若为了避免 OOM/极端慢必须开 cap：
  - 在图中必须显式标注该点为 “enum truncated/failed”，并从 oracle(min) 的分母中剔除该点（否则结论不可信）

### 9.2 ok_rate 与 correctness

EXP-6 默认不做重 correctness/uniformity 检验（那是 EXP-1 的任务）。
 EXP-6 的 `ok_rate` 更接近“该配置点是否正常跑完并输出结构合法样本”。

**建议在论文文本中明确**：

- “EXP-6 仅比较性能与自适应效果；正确性与采样质量已由 EXP-1 系统验证。”

------

## 10. 如何运行与复现实验（简明版）

在仓库根目录：

```
bash run/run_exp6.sh
```

常用覆盖参数（只举最常见的）：

```
# 提升统计稳定性
EXP6_REPEATS=5 bash run/run_exp6.sh

# 扫更多方法（更慢）
EXP6_METHODS="ours,kd_tree,r_tree" bash run/run_exp6.sh

# 调整自适应阈值（做敏感性）
EXP6_J_STAR=500000 bash run/run_exp6.sh
EXP6_J_STAR=2000000 bash run/run_exp6.sh

# 改规模（注意同时考虑 j_star 的缩放）
EXP6_N=200000 EXP6_J_STAR=2000000 bash run/run_exp6.sh
```

------

## 11. 论文写作建议：EXP-6 应该给出的“最小但有力结论”

主文建议输出两张图（phase + ratio）即可。配套文字建议包含：

1. 固定数据：stripe 控制 $|J|=\alpha(n_r+n_s)$，保证 $\alpha$ 可解释
2. 自适应策略：说明 $J_\star$ 与 $\alpha_\star\approx J_\star/(n_r+n_s)$
3. 结果：
   - phase 图展示“哪个区域哪种模式赢”
   - ratio 图展示“adaptive 多接近 oracle(min)”
4. 若 ratio>1 的区域明显：用 “adaptive_branch 图/或 raw 中 pilot_pairs” 解释为 pilot 固定成本，并在敏感性实验中展示调 $J_\star$ 后该区域如何收缩（这会显著提升说服力）

------

## 12. 一页 checklist（EXP-6 出图前自检）

-  `enum_cap=0`（或已在图中标注并剔除 truncated 点）
-  `threads=1` 且（最好）环境变量强制单线程（OMP/MKL 等）
-  repeats ≥ 5（论文最终图）
-  记录并附上 manifest + config + git sha
-  phase 图和 ratio 图都只用 ok_rate=1 的点（失败点要显式显示 NA）
-  若要解释拐点：额外导出 `adaptive_branch` 热图（推荐附录）