# EXP‑3 阅读文档：Runtime vs 密度参数 α（对应 RQ3）

本文档用于**阅读、复现、解释**实验 EXP‑3（RQ3：Sensitivity to density α）的结果。目标是让读者在不读代码的情况下，也能理解：实验做了什么、曲线应该怎么看、异常点怎么处理、结论如何写得“paper‑friendly”。

> **推荐复现入口（与最新 `run/run_exp3.sh` 对齐）**：  
> `bash run/run_exp3.sh`  
> 该脚本会：自动编译 → 生成 effective config → 跑 sweep → 导出 adaptive 分支比例 →（可选）画图 → 同步产物到 `results/raw/exp3`（或按时间戳保存）。

------

## 1. EXP‑3 在整篇论文中的位置

- **RQ1** 负责证明“输出语义正确：对 $J$ 的 i.i.d. 均匀有放回抽样”
- **RQ2** 负责证明“随样本量 $t$ 增大，时间如何增长”
- **RQ3（本实验）** 负责证明“随连接密度 α 增大，方法是否稳定/是否退化/自适应是否会换挡”

EXP‑3 的叙事重点不是“采样质量”（那是 RQ1），而是**性能曲线在稀疏→中等→稠密压力下的退化机制**与**自适应策略的拐点行为**。

------

## 2. 任务定义与密度参数（必须在论文写清）

### 2.1 任务定义（与全篇一致）

给定二维轴对齐矩形集合 $R,S$（**half‑open** 语义）：
$$
r=[L_x(r),R_x(r))\times[L_y(r),R_y(r)),\quad
s=[L_x(s),R_x(s))\times[L_y(s),R_y(s)).
$$
Join 结果：
$$
J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq\emptyset\},
$$
其中 **边界贴合不算相交**（half‑open：例如 $R_x(r)=L_x(s)$ 不相交）。

实验输出：从 $J$ 中抽取 $t$ 个 **i.i.d. 均匀、有放回**样本对 $(r,s)$。

### 2.2 密度参数 α（EXP‑3 的核心定义）

本文使用：
$$
\alpha = \frac{|J|}{n_r+n_s}
$$
而不是 $|J|/(|R||S|)$。直观含义：**每个对象平均匹配量级**，便于跨规模对齐与压力测试。

> 注意：合成数据生成器通常用  
> $k=\mathrm{round}(\alpha\cdot(n_r+n_s))$ 生成目标 $|J|=k$，  
> 因此实际 $\alpha$ 可能是离散近似（论文里建议写 “α_target/α_achieved” 或说明 rounding）。

------

## 3. RQ3 要回答什么（Claims / 读图的“问法”）

EXP‑3 的核心输出应回答以下问题（建议论文按这套问题写解释）：

1. **稳定性**：哪些方法对 α 不敏感（曲线近乎平坦）？
2. **退化分界**：哪些方法随 α 明显变慢？退化原因是什么（Build/Count/Sample 哪个阶段爆炸）？
3. **极端稠密**：α 很大时是否出现失败/截断（ok_rate<1、enum_truncated、OOM/超时）？
4. **自适应拐点**：Adaptive 是否在合适的 α 区间从枚举切换到采样？拐点是否与阈值 $J_\star$ 一致？
5. **拒绝采样类（如 Rejection）**：α 增大时 acceptance 是升还是降？曲线是否符合该机制？

------

## 4. 实验设置与计时口径（强烈建议固定写法）

### 4.1 计时口径（论文必须写清）

EXP‑3 比较聚焦“算法本体”——建议以**runner 内部 wall time**作为主口径：

- 数据生成 / 文件读写 **不计入**主要计时
- 计时覆盖：Build / Count / Enumerate（如适用）/ Sample
- 输出记录阶段分解（用于解释退化原因）

> 读结果时：主图建议用 **端到端 wall time**（median + p95），并可用 `phases_json` 做阶段解释图（EXP‑7 风格，但 EXP‑3 也推荐附一个小图/附录）。

### 4.2 随机性与重复

- 每个点重复 $k$ 次（推荐 3 或 5）
- 每次重复 seed 不同（base seed + rep_id 偏移）
- 主结果报告 median，误差条报告 p95（或 stdev）

### 4.3 线程与公平性

主文默认 **单线程**（`sys.threads=1`），确保公平比较；多线程作为补充实验放附录。

> 与脚本对齐：`run/run_exp3.sh` 默认 `EXP3_THREADS=1`；若你设置 `EXP3_THREADS>1`，脚本也会同步设置常见的线程池上限（如 `OMP_NUM_THREADS` 等）到同一值，避免“明明指定多线程但被外部线程池强制 1”的歧义。

------

## 5. 数据与工作负载（EXP‑3 专用）

### 5.1 数据生成器：SYN‑CTRL（可控密度条带）

EXP‑3 使用 **可控密度条带构造**（stripe 控制 $|J|$ 从而控制 α）：

- 固定 $n_r,n_s$
- 扫描 α
- 生成保证 $|J| \approx \alpha(n_r+n_s)$（取整误差）

> 为什么它对 EXP‑3 “paper‑friendly”：  
> 因为 $|J|$ 是你可控的，所以你能把退化归因于“密度压力”，而不是数据集偶然性。

### 5.2 固定参数（典型默认，建议先跑通）

- $n_r=n_s=100000$
- $t=100000$
- repeats=3
- 单线程

> 注意：上述是“典型默认”。实际跑出来的默认值以  
> `config/sweeps/sweep_alpha.json` 与脚本生成的 `run/temp/exp3/sweep_exp3_effective.json` 为准。

------

## 6. 方法集合与运行模式（EXP‑3 的比较口径）

### 6.1 三种模式（统一语义）

- **Sampling**：不物化 $J$，直接输出 $t$ 个样本
- **Enum+Sampling**：先枚举 $J$ 再抽样（通常在大 α 会失败/截断）
- **Adaptive**：$|J|$ 小走枚举，$|J|$ 大回退到采样（通过阈值或 pilot 决策）

> 论文写作建议：EXP‑3 主图一般只画 **Sampling + Adaptive**；  
> Enum+Sampling 可在小 α 区间出现，超出后用失败标注/截断，不要强行跑满 0–300。

### 6.2 特殊说明：Rejection 方法的变体选择

仓库实现中通常会出现“Rejection + Sampling 被禁用/跳过”的情况（为了避免语义/工程不合理组合）。因此：

- **Rejection 建议只展示 Adaptive**（并在图注注明 sampling 被跳过）
- 画图时要过滤 `ok_rate=0` 或 `wall_*=-1` 的记录，避免 log 轴崩溃

------

## 7. 如何运行 EXP‑3（与 `run/run_exp3.sh` 对齐的推荐路径）

### 7.1 一键运行（推荐）

在仓库根目录执行：

```bash
bash run/run_exp3.sh
```

脚本会在 `run/temp/exp3` 产生中间产物，然后同步到 `results/raw/exp3`（或按时间戳保存）。

#### 常用环境变量（脚本级）

- **构建与测试**
  - `EXP3_BUILD_TYPE=Release|Debug|RelWithDebInfo|MinSizeRel`
  - `EXP3_CLEAN_BUILD=1`：清理 build 目录后重建
  - `EXP3_JOBS=16`：编译并行度
  - `EXP3_RUN_TESTS=1`：build 后跑 `ctest`

- **实验控制（覆盖 sweep JSON 的关键参数）**
  - `EXP3_NR`, `EXP3_NS`
  - `EXP3_T`：采样次数 t
  - `EXP3_REPEATS`
  - `EXP3_JSTAR`：自适应阈值（决定理论换挡点）
  - `EXP3_ENUM_CAP`：枚举上限（防止 OOM）
  - `EXP3_THREADS`：`sys.threads`（默认 1）
  - `EXP3_METHODS` / `EXP3_VARIANTS`：只跑部分方法/模式（开发时可用）

- **α 列表（让“换挡点证据”更硬）**
  - `EXP3_ALPHA_LIST="0,0.03,0.1,0.3,1,3,4,5,6,10,30,100,300"`：显式指定点位
  - `EXP3_AUTO_ALPHA_LIST=1`（默认）：若未指定 `EXP3_ALPHA_LIST`，脚本会采用更 paper‑friendly 的默认列表；若能从 `EXP3_NR/EXP3_NS/EXP3_JSTAR` 推出 $\alpha_\star$，则会在 $\alpha_\star$ 附近自动加密整数点位。

- **画图与输出**
  - `EXP3_PLOT=0`：跳过画图
  - `EXP3_PLOT_ENUM=1`：额外生成包含 enum_sampling 的图（用于附录）
  - `EXP3_KEEP_HISTORY=1`：结果写入 `results/raw/exp3/<timestamp>/`，并更新 `results/raw/exp3/latest`（软链，best effort）

#### 推荐命令示例

1) **主文单线程、默认点位**（推荐）：

```bash
bash run/run_exp3.sh
```

2) **明确加密换挡点附近 α（更像论文）**：

```bash
EXP3_ALPHA_LIST="0,0.03,0.1,0.3,1,3,4,5,6,10,30,100,300" \
bash run/run_exp3.sh
```

3) **保留历史结果 + 生成附录 enum 图**：

```bash
EXP3_KEEP_HISTORY=1 EXP3_PLOT_ENUM=1 bash run/run_exp3.sh
```

4) **多线程附录实验（注意：主文仍建议单线程）**：

```bash
EXP3_KEEP_HISTORY=1 EXP3_THREADS=8 bash run/run_exp3.sh
```

### 7.2 手动运行（开发调试）

如需跳过脚本、直接跑 sweep（例如调参/单点 debug）：

```bash
./sjs_sweep --config=../config/sweeps/sweep_alpha.json
```

> 若不想覆盖旧结果：复制一份 JSON，修改 `base.output.out_dir`。

### 7.3（可选）先跑一次 sanity（质量与几何判定）

```bash
./sjs_verify \
  --dataset_source=synthetic --gen=stripe --dataset=verify \
  --n_r=2000 --n_s=2000 --alpha=0.3 --gen_seed=1 \
  --method=ours --variant=sampling --t=20000 --seed=1 --repeats=1
```

------

## 8. 输出文件与字段解读（读结果必看）

### 8.1 输出目录结构（脚本产物）

脚本会先写到：

- `run/temp/exp3/`

然后同步到：

- 默认覆盖模式：`results/raw/exp3/`
- 若 `EXP3_KEEP_HISTORY=1`：`results/raw/exp3/<timestamp>/`，并尝试更新 `results/raw/exp3/latest`

典型文件：

- `sweep_raw.csv`：每次 repeat 一行（最细粒度）
- `sweep_summary.csv`：按 (alpha, method, variant, …) 聚合（画主图首选）
- `derived/adaptive_branch_ratio.csv`：从 raw 导出的自适应分支比例
- `meta/manifest.txt`：运行环境、git 信息、参数覆盖等（artifact 友好）
- `sweep_exp3_effective.json`：本次实验真正使用的 effective config（可复现）
- `plots/`：图像输出目录（见第 9 节）

### 8.2 `sweep_summary.csv` 常用字段（主图用）

- `alpha`
- `method`, `variant`
- `wall_median_ms`, `wall_p95_ms`
- `ok_rate`：成功完成的比例（非常关键）

> **读图规则 1（强制）**：`ok_rate==0` 或 `wall_median_ms<=0` 的点不要上 log‑y 图；应当标注为失败/跳过。

### 8.3 `sweep_raw.csv` 常用字段（解释与分支比例）

- `ok`：单次运行是否成功
- `adaptive_branch`：`enumerate_all` / `fallback_sampling` / `fallback_sampling_no_pilot`
- `enum_truncated` / `enum_cap`：枚举是否触顶
- `phases_json`：阶段分解（Build/Count/Sample 等）

------

## 9. 作图规范（与脚本输出目录对齐）

### 9.1 主图：runtime vs α（用 summary）

- x：α（建议 `symlog`，因为 α=0 不能 log）
- y：`wall_median_ms`（通常用 log‑y；并可加 p95 误差条）
- 线：method×variant（常用做法：每个 method 一张图，里面画 sampling/adaptive 两条）

**必须过滤：**

- `ok_rate==0`
- `wall_median_ms<=0`（如跳过行会是 -1）

### 9.2 副图：自适应分支比例 vs α（用 raw）

过滤 `variant==adaptive` 且 `ok==1`：

- `enumerate_frac = count(adaptive_branch==enumerate_all)/runs`
- `fallback_frac = 1 - enumerate_frac`（可把 no_pilot 合并进 fallback）

这个图用于回答：“adaptive 是否在合理 α 区间换挡”。

### 9.3（强烈建议）解释图：phase breakdown vs α（用 phases_json）

从 raw 的 `phases_json` 提取并聚合（median/p95），画：

- `Build_ms`、`Count_ms`、`Sample_ms` 随 α 的变化（堆叠柱/多折线均可）

用于解释：某方法为何在大 α 爆炸（通常是 Count 或 Sample 阶段变慢）。

### 9.4 与脚本生成图的对应关系

- 默认主图输出：`plots/main/`（不含 enum_sampling；主文更干净）
- 若 `EXP3_PLOT_ENUM=1`：额外输出 `plots/with_enum/`（含 enum_sampling；适合附录/小 α 子图）

------

## 10. 如何解释曲线（论文叙事模板）

建议把 α 轴分成三段来写（读者更容易理解）：

### 10.1 稀疏区（α 很小）

- $|J|$ 很小，枚举成本低
- Adaptive 应更倾向于 `enumerate_all`
- Sampling 类方法通常会有较高固定开销（build/count），但总体仍可接受

### 10.2 中等区（α 接近换挡点）

理论换挡点：
$$
\alpha_\star \approx \frac{J_\star}{n_r+n_s}
$$
因此如果 `j_star=10^6, n_r=n_s=10^5`，则 $\alpha_\star\approx 5$。

**写作建议：**

- sweep 点位要覆盖 $\alpha_\star$ 附近（建议加 α=4/5/6；脚本默认会提供 paper‑friendly 列表）
- 用 branch ratio 图证明换挡发生在接近 $\alpha_\star$ 的位置

### 10.3 稠密区（α 很大）

- 枚举类方法可能出现：超时 / 内存爆炸 / enum_truncated
- 索引/采样类方法应更稳定（曲线更平）

**对拒绝采样类（Rejection）的正确解释：**

- 不要先入为主写“α 大 acceptance 下降”。在很多构造下，**join 越密 acceptance 反而越高**，于是 runtime 可能下降。
- 正确做法：用 `phases_json` 或接受率相关统计（如果有）来佐证它的机制。

------

## 11. 失败点与公平性（图里必须显式处理）

- 枚举失败/截断：
  - 由 `enum_truncated`/`enum_cap` 标识
  - 图上用 “× / 空心点 / 标注 fail” 显示，不要静默删掉
- 被实现禁用的组合（如某些 `method+variant` 被 SKIPPED）：
  - 仍可在表中说明（ok_rate=0），图里跳过并在 caption 写原因

------

## 12. EXP‑3 的复现检查清单（跑完后必做）

1. `sweep_summary.csv` 中每个点是否 `ok_rate≈1`（除明确 SKIPPED/失败）
2. `sweep_raw.csv` 中 adaptive 的 `adaptive_branch` 是否出现从 enumerate 到 fallback 的切换
3. 画主图前是否过滤了 `wall<=0` 的行（避免 log 轴崩）
4. 至少对 2–3 个代表点（稀疏/中等/稠密）抽 `phases_json` 做解释（哪一阶段导致差异）
5. 图注中是否写明：α 定义、half‑open 语义、计时不含数据生成/读写、单线程

------

## 13. 推荐的“最小增强”（让 EXP‑3 更像论文）

- 把 α 列表在换挡点附近加密：例如  
  `0,0.03,0.1,0.3,1,3,4,5,6,10,30,100,300`  
  （脚本默认会提供 paper‑friendly 列表；你也可以用 `EXP3_ALPHA_LIST` 显式指定。）
- 主图只放 sampling/adaptive；enum+sampling 放一个“小 α 子图或附录”
- 增加一张小图：phase breakdown（Build/Count/Sample）随 α 的变化（解释力极强）

------

### 一句话总结（给论文的 caption 可用）

> **EXP‑3（RQ3）**：在固定规模与固定样本量下扫描密度参数  
> $\alpha=|J|/(n_r+n_s)$，比较不同方法在稀疏→稠密 join 压力下的端到端采样时间，并用自适应分支比例与阶段分解解释退化分界与换挡行为（half‑open 相交语义，计时不含数据生成/读写，单线程，结果报告 median+p95）。
