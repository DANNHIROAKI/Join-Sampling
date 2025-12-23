# EXP‑4 阅读文档（Scalability vs N）

## 0. EXP‑4 的定位与目标（对应总纲）

### 0.1 对应研究问题（RQ4）

**RQ4：Scalability vs $N$**
 在固定采样任务语义下，随规模 $n_r=n_s=N$ 增大，比较不同方法的：

- **端到端 wall time–N**（主图）
- **Build time–N**（解释性强）
- **RSS 峰值**（内存表）
- **阶段瓶颈解释**（Build / Count / Enumerate / Sample 的贡献）

### 0.2 EXP‑4 的“最重要叙事”

EXP‑4 的核心不是“谁在所有设置永远最快”，而是回答：

1. **随 $N$ 增长，各方法是否保持可扩展（时间/内存趋势）？**
2. **瓶颈阶段在哪里（Build/Count/Sample/Enumerate）？**
3. **是否存在“结构性 regime”导致方法天然占优？**
   - 输出敏感（output‑sensitive）方法在稀疏 regime 通常占优
   - 输出依赖（output‑dependent）枚举/分区类方法在高匹配压力 regime 可能退化或失败
   - 输出无关（output‑oblivious）采样/计数结构在高匹配压力 regime 更稳定

> 重要提醒（对齐总纲）：EXP‑4 是“规模扩展性实验”。
>  输出爆炸/密度压力的主战场应由 **EXP‑3（$\alpha$ sweep）**承担；EXP‑4 只在固定 $\alpha$ 下观察随 $N$ 的增长趋势。

------

## 1. 任务定义（必须与总纲一致）

### 1.1 输入/输出语义

- 输入：二维 axis-aligned rectangles 集合 $R,S$，均为 **half‑open** 语义

- 空间连接结果：
  $$
  J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq\emptyset\}
  $$

- 输出目标：从 $J$ 中抽取 **$t$ 个 i.i.d. 均匀、有放回样本**（每个样本是 pair 的标识）

### 1.2 half‑open 语义强调

- **边界贴合不算相交**（严格 `<`）
- 论文中必须明确：所有方法共享相同谓词（同一个 half‑open 相交判定）

------

## 2. 指标与计时口径（对齐总纲，避免审稿人挑刺）

### 2.1 主指标：runtime–N

- 使用程序内部输出的 **`wall_ms`** 作为端到端算法时间
- 汇总：对每个点（N, method, variant）统计 **median**，误差条用 **p95** 或 stdev

> 关键口径：**不自己“拼 phases 的和”**，主图用 `wall_ms`。

### 2.2 解释性指标：build time–N

- 从 `run.csv` 里的 `phases_json` 提取 `run_build_ms`（或对应 key）
- 同样做 median/p95 汇总

### 2.3 内存指标：RSS 峰值

- 使用外部工具（推荐 GNU `time -v`）抓取 **Maximum resident set size**
- RSS **只用于内存表**，不用于 runtime 曲线

### 2.4 明确“不计入”的部分

- 数据生成、数据加载、文件读写（尤其是 samples 写盘）**不计入**主要 runtime
- 论文中建议一句写死：
   “We start internal timing after dataset loading; thus wall_ms excludes I/O by construction.”

------

## 3. 变量与固定参数（EXP‑4 的实验设计）

### 3.1 变量：规模 $N$

- 扫描：$n_r=n_s=N$

- 主文推荐：
  $$
  N\in\{50k,100k,200k,400k,800k\}
  $$

- 若资源允许：附录扩到 1.6M（更能拉开趋势）

### 3.2 固定：样本量 $t$

两种常见设定（任选其一作为主文，另一种做补充）：

- `t = 1e5`：更突出 build/count 的扩展性
- `t = 1e6`：更放大 sample 阶段的增长项

**建议主文固定一个 t**，避免不同 N 之间口径变复杂。

### 3.3 固定：密度参数 $\alpha$（双 regime）

对齐总纲的密度定义：
$$
\alpha=\frac{|J|}{n_r+n_s}
$$
EXP‑4 推荐一次跑两套（双 regime）：

- **Regime‑S（sparse / low‑α）**：`α = 5`
- **Regime‑H（high match pressure / high‑α）**：`α = 100`
   （如交叉点不明显，可在 quick check 后提高到 150/200，这在总纲允许范围内）

在 SYN‑CTRL 下 $n_r=n_s=N$ 时：
$$
|J|\approx 2\alpha N
$$
所以 α=5 是线性稀疏匹配压力，α=100 是高匹配压力但仍为线性规模（不是 $|R||S|$ 级爆炸）。

### 3.4 repeats 与随机性控制

- 每点重复 `k=3` 或 `k=5`（主文若画 p95 更建议 5）
- 固定数据生成种子 `gen_seed=1`
- 固定 base 运行种子 `seed=1`，repeat 内按编号偏移（对齐实现）

### 3.5 线程设置（公平性）

- 主文统一单线程：
  - 环境：`OMP_NUM_THREADS=1` 等
  - 程序参数：**显式传 `--threads=1`**（即使默认是 1，也建议写死）

------

## 4. 数据与工作负载（EXP‑4 默认采用 SYN‑CTRL）

### 4.1 默认数据：SYN‑CTRL（可控密度条带构造）

使用 SYN‑CTRL 的原因（对齐总纲）：

- 能精确控制 $|J|$ → 精确控制 $\alpha$
- 不同 $N$ 下匹配压力可比（$|J|\propto N$）

### 4.2 推荐做法：离线二进制数据

为了严格满足“数据生成/加载不计入 runtime”的口径：

1. 使用 `sjs_gen_dataset` 离线生成 `*.bin`
2. EXP‑4 运行时使用 `--dataset_source=binary`

------

## 5. 方法集合与运行模式（对齐总纲）

### 5.1 方法集合（示例）

覆盖多类别基线以增强结论普适性：

- 本文方法：`Ours`
- 索引采样类：AABB / IntervalTree / KDTree / RTree / RangeTree
- 分区/排序 join 类：PBSM / TLSOP
- 范围查询驱动采样：SIRS
- 拒绝采样：Rejection
- 映射/索引类：Tsunami（若实现/可用）

### 5.2 三种运行模式（每个 baseline 都跑）

对齐总纲固定写法：

1. **Sampling**
2. **Enum+Sampling**
3. **Adaptive**

> 重要现实提醒（建议写进 paper）：
>  Adaptive 目标是接近
> $$
> \min\{T(\text{sampling}),T(\text{enum+sampling})\}
> $$
> 但通常略高于 min（pilot/决策开销）。不要预期 “Adaptive 永远全胜”。

------

## 6. 公平性机制与失败处理（EXP‑4 必须写清）

### 6.1 统一语义

- 所有方法严格遵守同一 half‑open 相交判定
- 目标分布完全一致：对 $J$ 的 i.i.d. 均匀、有放回抽样

### 6.2 枚举风险：enum_cap 与“截断点”的处理规则

枚举类（enum_sampling）可能在大 $|J|$ 下内存爆炸或时间爆炸，因此需要统一安全策略：

- 统一枚举上限：`enum_cap = C`
- 超过即视为 **枚举截断**（truncated）

**关键：截断点不能当成正常成功点画进主图**（否则样本语义可能被破坏，审稿人会抓住）。
 建议规则：

- 对 **enum_sampling**：若 `run.csv` 任一 repeat 出现 `enum_truncated==1` → **该点视为失败/无效**（图上用空心点/叉号/N/A 标注）
- 对 **adaptive**：pilot 的 `enum_truncated==1` 不等价于失败（可能 fallback 到 sampling，语义仍正确），需要结合 `phases_json` 的 branch 信息判定

> 建议你们把这一条写进论文“failure handling”段落里，非常加分。

### 6.3 unsupported 组合

- baseline×variant 若不支持：
  - 在 status.csv 记录失败码
  - 图中标注 N/A / unsupported
  - 不可静默忽略

### 6.4 timeout（强烈建议）

为了避免某个点“卡死拖垮整套 sweep”，建议统一设置超时策略：

- `TIMEOUT_SEC`（例如 1h/点）
- 超时点作为失败点记录并标注

------

## 7. 推荐运行方式（双 regime 一键跑）

### 7.1 一键运行（推荐）

```
chmod +x run/run_exp4.sh
bash run/run_exp4.sh
```

默认：

- `N_LIST="50000 100000 200000 400000 800000"`
- `T=100000`
- `REPEATS=3`
- `REGIMES="alpha5:5 alpha100:100"`

### 7.2 常用覆盖参数（建议写进复现说明）

- 只跑一个 regime：

```
REGIMES="alpha5:5" bash run/run_exp4.sh
```

- 改 t / repeats：

```
T=1000000 REPEATS=5 bash run/run_exp4.sh
```

- 只跑部分方法/变体：

```
METHODS="ours pbsm" VARIANTS="sampling adaptive" bash run/run_exp4.sh
```

- 设置 enum_cap（强烈建议）：

```
EXTRA_RUN_ARGS="--enum_cap=50000000 --threads=1" bash run/run_exp4.sh
```

- 复跑时跳过 build / gen：

```
SKIP_BUILD=1 SKIP_GEN=1 bash run/run_exp4.sh
```

------

## 8. 输出结构与“去哪找结果”

每个 regime 输出独立目录，避免覆盖：

- `results/raw/exp4/<TAG>/<reg_dir>/alpha<alpha>_t<t>/`

其中每个点 `(N, method, variant)`：

- `N<N>/<method>/<variant>/run.csv`（每个 repeat 一行，含 `wall_ms` 与 `phases_json`）
- `stdout.log / stderr.log`
- `mem/*.timev.log`（GNU time -v 输出，用于 RSS）

每个 regime 还会输出汇总 CSV：

- `exp4_status.csv`：exit_code、out_dir、stderr 路径等
- `exp4_rss_peak_kb.csv`：RSS 峰值（KB）

------

## 9. 从 run.csv 得到图表与阶段分解（最常踩坑的地方）

### 9.1 runtime–N 主图

- y：`wall_ms`（median），误差条 p95/stdev
- x：N
- 建议分两张图：
  - α=5
  - α=100

> 不要用外部 time 的 elapsed 做 runtime，因为它包含 I/O。

### 9.2 build time–N

- y：`phases_json["run_build_ms"]`（median）
- x：N

这张图通常能解释：

- 索引类 baseline 的 build 是否 dominated
- 本文方法是否近线性×log 增长

### 9.3 RSS 峰值表

- 来源：`exp4_rss_peak_kb.csv`（转 MB）
- 推荐主表至少包含 N=800k（最有说服力）
- 可再补 50k/200k/800k 三点更清晰

### 9.4 wall_ms 与 phases sum 不一致怎么办？

正常现象：析构/allocator 回收/日志开销会导致 `wall_ms` 略大。
 论文写法建议：

- “wall_ms is end‑to‑end algorithm time; phase timers cover major steps; remaining time is categorized as other/cleanup.”

------

## 10. 如何解读 EXP‑4（你希望读者得出的结论）

### 10.1 在 α=5（sparse）看到 PBSM 更快是“合理现象”

稀疏匹配压力下，输出敏感 join 的线性扫描常数小，领先并不意味着本文方法失败。
 正确叙事（可写进 paper）：

> “In low‑α regimes, output‑sensitive baselines remain competitive due to low join output and small constants.”

### 10.2 在 α=100（high match pressure）才是关键：看交叉点/斜率变化

主要观察：

1. PBSM runtime–N 斜率是否变陡（受 $|J|$ 增大影响）
2. 本文方法是否更稳定（更接近只随 $N$ 和 $\log N$ 变化）
3. 是否出现内存风险点（尤其 enum_sampling）
4. Adaptive 是否能避开枚举风险并接近 min 分支

若 α=100 下没有明显趋势/交叉点：

- 合规做法：提高第二个 regime 到 α=150/200（先 quick check 再决定）

### 10.3 Adaptive 的合理读法

不要期待 Adaptive 永远最快。应看：

- 是否接近 $\min\{T(\text{sampling}),T(\text{enum+sampling})\}$
- 是否能在大 N / 高 α 下安全避开 enum 风险（不崩溃、不截断）
- pilot 开销占比是否可接受（通过 phases 分解解释）

------

## 11. 常见坑与跑前排雷清单（EXP‑4 必看）

1. **单线程一定要写死**：环境变量 + `--threads=1`
2. **enum_sampling 必须有截断/失败标注策略**：
   - 不能只看 exit_code，需要读 `run.csv` 的 `enum_truncated`
3. **强烈建议加 timeout**：避免单点拖垮 sweep
4. **两套 regime 的数据与结果必须分目录**：避免覆盖（脚本已做到）
5. **Tsunami/realdata_stub 若是占位**：主图不要作为强结论来源，需在图注说明

------

## 12. 论文写法模板（可直接拷贝）

### Setting（EXP‑4）

- “We evaluate scalability by sweeping $n_r=n_s=N$ from 50k to 800k.”
- “We fix $t=\dots$ and control join density via $\alpha=|J|/(n_r+n_s)$ using SYN‑CTRL.”
- “We report median wall‑clock time across $k$ repeats, with p95 error bars.”
- “We exclude dataset generation/loading from timing by construction, and separately report peak RSS.”

### Figures / Tables

- Fig X: runtime–N (α=5)
- Fig Y: runtime–N (α=100)
- Fig Z: build time–N
- Table T: peak RSS (MB) at N=800k（可选再加 50k/200k）

### Failure handling（必须写）

- “Enum+Sampling points that hit `enum_cap` are reported as truncated and excluded/marked, as they do not represent full uniform sampling over $J$.”

------

## 13. 最小 sanity（强烈建议正式扫 N 前做一次）

先跑：

- N=50k、α=5、t=10k，只跑 ours/pbsm sampling

检查：

- `run.csv` 存在
- `count_exact==1`
- 对 SYN‑CTRL 预期：$|J|\approx 2\alpha N$，检查 `count_value` 是否一致
- 若你启用 verify/quality：`missing-in-universe==0`

------

# 一段可直接发给组员的“EXP‑4 阅读总结”

> EXP‑4 evaluates scalability by sweeping $N$ under fixed $t$ and fixed density $\alpha$.
>  We use SYN‑CTRL to control $\alpha=|J|/(n_r+n_s)$ and run two regimes: α=5 (low‑α) and α=100 (high match pressure), where $|J|\approx 2\alpha N$.
>  We report median algorithm wall time `wall_ms` (excluding I/O by construction), build time from `phases_json`, and peak RSS from GNU `time -v`.
>  Enum+Sampling may hit `enum_cap`; truncated points must be explicitly marked/excluded.
>  The key interpretation is trend/phase bottlenecks and regime‑dependent dominance rather than “one method always wins”.