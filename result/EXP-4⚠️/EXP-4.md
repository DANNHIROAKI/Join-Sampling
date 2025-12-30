# EXP‑4 阅读文档（Scalability vs N）

---

## 0. EXP‑4 的定位与目标

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
3. **是否存在结构性 regime 导致方法天然占优？**
   - 输出敏感（output‑sensitive）方法在稀疏 regime 通常占优
   - 输出依赖（output‑dependent）枚举/分区类方法在高匹配压力 regime 可能退化或失败
   - 输出无关（output‑oblivious）采样/计数结构在高匹配压力 regime 更稳定

> 重要提醒（对齐总纲）：EXP‑4 是“规模扩展性实验”。  
> 输出爆炸/密度压力的主战场应由 **EXP‑3（$\alpha$ sweep）**承担；EXP‑4 只在固定 $\alpha$ 下观察随 $N$ 的增长趋势。

---

## 1. 任务定义

### 1.1 输入/输出语义

- 输入：二维 axis‑aligned rectangles 集合 $R,S$，均为 **half‑open** 语义
- 空间连接结果：
  $$
  J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq\emptyset\}
  $$
- 输出目标：从 $J$ 中抽取 **$t$ 个 i.i.d. 均匀、有放回样本**（每个样本是 pair 的标识）

### 1.2 half‑open 语义强调

- **边界贴合不算相交**（严格 `<`）
- 论文中必须明确：所有方法共享相同谓词（同一个 half‑open 相交判定）

---

## 2. 指标与计时口径

### 2.1 主指标：runtime–N

- 使用程序内部输出的 **`wall_ms`** 作为端到端算法时间
- 汇总：对每个点（N, method, variant）统计 **median**，误差条用 **p95**（建议主文统一 p95）

> 关键口径：**不自己“拼 phases 的和”**，主图用 `wall_ms`。

### 2.2 解释性指标：build time–N

- 从 `run.csv` 里的 `phases_json` 提取 `run_build_ms`（或对应 key）
- 同样做 median/p95 汇总

### 2.3 内存指标：RSS 峰值

- 使用外部工具抓取 **Maximum resident set size**
- 最新脚本要求 `time` 支持 `-v`（推荐 GNU time）
- RSS **只用于内存表**，不用于 runtime 曲线

### 2.4 明确“不计入”的部分

- 数据生成、数据加载、文件读写（尤其是 samples 写盘）**不计入**主要 runtime
- 论文中建议一句写死：  
  > “We start internal timing after dataset loading; thus wall_ms excludes I/O by construction.”

---

## 3. 变量与固定参数（EXP‑4 的实验设计）

### 3.1 变量：规模 $N$

- 扫描：$n_r=n_s=N$
- 主文推荐：
  $$
  N\in\{50k,100k,200k,400k,800k\}
  $$
- 若资源允许：附录扩到 1.6M（更能拉开趋势）

### 3.2 固定：样本量 $t$（脚本支持单个或列表）

两种常见设定（任选其一作为主文，另一种做补充）：

- `t = 1e5`：更突出 build/count 的扩展性
- `t = 1e6`：更放大 sample 阶段的增长项

**主文建议固定一个 t**。  
但最新脚本为方便 quick‑check / 生成补充图，支持：

- `T=<single>`（只跑一个 t）
- `T_LIST="100000 1000000"`（一次跑多个 t，输出目录自动按 `alpha*_t*` 分开）

### 3.3 固定：密度参数 $\alpha$（双 regime）

对齐总纲的密度定义：
$$
\alpha=\frac{|J|}{n_r+n_s}
$$

在 SYN‑CTRL 下 $n_r=n_s=N$ 时：
$$
|J|\approx 2\alpha N
$$
所以 α=5 是线性稀疏匹配压力；α=100/200 是更高匹配压力但仍为线性规模（不是 $|R||S|$ 级爆炸）。

**脚本支持两种“推荐默认”**（通过 `RUN_PROFILE` 控制）：

- `RUN_PROFILE=balanced`（论文默认口径）：`α = 5` 与 `α = 100`
- `RUN_PROFILE=ours`（更强的高压力 check，用于更清晰地观察交叉/斜率）：`α = 5` 与 `α = 200`  
  > 这符合总纲/EXP‑4 的合规做法：当 α=100 交叉点不明显，可提高到 150/200 后再正式扫 N。

你也可以显式覆盖：
- `REGIMES="alpha5:5 alpha100:100"`
- `REGIMES="alpha5:5 alpha200:200"`

### 3.4 repeats 与随机性控制

- 每点重复 `k=3` 或 `k=5`（主文若画 p95 更建议 5）
- 固定数据生成种子 `gen_seed=1`
- 固定 base 运行种子 `seed=1`，repeat 内按编号偏移（对齐实现）

最新脚本默认：
- `RUN_PROFILE=balanced`：`REPEATS=3`
- `RUN_PROFILE=ours`：`REPEATS=5`

### 3.5 线程设置（公平性）

- 主文统一单线程：
  - 环境：`OMP_NUM_THREADS=1` 等
  - 程序参数：**显式传 `--threads=1`**
- 最新脚本会 **强制** 将 `OMP/MKL/OPENBLAS/...` 的线程上限设置为 `THREADS`，并将 `--threads=$THREADS` 传入程序。

---

## 4. 数据与工作负载（EXP‑4 默认采用 SYN‑CTRL）

### 4.1 默认数据：SYN‑CTRL（可控密度条带构造）

使用 SYN‑CTRL 的原因（对齐总纲）：

- 能精确控制 $|J|$ → 精确控制 $\alpha$
- 不同 $N$ 下匹配压力可比（$|J|\propto N$）

### 4.2 推荐做法：离线二进制数据

为了严格满足“数据生成/加载不计入 runtime”的口径：

1. 使用 `sjs_gen_dataset` 离线生成 `*.bin`
2. EXP‑4 运行时使用 `--dataset_source=binary`

**脚本实际存放位置：**
- 离线数据在：`data/synthetic/exp4/<reg_dir>/`
- 为了集中查看，脚本会在：`run/temp/exp4/datasets/<reg_dir>/` 创建 symlink（可关闭 `LINK_DATA_IN_TEMP=0`）

---

## 5. 方法集合与运行模式（对齐总纲）

### 5.1 方法集合（默认覆盖多类别）

脚本默认会跑（可用 `METHODS=` 覆盖）：

- 本文方法：`ours`
- 索引采样类：`aabb interval_tree kd_tree r_tree range_tree`
- 分区/排序 join 类：`pbsm tlsop`
- 范围查询驱动采样：`sirs`
- 拒绝采样：`rejection`
- 映射/索引类：`tsunami`（若实现/可用）

### 5.2 三种运行模式（每个 baseline 都跑）

对齐总纲固定写法（脚本默认）：  
1) **sampling**  
2) **enum_sampling**（Enum+Sampling）  
3) **adaptive**

> 现实提醒（建议写进 paper）：  
> Adaptive 目标是接近
> $$
> \min\{T(\text{sampling}),T(\text{enum+sampling})\}
> $$
> 但通常略高于 min（pilot/决策开销）。不要预期 “Adaptive 永远全胜”。

---

## 6. 公平性机制与失败处理（EXP‑4 必须写清）

### 6.1 统一语义

- 所有方法严格遵守同一 half‑open 相交判定
- 目标分布完全一致：对 $J$ 的 i.i.d. 均匀、有放回抽样

### 6.2 枚举风险：enum_cap 与“截断点”的处理规则

枚举类（enum_sampling）可能在大 $|J|$ 下内存爆炸或时间爆炸，因此需要统一安全策略：

- 统一枚举上限：`enum_cap = C`
- 超过即视为 **枚举截断**（truncated）

**关键：截断点不能当成正常成功点画进主图**。建议规则：

- 对 **enum_sampling**：若 `run.csv` 任一 repeat 出现 `enum_truncated==1` → 该点视为失败/无效（图上用空心点/叉号/N/A 标注）
- 对 **adaptive**：pilot 的 `enum_truncated==1` 不等价于失败（可能 fallback 到 sampling，语义仍正确），需要结合 `phases_json` 的 branch 信息判定

> 最新脚本默认启用：`EXTRA_RUN_ARGS="--enum_cap=50000000"`（你仍可覆盖或清空它）。

### 6.3 unsupported 组合

- baseline×variant 若不支持：
  - 脚本继续跑（不中断 sweep）
  - 在 `exp4_status.csv` 里记录失败码与 stderr 路径
  - 图中标注 N/A / unsupported（不可静默忽略）

### 6.4 timeout（强烈建议）

为了避免某个点卡死拖垮 sweep：

- 脚本默认 `TIMEOUT_SEC=3600`（1h/点）
- **注意：timeout 为 best‑effort**，仅当系统存在 `timeout` 命令时生效（脚本会在开头打印 “timeout cmd: yes/no”）

---

## 7. 推荐运行方式（与脚本默认一致）

### 7.1 一键运行（当前脚本默认）

```
chmod +x run/run_exp4.sh
bash run/run_exp4.sh
```

**注意：默认 RUN_PROFILE=ours**，因此默认行为通常是：

- `N_LIST="50000 100000 200000 400000 800000"`
- `REGIMES="alpha5:5 alpha200:200"`
- `T_LIST="100000 1000000"`（一次跑两个 t）
- `REPEATS=5`
- `EXTRA_RUN_ARGS="--enum_cap=50000000"`
- `TIMEOUT_SEC=3600`

### 7.2 论文口径（与旧文档默认一致）

如果你要跑“论文默认设置”（α=5/100、t=1e5、repeats=3）：

```
RUN_PROFILE=balanced bash run/run_exp4.sh
```

### 7.3 常用覆盖参数

- 只跑一个 regime：

```
REGIMES="alpha5:5" bash run/run_exp4.sh
```

- 只跑一个 t（避免输出两套 t）：

```
T=100000 T_LIST="" bash run/run_exp4.sh
```

或只跑大样本：

```
T=1000000 T_LIST="" bash run/run_exp4.sh
```

- 只跑部分方法/变体：

```
METHODS="ours pbsm kd_tree" VARIANTS="sampling adaptive" bash run/run_exp4.sh
```

- 自己设 enum_cap / timeout：

```
EXTRA_RUN_ARGS="--enum_cap=80000000" TIMEOUT_SEC=5400 bash run/run_exp4.sh
```

- 复跑时跳过 build / gen：

```
SKIP_BUILD=1 SKIP_GEN=1 bash run/run_exp4.sh
```

---

## 8. 输出结构与“去哪找结果”（对齐最新脚本）

> **重要：脚本每次运行会覆盖 `results/raw/exp4/`**。  
> 若要保留历史，请在运行后手动复制：  
> `cp -a results/raw/exp4 results/raw/exp4_<TAG>`（TAG 可用 meta/manifest.txt 里的 RUN_ID）。

### 8.1 结果目录（覆盖写）

- 最终目录：`results/raw/exp4/`
- 临时目录：`run/temp/exp4/`（运行结束后已拷贝到 results）

### 8.2 每个 regime × t 的输出

- `results/raw/exp4/<reg_dir>/alpha<alpha>_t<t>/`

其中每个点 `(N, method, variant)`：

- `N<N>/<method>/<variant>/run.csv`（每个 repeat 一行，含 `wall_ms` 与 `phases_json`）
- `stdout.log / stderr.log`
- `mem/*.timev.log`（GNU time -v 输出，用于 RSS）

每个 `(regime, t)` 还会输出汇总 CSV：

- `exp4_status.csv`：exit_code、out_dir、stderr 路径、enum_truncated_any
- `exp4_rss_peak_kb.csv`：RSS 峰值（KB），enum_truncated_any

### 8.3 meta 信息（强烈建议保留）

- `results/raw/exp4/meta/manifest.txt`：RUN_ID、profile、参数快照
- `results/raw/exp4/meta/sysinfo.txt`：机器信息/编译器版本等
- `results/raw/exp4/commands.log`：全局命令日志
- `results/raw/exp4/meta/gen_reports/`：数据生成 JSON 报告（若有）

---

## 9. 从 run.csv 得到图表与阶段分解（最常踩坑）

### 9.1 runtime–N 主图

- y：`wall_ms`（median），误差条 p95
- x：N
- 建议按 regime 分图：α=5 一张、α=100/200 一张
- 如果脚本跑了多个 t（T_LIST），建议再按 t 分组（或挑一个 t 做主文）

> 不要用外部 time 的 elapsed 做 runtime，因为它包含 I/O。

### 9.2 build time–N

- y：`phases_json["run_build_ms"]`（median）
- x：N

### 9.3 RSS 峰值表

- 来源：`exp4_rss_peak_kb.csv`（转 MB）
- 主表建议至少包含 N=800k（最有说服力）

### 9.4 wall_ms 与 phases sum 不一致怎么办？

正常现象：析构/allocator 回收/日志开销会导致 `wall_ms` 略大。论文写法建议：

- “wall_ms is end‑to‑end algorithm time; phase timers cover major steps; remaining time is categorized as other/cleanup.”

---

## 10. 如何解读 EXP‑4（你希望读者得出的结论）

### 10.1 在 α=5（sparse）看到 PBSM 更快是合理现象

稀疏匹配压力下，输出敏感 join 的线性扫描常数小，领先并不意味着本文方法失败。可用叙事：

> “In low‑α regimes, output‑sensitive baselines remain competitive due to low join output and small constants.”

### 10.2 在高 match pressure（α=100/200）才是关键：看交叉点/斜率变化

主要观察：

1. PBSM runtime–N 斜率是否变陡（受 $|J|$ 增大影响）
2. 本文方法是否更稳定（更接近只随 $N$ 与 $\log N$ 变化）
3. enum_sampling 是否出现截断点（必须标注/排除）
4. Adaptive 是否能避开枚举风险并接近 min 分支

如果 α=100 下趋势不明显：提高到 α=150/200（先 quick check 再正式扫 N）。

### 10.3 Adaptive 的合理读法

不要期待 Adaptive 永远最快。应看：

- 是否接近 $\min\{T(\text{sampling}),T(\text{enum+sampling})\}$
- 是否能在大 N / 高 α 下安全避开 enum 风险（不崩溃、不截断）
- pilot 开销是否可接受（通过 phases 分解解释）

---

## 11. 跑前排雷清单（EXP‑4 必看）

1. **单线程写死**：`THREADS=1`（脚本会强制设线程环境变量）
2. **确保有 GNU time**：`time -v` 必须可用（macOS 建议安装 gnu-time 并使用 `gtime`）
3. **enum_sampling 截断处理**：画图必须读 `enum_truncated` / `enum_truncated_any`
4. **timeout**：脚本默认 1h/点，但仅在系统有 `timeout` 命令时生效
5. **results/raw/exp4 会被覆盖**：要保留历史请手动 copy

---

## 12. 论文写法模板（可直接拷贝）

### Setting（EXP‑4）

- “We evaluate scalability by sweeping $n_r=n_s=N$ from 50k to 800k.”
- “We fix $t=\dots$ and control join density via $\alpha=|J|/(n_r+n_s)$ using SYN‑CTRL.”
- “We report median wall‑clock time across $k$ repeats, with p95 error bars.”
- “We exclude dataset generation/loading from timing by construction, and separately report peak RSS.”
- “Enum+Sampling points that hit `enum_cap` are reported as truncated and excluded/marked.”

### Figures / Tables

- Fig X: runtime–N (α=5)
- Fig Y: runtime–N (α=100 或 α=200；请在 caption 写清)
- Fig Z: build time–N
- Table T: peak RSS (MB) at N=800k（可选再加 50k/200k）

---

## 13. 最小 sanity（正式扫 N 前强烈建议）

先跑一小点：

- `N=50k、α=5、t=10k`，只跑 ours/pbsm sampling

检查：

- `run.csv` 存在
- `count_exact==1`
- 对 SYN‑CTRL 预期：$|J|\approx 2\alpha N$，检查 `count_value` 是否一致
- 若启用 verify/quality：`missing-in-universe==0`

---

# 给组员的一段“EXP‑4 阅读总结”（可直接转发）

> EXP‑4 evaluates scalability by sweeping $N$ under fixed $t$ and fixed density $\alpha$.  
> We use SYN‑CTRL to control $\alpha=|J|/(n_r+n_s)$, and report median algorithm wall time `wall_ms` (excluding I/O by construction), build time from `phases_json`, and peak RSS from GNU `time -v`.  
> Enum+Sampling may hit `enum_cap`; truncated points must be explicitly marked/excluded.  
> The key interpretation is trend/phase bottlenecks and regime‑dependent dominance rather than “one method always wins”.
