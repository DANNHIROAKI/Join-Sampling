# EXP-2 阅读文档：Runtime vs t

> 本文档用于复现实验 **EXP-2（RQ2: Runtime vs t）**，并且**严格与仓库的一键 runner：`run/run_exp2.sh` 的实际行为一致**（包括：默认 t-profile、默认 repeats/warmup、默认 `j_star`、输出目录结构与覆盖策略、绘图入口与参数）。

---

## 1. 实验目的与要回答的问题（RQ2）

**EXP-2** 用来回答 RQ2：

> 在固定数据集（固定 $R,S$）与固定密度参数 $\alpha$ 下，随着采样数量 $t$ 从 1k 增加到 1M（以及更大扩展区间），各方法的端到端 wall-clock 时间如何变化？是否呈现“采样型随 $t$ 增长；枚举型对 $t$ 的敏感度取决于其额外步骤（如 rank-sampling 排序）”的典型分化？

EXP-2 的核心价值是：**把“生成更多样本”这一维的代价分离出来**，让读者看到方法随 $t$ 的扩展性，并为后续 RQ3/RQ4/RQ6 的趋势解释提供 baseline。

---

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

---

## 3. EXP-2 的实验设计（“固定什么、扫什么”）

### 3.1 固定项（runner 默认已帮你做的）

- **固定数据集**：同一组 $R,S$ 在 sweep 全程不变（由 sweep JSON 决定）。
- **固定密度参数 $\alpha$**：EXP-2 不扫 $\alpha$（由 sweep JSON 决定）。
- **固定线程数**：默认单线程（`--threads 1`）。
- **固定是否写样本**：默认关闭写样本（`--write_samples 0`），避免 I/O 噪声。

> 以上运行侧固定项都可在 `run/run_exp2.sh` 里通过参数覆盖。

### 3.2 扫描项：只 sweep $t$（profiles）

runner 支持三套 t-range（称为 **profiles**）：

#### (A) paper profile（默认 t-list）
$$
t\in \{10^3, 3\cdot10^3, 10^4, 3\cdot10^4, 10^5, 3\cdot10^5, 10^6\}
$$

对应默认参数：
- `--t_list_paper "1000,3000,10000,30000,100000,300000,1000000"`

#### (B) extended profile（默认 t-list）
用于更大 $t$，更容易暴露渐近项主导/潜在 crossover：
- `--t_list_ext "1000,1000000,3000000,10000000,30000000"`

#### (C) full profile（**runner 默认 profile**）
**full** 会在同一次 sweep 中跑一条“全范围曲线”，默认取：
- `union(paper, extended)` 的 **去重 + 升序**集合

因此默认 full 的 $t$ 点为：
$$
t\in \{10^3, 3\cdot10^3, 10^4, 3\cdot10^4, 10^5, 3\cdot10^5, 10^6, 3\cdot10^6, 10^7, 3\cdot10^7\}.
$$

对应参数：
- 若不显式指定 `--t_list_full`，runner 会自动合并 `t_list_paper` 与 `t_list_ext` 生成 full 列表；
- 也可以手动覆盖：`--t_list_full "<csv>"`。

#### profile 选择：`--t_profile`
- `full`：只跑 full（**默认**）
- `paper`：只跑 paper（最贴合主文图的范围）
- `extended`：只跑 extended（更像 appendix/supplement 图）
- `both`：分别跑 paper 与 extended（两个目录，两套图）

> 注意：`full` 与 `both` 的区别是：  
> `full` 会把所有点画在同一条曲线里（一个目录）；`both` 会输出两套独立 profile（两个目录）。

### 3.3 方法与三种运行模式（论文推荐固定写法）

每个方法都以三种模式评估（由 sweep 配置决定 method×variant 的组合）：

1. **Sampling**：不物化全部 $J$，直接输出 $t$ 个样本。  
2. **Enum+Sampling（枚举基线）**：先枚举 $J$，再从枚举流中做均匀有放回抽样。  
3. **Adaptive**：阈值切换 “枚举” 与 “采样”，目标接近
   $$
   \min\{T(\text{enumerate}),T(\text{sampling})\}.
   $$

runner 默认会把 `base.run.j_star` 设为 **10000**（可用 `--j_star` 覆盖）。

---

## 4. 计时口径与统计口径（读图前必须明确）

### 4.1 计时覆盖范围

EXP-2 的主指标是 **端到端 wall-clock 时间**，覆盖方法内部典型阶段（随实现不同会合并/拆分）：

- Build：索引/预处理  
- Count：计数或权重计算  
- Enumerate（如适用）：枚举 join 流  
- Sample：生成 $t$ 个样本（含二次扫描/定位）  

> 实际落盘字段以 `sweep_raw.csv` 与 `phases_json` 为准。

### 4.2 repeats / warmup

runner 采用 “**warmup + 有效 repeats**” 的策略减少 cold-start 偏差；每个 profile 的默认值如下：

- **paper profile 默认：**
  - 有效 repeats：`--repeats_paper 10`
  - warmup（不参与绘图/统计）：`--warmup_paper 2`
  - 总 runs（传给 sweep）：12

- **extended profile 默认：**
  - 有效 repeats：`--repeats_ext 5`
  - warmup：`--warmup_ext 1`
  - 总 runs：6

- **full profile 默认：**
  - 有效 repeats：`--repeats_full 5`
  - warmup：`--warmup_full 1`
  - 总 runs：6

绘图时 runner 会把 `warmup_reps` 传给 plotter，从而自动排除 warmup 轮；同时传 `--min_repeats=<effective_repeats>`，从而只使用 fully-successful 的点。

> 快速 sanity 可以把 repeats 调小，但论文主图建议保持 paper 默认（10 个有效 repeats + p95）。

---

## 5. 配置文件与运行入口（你应该看哪里）

### 5.1 Sweep 配置文件（输入定义）

默认使用：
- `config/sweeps/sweep_t.json`

runner 会自动复制一份并做 **patch**（线程数、repeats、write_samples、j_star、t_list 等），所以你通常不需要手动改 JSON。

### 5.2 一键运行脚本

对应文件：`run/run_exp2.sh`

最常用（**默认 full profile**）：

```bash
bash run/run_exp2.sh
```

只跑 paper t-range（更贴合论文主文 t 范围）：

```bash
bash run/run_exp2.sh --t_profile paper
```

只跑 extended（用于观察大 t 渐近趋势）：

```bash
bash run/run_exp2.sh --t_profile extended
```

分别跑 paper + extended 两套（两个目录）：

```bash
bash run/run_exp2.sh --t_profile both
```

覆盖线程数 / 写样本 / adaptive 阈值：

```bash
bash run/run_exp2.sh --threads 1 --write_samples 0 --j_star 10000
```

覆盖 t-list（逗号分隔整数）：

```bash
# full（单条曲线）
bash run/run_exp2.sh --t_profile full \
  --t_list_full "1000,3000,10000,30000,100000,300000,1000000,3000000,10000000,30000000"

# paper + extended（两套输出）
bash run/run_exp2.sh --t_profile both \
  --t_list_paper "1000,3000,10000,30000,100000,300000,1000000" \
  --t_list_ext   "1000,1000000,3000000,10000000,30000000"
```

构建控制：

```bash
# 重新干净构建
bash run/run_exp2.sh --clean

# 跳过构建（已编译好时）
bash run/run_exp2.sh --no-build

# Debug/RelWithDebInfo 等
bash run/run_exp2.sh --build_type RelWithDebInfo
```

绘图控制：

```bash
# 只跑数据，不画图
bash run/run_exp2.sh --no-plot

# 改 Δruntime 的基线点 t0（注意：t0 必须在 t-list 中，runner 会强制检查）
bash run/run_exp2.sh --plot_t0 1000
```

> 依赖：runner 会检查 `cmake/jq/awk/sort/find/tee` 等；绘图需要 `python3`（若找不到 python3 会自动跳过绘图）。

---

## 6. 输出目录结构（与 runner 完全一致）

### 6.1 输出位置与覆盖策略（重要）

runner 的输出路径是固定的：

- **临时目录（每次运行都会覆盖）：**  
  `run/temp/exp2/`

- **最终目录（运行成功后覆盖）：**  
  `results/raw/exp2/`

> ⚠️ 注意：`results/raw/exp2/` 会在成功时被 `rm -rf` 后整目录覆盖写入。  
> 若要保留多次运行结果，请在每次跑完后手动拷贝：
>
> ```bash
> cp -a results/raw/exp2 results/raw/exp2_$(date +%Y%m%d_%H%M%S)
> ```

### 6.2 每个 profile 的目录内容

最终目录 `results/raw/exp2/` 下会按 profile 分子目录：

- `full/`（默认）
- `paper/`（当 `--t_profile paper` 或 `both`）
- `extended/`（当 `--t_profile extended` 或 `both`）

每个 profile 目录下包含：

- `sweep_raw.csv`：每次 repeat 一行（含 ok、wall、phases_json 等）
- `sweep_summary.csv`：按 method×variant×t 聚合后的统计表（median/p95/stdev/ok_rate 等）
- `sweep_original.json`：runner 拷贝的原始 sweep 配置
- `sweep_used.json`：runner patch 后实际使用的配置（包含 `meta.patch`）
- `MANIFEST.txt`：本次运行关键参数（threads、t_list、repeats、plot 选项等）
- `logs/`：
  - `sjs_sweep.log`
  - `plot_exp2.log`（若启用绘图）

此外，在 `results/raw/exp2/logs/` 还有全局日志/环境信息：

- `env.txt`：包含 timestamp、硬件/编译器信息、git 版本（若在 git repo 内）等
- `cmake_configure.log`、`cmake_build.log`：构建日志

---

## 7. 绘图与读图（由 runner 自动调用 plotter）

### 7.1 绘图脚本位置（runner 的实际选择逻辑）

runner 会优先使用：

- `run/include/exp2_plot.py`

若不存在，则 fallback 到：

- `run/plot_exp2.py`

你一般不需要手动调用 plotter；runner 会在每个 profile 运行后自动调用（除非 `--no-plot`）。

### 7.2 绘图输入与关键参数（runner 已固定/可覆盖）

runner 调用 plotter 时的关键参数包括：

- `--out_dir <profile_dir>`（例如 `.../run/temp/exp2/full`）
- `--t0 <plot_t0>`：Δruntime 基线点，默认 1000（runner 会检查 t0 必须存在于 t-list）
- `--warmup_reps <warmup>`：自动从 raw 中排除 warmup
- `--min_repeats <effective_repeats>`：只使用完整成功的 repeats 点
- `--mode paper|full`：默认 `paper`
- paper 图方法选择：
  - `--topk <int>`（默认 6）
  - `--always_include <csv>`（默认 `ours`）
  - `--paper_methods <csv>`（若指定，则精确使用此列表 + always_include）
  - `--exclude_methods <csv>`（可排除某些方法，如 rejection）

> 说明：**raw/summary 永远保留全量**。plot 的 method subset 只影响“主图展示集合”，不会删原始结果。

---

## 8. 理论预期与“看到什么算正常”

### 8.1 Sampling 模式：通常 ~线性增长

Sampling 需要真实产出 $t$ 个样本，因此随 $t$ 增长是必然的：

- 预期：`Δruntime`、`sample-phase` 近似线性；
- 解释：Sample 阶段每个样本通常 $O(1)$ 摊还，整体 $O(t)$。

### 8.2 Enum+Sampling 模式：可能出现 $t\log t$ 项

Enum+Sampling 通常包含：

- 生成 $t$ 个 ranks 并排序；
- 扫描 join 流定位这些 ranks；

因此：

- $t$ 很大时可能看到 **$t\log t$** 的增长；
- $|J|$ 很大时枚举扫描成本主导，曲线可能对 $t$ “不那么敏感”。

### 8.3 Adaptive：是否“切换”取决于 $J_\star$ 与 $|J|$

- 若 $|J| < J_\star$：可能走枚举分支；
- 若 $|J| > J_\star$：会切到 sampling 分支。

> 读图时应同时报告 raw 中的分支统计（例如 pilot 规模 / 是否走枚举），否则读者会问“为什么 adaptive 不切换/为什么和某个模式一模一样”。

---

## 9. 公平性与失败处理（与 runner 的行为一致）

- runner **不会删除或隐藏**任何 raw/summary 结果。
- 绘图时可能会：
  - 排除 warmup repeats；
  - 仅使用满足 `--min_repeats` 的点；
  - 过滤失败点（ok_rate < 1）；
  - 对 paper mode 只展示 top-k 或指定 method 子集（但 raw 仍保留全量）。

若启用了枚举上限（enum_cap 等），应在论文说明截断点/失败点，不要静默消失。

---

## 10. 推荐写入论文的 EXP-2 小节结构（模板）

（可直接复用）

1. **Setup**：固定数据集（给出 $n_r,n_s,\alpha$）、单线程、paper profile t-range、有效 repeats=10（+2 warmup）、报告 median+p95。  
2. **Result (Sampling)**：runtime 随 $t$ 增长；Δruntime/sample-phase 近线性；给出 ns/sample 量级。  
3. **Result (Enum+Sampling)**：解释 $t\log t$（rank sort）与 $|J|$ 枚举项的主导区间。  
4. **Result (Adaptive)**：给出不同 t 区间/不同方法的分支选择现象；讨论其接近 $\min\{\cdot\}$ 的程度。  
5. **Takeaway**：一句话总结“谁的斜率更小/谁的固定开销更大/谁对 t 更稳”。

---

## 11. 常见坑与排障清单（跑完 EXP-2 必查）

1. 固定数据集下，`count_mean` 是否随 t 基本不变？（summary 中可快速 sanity）  
2. `ok_rate` 是否为 1？若不是，查看 `logs/sjs_sweep.log` 与 raw。  
3. 主图看起来不随 t 变：是否常数项过大？请看 Δruntime 与 sample-phase。  
4. repeats 太少却画 p95：建议增加 repeats（paper 默认已是 10 有效 repeats）。  
5. 你改了 `--plot_t0` 却报错：t0 必须出现在 t-list 中（runner 会强制检查）。  
