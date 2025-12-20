# EXP-0 阅读文档：基础 Sanity（单元级 + 端到端 Smoke）

> 本文档用于“读懂并复现” EXP-0。它不是性能实验，也不用于支撑论文主结论；它的唯一目标是：在进入 EXP-1～EXP-N 之前，确认**环境、几何语义、join 流程、采样输出链路**都正确可用。

------

## 1. EXP-0 在整套实验中的位置

在你给出的实验总纲里，所有实验组遵循统一的流程范式：
 “生成数据 → 运行方法（3 模式）→ 记录 raw → 汇总 summary → 作图解释”。

EXP-0 是这个流程的“启动前自检”版本：

- 只做**最小但覆盖面足够**的 sanity：编译、单测、端到端 smoke、oracle 校验与基本质量指标输出。
- 通过后，才值得开始 EXP-1（Correctness & Sample Quality：多 seed、小规模真值、系统统计检验）以及后续 runtime / scalability / robustness 等实验。

------

## 2. 背景语义：EXP-0 到底在验证什么？

### 2.1 任务定义（必须一致，否则后面所有图都不可信）

给定两集合 $R,S$，元素是二维轴对齐矩形盒，采用 **half-open** 语义 $[lo, hi)$。空间连接结果为
$$
J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq\emptyset\}
$$
其中“相交”严格遵循 half-open（边界贴合不算相交）。

### 2.2 实验总目标（采样语义）

最终要从 $J$ 中抽取 $t$ 个 **i.i.d. 均匀、有放回样本**。
 因此：任何方法只要输出了不属于 $J$ 的 pair，都属于**正确性失败**（后续 EXP-1/RQ1 会更系统地检验）。

### 2.3 三种运行模式（后续论文统一口径）

总纲建议每个方法都用三种模式评估：

1. Sampling：不物化 $J$，直接输出 $t$ 个均匀样本
2. Enum+Sampling：先枚举 $J$，再从枚举结果均匀有放回抽样
3. Adaptive：在 $J$ 小时走枚举更快，$J$ 大时走采样更稳 

EXP-0 不一定把三种模式都全覆盖（它更像“冒烟测试”），但至少要保证：

- pipeline 能跑通
- 几何/相交判定一致
- join 流程与采样输出能与 oracle 对齐 

------

## 3. EXP-0 的一键入口

仓库提供的一键脚本是：

```
bash run/run_exp0.sh
```

脚本的设计意图就是“EXP‑0：基础 sanity（单元级 + 端到端 smoke）”，包含编译、ctest、生成数据、sjs_run smoke、sjs_verify 校验等步骤。

------

## 4. EXP-0 的完整流程拆解（对应脚本的 Step 1～5）

> 你阅读日志时，可以按下面的结构逐段对照定位问题。

### Step 1/5：Clean build + compile（构建链路 sanity）

**目的**：确认环境、依赖、编译选项、目标可执行文件都正常生成。

**产物（重要）**：

- `results/raw/exp0/logs/cmake_configure.log`
- `results/raw/exp0/logs/cmake_build.log`
- `results/raw/exp0/logs/env.txt`（记录 cmake/ctest/uname 等）

**必备可执行文件**（至少会用到）：

- `sjs_gen_dataset`：生成合成数据并导出 bin/csv
- `sjs_run`：跑一次实验并输出 run.csv / samples
- `sjs_verify`：小规模 oracle + 采样质量校验

------

### Step 2/5：ctest（单元级 sanity）

**目的**：快速覆盖几何判定、事件扫描、oracle join、采样质量检查、baseline smoke 等单元/集成测试，确保“底座逻辑”没问题。

**产物**：

- `results/raw/exp0/logs/ctest.log` 

**通过标准**：

- `ctest` 所有用例必须 100% 通过（否则不要进入 Step 3 以后）。

------

### Step 3/5：生成一个小数据集到磁盘（验证数据生成 + IO）

**目的**：

- 验证合成数据生成器可用
- 验证二进制格式写出可用（给 Step 4 的 binary-path smoke 用）

脚本使用 `--gen=stripe`（stripe 是 `stripe_ctrl_alpha` 的别名）。
 该类生成器在实验设计里属于“可控密度条带构造”，用于精确控制 $|J|$ / $\alpha$（paper friendly，后续 RQ2/RQ3/RQ6 常用）。实验大纲

**产物（Step 3 必须看到这些文件）**：

- `data/synthetic/exp0/exp0_R.bin`
- `data/synthetic/exp0/exp0_S.bin`
- （可选）`exp0_R.csv`, `exp0_S.csv`
- `exp0_gen_report.json`（记录生成器报告）

------

### Step 4/5：sjs_run 端到端 smoke（synthetic + binary）

**目的**：验证“从数据源到输出文件”的端到端链路：

1. **Synthetic on-the-fly**：在线生成数据并跑一次
2. **Binary input path**：从 `.bin` 加载数据并跑一次（专门验证 IO + 加载一致性）

脚本对一个代表性方法（通常是 `ours`）跑三种 variant：`sampling / enum_sampling / adaptive`。

**你应该检查的关键输出**：
 每个 run 的输出目录中，至少要有：

- `<out_dir>/run.csv`（每次 repeat 一行）
- `<out_dir>/samples/*.tsv`（如果启用 `--write_samples=1`）

> 这一步的意义不是比快，而是确认：
>
> - 不同 dataset_source 路径都通；
> - 不同 variant 的关键分支都通；
> - 能正确写出结果文件；
> - 不 crash、不异常退出。

------

### Step 5/5：sjs_verify（oracle + sample-quality 的“最低限度校验”）

**目的**：把“正确性”这一条在 EXP‑0 就先硬钉住：

- 先用 oracle（小规模可算真值）得到 $|J|$
- 收集 universe pairs（全体真值对）
- 检查采样输出是否都落在 universe 内
- 并输出 χ²/KS/自相关等质量指标（主要用于 sanity 观测；严格的系统检验在 EXP‑1 做）

总纲强调：

- **missing-in-universe 必须为 0**，否则属于正确性失败。
- 统计检验（χ²/KS p-value、自相关）用于验证均匀性/独立性 sanity。

脚本在 EXP‑0 中会对一组方法（ours/aabb/interval_tree/kd_tree/r_tree/range_tree/pbsm/tlsop/sirs/rejection/tsunami）执行 `sjs_verify` 的 sampling 模式验证。

------

## 5. EXP-0 的“通过判据”建议（读日志时按这个判断）

### 必须通过（硬门槛）

1. 编译成功，关键可执行文件存在（`sjs_gen_dataset/sjs_run/sjs_verify`）
2. `ctest` 全部通过（100%）
3. `sjs_run` 每个 smoke 的输出目录都有 `run.csv`（以及可选 samples）
4. `sjs_verify` 中 `missing_in_universe = 0`（任何方法出现非零都应判为 FAIL）

### 允许观察但不建议在 EXP-0 “卡死”的项

- χ²/KS p-value、自相关：EXP‑0 可以只做“是否能算出来、指标量级看起来不离谱”的 sanity；真正“系统统计意义上的不拒绝”应放在 EXP‑1 多 seed/多配置点里做。

> 经验上：如果你像目前一样在 EXP‑0 就能对多方法跑出 oracle+质量指标，且 missing=0，那么后续 EXP‑1/2/3 才有意义。

------

## 6. 产物清单（跑完 EXP‑0 你应该看到什么）

以脚本默认目录为例：

### 日志

- `results/raw/exp0/logs/env.txt`
- `results/raw/exp0/logs/cmake_configure.log`
- `results/raw/exp0/logs/cmake_build.log`
- `results/raw/exp0/logs/ctest.log`
- `results/raw/exp0/logs/gen_dataset.log`
- `results/raw/exp0/logs/run_*.log`
- `results/raw/exp0/logs/verify_*.log`

### 数据（落盘）

- `data/synthetic/exp0/exp0_R.bin`
- `data/synthetic/exp0/exp0_S.bin`
- `data/synthetic/exp0/exp0_R.csv`（可选）
- `data/synthetic/exp0/exp0_S.csv`（可选）
- `data/synthetic/exp0/exp0_gen_report.json`

### 运行结果（raw）

每个 run 的 `<out_dir>` 里：

- `run.csv`
- `samples/*.tsv`（如果写样本）

------

## 7. 常见失败形态与“读日志定位点”

> 这部分用于你带新人/合作者快速排障：看到哪类错误该先看哪里。

1. **编译失败 / 链接失败**
   - 先看 `cmake_configure.log` / `cmake_build.log`，确认依赖与编译器版本
   - 然后确认目标是否启用 `SJS_BUILD_ROOT_APPS / SJS_BUILD_TESTS`（脚本默认启用）
2. **ctest 不通过（几何/事件流/Join oracle）**
   - 先不要跑后续 Step 3～5
   - 优先修复单测：因为后续实验的所有结论都会建立在这些语义之上（尤其 half-open 判定）
3. **sjs_run crash 或 run.csv 缺失**
   - 先看对应 `run_*.log`，确认是否参数/路径错误
   - 再看输出目录权限/磁盘空间
   - 若是内存错误，建议用 sanitizer 构建复现（这属于工程保障，不是 EXP‑0 主线，但非常有效）
4. **sjs_verify 出现 missing_in_universe != 0**
   - 直接判 FAIL（说明有样本不在 $J$ 内）
   - 优先怀疑：几何相交判定不一致、join 流程漏判、样本对构造/回填错误等
   - 总纲明确这一条是正确性硬门槛。