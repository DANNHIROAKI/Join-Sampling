# Run Test All（全仓库摔打性测试）说明文档

> 目标：对仓库进行一次**全面且严苛**的“摔打性测试（torture test）”，尽可能早地发现：
> - 语义/正确性错误（half-open 相交、join 语义、采样语义）
> - 算法实现偏差（oracle 不一致、枚举/自适应分支异常）
> - 工程性致命问题（崩溃、未定义行为、内存越界、潜在数据竞争）
> - 跑实验脚本的规范违背（build/results/temp/include 目录策略、内嵌 Python 等）

本说明文档对应脚本：

- `run/run_test_all.sh`：仓库级严苛测试入口（主脚本）
- `run/include/csv_check.py`：CSV 检查小工具（避免 bash 用 cut/awk 解析带引号 JSON 的 CSV 时误判）

---

## 0. 目录规范与“强制对齐点”

本测试脚本**同时也是一个“仓库规范守门人”**：它不仅跑程序，还会对 `run/*.sh` 做策略检查。

### 0.1 统一目录（与仓库约束对齐）

- **Build**：统一在根目录 `build/` 下（本脚本会用到）
  - `build/release/`
  - `build/debug/`
  - `build/asan/`（可选）
  - `build/tsan/`（可选）
- **结果输出**：统一在根目录 `results/raw/` 下（本脚本固定为）
  - `results/raw/test_all/`（包含 report、failures、logs、artifacts）
- **运行期临时文件**：统一在根目录 `run/temp/` 下（本脚本固定为）
  - `run/temp/test_all/`（生成的 configs、临时数据集等）
- **bash 内嵌 Python**：禁止（策略检查会抓），应剥离到
  - `run/include/*.py`（本脚本依赖 `run/include/csv_check.py`）

### 0.2 run 脚本策略检查（Step 0）

脚本会扫描 `run/*.sh`，默认**严格模式**（`POLICY_STRICT=1`）下，发现以下任一情况会直接 FAIL：

- 发现 **bash 内嵌 Python heredoc**（例如 `python3 - <<'PY' ...`）
- 发现 `build_*` 这类 **旧式 build 目录命名**（要求统一进 `build/`）
- 发现 `results/sweeps/...` 或 `results/expX/...` 这类 **旧式输出路径**
  - 必须统一用 `results/raw/...`

如果你在重构 run 脚本时违反了上述规范，`run_test_all.sh` 会第一时间把它拦下来。

---

## 1. 快速上手

### 1.1 放置文件

请确保仓库中存在以下两个文件（路径必须一致）：

```
run/run_test_all.sh
run/include/csv_check.py
```

并赋予执行权限：

```bash
chmod +x run/run_test_all.sh
```

### 1.2 一键启动（默认 deep 模式）

```bash
bash run/run_test_all.sh
```

deep 模式目标是“尽可能摔打”，覆盖构建矩阵、sanitizer、端到端 smoke、oracle gate、边界条件、随机 fuzz 等。

### 1.3 快速模式（更少 case、更少开销）

```bash
MODE=fast bash run/run_test_all.sh
```

fast 模式默认会减少 fuzz case，并默认关闭 ASAN（可按需手动打开）。

---

## 2. 输出与如何读报告

### 2.1 输出目录

测试完成后，核心产物在：

```
results/raw/test_all/
  ├─ REPORT.txt          # 总报告：每一步的命令、PASS/FAIL/SKIP
  ├─ FAILURES.txt        # 失败摘要（若有）
  ├─ env.txt             # 运行环境信息（uname/cmake/ctest/python/git 等）
  ├─ logs/               # 每个步骤一个 log，按 tag 命名
  └─ artifacts/          # 预留：可存放额外产物（目前主要用于扩展）
```

临时文件在：

```
run/temp/test_all/
  ├─ datasets/           # 生成的测试数据（bin/csv）
  └─ sweep_smoke.json    # sjs_sweep 的烟雾测试配置
```

build 产物在：

```
build/release/
build/debug/
build/asan/   (可选)
build/tsan/   (可选)
```

### 2.2 结果含义：PASS / FAIL / SKIP

- **PASS**：该项验证成功
- **FAIL**：该项验证失败（默认直接退出；若设置 `CONTINUE_ON_FAIL=1` 则继续跑完并汇总）
- **SKIP**：允许的跳过（常见原因：某方法/变体不支持；某工具未安装；某可执行不支持某 flag；或 sample 文件未生成导致 determinism 无法检验）

你应该先看：

1) `results/raw/test_all/REPORT.txt`  
2) 如有失败，再看 `results/raw/test_all/FAILURES.txt`  
3) 然后打开对应的 `results/raw/test_all/logs/<tag>.log`

---

## 3. 运行参数（环境变量）

> 下面所有参数都可以通过环境变量覆盖，适合 CI 或不同机器环境。

### 3.1 运行模式与失败策略

- `MODE=deep|fast`（默认 `deep`）
- `CLEAN=1|0`（默认 `1`，先删除旧的 `results/raw/test_all` 与 `run/temp/test_all`）
- `CONTINUE_ON_FAIL=0|1`（默认 `0`，一旦失败就立刻退出）
- `POLICY_STRICT=1|0`（默认 `1`，策略检查命中即 FAIL；设为 0 则只 warning）

### 3.2 资源控制

- `JOBS=<int>`：cmake build 并行度（默认自动检测）
- `TIMEOUT_SEC=<int>`：每个命令的超时（默认 600 秒；若系统无 `timeout` 命令则不启用）

### 3.3 Sanitizer / 静态检查

- `RUN_ASAN=1|0`：
  - deep 默认开启
  - fast 默认关闭
- `RUN_TSAN=1|0`：默认关闭（需要 clang/clang++ 且开销较大）
- `RUN_STATIC_ANALYSIS=1|0`：默认开启（缺工具会 SKIP，不强制安装）

### 3.4 数据与实验旋钮（用于 correctness 与 edge case）

- 数据规模：
  - `NR=2000 NS=2000 ALPHA=1 GEN_SEED=1`
- 采样规模：
  - `T=20000 SEED=1 REPEATS=2`
- oracle 上限：
  - `ORACLE_MAX_CHECKS=50000000`（用于 sjs_verify 的 oracle 检查预算）
- edge case：
  - `J_STAR_SMALL=1`（强制 adaptive 走 fallback 的 j_star）
  - `J_STAR_LARGE=1000000000`（强制 adaptive 走 enumerate 的 j_star）
  - `ENUM_CAP_TEST=1000`（enum_cap 截断测试用）

### 3.5 覆盖方法/变体集合

- `METHODS="ours aabb ..."`（空格分隔）
- `VARIANTS="sampling enum_sampling adaptive"`

---

## 4. 测试覆盖面总览（逐步说明）

脚本按“从规范→构建→端到端→语义正确性→边界→随机摔打→静态检查”的顺序执行。

### Step 0：run 脚本策略/卫生检查（Policy lint）

目标：防止 run 脚本逐渐“失控”，出现以下会直接破坏复现与对齐的情况：

- bash 内嵌 Python（heredoc）
- build 目录不在 `build/` 下
- results 不在 `results/raw/` 下

默认严格模式下任何命中都会 FAIL。

---

### Step 1：构建矩阵 + 单元测试（Release / Debug / ASAN / TSAN）

目标：最大化捕捉工程性问题。

- **Release 构建 + ctest**
- **Debug 构建 + ctest**
- **ASAN/UBSAN 构建 + ctest**（可选，deep 默认开启）
  - 用于抓：越界、UAF、UB、泄漏（取决于 ASAN_OPTIONS）
- **TSAN 构建 + ctest**（可选，默认关闭）
  - 用于抓：数据竞争（需要 clang）

> 通过标准：对应 build + ctest 必须全部通过（除非你主动关闭对应 sanitizer）。

---

### Step 2：定位核心可执行文件（Release）

目标：确保 build 产物可用，且路径不依赖固定 layout。

会在 `build/release` 下解析并定位：

- `sjs_gen_dataset`
- `sjs_run`
- `sjs_verify`
- `sjs_sweep`

缺任一都会 FAIL。

---

### Step 3：CLI smoke（--help 必须可用）

目标：最基础的“接口可调用性”检查（尤其适合 CI）。

对上述四个可执行运行 `--help`。

---

### Step 4：端到端 smoke（synthetic + binary）

目标：验证“生成→写盘→读取→跑算法→落盘 run.csv/samples”全链路可用。

包含：

1. 用 `sjs_gen_dataset` 生成离线数据（`.bin` + `.csv`）
2. `sjs_run` 跑 `ours` 的三种变体：
   - `sampling`
   - `enum_sampling`
   - `adaptive`
3. `sjs_run` 跑 binary 输入路径（保证 IO/解析一致）
4. baseline crash-guard：对 `r_tree` / `pbsm` 做一次 sampling 冒烟（允许 unsupported -> SKIP）

同时用 `csv_check.py` 检查 `run.csv` 是否非空（避免“跑了但没产出”的假阳性）。

---

### Step 4b：确定性检查（Determinism）

目标：同一数据、同一 seed、单线程下，样本输出应稳定一致。

- 连续跑两次 `ours/sampling`（`t=5000 seed=123 threads=1`）
- 对 samples 目录下的 `.tsv` 做 hash 比较
- hash 不一致会 FAIL

> 若 run 没写 samples 文件（或 samples 路径改变），该检查会 SKIP，并提示你查看输出格式。

---

### Step 5：正确性电池（Correctness battery，oracle gate）

目标：对每个 `method × variant` 做严格 correctness gate。

执行方式：调用 `sjs_verify`，并做以下硬性检查：

- **Gate 1（必须）**：`missing_in_universe == 0`
  - 任何非 0 直接 FAIL（这意味着输出样本不属于真实 join 结果集 J）
- **Gate 2（exact count）**：若日志行标记 `(exact)`
  - 必须 `count == oracle` 且 `rel_err == 0`
- 对 `(est)`：默认只 warning（可通过环境变量升级为 gate）
  - `CHECK_EST_COUNT=1` 且 `EST_REL_ERR_MAX=<阈值>`

此外会尝试跑 `sjs_verify` 的 binary-path（若 `--help` 支持 `--path_r`）。

---

### Step 6：边界条件测试（Edge-case）

目标：逼出“只在极端条件下出现”的逻辑错误。

包含：

1) **adaptive 分支强制切换**（如果 `sjs_run` 支持 `--j_star`）
- `j_star=J_STAR_SMALL`：期望走 fallback
- `j_star=J_STAR_LARGE`：期望走 enumerate  
并尝试从 `run.csv` 的 `adaptive_branch` 字段验证分支 token（字段缺失会 SKIP 并提示你更新检查条件）。

2) **enum_cap 截断显式性**（如果 `sjs_run` 支持 `--enum_cap`）
- 在较大 alpha 下跑 `ours/enum_sampling` 且 `enum_cap=ENUM_CAP_TEST`
- 期望出现 `enum_truncated==1` 或 `ok==0`（否则说明“截断未被显式记录”，会提示人工复核）

3) **空 join（alpha=0）处理**
- 期望该命令**失败**，或成功但在 `run.csv` 中明确标记 `ok=0`
- 若“成功但没有任何失败信号”，视为危险行为，直接 FAIL

---

### Step 7：sjs_sweep 冒烟（tiny config）

目标：验证 sweep 引擎能最小跑通并产出：

- `sweep_raw.csv`
- `sweep_summary.csv`

脚本会在 `run/temp/test_all/sweep_smoke.json` 写入一个极小 sweep（alpha×t×method×variant），并执行 `sjs_sweep --config=...`，随后用 `csv_check.py` 检查 CSV 非空。

---

### Step 8：随机 fuzz 电池（deep 模式）

目标：用大量随机小 case 覆盖“组合爆炸”的输入空间，尽可能触发罕见 bug。

- 随机抽取：
  - `n_r, n_s`（小规模）
  - `alpha`（从若干离散集合中选）
  - `generator`（stripe/uniform/clustered/hetero_sizes）
  - `t, seed`
- 调 `sjs_verify` 跑 `ours/sampling` 并检查 `missing_in_universe==0`  
  （unsupported -> SKIP）

fast 模式默认跳过 fuzz。

---

### Step 9：可选静态检查（best-effort）

目标：在不强制安装全部工具的前提下，尽可能覆盖代码质量问题。

- `shellcheck run/*.sh`（若安装）
- `python3 -m py_compile run/include/*.py`（若存在）
- `clang-tidy`：对 `compile_commands.json` 中前 20 个 TU 做检查（**warning-only，不作为 gate**）
- `cppcheck`：若安装则作为 gate（发现问题会 FAIL）

---

## 5. 常见失败与修复建议

### 5.1 Step 0 策略检查失败（最常见）

- 发现 run 脚本里有 `python3 - <<'PY'` 或 `<<'PY'`：
  - 把 Python 代码迁移到 `run/include/*.py`
  - bash 只负责调用，不再内嵌
- 发现 `build_*` 目录：
  - 统一改为 `build/<name>/...`
- 发现 `results/sweeps` 或 `results/expX`：
  - 统一改为 `results/raw/<...>`，并确保“新结果覆盖旧结果”的策略一致

### 5.2 ASAN 失败

说明存在真实的内存/未定义行为问题（越界、UAF、UB 等）。建议：

- 先打开失败对应的 `logs/ctest_asan.log` 或该 case 的 log
- 在本地用相同命令单独复现
- 最小化 case 再修（fuzz case 也会给出具体参数）

### 5.3 missing_in_universe != 0（正确性红线）

这是“样本不在真实 join 结果集”——属于语义级错误，应优先修复：

- half-open 相交判定是否统一？
- join 输出 pair id 是否正确对应原对象？
- adaptive / enum / sampling 分支在输出层是否复用了不同判定？
- 若某 baseline 语义不支持，应在程序侧显式标记 unsupported 并让脚本 SKIP，而不是输出错误样本。

### 5.4 Determinism hash 不一致

常见原因：

- 算法内部使用了非确定性迭代顺序（unordered 容器）
- 多线程引入竞态（脚本已强制 threads=1，但你可能仍有内部并行）
- seed 未真正控制所有随机源（例如只控制采样但未控制生成）

修复建议：把影响样本输出顺序/内容的随机源全部绑定到 seed，并保证单线程 deterministic。

### 5.5 空 join（alpha=0）没有明确失败信号

这是“坏味道”：程序看似成功但其实语义无法成立（无法从空集采样）。建议：

- 直接返回非 0 exit code；或
- 写出 `run.csv` 并显式标记 `ok=0`，并在 `note` 字段说明原因。

---

## 6. 扩展建议（怎么把它变成长期的 CI 资产）

你可以把 `run_test_all.sh` 当作“仓库质量门禁”，后续建议：

1. 在 CI 中至少跑一次 `MODE=fast`  
2. 合并前（或发 artifact 前）跑一次 `MODE=deep`  
3. 当新增方法/变体时：
   - 把 method token 加入 `METHODS_DEFAULT`
   - 确认 Step 5 的 correctness gate 能覆盖到
4. 当 `run.csv`/`sweep_*.csv` 格式变化时：
   - 同步更新 `csv_check.py` 的检查逻辑
   - 或在脚本里对缺失字段使用 SKIP 并给出明确提示（避免误判）

---

## 7. 附：一条推荐的“最小命令行套餐”

```bash
# 1) 快速门禁（适合频繁跑）
MODE=fast bash run/run_test_all.sh

# 2) 全量摔打（适合里程碑前）
MODE=deep RUN_TSAN=0 RUN_STATIC_ANALYSIS=1 bash run/run_test_all.sh

# 3) 发现问题后：保留现场（不要清理），方便调试复现
CLEAN=0 CONTINUE_ON_FAIL=1 bash run/run_test_all.sh
```

> 提醒：如果机器较慢或 CI 限制严苛，可以适当调大 `TIMEOUT_SEC`，或降低 `FUZZ_CASES`，再逐步加严。

