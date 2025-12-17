# 00_overview.md

本仓库是一个 **SIGMOD 级别实验用**的 C++ 框架，用来系统评估：

> 给定两组轴对齐半开盒（目前 Dim=2，为矩形）集合 \(R, S\)，其空间 Join 结果集合  
> \[
J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq\emptyset\}
\]
如何从 \(J\) 中抽取 **t 个 i.i.d. 均匀（有放回）样本**，并在不同数据分布/规模下比较各方法的时间、内存、误差与可扩展性。

---

## 1. 关键约定（决定“正确性”）

### 1.1 半开区间（half-open）

几何统一采用 **半开盒**：

- 2D: \([x_\text{lo}, x_\text{hi}) \times [y_\text{lo}, y_\text{hi})\)
- dD: \(\prod_i [\text{lo}_i, \text{hi}_i)\)

相交判定（逐维）：
\[
\forall i,\ \max(\text{lo}_i(a),\text{lo}_i(b)) < \min(\text{hi}_i(a),\text{hi}_i(b)).
\]

半开语义的意义：

- 避免“贴边是否相交”的二义性；
- 扫描线事件 tie-break 可以严格与 half-open 对齐（同坐标先 END 后 START）。

### 1.2 扫描线事件排序

当前框架多种方法都依赖 **Plane Sweep**（在 Dim=2 上是 x 轴 sweep）：

- 每个 box 产生 `START`（在 lo）与 `END`（在 hi）事件。
- 排序规则（稳定）：
  1) 按 sweep 坐标升序  
  2) 同坐标：`END` 在 `START` 前（匹配 half-open）  
  3) 同坐标且都是 `START`：固定 side tie-break（默认 R 在 S 前）  
  4) 同坐标且都是 `END`：无关紧要

这样可保证：任意 join pair \((r,s)\) 被分配到**唯一**的“较晚出现的 START 事件”。

---

## 2. 三种运行版本（每个 baseline 都必须实现）

框架把每个方法统一成 3 个 variant（`--variant=`）：

### 2.1 Sampling（纯采样，不显式枚举全部 join）

目标：在理论/工程上尽量做到 **时间/空间不依赖 \(|J|\)**（因为 \(|J|\) 可能是 \(\Theta(n^2)\)）。

典型结构（很多方法共享）：

1. **Phase 1（第一次 sweep）**：对每个 START 事件计算权重 \(w_e\)（该事件对应的局部 join 块大小），并累加 \(W=\sum w_e = |J|\)（或估计值）。  
2. **Phase 2（分配 slot）**：对事件按 \(w_e/W\) 建 alias；为每个样本位置 \(j\) 抽取事件并记录“该位置属于哪个事件/子块”。  
3. **Phase 3（第二次 sweep）**：再次 sweep，在事件点上做局部均匀采样并回填输出。

### 2.2 Enumerate+Sampling（先全枚举 join，再数组采样）

基线版本：把所有 \((r,s)\in J\) 真实枚举到数组 `Pairs`，再从数组中均匀抽下标（有放回）。

优点：实现直观、正确性简单。  
缺点：时间/空间**线性依赖 \(|J|\)**，大 join 容易爆内存/爆时间。

### 2.3 Adaptive（阈值自适应）

工程上最常用：

- 若 \(|J|\) 小于阈值 `J*`：走枚举分支（快，少一次 sweep）
- 若 \(|J|\) 超过阈值：丢弃已枚举结果，转入 Sampling 的两遍流程（避免 OOM）

关键：切换不允许影响输出分布；因此 Phase1 必须先做 Count，再判断阈值。

---

## 3. 代码结构总览（模块如何拼起来）

### 3.1 include/sjs：公共头文件（强可复用）

- `core/`：基础类型、RNG、计时器、统计汇总、日志等
- `geometry/`：`Point<Dim>`, `Box<Dim>`，相交/包含等 predicate；`embedding` 预留 join→range 映射
- `io/`：`Dataset<Dim>` 容器；二进制/CSV 读写；真实数据导入 stub
- `data/synthetic/`：合成数据生成器（stripe / uniform / clustered / hetero_sizes）
- `join/`：事件生成、枚举器、真值 oracle
- `sampling/`：alias、rank sampling、质量检验（chi2/KS/自相关）
- `index/`：各类索引（interval tree / aabb tree / kd tree / r tree / range tree）
- `baselines/`：统一 Baseline API + 各 baseline 的三版本实现

### 3.2 src：把“header-only 的模板世界”落地成可编译单元

主要包含：

- Dim=2 的显式实例化（避免模板膨胀、缩短编译时间）
- 工厂（`method+variant -> Baseline<2>`，用于 CLI）
- 一些 runtime 封装（如 load/generate Dataset 的统一入口）

### 3.3 apps：实验可执行程序（全 C++）

- `sjs_gen_dataset`：生成合成数据并导出 bin/csv
- `sjs_run`：单次实验运行（dataset+method+variant+t）
- `sjs_sweep`：读取 `config/sweeps/*.json` 批量扫参
- `sjs_verify`：小规模正确性 / 采样质量评估（oracle）
- `sjs_profile`：phase breakdown / 计数器输出
- `sjs_convert_realdata`：预留：真实数据转换到内部 bin 格式

---

## 4. 一条典型的 SIGMOD 实验流水线

1. **生成/准备数据**  
   - 合成：用 `sjs_gen_dataset`  
   - 真实：未来用 `sjs_convert_realdata` 预处理成 `.bin`

2. **跑 baseline**  
   - 单跑：`sjs_run`
   - 批量：`sjs_sweep`（输出 raw + summary）

3. **正确性 + 采样质量**  
   - 小规模用 `sjs_verify` 做 oracle 对比
   - 大规模做统计一致性（chi2/KS 等）

4. **性能剖析与论文图表**  
   - `sjs_profile` 看瓶颈（build/count/sample）
   - 用 `sweep_summary.csv` 画误差条 / 速度曲线

---

## 5. 如何扩展

- 新增 baseline：实现三种 variant，并在工厂注册（见 `docs/03_baselines.md`）
- 扩展到高维：添加 Dim 的显式实例化与更高维的数据结构（见 `docs/04_extending_to_high_dim.md`）
- 加入真实数据：补全 `realdata_stub` 并增加转换工具（见 `docs/05_adding_real_datasets.md`）
