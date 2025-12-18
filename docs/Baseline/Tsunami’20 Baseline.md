# Tsunami‑JSS Baseline 设计报告（90+版）

## 0. Baseline 总览（给读者的 5 行结论）

我们要在 spatial intersection join 的结果集合 $J$ 上做 **exact i.i.d. uniform sampling with replacement**。Baseline 不利用几何结构（如 plane sweep/分解），而是把“盒子相交”转为“高维点的范围过滤”，并把 **Tsunami** 当作一个强力的 learned multi‑dimensional range filtering 引擎来加速过滤。随后用三种策略实现均匀抽样：

1. **Sampling（2Pass‑RankSample）**：两遍扫描 join stream，用全局 rank 定位样本；
2. **Enumerate+Sampling（ArraySample）**：枚举全部 $J$ 到数组后抽样；
3. **Adaptive+Sampling**：输出小则枚举，大则切换到 2Pass 方案，避免内存爆炸。
    三者都能保证 **exact i.i.d. uniform with replacement**，但最坏复杂度仍会依赖 $|J|$（这是 baseline 的局限性，用来凸显你们 Ours 的优势）。

------

# 1. 问题定义与分析

## 1.1 输入、Join 结果与采样目标

给定两类轴对齐半开盒（MBR）集合：

- $R_c=\{r_{c1},\dots,r_{c n_1}\}$
- $R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\}$

每个盒子
$$
r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\quad L_i(r)<R_i(r).
$$
只关心跨集合的相交对（intersection join / filter step）：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
半开区间相交判定（逐维）：
$$
r_c\cap r_{\bar c}\neq\varnothing \iff \forall i,\ \max(L_i(r_c),L_i(r_{\bar c}))<\min(R_i(r_c),R_i(r_{\bar c})).
$$
**采样目标**：输出 $t$ 个样本
$$
Z_1,\dots,Z_t\in J
$$
满足 **i.i.d. 均匀有放回**：
$$
\Pr(Z_j=P)=\frac1{|J|}\ \ (\forall P\in J),\qquad Z_1,\dots,Z_t\ \text{相互独立}.
$$

> 这里的目标与 “Spatial Join Sampling” 语境完全一致：抽样对象是 join pair，分布必须是对 $J$ 的**严格均匀分布**（而非近似/加权）。

------

## 1.2 Baseline 的核心思想：Box‑intersection → 2d‑dim point range filtering

Tsunami 原生索引的是 **点（records/points）**，而 spatial join 的对象是 **盒子**。Baseline 采用标准转换：

对每个被索引侧盒子 $r$ 编码为 $2d$ 维点：
$$
p(r)=\big(L_1(r),\dots,L_d(r),\ R_1(r),\dots,R_d(r)\big)\in\mathbb{R}^{2d}.
$$
对查询侧盒子 $q$，与 $r$ 相交当且仅当对所有维度 $i$：

- $L_i(r) < R_i(q)$
- $R_i(r) > L_i(q)$

因此可等价为：点 $p(r)$ 落入一个 $2d$ 维正交范围查询（其中部分维度是半无限区间 $(-\infty,\cdot)$ 或 $(\cdot,+\infty)$）。

> 直观上：相交要求 “$r$ 的左端点在 $q$ 右端点左侧” 且 “$r$ 的右端点在 $q$ 左端点右侧”。这正好是对 $(L,R)$ 端点坐标的二维不等式集合（每一维各两条）。

------

## 1.3 Baseline 的研究定位：为什么“合理且强”

- **强系统 baseline**：不使用你们方法里的几何分解/模式结构，而是借助 Tsunami 的 learned multi‑dim index 做高维过滤，属于系统/DB reviewer 常认可的对照思路。
- **抽样是 exact**：通过 rank‑based 的定位方案，可实现严格 i.i.d. uniform with replacement（不是近似）。
- **局限性也“baseline 化”**：在我们将 Tsunami 当作 *range filter engine* 并通过“遍历返回的匹配结果”来计数/枚举的接口模型下，要做到严格均匀往往不可避免地触及大量输出，最坏仍会依赖 $|J|$（这正是它与 Ours 的本质差异点）。

------

# 2. 引用 Tsunami’20/’21 的关键内容

> 这里不是复述 Tsunami 全文，而是明确：baseline 使用 Tsunami 的哪些机制/假设，哪些是我们需要写进论文才能“自洽可复现”的点。

## 2.1 论文信息（建议写进引用段）

- **标题**：*Tsunami: A Learned Multi-dimensional Index for Correlated Data and Skewed Workloads*
- **作者**：Jialin Ding, Vikram Nathan, Mohammad Alizadeh, Tim Kraska
- **发表**：PVLDB Vol.14 No.2 (2021)（你们文内可继续称 Tsunami’20 以对齐工程实现来源）

------

## 2.2 Tsunami 的任务与数据模型：多维范围过滤（point filtering）

Tsunami 面向单表点数据，把 records 视为 $d$ 维空间中的点，支持典型 SQL 形式的多维范围谓词过滤（range predicates）。

> 这恰好与我们的 Box→Point 转换对齐：我们把每个盒子映射为 $2d$ 维点，然后对每个查询盒子 $q$ 发起一个 $2d$ 维范围过滤。

------

## 2.3 Tsunami 的结构与查询流程：为什么它是“强 filter 引擎”

Tsunami 是一个 **clustered、in‑memory、read‑optimized** 的 learned multi‑dim index。其核心由两部分组成：

- **Grid Tree**：按工作负载将数据空间划分为多个不重叠区域，用于缓解 query skew；
- **Augmented Grid**：在每个区域内部构建网格索引，并用相关性建模（如 functional mappings、conditional CDFs）来更好处理维度相关性。

**查询流程（我们 baseline 依赖其“过滤工作流”）**：

1. 遍历 Grid Tree 找到与查询范围相交的 regions；
2. 在每个 region 内定位相交的网格 cells，并通过 lookup table 映射到物理存储的连续区间；
3. 扫描这些物理区间内的点并检查过滤谓词，产出匹配结果。

> Baseline 将 Tsunami 视为提供接口 `Query(Index, Q) -> stream of matched ids` 的过滤引擎，并依赖它“能返回所有满足范围谓词的点 id”。

------

## 2.4 Tsunami 的线性代价模型：用于解释 baseline 的性能趋势

Tsunami 在优化 Augmented Grid 时使用一个简单的线性代价模型：
$$
\text{Time}=w_0\cdot(\#\text{cell ranges})+w_1\cdot(\#\text{scanned points})\cdot(\#\text{filtered dims}).
$$
并解释了 cell range（物理上连续的 cell 扫描段）、scanned points、filtered dims 等概念。

在我们的 baseline 中，过滤维度固定为 $\#\text{filtered dims}=2d$，因此该模型也能解释 “维度翻倍” 对性能的影响：如果 scanned points 相同，单点检查的维度因子约翻倍。

------

## 2.5 Tsunami 的实现细节：64-bit 整数与浮点缩放（影响严格不等号处理）

Tsunami 的实现与评测中使用 **64-bit integer-valued attributes**；若是浮点属性，通常会限制小数位数并按最小的 $10^k$ 进行缩放转整数。

> 这直接影响我们 baseline 的“严格不等号/半开边界”如何落地：
>
> - 若数据确实是整数（或固定精度缩放且无误差），$\pm 1$ 方案很自然；
> - 若存在浮点误差或非固定精度，必须改用 rank/二级键方案来保证严格正确（见 §3.3）。

------

# 3. 核心数据结构

## 3.1 选择索引侧与查询侧

令

- 被索引侧（data side）：$B := R_{\bar c}$
- 查询侧（workload side）：$A := R_c$

也可互换；工程上通常把更大的一侧放入索引更合算（减少构建多次查询的成本）。

------

## 3.2 Box → Point 编码（索引数据）

对每个 $r\in B$，构造 $2d$ 维点：
$$
p(r)=\big(L_1(r),\dots,L_d(r),R_1(r),\dots,R_d(r)\big).
$$
Tsunami 的“单表点过滤”即对集合 $\{p(r)\mid r\in B\}$ 建索引。

------

## 3.3 查询范围 $Q(q)$ 的构造（严格不等号与半无限范围）

对每个查询盒子 $q\in A$，相交条件等价于对所有 $i$：
$$
L_i(r) < R_i(q)\quad\land\quad R_i(r) > L_i(q).
$$
因此在 $2d$ 维空间中，查询范围可写为：

- 第 $1..d$ 维（存 $L_i(r)$）：$(-\infty,\ R_i(q))$
- 第 $d+1..2d$ 维（存 $R_i(r)$）：$(L_i(q),\ +\infty)$

但 Tsunami 的谓词通常是“闭区间/可执行”的 $a\le x\le b$。因此必须把严格不等号与半无限边界 **落地为可执行的整数边界**。我们给出两种实现方案（二选一写进论文即可）。

------

### 方案 A：整数域（或固定精度缩放）+$\pm 1$ 边界（工程直觉最强）

若所有坐标已在整数域（或你能保证浮点已按固定小数位缩放为整数且无误差）：

- $L_i(r) < R_i(q)$ 变为
  $$
  L_i(r)\le R_i(q)-1.
  $$

- $R_i(r) > L_i(q)$ 变为
  $$
  R_i(r)\ge L_i(q)+1.
  $$

半无限边界可用该维度的全局最小/最大值替代（或安全的 int64 sentinel，但必须避免溢出）。

**必须写清的前提（否则 reviewer 会质疑正确性）**：

- 数据确实是整数或固定精度缩放；
- $R_i(q)-1$、$L_i(q)+1$ 不溢出且语义正确（建议用饱和处理：若 $R_i(q)$ 已是最小值则该维度约束为空，若 $L_i(q)$ 已是最大值同理）。

------

### 方案 B：坐标压缩（rank）+ 二级键（严格可证明正确，推荐）

对每个维度属性（共有 $2d$ 个属性轴）分别做坐标压缩。关键是处理重复值与严格不等号：

- 对该维度的所有值 $\{x\}$ 构造键 $(x,\text{id})$ 或 $(x,\text{tie})$ 排序，得到严格全序；
- 用 `lower_bound/upper_bound` 将严格不等号映射到 rank 区间端点。

以某个属性 $X$ 为例：

- 约束 $X < b$ 对应 rank 区间
  $$
  X\in [\text{minRank},\ \text{lb}(b)-1],
  $$
  其中 $\text{lb}(b)$ 是 `lower_bound(b)` 的返回位置。

- 约束 $X > a$ 对应
  $$
  X\in [\text{ub}(a),\ \text{maxRank}],
  $$
  其中 $\text{ub}(a)$ 是 `upper_bound(a)` 的返回位置。

**优点**：不依赖“$\pm 1$”和浮点缩放的正确性；对任意离散/浮点数据都能保证边界不漏不重。
 **建议**：如果你们希望 baseline 在论文里“严格到可证明”，选方案 B 最稳。

------

## 3.4 Tsunami 在 baseline 中的黑盒接口（必须明确）

为了让 baseline 可复现，我们只需要 Tsunami 提供以下抽象接口：

- `BuildTsunamiIndex(Points, WorkloadQueries) -> Index`
  - `Points = { p(r) | r ∈ B }`
  - `WorkloadQueries = { Q(q) | q ∈ A }`（可用全量或采样子集）
- `Query(Index, Q(q)) -> stream/list of matched r-ids`
   返回所有满足范围过滤谓词的 $r$ 的 id。

> Tsunami 的内部结构（Grid Tree + Augmented Grid）与三步查询流程，解释了为什么它是强力 filter 引擎，但 baseline 不需要复刻它的学习细节，只需使用其查询接口。

------

## 3.5 Sampling 需要的辅助数组/结构（3 个版本共用）

**通用**（Sampling / Enumerate / Adaptive 都用到）：

- `A[1..n1]`：查询侧对象顺序（建议按 id 升序，保证确定性）
- `deg[i]`：$\deg(A[i]) = |\text{Matches}(A[i])|$
- `W = Σ deg[i] = |J|`

**Sampling（2Pass‑RankSample）额外**：

- `Ranks = [(U_j, j)]`：全局 rank 与样本位置绑定
- `Ans[1..t]`：输出数组

**Enumerate+Sampling 额外**：

- `AllPairs`：显式存储所有 join pair（可能很大）

**Adaptive 额外**：

- 阈值 `J_star`
- `mode ∈ {ENUMERATE, COUNT_ONLY}`

------

## 3.6 确定性 join stream（对 2Pass‑RankSample 至关重要）

为了让“rank → pair”的映射稳定，我们把 join 输出定义为一个**确定性序列** $P_1,\dots,P_W$：

- 外层：按 `A[1],A[2],...,A[n1]` 固定顺序遍历查询盒子；
- 内层：对固定 $q$，按 `TsunamiQuery(Q(q))` 返回的顺序依次枚举匹配的 $r$。

这要求：

1. `A` 的顺序固定；
2. Tsunami 返回的匹配结果顺序在同一 index layout 下可重复（通常等于物理扫描顺序）；
3. 对重复 key/相同点坐标的记录，物理排序或 tie‑break 必须稳定（例如按 record id 稳定排序）。

> 这点必须写进论文/报告，否则 reviewer 会指出：如果结果顺序不确定，“rank 第 k 个 pair 是谁”就不稳定，进而采样实现无法被严格证明。

------

# 4. 算法详细流程（三个版本）

## 4.0 公共预处理（所有版本共享）

### Step 0：决定索引侧/查询侧

默认用 $B=R_{\bar c}$ 建索引、$A=R_c$ 做 workload 查询（可互换）。

### Step 1：构造点表与 query 工作负载

- 点表：`Points = { p(r) | r ∈ B }`（维度 $2d$）
- workload：`WorkloadQueries = { Q(q) | q ∈ A }`（维度 $2d$）

### Step 2：Build Tsunami Index

调用 `BuildTsunamiIndex(Points, WorkloadQueries)` 得到 `Index`。
 Tsunami 的优化逻辑会基于 workload 来优化结构（Grid Tree + Augmented Grid）并组织物理布局，这也是它在 skew/correlation 下表现强的原因。

> 实验中可说明 build 是一次性离线成本：评价 join sampling 通常关注 query/sampling runtime；但你们也可以在附录报告 build 时间以更透明。

------

## 4.1 版本 I：Sampling（TSUNAMI‑2Pass‑RankSample）

### 4.1.1 核心思想

把所有 join pair 视为确定性序列 $P_1,\dots,P_W$。
 生成 $t$ 个 i.i.d. 均匀 rank：$U_j \sim \mathrm{Unif}\{1,\dots,W\}$。
 第二遍扫描时定位并输出 $P_{U_j}$。

------

### 4.1.2 Pass 1：计数（得到 deg 与 W）

对每个 $q\in A$：执行 `TsunamiQuery(Q(q))` 并计数：
$$
\deg(q)=|\text{Matches}(q)|,\qquad W=\sum_{q\in A}\deg(q).
$$
若 $W=0$，返回空。

> 注意：此处“计数”是在我们 baseline 的接口模型里通过遍历 query 返回结果完成，因此至少会触及每个实际匹配一次。

------

### 4.1.3 生成 ranks

对 $j=1..t$：独立生成
$$
U_j\stackrel{i.i.d.}{\sim}\mathrm{Unif}\{1,\dots,W\}
$$
并构造 `Ranks = sort([(U_j, j)])`。

------

### 4.1.4 Pass 2：定位输出（跳过无关 query + 早停）

维护全局已扫过 pair 数 `g`。对每个 $q\in A$：

- 若当前最小 rank $U > g+\deg(q)$，则该 $q$ 的 block 内没有样本 rank，**跳过执行 TsunamiQuery**，仅令 `g += deg(q)`。
- 否则执行 TsunamiQuery(Q(q))，并在匹配流中用局部计数 `c` 找到需要的第 $u=U-g$ 个匹配，回填 `Ans[pos] = (q,r)`。
- 同一 rank 可出现多次，输出同一对多次，正是有放回抽样。

**可复现伪代码**（建议直接放入报告）：

```
TSUNAMI-2Pass-RankSample(A, Index, t):

# Pass 1: count degrees
for q in A (fixed order):
    deg[q] = 0
    for r in TsunamiQuery(Index, Q(q)):
        deg[q]++

W = sum_q deg[q]
if W == 0: return empty

# draw i.i.d. global ranks with replacement
for j = 1..t:
    U[j] ~ UniformInt(1, W)
Ranks = sort([(U[j], j)] by U)

# Pass 2: locate outputs
Ans[1..t]
g = 0          # global offset = number of pairs before current q-block
p = 1          # pointer in Ranks (sorted)

for q in A (same order):
    if p > t: break

    # skip if no rank falls inside (g, g+deg[q]]
    if Ranks[p].U > g + deg[q]:
        g += deg[q]
        continue

    # collect all local positions needed in this q-block
    Local = []
    while p <= t and Ranks[p].U <= g + deg[q]:
        u = Ranks[p].U - g              # u in [1..deg[q]]
        Local.append((u, Ranks[p].pos)) # pos = original sample index
        p++

    sort Local by u
    umax = Local.last.u

    # stream matches and stop early at umax-th match
    c = 0
    i = 1
    for r in TsunamiQuery(Index, Q(q)):
        c++
        while i <= |Local| and Local[i].u == c:
            Ans[Local[i].pos] = (q, r)
            i++
        if c == umax: break

    g += deg[q]

return Ans
```

**必须写清的正确性条件（工程注意点）**：

1. `A` 顺序固定；
2. TsunamiQuery 的匹配输出顺序可重复（建议等于物理扫描顺序）；
3. 严格不等号必须按 §3.3 的方案实现，否则边界 pair 会误算/漏算。

------

## 4.2 版本 II：Enumerate+Sampling（TSUNAMI‑Enumerate‑ArraySample）

### 4.2.1 核心思想

直接枚举全部 join 输出对 $J$ 到数组 `AllPairs`，然后在数组上做 i.i.d. 均匀下标采样。实现最直观，但空间 $\Theta(|J|)$。

### 4.2.2 伪代码

```
TSUNAMI-Enumerate-ArraySample(A, Index, t):

AllPairs = []
for q in A (fixed order):
    for r in TsunamiQuery(Index, Q(q)):
        AllPairs.append((q, r))

W = |AllPairs|
if W == 0: return empty

for j = 1..t:
    idx ~ UniformInt(0, W-1)
    Ans[j] = AllPairs[idx]

return Ans
```

------

## 4.3 版本 III：Adaptive+Sampling（TSUNAMI‑Adaptive）

### 4.3.1 核心思想

用阈值 $J_\star$ 自适应选择：

- 若 $|J|\le J_\star$：直接 Enumerate+Sampling（只需一遍扫描 + 数组采样）
- 若 $|J|> J_\star$：走 2Pass‑RankSample（节省内存；第二遍仅对命中 query 执行）

------

### 4.3.2 实现现实点：Tsunami 未必提供 CountOnly API

因此给出两种写法：

- **写法 A（理想）**：若 Tsunami 提供 `CountOnly(Q)`，可先 count 再决定是否枚举。
- **写法 B（通用、可落地）**：Query 流式返回时永远 count，同时只在未超阈值时追加到 `AllPairs`；一旦超过阈值立刻切换并释放 `AllPairs`，后续只计数。

下面给出写法 B（最通用）。

------

### 4.3.3 写法 B：一次遍历完成精确计数 +（可选）小规模枚举

**Phase 1：count + maybe enumerate until threshold**

- `mode = ENUMERATE`
- `AllPairs = []`（最多 $J_\star$）
- `W = 0`
- 同时维护 `deg[q]`

遍历 $q\in A$ 并流式扫描匹配 $r$：

- 永远：`deg[q]++ ; W++`
- 若 `mode==ENUMERATE`：
  - 若 `|AllPairs| < J_star`：append `(q,r)`
  - 否则：切换 `mode=COUNT_ONLY`，释放 `AllPairs`

**Phase 1 后分支**：

- 若 `mode==ENUMERATE`：说明 $W\le J_\star$，`AllPairs` 已完整枚举 $J$，直接数组采样
- 若 `mode==COUNT_ONLY`：说明 $W>J_\star$，转入 Sampling 分支：
  - 生成 ranks（用已知 $W$）
  - 执行“仅 Pass2 的定位输出”（deg 已在 Phase1 得到，无需再计数一遍）

------

### 4.3.4 伪代码（写法 B）

```
TSUNAMI-Adaptive(A, Index, t, J_star):

mode = ENUMERATE
AllPairs = []
W = 0

# Phase 1: exact count; enumerate only if small
for q in A:
    deg[q] = 0
    for r in TsunamiQuery(Index, Q(q)):
        deg[q]++
        W++

        if mode == ENUMERATE:
            if |AllPairs| < J_star:
                AllPairs.append((q, r))
            else:
                mode = COUNT_ONLY
                AllPairs.clear(); AllPairs.shrink_to_fit()

if W == 0: return empty

if mode == ENUMERATE:
    # AllPairs contains the full join result
    for j = 1..t:
        idx ~ UniformInt(0, W-1)
        Ans[j] = AllPairs[idx]
    return Ans

# mode == COUNT_ONLY: do rank-sampling branch using deg and W
return TSUNAMI-Pass2-LocateUsingDeg(A, Index, t, deg, W)
```

------

# 5. Adaptive 阈值 $J_\star$ 的选择策略

阈值本质上是在 “**省一遍扫描**” 与 “**避免存储爆炸**” 之间折中：

- 小 $|J|$：枚举 + 数组采样更快更简单；
- 大 $|J|$：必须切换到 rank‑sampling，避免 $\Theta(|J|)$ 内存。

给出三层策略（建议都写进论文，层次清晰）。

------

## 5.1 内存硬约束（必须写）

设给 `AllPairs` 的内存预算为 `MemBudget·ρ`（$0<\rho<1$），每个 pair 的存储开销为 `sizeof(Pair)` 字节，则：
$$
J_\star^{\text{mem}}
=\left\lfloor\frac{\rho\cdot \text{MemBudget}}{\text{sizeof(Pair)}}\right\rfloor,
\qquad
J_\star \le J_\star^{\text{mem}}.
$$
这保证不会 OOM。

------

## 5.2 时间权衡（建议写：避免大分支被额外写入拖垮）

在写法 B 中，一旦触发切换，最多额外“写入/尝试存储” $J_\star$ 个 pair 就会丢弃，因此大分支多了一个 $O(J_\star)$ 的写入开销。为了不让它压倒主要开销，常见建议是：
$$
J_\star^{\text{time}} = C\cdot t
\quad\text{或}\quad
J_\star^{\text{time}} = C_1\cdot t + C_2\cdot n_1.
$$
直觉：当只需要 $t$ 个样本时，不值得为省第二遍扫描而写入过多 pair。

最终可取：
$$
J_\star=\min\bigl(J_\star^{\text{mem}},\ J_\star^{\text{time}}\bigr).
$$

------

## 5.3 工程标定（最推荐：交叉点拟合）

做 benchmark 标定（与你们目前草稿一致，但这里把“怎么写进论文”说完整）：

1. 固定若干代表性配置 $(n_1,n_2,d,t)$ 与数据分布；
2. 分别测：纯 Enumerate+Sampling 与纯 2Pass‑RankSample 的耗时随 $|J|$ 变化；
3. 找到两者耗时相当的交叉点 $|J_{\text{cross}}|$；
4. 取 $J_\star \approx \alpha |J_{\text{cross}}|$（$\alpha\in(0,1)$ 略小），并受 $J_\star^{mem}$ 截断。

这样写的好处是：阈值不是拍脑袋，而是“与实现/硬件相关的经验最优点”。

------

# 6. 算法分析（正确性 + 复杂度；三版本都包含）

## 6.1 正确性基础：Box→Point 过滤的充要等价

**命题**：对任意半开盒 $q,r$，有
$$
q\cap r\neq\varnothing
\iff
\forall i,\ L_i(r) < R_i(q)\ \land\ R_i(r) > L_i(q).
$$
**证明（逐维）**：
 半开区间相交等价于
$$
\max(L_i(q),L_i(r))<\min(R_i(q),R_i(r)).
$$
这等价于同时满足：

- $L_i(r)<R_i(r)$（否则 $L_i(r)\ge R_i(q)\Rightarrow \max(\cdot)\ge R_i(q)\ge \min(\cdot)$ 不相交）
- $R_i(r)>L_i(q)$（否则 $R_i(r)\le L_i(q)\Rightarrow \min(\cdot)\le L_i(q)\le \max(\cdot)$ 不相交）
   反向同理。合并所有维度即可。∎

因此，用 Tsunami 对 $p(r)$ 做 $2d$ 维范围过滤得到的匹配集合，与真实相交集合一一对应。

------

## 6.2 版本 I：2Pass‑RankSample 的正确性（exact i.i.d. uniform）

把 join 输出定义为确定性序列 $P_1,\dots,P_W$。算法生成 $U_j\sim \mathrm{Unif}\{1..W\}$，并输出 $Z_j=P_{U_j}$。

对任意固定 pair $P_k\in J$：
$$
\Pr(Z_j=P_k)=\Pr(U_j=k)=\frac1W=\frac1{|J|}.
$$
且 $U_1,\dots,U_t$ 独立生成，因此 $Z_1,\dots,Z_t$ **i.i.d. 均匀有放回**。∎ 

**注意**：该证明依赖 §3.6 的“确定性 join stream”条件（外层顺序与 Tsunami 返回顺序可重复）。

------

## 6.3 版本 II：Enumerate+Sampling 的正确性（exact i.i.d. uniform）

`AllPairs` 是对 $J$ 的一次完整枚举（每条 `(q,r)` 出现一次）。在数组上对下标做 i.i.d. 均匀采样显然得到 $J$ 上的 i.i.d. 均匀有放回样本。∎ 

------

## 6.4 版本 III：Adaptive 的正确性（切换不改变分布）

Adaptive 只是在两种**已正确**的方法之间选择：

- 若未切换（$W\le J_\star$）：等价 Enumerate+Sampling；
- 若切换（$W> J_\star$）：丢弃已枚举结果但保留精确 `deg[]` 与 `W`，随后执行与 2Pass‑RankSample 相同的 rank‑定位输出（无需再做 Pass1 计数），因此与 Sampling 分支完全等价。∎ 

------

## 6.5 复杂度：记号、Tsunami 代价模型与三版本结论

### 6.5.1 记号

- $n_1=|A|$, $n_2=|B|$
- $\deg(q)=|\text{Matches}(q)|$
- $W=|J|=\sum_{q\in A}\deg(q)$

### 6.5.2 Tsunami 的单次查询代价模型（解释趋势用）

Tsunami 给出线性模型：
$$
\text{Time}(q)\approx w_0\cdot C_q + w_1\cdot S_q\cdot( \#\text{filtered dims}),
$$
其中 $C_q$ 是 cell ranges 数，$S_q$ 是扫描点数，$\#\text{filtered dims}$ 为过滤维度数。

在 baseline 中：$\#\text{filtered dims}=2d$。

------

### 6.5.3 版本 I：Sampling（2Pass‑RankSample）

**时间**

- Build：记为 $T_{\text{build}}$（离线优化+重排等一次性成本）

- Pass1 计数：
  $$
  T_1 \approx \sum_{q\in A}\big(w_0 C_q + w_1 S_q\cdot 2d\big).
  $$
  在本 baseline 的接口模型中，计数通过遍历返回结果完成，因此至少需要触及每个匹配一次，从而有 baseline‑特定下界：
  $$
  T_1=\Omega(W)
  \quad (\text{在“遍历 Query 输出计数”的实现模型下}).
  $$
  这不是问题本身的下界，而是本 baseline 的限制。

- ranks 排序：$O(t\log t)$

- Pass2 定位输出：只对命中的 $A'\subseteq A$ 执行 query；最坏仍可达到 $\Omega(W)$，但当 $t$ 小且命中稀疏时通常明显小于全量 query。

综合：
$$
T_{\text{Sampling}} = T_{\text{build}} + T_1 + O(t\log t) + T_2.
$$
**额外空间（不含 Tsunami index 本体）**

$$
S_{\text{Sampling,extra}}=O(n_1+t).
$$

---

### 6.5.4 版本 II：Enumerate+Sampling 

**时间**   枚举必须写出全部 \(W\) 个 pair：\(\Omega(W)\)，加上 Tsunami 扫描开销： 

$$
T_{\text{Enum}}\approx T_{\text{build}} + \sum_{q\in A}\big(w_0 C_q + w_1 S_q\cdot 2d\big) + \Theta(W) + O(t).
$$


**空间**
$$
S_{\text{Enum,extra}}=\Theta(W)+O(t).
$$
这是该版本最大风险（密集 join 时 OOM）。

------

### 6.5.5 版本 III：Adaptive+Sampling

分两种情况：

- **Case A（未切换，$W\le J_\star$）**：等价 Enumerate+Sampling，但 $W$ 被阈值限制，因此空间受控：
  $$
  S_{\text{Adaptive,peak}} = O(J_\star+t).
  $$
- **Case B（触发切换，$W>J_\star$）**：Phase1 仍要完成精确计数，同时最多额外写入 $O(J_\star)$ 个 pair（随后丢弃），再进入 rank‑sampling 的定位输出：
  $$
  T_{\text{Adaptive}} = T_{\text{build}} + T_1 + O(J_\star) + O(t\log t) + T_2,
  $$

  $$
  S_{\text{Adaptive,peak}} = O(n_1+t+J_\star).
  $$

  若按 §5 的建议取 $J_\star$（受内存与 $t$ 控制），即可保证大分支不会被额外写入拖垮。

------

# 附录（可选但强烈建议加）：实现检查清单（让 reviewer 一眼安心）

1. **严格不等号处理**：明确选方案 A 或 B，并写清前提。若用 A，说明数据是 int64 或固定精度缩放；否则用 B。
2. **确定性输出顺序**：固定 `A` 顺序；Tsunami 返回顺序需可重复；重复点需稳定 tie‑break。
3. **计数溢出**：`deg`、`W` 建议用 64-bit/128-bit 计数（$|J|$ 可能非常大）。
4. **随机数可复现**：固定 seed；`UniformInt(1,W)` 需支持大 $W$。
5. **Pass2 早停**：仅在你明确定义了“rank 对应匹配输出顺序”且 TsunamiQuery 流式可中断时使用；否则可以不早停以避免实现差异。