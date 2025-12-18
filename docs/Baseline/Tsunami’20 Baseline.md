# 1. 问题定义与分析（对齐 Spatial Join Sampling）+ 引用 Tsunami’20

## 1.1 输入、join 结果与采样目标

给定两个集合（两张表）的轴对齐半开盒子（MBR）：

- $R_c=\{r_{c1},\dots,r_{c n_1}\}$
- $R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\}$

每个盒子
$$
r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\quad L_i(r)<R_i(r).
$$
只关心跨集合相交对（intersection join / filter step）：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
半开区间相交判定为（逐维）：
$$
r_c\cap r_{\bar c}\neq\varnothing \iff \forall i,\ \max(L_i(r_c),L_i(r_{\bar c}))<\min(R_i(r_c),R_i(r_{\bar c})).
$$
（这与您给出的对齐描述一致。） 

**采样目标（你们研究的标准目标）**：输出 $t$ 个样本
$$
Z_1,\dots,Z_t\in J
$$
满足 **i.i.d. 均匀有放回**：
$$
\Pr(Z_j=P)=\frac1{|J|}\ (\forall P\in J),\qquad Z_1,\dots,Z_t\ \text{相互独立}.
$$


------

## 1.2 Baseline 的核心思路：把“盒子相交”转为“点的高维范围过滤”

Tsunami 本质索引的是**点记录**，而我们的对象是**盒子**。Baseline 采用标准转换：

对每个盒子 $r\in R_{\bar c}$ 编码为 $2d$ 维点
$$
p(r)=\big(L_1(r),\dots,L_d(r),\ R_1(r),\dots,R_d(r)\big)\in\mathbb{R}^{2d}.
$$
盒子 $q\in R_c$ 与 $r$ 相交当且仅当对所有维度 $i$：

- $L_i(r) < R_i(q)$
- $R_i(r) > L_i(q)$

这等价于对点 $p(r)$ 做一个 $2d$ 维正交范围查询（部分维度是 $(-\infty,\cdot)$、$(\cdot,+\infty)$ 形式）。 

因此 join 的“枚举/计数”可以通过“对每个 $q\in R_c$ 做一次 Tsunami range filter 得到匹配集合”完成。

------

## 1.3 引用 Tsunami’20 的关键内容（我们 baseline 依赖哪些“论文要点”）

**论文信息**：*Tsunami: A Learned Multi-dimensional Index for Correlated Data and Skewed Workloads*（Ding, Nathan, Alizadeh, Kraska）。 

Baseline 将 Tsunami 当作“range filtering 引擎”，主要继承以下内容（与实验/实现强相关）：

1. **目标查询形态是多维范围谓词过滤**（典型 SQL 形式的多维范围过滤）。
2. **结构与查询流程**：Tsunami 由 **Grid Tree + Augmented Grid** 组成；查询时先定位相交 region，再在 region 内定位 cell ranges，并通过 lookup table 映射到物理存储区间，最后扫描并检查过滤谓词。
3. **线性代价模型**（用于分析/解释性能趋势）：

$$
\text{Time}=w_0(\#\text{cell ranges})+w_1(\#\text{scanned points})(\#\text{filtered dims})
$$

其中 filtered dims 在我们 baseline 中固定为 $2d$。

4. **工程实现习惯**：实验实现使用 64-bit integer 属性；若是浮点属性会缩放为整数。这个点对我们处理“严格不等号/半开边界”很关键。

------

## 1.4 Baseline 的研究定位（为什么“合理且强”）

- 这是一个**强 DB/系统派 baseline**：它不引入你们论文的特殊几何结构，而是借助 Tsunami’20 的 learned multi‑dim index 做 range filtering（很多 reviewer 会认可）。
- 它能做到 **exact i.i.d. uniform with replacement**（不是近似、不是加权近似）。
- 其缺点也很“baseline 化”：为了保证严格均匀，仍不可避免要“触及/计数”大量输出，最坏仍 $\Omega(|J|)$。

------

# 2. 核心数据结构

## 2.1 点编码与查询范围（Box → Point）

令被索引侧 $B:=R_{\bar c}$，查询侧 $A:=R_c$（也可互换，通常把更大的一侧放进 Tsunami 更合算）。

**点编码：**
$$
p(r)=\big(L_1(r),\dots,L_d(r),R_1(r),\dots,R_d(r)\big)\in\mathbb{R}^{2d}.
$$


**查询范围（理论形式）：** 对 $q\in A$，相交条件等价于：
$$
\forall i:\quad L_i(r) < R_i(q)\ \land\ R_i(r) > L_i(q).
$$
对应 $2d$ 维范围：

- 第 1..d 维（存 $L_i(r)$）：$(-\infty,\ R_i(q))$
- 第 d+1..2d 维（存 $R_i(r)$）：$(L_i(q),\ +\infty)$



------

## 2.2 严格不等号/半开边界的“可执行化”（强烈建议写进实现细节）

Tsunami 在实现中偏好整数域（64-bit integer），浮点会缩放成整数。
 因此 baseline 需要把严格不等号转成可执行的整数边界。

给出两种实现方案（你们选其一写进论文即可）：

### 方案 A：缩放到整数 + $\pm 1$ 边界（最贴近 Tsunami’20 工程习惯）

若坐标已是整数、盒子为半开 $[L,R)$，则：

- $L_i(r) < R_i(q)$ 可写成 $L_i(r)\le R_i(q)-1$
- $R_i(r) > L_i(q)$ 可写成 $R_i(r)\ge L_i(q)+1$

浮点时先乘 $10^k$ 并取整，再用 $\pm 1$ 做“严格不等号”边界修正。

同时把 $(-\infty,\cdot)$、$(\cdot,+\infty)$ 用数据域全局最小/最大值（或 int64 最小/最大但要避免溢出）替代。

### 方案 B：坐标压缩（rank）+ 二级键（更稳健，避免“$\pm 1$”假设）

对每个维度的所有端点值做排序，并用二级键 $(\text{value},\text{id})$ 做 tie‑break，把严格 $<$、$>$ 映射为 rank 区间端点的 `lower_bound/upper_bound`。
 优点：对任意离散/浮点数据都不会因为缩放误差产生边界错判。

> 如果你们打算“严谨到可证明正确”，方案 B 往往更讨 reviewer 喜欢；若更偏工程可复现，方案 A 更直接。

------

## 2.3 Tsunami 索引（作为黑盒接口描述即可）

你们不需要复述 Tsunami 的所有 learned 细节，但**需要明确它在 baseline 中提供什么接口**：

- `BuildTsunamiIndex(Points, WorkloadQueries) -> Index`
  - Points：$\{p(r)\mid r\in B\}$
  - Workload：$\{Q(q)\mid q\in A\}$
     Tsunami 的 index creation 包含离线优化与数据重排（clustered layout）。
- `Query(Index, Q(q)) -> stream of matched r-ids`
   查询流程由 Grid Tree + Augmented Grid 定位 cell ranges 并扫描过滤。

> 实现上：把 Tsunami 当作“高维 range filter engine”即可。

------

## 2.4 Sampling 需要的辅助结构

对三个版本分别需要：

### 通用

- `A[1..n1]`：查询侧对象顺序（建议按对象 id 升序，保证确定性）。
- `deg[i]`：第 $i$ 个查询对象 $q=A[i]$ 的匹配数 $|\text{Matches}(q)|$。
- `W = Σ deg[i] = |J|`。

### Sampling 版本额外

- `Ranks = [(U_j, j)]`：全局秩（rank）与样本位置绑定，并按 rank 排序。
- `Ans[1..t]`：输出样本数组。

### Enumerate+Sampling 版本额外

- `AllPairs`：显式存储 join 输出对（可能很大）。

### Adaptive 版本额外

- `J_*`：阈值（见第 4 节）
- `mode ∈ {ENUMERATE, COUNT_ONLY}`

------

# 3. 算法详细流程（三个版本）

下面把 Tsunami‑JSS 扩成三个版本，并给出足够“可复现”的确定性细节与伪代码。

------

## 3.0 三个版本共享的预处理

### Step 0：选择哪一侧建 Tsunami

默认：

- $A := R_c$（查询侧 / workload）
- $B := R_{\bar c}$（被索引侧 / data）

通常把更大的一侧放到 Tsunami 更划算。

### Step 1：构造点表与工作负载查询集合

- 点表：$\{p(r)\mid r\in B\}\subset\mathbb{R}^{2d}$ 
- 工作负载：$\{Q(q)\mid q\in A\}$（由 2.1/2.2 定义）

### Step 2：建索引

用 Tsunami 的离线优化（Grid Tree + Augmented Grid）构建 clustered index。

------

## 3.1 版本 I：Sampling（TSUNAMI‑2Pass‑RankSample）

### 核心思想

先把全体 join 对看成一个确定性“流”（stream）：

- 外层：按固定顺序遍历 $q\in A$
- 内层：对每个 $q$，以 Tsunami Query 返回的固定扫描/返回顺序遍历 $\text{Matches}(q)$

从而得到序列 $P_1,\dots,P_W$（$W=|J|$）。
 Sampling 就是：抽 $t$ 个 i.i.d. 均匀 rank，再在第二遍扫描中“命中这些 rank”。

### 算法流程

**Pass 1（计数）：**
 对每个 $q\in A$ 运行一次 Tsunami range query，并计数输出数目得到：
$$
\deg(q)=|\text{Matches}(q)|,\qquad W=\sum_{q\in A}\deg(q).
$$
**生成 ranks：**
$$
U_j\stackrel{i.i.d.}{\sim}\mathrm{Unif}\{1,\dots,W\},\quad j=1..t
$$
并把 `(U_j, j)` 按 `U_j` 升序排序。

**Pass 2（定位输出）：**
 维护全局已扫过对数 `g`。遍历 $q\in A$：

- 若当前最小 rank $U$ 满足 $U>g+\deg(q)$：说明该 $q$ 的 block 内没有任何样本位置，**跳过执行 query**（这是工程上节省 pass2 的关键）。
- 否则执行 query，流式遍历匹配，同时用局部计数 `c` 命中需要的局部位置 $u = U-g$。

> 同一个 rank 出现多次会输出同一对多次，这正是“有放回”的体现。

### 伪代码（可直接放文稿）

（以下伪代码的“计数+排序+二遍命中”结构与你们已写版本一致，我把工程细节补齐了。）

```
TSUNAMI-2Pass-RankSample(A, B, t):

# Preprocess: build Tsunami index on points p(r), r in B, with workload queries from A.

# Pass 1: count
for q in A (fixed order):
    deg[q] = 0
    for r in TsunamiQuery(Q(q)):   # stream results
        deg[q]++

W = sum_q deg[q]
if W == 0: return empty

# draw global ranks with replacement
for j = 1..t:
    U[j] ~ UniformInt(1, W)
Ranks = sort([(U[j], j)] by U)

# Pass 2: locate outputs
Ans[1..t]
g = 0               # processed pairs so far
p = 1               # pointer in Ranks

for q in A (same order):
    if p > t: break

    # skip this query if no rank falls in (g, g+deg[q]]
    if Ranks[p].U > g + deg[q]:
        g += deg[q]
        continue

    # collect all needed local positions within this q-block
    Local = []   # list of (u, original_pos)
    while p <= t and Ranks[p].U <= g + deg[q]:
        u = Ranks[p].U - g    # 1..deg[q]
        Local.append((u, Ranks[p].pos))
        p++

    Local.sort by u
    umax = Local.last.u

    # run Tsunami query and stop early at umax
    c = 0
    i = 1   # pointer in Local
    for r in TsunamiQuery(Q(q)):
        c++
        while i <= |Local| and Local[i].u == c:
            Ans[ Local[i].pos ] = (q, r)
            i++
        if c == umax: break

    g += deg[q]

return Ans
```

### 关键实现注意点（保证严格正确）

1. **固定 A 的顺序**（比如按 id），否则 “stream 定义”不确定。
2. **Tsunami Query 的扫描/返回顺序要可重复**：通常用“物理存储顺序”即可（不要在返回结果上再做随机 shuffle）。
3. **半开/严格不等号**必须通过 2.2 的方案处理，否则会出现边界 pair 被误算/漏算。

------

## 3.2 版本 II：Enumerate+Sampling（TSUNAMI‑Enumerate‑ArraySample）

### 核心思想

直接用 Tsunami 把 join 输出对 $J$ **完整枚举到数组** `AllPairs`，然后在数组上做均匀 i.i.d. 采样（最朴素、常数小，但空间爆）。这与你们“Enumerate+Sampling”版本的精神一致：**先枚举，再随机取样**。

### 算法流程

1. 枚举：
   - 对每个 $q\in A$：执行 TsunamiQuery(Q(q))，对每个返回 $r\in B$：append `(q,r)` 到 `AllPairs`
2. 抽样：
   - 对 $j=1..t$：`idx ~ Unif(0..|AllPairs|-1)`，输出 `AllPairs[idx]`

### 伪代码

```
TSUNAMI-Enumerate-ArraySample(A, B, t):

AllPairs = []
for q in A (fixed order):
    for r in TsunamiQuery(Q(q)):
        AllPairs.append((q, r))

W = |AllPairs|
if W == 0: return empty

for j = 1..t:
    idx ~ UniformInt(0, W-1)
    Ans[j] = AllPairs[idx]

return Ans
```

------

## 3.3 版本 III：Adaptive+Sampling（TSUNAMI‑Adaptive）

### 核心思想

用阈值 $J_\star$ 自动选择：

- 若 $|J|\le J_\star$：走 **Enumerate+Sampling**（只需一遍遍历 + 数组采样）
- 若 $|J|>J_\star$：走 **Sampling（两遍 rank）**（节省内存；第二遍只对命中 query 再查）

> 这正是你们 Adaptive 框架的 baseline 化版本：小输出用枚举，大输出用两遍采样。

------

### 重要现实点：Tsunami 不一定天然提供“只 Count 不枚举对象 id”的接口

因此 Adaptive 的 Phase1 实现有两种写法：

- **写法 A（理想，若 Tsunami 能提供 CountOnly）**：先 Count 再决定是否枚举，保证枚举量严格 $\le J_\star$（和你们别的 Adaptive 文档一致的“先 Count 再判断”哲学）。
- **写法 B（更普适）**：Query 一次流式返回时同时 count；若还处于 ENUMERATE 且 `AllPairs.size < J_star` 就 append；一旦达到阈值立刻切换并清空 `AllPairs`。这样不会 OOM，且额外写入最多 $J_\star$ 条。

下面给出 **写法 B** 的“可落地 baseline”（最常用，也最省事）。

------

### TSUNAMI‑Adaptive（写法 B）流程

**Phase 1：一次遍历（永远 count，必要时 append 直到阈值）**

初始化：

- `mode = ENUMERATE`
- `AllPairs` 预留容量 `J_star`（可选）
- `W = 0`

对每个 $q\in A$：

- `deg[q]=0`
- 遍历 TsunamiQuery(Q(q)) 返回的每个 $r$：
  - `deg[q]++ ; W++`（永远做）
  - 若 `mode==ENUMERATE`：
    - 若 `AllPairs.size < J_star`：append `(q,r)`
    - 否则（达到阈值）：
       `mode = COUNT_ONLY`；清空并释放 `AllPairs`（例如 `clear+shrink_to_fit`）；后续不再 append

遍历结束后我们总能得到精确的 `deg[]` 与 $W=|J|$。

**Phase 1 后分支：**

- 若 `mode==ENUMERATE`（说明 $W \le J_\star$，且 `AllPairs` 已完整）：直接数组采样输出（等价 Enumerate+Sampling）
- 若 `mode==COUNT_ONLY`（说明 $W > J_\star$）：转入 Sampling 的 ranks + Pass2 输出

------

### 伪代码（写法 B）

```
TSUNAMI-Adaptive(A, B, t, J_star):

mode = ENUMERATE
AllPairs = []
W = 0

# Phase 1: count + maybe store until threshold
for q in A:
    deg[q] = 0
    for r in TsunamiQuery(Q(q)):
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
    # W <= J_star, AllPairs contains full J
    for j = 1..t:
        idx ~ UniformInt(0, W-1)
        Ans[j] = AllPairs[idx]
    return Ans

# else mode == COUNT_ONLY: do Sampling branch (ranks + Pass2)
return TSUNAMI-2Pass-RankSample-UsingDeg(A, B, t, deg, W)
```

> 其中 `TSUNAMI-2Pass-RankSample-UsingDeg` 就是 3.1 的 Pass2，但 Pass1 已经得到 deg 与 W 了，所以无需再跑一遍计数。

------

# 4. Adaptive 版本阈值 $J_\star$ 的选择策略（给出可写进论文的规则）

阈值的本质：**允许枚举存储多少 join pairs，来换取“少一遍查询扫描”**。

下面给出三个层次的策略（从必须到可选）。

------

## 4.1 内存硬约束（必须满足）

假设每个 pair 只存两个对象 id（或指针），`sizeof(Pair)` 字节；给 `AllPairs` 的预算为 `MemBudget·ρ`（$0<\rho<1$）：
$$
J_\star^{\text{mem}}
=\left\lfloor\frac{\rho\cdot \text{MemBudget}}{\text{sizeof(Pair)}}\right\rfloor.
$$
实际取值必须满足：
$$
J_\star \le J_\star^{\text{mem}}.
$$

> 这条规则可以防止 baseline 在密集 join 上直接 OOM（尤其 Tsunami 返回结果很大时）。这一条是“必须写”的。

------

## 4.2 时间权衡（建议写：让大分支别被额外写入拖垮）

Tsunami‑Adaptive（写法 B）在“大分支”里最多额外写入/尝试存储 $J_\star$ 个 pair，随后立即丢弃并转为 COUNT_ONLY。

因此大分支的额外开销是 $O(J_\star)$ 级别（主要是内存写入），你可以给出一个建议：
$$
J_\star^{\text{time}} = C_1\cdot t + C_2\cdot n_1
\quad\text{或更粗：}\quad
J_\star^{\text{time}}=C\cdot t,
$$
直觉是：当你只需要 $t$ 个样本时，没必要为了省第二遍而“写入过大规模的 AllPairs”。

最终取：
$$
J_\star=\min\bigl(J_\star^{\text{mem}},\ J_\star^{\text{time}}\bigr).
$$

------

## 4.3 工程标定（最推荐：交叉点拟合）

这与“通用 Adaptive 框架”的常见做法一致：

1. 在代表性数据分布与参数 $(n_1,n_2,d,t)$ 下跑：
   - “纯 Enumerate+Sampling”
   - “纯 Sampling（两遍 rank）”
2. 找到两者耗时相当的输出规模 $|J_{\text{cross}}|$
3. 令 $J_\star$ 略小于该交叉点（并受 $J_\star^{mem}$ 截断）

这样写进论文会非常自然：阈值来自 **benchmark 标定**，不是拍脑袋。

------

# 5. 算法分析（正确性、复杂度；三版本都包含）

## 5.1 正确性：Box→Point 范围过滤的等价性

对任意两盒子 $q,r$，在第 $i$ 维相交当且仅当：
$$
\max(L_i(q),L_i(r))<\min(R_i(q),R_i(r))
\iff L_i(r)<R_i(q)\ \land\ R_i(r)>L_i(q).
$$
将所有维度合并，就得到“点 $p(r)$”落入“查询范围 $Q(q)$”的充要条件。该等价转换是 Tsunami‑JSS 正确性的几何基础。

------

## 5.2 正确性：Sampling（两遍 rank）是 exact i.i.d. uniform

把 join 输出定义成确定性流：

- 先按固定顺序遍历 $q\in A$
- 再按 Tsunami query 内部扫描/返回顺序遍历 $\text{Matches}(q)$

得到 $P_1,\dots,P_W$，其中 $W=|J|=\sum_q \deg(q)$。

Sampling 版本生成 $U_j\sim \text{Unif}\{1..W\}$，并令 $Z_j=P_{U_j}$。于是对任意固定 pair $P_k\in J$：
$$
\Pr(Z_j=P_k)=\Pr(U_j=k)=\frac1W=\frac1{|J|}.
$$
因为 $U_1,\dots,U_t$ 独立生成，所以 $Z_1,\dots,Z_t$ i.i.d.。

------

## 5.3 正确性：Enumerate+Sampling 是 exact i.i.d. uniform

`AllPairs` 是对 $J$ 的一次枚举（每个 $q$ 只对应其匹配 $r$，不会跨 $q$ 重复），在数组上做独立均匀下标采样显然得到 $J$ 上 i.i.d. 均匀有放回样本。

------

## 5.4 正确性：Adaptive 不破坏正确性

Adaptive 只是根据 $W$（或触发阈值）在两种**已正确**的方法之间选择：

- 若未切换：等价于 Enumerate+Sampling
- 若切换：等价于 Sampling（两遍 rank）

切换行为不改变 $Q(q)$ 的定义、不改变 Tsunami 返回的匹配集合，因此不会引入偏差。

------

## 5.5 复杂度：通用记号与 Tsunami 代价模型

记：

- $n_1=|A|$，$n_2=|B|$
- $\deg(q)=|\text{Matches}(q)|$
- $W=|J|=\sum_q \deg(q)$

Tsunami 提供的可用于解释实验趋势的代价模型（单 query）：
$$
\text{Time}(q)\approx w_0\cdot C_q + w_1\cdot S_q\cdot (2d),
$$
其中 $C_q$ 是 cell ranges 数，$S_q$ 是扫描点数，过滤维度固定为 $2d$。

另外，一个必须写清楚的下界：因为每个真实相交对至少要被“产出并计数”一次，计数阶段有
$$
T_1=\Omega(W).
$$


------

## 5.6 版本 I：Sampling 的复杂度

### 时间

- **建索引**：Tsunami 的 build 包含离线优化与重排，通常记为
  $$
  T_{\text{build}}=T_{\text{opt}}+T_{\text{reorder}}.
  $$
  

- **Pass1（计数）**：
  $$
  T_1\approx \sum_{q\in A}\big(w_0 C_q + w_1 S_q(2d)\big)=\Omega(W).
  $$

- **排序 ranks**：$O(t\log t)$ 

- **Pass2（抽样输出）**：只对命中的 $A'\subseteq A$ 再执行 query，最坏仍可到 $\Omega(W)$，但 $t$ 小时常显著少于全体 query。

综合（按你们 baseline 写法）：
$$
T_{\text{Sampling}}
=
T_{\text{build}} + T_1 + O(t\log t) + T_2,
\quad T_2\text{ 最坏 }=\Omega(W).
$$


### 空间（不含 Tsunami 索引本体）

- `deg[q]`：$O(n_1)$
- `Ranks` 与 `Ans`：$O(t)$

因此：
$$
S_{\text{Sampling,extra}}=O(n_1+t).
$$


------

## 5.7 版本 II：Enumerate+Sampling 的复杂度

### 时间

- 枚举阶段：必须产出全部 $W$ 个 pair（至少 $\Omega(W)$），外加 Tsunami query 扫描成本（用上面的代价模型解释即可）。
- 采样阶段：$O(t)$。

### 空间

- `AllPairs`：$\Theta(W)$（主风险点）
- 输出：$O(t)$

所以：
$$
S_{\text{Enumerate}}=\Theta(W)+O(t).
$$

> 这正是该版本的 baseline 意义：小 $W$ 时常数最好；大 $W$ 时会内存爆。

------

## 5.8 版本 III：Adaptive+Sampling 的复杂度

分两种情况：

### Case A：未触发切换（$W\le J_\star$）

- 时间：一遍枚举 + 数组采样
- 空间：$\Theta(W)\le\Theta(J_\star)$

### Case B：触发切换（$W>J_\star$）

- Phase1：无论如何都要完成计数；额外写入/尝试存储最多 $J_\star$ 个 pair（随后丢弃）
- Phase2+Pass2：等价 Sampling 的 ranks + 第二遍输出

因此可写为：
$$
T_{\text{Adaptive}}
=
T_{\text{build}}
+
T_1
+
O(J_\star)
+
O(t\log t)
+
T_2,
$$
其中 $T_1=\Omega(W)$，$T_2$ 最坏也可 $\Omega(W)$。

峰值额外空间（不含索引）：
$$
S_{\text{Adaptive,peak}}
=
O(n_1+t+\min(W,J_\star)).
$$

------

# 结尾：你们在论文里怎么“摆放”这个 baseline（建议一句话版）

> **Tsunami‑JSS Baseline**：把盒子相交 join 转成 $2d$ 维点的正交范围过滤；用 Tsunami’20（Grid Tree + Augmented Grid 的 learned multi‑dim index）作为 range filter 引擎，对每个查询盒子返回所有相交对象并形成确定性 join stream；再通过（1）两遍 rank 定位或（2）显式枚举数组采样或（3）阈值自适应切换，实现 **exact i.i.d. uniform sampling with replacement**，用于对比你们“在最坏情况下不依赖 $|J|$”的采样方法。