# TS19-PBSM Baseline 系列算法设计报告（升级版）

（Sampling / Enumerate+Sampling / Adaptive+Sampling）

------

## 1. 问题定义与分析

### 1.1 输入与对象表示

给定两张空间对象表（或集合）：

- $R=\{r_1,\dots,r_{n_R}\}$
- $S=\{s_1,\dots,s_{n_S}\}$

每个对象用二维轴对齐矩形（MBR）表示。我们统一采用**半开区间**语义来消除“贴边是否算相交”的歧义：
$$
r=[x_l(r),x_u(r))\times[y_l(r),y_u(r)),\quad x_l<x_u,\ y_l<y_u.
$$

> 说明：论文背景中强调 spatial join 常用 “filter step 用 MBR 过滤 + refinement step 精确几何验证”；Tsitsigkos’19 重点研究 filter step（即 MBR 层面的 join）。

### 1.2 相交连接（Intersection Join）

我们关心的 join 结果集合：
$$
J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq\varnothing\}.
$$
在半开区间语义下，矩形相交当且仅当在两维投影都**严格重叠**：
$$
r\cap s\neq\varnothing 
\iff 
\max(x_l(r),x_l(s))<\min(x_u(r),x_u(s))
\ \wedge\
\max(y_l(r),y_l(s))<\min(y_u(r),y_u(s)).
$$

### 1.3 采样目标：在 $J$ 上精确 i.i.d. 均匀有放回

输出 $t$ 个样本对：
$$
Z_1,\dots,Z_t\in J
$$
要求 **精确 i.i.d.、精确均匀、有放回**：
$$
\Pr(Z_j=p)=\frac1{|J|}\ (\forall p\in J),
\quad
Z_1,\dots,Z_t\ \text{相互独立}.
$$

> 核心挑战：$W:=|J|$ 可能非常大（最坏 $\Theta(n_Rn_S)$），因此我们需要多个 baseline 版本权衡时间/空间。

### 1.4 三个版本的定位与取舍

设 $W=|J|$。

- **Enumerate+Sampling**：先枚举并存下全部 $J$，再在数组上均匀采样。
  - 优点：一遍 join，逻辑最简单；
  - 缺点：空间 $O(W)$，当 $W$ 大会爆内存。
- **Sampling（两遍枚举 + 索引定位）**：Pass1 精确计数 $W$，Pass2 用随机 rank 定位输出样本。
  - 优点：空间不依赖 $W$（仅需索引数组与输出）；
  - 缺点：需要二次 join（时间显式依赖 $W$ 且约翻倍）。
- **Adaptive+Sampling**：当 $W$ 小则走 Enumerate+Sampling（快），当 $W$ 大则切到两遍 Sampling（省内存）。
  - 优点：避免 OOM 且尽量减少二次扫描；
  - 缺点：需要阈值 $J_\star$ 的策略。

------

## 2. 引用文章内容与我们继承的关键点

### 2.1 文献来源（必须写在报告里）

**论文**：Dimitrios Tsitsigkos, Panagiotis Bouros, Nikos Mamoulis, Manolis Terrovitis,
 *Parallel In-Memory Evaluation of Spatial Joins*（扩展版；文中注明是 ACM SIGSPATIAL’19 同题扩展）。

该文聚焦：在现代多核、内存充足的背景下，重新设计与调优经典 partition-based spatial join（PBSM）框架，并给出 1D/2D 分区、去重、sweep 轴选择、经验调参等系统性结论。

------

### 2.2 我们 baseline 中“枚举 join 流”的来源点（TS19 → Baseline）

Baseline 的 join 枚举器（EnumerateUniquePairs）在逻辑上继承/复用 Tsitsigkos’19 的以下构件：

#### (A) PBSM（Partition-Based Spatial-Merge Join）= MASJ 框架

- 将空间划分为规则网格 tiles（2D）或 stripes（1D）。
- 每个对象被**多重分配（multi-assignment）**到所有与其相交的分区。
- 每个分区形成独立 join 任务，可并行处理。

> Tsitsigkos’19 在引言/背景部分明确指出 PBSM 是常用 MASJ 框架：先分区，再对每个 tile/stripe 做局部 join。

#### (B) 分区内 join：Forward-Scan Plane Sweep（Algorithm 1）

论文给出一种经典、易实现的 plane sweep 变体（forward scans）：
 先按某一维 lower endpoint 排序（例如 $x_l$），再扫描并对潜在 overlap 的候选对做另一维验证。

> 论文 Algorithm 1 中的关键形式是：对某矩形 $r$，向前扫描另一表中所有满足 $r.x_u \ge s'.x_l$ 的对象作为候选，然后检查 $y$ 是否相交。
>  我们 baseline 使用同等语义，但为了与**半开区间**一致，把“$\ge$”改为严格“$>$”（见 3.4 与 4.2 伪代码），避免“贴边被当成相交候选”。

#### (C) Duplicate elimination：reference point 去重（Eq. (1)）

PBSM 多重分配会使同一对 $(r,s)$ 在多个分区被发现。论文采用 Dittrich & Seeger 的 reference-point 思想：
 只在**某个参考点落入该分区**时报告输出，从而避免重复。论文给出 2D tiles 的判定式（Eq. (1)）以及 1D stripes 简化。

> 论文文字称“top-left corner”，但 Eq.(1) 采用 $\max(\cdot)$ 形式对应的是 intersection 区域的 **(max lower endpoints)** 角点。为避免坐标系歧义，本报告直接把 reference point 定义为 $(\max x_l,\max y_l)$（见 3.5），并用 Eq.(1) 表达式实现。

#### (D) sweep axis 选择模型：直方图估计候选数（Eq. (2)）

论文提出：对每个分区 $T$，估计若沿 x 扫描会产生多少“x 投影相交候选对”：
$$
I_T^x \approx \sum_{i=0}^{k-1} H_{R,T}^x[i]\cdot H_{S,T}^x[i]
$$
并取 $\min(I_T^x, I_T^y)$ 对应轴作为 sweep axis；同时给出直方图 bucket 数 $k$ 与采样参数 $\phi=100$ 的工程建议。

#### (E) 分区数 $K$ 的经验规则 + 1D 优于 2D 的结论

论文实验发现：

- 增加分区数会提升性能，但过大时复制（replication）成本会反噬；
- “经验规则”：令分区的尺度约为对象平均尺度的一个数量级（约 10 倍）通常表现最好；
- 1D stripes 往往比 2D grid 更快（更少复制等原因）。

------

### 2.3 我们 baseline 中“Sampling 三版本”属于在 PBSM 枚举器之上的严格抽样封装

Tsitsigkos’19 关注“高性能枚举 join”，并不直接讨论“精确均匀抽样”。因此：

- 我们 baseline 的核心思想是：把 TS19-PBSM 作为**唯一枚举流生成器**；
- 再用“数组抽样 / 两遍 rank 抽样 / 自适应切换”实现**精确 i.i.d. uniform（with replacement）**。

------

## 3. 核心数据结构

本章把 Baseline 的组件拆为两层：

1. **PBSM 唯一枚举器（Enumerator）**
2. **抽样器（Sampler）**：三种版本共享的抽样相关结构

### 3.1 基础对象结构

对每个矩形 $r$ 保存：

- `id`：唯一整数（建议 32-bit 或 64-bit，取决于规模）
- `xl, xu, yl, yu`：双精度或定点坐标（建议 double）
- 可选：原对象指针/索引（用于输出时回指）

**统一相交判定（半开）**：

```
Intersect(r,s) :=
    max(r.xl, s.xl) < min(r.xu, s.xu)
and max(r.yl, s.yl) < min(r.yu, s.yu)
```

> 这一定义要贯穿：分区覆盖判定、forward-scan 的 while 条件、另一维过滤、去重 reference point 逻辑一致。

------

### 3.2 分区结构（PBSM）

我们支持两种分区方案：

#### 3.2.1 1D stripes（推荐默认）

以 x 轴为例（垂直条带）：

- 先确定全局 join 空间范围（domain）：
  $$
  X_{\min}=\min_{r\in R\cup S} x_l(r),\quad
  X_{\max}=\max_{r\in R\cup S} x_u(r).
  $$

- 等分为 $K$ 个半开条带：
  $$
  T_j=[X_j,X_{j+1})\times[Y_{\min},Y_{\max}),\quad j=0,\dots,K-1.
  $$

每个条带维护两个列表：

- `R_T[j]`：所有与条带 **x 范围严格重叠** 的 $r\in R$（multi-assignment）
- `S_T[j]`：同理

> 这里的 “$\times(-\infty,\infty)$” 在理论上可写，但实现中应使用整体数据 y 范围 $[Y_{\min},Y_{\max})$ 或直接忽略 y（因为 stripes 不按 y 分割）。

#### 3.2.2 2D grid（可选）

- 将 $[X_{\min},X_{\max})\times[Y_{\min},Y_{\max})$ 分为 $K\times K$ tiles：
  $$
  T_{i,j}=[X_i,X_{i+1})\times[Y_j,Y_{j+1})
  $$

- 每 tile 存 `R_T[i][j]`, `S_T[i][j]`

------

### 3.3 分区内排序结构（为 forward-scan 服务）

对每个分区 $T$：

- 选择 sweep axis：`axis[T] ∈ {x, y}`
- 分别对 `R_T[T]` 与 `S_T[T]` 排序：

若 `axis[T]=x`：按 $(x_l,\ id)$ 升序
 若 `axis[T]=y`：按 $(y_l,\ id)$ 升序

> 强烈建议：**直接用 (key, id) 作为比较器形成全序**，这样即使底层 sort 不稳定也能保证确定性。

------

### 3.4 Forward-Scan 所需的访问函数

对某轴 `axis ∈ {x,y}`：

- `lower(axis, r)`：若 x 轴则 $x_l(r)$，若 y 轴则 $y_l(r)$
- `upper(axis, r)`：若 x 轴则 $x_u(r)$，若 y 轴则 $y_u(r)$
- `IntersectOtherAxis(r,s,axis)`：若 axis=x 则检查 y 轴相交（严格），反之检查 x

------

### 3.5 Duplicate elimination（reference point）结构与判定

为每对相交矩形定义 reference point：
$$
p(r,s) := \bigl(\max(x_l(r),x_l(s)),\ \max(y_l(r),y_l(s))\bigr).
$$

- **2D tile 去重（来自论文 Eq.(1)）**：对 tile $T$，当且仅当
  $$
  p_x(r,s)\ge T.x_l \ \wedge\ p_y(r,s)\ge T.y_l
  $$
  才允许该 tile “报告/计数”该 pair。

- **1D vertical stripes 简化**：对 stripe $T=[T.x_l,T.x_u)\times\mathbb R$，当且仅当
  $$
  p_x(r,s)\ge T.x_l
  $$
  才报告。

> 为什么只检查“$\ge T.x_l$”而不检查 “$”？
>  直观原因：若 pair 在分区 $T$ 中被发现（说明两矩形都被分配到 $T$），那么必然有 $p_x（否则至少有一个矩形不会与该分区重叠），因此上界由“能被发现”隐含保证。该判定正是 Dittrich & Seeger 去重的经典形式，论文给出同样的简化说明。

------

### 3.6 抽样相关结构（与三版本对应）

- `Pairs[]`：仅 Enumerate+Sampling / Adaptive 小分支使用。每个元素为 `(rid, sid)`
  - 若 id 用 32-bit，则理论最小 8B；但考虑容器开销建议按 16B 估算阈值更稳。
- `W`：`uint64` 计数器，记录 $|J|$
- `I[1..t]`：随机索引数组（`uint64`），每个均匀落在 $[0,W-1]$
- `IdxPos[]`：把 `(I_j, j)` 按 `I_j` 升序排序，Pass2 线性匹配输出
- `Ans[1..t]`：输出样本数组

> 随机数必须无偏生成：不要直接 `rand64 % W`（可能引入偏差）。建议 rejection sampling（见 4.5 工程细节）。

------

## 4. 算法详细流程（三个版本）

我们先定义一个公共的 **PBSM 唯一枚举器**，再在其上实现三个版本。

------

### 4.1 公共预处理：BuildPartitionsAndSort

**输入**：矩形集合 $R,S$，分区方案（1D/2D），分区数 $K$，是否启用 sweep axis 模型
 **输出**：分区列表 `R_T, S_T`（已分配并排序）、每分区 `axis[T]`

#### Step A：确定 join domain（全局边界）

$$
X_{\min}=\min x_l,\ X_{\max}=\max x_u,\ 
Y_{\min}=\min y_l,\ Y_{\max}=\max y_u.
$$

#### Step B：multi-assignment 分区分配

以 1D vertical stripes 为例，条带宽度：
$$
w = \frac{X_{\max}-X_{\min}}{K}.
$$
一个矩形 $r$ 被分配到条带 $T_j$ 的条件（严格重叠）为：
$$
[x_l(r),x_u(r)) \cap [X_j,X_{j+1}) \neq \varnothing
\iff x_l(r)<X_{j+1}\ \wedge\ x_u(r)>X_j.
$$
实现中通常通过计算覆盖的条带索引范围 $[j_{\min},j_{\max}]$ 并循环追加：

- $j_{\min}=\left\lfloor \frac{x_l(r)-X_{\min}}{w}\right\rfloor$（clamp 到 $[0,K-1]$）
- $j_{\max}=\left\lceil \frac{x_u(r)-X_{\min}}{w}\right\rceil-1$（clamp 到 $[0,K-1]$）

然后对所有 $j\in[j_{\min},j_{\max}]$ 执行 `R_T[j].push_back(r.id)`（或 `S_T`）。

> **浮点边界建议**：若担心 $x_u$ 正好落在分区边界造成误判，可用 `nextafter(xu, -∞)` 计算“左邻值”来实现半开语义的一致性。

2D grid 情况同理，分别算 x/y 维覆盖范围并做双循环。

#### Step C：为每分区选择 sweep axis（可选但建议）

**论文模型（Eq.(2)）**：对分区 $T$，构造直方图估计“沿 x 扫描的候选对数量”：
$$
I_T^x = \sum_{i=0}^{k-1} H_{R,T}^x[i]\cdot H_{S,T}^x[i],
\quad
I_T^y \text{ 同理}.
$$
取 `axis[T] = argmin(I_T^x, I_T^y)`。论文还给出：

- 大 tile 取 $k=1000$；小 tile 则按 tile 尺度与平均对象尺度比例设置；
- 直方图用采样构建，默认每 $\phi=100$ 个对象取 1 个用于统计，开销很小。

**实现建议（自包含，不依赖论文细节）**：

- 直方图 bins：对 x 轴把 $[T.x_l,T.x_u)$ 等分为 $k$ 个 bins；
- 对每个矩形投影 $[x_l,x_u)$，计算覆盖的 bin 区间 $[b_l,b_u]$，对 difference array 做区间加一，最后 prefix 得到每个 bin 覆盖计数；
- 分别对 $R_T, S_T$ 做上述统计即可得到 $H_{R,T}^x, H_{S,T}^x$。

**1D stripes 的简化经验**：论文实验观察：若分区轴为 x，则 sweep 轴通常应选 y（正交组合效果好），反之亦然；且若分区轴与 sweep 轴同轴（xx/yy），随 $K$ 增大 join 成本不明显下降。
 因此：若不做直方图模型，1D stripes 可以默认选正交轴。

#### Step D：排序

按 `axis[T]` 对 `R_T[T]`、`S_T[T]` 排序，key 为 `(lower(axis), id)`。

------

### 4.2 公共子过程：EnumerateUniquePairs（PBSM 唯一枚举器）

这个枚举器做三件事：

1. 遍历每个分区 $T$
2. 在分区内用 forward-scan plane sweep 找出相交对
3. 用 DuplicateTest 保证每个全局 pair 只输出一次

> 分区遍历顺序必须固定（1D：`j=0..K-1`；2D：row-major）以便复现实验。

#### Forward-Scan（半开语义版）伪代码

下面是对论文 Algorithm 1 的工程化整合写法：核心逻辑对应 forward scan；差别仅在把 `>=` 改为 `>` 以匹配半开“严格重叠”。

```
PROC EnumerateUniquePairs(Partitions, R, S, axis, callback):
    for each partition T in deterministic order:
        Rt = R_T[T]  // sorted by (lower(axis[T]), id)
        St = S_T[T]  // sorted similarly
        a  = axis[T]

        i = 0; j = 0
        while i < |Rt| and j < |St|:
            r = Rt[i]; s = St[j]

            if lower(a,r) < lower(a,s):
                j2 = j
                // half-open: need strict overlap => upper(a,r) > lower(a,St[j2])
                while j2 < |St| and upper(a,r) > lower(a,St[j2]):
                    s2 = St[j2]
                    if IntersectOtherAxis(r, s2, a):
                        if DuplicateTest(T, r, s2):
                            callback(r.id, s2.id)
                    j2++
                i++
            else:
                i2 = i
                while i2 < |Rt| and upper(a,s) > lower(a,Rt[i2]):
                    r2 = Rt[i2]
                    if IntersectOtherAxis(r2, s, a):
                        if DuplicateTest(T, r2, s):
                            callback(r2.id, s.id)
                    i2++
                j++
```

**IntersectOtherAxis（严格）**：

- 若 `a=x`：检查 $\max(y_l(r),y_l(s)) < \min(y_u(r),y_u(s))$
- 若 `a=y`：检查 $\max(x_l(r),x_l(s)) < \min(x_u(r),x_u(s))$

**DuplicateTest**：

- 2D tiles：$\max(x_l)\ge T.x_l \wedge \max(y_l)\ge T.y_l$（Eq.(1)）
- 1D vertical stripes：$\max(x_l)\ge T.x_l$（论文简化）

------

### 4.3 版本一：Enumerate+Sampling（存全量后抽样）

**适用**：$W=|J|$ 不大，内存足够。

#### 流程

1. 预处理：`BuildPartitionsAndSort(...)`
2. 枚举并存储所有 pair：

```
Pairs = []
EnumerateUniquePairs(..., callback = (rid,sid) => Pairs.push((rid,sid)))
W = |Pairs|
```

1. i.i.d. 有放回均匀采样：

```
if W == 0: return empty
for j in 1..t:
    idx ~ Unif{0..W-1}   // unbiased
    Ans[j] = Pairs[idx]
return Ans
```

#### 关键性质

- `Pairs` 是 $J$ 的一一枚举（由去重机制保证），故数组均匀取下标 ⇒ 精确均匀；独立取下标 ⇒ 精确独立。

------

### 4.4 版本二：Sampling（两遍枚举 + 随机 rank 定位）

**适用**：$W$ 可能很大，不允许存下全部 $J$，但仍需精确 i.i.d. uniform。

#### 流程概览

- Pass1：只计数 $W$
- 生成 $t$ 个独立均匀 rank：$I_1,\dots,I_t \sim \text{Unif}\{0,\dots,W-1\}$
- Pass2：再次枚举同一个“唯一 join 流”，当枚举计数器走到这些 rank 时输出对应 pair

#### 详细步骤

**Step 0：预处理一次并复用（推荐）**

- 预处理 `BuildPartitionsAndSort` 一次，保存所有分区列表与排序结果。
- Pass1/Pass2 直接复用，避免重复排序。

> **理论上**：正确性只要求 Pass2 是对 $J$ 的一一枚举；Pass1 仅用于得到 $W$。
>  **工程上**：为了复现实验与 debug，仍建议使用确定性的分区顺序与排序 key。

**Pass1：计数**

```
W = 0
EnumerateUniquePairs(..., callback = (rid,sid) => W++)
```

**生成随机索引（有放回）**
$$
I_1,\dots,I_t \stackrel{i.i.d.}{\sim} \mathrm{Unif}\{0,\dots,W-1\}.
$$
工程实现：

```
IdxPos = [(I_j, j) for j=1..t]
sort IdxPos by I_j ascending
```

**Pass2：定位输出**

```
Ans[1..t]
p = 1     // pointer in IdxPos
c = 0     // current rank in the unique join stream, 0..W-1

EnumerateUniquePairs(..., callback = (rid,sid) => {
    while p <= t and IdxPos[p].index == c:
        Ans[IdxPos[p].pos] = (rid, sid)
        p++
    c++
})
return Ans
```

- 若多个索引相同：会在同一 $c$ 上输出同一个 pair 多次 ⇒ **有放回**

------

### 4.5 版本三：Adaptive+Sampling（自动切换）

**目标**：

- 小 $W$：一遍枚举存全量 + 数组采样（快）
- 大 $W$：两遍 rank sampling（不爆内存）

引入阈值 $J_\star$。

#### Adaptive Phase1：一次扫描（始终计数，枚举是可选附加）

```
mode = ENUMERATE
Pairs = []        // capacity J_star
W = 0

EnumerateUniquePairs(..., callback = (rid,sid) => {
    W++
    if mode == ENUMERATE:
        if W <= J_star:
            Pairs.push((rid,sid))
        else:
            mode = COUNT_ONLY
            Pairs = empty_and_free()  // discard partial enumeration
})
```

结束后：

- 若 `mode==ENUMERATE` ⇒ $W\le J_\star$ 且 `Pairs` 恰好是全体 $J$
- 若 `mode==COUNT_ONLY` ⇒ $W>J_\star$ 且已丢弃枚举结果（只保留计数）

#### 分支 A：小 join（无需 Pass2）

```
if W == 0: return empty
for j in 1..t:
    idx ~ Unif{0..W-1}
    Ans[j] = Pairs[idx]
return Ans
```

#### 分支 B：大 join（执行两遍 Sampling 的 Pass2）

- 生成并排序 `IdxPos`（同 4.4）
- 再跑一次 `EnumerateUniquePairs` 按 rank 定位输出（同 4.4 Pass2）

------

## 5. Adaptive 阈值 $J_\star$ 的选择策略

目标：**不 OOM** + **尽量避免不必要的第二遍**。

### 5.1 内存硬约束（必须满足）

设可用内存预算为 `MemBudget`（bytes），允许给 `Pairs` 使用比例 $\rho\in(0,1)$。
 每个 pair 存储开销 `sizeof(pair)`：

- 理论最小：两个 32-bit id ⇒ 8B
- 实际更稳：考虑对齐/容器开销 ⇒ 建议按 16B 估算

则
$$
J_\star^{\text{mem}}
=\left\lfloor\frac{\rho\cdot \text{MemBudget}}{\text{sizeof(pair)}}\right\rfloor,
\quad
J_\star \le J_\star^{\text{mem}}.
$$

### 5.2 时间权衡（交叉点思想）

粗略对比：

- Enumerate+Sampling：
  $$
  T_{\text{enum}} \approx T_{\text{join}} + \text{WriteCost}\cdot W + O(t)
  $$

- 两遍 Sampling：
  $$
  T_{\text{2pass}} \approx 2\cdot T_{\text{join}} + O(t\log t)
  $$

交叉点 $W_{\text{cross}}$ 可用 micro-benchmark 拟合：
 当 $W \approx W_{\text{cross}}$ 时两者耗时相当，则
$$
J_\star^{\text{time}} \approx W_{\text{cross}}.
$$
最终建议：
$$
J_\star=\min\left(J_\star^{\text{mem}},\ J_\star^{\text{time}}\right).
$$

### 5.3 与 TS19 的参数调优协同（选择 $K$ 与 axis）

$T_{\text{join}}$ 强烈依赖 PBSM 参数（分区数、分区类型、sweep 轴选择）。论文给出强经验结论：

- **分区尺度约为对象平均尺度的一个数量级（~10倍）**常常最好；过细会因复制过度变慢。
- 1D stripes 往往优于 2D grid；
- 1D 情况下 partition 轴与 sweep 轴正交（xy/yx）随 $K$ 增大效果更好，而同轴（xx/yy）不明显改善。

因此推荐流程：

1. 先用 TS19 经验规则选 $K$ 与 axis（或直方图模型）做默认配置；
2. 在该配置下测得 $W_{\text{cross}}$，再决定 $J_\star$。

------

## 6. 算法分析（正确性 + 复杂度）

本章覆盖：枚举器正确性、三版本采样正确性、时间/空间复杂度。

------

### 6.1 枚举器正确性（EnumerateUniquePairs 输出的是 $J$ 的一一枚举）

#### 引理 1（覆盖性）：任意 $(r,s)\in J$ 至少会在某个分区中被“发现”

由于 PBSM multi-assignment：每个矩形被分配到所有与其空间范围相交的 tiles/stripes。对任意相交对 $(r,s)$，它们的交集区域 $r\cap s$ 非空，取任一点 $p\in r\cap s$。必存在唯一分区 $T$ 使得 $p\in T$（分区是空间的半开剖分）。
 因为 $p\in r$ 且 $p\in s$，故 $r$ 与 $T$ 相交、$s$ 与 $T$ 相交，于是二者都会被分配到 $T$，从而在 $T$ 的局部 join 中成为候选并被检测。

> 这正是 PBSM/MASJ “对象复制到所有重叠分区”确保不漏的基本原理。

#### 引理 2（分区内完整性）：Forward-Scan 能枚举出分区内所有相交对（未去重前）

在固定分区 $T$ 内，设 $R_T,S_T$ 已按 sweep 轴 lower endpoint 排序。forward-scan 的逻辑等价于论文 Algorithm 1：
 每当扫描遇到某个对象的 lower endpoint，就向前扫描另一列表中所有 lower 小于自身 upper 的对象作为候选，并在另一维做严格相交验证。该过程保证不会漏掉任何在 sweep 轴上可能重叠的配对，从而不漏掉真实相交对。

> 我们仅把候选条件由 “$\ge$” 改为 “$>$” 以匹配半开区间定义，这不会影响“覆盖所有真正相交”的性质，只会去掉“贴边候选”。

#### 引理 3（全局唯一性）：DuplicateTest 保证每个全局 pair 恰好被报告一次

对任意 $(r,s)\in J$，定义 reference point：
$$
p(r,s)=(\max x_l,\max y_l).
$$

- 由于分区是半开剖分，存在唯一分区 $T^\star$ 使得 $p(r,s)\in T^\star$。
- 论文 Eq.(1) 的 DuplicateTest 恰好刻画了 “$p(r,s)$ 不在分区左/下方” 的条件，结合 “pair 在该分区被发现” 的前提，可推出：
  - 若 $T=T^\star$，则 DuplicateTest 通过；
  - 若 $T\neq T^\star$，则至少在 x 或 y 上落在 $T$ 的左/下方，从而 DuplicateTest 失败。

因此每个 pair 恰好在唯一分区被输出一次。

#### 推论（枚举流是 $J$ 的一一对应）

由引理 1（不漏）+ 引理 3（不重且仍不漏）得：
 `EnumerateUniquePairs` 输出流与 $J$ 存在一一对应。

------

### 6.2 三版本采样正确性

#### 定理 1（Enumerate+Sampling 正确）

`Pairs` 是 $J$ 的一一枚举且长度 $W=|J|$。
 每个样本位置独立生成 `idx ~ Unif{0..W-1}` 并输出 `Pairs[idx]`，显然得到：

- $\Pr(Z_j=p)=1/|J|$
- 各 $Z_j$ 相互独立

#### 定理 2（两遍 Sampling 正确）

Pass1 得到精确 $W=|J|$。
 Pass2 枚举流是 $J$ 的一一枚举，给每个 pair 定义其在 Pass2 流中的 rank：$\mathrm{rank}(p)\in[0,W-1]$。
 因 $I_j\sim \mathrm{Unif}\{0,\dots,W-1\}$：
$$
\Pr(Z_j=p)=\Pr(I_j=\mathrm{rank}(p))=\frac{1}{W}=\frac{1}{|J|}.
$$
又 $I_j$ 独立 ⇒ $Z_j$ 独立。
 重复索引自然对应“有放回”。

> 注意：**正确性不要求 Pass1 与 Pass2 的枚举顺序一致**；只要求 Pass2 是一一枚举。
>  但为了可复现实验，建议仍固定分区遍历顺序与排序比较器，使 Pass2 顺序确定。

#### 定理 3（Adaptive 正确）

Adaptive 仅在两条均正确的路径间选择：

- 若 $W\le J_\star$：退化为 Enumerate+Sampling（定理 1）
- 若 $W>J_\star$：退化为两遍 Sampling（定理 2）

阈值仅影响性能与内存，不影响分布。

------

### 6.3 复杂度分析（统一符号）

记：

- $n_R=|R|,\ n_S=|S|,\ n=n_R+n_S$

- 分区数：1D stripes 为 $K$，2D grid 为 $K^2$

- multi-assignment 后总列表长度：
  $$
  n'=\sum_T |R_T|+\sum_T|S_T|
  $$

- 排序开销：
  $$
  \text{SortCost}=\sum_T\bigl(|R_T|\log|R_T|+|S_T|\log|S_T|\bigr)
  $$

- forward-scan 候选检查/扫描成本（数据相关）记为：
  $$
  \text{SweepCost}=\sum_T \mathrm{JoinCost}(T)
  $$
  （其中 $\mathrm{JoinCost}(T)$ 可理解为候选对数量级 + 指针扫描开销）

此外，DuplicateTest 是每个候选对 $O(1)$ 的常数开销（Eq.(1) 或 1D 简化）。

------

#### 6.3.1 预处理成本（所有版本共享）

- 分区分配：$O(n')$
- sweep axis 决策：论文指出模型开销远小于 join，可视为轻量统计；严格计入也近似为每分区线性于样本量 + $k$ 的成本。
- 排序：$\text{SortCost}$

------

#### 6.3.2 Enumerate+Sampling

- 时间：
  $$
  T_{\text{Enum+Samp}}=
  O\bigl(n' + \text{SortCost} + \text{SweepCost} + |J| + t\bigr)
  $$
  （$|J|$ 来自写入 `Pairs`）

- 空间：
  $$
  S_{\text{Enum+Samp}}=O(n'+|J|+t)
  $$

------

#### 6.3.3 Sampling（两遍）

若复用已排序的分区列表，则 Pass2 不再排序：

- 时间：
  $$
  T_{\text{Samp}}=
  O\bigl(n' + \text{SortCost} + 2\cdot\text{SweepCost} + 2|J| + t\log t\bigr)
  $$
  （两遍枚举都会“访问到”每个真实输出对一次，因此有 $2|J|$ 计数/推进）

- 空间：
  $$
  S_{\text{Samp}}=O(n'+t)
  $$

------

#### 6.3.4 Adaptive+Sampling

- 若 $W\le J_\star$：同 Enumerate+Sampling
- 若 $W>J_\star$：Phase1 最多额外存 $J_\star$ 个再丢弃，然后还要做 Pass2：

时间：
$$
T_{\text{Adaptive}}=
O\bigl(n'+\text{SortCost}+2\cdot\text{SweepCost}+2|J|+t\log t+J_\star\bigr)
$$
空间：
$$
S_{\text{Adaptive}}=O(n'+J_\star+t).
$$

------

## 附录：工程实现关键细节（强烈建议写进报告最后）

这些点不写，baseline 很容易“看似正确但实现踩坑”，尤其是你们强调“精确 i.i.d.”时，审稿人/复现者会抓。

### A1. 生成 $[0,W-1]$ 的无偏均匀整数

不要直接 `u % W`。推荐 rejection sampling：

- 取 64-bit 均匀随机数 $U\in[0,2^{64}-1]$
- 令 $M=\left\lfloor\frac{2^{64}}{W}\right\rfloor\cdot W$
- 若 $U\ge M$ 则丢弃重抽；否则输出 $U\bmod W$

这样保证严格均匀。

### A2. 半开语义的一致性（必须全链路一致）

- 相交判定用严格 `<`
- forward-scan 候选 while 用严格 `>`
- 分区覆盖条件用严格重叠（不要把贴边当 overlap）
- 这样才能保证定义、实现、报告一致

### A3. 确定性顺序（为了可复现，不是 correctness 必需）

- 分区遍历顺序固定（stripe id 升序 / tile row-major）
- 排序比较器用 `(lower, id)` 形成全序（不依赖 sort 稳定性）
- forward-scan 在相等 lower 的分支固定走 `else`（或明确 tie-break）

### A4. 为什么推荐 1D stripes 作为 baseline 默认

论文实验明确指出：

- 1D partitioning 往往比 2D 更快；
- 正交的 partition/sweep 组合（xy/yx）随 $K$ 增大更有效；
- 经验规则：分区尺度约为对象尺度的 ~10 倍。

这些都非常适合 baseline：实现简单、常数小、性能稳定。