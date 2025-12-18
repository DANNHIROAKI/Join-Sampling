# TS19-PBSM Baseline 系列算法文档

（Sampling / Enumerate+Sampling / Adaptive+Sampling）

## 1. 问题定义与分析 + 引用文章内容

### 1.1 问题定义（与我们对齐）

给定两张空间对象表（或集合）：

- $R=\{r_1,\dots,r_{n_R}\}$
- $S=\{s_1,\dots,s_{n_S}\}$

每个对象用 **2D 轴对齐矩形/MBR**（半开区间）表示：
$$
r=[x_l(r),x_u(r))\times[y_l(r),y_u(r)), \quad x_l<x_u,\ y_l<y_u.
$$
我们关心的空间连接是 **相交连接（intersection join）**：
$$
J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq\varnothing\}.
$$
矩形相交条件（半开区间）：
$$
r\cap s\neq\varnothing 
\iff 
\max(x_l(r),x_l(s))<\min(x_u(r),x_u(s))
\ \wedge\
\max(y_l(r),y_l(s))<\min(y_u(r),y_u(s)).
$$
Tsitsigkos’19 同样以“相交 join”为核心，并指出实际系统常采用 MBR 做 filter step（本文也可视作 filter step 的 join），与我们的任务对齐。 [arXiv](https://arxiv.org/pdf/1908.11740)

------

### 1.2 我们的采样目标（对齐你们的 Sampling 目标）

输出 $t$ 个样本对：
$$
Z_1,\dots,Z_t\in J
$$
要求 **i.i.d. 均匀有放回采样**：
$$
\Pr(Z_j=p)=\frac1{|J|}\ (\forall p\in J),
\quad
Z_1,\dots,Z_t\ \text{相互独立}.
$$

> 下面三种版本都会做到“**精确** uniform + **精确** independent（不是近似）”，区别只在时间/空间开销与是否需要二次扫描。

------

### 1.3 Baseline 的文献来源与可继承关键点（Tsitsigkos’19 → 我们的 Baseline）

**论文信息：**
 Tsitsigkos, Bouros, Mamoulis, Terrovitis, *Parallel In-Memory Evaluation of Spatial Joins*（arXiv 扩展版，亦说明是 SIGSPATIAL’19 同题扩展）。 [arXiv+1](https://arxiv.org/pdf/1908.11740)

我们继承/复用的“关键方法”如下（每条都能在论文中找到对应）：

1. **PBSM（Partition-Based Spatial-Merge Join）作为 MASJ 框架**
    把空间分区成 tiles 或 stripes；每个分区形成独立 join 任务，易并行；且可通过机制避免重复结果。 [arXiv+1](https://arxiv.org/pdf/1908.11740)
2. **分区内 join：Forward-Scan Plane Sweep（论文 Algorithm 1）**
    先按某轴 lower endpoint 排序，再 forward scan 产生候选对，并检查另一维是否相交。该算法在论文中作为“处理小 join 的典型 plane sweep”。 [arXiv](https://arxiv.org/pdf/1908.11740)
3. **Duplicate avoidance：Reference-point / 参考点去重（论文 Eq.(1)）**
    PBSM 多重分配会导致同一对 $(r,s)$ 在多个分区被发现。Tsitsigkos’19 采用 Dittrich & Seeger 的思想：只在**参考点**落入的分区报告该 pair，并给出 2D tile 的判定式(1)。 [arXiv](https://arxiv.org/pdf/1908.11740)
4. **1D stripes 的简化去重**
    当采用 1D 垂直条带时，duplicate test 只需单维比较（更易实现、常数更小）。 [arXiv](https://arxiv.org/pdf/1908.11740)
5. **分区内 sweeping axis 的选择模型（论文 Eq.(2) + 直方图估计）**
    论文用直方图估计在 x 或 y 上的“投影相交候选数”，选择候选更少的轴作为 sweep axis；并给出采样构建直方图的工程参数（如默认每 $\phi=100$ 个对象取 1 个参与统计）。 [arXiv+1](https://arxiv.org/pdf/1908.11740)
6. **分区数 $K$ 的经验规则（rule of thumb）**
    论文指出增加分区数会提升性能但到某个点后复制开销反而变差，并给出经验：选择 $K$ 使得分区尺度约比矩形尺度“大一个数量级”。 [arXiv](https://arxiv.org/pdf/1908.11740)

> 综上：我们的 Baseline 在“**枚举 join 流**”这一层完全站在 Tsitsigkos’19 的 PBSM 框架上；而 Sampling 的“抽样实现”是我们为了匹配你们的 i.i.d. uniform 目标而做的严格实现（两遍枚举或数组抽样）。

------

### 1.4 三个版本的定位（为什么需要三个）

设 $W:=|J|$。

- **Enumerate+Sampling（最直接 baseline）**：先枚举并存下全部 $J$，再在数组上均匀抽样。
   优点：一遍 join；采样简单。缺点：空间 $O(W)$。
- **Sampling（两遍枚举 baseline）**：不存下 $J$，Pass1 精确计数 $W$，Pass2 用随机索引定位样本。
   优点：空间不依赖 $W$。缺点：需要二次 join（时间依赖 $W$ 且约翻倍）。
- **Adaptive+Sampling（工程友好 baseline）**：当 $W$ 小则走 Enumerate+Sampling（快），当 $W$ 大则走 Sampling（省内存），自动切换。
   优点：小结果快、大结果不爆内存。缺点：需要阈值策略 $J_\star$。

------

## 2. 核心数据结构

下面的数据结构按“PBSM 枚举器 + 抽样器”拆分。

### 2.1 基础对象与索引

- 每个矩形保存：
  - $x_l,x_u,y_l,y_u$
  - 唯一 id（整数），用于稳定排序与输出 pair（只输出 id 对即可）。
- 统一使用半开区间 $[x_l,x_u)\times[y_l,y_u)$ 来避免边界相切的歧义（与 Tsitsigkos’19 的写法一致）。 [arXiv+1](https://arxiv.org/pdf/1908.11740)

------

### 2.2 空间分区结构（PBSM）

#### 2.2.1 1D stripes（推荐默认 baseline）

选定分区轴（例如沿 x 轴做垂直条带），分成 $K$ 个条带：
$$
T_j=[X_j,X_{j+1})\times(-\infty,\infty),\quad j=0,\dots,K-1.
$$
每个条带维护两个列表：

- `R_T[j]`: 所有与条带 x 范围有交叠的 $r\in R$（多重分配）
- `S_T[j]`: 所有与条带 x 范围有交叠的 $s\in S$

**经验选择：**
 论文对 1D/2D 分区做了系统调参，并给出经验：选 $K$ 使得分区尺度约比对象尺度“大一个数量级”，过细会导致复制过度。 [arXiv](https://arxiv.org/pdf/1908.11740)

#### 2.2.2 2D grid（可选）

将空间分成 $K\times K$ tiles。每 tile $T_{i,j}$ 维护 `R_T[i][j]`、`S_T[i][j]`。
 去重条件要用 2D 版本（Eq.(1)）。 [arXiv+1](https://arxiv.org/pdf/1908.11740)

------

### 2.3 分区内 plane sweep 结构（Forward-Scan）

对每个分区 $T$：

- 选择 sweep axis：x 或 y
- 将 `R_T` 与 `S_T` 按 sweep axis 的 lower endpoint 排序（稳定）：
  - 若 sweep=x：按 $(x_l,\text{id})$ 排
  - 若 sweep=y：按 $(y_l,\text{id})$ 排

> 稳定排序 + id tie-break 是为了后续“两遍枚举”能得到完全相同的枚举顺序（Sampling/Adaptive 分支需要）。

Forward-scan 本身不需要额外复杂结构，只需要排序数组 + 游标指针。算法原型见 Tsitsigkos’19 Algorithm 1。 [arXiv](https://arxiv.org/pdf/1908.11740)

------

### 2.4 Duplicate avoidance 判定（核心：保证“唯一枚举流”）

#### 2.4.1 2D tile 的 duplicate test（论文 Eq.(1)）

若 $(r,s)$ 在 tile $T$ 内被发现相交，则仅当满足：
$$
\max\{x_l(r),x_l(s)\}\ge T.x_l \ \wedge\ \max\{y_l(r),y_l(s)\}\ge T.y_l
$$
才“报告/计数”该 pair。 [arXiv](https://arxiv.org/pdf/1908.11740)

#### 2.4.2 1D 垂直条带的简化

若是垂直条带 $T=[T.x_l,T.x_u)\times\mathbb{R}$，则只需：
$$
\max\{x_l(r),x_l(s)\}\ge T.x_l
$$
即可决定是否由该条带报告结果。 [arXiv](https://arxiv.org/pdf/1908.11740)

------

### 2.5 Sweep axis 选择的直方图结构（可选优化，但建议加上）

对于分区 $T$，Tsitsigkos’19 用直方图估计在 x/y 上投影相交的候选数：

- 把 tile 在 x 或 y 上分成 $k$ 段
- 统计每段内 $R_T$ 与 $S_T$ 的投影覆盖数量形成直方图
- 估计在 x 上的候选对数：

$$
I_T^x \approx \sum_{i=0}^{k-1} H_{R,T}^x[i]\cdot H_{S,T}^x[i]
$$

选 $\min(I_T^x,I_T^y)$ 对应的轴做 sweep axis。 [arXiv+1](https://arxiv.org/pdf/1908.11740)

工程参数：论文也给出了 $k$ 的设置思路以及“用采样构建直方图”的默认 $\phi=100$。 [arXiv+1](https://arxiv.org/pdf/1908.11740)

------

### 2.6 采样相关结构（因版本不同而不同）

- `Pairs[]`：仅 Enumerate+Sampling / Adaptive 的小分支使用；存所有（或最多 $J_\star$ 个）pair 的 (rid, sid)。
- `W`：64-bit 计数器，记录 $|J|$。
- `I[1..t]`：Sampling/Adaptive 大分支用的随机索引（64-bit），每个均匀落在 $[0,W-1]$。
- `IdxPos[]`：把 $(I_j,j)$ 排序（按索引升序）用于 Pass2 线性扫描定位。
- `Ans[1..t]`：输出样本数组。

------

## 3. 算法详细流程（Sampling / Enumerate+Sampling / Adaptive+Sampling）

下面先定义一个**公共的“PBSM 唯一枚举器”接口**，然后在其上实现三个版本。

### 3.1 公共子过程：PBSM 唯一枚举器

#### 3.1.1 预处理：构建分区 + 每分区排序 + 记录 sweep axis

输入：矩形集合 $R,S$，分区方案（1D/2D），分区数 $K$。

输出：每个分区 $T$ 的 `R_T`、`S_T`（已排序）以及 `SweepAxis[T]`。

**步骤 A：分区（multi-assignment）**

- 对每个 $r\in R$：计算其覆盖的条带/tiles，把 $r.id$ 追加到对应 `R_T`。
- 对每个 $s\in S$：同理追加到 `S_T`。

**步骤 B：为每个分区选择 sweep axis（可选但推荐）**

- 用 Tsitsigkos’19 的直方图模型计算 $I_T^x,I_T^y$，令 `SweepAxis[T] = argmin`. [arXiv+1](https://arxiv.org/pdf/1908.11740)
- 若不做模型：1D stripes 时可默认“分区轴与 sweep 轴正交”，论文调参时也强调同轴并不能随 $K$ 下降、而正交（xy/yx）会随 $K$ 改善。 [arXiv](https://arxiv.org/pdf/1908.11740)

**步骤 C：对每个分区排序**

- 若 `SweepAxis[T]=x`：按 $(x_l,\text{id})$ 排序 `R_T,S_T`
- 若 `SweepAxis[T]=y`：按 $(y_l,\text{id})$ 排序 `R_T,S_T`

> 这一步完成后，多次 Pass 只要复用同样的分区与排序结果，就能保证枚举顺序一致。

------

#### 3.1.2 分区内 join：Forward-Scan Plane Sweep（改写自 Algorithm 1）

对每个分区 $T$，在 `SweepAxis[T]` 上做 forward-scan：

**核心思想（来自论文 Algorithm 1）：**
 把 `R_T` 与 `S_T` 按 sweep 轴的 lower endpoint 排序；扫描时，每当遇到某个矩形的 lower endpoint，就在另一边 forward scan 出所有满足 upper ≥ other.lower 的候选，再检查另一维相交后输出。 [arXiv](https://arxiv.org/pdf/1908.11740)

------

#### 3.1.3 去重：DuplicateTest

对每个相交候选对 $(r,s)$：

- 若 `DuplicateTest(T,r,s)` 为真，则认为这是该 pair 的**唯一归属分区**，可“输出/计数/存储”。

该 test 直接使用 Tsitsigkos’19 的 Eq.(1) 及 stripes 简化形式。 [arXiv](https://arxiv.org/pdf/1908.11740)

------

#### 3.1.4 唯一枚举器接口（伪代码）

> 说明：以下伪代码是对 PBSM + Algorithm1 + Eq(1) 的“工程化整合写法”；不逐行照抄论文，但语义对应。 [arXiv+1](https://arxiv.org/pdf/1908.11740)

```
PROC EnumerateUniquePairs(R, S, Partitions, SweepAxis, callback):
    for each partition T in deterministic order:
        Rt = R_T[T]  // already sorted by SweepAxis[T]
        St = S_T[T]  // already sorted
        axis = SweepAxis[T]

        // Forward-scan plane sweep (Tsitsigkos’19 Algorithm 1 style)
        i = 0; j = 0
        while i < |Rt| and j < |St|:
            r = Rt[i]; s = St[j]
            if lower(axis,r) < lower(axis,s):
                // scan forward in S from j while overlap in axis possible
                j2 = j
                while j2 < |St| and upper(axis,r) >= lower(axis,St[j2]):
                    if IntersectOtherAxis(r, St[j2], axis):
                        if DuplicateTest(T, r, St[j2]):
                            callback(r, St[j2])
                    j2++
                i++
            else:
                // symmetric: scan forward in R from i while overlap in axis possible
                i2 = i
                while i2 < |Rt| and upper(axis,s) >= lower(axis,Rt[i2]):
                    if IntersectOtherAxis(Rt[i2], s, axis):
                        if DuplicateTest(T, Rt[i2], s):
                            callback(Rt[i2], s)
                    i2++
                j++
```

- `IntersectOtherAxis`：若 sweep=x 则检查 y 维是否相交，反之检查 x。
- `DuplicateTest`：
  - 2D grid：用 Eq.(1)
  - 1D vertical stripes：只检查 $\max(x_l)\ge T.x_l$

------

### 3.2 版本一：Enumerate+Sampling（存下全部 $J$ 再采样）

**适用场景：** $|J|$ 不大、内存充足；想要“一遍 join 完事”。

#### 流程

1. **预处理**：分区 + sweep axis + 排序（3.1.1）。
2. **枚举并存储所有 pair**：

```
Pairs = []
EnumerateUniquePairs(..., callback = (r,s) => Pairs.append((r.id, s.id)))
W = |Pairs|
```

1. **i.i.d. 有放回抽样**：

```
if W == 0: return empty
for j in 1..t:
    idx ~ Unif{0..W-1}
    Ans[j] = Pairs[idx]
return Ans
```

#### 关键性质

- 采样严格均匀且独立，因为直接对数组做独立均匀下标采样。

------

### 3.3 版本二：Sampling（两遍 PBSM 枚举 + 随机索引定位，精确 i.i.d.）

**适用场景：** $|J|$ 可能很大，不允许存下 $J$，但仍想要“精确 i.i.d. uniform”。

#### 流程概览

- Pass1：只计数 $|J|=W$（不存 pair）
- 生成 $t$ 个独立均匀索引 $I_1,\dots,I_t$
- Pass2：重新枚举同一个“唯一 join 流”，当计数器走到这些索引位置时输出对应 pair（重复索引自然支持有放回）

#### 详细步骤

**Step 0：预处理一次**（分区 + sweep axis + 排序），并复用到两次 Pass（保证顺序一致）。

**Pass 1：计数**

```
W = 0
EnumerateUniquePairs(..., callback = (r,s) => W++)
```

> 要用 64-bit 计数器，避免溢出。

**生成随机索引（有放回）**
$$
I_1,\dots,I_t \stackrel{i.i.d.}{\sim} \mathrm{Unif}\{0,\dots,W-1\}.
$$
工程上做：

- 构造 `IdxPos = [(I_j, j)]`
- 按 `I_j` 升序排序（便于 Pass2 一次线性扫描填完全部样本）

**Pass 2：定位输出**

```
Ans[1..t]
p = 1        // pointer in IdxPos
c = 0        // global rank in the unique join stream

EnumerateUniquePairs(..., callback = (r,s) => {
    while p <= t and IdxPos[p].index == c:
        Ans[IdxPos[p].pos] = (r.id, s.id)
        p++
    c++
})
return Ans
```

- 若多个索引相同，会在同一 $c$ 上输出同一个 pair 多次 → **有放回**。

------

### 3.4 版本三：Adaptive+Sampling（自动在“存全量”与“两遍采样”之间切换）

**目标：**
 小 $|J|$：走 Enumerate+Sampling（快且一遍）
 大 $|J|$：走 Sampling（不爆内存）

我们引入阈值 $J_\star$。

#### Adaptive Phase1：一次扫描（计数 + 尝试枚举到阈值）

> 注意：和很多“边枚举边切换”的实现不同，这里设计成**“永远计数，枚举只做可选附加动作”**，保证逻辑干净。

```
mode = ENUMERATE
Pairs = []     // capacity J_star
W = 0

EnumerateUniquePairs(..., callback = (r,s) => {
    W++
    if mode == ENUMERATE:
        if W <= J_star:
            Pairs.append((r.id, s.id))
        else:
            mode = COUNT_ONLY
            Pairs.clear()   // discard partial enumeration
})
```

结束后：

- 若 `mode == ENUMERATE` 则必有 $W \le J_\star$，且 `Pairs` 包含全体 $J$
- 若 `mode == COUNT_ONLY` 则 $W > J_\star$，且 `Pairs` 已丢弃

#### 分支 A（小 join）：直接数组采样（无需 Pass2）

若 `mode==ENUMERATE`：

```
if W == 0: return empty
for j in 1..t:
    idx ~ Unif{0..W-1}
    Ans[j] = Pairs[idx]
return Ans
```

#### 分支 B（大 join）：走两遍 Sampling 的 Pass2

若 `mode==COUNT_ONLY`：

1. 生成 $t$ 个随机索引并排序（同 3.3）
2. 再跑一遍 `EnumerateUniquePairs`，按索引输出（同 3.3 Pass2）

------

## 4. Adaptive 版本阈值 $J_\star$ 的选择策略

阈值的设计目标是：**既不 OOM，又尽可能减少“不必要的第二遍”。**

### 4.1 内存硬约束（必须满足）

设可用内存预算为 `MemBudget`（字节），允许给 `Pairs` 使用比例 $\rho\in(0,1)$。
 若每个 pair 存储为两个 32-bit id，则 `sizeof(pair)≈8 bytes`（实际还要考虑 vector 扩容/对齐，工程上常按 16 bytes 估更稳）。
$$
J_\star^{\text{mem}}
=\left\lfloor\frac{\rho\cdot \text{MemBudget}}{\text{sizeof(pair)}}\right\rfloor.
$$
必须取：
$$
J_\star \le J_\star^{\text{mem}}.
$$

------

### 4.2 时间权衡（“一遍存全量” vs “两遍不存”）

粗略比较：

- Enumerate+Sampling 时间：
   $T_{\text{enum}}\approx T_{\text{join}} + \text{WriteCost}\cdot W + O(t)$
- Sampling 两遍时间：
   $T_{\text{2pass}}\approx 2\cdot T_{\text{join}} + O(t\log t)$

其中 `WriteCost·W` 来自把 $W$ 个 pair 写入内存（对大 $W$ 可能很贵）。

**交叉点** $J_\star^{\text{time}}$ 可以用基准测试拟合：
 在同一机器上跑若干组数据，测得在 $W=W_{\text{cross}}$ 时两者耗时相当，则：
$$
J_\star^{\text{time}} \approx W_{\text{cross}}.
$$
最终建议：
$$
J_\star = \min\left(J_\star^{\text{mem}},\ J_\star^{\text{time}}\right).
$$

------

### 4.3 与 Tsitsigkos’19 的调参思想协同（建议）

Tsitsigkos’19 明确指出：分区数增加会提升性能但过大时复制过度会变差，并给出“分区尺度比对象尺度大一个数量级”的经验规则。你们可以把这条经验用于选择 $K$，从而影响 $T_{\text{join}}$ 的常数；这会进一步影响“是否值得两遍”。 [arXiv](https://arxiv.org/pdf/1908.11740)

此外，1D 分区的调参中也强调“分区轴与 sweep 轴正交”的组合（xy/yx）随 $K$ 增大性能提升更明显，而同轴（xx/yy）不太受益。对 baseline 性能很关键。 [arXiv](https://arxiv.org/pdf/1908.11740)

------

## 5. 算法分析（正确性 + 复杂度；三个版本都包含）

### 5.1 正确性分析

#### 引理 1（分区内 Forward-Scan 能枚举出分区内所有相交对）

对给定分区 $T$ 的两组列表 $R_T,S_T$，Forward-Scan plane sweep 会输出所有在该分区内相交的 pair（未去重前）。该逻辑对应 Tsitsigkos’19 的 Algorithm 1：按 lower endpoint 排序并 forward scan，保证所有在 sweep 轴上可能重叠者都会被检查，另一维检测过滤真实相交。 [arXiv](https://arxiv.org/pdf/1908.11740)

> 我们实现时只要保持与 Algorithm 1 等价的扫描逻辑即可。

------

#### 引理 2（DuplicateTest 保证每个全局 pair 恰好被“报告一次”）

PBSM multi-assignment 会让同一对 $(r,s)$ 出现在多个分区 join 任务中。Tsitsigkos’19 使用 reference-point 的 duplicate elimination：只在满足 Eq.(1) 的 tile 报告结果；在 1D stripes 情况下进一步简化为单维比较。 [arXiv](https://arxiv.org/pdf/1908.11740)

因此，对任意 $(r,s)\in J$，在所有可能发现它的分区中，**只有唯一一个分区**会通过 duplicate test，从而该 pair 被输出一次且仅一次。

------

#### 推论 1（`EnumerateUniquePairs` 枚举得到的是 $J$ 的一一对应流）

由引理 1（分区内不漏）与引理 2（跨分区不重且不漏）可得：
 `EnumerateUniquePairs` 在所有分区上输出的 pair 流，等价于对 $J$ 的一次完整枚举（每个元素恰好一次）。

------

#### 定理 1（Enumerate+Sampling 输出 i.i.d. 均匀有放回样本）

`Pairs` 数组是 $J$ 的一一枚举，长度 $W=|J|$。
 每次采样独立生成 $\text{idx}\sim \mathrm{Unif}\{0,\dots,W-1\}$，输出 `Pairs[idx]`。
 显然得到：
$$
\Pr(Z_j=p)=\frac1{|J|}
\quad\text{且}\quad Z_1,\dots,Z_t\ \text{独立}.
$$

------

#### 定理 2（Sampling 两遍索引定位 输出 i.i.d. 均匀有放回样本）

Pass1 得到 $W=|J|$ 精确值；生成的索引 $I_j$ 独立均匀。
 Pass2 的枚举流是 $J$ 的固定一一枚举（推论 1）。设 pair $p$ 在该枚举流中的 rank 为 $\mathrm{rank}(p)\in[0,W-1]$，则：
$$
\Pr(Z_j=p)=\Pr(I_j=\mathrm{rank}(p))=\frac1W=\frac1{|J|}.
$$
且 $I_j$ 独立 ⇒ 样本独立。

> 关键实现条件：两遍枚举必须产生完全一致的“唯一枚举流顺序”。这就是为什么我们要求“分区顺序固定 + 排序稳定 + tie-break by id”。

------

#### 定理 3（Adaptive+Sampling 不会破坏正确性）

Adaptive 只是在两条**各自正确**的路径间选一条：

- 若 $W\le J_\star$：走 Enumerate+Sampling（定理1正确）
- 若 $W>J_\star$：走 Sampling 两遍（定理2正确）

阈值只影响性能与内存，不影响样本分布的正确性。

------

### 5.2 复杂度分析（统一符号）

记：

- $n_R=|R|,\ n_S=|S|,\ n=n_R+n_S$
- 分区数：1D stripes 为 $K$，2D grid 为 $K^2$
- 复制后的总元素数（多重分配后落入分区列表的总长度）：

$$
n'=\sum_T |R_T|+\sum_T |S_T|.
$$

- 分区内排序开销：

$$
\text{SortCost}=\sum_T\bigl(|R_T|\log|R_T|+|S_T|\log|S_T|\bigr).
$$

- 分区内扫描/候选检查成本（与数据分布相关），记为：

$$
\text{SweepCost}=\sum_T \mathrm{JoinCost}(T).
$$

- 去重检查是 $O(1)$ 每对候选的常数开销（Eq.(1) 或其简化）。 [arXiv](https://arxiv.org/pdf/1908.11740)

------

### 5.3 预处理成本（所有版本共享）

- 分区（multi-assignment）：$O(n')$
- sweep axis 直方图估计（可选）：论文称该决策开销相对 join 很小，但严格计入也可写成每分区线性于样本量的轻量统计。 [arXiv+1](https://arxiv.org/pdf/1908.11740)
- 排序：$\text{SortCost}$

------

### 5.4 各版本时间复杂度

#### 5.4.1 Enumerate+Sampling

- 枚举（一次）：
  $$
  O\bigl(n' + \text{SortCost} + \text{SweepCost} + |J|\bigr)
  $$
  （输出 $|J|$ 个 pair 并写入数组）

- 采样：$O(t)$

综合：
$$
T_{\text{Enum+Samp}}=
O\bigl(n'+\text{SortCost}+\text{SweepCost}+|J|+t\bigr).
$$

------

#### 5.4.2 Sampling（两遍）

- Pass1 枚举计数：
   $O(n' + \text{SortCost} + \text{SweepCost} + |J|)$

- 索引生成：$O(t)$，索引排序：$O(t\log t)$

- Pass2 再枚举：
   $O(n' + \text{SweepCost} + |J|)$

  > 若复用已排序列表，则 Pass2 不必重复排序；若不复用则再加一次 $\text{SortCost}$。

综合（复用排序的常见实现）：
$$
T_{\text{Samp}}=
O\bigl(n'+\text{SortCost}+2\cdot\text{SweepCost}+2|J|+t\log t\bigr).
$$

> 你们写论文时可以强调：这是 baseline，时间仍然显式依赖 $|J|$（而你们方法目标是最坏情况不依赖 $|J|$）。

------

#### 5.4.3 Adaptive+Sampling

分两种情况：

- 若 $|J|\le J_\star$：走 Enumerate+Sampling
  $$
  T=O\bigl(n'+\text{SortCost}+\text{SweepCost}+|J|+t\bigr)
  $$

- 若 $|J|> J_\star$：走两遍 Sampling
   Phase1 中最多额外写入 $J_\star$ 个 pair 后丢弃（上限可控），其余只是计数；再做 Pass2。
  $$
  T=O\bigl(n'+\text{SortCost}+2\cdot\text{SweepCost}+2|J|+t\log t + J_\star\bigr).
  $$

$J_\star$ 只带来可控的“试存”开销。

------

### 5.5 各版本空间复杂度

#### 5.5.1 Enumerate+Sampling

- 分区列表：$O(n')$
- 存储全体 pair：$O(|J|)$
- 输出：$O(t)$

$$
S_{\text{Enum+Samp}}=O(n'+|J|+t).
$$

#### 5.5.2 Sampling（两遍）

- 分区列表：$O(n')$
- 索引 + 输出：$O(t)$
- 不存全体 pair

$$
S_{\text{Samp}}=O(n'+t).
$$

#### 5.5.3 Adaptive+Sampling

- 分区列表：$O(n')$
- 最多暂存 $J_\star$ 个 pair（小分支保留，大分支最终会清空）
- 索引 + 输出：$O(t)$

$$
S_{\text{Adaptive}}=O(n'+J_\star+t).
$$

------

## 额外工程建议（保证可复现实验 + “严格 i.i.d.” 不出坑）

1. **两遍一致性（Sampling / Adaptive 大分支必需）**
   - 分区遍历顺序固定（stripe id 升序，或 grid row-major）
   - 每分区排序稳定 + id tie-break
   - Forward-scan 的扫描推进顺序固定
2. **计数器与随机数**
   - $|J|$ 可能很大：用 64-bit
   - 生成 $[0,W-1]$ 的均匀整数：要用无偏的 64-bit uniform（不要 `% W` 直接取模导致偏差，除非使用 rejection sampling）
3. **分区参数 $K$**
   - 推荐用论文经验规则作为默认：分区尺度比对象尺度大约一数量级；并在实验中扫描少量 $K$ 值验证。 [arXiv](https://arxiv.org/pdf/1908.11740)
4. **优先 1D stripes**
   - 论文的调参与结论指出 1D partitioning 往往更高效，并强调正交的 partition/sweep 组合（xy 或 yx）会随 $K$ 改善，而同轴组合（xx/yy）不太受益。 