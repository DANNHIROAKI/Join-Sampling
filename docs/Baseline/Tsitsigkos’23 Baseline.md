# 1. 问题定义与分析

## 1.1 数据模型与 join 语义（filter step）

给定两组空间对象，仅考虑其 MBR（Minimum Bounding Rectangle）作为过滤步（filter step）的输入：

- $R=\{r_1,\dots,r_{n_R}\}$
- $S=\{s_1,\dots,s_{n_S}\}$

每个对象用二维轴对齐矩形表示。为避免“贴边是否算相交”的歧义，本 baseline 统一采用**半开区间**表示（你们当前草稿也如此）： 
$$
r=[x_l(r),x_u(r))\times[y_l(r),y_u(r)),\qquad
s=[x_l(s),x_u(s))\times[y_l(s),y_u(s)).
$$
我们要求每个矩形满足 $x_l < x_u$ 且 $y_l < y_u$（否则面积为 0，按过滤步通常可忽略）。

**相交连接（filter step）结果集合：**
$$
J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq \varnothing\}.
$$
在半开语义下，矩形相交等价于两个维度投影都严格相交：
$$
\max(x_l(r),x_l(s))<\min(x_u(r),x_u(s))
\ \land\
\max(y_l(r),y_l(s))<\min(y_u(r),y_u(s)).
$$

> 备注（写进论文更稳）：如果你们后续实验需要对齐 GIS 里“边界接触也算相交（闭区间）”的语义，则需要改不等号与 tie‑break；本 baseline 固定采用半开语义以便描述清晰、实现确定。

------

## 1.2 采样目标：在 $J$ 上 i.i.d. 均匀有放回

输出 $t$ 个样本对 $Z_1,\dots,Z_t\in J$，满足：

- **均匀性**：$\Pr(Z_j=P)=1/|J|$（任意 $P\in J$）
- **独立性**：$Z_1,\dots,Z_t$ 相互独立（i.i.d.）
- **有放回**：允许同一 pair 被多次抽到（with replacement）

这是 spatial join sampling 的标准目标。 

------

## 1.3 为什么 SOP（网格 + multi‑assignment）会产生 duplicates？

在 SOP（space‑oriented partitioning，例如规则网格）中，空间被切成**不重叠**的 tiles。非点对象（MBR）可能跨越多个 tile，因此采用 **multi‑assignment（多重分配）**：一个矩形会被复制/登记到所有与其相交的 tiles。

问题是：当我们在每个 tile 内把 $R$ 与 $S$ 的对象子集做局部 join 时，**同一对相交矩形**可能在多个 tiles 内同时出现，从而产生 duplicates。

传统去重方式包括：哈希去重、reference‑point 去重（对每个候选结果计算交集的参考点，只在包含参考点的 tile 输出）。Tsitsigkos’23 的核心动机就是：reference‑point 虽然比 hash 快，但仍会**生成重复候选**并付出额外判断代价。 

------

# 2. 引用文章内容（Tsitsigkos’23）与本 baseline 继承点

## 2.1 文章信息

- **标题**：*Two-layer Space-oriented Partitioning for Non-point Data*
- **作者**：Dimitrios Tsitsigkos, Panagiotis Bouros, Konstantinos Lampropoulos, Nikos Mamoulis, Manolis Terrovitis
- **年份/版本**：arXiv:2307.09256v1（2023‑07‑18） 

论文研究 main‑memory 场景下，用 SOP（特别是 grid）处理非点对象 MBR 的范围查询与相交 join，并提出 “secondary partitioning（二层次分区）” 来**避免 duplicates 的生成**。 

------

## 2.2 Tsitsigkos’23 的 two‑layer SOP：A/B/C/D 次分区

在每个 tile $T$ 内，把分配到 $T$ 的矩形按其 “begin（start）是否在 tile 起点之后”分类为四类：A/B/C/D。论文在 Section 3 给出了这一 secondary partitioning，并配图说明。 

设 tile $T$ 的起点（左下角）为 $(T.x_l, T.y_l)$。对矩形 $r$：

- **A 类**：$x_l(r)\ge T.x_l$ 且 $y_l(r)\ge T.y_l$
- **B 类**：$x_l(r)\ge T.x_l$ 且 $y_l(r)< T.y_l$
- **C 类**：$x_l(r)< T.x_l$ 且 $y_l(r)\ge T.y_l$
- **D 类**：$x_l(r)< T.x_l$ 且 $y_l(r)< T.y_l$

直观解释：
 A “从 tile 内开始”；B “x 内开始、y 早于 tile”；C “y 内开始、x 早于 tile”；D “两维都早于 tile”。 

------

## 2.3 Tsitsigkos’23 的 join 关键结论：16 个 mini‑join 中仅需评估 9 个即可无重复

对同一 tile $T$，把 $R$ 分成 $R_T^A,R_T^B,R_T^C,R_T^D$，把 $S$ 分成 $S_T^A,S_T^B,S_T^C,S_T^D$。那么该 tile 内的 join 可以分解为 $4\times 4=16$ 个 class‑to‑class 的 **mini‑joins**。论文在 Section 5.1.1 指出：其中 **7 个 mini‑join 只会产生 duplicates**，因此可以完全跳过；只评估剩余 **9 个 mini‑join** 即可得到全局无重复 join 结果，无需 reference‑point 去重。 

论文给出直观例子：若在 tile $T$ 内两条 **B 类**矩形相交（即结果属于 $R_T^B \Join S_T^B$），那么它们必然也会在 $T$ 上方的某个 tile 内同时出现并被报告，因此当前 tile 的输出必定是重复。 

------

## 2.4 Tsitsigkos’23 的 mini‑join 枚举器（plane‑sweep）与优化

论文在 Section 5.1.2 给了用于 mini‑join 的 plane‑sweep 框架（Algorithm 2），以及两类优化： 

1. **Reduced plane‑sweep**（Algorithm 3）：在某些 class 组合下，由定义可知一边的起点必在另一边之前，从而可以省掉一侧排序，并用简化扫描；论文还指出总体上可不必排序 $R_T^C,R_T^D,S_T^C,S_T^D$ 等。 
2. **Batch output / Avoid redundant comparisons**（Algorithm 4）：进一步减少冗余比较，通过按 $x_u$ 排序并批量确认一段候选。 

> 本 baseline 的原则：**优先保证正确性 + 可复现**；优化作为可选增强，写清楚“实现最低配也能工作”。

------

# 3. 核心数据结构

## 3.1 主分区：规则网格 tiles（SOP）

设空间域（或数据包围盒/归一化后的域）为：
$$
\mathcal{D}=[X_{\min},X_{\max})\times[Y_{\min},Y_{\max}).
$$
选网格分辨率 $(n_x,n_y)$，tile 宽高为：
$$
w=\frac{X_{\max}-X_{\min}}{n_x},\qquad
h=\frac{Y_{\max}-Y_{\min}}{n_y}.
$$
第 $(i,j)$ 个 tile 定义为：
$$
T_{i,j}=[X_{\min}+iw,\ X_{\min}+(i+1)w)\times [Y_{\min}+jh,\ Y_{\min}+(j+1)h),
\quad 0\le i<n_x,\ 0\le j<n_y.
$$
**tile id（用于确定性顺序）**：固定采用 row‑major：
$$
\text{id}(i,j)=j\cdot n_x+i.
$$

------

## 3.2 Multi‑assignment：对象到 tiles 的分配

矩形 $r$ 分配到所有与其相交的 tiles。工程实现通常通过计算覆盖的 index 范围来完成：
$$
i_{\min}=\left\lfloor \frac{x_l(r)-X_{\min}}{w}\right\rfloor,\quad
i_{\max}=\left\lfloor \frac{x_u(r)-X_{\min}-\varepsilon}{w}\right\rfloor,
$$
其中 $\varepsilon$ 是一个“半开区间”意义上的微小量（实现上可用 `nextafter` 或“若坐标离散则直接用整数边界”来避免边界落在上界导致多分配）。然后对所有 $i\in[i_{\min},i_{\max}], j\in[j_{\min},j_{\max}]$ 做插入。

> 如果你们数据坐标本身是整数网格（常见于栅格化/归一化），可以把 $\varepsilon$ 的问题转化为“上界是否包含”的整数公式，从而完全确定。

------

## 3.3 次分区：tile 内 A/B/C/D（Two‑layer SOP）

每个 tile $T$ 存 8 个列表：
$$
R_T^A,R_T^B,R_T^C,R_T^D;\quad
S_T^A,S_T^B,S_T^C,S_T^D.
$$
每个列表元素建议存储（MBR + id）：
$$
\langle id,\ x_l,\ x_u,\ y_l,\ y_u\rangle.
$$
**可选空间优化（来自论文）**：在固定以 $x$ 为 sweep 维度时，某些类在 $y_l$ 或 $x_l$ 上可省存（论文 Section 5.2.3 的讨论）。baseline 允许但不强制实现。 

------

## 3.4 9 个 mini‑join 列表（固定顺序）

对每个 tile $T$，只评估以下 9 个 mini‑join（固定顺序，用于 JoinStream 的确定性输出）：

1. $R_T^A \Join S_T^A$
2. $R_T^A \Join S_T^B$
3. $R_T^A \Join S_T^C$
4. $R_T^A \Join S_T^D$
5. $R_T^B \Join S_T^A$
6. $R_T^B \Join S_T^C$
7. $R_T^C \Join S_T^A$
8. $R_T^C \Join S_T^B$
9. $R_T^D \Join S_T^A$

这正对应 Tsitsigkos’23 所述“16 个中跳过 7 个 shaded mini‑joins，保留 9 个即可无重复”。 

------

## 3.5 JoinStream：确定性 join 输出流（核心接口）

我们把“全局枚举 join pairs”抽象为一个生成器：
$$
\texttt{JOIN\_STREAM(Index)} \to (r,s)\in J
$$
要求：

- **完备**：输出集合等于 $J$
- **无重复**：每个 pair 恰好输出一次
- **确定性顺序**：相同输入与参数下，两次运行输出顺序完全一致

确定性对 2‑pass index sampling 和 Adaptive 的大分支至关重要（否则第二遍索引定位会错，导致采样失真）。你们草稿也强调了这一点，这里把它上升为 baseline 的硬约束。 

------

## 3.6 Sampling 相关辅助结构（按版本）

- Sampling（2‑pass Index）：
  - `K[1..t]` 随机索引
  - `K_sorted`（按索引排序的 (value, pos)）
  - `Ans[1..t]`
- Enumerate+Sampling：
  - `AllPairs[0..N-1]`，其中 $N=|J|$
- Adaptive：
  - 阈值 $J_\star$
  - `mode ∈ {ENUMERATE, COUNT_ONLY}`
  - `AllPairs`（最多临时存到 $J_\star$）

------

# 4. 算法详细流程（Sampling / Enumerate+Sampling / Adaptive+Sampling）

本章给出“可直接实现”的流程与伪代码，读者无需阅读原论文即可复现。

------

## 4.1 BuildIndex：构建 two‑layer grid 索引

### 4.1.1 输入输出

- 输入：集合 $R,S$，网格参数 $(n_x,n_y)$，域 $\mathcal{D}$
- 输出：Index（包含所有非空 tiles 及其 8 个列表）

### 4.1.2 伪代码（推荐两遍分配，可选）

```
BuildIndex(R, S, nx, ny, Domain):
    init tiles[0 .. nx*ny-1] each with 8 empty lists

    // assign R
    for each r in R:
        for each tile T overlapped by r (multi-assignment):
            X = Classify(r, T)  // A/B/C/D
            append r to tiles[T].R[X]

    // assign S
    for each s in S:
        for each tile T overlapped by s (multi-assignment):
            Y = Classify(s, T)
            append s to tiles[T].S[Y]

    // determinism: sort each list once (stable)
    for each non-empty tile T in increasing tile_id:
        for X in {A,B,C,D}:
            stable_sort(tiles[T].R[X], key = (x_l, y_l, x_u, y_u, id))
            stable_sort(tiles[T].S[X], key = (x_l, y_l, x_u, y_u, id))

    return tiles
```

> 说明：
>
> - 排序键你们也可以简化为 `(x_l, id)`，但建议把 `(y_l, x_u, y_u)` 纳入 key 以减少浮点 tie 时的不确定性（尤其不同平台/编译器）。
> - 若使用 sparse 存储（如只存 non-empty tiles），必须保证 tile 遍历顺序仍然是稳定的（如把 tile_id 放入 vector 后排序）。

------

## 4.2 JoinStream：全局无重复枚举器（确定性输出流）

### 4.2.1 JoinStream 的全局顺序（必须固定）

1. tile 顺序：tile_id 升序（row‑major）
2. tile 内 mini‑join 顺序：固定为上文 9 个列表
3. mini‑join 内输出顺序：由 mini‑join 枚举器的排序与扫描规则决定，必须确定

这些约束的目的，是让 `JOIN_STREAM()` 在不同 pass 中输出完全一致的序列。 

### 4.2.2 伪代码

```
JOIN_STREAM(Index):
    for tile_id from 0 to nx*ny-1:
        T = Index[tile_id]
        if T is empty: continue

        // fixed mini-join order
        yield all pairs from MiniJoin(T.R[A], T.S[A])
        yield all pairs from MiniJoin(T.R[A], T.S[B])
        yield all pairs from MiniJoin(T.R[A], T.S[C])
        yield all pairs from MiniJoin(T.R[A], T.S[D])

        yield all pairs from MiniJoin(T.R[B], T.S[A])
        yield all pairs from MiniJoin(T.R[B], T.S[C])

        yield all pairs from MiniJoin(T.R[C], T.S[A])
        yield all pairs from MiniJoin(T.R[C], T.S[B])

        yield all pairs from MiniJoin(T.R[D], T.S[A])
```

------

## 4.3 mini‑join 枚举器：Plane‑Sweep（deterministic）

> 这里我们给一个“baseline 级可复现”的 plane‑sweep 实现框架。你们可以实现为 Tsitsigkos’23 Algorithm 2 的等价版本，并可选加入 Algorithm 3/4 的优化。 Tsitsigkos’23

### 4.3.1 最低配 MiniJoin（正确性优先）

输入：两个已排序列表 $A,B$（按稳定键排序）
 输出：所有相交 pair $(a,b)$

**关键点**：必须完全确定性（同输入同输出顺序）。

一个容易实现且确定的做法（“双指针 + 活跃窗口”）如下：

- 列表均按 $x_l$ 升序
- 对每个 $a\in A$，维护 $B$ 中所有满足 $x_l(b) < x_u(a)$ 的候选窗口（用指针推进），再在窗口内筛掉 $x_u(b)\le x_l(a)$ 的元素（可用队列/指针）
- 对候选做 y 相交检查并输出

伪代码（简化版）：

```
MiniJoin(A, B):
    // A, B are sorted by (x_l, y_l, x_u, y_u, id)
    j_end = 0
    ActiveB = empty list  // store indices of B in deterministic order

    for each a in A:
        // extend ActiveB with new b whose x_l < a.x_u
        while j_end < |B| and B[j_end].x_l < a.x_u:
            ActiveB.push_back(j_end)
            j_end++

        // remove b that cannot overlap in x anymore: b.x_u <= a.x_l
        while ActiveB not empty and B[ActiveB.front].x_u <= a.x_l:
            pop_front(ActiveB)

        // now for all active b, x-overlap holds; check y-overlap
        for each idx in ActiveB in order:
            b = B[idx]
            if Intersect2D(a, b):   // strict by half-open
                yield (a, b)
```

> 复杂度：最坏情况下（大量重叠）会接近输出规模；但 baseline 重点是可复现与正确性。若你们希望更贴近论文性能，可替换为更经典的 plane sweep（带平衡树/区间树）或直接实现论文的 reduced/batch 版本。 

------

## 4.4 版本一：Sampling（Two‑Pass Index Sampling）

**名称**：TS23‑TLSOP‑2PassIndex
 **思想**：

- Pass1：只计数 $N=|J|$
- 抽 $t$ 个独立均匀索引 $K_j\sim\text{Unif}\{1..N\}$（有放回）
- Pass2：再次按同一 JoinStream 枚举，命中索引时输出

### 4.4.1 伪代码

```
TS23_TLSOP_2PassIndex(R, S, t, nx, ny, Domain):

1) Index = BuildIndex(R, S, nx, ny, Domain)

2) // Pass-1: count
   N = 0
   for pair in JOIN_STREAM(Index):
       N++

   if N == 0 or t == 0:
       return empty

3) // draw indices with replacement
   for j = 1..t:
       K[j] ~ UniformInt(1, N)
   K_sorted = sort pairs (K[j], j) by K ascending, tie by j

4) // Pass-2: select by index
   idx = 0
   p = 1
   Ans[1..t]
   for pair in JOIN_STREAM(Index):
       idx++
       while p <= t and K_sorted[p].value == idx:
           Ans[K_sorted[p].pos] = pair
           p++
       if p > t: break

5) return Ans
```

### 4.4.2 实现注意

- `K_sorted` 的 tie by j 可以保证完全确定性（即使多个样本抽到同一 idx，也按样本位置递增回填）。
- 若 `JOIN_STREAM` 两遍输出顺序不一致，则该算法会**直接错误**（不是“变慢”，而是输出分布不对）。因此必须满足第 3.5 节的确定性约束。 

------

## 4.5 版本二：Enumerate+Sampling（显式枚举后数组采样）

**名称**：TS23‑TLSOP‑EnumerateThenSample
 **思想**：一次 JoinStream 把全部 pairs 物化到数组 `AllPairs`，然后在数组上 i.i.d. 均匀采样。

### 4.5.1 伪代码

```
TS23_TLSOP_EnumerateThenSample(R, S, t, nx, ny, Domain):

1) Index = BuildIndex(R, S, nx, ny, Domain)

2) AllPairs = []
   for pair in JOIN_STREAM(Index):
       AllPairs.append(pair)

   N = |AllPairs|
   if N == 0 or t == 0:
       return empty

3) Ans[1..t]
   for j = 1..t:
       idx ~ UniformInt(0, N-1)
       Ans[j] = AllPairs[idx]

4) return Ans
```

> 优点：只做 1 次 join stream，通常最快。
>  缺点：空间 $\Theta(|J|)$ 可能爆内存。

------

## 4.6 版本三：Adaptive+Sampling（阈值自适应）

**名称**：TS23‑TLSOP‑Adaptive
 **目标**：

- 若 $|J|\le J_\star$：退化为 Enumerate+Sampling（1 pass，快）
- 若 $|J|> J_\star$：退化为 2‑PassIndex（2 pass，省内存）

### 4.6.1 推荐的切换逻辑（准确计数 + 最多白存 $J_\star$）

第一遍扫描 JoinStream 时始终做计数 $N++$，并在 $N\le J_\star$ 时暂存 pairs；一旦 $N>J_\star$，立即清空并停止存储，后续只计数直到结束。这样可得到精确 $N=|J|$。

### 4.6.2 伪代码

```
TS23_TLSOP_Adaptive(R, S, t, J_star, nx, ny, Domain):

1) Index = BuildIndex(R, S, nx, ny, Domain)

2) mode = ENUMERATE
   AllPairs = []
   N = 0

   // Pass-1: count (and maybe enumerate)
   for pair in JOIN_STREAM(Index):
       N++
       if mode == ENUMERATE:
           if N <= J_star:
               AllPairs.append(pair)
           else:
               mode = COUNT_ONLY
               AllPairs.clear()    // release memory

3) if N == 0 or t == 0:
       return empty

4) if mode == ENUMERATE:
       // N <= J_star, AllPairs is full J
       for j = 1..t:
           idx ~ UniformInt(0, N-1)
           Ans[j] = AllPairs[idx]
       return Ans

5) // mode == COUNT_ONLY: do 2-pass index sampling
   for j = 1..t:
       K[j] ~ UniformInt(1, N)
   K_sorted = sort (K[j], j) by K ascending, tie by j

   idx = 0
   p = 1
   for pair in JOIN_STREAM(Index):
       idx++
       while p <= t and K_sorted[p].value == idx:
           Ans[K_sorted[p].pos] = pair
           p++
       if p > t: break

6) return Ans
```

------

## 4.7 复现与确定性规则（必须写进实验/实现说明）

为保证两遍流一致（Sampling / Adaptive 大分支），你们必须明确并实现以下规则： 

1. **tile 顺序固定**：tile_id 升序（row‑major），不要依赖 unordered_map 的遍历顺序
2. **mini‑join 顺序固定**：严格按 9 个组合固定列表
3. **类内排序稳定**：稳定排序，tie‑break 必须包含 id（建议用多字段 key）
4. **mini‑join 输出顺序固定**：不能用不确定迭代容器；扫描逻辑必须完全确定
5. **浮点比较一致**：若跨平台复现，建议把输入坐标归一化为整数/定点或明确比较策略（例如使用 IEEE `nextafter` 处理半开边界）

------

# 5. Adaptive 阈值 $J_\star$ 的选择策略

阈值设计要满足：**不 OOM（硬约束）** + **尽可能快（软约束）**。

------

## 5.1 内存硬约束（必须满足）

设给 `AllPairs` 的可用内存预算为 `MemBudgetPairs` 字节。每个 pair 存储开销近似：

- 两个 32‑bit id：8 字节，或
- 两个 64‑bit id/指针：16 字节

再乘容器与对齐开销系数 $c_{\text{over}}\in[1.2,2]$。

定义：
$$
J_\star^{\text{mem}}
=
\left\lfloor
\frac{\text{MemBudgetPairs}}{\text{PairBytes}\cdot c_{\text{over}}}
\right\rfloor.
$$
推荐取保守值（预留索引结构、排序缓冲、线程栈等）：
$$
J_\star = \left\lfloor 0.7\cdot J_\star^{\text{mem}}\right\rfloor.
$$

------

## 5.2 时间软约束（可选：离线标定 / 交叉点）

对固定硬件与实现：

- Enumerate+Sampling：1 次 JoinStream + 写 `AllPairs` + $O(t)$ 采样
- Sampling（2‑pass）：2 次 JoinStream + $O(t\log t)$ 排序索引

实际中，“写 AllPairs 的内存带宽 + cache miss”可能很贵，因此当 $|J|$ 较大时两遍反而更快。

推荐在 2–3 个代表性数据集上做 pilot：

1. 统计不同 $|J|$ 下 Enumerate 与 2‑pass 的耗时；
2. 找到二者耗时交叉点 $J_\star^{\text{time}}$；
3. 取最终阈值：

$$
J_\star = \min(J_\star^{\text{mem}},\ J_\star^{\text{time}}).
$$

------

## 5.3 实用默认（论文写法最稳）

如果你们不想引入时间标定，直接写：

> “$J_\star$ 由 pair 物化的内存预算决定；若超出阈值则切换到 2‑pass sampling 以避免 OOM。”

这在 baseline 报告里足够合理、可复现。

------

# 6. 算法分析（正确性 + 复杂度；三版本都包含）

------

## 6.1 正确性（一）：为什么只做 9 个 mini‑join 可以无重复枚举？

这里给一个**读者无需看原论文也能理解**的证明框架（与论文直觉一致）。 

### 定义：pair 的“参考点”（不用实际计算，只用于证明）

对任意相交 pair $(r,s)\in J$，定义其重叠区域的“左下参考点”：
$$
p(r,s)=\big(\max(x_l(r),x_l(s)),\ \max(y_l(r),y_l(s))\big).
$$
由于 $r,s$ 相交，参考点 $p(r,s)$ 一定位于二者重叠区域内。

因为 tiles 彼此不相交且采用半开划分，域内任意点都属于唯一一个 tile，所以存在唯一 tile $T^*$ 使得：
$$
p(r,s)\in T^*.
$$
我们称 $T^*$ 为该 pair 的**归属 tile**。

------

### 关键引理 A：在归属 tile $T^*$ 中，$(r,s)$ 必落在 9 个 mini‑join 之一

因为 $p_x=\max(x_l(r),x_l(s))\in[T^*.x_l,T^*.x_u)$，所以至少有一个对象（$r$ 或 $s$）满足其 $x_l\ge T^*.x_l$。换句话说，至少一方在 $x$ 维属于 “A 或 B”（即 begin inside x）。同理，至少一方在 $y$ 维属于 “A 或 C”。

枚举所有可能组合，你会发现 $(r,s)$ 在 $T^*$ 中不可能落入被跳过的 7 种组合：
$$
(B,B),(B,D),(C,C),(C,D),(D,B),(D,C),(D,D)
$$
因为这些组合都意味着两者在某一维都 “begin before tile”，从而 $\max(\cdot)$ 会落在 tile 起点之前，矛盾于 $p(r,s)\in T^*$。因此在 $T^*$ 中，该 pair 必落入保留的 9 个 mini‑join 之一，并会被枚举输出。 

------

### 关键引理 B：在任何非归属 tile $T\ne T^*$ 中，$(r,s)$ 若同时出现，则必落入那 7 个被跳过的 mini‑join 之一

直观理由：若 $T\ne T^*$ 但两者都被分配到 $T$，那么 $T$ 必然是某个方向上的“后继 tile”，即其起点在 $x$ 或 $y$ 维上不小于 $T^*$ 的起点。此时参考点 $p(r,s)$ 不在 $T$ 中，说明在至少一个维度上，两者的 begin 均在 $T$ 起点之前（否则参考点会被拉进 $T$）。这恰好对应那 7 个会产生 duplicates 的组合；因此在 $T$ 内该 pair 会被**刻意跳过**。 

------

### 结论（JoinStream 无重不漏）

由引理 A：每个 $(r,s)\in J$ 至少在归属 tile $T^*$ 被输出一次；
 由引理 B：在所有 $T\ne T^*$ 的 tiles 中不会输出它；
 因此全局 JoinStream **不重不漏**，每个 pair 恰好输出一次。

这就是 “9 mini‑join 无重复枚举” 的正确性核心。

------

## 6.2 正确性（二）：三种采样版本均输出 i.i.d. 均匀有放回样本

### 引理 1：2‑PassIndex 是 i.i.d. 均匀

设 JoinStream 输出序列为 $P_1,\dots,P_N$（这是 $J$ 的一个固定编号，且 $N=|J|$）。算法抽 $K_j\sim\mathrm{Unif}\{1..N\}$，输出 $Z_j=P_{K_j}$。

则对任意 $P_i$：
$$
\Pr(Z_j=P_i)=\Pr(K_j=i)=\frac{1}{N}.
$$
又因为 $K_1,\dots,K_t$ 独立，故 $Z_1,\dots,Z_t$ i.i.d.。

------

### 引理 2：Enumerate+Sampling 是 i.i.d. 均匀

若 `AllPairs` 恰好包含每个 pair 一次，则在数组下标上做独立均匀采样显然得到 i.i.d. 均匀。

------

### 定理：Adaptive 切换不影响正确性

Adaptive 只有两种终态：

- 未切换：等价 Enumerate+Sampling
- 切换：等价 2‑PassIndex（第二遍按索引定位）

切换只改变“是否保留 AllPairs”，不改变 JoinStream，也不改变最终计数 $N=|J|$ 与索引分布，因此输出仍为 $J$ 上 i.i.d. 均匀有放回样本。 

------

## 6.3 复杂度分析（构建 + 一次 JoinStream + 三版本总成本）

### 6.3.1 记号

- $n=n_R+n_S$
- replication factor（multi‑assignment 复制系数）：$\rho_R,\rho_S$
- 索引规模：$n' = n_R\rho_R+n_S\rho_S$
- 9 个 mini‑join 集合：$\mathcal{M}$
- tile $T$ 内类规模：$|R_T^X|,|S_T^Y|$

------

### 6.3.2 BuildIndex（two‑layer grid）

- **分配时间**：$O(n_R\rho_R+n_S\rho_S)=O(n')$
- **排序时间**（若对每个 tile/class 只排序一次）：

$$
O\left(\sum_T \sum_{X\in\{A,B,C,D\}} |R_T^X|\log|R_T^X|
      +\sum_T \sum_{Y\in\{A,B,C,D\}} |S_T^Y|\log|S_T^Y|
\right).
$$

- **空间**：$O(n')$（存所有复制后的 (MBR,id) 记录）

> 你们草稿里 BuildIndex 空间公式缺一个右括号，这里已修正为 $O(n')$。 

------

### 6.3.3 一次 JoinStream 的代价 $T_{\text{stream}}$

把每个 tile 的 9 个 mini‑joins 的枚举成本加总：
$$
T_{\text{stream}}=
\sum_T \sum_{(X,Y)\in\mathcal{M}}
\Big(
T_{\text{MiniJoin}}(R_T^X,S_T^Y)
\Big).
$$
若 mini‑join 使用“plane‑sweep + 输出”，通常可写成：
$$
T_{\text{MiniJoin}}(A,B)=\text{sweep\_scan}(A,B)+O(\text{out}_{A,B}),
$$
其中 $\text{out}_{A,B}$ 是该 mini‑join 输出的 pair 数。总输出规模为 $|J|$，因此任何正确枚举算法至少需要 $\Omega(|J|)$ 的输出写入成本（这是不可避免的）。

------

### 6.3.4 Sampling（2‑pass Index）

- **时间**：

$$
T_{\text{Sampling}} = 2\cdot T_{\text{stream}} + O(t\log t).
$$

- **空间**：

$$
S_{\text{Sampling}} = O(n' + t).
$$

------

### 6.3.5 Enumerate+Sampling

- **时间**：

$$
T_{\text{Enum+Sampling}} = T_{\text{stream}} + O(|J|) + O(t).
$$

（其中 $O(|J|)$ 是把 pairs 写入 `AllPairs` 的物化成本。）

- **空间**：

$$
S_{\text{Enum+Sampling}} = O(n' + |J| + t).
$$

------

### 6.3.6 Adaptive+Sampling

- 若 $|J|\le J_\star$（未切换）：

$$
T = T_{\text{stream}} + O(|J|) + O(t),\qquad
S = O(n' + |J| + t)\le O(n' + J_\star + t).
$$

- 若 $|J|> J_\star$（切换）：

第一遍最多白存 $J_\star$ 个 pair 后停止物化，因此额外开销 $O(J_\star)$：
$$
T = 2\cdot T_{\text{stream}} + O(J_\star) + O(t\log t),\qquad
S_{\max} = O(n' + J_\star + t).
$$

------

## 6.4 baseline 的“适用性边界”说明（建议写进报告）

- 该 baseline 本质是“枚举驱动”的：即使在 2‑pass/Adaptive 大分支里，也需要完整扫描 JoinStream 至少一次来得到 $N=|J|$。因此当 $|J|$ 极大时（接近 $n^2$），它在时间上不可避免地重。
- 这正是它作为 baseline 的意义：它代表“工程上直观、但最坏情况会被输出规模拖垮”的方案；与 Ours（最坏情况不依赖 $|J|$ 的采样）形成清晰对比。