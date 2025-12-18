# SIRS‑JS Baseline 设计报告

**目标**：用 SIGMOD’21 SIRS 作为“range i.i.d. sampling 黑盒”，构造 Spatial Join Sampling 的三个 baseline：

- Sampling（SIRS‑JS‑Sampling）
- Enumerate+Sampling（SIRS‑JS‑Enumerate）
- Adaptive+Sampling（SIRS‑JS‑Adaptive）

并保证：输出 $t$ 个样本 **i.i.d. 且在真实 join 结果集合 $J$ 上均匀**（with replacement）。

------

## 1. 问题定义与分析

### 1.1 输入对象：两类半开盒（half-open boxes）

在 $d\ge 2$ 维欧氏空间 $\mathbb{R}^d$ 中，给定两类**轴对齐半开盒**集合：
$$
R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\}.
$$
每个盒子（半开）定义为：
$$
r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\qquad L_i(r)<R_i(r).
$$
我们只关心**跨集合相交对**（cross-set intersection pairs）：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$

### 1.2 输出目标：在 $J$ 上 i.i.d. 均匀有放回采样

输出 $t$ 个样本 $Z_1,\dots,Z_t\in J$，要求：

- **均匀性（uniform）**：$\Pr(Z_j=P)=1/|J|$，$\forall P\in J$
- **独立性（i.i.d. with replacement）**：$Z_1,\dots,Z_t$ 相互独立，允许重复

> 这对应你们整体论文对 Sampling/Adaptive baseline 的要求。

### 1.3 半开相交语义（严格不等号必须写清楚）

对半开盒 $a,b$：
$$
a\cap b\neq\varnothing
\iff
\forall i\in[d]:
\big(L_i(a)<R_i(b)\big)\wedge \big(L_i(b)<R_i(a)\big).
$$
等价改写（更贴合后续嵌入）：
$$
\forall i:
\big(L_i(b)<R_i(a)\big)\wedge \big(R_i(b)>L_i(a)\big).
$$
**贴边不算相交**：例如 $R_i(b)=L_i(a)$ 时在该维没有交叠，因此整体不相交。

### 1.4 baseline 的核心采样分解：两层采样 = 全局均匀

对任意 outer 盒子 $a$ 定义其 join 度数：
$$
\deg(a)=|\{b\in \text{inner}:(a,b)\in J\}|.
$$
则：
$$
|J|=\sum_{a\in \text{outer}} \deg(a).
$$
如果我们能做到：

1. 按 $\Pr(\text{pick }a)=\deg(a)/|J|$ 抽 outer；
2. 在给定 $a$ 的条件下，按 $\Pr(\text{pick }b\mid a)=1/\deg(a)$ 在匹配的 inner 中均匀抽 $b$；

则对任意 $(a,b)\in J$：
$$
\Pr((a,b))=\frac{\deg(a)}{|J|}\cdot \frac1{\deg(a)}=\frac1{|J|},
$$
从而得到**单次样本全局均匀**；若两层抽样都独立重复，则得到 **i.i.d.**。

baseline 的全部工作，就是让这两层抽样“可实现且可复现”，其中第二层由 SIRS 提供。

------

## 2. 引用文章内容：SIRS’21 的核心定义与我们复用的能力

### 2.1 文章信息（baseline 依赖的来源）

- 标题：**Spatial Independent Range Sampling**
- 作者：Dong Xie, Jeff M. Phillips, Michael Matheny, Feifei Li
- 会议：SIGMOD 2021
- DOI：**10.1145/3448016.3452806**（论文首页）

### 2.2 SIRS 的问题定义：Uniform SIRS（Definition 1）

SIRS 将问题形式化为：给定点集 $P\subset\mathbb{R}^d$、查询 MBR $R$、整数 $k$，返回 $k$ 个来自 $R\cap P$ 的**独立**均匀样本，每个点的概率为 $1/|R\cap P|$。

> baseline 复用的“黑盒能力”就是：
>  **对点集 $P$ 的任意 MBR 查询 $R$，返回 $k$ 个 i.i.d. uniform samples（with replacement）。**

### 2.3 SIRS 对“跨 query 独立性”的明确强调（对 join baseline 非常关键）

SIRS 论文在介绍中明确指出：独立性不仅要求同一 query 内的样本彼此独立，也要对**不同 queries 返回的样本**仍然独立。

这点对 join baseline 很关键，因为 join 的每个样本可能来自不同的 outer 盒子 $a$，也就是不同的查询 $Q(a)$。

### 2.4 SIRS 的通用采样框架与 Lemma 1（复杂度抽象）

SIRS 的核心思想是：把“在 MBR 内抽样”转化为“在一维线性布局上的若干连续索引区间里抽样”，对这些区间建一个 top-level alias，然后按 alias 选区间、再在区间内均匀抽点。

**Lemma 1（SIRS 论文）**给出抽象复杂度：
 如果能把查询 MBR $R$ 映射成 $I(n,R)$ 个“连续索引区间”，映射耗时 $M(n,R)$，则可在
$$
O(M(n,R)+I(n,R)+k)
$$
时间内返回 $k$ 个独立均匀样本。

baseline 在复杂度分析里引用 Lemma 1 的这种“抽象符号”写法最审稿友好。

### 2.5 Walker alias 方法（Background）：离散按权采样的 $O(1)$ 组件

SIRS 在 Background 部分回顾：对离散对象集合按权采样，Walker alias 能做到**线性建表 $O(n)$**、**单次采样 $O(1)$**，且每次采样独立。

baseline 里我们用 alias 来实现外层按 $\deg(a)$ 加权抽 outer（以及 SIRS 内部也用 alias 抽区间）。

### 2.6 SIRS 的具体 instantiation（你可选，但 baseline 建议默认 KD-tree）

SIRS 给出多种 instantiation，例如：

- **Z-order / ZV-tree**：查询映射到若干 z-value 连续区间；论文给出相应的复杂度讨论（Corollary 1）。
- **KD-tree (KDS)**：把查询映射到 $O(\sqrt{n})$ 个节点/区间，并给出 $O(\sqrt{n}+k)$ 的期望复杂度结论（Corollary 2，固定展示 $d=2$）。
- 也讨论了可推广到 R-tree 等层次索引（Sec. 3.4）。

同时 SIRS 论文说明：文中固定 $d=2$ 讲解，但预期可用于低维（文中提到 $d<10$，评估到 $d=7$）。

> **baseline 需要诚实写清楚的点**：我们做 box→point 嵌入后，维度变为 **$2d$**。因此当原始 $d$ 较大时，SIRS 的实际性能可能显著下降；这不是 correctness 问题，但会影响 baseline 的可扩展性解释。

------

## 3. 核心数据结构（可复现 API + 工程细节）

本节给出 baseline 的“实现接口”和“实现要点”。你可以把它当作一个独立模块来实现。

### 3.1 数据对象与映射

**Box**

- `id`: 全局唯一（建议 64-bit）
- `L[1..d], R[1..d]`：半开盒坐标（保证 $L_i）

**EmbeddedPoint（inner 盒子的 $2d$ 维嵌入点）**

- `coord[1..2d]`：$\phi(b)$
- `inner_id`：指回原盒子 $b$

**映射表**

- `inner_id -> Box`
- `outer_id -> Box`

### 3.2 关键对齐：把 box 相交 join 转成点集 MBR 查询

#### 3.2.1 嵌入：inner 盒子 $\to$ 点

对每个 inner 盒子 $b$，定义：
$$
\phi(b)=\big(L_1(b),\dots,L_d(b),\ -R_1(b),\dots,-R_d(b)\big)\in\mathbb{R}^{2d}.
$$
令点集：
$$
P=\{\phi(b)\mid b\in \text{InnerSet}\}.
$$

#### 3.2.2 查询：outer 盒子 $\to$ $2d$ 维 MBR

对 outer 盒子 $a$，构造：
$$
Q(a)= \prod_{i=1}^d [\underline{L}_i,\ R_i(a)) \ \times\ \prod_{i=1}^d [\underline{U}_i,\ -L_i(a)).
$$
其中下界只要“足够小”覆盖所有点即可，例如：
$$
\underline{L}_i := \min_{b\in \text{Inner}} L_i(b),\qquad
\underline{U}_i := -\max_{b\in \text{Inner}} R_i(b).
$$

#### 3.2.3 等价性（核心引理）

对任意 outer $a$ 和 inner $b$：
$$
a\cap b\neq\varnothing\quad\Longleftrightarrow\quad \phi(b)\in Q(a).
$$
解释（逐维）：

- $\phi(b)$ 的第 $i$ 维是 $L_i(b)$，落入 $[\,\underline L_i, R_i(a)\,)$ 等价于 $L_i(b)；
- $\phi(b)$ 的第 $d+i$ 维是 $-R_i(b)$，落入 $[\,\underline U_i, -L_i(a)\,)$ 等价于 $-R_i(b)<-L_i(a)\iff R_i(b)>L_i(a)$。

因此：
$$
\deg(a)=|P\cap Q(a)|,
$$
并且在 $P\cap Q(a)$ 上均匀采样一个点，就等价于在与 $a$ 相交的 $b$ 上均匀采样一个。

> 这部分是 baseline 正确性的“硬基石”。

------

### 3.3 严格不等号/半开语义的工程化处理（必须写进 baseline）

如果边界语义处理错（比如把 “$<$” 实现成 “$\le$”），你采样的就不是“真实的 $J$”了，直接破坏定义。

推荐两种方案（**优先方案 A**）：

#### 方案 A：坐标压缩（rank） + 半开上界（最干净）

对每个维度分别压缩坐标，使查询天然用半开上界表达严格不等号。

- 对前半部分维度 $i\in[1..d]$：
   收集所有 inner 的 $L_i(b)$ 与 outer 的 $R_i(a)$，排序去重得数组 `V_i`。
- 对后半部分维度 $i\in[1..d]$：
   收集所有 inner 的 $-R_i(b)$ 与 outer 的 $-L_i(a)$，排序去重得数组 `W_i`。

定义 `rank(x)` 为其在对应数组中的下标（或稳定映射）。
 对查询上界 `upper`，用 `hi = lower_bound(array, upper)`，则 “$x” 等价于 “rank(x) ∈ [0, hi)”。

这样你在嵌入空间中的查询区间全部是标准的半开 MBR：
$$
[0, hi_1)\times\cdots\times[0, hi_{2d}).
$$

#### 方案 B：浮点转整数 + 单位化边界

若输入为浮点且你不想做 rank：

- 统一放大 $10^k$ 并取整，转成 int
- 用整数严格不等号实现：
  - $x 等价于 $x\le y-1$
  - $x>y$ 等价于 $x\ge y+1$

> baseline 建议明确写：**默认使用方案 A**，因为更稳定、无 epsilon、实现更一致。

------

### 3.4 SIRSIndex：作为“range i.i.d. sampling 黑盒”，但要补齐 Count/Report

我们在点集 $P$ 上建立 `SIRSIndex`，其核心承诺来自 SIRS Definition 1：

- `Sample(Q, k)`：返回 $k$ 个来自 $P\cap Q$ 的 i.i.d. uniform 样本，并且跨 query 独立。

为了让 join baseline **严格无偏**，我们额外要求：

- `Count(Q)`：精确返回 $|P\cap Q|$
- `Report(Q)`：枚举 $P\cap Q$（Enumerate/Adaptive 分支需要）

> 重要澄清：SIRS 论文主目标是“采样”，baseline 要求 Count/Report 是为了 join 的外层权重 $\deg(a)$ 必须精确，否则就会引入系统性偏差。

下面给一个**可落地的实现建议**：
 选择 SIRS 的 **KD-tree instantiation (KDS)** 作为默认（论文 Sec. 3.3）。

------

### 3.5 推荐的可复现实现：KD-tree SIRSIndex（KDS 风格）+ 精确 Count/Report

> 你们论文里 baseline 可以写：
>  “我们用 SIRS’21 的 KD-tree sampling instantiation 建索引，并在同一棵 KD-tree 上实现精确 Count/Report。”

#### 3.5.1 结点结构（Node）

每个结点 `u` 存：

- `bbox_u`: 该节点子树点集的 MBR（在 $2d$ 维嵌入空间）
- `start_u, end_u`: 该子树点在存储数组中的连续区间（闭区间或半开均可，建议半开 `[start,end)`）
- `size_u = end_u - start_u`
- `left_u, right_u`（二叉 KD-tree）

并设置叶子容量 `B`（例如 256，SIRS 实验也常用类似叶子大小）。

> SIRS 的 KDS 核心要求是：每个节点对应一个连续存储区间，子节点区间被父节点区间覆盖。

#### 3.5.2 Build：构建 KD-tree 与连续布局（离线构建）

伪代码（示意）：

```
BuildKD(points array A, dim D=2d, leaf_size B):
    return BuildRec(A, l=0, r=|A|, depth=0)

BuildRec(A, l, r, depth):
    u = new Node()
    u.start = l; u.end = r
    u.bbox = MBR(A[l:r])
    if r-l <= B: u.is_leaf = true; return u
    split_dim = depth mod D
    mid = (l+r)/2  # median
    nth_element(A[l:r], mid, key=A[*].coord[split_dim])  # 使 A[l:mid) <= A[mid:r)
    u.left  = BuildRec(A, l, mid, depth+1)
    u.right = BuildRec(A, mid, r, depth+1)
    return u
```

> 这样建出来的树，天然满足“每个节点在数组上对应连续区间”的性质，从而可按 SIRS 的 Lemma 1 走“区间分解 + alias”的采样框架。

#### 3.5.3 QueryDecompose：把查询 MBR 分解为“完全包含结点 + 边界叶子”

我们需要一个统一的分解过程，供 Sample / Count / Report 复用：

```
Decompose(u, Q, full_nodes, boundary_leaves):
    if bbox_u ∩ Q == ∅: return
    if bbox_u ⊆ Q:
        full_nodes.append(u)
        return
    if u.is_leaf:
        boundary_leaves.append(u)
        return
    Decompose(u.left, Q, full_nodes, boundary_leaves)
    Decompose(u.right, Q, full_nodes, boundary_leaves)
```

- `full_nodes`：bbox 完全被 Q 包含的节点，节点区间 `[start,end)` 可直接用；
- `boundary_leaves`：与 Q 相交但不完全包含的叶子节点，需要扫描叶子点来做精确判断（用于精确 Count/Report，也可用于“无拒绝采样”的增强版 Sample）。

#### 3.5.4 Sample(Q,k)：SIRS 采样（黑盒语义 + 可复现实现轮廓）

SIRS 的一般框架是：

1. 用 `Decompose` 得到若干候选区间（来自 full_nodes，以及必要时的 boundary leaves）；
2. 对候选区间建立 top-level alias（权重是“该区间内属于 Q 的点数”）；
3. 重复 $k$ 次：alias 选区间 + 区间内均匀选点（必要时拒绝/或扫描后无拒绝）。

- SIRS 论文讨论了“边界叶子导致的 rejection”，以及“扫描边界叶子避免 rejection”的折中与开销。
- 对 baseline 而言，只要 `Sample` 满足 Definition 1（uniform + independence），我们即可把它当黑盒调用；但为了“可复现”，你可以在附录写上述轮廓与伪代码。

#### 3.5.5 Count(Q)：精确范围计数（baseline 必须补强的关键）

对 join 的正确性而言，`Count(Q(a))` 必须等于 $|P\cap Q(a)|$ 的**精确值**，不能是估计。

实现：

```
Count(Q):
    full_nodes, boundary_leaves = Decompose(root, Q)
    cnt = sum(u.size for u in full_nodes)  # full nodes 全部在 Q 内，直接加
    for leaf in boundary_leaves:
        for p in A[leaf.start:leaf.end]:
            if p.coord ∈ Q: cnt += 1
    return cnt
```

- 为什么精确：
  - full_nodes 的 bbox ⊆ Q ⇒ 该子树所有点都在 Q 内；
  - boundary_leaves 必须逐点检查，避免把“bbox 相交但点不在 Q 内”的点误计。

> 这一步就是你们 baseline 之前“需要补强但没写到可落地”的点：**精确 Count 不是 SIRS Lemma 1 的采样复杂度能自动替代的，而是一个明确的范围计数过程**。

#### 3.5.6 Report(Q)：精确范围报告（Enumerate/Adaptive 分支需要）

实现同理：

```
Report(Q):
    full_nodes, boundary_leaves = Decompose(root, Q)
    out = []
    for u in full_nodes:
        for p in A[u.start:u.end]:
            out.append(p)
    for leaf in boundary_leaves:
        for p in A[leaf.start:leaf.end]:
            if p.coord ∈ Q: out.append(p)
    return out
```

若考虑内存与流式输出，可提供 `ReportStream(Q)` 迭代器：逐点 yield。

------

### 3.6 AliasOuter：按 $\deg(a)$ 加权选择 outer

构建离散分布：
$$
\Pr(\text{AliasOuter.sample()}=a)=\frac{\deg(a)}{|J|},
$$
只对 $\deg(a)>0$ 的 $a$ 建表。

alias 的理论性质（线性建表、常数采样、独立）来自 SIRS Background 对 Walker alias 的回顾。

------

### 3.7 Enumerate/Adaptive 的结果容器

- `AllPairs`: 动态数组（vector）存储所有 join 对 $(a.id,b.id)$
- `sizeof(Pair)`：建议在报告里写出工程估计（例如两 64-bit id 最少 16B，但 vector/对齐常会到 16~24B+）
- `deg[a]`: 存储每个 outer 的精确度数（至少 64-bit）

------

## 4. 算法详细流程（三个版本）

为方便表述，定义：

- `OuterSet`, `InnerSet`：执行 baseline 时选定的 outer/inner（可交换，见 4.0）
- $P=\{\phi(b)\mid b\in \text{InnerSet}\}$
- `BuildQueryMBR(a)`：返回 $Q(a)$

### 4.0 所有版本共享的预处理

输入：$R_c, R_{\bar c}, t$

#### Step 0：可选 outer/inner 交换（只影响性能，不影响正确性）

可用简单 heuristic：

- 若 $n_1\gg n_2$，选更小的一侧做 outer，以减少 `Count(Q(a))` 的次数
   交换后要保证输出仍是 $(r_c,r_{\bar c})$ 的顺序。

#### Step 1：构建 inner 点集 $P$

对每个 inner 盒子 $b$：

- 计算嵌入点 $\phi(b)$
- 存入 `EmbeddedPoint(coord=phi(b), inner_id=b.id)`

#### Step 2：严格不等号处理（推荐 rank）

按 3.3 的方案 A 做 rank 压缩，得到整数域的点与查询。

#### Step 3：构建 `SIRSIndex`（KD-tree instantiation）

`SIRSIndex.Build(P)`
 SIRS 论文中 KD-tree instantiation 的整体构造与采样框架在 Sec. 3.3。

------

### 4.1 版本 A：Sampling（SIRS‑JS‑Sampling）

#### 4.1.1 思路

- Phase A：对每个 outer $a$ 做精确 `Count(Q(a))` 得 $\deg(a)$，建 `AliasOuter`
- Phase B：重复 $t$ 次：
   outer 用 alias 抽 $a$，inner 用 `SIRSIndex.Sample(Q(a),1)` 抽一个 $b$，输出 $(a,b)$

#### 4.1.2 伪代码（基础版：每次 Sample 1 个）

```
Algorithm SIRS-JS-Sampling(OuterSet, InnerSet, t):

# Preprocess: build P={phi(b)}, rank compress, build SIRSIndex

# Phase A: exact degrees
for each a in OuterSet:
    Qa = BuildQueryMBR(a)
    deg[a] = SIRSIndex.Count(Qa)      # MUST be exact

J = sum_a deg[a]
if J == 0: return []

AliasOuter = BuildAlias({a:deg[a] | deg[a] > 0})

# Phase B: i.i.d. samples
Ans = []
for j in 1..t:
    a = AliasOuter.sample()
    Qa = BuildQueryMBR(a)
    p = SIRSIndex.Sample(Qa, 1)[1]    # 1 uniform point in P∩Qa
    b = p.inner_id
    Ans.append(OutputAsCrossSetPair(a.id, b))

return Ans
```

#### 4.1.3 可选但推荐的工程优化：按 outer 分组批量 Sample

为了减少 “每次 Sample(Q,1) 都要构建 top-level alias” 的开销（SIRS Lemma 1 里的 $M+I$ 部分），可以：

1. 先抽 $t$ 次 outer：得到序列 $a_1,\dots,a_t$
2. 统计每个 outer 的出现次数 $k_a$
3. 对每个 $k_a>0$ 的 outer，只调用一次 `Sample(Q(a), k_a)` 返回 $k_a$ 个 i.i.d. 样本

由于 SIRS Definition 1 本身就要求同一 query 内返回 $k$ 个样本相互独立，并且跨 query 独立性也被强调，因此这种“批量调用”不改变分布，只提升性能。

------

### 4.2 版本 B：Enumerate+Sampling（SIRS‑JS‑Enumerate）

#### 4.2.1 思路

先枚举全部 $J$ 到数组 `AllPairs`，再对数组做独立均匀下标采样。

#### 4.2.2 伪代码

```
Algorithm SIRS-JS-Enumerate(OuterSet, InnerSet, t):

# Preprocess: build P, rank compress, build SIRSIndex

AllPairs = []
for each a in OuterSet:
    Qa = BuildQueryMBR(a)
    List = SIRSIndex.Report(Qa)     # exact report of P∩Qa
    for p in List:
        b = p.inner_id
        AllPairs.append(OutputAsCrossSetPair(a.id, b))

W = |AllPairs|
if W == 0: return []

Ans = []
for j in 1..t:
    idx = UniformInt(0, W-1)
    Ans.append(AllPairs[idx])

return Ans
```

------

### 4.3 版本 C：Adaptive+Sampling（SIRS‑JS‑Adaptive）

#### 4.3.1 设计目标

给定阈值 $J_\star$：

- 若 $|J|\le J_\star$：走 Enumerate+Sampling（快）
- 若 $|J|> J_\star$：走 Sampling（省内存/省枚举）

#### 4.3.2 无偏切换原则（必须强调）

1. **不能混合**“部分枚举的 AllPairs”与“后续 Sampling”共同输出
    → 一旦决定走 Sampling，必须丢弃 AllPairs
2. Sampling 分支必须有所有 outer 的**精确** $\deg(a)$ 以建 alias

#### 4.3.3 推荐实现：流式枚举尝试 + 越界时用 Count 补齐当前 outer

关键优化：越界发生时，不要继续枚举剩余 pair（否则在 $|J|$ 巨大时仍然可能被拖成 $\Omega(|J|)$）。做法是：

- 当处理某个 outer $a$ 的 Report 流时首次超过 $J_\star$：
  - `partial = deg_partial(a)`（已枚举的命中数）
  - `total = Count(Q(a))`（精确总数）
  - 用 `total - partial` 补齐累计计数
  - 立即停止 Report 当前 outer 的剩余输出
  - 清空 AllPairs，切换到 COUNT_ONLY

#### 4.3.4 伪代码

```
Algorithm SIRS-JS-Adaptive(OuterSet, InnerSet, t, J_star):

# Preprocess: build P, rank compress, build SIRSIndex

mode = ENUMERATE
AllPairs = []
W = 0
for each a in OuterSet: deg[a] = 0

for each a in OuterSet:
    Qa = BuildQueryMBR(a)

    if mode == ENUMERATE:
        for p in SIRSIndex.ReportStream(Qa):
            deg[a] += 1
            W += 1
            if W <= J_star:
                AllPairs.append(OutputAsCrossSetPair(a.id, p.inner_id))
            else:
                partial = deg[a]
                total = SIRSIndex.Count(Qa)        # exact total
                deg[a] = total
                W += (total - partial)             # now W includes full deg[a]
                AllPairs.clear()
                mode = COUNT_ONLY
                break

    else:  # COUNT_ONLY
        deg[a] = SIRSIndex.Count(Qa)               # exact
        W += deg[a]

if W == 0: return []

if mode == ENUMERATE:
    return UniformArraySample(AllPairs, t)

AliasOuter = BuildAlias({a:deg[a] | deg[a] > 0})
return OuterAliasPlusSIRS(AliasOuter, SIRSIndex, t)
```

------

## 5. Adaptive 阈值 $J_\star$ 的选择策略（可写进论文 + 可落地）

建议写成：
$$
J_\star=\min\Big(J_\star^{\text{mem}},\ J_\star^{\text{time}}\Big).
$$

### 5.1 内存硬约束（必须满足，否则 OOM）

设：

- 可用内存预算 `MemBudget`（字节）
- 允许给 `AllPairs` 的比例 $\rho\in(0,1)$（建议 0.2～0.5）
- pair 存储开销 `sizeof(Pair)`（至少 16B，两 64-bit id）

则：
$$
J_\star^{\text{mem}}
=
\left\lfloor
\frac{\rho\cdot \text{MemBudget}}{\text{sizeof(Pair)}}
\right\rfloor.
$$

### 5.2 时间权衡（让 Adaptive “值回票价”）

核心比较：

- 枚举版至少 $\Omega(|J|)$ 输出成本
- 采样版主要随 $t$ 增长（每个样本一次或批量几次 SIRS sampling）

推荐两种写法：

#### 写法 A（最审稿友好）：microbenchmark 标定交叉点

抽若干 query 做标定：

- 单次 `Sample(Q,1)` 平均耗时 $\hat c_{\text{samp}}$
- 枚举并写入 1 条 pair 的摊销耗时 $\hat c_{\text{pair}}$

估计交叉点：
$$
|J_{\text{cross}}|\cdot \hat c_{\text{pair}}
\approx
t\cdot \hat c_{\text{samp}}.
$$
取保守系数 $\gamma\in[0.7,0.9]$：
$$
J_\star^{\text{time}}=\lfloor \gamma\cdot |J_{\text{cross}}|\rfloor.
$$

#### 写法 B（不做 benchmark 的可解释线性阈值）

$$
J_\star^{\text{time}} = C_1\cdot t + C_2\cdot n_{\text{out}},
$$

其中 $C_1,C_2$ 是经验常数（在实验部分可固定）。

------

## 6. 算法分析（正确性 + 复杂度，三版本都给）

记：

- $n_{\text{out}}=|\text{OuterSet}|$，$n_{\text{in}}=|\text{InnerSet}|$
- $P=\{\phi(b)\}$
- $\deg(a)=|P\cap Q(a)|$
- $|J|=\sum_a \deg(a)$

### 6.1 正确性

#### 6.1.1 引理：相交 $\Leftrightarrow$ 点落入 $Q(a)$

已在 3.2.3 给出逐维等价证明，因此：

- `Count(Q(a)) = deg(a)`（精确）
- `Report(Q(a))` 枚举的 inner 集合等于所有与 $a$ 相交的 $b$
- `Sample(Q(a),k)` 在这些 $b$ 上均匀采样（SIRS 定义保证）

#### 6.1.2 Sampling 版本（SIRS‑JS‑Sampling）的均匀性

对任意固定 $(a,b)\in J$：
$$
\Pr(\text{pick }a)=\frac{\deg(a)}{|J|},\qquad
\Pr(\text{pick }b\mid a)=\frac{1}{\deg(a)}.
$$
因此：
$$
\Pr(\text{output }(a,b))=\frac{1}{|J|}.
$$

#### 6.1.3 Sampling 版本的独立性（i.i.d.）

- `AliasOuter.sample()` 每次调用独立（alias 方法要求不同 query 返回独立索引）
- `SIRSIndex.Sample(Q,1)` 返回独立样本，且 SIRS 强调跨 query 也独立 

因此 $Z_1,\dots,Z_t$ i.i.d. 且均匀。

#### 6.1.4 Enumerate+Sampling 的正确性

- 枚举阶段：对每个 $a$，Report(Q(a)) 得到所有相交 $b$，因此 `AllPairs` 不重不漏地包含 $J$
- 采样阶段：对数组下标独立均匀采样，自然得到 i.i.d. uniform with replacement

#### 6.1.5 Adaptive 的正确性

- 若未越界：等价于 Enumerate+Sampling
- 若越界：丢弃 AllPairs，依赖精确 $\deg(a)$ 建 alias，再走 Sampling
   切换只影响“是否保存枚举结果”，不会改变概率分布，因此无偏。

------

### 6.2 复杂度分析（建议论文中使用 SIRS Lemma 1 的抽象符号）

#### 6.2.1 SIRS 的抽象采样复杂度符号（直接引用 Lemma 1）

记：
$$
\text{SIRSCost}(R,k)=O\big(M(n,R)+I(n,R)+k\big)
$$
其中 $M,I$ 来自 SIRS Lemma 1。

并定义 baseline 额外需要的：

- $\text{CountCost}(R)$：精确 Count 的成本
- $\text{ReportCost}(R)$：精确 Report 的成本（至少输出量 $|P\cap R|$）

> 若你采用我在 3.5 的 KD-tree 实现，Count/Report 的过程是标准范围查询：全包含结点直接加子树规模；边界叶逐点检查。其时间由“访问节点数 + 边界叶扫描点数”构成，叶容量固定时可视作 $O(\#\text{visited nodes})$ 的主导项。

#### 6.2.2 Sampling（SIRS‑JS‑Sampling）

**时间**
$$
T_{\text{samp}}
=
\text{Build\_SIRS}(n_{\text{in}})
+
\sum_{a\in \text{Outer}} \text{CountCost}(Q(a))
+
O(n_{\text{out}})
+
\sum_{j=1}^{t} \text{SIRSCost}(Q(a_j),1).
$$
（若使用“按 outer 分组批量 Sample”，则最后一项替换为 $\sum_{a: k_a>0}\text{SIRSCost}(Q(a),k_a)$，并且 $\sum_a k_a=t$。）

**空间**
$$
S_{\text{samp}}=O(n_{\text{in}}+n_{\text{out}}+t).
$$

#### 6.2.3 Enumerate+Sampling（SIRS‑JS‑Enumerate）

**时间**
$$
T_{\text{enum}}
=
\text{Build\_SIRS}(n_{\text{in}})
+
\sum_a \text{ReportCost}(Q(a))
+
O(t),
$$
且 $\sum_a \text{ReportCost}(Q(a))$ 至少包含 $\Omega(|J|)$ 的输出成本。

**空间**
$$
S_{\text{enum}}
=
O(n_{\text{in}}+n_{\text{out}}+|J|+t).
$$

#### 6.2.4 Adaptive+Sampling（SIRS‑JS‑Adaptive）

- 若 $|J|\le J_\star$：时间/空间同 Enumerate，但 $|J|$ 被 $J_\star$ 截断
- 若 $|J|>J_\star$：最多枚举写入 $O(J_\star)$ 条后切换，之后同 Sampling

因此可写为：
$$
T_{\text{adp}}
=
\text{Build\_SIRS}(n_{\text{in}})
+
\sum_a \text{CountCost}(Q(a))
+
O(n_{\text{out}})
+
\sum_{j=1}^{t} \text{SIRSCost}(Q(a_j),1)
+
O(J_\star),
$$

------

## 你在论文里建议额外加的一段“适用范围说明”（避免 reviewer 卡你）

由于嵌入把原 $d$ 维盒子变成 $2d$ 维点，baseline 实际运行在 $2d$ 维空间；而 SIRS 论文固定展示 $d=2$，并讨论面向低维（文中提到 $d<10$，评估到 $d=7$）。
 因此建议在 baseline 小节加一句“我们采用 SIRS 的低维假设；当 $d$ 较大时 baseline 可能退化，这属于 baseline 局限而非正确性问题”。