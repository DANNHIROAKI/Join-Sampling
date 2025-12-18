# 1. 问题定义与分析 + 引用文章内容

## 1.1 Spatial Join Sampling 任务定义（与你们整体论文对齐）

在 $d\ge 2$ 维欧氏空间 $\mathbb{R}^d$ 中，给定两类**轴对齐半开盒（half-open boxes）**集合：
$$
R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\}.
$$
每个盒子（半开）定义为：
$$
r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\qquad L_i(r)<R_i(r).
$$
只关心**跨集合相交对**（cross-set intersection pairs）：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
输出 $t$ 个样本 $Z_1,\dots,Z_t\in J$，要求：

- **均匀性（uniform）**：$\Pr(Z_j=P)=1/|J|$，$\forall P\in J$
- **独立性（i.i.d. with replacement）**：$Z_1,\dots,Z_t$ 相互独立，允许重复

------

## 1.2 引用文章：SIRS’21 的核心定义与“我们复用的能力”

**文章信息（baseline 依赖的“黑盒能力”来源）**

- 标题：*Spatial Independent Range Sampling*
- 作者：Dong Xie, Jeff M. Phillips, Michael Matheny, Feifei Li
- 会议：SIGMOD 2021
- DOI：10.1145/3448016.3452806 [Pure at Penn State+1](https://pure.psu.edu/en/publications/spatial-independent-range-sampling?utm_source=chatgpt.com)

**SIRS 定义（Uniform SIRS）**
 SIRS 将问题形式化为：给定点集 $P\subset\mathbb{R}^d$、查询 MBR $R$、整数 $k$，返回 $k$ 个来自 $R\cap P$ 的**独立**均匀样本，每个点的概率为 $1/|R\cap P|$。这对应论文 Definition 1。 [users.cs.utah.edu](https://users.cs.utah.edu/~jeffp/papers/sirs-sigmod.pdf)

**SIRS 对独立性的关键强调：跨 query 也必须独立**
 你草稿里写的“跨 query 独立性”并不是你们额外加的假设——SIRS 论文明确写道独立性不仅针对同一 query 内样本，也应对不同 queries 返回的样本成立。 [users.cs.utah.edu](https://users.cs.utah.edu/~jeffp/papers/sirs-sigmod.pdf)
 （这点对 join sampling 非常关键，因为每个样本可能来自不同的 outer 盒子 $a$，即不同的 $Q(a)$。）

**SIRS 的通用框架 Lemma 1（复杂度抽象）**
 只要能把查询 MBR $R$ 映射成 $I(n,R)$ 个“连续索引区间”，且映射耗时 $M(n,R)$，就能在
$$
O(M(n,R)+I(n,R)+k)
$$
时间内返回 $k$ 个独立均匀样本。 [users.cs.utah.edu](https://users.cs.utah.edu/~jeffp/papers/sirs-sigmod.pdf)
 这也是你们 baseline 分析里最“可引用且审稿友好”的复杂度刻画方式。

**SIRS 中的 Walker alias（离散按权采样的 O(1) 查询组件）**
 SIRS 论文在 Background 里专门回顾了 Walker alias 方法：线性预处理 + $O(1)$ 抽样，用于按权重选择离散对象（例如区间、节点等）。 [users.cs.utah.edu](https://users.cs.utah.edu/~jeffp/papers/sirs-sigmod.pdf)
 我们在 join sampling 外层“按 $deg(a)$ 加权选 $a$”也正好复用这一组件。

> 一句话总结：SIRS 给我们一个可以当黑盒调用的能力：
>  **对点集 $P$ 的任意 MBR 查询 $R$，返回 $k$ 个独立均匀样本（且跨 query 独立）。** [users.cs.utah.edu+1](https://users.cs.utah.edu/~jeffp/papers/sirs-sigmod.pdf)

------

## 1.3 核心对齐：把“盒子相交 join”转成“点集 MBR 查询 + range sampling”

SIRS 的输入是点集，但 join 的 inner 侧是盒子。我们采用标准等价嵌入，把 **“盒子相交”** 转成 **“点落入 MBR”**：

### 1.3.1 半开盒相交的严格条件（必须说清楚）

对半开盒 $a,b$，有：
$$
a\cap b\neq\varnothing
\quad\Longleftrightarrow\quad
\forall i\in[d]:
\big(L_i(a)<R_i(b)\big)\ \wedge\ \big(L_i(b)<R_i(a)\big).
$$
等价写法（你草稿里用的形式）：
$$
\forall i:
\big(L_i(b)<R_i(a)\big)\ \wedge\ \big(R_i(b)>L_i(a)\big).
$$
两种写法完全等价（都是严格不等号），并且与半开语义一致：**贴边（如 $R_i(b)=L_i(a)$）不算相交**。

### 1.3.2 2d 维嵌入：inner 盒子 → 点

对每个 inner 盒子 $b\in R_{\bar c}$，定义点：
$$
\phi(b)=\big(L_1(b),\dots,L_d(b),\ -R_1(b),\dots,-R_d(b)\big)\in\mathbb{R}^{2d}.
$$

### 1.3.3 outer 盒子 → 查询 MBR $Q(a)$

对 outer 盒子 $a\in R_c$，构造 $2d$ 维查询 MBR：
$$
Q(a)= \prod_{i=1}^d [\underline{L}_i,\ R_i(a)) \ \times\ \prod_{i=1}^d [\underline{U}_i,\ -L_i(a)),
$$
其中下界常量只需保证“足够小”覆盖所有 $\phi(b)$ 的坐标，例如：
$$
\underline{L}_i := \min_{b\in R_{\bar c}} L_i(b),\qquad
\underline{U}_i := \min_{b\in R_{\bar c}} (-R_i(b)) = -\max_{b\in R_{\bar c}} R_i(b).
$$

### 1.3.4 等价性（join ↔ range query）：正确性关键性质

对任意 $a\in R_c, b\in R_{\bar c}$：
$$
a\cap b\neq\varnothing
\quad\Longleftrightarrow\quad
\phi(b)\in Q(a).
$$
理由：

- $\phi(b)$ 前 $d$ 维落入 $[-\infty, R_i(a))$ 等价于 $L_i(b)
- $\phi(b)$ 后 $d$ 维落入 $[-\infty, -L_i(a))$ 等价于 $-R_i(b)<-L_i(a)\iff R_i(b)>L_i(a)$

因此：

- $deg(a):=|\{b:(a,b)\in J\}| = |P\cap Q(a)|$，其中 $P=\{\phi(b)\}$
- 在 $P\cap Q(a)$ 上均匀采一个点，就等价于在与 $a$ 相交的 $b$ 上均匀采一个

------

# 2. 核心数据结构

本节把 baseline 落地所需结构写成“可复现 API”，避免停留在概念层。

## 2.1 数据对象与基本映射

- **Box**
  - `id`：全局唯一
  - `L[1..d], R[1..d]`：半开盒坐标
- **EmbeddedPoint**（inner 盒子的 2d 维嵌入点）
  - `coord[1..2d]`：即 $\phi(b)$
  - `inner_id`：指回原盒子 $b$ 的 id
- 映射表：
  - `inner_id -> Box b`
  - `outer_id -> Box a`

------

## 2.2 严格不等号/半开语义的工程化处理（必须写进可复现细节）

你的草稿已经意识到这点，我这里把它提升为“baseline 的硬性实现要求”：

> 必须严格实现
>  $\;L_i(b)<R_i(a)$ 与 $\;R_i(b)>L_i(a)\;$
>  否则会把“贴边不相交”的 pair 错当成相交，进而 **改变 $J$**，破坏“uniform over true $J$”的定义。

推荐两种等价方案（二选一）：

### 方案 A：坐标压缩（rank）+ 半开区间（推荐，最干净）

对每个维度分别构建排序数组：

- 对前半部分维度 $i\in[1..d]$：压缩所有 $L_i(b)$ 与所有 $R_i(a)$
- 对后半部分维度 $i\in[1..d]$：压缩所有 $-R_i(b)$ 与所有 $-L_i(a)$

实现时：

- 点坐标用 `rank(x)` 映射到整数
- 查询上界用 `lower_bound` 找到第一个 $\ge$ 上界的位置作为 `hi`
   则 “< upper” 等价于 rank in `[0, hi)`

这样天然保留严格不等号（因为上界是半开）。

### 方案 B：浮点转整数 + 边界单位化

若输入为浮点：

- 乘 $10^k$ 并四舍五入成 int
- 用整数域表示严格不等号
  - “$x” 可实现为 “$x\le y-1$”
  - “$x>y$” 可实现为 “$x\ge y+1$”

> 方案 A 更推荐：不需要拍 epsilon，行为最稳定。

------

## 2.3 SIRSIndex：我们把 SIRS 当“可调用黑盒”

建立在 $P=\{\phi(b)\}$ 上的结构 `SIRSIndex`，提供：

- `Build(P)`：构建 SIRS 索引（SIRS 论文给出多种 instantiation，如 KD-tree / z-order 等；论文固定 $d=2$ 讲解，并预计适用于低维 $d<10$ 的场景，评估到 $d=7$。 [users.cs.utah.edu](https://users.cs.utah.edu/~jeffp/papers/sirs-sigmod.pdf)）
- `Sample(Q, k)`：返回 $k$ 个来自 $P\cap Q$ 的**i.i.d. uniform**样本（with replacement），并且跨 query 仍保持独立性（SIRS 要求）。 [users.cs.utah.edu+1](https://users.cs.utah.edu/~jeffp/papers/sirs-sigmod.pdf)

为了让 join baseline 自洽，我们还需要 **精确度数**：

- `Count(Q)`：精确返回 $|P\cap Q|$
- （Enumerate/Adaptive 分支需要）`Report(Q)`：枚举所有 $P\cap Q$

> 重要：SIRS 本身解决的是“采样”。你们 baseline 额外要求 Count/Report，是为了保证 join 上的严格均匀性（外层必须按精确 $deg(a)$ 加权）。
>  这可以通过**同一棵 KD-tree / R-tree**实现范围计数与报告（bbox containment → 加子树规模；边界 → 递归或逐点检查）。

------

## 2.4 Outer alias：按 join 度加权选择 outer 盒子

定义 outer 侧每个盒子 $a$ 的 join 度：
$$
deg(a):=|P\cap Q(a)|.
$$
则：
$$
|J|=\sum_{a\in R_c}deg(a).
$$
构建 Walker alias 结构 `AliasOuter`，只对 $deg(a)>0$ 的 $a$ 建表：
$$
\Pr(\text{AliasOuter.sample()}=a)=\frac{deg(a)}{|J|}.
$$
alias 的 $O(1)$ 抽样、$O(n)$ 建表特性是 SIRS 论文 Background 中明确回顾的组件。 [users.cs.utah.edu](https://users.cs.utah.edu/~jeffp/papers/sirs-sigmod.pdf)

------

## 2.5 Enumerate/Adaptive 用的结果容器

- `AllPairs`：动态数组，存储所有 join pair $(a.id, b.id)$
- `sizeof(Pair)`：通常是两个 32/64-bit id（工程上建议明确）
- 枚举版的空间是 $\Theta(|J|)$，Adaptive 用阈值 $J_\star$ 控制

------

# 3. 算法详细流程（三个版本）

## 3.0 共同预处理（所有版本共享）

输入：$R_c, R_{\bar c}, t$

### Step 0：可选的 outer/inner 交换（只影响性能，不影响正确性）

默认 outer=$R_c$，inner=$R_{\bar c}$。你可以加入一个简单 heuristic：

- 若 $n_1\gg n_2$，把 outer 设为更小的一侧，以减少要做多少次 `Count(Q(a))`
- 交换后，所有定义保持一致：outer 变成 $R_{\bar c}$，inner 变成 $R_{c}$，只要统一使用 $\phi(\cdot)$ 嵌入 inner 即可

**正确性不受影响**：因为 join 是对称关系，外层只是“抽样的组织方式”。

### Step 1：构建 inner 点集 $P$

对每个 $b\in$ inner：

- 计算嵌入点 $\phi(b)$
- 写入 `EmbeddedPoint(coord=phi(b), inner_id=b.id)` 到 $P$

### Step 2：准备查询构造所需常量（可用 rank 压缩替代）

计算每维的 $\underline{L}_i, \underline{U}_i$（或使用 rank 后下界就是 0）。

### Step 3：构建 SIRSIndex

`SIRSIndex.Build(P)`
 （你可以在实现里选择 SIRS 的 KD-tree instantiation 或其它；baseline 文稿只需要把它当满足 Definition 1 的黑盒即可。 [users.cs.utah.edu+1](https://users.cs.utah.edu/~jeffp/papers/sirs-sigmod.pdf)）

------

## 3.1 版本 A：Sampling（SIRS‑JS‑Sampling）

### 3.1.1 思想

- 外层：按 $deg(a)$ 加权选 $a$
- 内层：在 $P\cap Q(a)$ 上调用 SIRS 采 1 个点，对应 $b$
- 输出 $(a,b)$

### 3.1.2 详细流程

**Phase A：计算 join 度 & 建 alias**

1. 对每个 outer 盒子 $a$：
   - 构造 $Q(a)$
   - $deg(a)\leftarrow \texttt{SIRSIndex.Count}(Q(a))$（必须精确）
2. 计算 $|J|=\sum_a deg(a)$。若 $|J|=0$ 返回空。
3. 建 alias：`AliasOuter.Build({a:deg(a) | deg(a)>0})`

**Phase B：采样 $t$ 次**

重复 $t$ 次：

1. `a ← AliasOuter.sample()`
2. 构造 $Q(a)$
3. `p ← SIRSIndex.Sample(Q(a), 1)`（返回一个点）
4. `b ← p.inner_id`
5. 输出 pair $(a.id, b)$

### 3.1.3 伪代码

```
Algorithm SIRS-JS-Sampling(OuterSet, InnerSet, t):

# Preprocess: build P={phi(b)} and SIRSIndex on P

# Phase A: exact degrees
for each a in OuterSet:
    Qa = BuildQueryMBR(a)
    deg[a] = SIRSIndex.Count(Qa)   # exact

J = sum_a deg[a]
if J == 0: return []

AliasOuter = BuildAlias({a: deg[a] | deg[a] > 0})

# Phase B: i.i.d. samples
Ans = []
for j = 1..t:
    a = AliasOuter.sample()
    Qa = BuildQueryMBR(a)
    p = SIRSIndex.Sample(Qa, 1)[1]
    b = p.inner_id
    Ans.append((a.id, b))

return Ans
```

------

## 3.2 版本 B：Enumerate + Sampling（SIRS‑JS‑Enumerate）

### 3.2.1 思想

先把全部 $J$ 枚举进数组 `AllPairs`，再在数组上做 i.i.d. 均匀下标采样。

- 优点：当 $|J|$ 小且 $t$ 很大，采样阶段极快
- 缺点：空间 $\Theta(|J|)$

### 3.2.2 流程

1. `AllPairs ← []`
2. 对每个 outer 盒子 $a$：
   - 构造 $Q(a)$
   - `List ← SIRSIndex.Report(Q(a))`（枚举全部命中点）
   - 对每个点 $p\in List$：`AllPairs.append((a.id, p.inner_id))`
3. $W=|AllPairs|$。若 $W=0$ 返回空
4. 对 $j=1..t$：采样 `idx ~ Uniform{0..W-1}`，输出 `AllPairs[idx]`

### 3.2.3 伪代码

```
Algorithm SIRS-JS-Enumerate(OuterSet, InnerSet, t):

# Preprocess: build P and SIRSIndex

AllPairs = []
for each a in OuterSet:
    Qa = BuildQueryMBR(a)
    List = SIRSIndex.Report(Qa)
    for p in List:
        AllPairs.append((a.id, p.inner_id))

W = |AllPairs|
if W == 0: return []

Ans = []
for j = 1..t:
    idx = UniformInt(0, W-1)
    Ans.append(AllPairs[idx])

return Ans
```

------

## 3.3 版本 C：Adaptive + Sampling（SIRS‑JS‑Adaptive）

### 3.3.1 目标与无偏切换原则

给定阈值 $J_\star$：

- 若 $|J|\le J_\star$：走 Enumerate+Sampling
- 若 $|J|>J_\star$：走 Sampling

**无偏的关键原则**：

1. **不能混合**“部分枚举得到的 AllPairs”与“后续 Sampling”一起输出
    → 一旦决定走 Sampling，必须丢弃 `AllPairs`，否则会破坏均匀性（某些 pair 被“过度代表”）。
2. Sampling 分支必须拥有所有 $deg(a)$ 的精确值（用于 alias 权重）。

### 3.3.2 推荐实现：流式枚举尝试 + 越界后用 Count 补齐当前 outer 的剩余

你原稿的“继续扫描 stream 只计数不存储”在 $|J|$ 巨大时会退化为 $\Omega(|J|)$ 的无谓代价。更好的做法：

- 在某个 $a$ 的 `Report(Q(a))` 过程中第一次超过 $J_\star$ 时：
  1. 记录已经枚举的计数 `partial = deg_partial(a)`
  2. 立刻调用一次 `Count(Q(a))` 得到 `total = deg(a)`
  3. 把 $W$ 增加 `total - partial`（补齐真实总数）
  4. **停止该 outer 的后续 Report**（不再扫描剩余输出）
  5. 丢弃 AllPairs，切换到 Count-only 处理后续 outer

这样可保证：

- 额外枚举输出最多 $J_\star$ 条
- 越界后不会再被 $|J|$ 拖死

### 3.3.3 流程（单 pass + 早停）

初始化：

- `mode = ENUMERATE`
- `AllPairs=[]`
- `W=0`
- `deg[a]=0`

遍历每个 outer $a$：

- 若 `mode == ENUMERATE`：
  - 开始 `Report(Q(a))` 流式输出
  - 每产出一个命中点 $p$：
     `deg[a]++ ; W++`
    - 若 `W <= J_star`：存到 `AllPairs`
    - 若 `W == J_star + 1`（首次越界）：
      1. `partial = deg[a]`
      2. `total = Count(Q(a))`（精确）
      3. `deg[a] = total`
      4. `W += (total - partial)`（补齐）
      5. `AllPairs.clear()`，`mode=COUNT_ONLY`
      6. **break**（停止继续 Report 该 $a$）
- 若 `mode == COUNT_ONLY`：
  - `deg[a] = Count(Q(a))`
  - `W += deg[a]`

遍历结束：

- 若 `W==0` 返回空
- 若 `mode==ENUMERATE`：AllPairs 即完整 $J$，走数组采样
- 否则：构建 alias，走 Sampling 输出

### 3.3.4 伪代码

```
Algorithm SIRS-JS-Adaptive(OuterSet, InnerSet, t, J_star):

# Preprocess: build P and SIRSIndex

mode = ENUMERATE
AllPairs = []
W = 0
for each a in OuterSet:
    deg[a] = 0

for each a in OuterSet:
    Qa = BuildQueryMBR(a)

    if mode == ENUMERATE:
        # stream report
        for p in SIRSIndex.ReportStream(Qa):
            deg[a] += 1
            W += 1
            if W <= J_star:
                AllPairs.append((a.id, p.inner_id))
            else:
                # first overflow: switch
                partial = deg[a]
                total = SIRSIndex.Count(Qa)     # exact total
                deg[a] = total
                W += (total - partial)
                AllPairs.clear()
                mode = COUNT_ONLY
                break  # stop reporting remaining for this a

    else:
        cnt = SIRSIndex.Count(Qa)
        deg[a] = cnt
        W += cnt

if W == 0: return []

if mode == ENUMERATE:
    return UniformArraySample(AllPairs, t)

AliasOuter = BuildAlias({a: deg[a] | deg[a] > 0})
return OuterAliasPlusSIRS(AliasOuter, SIRSIndex, t)
```

------

# 4. Adaptive 版本阈值 $J_\star$ 的选择策略（能写进论文 + 能落地）

阈值要同时满足**内存硬约束**与**时间权衡**。推荐写成两部分，然后取最小值：
$$
J_\star=\min\Big(J_\star^{\text{mem}},\ J_\star^{\text{time}}\Big).
$$

------

## 4.1 内存硬约束（必须满足，否则枚举分支可能 OOM）

设：

- 可用内存预算：`MemBudget`（字节）
- 允许给 `AllPairs` 的比例：$\rho\in(0,1)$（建议 0.2～0.5）
- 每个 pair 的存储开销：`sizeof(Pair)`（两 id + 容器开销）

则：
$$
J_\star^{\text{mem}}
=
\left\lfloor
\frac{\rho\cdot \text{MemBudget}}{\text{sizeof(Pair)}}
\right\rfloor.
$$
工程上建议把 `sizeof(Pair)` 写成：

- 若仅存两 64-bit id：至少 16B，但向量/对象对齐后常常是 16~24B
- 如果还存指针/额外字段，请按真实结构体大小估

------

## 4.2 时间权衡（让 Adaptive “值回票价”）

核心比较：

- 枚举版开销至少 $\Omega(|J|)$（因为要输出所有 pair）
- 采样版开销与 $t$ 成正比（反复 `Sample(Q,1)`）

推荐一个可写入论文的、足够稳健的模型：

### 4.2.1 经验交叉点（最推荐：microbenchmark 标定）

做一次小规模标定：

1. 随机抽若干 outer $a$，测：
   - 单次 `Sample(Q(a),1)` 平均耗时 $\hat{c}_{\text{samp}}$
   - 生成一条 pair（即 report 输出一条命中并写入 AllPairs）的平均摊销耗时 $\hat{c}_{\text{pair}}$
2. 估计交叉点：

$$
|J_{\text{cross}}|\cdot \hat{c}_{\text{pair}}
\approx
t\cdot \hat{c}_{\text{samp}}.
$$

1. 设：

$$
J_\star^{\text{time}}=\left\lfloor \gamma\cdot |J_{\text{cross}}| \right\rfloor
$$

其中 $\gamma\in[0.7,0.9]$（保守一点更安全）。

### 4.2.2 若不做 benchmark：给一个“可解释”的线性阈值

$$
J_\star^{\text{time}} = C_1\cdot t + C_2\cdot n_{\text{outer}},
$$

其中：

- $C_1$：每个样本等价的枚举输出条数（经验常数）
- $C_2$：外层固定开销折算

这比拍一个固定阈值更审稿友好：你说明了“随 $t$ 增大，枚举更可能划算”。

------

# 5. 算法分析（正确性 + 复杂度：三版本都包含）

记：

- $n_{\text{out}}=|{\text{OuterSet}}|$，$n_{\text{in}}=|{\text{InnerSet}}|$
- $P=\{\phi(b)\mid b\in\text{InnerSet}\}$
- $deg(a)=|P\cap Q(a)|$
- $|J|=\sum_a deg(a)$

------

## 5.1 关键引理：相交 $\Leftrightarrow$ 点落入 $Q(a)$

对任意 $a$（outer）和 $b$（inner）：
$$
(a,b)\in J \iff a\cap b\neq\varnothing \iff \phi(b)\in Q(a).
$$
因此：

- `Count(Q(a))` 精确等于 $deg(a)$
- `Report(Q(a))` 枚举的 inner 集合精确等于 $\{b:(a,b)\in J\}$
- `Sample(Q(a),k)` 均匀采样点 $\phi(b)$ 等价于均匀采样匹配盒子 $b$

------

## 5.2 Sampling 版本的正确性

### 5.2.1 均匀性（单次样本）

任意固定 $(a,b)\in J$：

- 外层选中 $a$：

$$
\Pr(\text{pick }a)=\frac{deg(a)}{|J|}.
$$

- 条件在 $a$ 上，SIRS 在 $P\cap Q(a)$ 上均匀采样 1 个点：

$$
\Pr(\text{pick }b\mid a)=\frac{1}{deg(a)}.
$$

所以：
$$
\Pr(\text{output }(a,b))
=
\frac{deg(a)}{|J|}\cdot\frac{1}{deg(a)}
=
\frac{1}{|J|}.
$$

### 5.2.2 独立性（i.i.d.）

每次输出由两次独立随机过程组成：

1. `AliasOuter.sample()`（每次调用独立）
2. `SIRSIndex.Sample(Q,1)`（SIRS 定义要求样本独立，且跨 query 也独立） [users.cs.utah.edu+1](https://users.cs.utah.edu/~jeffp/papers/sirs-sigmod.pdf)

因此 $Z_1,\dots,Z_t$ i.i.d. 且均匀。

------

## 5.3 Enumerate+Sampling 版本的正确性

1. 枚举阶段：由 §5.1 引理，`AllPairs` 不重不漏地包含所有 $J$
2. 采样阶段：在长度 $|J|$ 的数组上做独立均匀下标采样 → 显然 i.i.d. uniform with replacement

------

## 5.4 Adaptive 版本的正确性

Adaptive 有两种输出模式：

- 若 $|J|\le J_\star$：等价于 Enumerate+Sampling → 正确
- 若 $|J|>J_\star$：等价于 Sampling → 正确

关键在于：分支选择只依赖确定的输入统计（最终的 $|J|$ 或越界判定），且**一旦越界丢弃 AllPairs**，不会引入任何“部分枚举导致的偏差”。

------

## 5.5 复杂度分析（用 SIRS Lemma 1 的抽象记号最审稿友好）

### 预处理

- 构建点集：$O(n_{\text{in}})$
- 建 SIRSIndex：记为 `Build_SIRS(n_in)`（SIRS 论文讨论了线性空间、并给出 KD-tree 等实现；这里 baseline 可以引用其结论/实现） [users.cs.utah.edu](https://users.cs.utah.edu/~jeffp/papers/sirs-sigmod.pdf)

------

### 5.5.1 Sampling（SIRS‑JS‑Sampling）

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
其中 `SIRSCost` 可直接引用 SIRS Lemma 1：$\text{SIRSCost}(R,k)=O\big(M(n,R)+I(n,R)+k\big)$

 **空间** $S_{\text{samp}} = O(n_{\text{in}} + n_{\text{out}} + t).$

------

### 5.5.2 Enumerate+Sampling（SIRS‑JS‑Enumerate）

**时间**
$$
T_{\text{enum}}
=
\text{Build\_SIRS}(n_{\text{in}})
+
\sum_{a} \text{ReportCost}(Q(a))
+
O(t).
$$
并且 $\sum_a \text{ReportCost}(Q(a))$ 至少包含 $\Omega(|J|)$ 的输出成本。

**空间**
$$
S_{\text{enum}} = O(n_{\text{in}} + n_{\text{out}} + |J| + t),
$$
其中 $|J|$ 是主瓶颈。

------

### 5.5.3 Adaptive+Sampling（SIRS‑JS‑Adaptive）

- 若未越界（$|J|\le J_\star$）：
  - 时间/空间同 Enumerate+Sampling，但 $|J|\le J_\star$
- 若越界（$|J|>J_\star$）：
  - 额外枚举写入最多 $J_\star$ 条（然后丢弃）
  - 之后转为 Sampling

因此可写为：
$$
T_{\text{adp}}
=
\text{Build\_SIRS}(n_{\text{in}})
+
\sum_{a} \text{CountCost}(Q(a))
+
O(n_{\text{out}})
+
\sum_{j=1}^{t} \text{SIRSCost}(Q(a_j),1)
+
O(J_\star).
$$
空间：
$$
S_{\text{adp}}
=
O(n_{\text{in}} + n_{\text{out}} + t + J_\star),
$$
并且越界后 `AllPairs` 被清空，实际峰值受 $J_\star$ 控制。