# Plane Sweep + Dynamic AABB Tree：Spatial Join 的 i.i.d. 均匀抽样（统一算法文档·增强版）

## 目录

1. 问题定义与分析
2. 核心数据结构
3. 算法详细流程（Sampling / Enumerate+Sampling / Adaptive+Sampling）
4. Adaptive 阈值 $J_\star$ 的选择策略
5. 算法分析（正确性与复杂度；三版本齐全）

------

## 1. 问题定义与分析

### 1.1 输入：两类 $d$ 维半开 AABB 与跨集合相交 Join

给定 $d\ge 2$（通常视为常数）维欧氏空间 $\mathbb{R}^d$ 中两类轴对齐半开盒集合
$$
R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\},\quad n=n_1+n_2.
$$
每个盒子（AABB）为半开形式：
$$
r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\qquad L_i(r)<R_i(r).
$$
我们只关心跨集合的相交对（spatial join result）：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
**半开相交判定（精确语义）**
 两个半开盒 $A=\prod_i[a_i,b_i)$、$B=\prod_i[c_i,d_i)$ 相交当且仅当逐维严格相交：
$$
A\cap B\neq\varnothing
\iff
\forall i,\ \max(a_i,c_i)<\min(b_i,d_i).
$$
等价写法：$\forall i, a_i<d_i \wedge c_i<b_i$，但实现时推荐用 $\max<\min$ 保持一致。

> 备注（工程实现非常重要）：使用浮点坐标时，请**统一比较策略**（不引入 epsilon 或者全程固定 epsilon 且保证一致性），否则“半开 + END-before-START”的边界语义可能被破坏。

------

### 1.2 输出：在 $J$ 上 i.i.d. 均匀（有放回）采样

目标输出 $t$ 个样本：
$$
Z_1,\dots,Z_t\in J,\quad Z_j=(r_c,r_{\bar c}),
$$
要求严格 **i.i.d. 均匀（with replacement）**：
$$
\Pr(Z_j=P)=\frac{1}{|J|}\ (\forall P\in J),\qquad
Z_1,\dots,Z_t\ \text{相互独立}.
$$
统一约定输出 pair 顺序始终为 $(r_c,r_{\bar c})$：若某次 START 来自 $R_{\bar c}$，则输出时交换顺序。

------

### 1.3 第 1 维 Plane Sweep：事件系统、活跃集与“唯一归属”分块

选择第 1 维 $x_1$ 扫描。每个盒子 $r$ 产生两个事件：

- $\texttt{START}(r)$ 发生在 $x_1=L_1(r)$
- $\texttt{END}(r)$ 发生在 $x_1=R_1(r)$

**事件排序必须是稳定全序**（用于唯一归属）：

1. 按 $x_1$ 升序；
2. 同 $x_1$：**END 在 START 之前**（匹配半开 $[L,R)$ 的“右端点不包含”语义）；
3. 同 $x_1$ 且同为 START：固定类优先级（例如先 $R_c$ 再 $R_{\bar c}$）；若还相同，可再按 `id` 升序，保证严格全序；
4. 同 $x_1$ 且同为 END：任意（可保持输入稳定性）。

扫描到位置 $x_0$ 时，定义两侧活跃集：
$$
\mathrm{Active}(R_c)=\{r_c\in R_c\mid x_0\in[L_1(r_c),R_1(r_c))\},
\quad
\mathrm{Active}(R_{\bar c})\ \text{同理}.
$$
对每个 START 事件 $e$，令新盒子 $q=r(e)$。定义局部集合：

- 若 $q\in R_c$：
  $$
  K_e=\{r_{\bar c}\in \mathrm{Active}(R_{\bar c})\mid q\cap r_{\bar c}\neq\varnothing\},\quad
  J_e=\{(q,r_{\bar c})\mid r_{\bar c}\in K_e\},\quad
  w_e=|K_e|.
  $$

- 若 $q\in R_{\bar c}$（对称）：
  $$
  K_e=\{r_{c}\in \mathrm{Active}(R_{c})\mid q\cap r_{c}\neq\varnothing\},\quad
  J_e=\{(r_{c},q)\mid r_{c}\in K_e\},\quad
  w_e=|K_e|.
  $$

**命题 1（事件分块的不交并分解）**
 设 $E$ 为所有 START 事件集合，则在“半开 + END-before-START + START tie-break”的排序规则下：
$$
J=\biguplus_{e\in E} J_e,\qquad |J|=\sum_{e\in E} w_e.
$$
**证明要点（建议直接放进正文）**
 对任意相交对 $(a,b)\in J$，考虑两者 START 事件在全序中的较晚者 $e^\star$。在处理 $e^\star$ 时：

- 较早 START 的那个盒子已开始且尚未 END（因为二者在第 1 维相交，且 END 在 START 前处理），故它属于对立活跃集；
- 因此该对必出现在 $J_{e^\star}$；
- 反过来，在更早 START 时刻，较晚 START 的盒子尚未活跃，所以该对不可能被计入更早的任何 $J_e$。

故不重不漏且归属唯一。

------

### 1.4 关键降维：处理 START($q$) 时只需看维度 $2..d$

定义投影到维度 $2..d$ 的 $(d-1)$ 维半开盒：
$$
q^\star=\prod_{i=2}^{d}[L_i(q),R_i(q)),\quad
r^\star=\prod_{i=2}^{d}[L_i(r),R_i(r)).
$$
**引理 2（投影相交等价）**
 当处理事件 $\texttt{START}(q)$ 且扫描位置 $x_0=L_1(q)$ 时，对于任意对立活跃盒 $r\in \mathrm{Active}(\text{Other})$：
$$
q\cap r\neq\varnothing\quad \Longleftrightarrow\quad q^\star\cap r^\star\neq\varnothing.
$$
**原因**：活跃性给出 $x_0\in[L_1(r),R_1(r))$ 且显然 $x_0\in[L_1(q),R_1(q))$，故第 1 维必相交；其余维度完全由投影决定。

> 因此：在扫线过程中，我们只需维护对立活跃集的投影集合 $\{r^\star\}$，支持
>
> - `CountIntersect(q^\star)` 精确得到 $w_e$
> - `SampleIntersect(q^\star,t')` 在 $K_e$ 上精确 i.i.d. 均匀采样（有放回）

这两项由动态 AABB Tree（动态 BVH）完成。

------

## 2. 核心数据结构

令 $m=d-1$。对任意盒 $r$，记其投影为
$$
r^\star=\prod_{k=1}^{m}[a_k(r),b_k(r)),\qquad a_k(r)=L_{k+1}(r),\ b_k(r)=R_{k+1}(r).
$$

### 2.1 事件数组与索引字段

事件记录结构（建议实现字段）：

- `x`：事件坐标 $x_1$
- `type`：`START` / `END`
- `cls`：`c` / `barc`
- `id`：盒子唯一 id
- `start_idx`：若为 START，则为该 START 在排序后得到的连续编号 $e\in\{1,\dots,|E|\}$；否则为空

> 关键：Phase1/Phase3 都需要从事件记录中**O(1)** 找到 $e$ 来读写 $w_e$ 与 slot。

------

### 2.2 两棵动态 AABB Tree（维护对立活跃投影）

维护两棵动态 BVH：

- $\mathcal{T}^{(c)}$：存当前活跃的 $\{r_c^\star\}$
- $\mathcal{T}^{(\bar c)}$：存当前活跃的 $\{r_{\bar c}^\star\}$

每个节点包含：

- `AABB(node)`：该子树所有叶条目的包围盒（可能是 **conservative**，例如 fat AABB）
- `size(node)`：子树叶数（用于 full-accept 与子树内均匀抽样）
- `left/right`：子指针（若实现为多叉，也可扩展到 children 列表）
- 叶节点存储：`(id, true_box = r^\star)`，可额外存 `fat_box` 用于插入/更新策略

#### 2.2.1 必须满足的 BVH 不变量（写清=加分关键点）

为保证 Count/Sample 正确性，BVH 只需满足以下 **最小正确性不变量**（与具体旋转/启发式无关）：

- **(I1) 包围正确（conservative）**：对任意节点 `node`，`AABB(node)` 是其所有后代叶 `true_box` 的超集（允许 fat 扩张，但不允许漏包含）。
- **(I2) size 正确**：`size(node)` 等于子树叶子条目数；插入/删除/旋转后必须维护。
- **(I3) 句柄可删除**：插入返回 `handle`，在对象仍活跃期间可用于 O(1) 定位并删除该叶（随后沿路径修复 AABB 与 size）。
- **(I4) 查询期间结构不变**：一次 `SampleIntersect(Q,t')` 期间，BVH 不发生插删（本算法天然满足：对某个 START 的查询发生在“插入 q 之前”。）

只要满足这些不变量，即使 BVH 最坏退化，**正确性仍严格成立**，差异仅体现在访问节点数和运行时间。

------

### 2.3 半开相交与包含判定（精确/安全）

#### 2.3.1 半开相交（叶子必须精确）

对两个 $m$ 维半开盒 $A=\prod_k[\alpha_k,\beta_k)$、$B=\prod_k[\gamma_k,\delta_k)$：
$$
\texttt{Intersect}(A,B)\iff \forall k:\ \max(\alpha_k,\gamma_k)<\min(\beta_k,\delta_k).
$$
**叶节点必须用 `true_box` 做上述精确判定**，用于严格遵守半开语义与防止贴边误计。

#### 2.3.2 半开包含（用于 full-accept，必须安全）

对节点包围盒 $B$ 与查询盒 $Q$：
$$
\texttt{Contained}(B,Q)\iff \forall k:\ L_k(B)\ge L_k(Q)\ \land\ R_k(B)\le R_k(Q).
$$
若成立，则整子树所有叶盒都 $\subseteq B\subseteq Q$，因此必与 $Q$ 相交（因为叶盒非空）。
 若 BVH 用 fat AABB：`Contained(fat, Q)` 仍安全（fat 在 Q 内 $\Rightarrow$ true 也在 Q 内）。

------

### 2.4 `CountIntersect(Q)`：精确计数

递归（或迭代栈）实现：

```
Count(node, Q):
  if node == null: return 0
  if not Intersect(AABB(node), Q): return 0
  if Contained(AABB(node), Q): return size(node)
  if node is leaf:
       return 1 if Intersect(node.true_box, Q) else 0
  return Count(node.left, Q) + Count(node.right, Q)
```

输出：
$$
\mathrm{CountIntersect}(Q)=|\{r^\star\in \mathrm{Active}:\ r^\star\cap Q\neq\varnothing\}|.
$$

> 说明：内部节点用 `AABB(node)` 做剪枝允许 false positive，但不能 false negative（由不变量 I1 保证），因此不会漏计。

------

### 2.5 `SampleIntersect(Q,t')`：在命中集合上严格均匀 i.i.d. 采样

AABB Tree 通常支持 report，但不保证“从 query 结果中均匀采样”。我们采用你文稿里的思路，并把关键引理与实现不变量补齐：**BuildGuide（命中计数结构）+ 按 cnt 加权随机下降**。

#### 2.5.1 `BuildGuide(Q)`：一次遍历生成“命中块”的不交分解

GuideNode 字段：

- `type ∈ {FULL, LEAF, PARTIAL}`
- `cnt`：该 guide 覆盖的命中叶数
- 若 `FULL`：保存指向 BVH 子树根 `node*`
- 若 `LEAF`：保存叶条目 `id`
- 若 `PARTIAL`：保存子 guide 列表（常见二叉即左右）

构造过程：

```
BuildGuide(node, Q):
  if node == null: return null (cnt=0)
  if not Intersect(AABB(node), Q): return null
  if Contained(AABB(node), Q):
        return Guide(FULL, cnt=size(node), ptr=node)
  if node is leaf:
        if Intersect(node.true_box, Q):
             return Guide(LEAF, cnt=1, id=node.id)
        else return null
  gL = BuildGuide(node.left, Q)
  gR = BuildGuide(node.right, Q)
  cnt = cnt(gL)+cnt(gR)   (treat null as 0)
  if cnt==0: return null
  return Guide(PARTIAL, cnt=cnt, children=[gL,gR] filtered non-null)
```

**引理 3（Guide 的不交分解性质）**
 `BuildGuide(Q)` 生成的 guide 结构把命中叶集合 $K(Q)$ 分解成若干不相交块（FULL 子树或命中 LEAF），并且根 `cnt` 恰为 $|K(Q)|$。
 证明可按递归结构：FULL 子树覆盖所有叶且与其他块不重；PARTIAL 由左右不交子树合并；剪枝只移除必不相交部分，不影响命中集合。

#### 2.5.2 `SampleSubtreeUniform(node)`：在子树叶上均匀抽 1 个

对 BVH 的任意节点 `node`（非 guide）：

```
SampleSubtreeUniform(node):
  while node is not leaf:
       L=node.left, R=node.right
       p = size(L)/size(node)
       node = L with prob p else R
  return node.id
```

**引理 4（子树均匀性）**
 若 (I2) `size` 正确，则 `SampleSubtreeUniform(node)` 对该子树每个叶条目输出概率均为 $1/\mathrm{size}(node)$。
 证明：对子树高度归纳即可。

> 若你的 BVH 是多叉节点：把“走子节点”改成按 `size(child)` 加权在 children 中采样即可，引理仍成立。

#### 2.5.3 `SampleOne(GuideRoot)`：对命中集合 $K(Q)$ 均匀抽 1 个

```
SampleOne(g):
  if g.type == FULL:
        return SampleSubtreeUniform(g.ptr)
  if g.type == LEAF:
        return g.id
  // PARTIAL
  choose child gi with prob cnt(gi)/cnt(g)
  return SampleOne(gi)
```

**引理 5（命中集合均匀性）**
 对任意命中叶条目 $x\in K(Q)$，`SampleOne` 返回 $x$ 的概率为 $1/|K(Q)|$。
 证明：对 guide 结构递归归纳——PARTIAL 先按块大小选块，再在块内均匀；FULL 块由引理 4 保证块内均匀。

#### 2.5.4 `SampleIntersect(Q,t')`：返回 $t'$ 个 i.i.d. 均匀样本（有放回）

```
SampleIntersect(Q, t'):
  g = BuildGuide(root, Q)
  W = g.cnt
  if W==0: return []
  for u in 1..t':
       out[u] = SampleOne(g)   // independent RNG draws
  return out
```

**引理 6（i.i.d.）**
 若每次 `SampleOne` 使用独立随机源，则输出序列为在 $K(Q)$ 上 i.i.d. 均匀（有放回）。

------

### 2.6 `IntersectReport(Q)`：枚举命中（用于 Baseline 或 Adaptive 小分支）

```
Report(node, Q, out):
  if node==null: return
  if not Intersect(AABB(node), Q): return
  if node is leaf:
        if Intersect(node.true_box, Q): out.append(node.id)
        return
  Report(node.left, Q, out)
  Report(node.right, Q, out)
```

> 说明：如果你的 BVH 实现内部就会对叶做精确判定，则 out 直接是 exact；
>  如果实现返回候选，需要在外部再做一次叶级精确过滤（但本增强版建议直接把精确过滤写入 report）。

------

### 2.7 Alias 表与 slot 容器（Sampling / Adaptive 大分支必需）

- **Alias 表**：在 START 事件集合 $E$ 上按权重 $w_e$ 构建分布
  $$
  p_e=\frac{w_e}{W},\quad W=\sum_{e\in E} w_e=|J|.
  $$
  采样 O(1)，构建 O(|E|)=O(n)。

- **零权重事件处理（必须写清）**：只对 $w_e>0$ 的 START 事件构建 alias；否则 alias 可能被零权重破坏（或引入无意义 slot）。

- **slot**：为每个 START 事件维护一个动态数组 $S_e\subseteq\{1,\dots,t\}$，表示哪些输出位置分配给该事件。Phase3 一次性采 $t_e=|S_e|$ 个并回填。

------

## 3. 算法详细流程（三个版本）

### 3.0 公共预处理（所有版本共享）

**输入**：$R_c, R_{\bar c}$，样本数 $t$，（Adaptive 额外给阈值 $J_\star$）。

1. 构造 `Events[1..2n]`：每个盒子生成 START/END 事件；
2. 按 1.3 的稳定全序排序；
3. 扫描 `Events`，为每个 START 事件赋 `start_idx = ++cnt_start`；
4. 为每个盒子准备 `handle` 存储位（用于删除）。

> 注意：Sampling/Adaptive 的 Phase3 需要“重新扫一遍”，所以 `Events` 必须可复用且排序规则固定。

------

### 3.1 版本一：Sampling（两遍扫描，不枚举）

**目标**：不显式存储 $J$，仍输出在 $J$ 上 i.i.d. 均匀样本。

#### Phase 1：第一次扫描（只计数得到 $w_e$ 与 $W=|J|$）

初始化：

- $\mathcal{T}^{(c)}\gets\emptyset,\ \mathcal{T}^{(\bar c)}\gets\emptyset$
- `w[1..|E|] = 0`
- `W = 0`

伪代码：

```
for event in Events:
  if event.type == END:
      if event.cls==c:      T^(c).Delete(handle[event.id])
      else:                 T^(barc).Delete(handle[event.id])

  else: // START
      e = event.start_idx
      q = box[event.id]
      Q = q^*   // projection to dims 2..d

      if q in R_c:
           w[e] = T^(barc).CountIntersect(Q)
           W += w[e]
           handle[q] = T^(c).Insert(q^*)
      else: // q in R_barc
           w[e] = T^(c).CountIntersect(Q)
           W += w[e]
           handle[q] = T^(barc).Insert(q^*)
```

若 `W==0` 则 $J=\varnothing$，直接返回空输出。

#### Phase 2：事件级 alias + slot 分配

1. 用所有 $w_e>0$ 的事件建立 alias（权重 $w_e$）。
2. 初始化 `S[e]` 为空（可用 vector/list）。

```
for j in 1..t:
    e = AliasSample()
    S[e].append(j)
```

令 $t_e=|S_e|$，有 $\sum_e t_e=t$。

#### Phase 3：第二次扫描（按 slot 批量局部采样并回填）

重置：

- 清空两棵树（活跃集为空）
- `Ans[1..t]`

```
reset T^(c), T^(barc)

for event in Events:
  if event.type == END:
      delete as in Phase1
  else:
      e = event.start_idx
      q = box[event.id]
      Q = q^*

      te = |S[e]|
      if te > 0:
          if q in R_c:
              List = T^(barc).SampleIntersect(Q, te) // ids in R_barc
              for u in 1..te:
                   j = S[e][u]
                   Ans[j] = (q, List[u])
          else:
              List = T^(c).SampleIntersect(Q, te) // ids in R_c
              for u in 1..te:
                   j = S[e][u]
                   Ans[j] = (List[u], q)          // keep (r_c, r_barc)

      insert q into its side tree (as in Phase1)
return Ans
```

------

### 3.2 版本二：Enumerate+Sampling（Baseline：枚举 $J$ 再数组采样）

**目标**：最直观基线；实现简单，但时间/空间依赖 $|J|$。

#### Step 1：一次扫描枚举所有 pair

初始化：

- 两棵树为空
- `Pairs = []`

```
for event in Events:
  if END: delete
  else: // START
      q = box[event.id]
      Q = q^*
      if q in R_c:
           ids = T^(barc).IntersectReport(Q)
           for id in ids: Pairs.append( (q, id) )
           handle[q] = T^(c).Insert(q^*)
      else:
           ids = T^(c).IntersectReport(Q)
           for id in ids: Pairs.append( (id, q) ) // keep order
           handle[q] = T^(barc).Insert(q^*)
W = |Pairs| = |J|
```

若 `W==0` 返回空；否则进入采样。

#### Step 2：在数组上独立均匀采样 $t$ 次

```
for j in 1..t:
    idx ~ UniformInt(0, W-1)
    Ans[j] = Pairs[idx]
return Ans
```

------

### 3.3 版本三：Adaptive+Sampling（阈值自适应）

**目标**：当 $|J|$ 小时采用 Baseline（省掉第二遍扫线）；当 $|J|$ 大时自动切到 Sampling（避免爆内存/爆枚举开销），并严格保持正确分布。

#### Phase 1：一次扫描（永远 Count；必要时 Report，但严格受阈值控制）

输入阈值 $J_\star$。

初始化：

- `mode = ENUMERATE`
- `AllPairs = []`（最多允许增长到 $J_\star$）
- 两棵树为空
- `w[e]` 数组
- `W=0`

伪代码（关键逻辑：先 Count，再判断阈值；且一旦超阈值，本事件不再枚举，并立即丢弃旧枚举结果）：

```
for event in Events:
  if END: delete

  else: // START
      e = event.start_idx
      q = box[event.id]
      Q = q^*

      // A: always count
      if q in R_c:      w[e]=T^(barc).CountIntersect(Q)
      else:             w[e]=T^(c).CountIntersect(Q)
      W += w[e]

      // B: optional enumeration, but only if still safe
      if mode==ENUMERATE and W <= J_star:
            if q in R_c:
                 ids = T^(barc).IntersectReport(Q)
                 for id in ids: AllPairs.append( (q,id) )
            else:
                 ids = T^(c).IntersectReport(Q)
                 for id in ids: AllPairs.append( (id,q) )
      else if mode==ENUMERATE and W > J_star:
            mode = COUNT_ONLY
            AllPairs.clear();   // release memory (optional shrink_to_fit)

      // C: insert q
      insert q into its tree
```

扫描结束后：得到精确 $W=|J|$ 与所有 $w_e$。

#### 分支 A：未切换（$W\le J_\star$）

此时 `AllPairs` 完整等于 $J$。在 `AllPairs` 上数组采样输出（同 baseline）。

#### 分支 B：已切换（$W>J_\star$）

直接执行 Sampling 的 Phase2+Phase3（事件 alias + slot + 第二遍扫线采样）。注意：无需再做 Phase1，因为计数 $w_e$ 已经保存。

------

## 4. Adaptive 阈值 $J_\star$ 的选择策略（更“可落地”的写法）

AABB-tree/BVH 的运行时间强依赖数据分布，理论上最坏可退化；因此阈值不宜仅凭理论式，而应结合 **硬内存约束 + 可测交叉点**。

### 4.1 内存硬约束（必须满足）

设内存预算 `MemBudget`（字节），允许 `AllPairs` 使用比例 $\rho\in(0,1)$。每条 pair 存储开销 `sizeof(Pair)`（建议写清：两个 32/64-bit id，或两个指针等）。则
$$
J_\star^{\text{mem}}=
\left\lfloor\frac{\text{MemBudget}\cdot \rho}{\text{sizeof(Pair)}}\right\rfloor,
\qquad
J_\star\le J_\star^{\text{mem}}.
$$

> 额外建议：若 `AllPairs` 采用动态扩容数组，请在达到阈值前预留容量 `reserve(J_star)`，避免扩容抖动。

------

### 4.2 时间权衡：交叉点 $J_{\text{cross}}$（推荐做法）

在发生切换时，Adaptive 的 Phase1 至多枚举 $\le J_\star$ 条 pair（因为超阈值立即停止 report 并丢弃），因此大分支“额外枚举开销”约为 $O(J_\star)$。

**推荐流程**：

1. 在代表性数据上 benchmark：
   - `Baseline`: Enumerate+Sampling（一次扫线枚举 + 数组采样）
   - `Sampling`: 两遍扫线（只计数 + 事件 alias + 局部采样）
2. 用 $|J|$ 作为横轴，找两者耗时大致相等的交叉点 $J_{\text{cross}}$。
3. 取

$$
J_\star^{\text{time}} \approx \alpha\cdot J_{\text{cross}},\qquad \alpha\in[0.5,0.8].
$$

1. 最终阈值：

$$
J_\star=\min(J_\star^{\text{mem}},\ J_\star^{\text{time}}).
$$

> $\alpha$ 取 0.8 的直觉：尽量让“该枚举的都枚举”，但要给分布波动留安全余量，避免切换太晚导致浪费。

------

### 4.3 在线自适应（可选增强，写进附录也很加分）

如果希望“部署时自动适配不同数据分布”，可在扫线前段做轻量 profiling：

- 统计若干次 `CountIntersect` 的平均访问节点数 $\overline{V}$ 与平均耗时；
- 统计若干次 `IntersectReport` 的平均输出密度；
- 估计 “继续枚举” 与 “切换两遍采样” 的边际成本，动态调整 $J_\star$。

这不影响算法正确性，仅优化工程性能。

------

## 5. 算法分析（正确性与复杂度；三版本均包含）

### 5.1 正确性公共基石

#### 引理 A（事件分块不交并）：$J=\biguplus_e J_e$

这就是 1.3 的命题 1，不再重复。其作用是：把“全局均匀”转化为“先选事件、再在事件内均匀”。

#### 引理 B（投影等价）：在 START 时只需判断 $q^\star\cap r^\star$

这就是 1.4 引理 2。

#### 引理 C（CountIntersect 精确）与 引理 D（SampleIntersect 严格均匀 i.i.d.）

来自 2.4 与 2.5 的引理 3–6。其核心依赖 BVH 不变量 (I1)-(I4)。

------

### 5.2 Sampling（两遍扫描）的正确性

记 Phase1 得到 $w_e=|J_e|$，且 $W=\sum_e w_e=|J|$。

**定理 1（单样本全局均匀）**
 对任意 $P\in J$，存在唯一事件 $e^\star$ 使 $P\in J_{e^\star}$。Sampling 输出的任意位置 $j$ 满足
$$
\Pr(Z_j=P)=\frac{1}{|J|}.
$$
**证明**
 Phase2 选事件：$\Pr(E_j=e)=w_e/W$。
 Phase3 在事件 $e$ 内由 `SampleIntersect` 对 $K_e$ 均匀：$\Pr(P\mid E_j=e)=1/w_e$。
 故
$$
\Pr(Z_j=P)=\Pr(E_j=e^\star)\Pr(P\mid E_j=e^\star)
=\frac{w_{e^\star}}{W}\cdot \frac{1}{w_{e^\star}}
=\frac{1}{|J|}.
$$
证毕。

**定理 2（i.i.d.）**
 $Z_1,\dots,Z_t$ 相互独立，且边缘分布均为 $J$ 上均匀。

**证明要点**

- Phase2 中对每个位置 $j$ 的事件抽样使用独立随机数，故 $(E_j)$ 独立；
- 对固定事件 $e$，`SampleIntersect(Q,t_e)` 内部重复独立 `SampleOne`，得到 $t_e$ 个 i.i.d.；
- slot 回填是确定性的索引写入，不引入额外相关性。

------

### 5.3 Enumerate+Sampling（Baseline）的正确性

#### 定理 3（枚举不重不漏：`Pairs` 恰为 $J$）

**不漏**：任意 $(a,b)\in J$，由引理 A 存在唯一较晚 START 事件 $e^\star$，处理该 START 时另一侧盒活跃；由引理 B，投影相交，故 `IntersectReport` 必会访问到对应命中叶并输出该 id；于是 pair 被追加。

**不重**：同一 pair 不可能在更早 START 时出现（较晚 START 盒尚未活跃），也不会在另一 START 再出现（唯一归属）。

故 `Pairs` 一一枚举 $J$。

#### 定理 4（数组采样 i.i.d.）

在 `Pairs[0..|J|-1]` 上独立均匀抽下标，显然得到 $J$ 上 i.i.d. 均匀（有放回）样本。

------

### 5.4 Adaptive+Sampling 的正确性（切换不改变分布）

#### 定理 5（不切换时等价于 Baseline）

若最终 $W\le J_\star$ 且 `mode==ENUMERATE`，则 Phase1 全程枚举且 `AllPairs` 完整等于 $J$，之后数组采样即为 Baseline，正确。

#### 定理 6（切换后等价于 Sampling）

一旦 $W>J_\star$，算法：

- 立即停止 report，不再保留任何 `AllPairs`（避免混合分布难题）；
- 但仍保留 Phase1 对所有事件的精确计数 $w_e$；
- 后续执行 Sampling 的 Phase2+Phase3，仅依赖 $w_e$ 与事件序，不依赖任何已枚举 pair。

因此切换后输出分布与“从一开始就运行 Sampling”完全一致，仍为 $J$ 上 i.i.d. 均匀。

------

### 5.5 复杂度分析（参数化 + 最坏退化说明；三版本分别给出）

由于 BVH 是工程索引，最坏情形可能退化；因此用“访问节点数”参数化描述，并给出常见情况近似。

#### 5.5.1 定义参数（一次查询的代价刻画）

- $h$：BVH 高度（平衡时接近 $O(\log n)$，最坏可到 $O(n)$）
- 对查询盒 $Q$，定义：
  - $V(Q)$：访问到的内部节点数（未被剪枝）
  - $L(Q)$：访问到并进行叶级精确判定的叶数量
  - $k(Q)$：真实输出大小（只对 report 有意义）

在满足 (I1) conservative 的前提下，典型实现满足：

- `CountIntersect(Q)`：$\ \ O(V(Q)+L(Q))$

- `IntersectReport(Q)`：$\ O(V(Q)+L(Q)+k(Q))$（若 out 仅对真命中追加）

- `BuildGuide(Q)`：$\ O(V(Q)+L(Q))$

- `SampleIntersect(Q,t')`：
  $$
  O\bigl(V(Q)+L(Q)+t'\cdot h\bigr)
  $$

插入/删除成本记为 `UpdateCost`（工程上常见均摊 $O(h)$，最坏可退化）。

------

#### 5.5.2 Sampling（两遍扫描）

- **时间**

Phase1：
$$
T_1 = O(n\log n) + O(\text{总更新}) + \sum_{e\in E} O\bigl(V(q_e^\star)+L(q_e^\star)\bigr).
$$
Phase2：
$$
T_2 = O(|E|+t)=O(n+t).
$$
Phase3：
$$
T_3 = O(\text{总更新})
+\sum_{e\in E:t_e>0} O\bigl(V(q_e^\star)+L(q_e^\star)+t_e\cdot h\bigr).
$$
且 $\sum_e t_e=t$，所以采样下降总代价为 $O(t\cdot h)$。

合并：
$$
T_{\text{Sampling}}
=
O(n\log n) + O(\text{两遍总更新})
+ \sum_{e\in E}O(V+L)
+ \sum_{e:t_e>0}O(V+L)
+ O(t\cdot h).
$$

- **常见情形近似**：若 $h=O(\log n)$，且大部分查询 $V(Q)+L(Q)\approx O(\log n)$，则经验上

$$
T_{\text{Sampling}}\approx O((n+t)\log n).
$$

- **最坏退化（必须写清）**：若 AABB 大量重叠导致剪枝失效，可能出现 $V(Q)=\Theta(n)$，从而两遍扫描总体可能退化到 $\Theta(n^2)$。这不影响正确性，只影响性能。
- **空间**
  - 两棵 BVH：$O(n)$
  - `Events,w,alias`：$O(n)$
  - slots + 输出：$O(t)$
  - guide：单次查询临时结构，峰值 $O(\max_Q(V(Q)+L(Q)))\le O(n)$

故峰值空间：
$$
S_{\text{Sampling}}=O(n+t).
$$

------

#### 5.5.3 Enumerate+Sampling（Baseline）

- **时间**
  $$
  T_{\text{Enum}}=
  O(n\log n) + O(\text{总更新})
  + \sum_{e\in E} O\bigl(V(Q_e)+L(Q_e)+k(Q_e)\bigr)
  + O(t).
  $$
  注意 $\sum_e k(Q_e)=|J|$，因此必然含 $|J|$ 级输出代价。

- **空间**
  $$
  S_{\text{Enum}}=O(n+t+|J|),
  $$
  最坏 $|J|=\Theta(n^2)$ 导致爆内存是 baseline 的必然缺陷。

------

#### 5.5.4 Adaptive+Sampling

- **Case A（未切换，$W\le J_\star$）**：等价 baseline，但 $|J|\le J_\star$ 被限制
  $$
  T_{\text{Adap-A}} \approx T_{\text{Enum}} \ \text{且}\ |J|\le J_\star,\qquad
  S_{\text{Adap-A}}=O(n+t+|J|)\le O(n+t+J_\star).
  $$

- **Case B（切换，$W>J_\star$）**：

  - Phase1 枚举代价被强制 $\le O(J_\star)$（因为超阈值立即停止 report 并丢弃）；
  - 后续 Phase2+Phase3 与 Sampling 同阶。

因此（参数化形式）：
$$
T_{\text{Adap-B}}
=
T_{\text{Sampling}}\ +\ O(J_\star),
\qquad
S_{\text{Adap-B}}=O\bigl(n+\max(t,J_\star)\bigr).
$$
