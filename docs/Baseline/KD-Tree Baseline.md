# KD‑Join i.i.d. Uniform Sampling Baseline（强化版）

**关键词：Plane Sweep / Active Set / rank‑embedding / static balanced KD‑tree / orthogonal range query / alias / slot / i.i.d. uniform sampling with replacement**

------

## 记号与约定

- 维度：$d\ge 2$。令 $m=d-1$，嵌入维数 $k=2m=2(d-1)$。

- 两类半开轴对齐盒：
  $$
  R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
  R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\},\quad
  n=n_1+n_2.
  $$

- 每个盒子 $r$：
  $$
  r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\quad L_i(r)<R_i(r).
  $$

- 统一输出 pair 顺序：任何输出都写为 $(r_c,r_{\bar c})$（若 START 来自 $R_{\bar c}$，则输出时交换顺序）。

- 随机数：假设调用的 RNG（如 `rand()` / `mt19937`）能产生相互独立的均匀随机数；alias、KD.Sample、数组采样都使用同一 RNG 流即可（只要调用序列确定）。

------

# 1. 问题定义与分析

## 1.1 输入：两类 $d$ 维半开盒

输入为两类盒集合 $R_c,R_{\bar c}$。半开区间的好处是消除“贴边是否算相交”的歧义，并能和事件顺序（同坐标 END 在 START 前）形成严格一致的语义闭环。

**一维半开区间相交判据：**
$$
[a,b)\cap[c,d)\ne\varnothing \iff \max(a,c)<\min(b,d).
$$

------

## 1.2 Join 结果与采样目标

只关心跨集合相交对：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\ne\varnothing\}.
$$
目标：输出 $t$ 个样本 $Z_1,\dots,Z_t\in J$，要求 **i.i.d. 均匀有放回**：

- **均匀性**：对任意 $P\in J$，$\Pr(Z_j=P)=\frac1{|J|}$；
- **独立性**：$Z_1,\dots,Z_t$ 相互独立。

------

## 1.3 第 1 维 Plane Sweep 与事件系统（“唯一归属”的关键前提）

选择第 1 维 $x_1$ 做扫描。每个盒子 $r$ 产生两事件：

- $\texttt{START}(r)$ at $x_1=L_1(r)$
- $\texttt{END}(r)$ at $x_1=R_1(r)$

### 事件全序（必须稳定排序）

对所有事件定义一个**全序**（stably sorted total order）：

1. 按 $x_1$ 升序；
2. 同 $x_1$：**END 在 START 之前**（与半开一致，避免边界贴合误算相交）；
3. 同 $x_1$ 且同为 START：按固定优先级（例如先 $R_c$ 后 $R_{\bar c}$），再按全局唯一 `id` 打破平局；
4. 同 $x_1$ 且同为 END：相对顺序无关（可按 id）。

### 扫描时序约定（必须写死）

- 处理 $\texttt{START}(q)$：**先查询对立活跃集，再激活 $q$**；
- 处理 $\texttt{END}(r)$：先失活 $r$。

> 这两条 + 半开 + END-before-START + START tie-break，合起来才保证后续“事件块分解不重不漏”。

### 活跃集定义

当扫描线位于 $x_0$ 时：
$$
\mathrm{Active}(R_c)=\{r_c\in R_c: x_0\in[L_1(r_c),R_1(r_c))\},
\quad \mathrm{Active}(R_{\bar c})\ \text{同理}.
$$

------

## 1.4 事件诱导的块分解（全局 Join 的不交并分解）

令 $E$ 为所有 START 事件集合。对任一 START 事件 $e\in E$，设新盒子 $q=r(e)$。

定义对立活跃集（Other side）：
$$
\mathrm{Other}(q)=
\begin{cases}
R_{\bar c}, & q\in R_c,\\
R_c, & q\in R_{\bar c}.
\end{cases}
$$
在处理 $\texttt{START}(q)$ 时，定义：
$$
K_e := \{r\in \mathrm{Active}(\mathrm{Other}(q)) : q\cap r\ne\varnothing\},\qquad
w_e:=|K_e|.
$$
定义局部 join 块 $J_e$（统一输出顺序为 $(r_c,r_{\bar c})$）：

- 若 $q\in R_c$：$\displaystyle J_e=\{(q,r_{\bar c})\mid r_{\bar c}\in K_e\}$
- 若 $q\in R_{\bar c}$：$\displaystyle J_e=\{(r_c,q)\mid r_c\in K_e\}$

### 引理 1（事件分块不交分解）

在 §1.3 的事件全序 + 半开语义 + “START 先查后插”下：
$$
J=\biguplus_{e\in E}J_e,\qquad |J|=\sum_{e\in E}w_e.
$$
**证明要点（写成可审稿的版本）：**

对任意相交对 $(r_c,r_{\bar c})\in J$，考虑它们各自的 START 事件（在事件全序下）$s_c=\texttt{START}(r_c)$, $s_{\bar c}=\texttt{START}(r_{\bar c})$。由于全序，二者必有唯一“更靠后”的 START：
$$
e^\star := \max_{\prec}\{s_c,s_{\bar c}\}.
$$

- **(不漏)** 在处理 $e^\star$ 时，另一个盒子 $r'$ 的 START 早于 $e^\star$，且它尚未 END（否则二者在第 1 维不可能相交，或者同坐标下已被 END-before-START 删除），因此 $r'$ 必在对立活跃集中；又因为两盒确实相交，所以 $r'\in K_{e^\star}$，因此该 pair 会被输出到 $J_{e^\star}$。
- **(不重)** 在处理较早的 START 时，另一个盒子尚未被激活（同坐标 START 的先后由 tie-break 决定，确保只有“后处理者”能看到“先处理者”），因此该 pair 不会在较早 START 的 $J_e$ 中出现。不同事件集合显然不交。

从而得到不交并分解与计数等式。

------

## 1.5 第 1 维相交的“自动满足”与降维过滤（只需检查维度 2..d）

在处理 $\texttt{START}(q)$ 时，对立活跃盒 $r\in \mathrm{Active}(\mathrm{Other}(q))$ 必满足：
$$
L_1(r)\le L_1(q) < R_1(r)
$$
这是由“活跃定义 + END-before-START”直接推出的，因此二者在第 1 维一定真重叠（非贴边）。

于是：
$$
q\cap r\ne\varnothing
\iff
\text{在维度 }2..d\text{ 的投影相交}.
$$

------

## 1.6 维度 2..d 的相交等价：点嵌入 + 正交范围查询

令 $m=d-1$，只取维度 $2..d$ 的端点向量：
$$
\mathbf L(r)=(L_2(r),\dots,L_d(r))\in\mathbb R^{m},\quad
\mathbf R(r)=(R_2(r),\dots,R_d(r))\in\mathbb R^{m}.
$$
定义 $k=2m$ 维点嵌入：
$$
p(r)=(L_2(r),\dots,L_d(r),\ R_2(r),\dots,R_d(r))\in\mathbb R^{2m}.
$$
对查询盒 $q$，对任意 $i\ge 2$：
$$
[L_i(r),R_i(r))\cap[L_i(q),R_i(q))\ne\varnothing
\iff
L_i(r)<R_i(q)\ \land\ R_i(r)>L_i(q).
$$
因此定义半无限正交范围：
$$
\mathcal Q(q)=
\prod_{i=2}^{d}(-\infty,R_i(q))
\times
\prod_{i=2}^{d}(L_i(q),+\infty).
$$

### 引理 2（嵌入等价性）

$$
r \text{ 与 } q \text{ 在维度 }2..d\text{ 相交}
\iff
p(r)\in \mathcal Q(q).
$$

------

## 1.7 三个版本的定位（Baseline 三件套）

- **Sampling（两遍扫描）**：不显式枚举 $|J|$，但依赖 KD 范围查询性能；最坏可能退化。
- **Enumerate+Sampling（一次扫描枚举 + 数组采样）**：实现最直观，常数小，但时间/空间线性依赖 $|J|$。
- **Adaptive+Sampling（阈值切换）**：$|J|$ 小时只扫一遍更快；$|J|$ 大时避免 $|J|$ 级枚举存储，切到 Sampling。

------

# 2. 核心数据结构

Baseline 需要四个组件：

1. 事件数组 `Events`
2. 严格不等号的 **rank 化**（将 $\mathcal Q(q)$ 变成整数闭区间盒 $Q(q)$）
3. 两棵静态平衡 KD‑tree（分别对应 $R_c$ 与 $R_{\bar c}$ 的点集），扫描时只做 active/cnt 更新
4. alias + slot（事件级采样规划；KD.Sample 内部也可用 piece‑slot）

------

## 2.1 事件数组 `Events`

每个事件记录：
$$
(x_1,\ \mathrm{type}\in\{\texttt{START},\texttt{END}\},\ \mathrm{class}\in\{c,\bar c\},\ \mathrm{id}).
$$

- `id` 为全局唯一盒子编号（1..n）
- `class` 标记盒子属于哪一侧
- 按 §1.3 的全序稳定排序

> **工程建议：** 用一个数组保存排序后的事件；另外保存 `start_index_of_box[id]` 方便 debug，但不是必须。

------

## 2.2 严格不等号 rank 化（二级键，保证“严格 < / >”可实现）

KD‑tree 通常在整数/可比较坐标上做闭区间查询最方便。为了把严格不等号稳定映射到闭区间，用二级键 $(\text{value},\text{id})$ 进行坐标压缩。

### 2.2.1 全局唯一 id 与哨兵

每个盒子分配全局唯一：
$$
\mathrm{id}(r)\in\{1,2,\dots,n\}.
$$
设置哨兵：
$$
\mathrm{LOW}=0,\qquad \mathrm{HIGH}=n+1.
$$

### 2.2.2 两类排序数组（每个维度各两份）

对每个维度 $i\in\{2,\dots,d\}$：

- $A_i$：排序 $\{(L_i(r),\mathrm{id}(r))\mid r\in R_c\cup R_{\bar c}\}$
- $B_i$：排序 $\{(R_i(r),\mathrm{id}(r))\mid r\in R_c\cup R_{\bar c}\}$

二级键按字典序比较。

定义 rank（0-based）：
$$
\mathrm{rank}_{A_i}(r)=\text{pair }(L_i(r),\mathrm{id}(r))\text{ 在 }A_i\text{ 中的位置},
\quad
\mathrm{rank}_{B_i}(r)=\text{pair }(R_i(r),\mathrm{id}(r))\text{ 在 }B_i\text{ 中的位置}.
$$
于是嵌入点变为整数点：
$$
p(r)=(
\mathrm{rank}_{A_2}(r),\dots,\mathrm{rank}_{A_d}(r),\
\mathrm{rank}_{B_2}(r),\dots,\mathrm{rank}_{B_d}(r)
)\in\mathbb Z^{k}.
$$

### 2.2.3 查询范围从 $\mathcal Q(q)$ 到整数闭区间盒 $Q(q)$

对任意维 $i\ge 2$：

- 约束 $L_i(r)<R_i(q)$（严格 <）：

  用
  $$
  U_i=\mathrm{lower\_bound}_{A_i}\big((R_i(q),\mathrm{LOW})\big)-1
  $$
  得到 rank 区间：
  $$
  \mathrm{rank}_{A_i}(r)\in[0,U_i].
  $$
  **解释：** `lower_bound((R,LOW))` 返回第一个 $\ge (R,LOW)$ 的位置；由于所有 id $\ge 1>\mathrm{LOW}$，所以任何 $(R,\mathrm{id})$ 都 $\ge (R,\mathrm{LOW})$，从而 $\le U_i$ 恰好等价于 value < R。

- 约束 $R_i(r)>L_i(q)$（严格 >）：

  用
  $$
  L'_i=\mathrm{upper\_bound}_{B_i}\big((L_i(q),\mathrm{HIGH})\big)
  $$
  得到 rank 区间：
  $$
  \mathrm{rank}_{B_i}(r)\in[L'_i, |B_i|-1].
  $$
  **解释：** 由于所有 id $\le n < \mathrm{HIGH}$，所以所有 $(L,\mathrm{id})\le (L,\mathrm{HIGH})$。`upper_bound((L,HIGH))` 返回第一个 $>(L,HIGH)$，即 value > L 的第一个位置，从而 $\ge L'_i$ 等价于 value > L。

若存在某维 $U_i<0$ 或 $L'_i>|B_i|-1$，则 $Q(q)$ 为空，命中集合为空。

最终闭区间盒：
$$
Q(q)=
\prod_{i=2}^{d}[0,U_i]
\times
\prod_{i=2}^{d}[L'_i,|B_i|-1]
\subseteq \mathbb Z^{k}.
$$

### 引理 2′（严格不等号在 rank 空间被精确实现）

对任意盒 $r$ 与查询 $q$：
$$
p(r)\in \mathcal Q(q)\quad\Longleftrightarrow\quad p(r)\in Q(q).
$$

------

## 2.3 两棵静态平衡 KD‑tree：骨架静态，Active 动态

构建两棵静态平衡 KD‑tree（典型用 median split）：

- $\mathcal{KD}^{(c)}$：点集 $\{p(r_c)\mid r_c\in R_c\}$
- $\mathcal{KD}^{(\bar c)}$：点集 $\{p(r_{\bar c})\mid r_{\bar c}\in R_{\bar c}\}$

扫描阶段**不做结构性插删**，只做：

- `active(v) ∈ {0,1}`：结点对应点是否活跃
- `cnt(v)`：子树内活跃点数（包含结点自身点）
- `bbox(v)`：子树点集的包围盒（静态预计算；每维 min/max，闭区间）

并维护 `node_ptr[id] -> node`，用于 O(height) 找到结点并沿 parent 更新 cnt。

### 2.3.1 结点字段（建议实现）

每个结点 $v$：

- `pt(v) ∈ Z^k`：点坐标
- `id(v)`：对应盒子 id
- `left(v), right(v), parent(v)`
- `split_dim(v)`：该结点的切分维（常用 depth mod k）
- `bbox(v)`：静态 bounding box（k 维闭区间）
- `active(v)`：动态 0/1
- `cnt(v)`：动态整数

**可选增强：Bucket（加速子树内均匀抽样与整段报告）**

- `Bucket(v)`：存“当前子树内活跃点的 id 列表”，并维护 id→位置句柄以支持 O(1) 删除（swap-delete）
- 若对每个点沿 ancestor 更新 Bucket，则总空间 $\Theta(nh)=\Theta(n\log n)$，能把子树内均匀抽样降到 O(1)

> 你们原稿已提 Bucket，这里把“如何维护/代价是什么/能换来什么”写清楚。KD-Tree Baseline

### 2.3.2 Activate / Deactivate（动态标记）

**Activate(id)**

1. $v\leftarrow node\_ptr[id]$
2. 若 `active(v)=0`：置 1；并沿 parent 向上对 `cnt(u)` 加 1（启用 Bucket 则同步插入）

**Deactivate(id)** 对称（置 0，沿 parent 向上减 1；Bucket 删除）

**高度假设：** 若 KD‑tree 用 median split（或等价平衡策略）构建，则高度 $h=O(\log n)$（平均/典型）。我们在复杂度里将 Update 记为 $U(n)=O(h)$。

------

## 2.4 KD‑tree 查询接口：Count / Report / Sample

对闭区间 query box $Q=\prod_{j=1}^{k}[L_j,U_j]\subseteq\mathbb Z^k$，定义命中集合：
$$
\mathrm{Hit}(Q)=\{r\ \text{active}:\ p(r)\in Q\}.
$$
我们需要支持：

- `Count(Q)`：返回 $|\mathrm{Hit}(Q)|$
- `Report(Q)`：枚举 $\mathrm{Hit}(Q)$ 的所有 id
- `Sample(Q,t′)`：从 $\mathrm{Hit}(Q)$ 上 **i.i.d. 均匀有放回**采样 $t′$ 次

------

### 2.4.1 Count(Q)：精确计数

递归 `Count(v,Q)`：

- 若 $v=\varnothing$ 或 `cnt(v)=0`：返回 0

- 若 `bbox(v)` 与 $Q$ 不相交：返回 0

- 若 `bbox(v) ⊆ Q`：返回 `cnt(v)`

- 否则：
  $$
  ans=\mathbf 1[\texttt{active}(v)=1\land pt(v)\in Q]+Count(left,Q)+Count(right,Q).
  $$

**剪枝正确性理由：**

- `bbox ∩ Q = ∅` ⇒ 子树无点命中；
- `bbox ⊆ Q` ⇒ 子树所有活跃点都命中，因此计数等于 `cnt`。

------

### 2.4.2 Report(Q)：精确枚举

递归 `Report(v,Q)`：

- 若 $v=\varnothing$ 或 `cnt(v)=0` 或 `bbox(v)∩Q=∅`：返回空
- 若 `bbox(v)⊆Q`：
  - 若启用 Bucket：直接输出 `Bucket(v)`（全部活跃点 id）
  - 否则：做一次“子树遍历输出活跃点”（代价 $\Theta(\mathrm{cnt}(v))$）
- 否则（部分覆盖）：
  - 若 `active(v)=1` 且 `pt(v)∈Q`：输出 `id(v)`
  - 递归左右子树

------

### 2.4.3 Sample(Q,t′)：canonical pieces + piece‑slot（保证“有序 i.i.d.”）

#### 目标

输出序列 $(X_1,\dots,X_{t′})$，满足：

- 每个 $X_u$ 在 $\mathrm{Hit}(Q)$ 上均匀（概率 $1/|\mathrm{Hit}(Q)|$）
- 且 $X_1,\dots,X_{t′}$ 相互独立（有放回）

#### Step 1：构造 canonical pieces（不交并分解）

递归 `CollectPieces(v,Q)` 生成 piece 列表 $\mathcal P$。piece 有两种：

- `SubtreePiece(v)`：表示“结点 $v$ 子树内所有活跃点”（前提 `bbox(v)⊆Q`）
- `PointPiece(v)`：表示“单点 {id(v)}”（前提 `active(v)=1` 且 `pt(v)∈Q` 且当前是部分覆盖场景）

伪代码：

```
CollectPieces(v,Q):
  if v==null or cnt(v)==0: return
  if bbox(v) ∩ Q == ∅: return
  if bbox(v) ⊆ Q:
      add SubtreePiece(v) with weight w = cnt(v)
      return
  else:
      if active(v)==1 and pt(v)∈Q:
          add PointPiece(v) with weight 1
      CollectPieces(left(v),Q)
      CollectPieces(right(v),Q)
```

令每个 piece $P$ 对应集合 $S_P$：

- 若 $P=\text{PointPiece}(v)$，则 $S_P=\{id(v)\}$
- 若 $P=\text{SubtreePiece}(v)$，则 $S_P=\{\text{active points in subtree}(v)\}$

定义权重 $w(P)=|S_P|$，以及总权重：
$$
W=\sum_{P\in\mathcal P} w(P).
$$
若 $W=0$ 直接返回空序列（无命中点）。

##### 引理 3（pieces 构成不交并分解）

$$
\mathrm{Hit}(Q)=\biguplus_{P\in\mathcal P} S_P,\quad \text{且 } w(P)=|S_P|.
$$

**证明要点：**

- `bbox⊆Q` 时停止展开，确保每个命中点落入“最高层”的某个 SubtreePiece（不会被子节点再覆盖）。
- 部分覆盖时只把当前结点自身点（若命中）作为 PointPiece，子树点交给递归，保证不漏。
- 由于 KD‑tree 是“点唯一存储”（每点只在一个结点），PointPiece 之间不重；SubtreePiece 之间也不重（它们来自不同子树或不同停止点）；且 SubtreePiece 不会与其后代 piece 重叠（因为停止展开）。

#### Step 2：在 pieces 上按权重抽样，并用 slot 归档

在 $\mathcal P$ 上构建离散分布：
$$
\Pr(P)=\frac{w(P)}{W}.
$$
对 $u=1..t′$：

- 独立抽取 piece $P_u\sim \Pr(\cdot)$
- 把位置 $u$ 记录到 `Slot[P_u]`

> 这里 slot 的作用是：把“逐次抽 piece”改成“按 piece 分组批量生成”，但**保持每个位置的抽样决策独立**。

#### Step 3：piece 内均匀采样并回填

对每个 piece $P$：

- 若 $P=\text{PointPiece}(v)$：对 `Slot[P]` 内所有位置直接填 `id(v)`（允许重复）
- 若 $P=\text{SubtreePiece}(v)$：需要从“结点 $v$ 子树活跃点集合”上均匀采样 $|Slot[P]|$ 次

提供两种实现：

**实现 C‑1（无需 Bucket，按 cnt 权重随机下降；单次 O(h)）**

```
SampleSubtreeUniform(v):
  // assume cnt(v)>0
  w_self = 1 if active(v)==1 else 0
  w_L = cnt(left(v))  if left!=null else 0
  w_R = cnt(right(v)) if right!=null else 0
  T = w_self + w_L + w_R = cnt(v)

  draw u uniformly in {1..T}
  if u <= w_self: return id(v)
  else if u <= w_self + w_L: return SampleSubtreeUniform(left(v))
  else: return SampleSubtreeUniform(right(v))
```

对 `Slot[SubtreePiece(v)]` 中每个位置调用一次。

**实现 C‑2（启用 Bucket；单次 O(1)）**

- `Bucket(v)` 保存子树活跃点 id 的紧凑数组
- 均匀抽样：随机下标访问 `Bucket(v)[rand % cnt(v)]`

#### 引理 4（KD.Sample(Q,t′) 输出 i.i.d. 均匀）

对任一命中点 $x\in\mathrm{Hit}(Q)$，单次输出满足：
$$
\Pr(\text{output}=x)=\frac1{|\mathrm{Hit}(Q)|}.
$$
并且 $t′$ 次输出相互独立（有放回）。

**证明要点：**

- 由引理 3，$x$ 属于唯一 $P^\star$ 且 $w(P^\star)=|S_{P^\star}|$。
- 抽到 piece 的概率：$w(P^\star)/W$；
- 给定 piece，内部均匀抽到 $x$ 的概率：$1/w(P^\star)$；
- 乘积得 $1/W$；且 slot 只是把每次独立抽 piece 的结果按组汇总回填，不改变每个位置的独立性。

------

## 2.5 alias + slot：事件级采样规划（Sampling / Adaptive 大分支使用）

- 事件权重：每个 START 事件 $e$ 的 $w_e=|K_e|=|J_e|$

- 构建事件级 alias（或前缀和）：
  $$
  \Pr(E=e)=\frac{w_e}{|J|},\quad |J|=\sum_{e\in E}w_e.
  $$

- slot：对每个事件 $e$ 维护 $S_e\subseteq\{1,\dots,t\}$ 表示哪些输出位置分配给该事件（由 Phase2 随机产生）

> **工程建议：** 构建 alias 时可直接过滤掉 $w_e=0$ 的事件（省空间，也避免边界处理）。

------

# 3. 算法详细流程（三版本）

下面给出三版本都能直接实现的流程（含关键 corner cases）。

------

## 3.0 公共预处理（所有版本共享）

输入：$R_c,R_{\bar c}$，样本数 $t$（Adaptive 另有 $J_\star$）。

**Preprocess：**

1. 为每个盒子分配全局唯一 id（1..n），保存 `id -> box` 映射。
2. 构造 `Events`：每个盒子产生 START/END 各一个。
3. 按 §1.3 全序稳定排序 `Events`。
4. 对每个维度 $i=2..d$ 构建排序数组 $A_i,B_i$，并为每个盒子计算整数点 $p(r)\in\mathbb Z^{k}$。
5. 构建两棵平衡 KD‑tree：
   - $\mathcal{KD}^{(c)}$ 存 $R_c$ 的点
   - $\mathcal{KD}^{(\bar c)}$ 存 $R_{\bar c}$ 的点
      同时填好 `node_ptr[id]`、静态 `bbox`、父指针等。
6. 实现 `BuildRange(q)`：按 §2.2 将 $q$ 映射为闭区间盒 $Q(q)$。若空则返回 `EMPTY`。

初始化：两棵 KD‑tree 全部 `active=0`，`cnt(root)=0`（Bucket 若启用也为空）。

------

## 3.1 版本一：Sampling（两遍扫描 + 事件 alias+slot + 局部 KD.Sample）

### 适用场景

- 不希望时间/空间线性依赖 $|J|$
- 但接受 KD‑tree 在高维/坏分布下可能退化（见 §5.3）

------

### Phase 1：第一次扫描（只计数，得到所有 $w_e$ 与 $W=|J|$）

初始化：清空两棵 KD 活跃状态；数组 `w_e`（只对 START 事件存）。`W=0`。

按 `Events` 顺序扫描：

- 若事件为 `END(r)`：对所属类 KD 执行 `Deactivate(id(r))`
- 若事件为 `START(q)`（该 START 对应事件 $e$）：
  1. $Q\leftarrow \texttt{BuildRange}(q)$。若 `EMPTY`，置 $w_e=0$。
  2. 若 $q\in R_c$：$w_e \leftarrow \mathcal{KD}^{(\bar c)}.\texttt{Count}(Q)$
      若 $q\in R_{\bar c}$：$w_e \leftarrow \mathcal{KD}^{(c)}.\texttt{Count}(Q)$
  3. `W += w_e`
  4. `Activate(id(q))` 到本侧 KD（注意：**查询在激活之前**）

扫描结束：由引理 1 得 $W=\sum_e w_e=|J|$。若 $W=0$ 直接返回空。

------

### Phase 2：事件级 alias + slot 分配（生成“采样计划”）

1. 在 START 事件集合上，过滤 $w_e>0$ 的事件，构建 alias（或前缀和）。
2. 初始化每个事件的 slot 容器：$S_e\gets\emptyset$。
3. 对每个输出位置 $j=1..t$：
    抽事件 $E_j\sim (w_e/W)$，并把 $j$ push 到 $S_{E_j}$。

记 $t_e=|S_e|$，则 $\sum_e t_e=t$。

> **实现细节：**
>  因为 $j$ 从小到大遍历并 push，故每个 $S_e$ 内天然是递增序列，无需额外排序。

------

### Phase 3：第二次扫描（按事件批量 KD.Sample + 回填）

1. 再次清空两棵 KD 的活跃状态（结构与 bbox 不变）。
2. 初始化输出数组 `Ans[1..t]`。
3. 再按同一 `Events` 顺序扫描：

- `END(r)`：Deactivate
- `START(q)`（事件 $e$）：
  - 若 $t_e=0$：跳过采样，仅 `Activate(id(q))`
  - 若 $t_e>0$：
    1. $Q\leftarrow BuildRange(q)$。理论上 $Q$ 不可能为空（因为 $w_e>0$ 才会被采到），但工程上可 `assert`。
    2. 若 $q\in R_c$：`List ← KD^(barc).Sample(Q, t_e)`
        若 $q\in R_{\bar c}$：`List ← KD^(c).Sample(Q, t_e)`
    3. 设 $S_e=(j_1,\dots,j_{t_e})$，回填：
       - 若 $q\in R_c$：`Ans[j_u] = (q, box(List[u]))`
       - 若 $q\in R_{\bar c}$：`Ans[j_u] = (box(List[u]), q)`
  - 最后 `Activate(id(q))`

输出 `Ans`。

------

## 3.2 版本二：Enumerate+Sampling（一次扫描枚举全部 $J$ + 数组采样）

### 适用场景

- $|J|$ 较小（或内存足够）
- 希望实现简单、常数小、调试方便

------

### Step 1：一次扫描枚举全部 $J$

初始化：`Pairs=[]`；两棵 KD 全 inactive。

扫描 `Events`：

- `END(r)`：Deactivate
- `START(q)`：
  1. $Q\leftarrow BuildRange(q)$。若空，则 `List=[]`
  2. 若 $q\in R_c$：`List ← KD^(barc).Report(Q)`
      若 $q\in R_{\bar c}$：`List ← KD^(c).Report(Q)`
  3. 对每个命中 id 形成 pair（统一顺序）：
     - 若 $q\in R_c$：`Pairs.append((q, r_bar))`
     - 若 $q\in R_{\bar c}$：`Pairs.append((r_c, q))`
  4. `Activate(id(q))`

扫描结束后，`|Pairs|=|J|`（由引理 1 保证不重不漏）。

------

### Step 2：数组 i.i.d. 下标采样（有放回）

- 若 `Pairs` 为空：返回空；
- 否则对 $j=1..t$：
   `idx_j ~ Unif{0..|Pairs|-1}`，`Ans[j]=Pairs[idx_j]`。

------

## 3.3 版本三：Adaptive+Sampling（阈值切换）

### 目标

- 当 $|J|$ 小：只扫一遍并枚举（更快）
- 当 $|J|$ 大：避免枚举爆内存，切换到 Sampling（两遍扫描）

------

### Phase 1：一次扫描（永远 Count；在阈值内可选枚举）

输入阈值 $J_\star$。

初始化：

- `mode = ENUMERATE`
- `AllPairs=[]`
- `W=0`
- 保存所有 START 的 `w_e`
- 两棵 KD 全 inactive

扫描 `Events`：

- `END(r)`：Deactivate
- `START(q)`（事件 $e$）：
  1. **先计数（无论是否枚举都做）**
     - $Q\leftarrow BuildRange(q)$，若空则 $w_e=0$
     - 否则：`w_e ← KD_other.Count(Q)`
     - `W += w_e`，保存 `w_e`
  2. **阈值判断 + 可选枚举**
     - 若 `mode==ENUMERATE` 且 `W ≤ J_star`：
       - `List ← KD_other.Report(Q)`（若 Q 空则 List 空）
       - 把 List 转成 pairs append 到 `AllPairs`（统一顺序）
     - 若 `mode==ENUMERATE` 且 `W > J_star`：
       - `mode ← COUNT_ONLY`
       - 释放/丢弃 `AllPairs`（防止继续占内存）
       - 后续不再 Report
  3. `Activate(id(q))`

扫描结束：得到精确 $W=|J|$。若 $W=0$ 返回空。

------

### 分支 A：未切换（$W\le J_\star$）

此时 `AllPairs` 完整枚举 $J$，直接做数组采样（同 §3.2 Step2）。

------

### 分支 B：已切换（$W>J_\star$）

执行 Sampling 的 Phase2 + Phase3：

- Phase2：用 Phase1 保存的 $w_e$ 与 $W$ 构建事件 alias + slot
- Phase3：二遍扫描调用 `KD.Sample` 生成样本

------

# 4. Adaptive 阈值 $J_\star$ 的选择策略（更可落地）

阈值设计的本质是：
 **枚举分支的成本 $\approx$ 写出 $|J|$ 个 pair 的成本**
 vs
 **切换分支的成本 $\approx$ 多做一次完整扫描 + 为 $t$ 个样本做局部 KD.Sample 的成本**

我们给出从“硬约束”到“可标定”的三层策略。

------

## 4.1 内存硬约束（必须满足）

设可用内存预算 `MemBudget` 字节，允许 `AllPairs` 使用比例 $\rho\in(0,1)$。每个 pair 存储开销：

- 理想情况：两个 32-bit id ⇒ 8B；两个 64-bit ⇒ 16B
- 真实情况：vector/结构体开销、对齐 padding、allocator overhead 等通常更大

建议用保守估计 `sizeof(Pair)`（例如 16B 或 24B）。

则：
$$
J_\star^{\mathrm{mem}}=
\left\lfloor \frac{\rho\cdot \mathrm{MemBudget}}{\mathrm{sizeof(Pair)}}\right\rfloor,\qquad
J_\star\le J_\star^{\mathrm{mem}}.
$$

> 工程建议：取 $\rho\in[0.3,0.6]$ 给系统与其他结构留余量。

------

## 4.2 时间权衡：用“多一遍扫描”换“不写 $|J|$”

记：

- 更新代价（Activate/Deactivate）：$U(n)\approx O(h)$
- KD 范围遍历代价：用 $Q_{\mathrm{KD}}(n,k)$ 表示（见 §5.3）

一次“只计数扫描”的近似：
$$
T_{\mathrm{sweep}}
\approx 2n\cdot U(n) + n\cdot Q_{\mathrm{KD}}(n,k).
$$
切换后会多一遍扫描（Phase3），因此多出来的“固定成本”约为 $T_{\mathrm{sweep}}$ 级别；同时还要生成 $t$ 个样本（总体 $O(t\cdot H)$，其中 $H=h$ 或 1 取决于 Bucket）。

因此一个可落地的经验阈值：
$$
J_\star^{\mathrm{time}}
= C_1\cdot T_{\mathrm{sweep}} + C_2\cdot t\cdot H,
\qquad
J_\star=\min(J_\star^{\mathrm{mem}},\ J_\star^{\mathrm{time}}),
$$
其中 $C_1,C_2$ 用 benchmark 标定（一般是 0.5~3 的量级）。

------

## 4.3 交叉点拟合（推荐、最稳健）

1. 固定数据分布、维数 $d$、样本数 $t$，对不同 $n$ 运行：

   - Enumerate+Sampling
   - Sampling

2. 找到耗时交叉点对应的 $|J_{\mathrm{cross}}|$。

3. 取
   $$
   J_\star \approx \alpha\cdot |J_{\mathrm{cross}}|,\quad \alpha\in[0.6,0.9]
   $$
   并用 $J_\star^{\mathrm{mem}}$ 截断。

------

# 5. 算法分析（正确性 + 复杂度；三版本都覆盖）

------

## 5.1 正确性：基础引理

### 引理 1（事件分块不交分解）

已在 §1.4 给出：
$$
J=\biguplus_{e\in E}J_e,\qquad |J|=\sum_{e\in E}w_e.
$$

### 引理 2（点嵌入与范围查询等价）

已在 §1.6 给出：
$$
r \text{ 与 } q \text{ 在维度 }2..d\text{ 相交}\iff p(r)\in\mathcal Q(q).
$$

### 引理 2′（rank 化精确实现严格不等号）

已在 §2.2 给出：
$$
p(r)\in \mathcal Q(q)\iff p(r)\in Q(q).
$$

### 引理 3（KD.Count / KD.Report 语义正确）

`bbox` 剪枝与 `bbox⊆Q` 整段接受保证 Count 精确、Report 精确枚举（只统计 active 点，正确性直接来自包围盒包含关系）。

### 引理 4（KD.Sample 输出 i.i.d. 均匀）

已在 §2.4.3 给出：canonical pieces 不交并 + 加权抽 piece + piece 内均匀抽样 ⇒ i.i.d. uniform with replacement。

------

## 5.2 三个版本正确性结论

### 定理 A（Enumerate+Sampling 的 i.i.d. 均匀性）

- 枚举阶段：由引理 1（事件唯一归属），每个 pair 恰好被加入一次，故 `Pairs` 是 $J$ 的一一枚举。
- 采样阶段：对 `Pairs` 做独立均匀下标采样，得到 $J$ 上 i.i.d. 均匀有放回样本。

------

### 定理 B（Sampling 的 i.i.d. 均匀性）

对任意 $P\in J$，存在唯一事件 $e^\star$ 使 $P\in J_{e^\star}$ 且 $|J_{e^\star}|=w_{e^\star}$。

对任一输出位置 $j$：

- $\Pr(E_j=e)=w_e/W$
- 给定 $E_j=e^\star$，局部 `KD.Sample(Q(q),1)` 在 $K_{e^\star}$ 上均匀，故 $\Pr(P\mid E_j=e^\star)=1/w_{e^\star}$

因此：
$$
\Pr(Z_j=P)=\frac{w_{e^\star}}{W}\cdot\frac{1}{w_{e^\star}}
=\frac1W=\frac1{|J|}.
$$
独立性来自：事件抽样对不同 $j$ 独立；且每个事件内 KD.Sample 产生 i.i.d. 序列，slot 回填只是位置映射，不引入相关。

------

### 定理 C（Adaptive 的正确性）

- 未切换：等价 Enumerate+Sampling ⇒ 正确；
- 已切换：Phase1 已保存精确 $w_e,W$，后续完全等价 Sampling ⇒ 正确；
- 切换只影响“是否保留枚举结果”，不改变任何计数与随机过程 ⇒ 不会引入分布偏差。

------

## 5.3 复杂度（参数化表达 + 退化说明）

### 5.3.1 记号与关键参数

- KD‑tree 高度：$h$（平衡构造通常 $h=O(\log n)$）

- Update 代价：$U(n)=O(h)$

- 范围遍历代价：用
  $$
  Q_{\mathrm{KD}}(n,k)
  $$
  表示一次 `Count/CollectPieces/Report` 访问的结点/检查成本（与分布、维数 $k=2(d-1)$ 强相关）

- 子树内均匀抽样单次代价：

  - 无 Bucket：$H=h$
  - 有 Bucket：$H=1$

------

### 5.3.2 预处理成本

- 构事件并排序：$O(n\log n)$
- 构 $A_i,B_i$ 并计算 rank：$O((d-1)\,n\log n)$（$d$ 常数 ⇒ $O(n\log n)$）
- 构建两棵 KD‑tree：
  - 构造（median split）：典型 $O(n\log n)$
  - `bbox` 预计算：$O(nk)$（$k$ 常数 ⇒ $O(n)$）

空间：

- KD‑tree（无 Bucket）：$O(n)$
- Bucket（若启用）：$\Theta(nh)=\Theta(n\log n)$

------

### 5.3.3 Enumerate+Sampling

**时间：**

- 扫描更新：$2n\cdot U(n)=O(nh)$

- 每个 START 一次 `Report(Q)`：
  $$
  O(Q_{\mathrm{KD}}(n,k)+k_{\text{out}}(e))
  $$

- $\sum_e k_{\text{out}}(e)=|J|$

- 数组采样：$O(t)$

因此：
$$
T_{\mathrm{enum}}
=
O\bigl(n\log n + nh + n\cdot Q_{\mathrm{KD}}(n,k) + |J| + t\bigr).
$$
**空间：**
$$
S_{\mathrm{enum}}
=
O(\mathrm{KD\text{-}space}+|J|+t).
$$

------

### 5.3.4 Sampling

Phase1（Count）：
$$
T_1=O(nh)+n\cdot O(Q_{\mathrm{KD}}(n,k)).
$$
Phase2（alias + slot）：
$$
T_2=O(n+t),\qquad S_2=O(n+t).
$$
Phase3（第二遍扫描 + KD.Sample）：

- 扫描更新：$O(nh)$
- 设 $M=\#\{e: t_e>0\}\le \min(n,t)$
- 每个这样的事件一次 `KD.Sample(Q,t_e)`：
  - pieces 构造：$O(Q_{\mathrm{KD}}(n,k))$
  - 生成 $t_e$ 个样本：$O(t_e\cdot H)$

汇总 $\sum_e t_e=t$：
$$
T_3
=
O(nh) + M\cdot O(Q_{\mathrm{KD}}(n,k)) + O(t\cdot H).
$$
合并：
$$
T_{\mathrm{samp}}
=
O\Bigl(
n\log n + nh + (n+M)\cdot Q_{\mathrm{KD}}(n,k) + t\cdot H + (n+t)
\Bigr).
$$
空间：
$$
S_{\mathrm{samp}}=O(\mathrm{KD\text{-}space}+n+t).
$$

------

### 5.3.5 Adaptive+Sampling

- 未切换（$W=|J|\le J_\star$）：同 Enumerate+Sampling，但 $|J|\le J_\star$
- 已切换（$W>J_\star$）：
  - Phase1 最多枚举写入 $O(J_\star)$ 个 pair（越界即停止），其余只 Count
  - 后续同 Sampling

$$
T_{\mathrm{adaptive,\ switch}}
\approx T_{\mathrm{samp}} + O(J_\star),
\qquad
S_{\mathrm{adaptive,\ switch}}
=O(\mathrm{KD\text{-}space}+\max(J_\star,t)+n).
$$

------

### 5.3.6 最坏情况退化说明（Baseline 的“诚实边界”）

由于 KD‑tree 在高维正交范围查询上可能出现严重退化（特别是 $k=2(d-1)$ 增大、数据高度重叠、或分裂维度与分布不匹配时），可能出现：
$$
Q_{\mathrm{KD}}(n,k)=\Theta(n),
$$
从而 Sampling 的总时间在极端情况下可达 $\Theta(n^2)$ 量级（比如 Phase1/Phase3 中有 $\Theta(n)$ 次范围计数/遍历，每次成本 $\Theta(n)$）。这也是 KD‑tree 作为 baseline 的典型性质。

------

## （可选）附录：实现清单（强烈建议写在代码/论文补充材料里）

1. **事件排序必须稳定**；同 $x_1$ 时 END-before-START；START 再按 class 再按 id。
2. START 事件处理必须“先 query 再 activate”。
3. rank 化必须用二级键（value,id）并配 LOW/HIGH 哨兵，严格实现 `<`/`>`。
4. alias 构建时过滤 $w=0$。
5. KD.Sample 内部 pieces 不交并：`bbox⊆Q` 必须立即停止展开（否则会重复）。
6. 所有 sampling 都必须“有放回”，因此允许同 id 重复出现。
7. Adaptive 切换时必须丢弃已枚举 `AllPairs`，但**不能**丢弃 `w_e` 与 `W`。