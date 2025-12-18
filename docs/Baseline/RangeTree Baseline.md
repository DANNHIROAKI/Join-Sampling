# 纯 Range Tree 的 Spatial Join i.i.d. 均匀抽样统一算法（RT‑Sampling / RT‑Enumerate / Adaptive‑RT）— 升级补强版

## 记号与约定

- 维度：$d\ge 2$，且假设 $d$ 为常数。令
  $$
  m:=d-1,\qquad k:=2m=2(d-1).
  $$
  （因此 $k$ 亦为常数，但注意 $k$ 可能偏大，baseline 的常数因子会显著。）

- 两类半开轴对齐盒子集合：
  $$
  R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
  R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\},\quad n=n_1+n_2.
  $$
  每个盒子
  $$
  r=\prod_{i=1}^d [L_i(r),R_i(r)),\qquad L_i(r)<R_i(r).
  $$

- **统一输出顺序**：所有输出 pair 均写成 $(r_c,r_{\bar c})$。当 START 来自 $R_{\bar c}$ 时，必须交换顺序写入输出。

- 半开区间相交判定（严格不等号）：
  $$
  [a,b)\cap[c,d)\neq\varnothing \iff \max(a,c)<\min(b,d).
  $$

- 每个盒子有唯一 id：$\mathrm{id}(r)\in\{1,\dots,n\}$。id 仅用于打破相等与做 rank 边界处理。

------

## 1. 问题定义与分析

### 1.1 Spatial Join 结果与采样目标

我们只关心跨集合相交对：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
目标输出 $t$ 个样本 $Z_1,\dots,Z_t\in J$，要求 **i.i.d. 均匀有放回**：
$$
\Pr(Z_j=P)=\frac{1}{|J|}\ (\forall P\in J),\qquad Z_1,\dots,Z_t\ \text{相互独立}.
$$

------

### 1.2 第 1 维 plane sweep 与事件全序（保证“唯一归属”）

选择第 1 维 $x_1$ 做扫描。每个盒子 $r$ 产生两事件：

- $\texttt{START}(r)$ at $x_1=L_1(r)$
- $\texttt{END}(r)$ at $x_1=R_1(r)$

定义事件结构（便于写成严格全序）：
$$
e=(x,\ \text{type},\ \text{cls},\ \text{id}),
$$
其中 $\text{type}\in\{\texttt{END},\texttt{START}\}$，$\text{cls}\in\{c,\bar c\}$。

**事件排序（稳定，必须与半开边界一致）**：

1. 按 $x$ 升序；
2. 同 $x$：**END 在 START 前**；
3. 同 $x$ 且同为 START：固定类优先级（例如 $c$ 优先于 $\bar c$），再按 id 升序；
4. 同 $x$ 且同为 END：顺序不影响正确性（可按 id 升序固定化）。

扫描位置为 $x_0$ 时定义活跃集（由数据结构隐式维护）：
$$
\mathrm{Active}(R_c)=\{r_c\in R_c:\ x_0\in[L_1(r_c),R_1(r_c))\},\quad
\mathrm{Active}(R_{\bar c})\ \text{同理}.
$$
**关键观察（START 时第 1 维必然重叠）**：处理 $\texttt{START}(q)$ 时，$x_0=L_1(q)$。对任何对立活跃盒子 $r\in\mathrm{Active}(\text{Other})$，有 $x_0\in[L_1(r),R_1(r))$。由于半开+“END-before-START”，若仅贴边（例如 $R_1(r)=L_1(q)$）则 $r$ 会先 END 被移出，不会误计。于是：**在 START 时刻，只需判断维度 $2..d$ 相交，即可推出整盒相交**。 

------

### 1.3 START 事件诱导的块分解：$J=\biguplus_e J_e$

令 $E$ 为所有 START 事件集合。对 START 事件 $e$（新盒子 $q=r(e)$），令“另一侧”集合为 $\text{Other}(e)\in\{R_c,R_{\bar c}\}\setminus\{\text{cls}(q)\}$。定义：
$$
K_e=\{r\in \mathrm{Active}(\text{Other}(e))\mid q\cap r\neq\varnothing\},\qquad w_e:=|K_e|.
$$
并定义局部相交对集合（输出统一为 $(r_c,r_{\bar c})$）：

- 若 $q\in R_c$：$J_e=\{(q,r_{\bar c})\mid r_{\bar c}\in K_e\}$
- 若 $q\in R_{\bar c}$：$J_e=\{(r_c,q)\mid r_c\in K_e\}$

**命题（事件分块不交分解）**
$$
J=\biguplus_{e\in E} J_e,\qquad |J|=\sum_{e\in E}w_e.
$$
**证明要点（写到论文也够用）**：任取 $(r_c,r_{\bar c})\in J$。考虑二者 START 事件在上述全序中的较晚者 $e^\star$。在处理 $e^\star$ 时，另一盒必已 START 并尚未 END（否则第 1 维不可能相交），故它在对立活跃集中；又由于整盒相交，因此维度 $2..d$ 也相交，于是该对属于 $J_{e^\star}$。而在较早 START 时刻，较晚者尚未插入活跃结构，因此该对不可能更早被输出。故每对恰归入唯一 $J_e$，得到不交并与计数恒等式。

> 这条分解是 RT‑Sampling / Adaptive‑RT 中“按事件权重抽样”的核心。 

------

### 1.4 三个版本的定位（同一事件系统 + 同一 DS）

- **RT‑Enumerate（Enumerate+Sampling）**：一次 sweep 用 `Report` 显式枚举全部 $J$ 到数组，再数组均匀采样。直观但依赖 $|J|$。
- **RT‑Sampling（Sampling，不枚举 $J$）**：两遍 sweep：第一遍只 Count 得到全部 $w_e$ 与 $W=|J|$；第二遍按 slot 批量调用 `Sample` 回填答案。最坏时间/空间不依赖 $|J|$。
- **Adaptive‑RT（Adaptive+Sampling）**：先尝试枚举，累计超过阈值 $J_\star$ 立即切换到 “Count‑only + RT‑Sampling”，兼顾小 $|J|$ 的常数优势与大 $|J|$ 的稳健性。 

------

## 2. 核心数据结构

RT baseline 的核心是：在每个 START 时刻，我们需要对“对立活跃集”做一次“维度 $2..d$ 相交”的判定/计数/采样。Range Tree 做法是把该条件一次性写成 **$k=2(d-1)$ 维正交范围查询**。

------

### 2.1 从盒子相交到 $k$ 维正交范围查询

对任意盒子 $r$，取其在维度 $2..d$ 的端点向量：
$$
\mathbf{L}(r)=(L_2(r),\dots,L_d(r))\in\mathbb{R}^m,\quad
\mathbf{R}(r)=(R_2(r),\dots,R_d(r))\in\mathbb{R}^m.
$$
定义嵌入点：
$$
p(r)=(\mathbf{L}(r),\mathbf{R}(r))=
(L_2(r),\dots,L_d(r),R_2(r),\dots,R_d(r))\in\mathbb{R}^{2m}.
$$
对查询盒子 $q$，维度 $i\ge2$ 的相交条件为：
$$
[L_i(r),R_i(r))\cap[L_i(q),R_i(q))\neq\varnothing
\iff
L_i(r)<R_i(q)\ \land\ R_i(r)>L_i(q).
$$
因此定义正交范围（带无穷边界、严格不等号）：
$$
\mathcal{Q}(q)=
\prod_{i=2}^{d}(-\infty, R_i(q))
\ \times\
\prod_{i=2}^{d}(L_i(q), +\infty).
$$
**引理（等价性）**
$$
r\text{ 与 }q\text{ 在维度 }2..d\text{ 相交}\iff p(r)\in \mathcal{Q}(q).
$$
结合 §1.2 的关键观察：在处理 START(q) 时，第 1 维已保证重叠，因此只需判断 $p(r)\in\mathcal{Q}(q)$ 即可决定整盒相交。 

------

### 2.2 严格不等号的可实现化：二级键 + rank 区间（补强到可直接实现）

难点：$\mathcal{Q}(q)$ 中包含严格不等号 “$<$” 与 “$>$”。为了在离散 rank 域上稳定实现，采用**二级键**：
$$
\kappa(x,r):=(x,\mathrm{id}(r))
$$
并按字典序排序。因为 id 唯一，所以对同一轴上，所有点坐标 $(\text{value},\text{id})$ 都互异。

为处理严格边界，引入两个哨兵 id：
$$
\text{min sentinel }0 < 1,\dots,n,\qquad
\text{max sentinel }n < n+1.
$$
下面给出两条可直接落地的等价关系（这是全文最关键的“off‑by‑one 防雷”补强）：

**性质 A（严格小于）**
 对任意真实坐标值 $x$，对任意点 $r$：
$$
L_i(r)<x
\iff
(L_i(r),\mathrm{id}(r)) < (x,0).
$$
**性质 B（严格大于）**
 对任意真实坐标值 $x$，对任意点 $r$：
$$
R_i(r)>x
\iff
(R_i(r),\mathrm{id}(r)) > (x,n).
$$

> 直觉：$(x,0)$ 比任何 $(x,\mathrm{id})$ 都小；$(x,n)$ 比任何 $(x,\mathrm{id})$ 都大。

------

#### 2.2.1 预处理：为每个轴建立“二级键排序表”

对每个维 $i\in\{2,\dots,d\}$，分别构建两张有序表：

- $A_i$：所有 $(L_i(r),\mathrm{id}(r))$ 的升序数组（长度 $n$）
- $B_i$：所有 $(R_i(r),\mathrm{id}(r))$ 的升序数组（长度 $n$）

并记：

- $\mathrm{rankA}_i(r)$ 为 $(L_i(r),\mathrm{id}(r))$ 在 $A_i$ 中的位置（0-index）
- $\mathrm{rankB}_i(r)$ 为 $(R_i(r),\mathrm{id}(r))$ 在 $B_i$ 中的位置（0-index）

于是每个盒子映射为 rank 点：
$$
p^\#(r)=\bigl(\mathrm{rankA}_2(r),\dots,\mathrm{rankA}_d(r),\mathrm{rankB}_2(r),\dots,\mathrm{rankB}_d(r)\bigr)\in[0,n)^{2m}.
$$

------

#### 2.2.2 查询：把 $\mathcal{Q}(q)$ 变成 rank 空间中的半开区间笛卡尔积

对每个 $i\ge2$：

- 约束 $L_i(r)<R_i(q)$
   令
  $$
  \mathrm{hiA}_i(q):=\mathrm{lower\_bound}\bigl(A_i,(R_i(q),0)\bigr)\in[0,n],
  $$
  则满足该约束的 rank 为：
  $$
  \mathrm{rankA}_i(r)\in [0,\ \mathrm{hiA}_i(q)).
  $$

- 约束 $R_i(r)>L_i(q)$
   令
  $$
  \mathrm{loB}_i(q):=\mathrm{upper\_bound}\bigl(B_i,(L_i(q),n)\bigr)\in[0,n],
  $$
  则满足该约束的 rank 为：
  $$
  \mathrm{rankB}_i(r)\in [\ \mathrm{loB}_i(q),\ n).
  $$

综上，rank 空间中的查询盒为：
$$
\mathcal{Q}^\#(q)=
\prod_{i=2}^{d}[0,\mathrm{hiA}_i(q))
\ \times\
\prod_{i=2}^{d}[\mathrm{loB}_i(q),n).
$$
**引理（严格不等号的 rank 等价性）**
$$
p(r)\in\mathcal{Q}(q)\iff p^\#(r)\in\mathcal{Q}^\#(q).
$$
**证明要点**：逐维用性质 A/B 与二分边界定义即可；“$<$” 对应 lower_bound，“$>$” 对应 upper_bound，且统一用半开区间 $[\,\cdot,\cdot\,)$ 表示，避免歧义。

> 这一节补强的目的：让 reviewer 无法用 off‑by‑one 或 “严格不等号怎么做” 卡你。 

------

### 2.3 动态 $k$ 维 Range Tree（作为黑盒 DS 的接口与假设）

维护两棵动态 $k$-维 Range Tree：

- $\mathcal{RT}^{(c)}$：维护当前 $\mathrm{Active}(R_c)$ 的点集 $\{p^\#(r_c)\}$
- $\mathcal{RT}^{(\bar c)}$：维护当前 $\mathrm{Active}(R_{\bar c})$ 的点集 $\{p^\#(r_{\bar c})\}$

每个点携带 box-id（或指针）以便从点回溯盒子对象。

#### 2.3.1 需要支持的接口（统一定义）

对任一 $\mathcal{RT}$，支持：

- `Insert(r)` / `Delete(r)`：动态维护活跃点集
- `Count(Q#)`：返回 $\#\{r:\ p^\#(r)\in Q^\#\}$
- `Report(Q#)`：枚举所有 $r$ 使 $p^\#(r)\in Q^\#$
- `Sample(Q#, t')`：从集合 $\{r:\ p^\#(r)\in Q^\#\}$ 中 **i.i.d. 均匀有放回**采样 $t'$ 次（返回长度为 $t'$ 的列表）

#### 2.3.2 复杂度“baseline 假设”（写清前提，避免被挑刺）

在常数维 $k$ 下，假设动态 Range Tree（含关联结构）可达：
$$
\texttt{Insert/Delete/Count}=O((\log n)^k),\quad
\texttt{Report}=O((\log n)^k+k_{\text{out}}),\quad
\text{Space}=O(n(\log n)^{k-1}).
$$
并且我们将进一步说明如何在该框架上实现 `Sample(Q#,t')` 达到：
$$
\texttt{Sample}(Q^\#,t')=O((\log n)^k+t').
$$

> 注：这里是 baseline 的“理论 DS 假设”；论文/报告里建议加一句：实现可用递归关联结构 + 每层平衡 BST/segment tree 完成；常数维下该假设是标准范围树量级。 

------

### 2.4 在 Range Tree 上实现 `Sample(Q#, t')`：规范分解 + 子块加权 + slot 批量化（重点补强）

这一节把“Range Tree 原生无采样接口”的疑点彻底补上：给出**可复现的抽样构造**与**为什么仍 i.i.d.**。

#### 2.4.1 规范分解（canonical decomposition）回顾

在最外层坐标轴（第 1 轴）上，对任意一维区间 $[a,b)$，Range Tree（或等价的平衡 segment tree）会返回一组节点集合：
$$
\mathcal{C}=\{v_1,\dots,v_s\},\quad s=O(\log n),
$$
满足：

- 各节点代表的子树点集两两不交；
- 它们的并集恰好覆盖所有满足该轴约束的点；
- 对每个节点 $v$，其关联结构 $\mathcal{RT}_v$ 是一个 $(k-1)$ 维结构，维护该节点子树点集在剩余坐标上的投影（动态版本下同样成立，只是更新时会维护这些关联结构）。

对半无限区间（如 $[0,h)$、$[l,n)$）也同理，仍是 $O(\log n)$ 个规范节点。

------

#### 2.4.2 单次采样（1 个样本）：按子块大小加权的递归

给定 $k$ 维查询盒 $Q^\#$。把它写成：
$$
Q^\# = I \times Q^\#_{\text{tail}},
$$
其中 $I$ 是最外层轴的一维区间，$Q^\#_{\text{tail}}$ 是剩余 $k-1$ 轴的笛卡尔积约束。

令 $\mathcal{C}$ 是 $I$ 的规范分解节点集合。对每个 $v\in\mathcal{C}$，定义
$$
w_v := \texttt{Count}_{\mathcal{RT}_v}(Q^\#_{\text{tail}}),\qquad
W:=\sum_{v\in\mathcal{C}}w_v.
$$
**单次 `Sample` 过程**：

1. 以概率 $w_v/W$ 选择一个节点 $v\in\mathcal{C}$；
2. 在 $\mathcal{RT}_v$ 上递归调用 `Sample(Q_tail#, 1)`，得到一个点（最终回到盒子 id）。

**为什么均匀**：任何候选点 $x\in Q^\#$ 只属于唯一一个节点 $v(x)\in\mathcal{C}$。于是
$$
\Pr(\text{输出 }x)=\Pr(v=v(x))\cdot \Pr(x\mid v=v(x))
=\frac{w_{v(x)}}{W}\cdot \frac{1}{w_{v(x)}}=\frac{1}{W}.
$$

------

#### 2.4.3 批量采样（$t'$ 个样本）：slot 聚合把 $t'(\log n)^k$ 降为 $(\log n)^k+t'$

若直接重复单次过程 $t'$ 次，每次都要重新做规范分解与多次 Count，会变成 $t'(\log n)^k$。为得到 $(\log n)^k+t'$，做法是：

1. **一次性计算权重**：计算所有 $w_v$（共 $O(\log n)$ 个），并建 alias（或前缀和）用于在 $\mathcal{C}$ 上按 $w_v/W$ 采样；

2. **先抽节点，再统计次数**：独立抽 $t'$ 次节点 $V_1,\dots,V_{t'}$，统计每个节点被抽到次数 $t_v$；

3. **每个节点只递归一次**：对每个 $v$ 调用一次
   $$
   \mathcal{RT}_v.\texttt{Sample}(Q^\#_{\text{tail}}, t_v),
   $$
   得到 $t_v$ 个 i.i.d. 样本，回填到对应的 slot 位置。

**i.i.d. 不会被破坏（关键论证）**：上述过程与“逐次独立地：先抽节点 v，再在该节点内均匀抽点”的生成过程完全等价；slot 只是把相同 $v$ 的试验集中执行，因此联合分布不变。

**复杂度（归纳界）**：记 $T_k(t')$ 为 $k$ 维 `Sample` 的时间。则
$$
T_k(t') = O((\log n)^k) + \sum_{v\in\mathcal{C}:t_v>0} T_{k-1}(t_v).
$$
而 $|\mathcal{C}|=O(\log n)$、$\sum_v t_v=t'$。用归纳假设 $T_{k-1}(t_v)=O((\log n)^{k-1}+t_v)$，得
$$
T_k(t') = O((\log n)^k) + O(|\mathcal{C}|(\log n)^{k-1}) + O(t')
= O((\log n)^k + t').
$$

> 这一节等价于把你原稿 §2.4 的想法补成“可直接写进论文的严谨版本”。 

------

### 2.5 Alias 表（事件级 / 规范节点级）实现细节（可选但建议写）

对任意离散分布 $\{p_i\}_{i=1}^M$，alias 构建 $O(M)$，单次采样 $O(1)$。本文中会用两次：

- **事件级 alias**：在 START 事件集合上按 $w_e/W$；
- **范围树节点级 alias**：在规范分解集合 $\mathcal{C}$ 上按 $w_v/W$。

工程建议：若 $M=O(\log n)$ 很小，也可用前缀和 + 二分（$O(\log\log n)$ 近似常数）；但 alias 写起来最统一。

------

## 3. 算法详细流程（三个版本）

### 3.0 公共预处理（所有版本共享）

输入：$R_c,R_{\bar c}$，样本数 $t$；Adaptive 还需阈值 $J_\star$。

**步骤 P0：分配 id**

- 为每个盒子分配唯一 $\mathrm{id}(r)\in[1..n]$。

**步骤 P1：构造事件数组**

- 每个盒子 $r$ 生成两条事件：
  - $(L_1(r),\texttt{START},\text{cls}(r),\mathrm{id}(r))$
  - $(R_1(r),\texttt{END},\text{cls}(r),\mathrm{id}(r))$

**步骤 P2：稳定排序**

- 按 §1.2 的事件排序规则排序，得到 `Events[1..2n]`。

**步骤 P3：二级键压缩（维度 2..d）**

- 对每个 $i\in\{2,\dots,d\}$：
  - 构建 $A_i=\mathrm{sort}\{(L_i(r),\mathrm{id}(r))\}$
  - 构建 $B_i=\mathrm{sort}\{(R_i(r),\mathrm{id}(r))\}$
- 为每个盒子 $r$ 预计算点坐标 $p^\#(r)\in[0,n)^{2m}$。

**步骤 P4：初始化两棵空 Range Tree**

- $\mathcal{RT}^{(c)}\leftarrow\emptyset$，$\mathcal{RT}^{(\bar c)}\leftarrow\emptyset$。

**辅助函数：构造查询盒 $\mathcal{Q}^\#(q)$**

- 对每个 $i\ge2$：

  - $\mathrm{hiA}_i(q)\leftarrow \mathrm{lower\_bound}(A_i,(R_i(q),0))$
  - $\mathrm{loB}_i(q)\leftarrow \mathrm{upper\_bound}(B_i,(L_i(q),n))$

- 返回
  $$
  Q^\#(q)=\prod_{i=2}^{d}[0,\mathrm{hiA}_i(q))\times\prod_{i=2}^{d}[\mathrm{loB}_i(q),n).
  $$

------

### 3.1 版本一：RT‑Enumerate（Enumerate+Sampling）

**目标**：显式构造全部 $J$，然后数组均匀采样得到 i.i.d.

**核心流程（一次 sweep + Report）**：

```
RT-Enumerate(R_c, R_barc, t):
  Pairs = empty dynamic array
  init RT^(c), RT^(barc) as empty

  for e in Events in order:
    if e is END(r):
      if r in R_c:    RT^(c).Delete(r)
      else:           RT^(barc).Delete(r)

    else e is START(q):
      Q = BuildQueryBox(q)  // Q#(q)
      if q in R_c:
        List = RT^(barc).Report(Q)
        for r in List: Pairs.push( (q, r) )
        RT^(c).Insert(q)
      else:  // q in R_barc
        List = RT^(c).Report(Q)
        for r in List: Pairs.push( (r, q) )
        RT^(barc).Insert(q)

  W = |Pairs|
  if W == 0: return empty

  Ans[1..t]
  for j=1..t:
    idx ~ Unif{0..W-1}
    Ans[j] = Pairs[idx]
  return Ans
```

**要点补强**：

- `Report(Q)` 返回的是满足 $p^\#(r)\in Q^\#(q)$ 的活跃点（即对立活跃盒子中在维度 $2..d$ 相交者），由 §1.2 可推出整盒相交。
- 不重不漏来自 §1.3 的“事件唯一归属”。

------

### 3.2 版本二：RT‑Sampling（不显式枚举 $J$）

三阶段：**Count → 事件抽样规划 → 第二遍 sweep 批量 Sample 回填**。

#### Phase 1：第一次扫描（只计数，得到 $w_e$ 与 $W=|J|$）

```
Phase1-Count():
  init RT^(c), RT^(barc) as empty
  for each START event e: w[e] = 0
  W = 0

  for event e in Events:
    if END(r):
      delete r from its side RT
    else START(q):
      Q = BuildQueryBox(q)
      if q in R_c:     w[e] = RT^(barc).Count(Q)
      else:            w[e] = RT^(c).Count(Q)
      W += w[e]
      insert q into its side RT

  return (W, w[·])
```

结束后由 §1.3 命题得 $W=\sum_e w_e=|J|$。若 $W=0$ 直接返回空。

------

#### Phase 2：事件级 alias + slot 分配

```
Phase2-Plan(t, W, w[·]):
  build alias over START events with w[e] > 0 using prob w[e]/W

  for each START event e: Slots[e] = empty list
  for j=1..t:
    e = AliasSampleEvent()
    Slots[e].push(j)

  return Slots[·]
```

记 $t_e:=|Slots[e]|$，则 $\sum_e t_e=t$。

------

#### Phase 3：第二次扫描（按 slot 调用 Sample 并回填）

```
Phase3-Generate(Slots[·]):
  clear RT^(c), RT^(barc)
  Ans[1..t]

  for event e in Events:
    if END(r):
      delete r from its side RT
    else START(q):
      Q = BuildQueryBox(q)
      te = |Slots[e]|
      if te > 0:
        if q in R_c:
          List = RT^(barc).Sample(Q, te)   // i.i.d uniform with replacement
          fill Ans at positions Slots[e] with (q, List[u])
        else:
          List = RT^(c).Sample(Q, te)
          fill Ans at positions Slots[e] with (List[u], q)

      insert q into its side RT

  return Ans
```

**为什么 Phase3 总复杂度是 $O(n(\log n)^k + t)$**：最多只有 $n$ 个 START 事件，每个事件最多调用一次 `Sample(Q, t_e)`，其代价 $O((\log n)^k + t_e)$。求和得 $O(n(\log n)^k + \sum_e t_e)=O(n(\log n)^k+t)$。

------

### 3.3 版本三：Adaptive‑RT（阈值自适应切换）

目标：小 $J$ 时省掉第二遍 sweep；大 $J$ 时不被 $|J|$ 拖垮。

#### Adaptive Phase1：一次扫描（永远 Count；必要时 Enumerate，超过阈值立即停）

```
Adaptive-Phase1(J_star):
  mode = ENUMERATE
  AllPairs = empty array     // only used if mode stays ENUMERATE
  init RT^(c), RT^(barc) empty
  for each START event e: w[e] = 0
  W = 0

  for event e in Events:
    if END(r):
      delete r from its side RT

    else START(q):
      Q = BuildQueryBox(q)

      // (1) Always Count
      if q in R_c:   w[e] = RT^(barc).Count(Q)
      else:          w[e] = RT^(c).Count(Q)
      W += w[e]

      // (2) Optionally Enumerate until threshold exceeded
      if mode == ENUMERATE:
        if W <= J_star:
          if q in R_c:
            List = RT^(barc).Report(Q)
            for r in List: AllPairs.push((q, r))
          else:
            List = RT^(c).Report(Q)
            for r in List: AllPairs.push((r, q))
        else:
          mode = COUNT_ONLY
          free(AllPairs)   // discard enumerated pairs

      // (3) Insert
      insert q into its side RT

  return (mode, W, w[·], AllPairs if kept)
```

#### 分支 A：未切换（$W=|J|\le J_\star$）

- `AllPairs` 已完整枚举 $J$，直接数组均匀采样 $t$ 次（同 RT‑Enumerate 的采样阶段）。

#### 分支 B：已切换（$W>|J_\star|$）

- 执行 RT‑Sampling 的 Phase2 + Phase3：
  - Phase2 用已保存的精确 $w[e]$ 做事件 alias + slot；
  - Phase3 第二遍 sweep 生成样本回填。

> 关键点：切换只影响“是否保留枚举结果”，不影响 `Count` 的精确性与事件序，因此不会引入偏差。 

------

## 4. Adaptive 阈值 $J_\star$ 的选择策略（更细、更可落地）

阈值决定：枚举到多大就停止并切换到 Sampling。建议用 **硬约束（内存）+ 软约束（时间）+ 标定（交叉点）** 三层策略。

### 4.1 内存硬约束（必须满足，避免 OOM）

设可用内存预算 `MemBudget`（字节），允许给 `AllPairs` 使用比例 $\rho\in(0,1)$。若每个 pair 存两条 32-bit id，则
$$
\mathrm{sizeof(Pair)}\approx 8\ \text{bytes};
$$
若用 64-bit id，则约 16 bytes（还不含 vector/allocator 开销）。

建议用保守估计 $\mathrm{sizeof(Pair)}\gets 16$ 或 $24$ bytes（包含容器开销）。则：
$$
J_\star^{\text{mem}}=\left\lfloor\frac{\rho\cdot \text{MemBudget}}{\text{sizeof(Pair)}}\right\rfloor,
\qquad
J_\star\le J_\star^{\text{mem}}.
$$
**工程建议**：`AllPairs` 用“结构体数组 + 预留容量 + 只存 id（不存大对象）”，可以显著降低真实常数。

------

### 4.2 时间权衡（理论建议，写清“为何这样选”）

发生切换时，Adaptive Phase1 最多额外枚举 $O(J_\star)$ 个 pair，之后只 Count。大分支总时间上界：
$$
T_{\text{Adaptive-B}}
=O\bigl(n(\log n)^k + J_\star + t\bigr).
$$
若希望大分支仍与 RT‑Sampling 同阶 $O(n(\log n)^k+t)$，应选：
$$
J_\star=O\bigl(n(\log n)^k+t\bigr).
$$
工程写法：
$$
J_\star^{\text{time}}=C_1\,n(\log n)^k+C_2\,t,\qquad
J_\star=\min(J_\star^{\text{mem}},J_\star^{\text{time}}).
$$

> 注：因为本 baseline 的 $k=2(d-1)$ 可能偏大，$(\log n)^k$ 的常数可能很高，很多情况下会倾向于取更小的 $J_\star$，让大分支尽快进入 Sampling 路径以避免枚举拖慢。

------

### 4.3 工程标定（推荐：用交叉点拟合常数）

1. 在代表性数据上跑：

   - RT‑Enumerate（完整枚举 + 数组采样）
   - RT‑Sampling（两遍扫描采样）

2. 找到两者耗时交叉点 $J_{\text{cross}}$（以真实 $|J|$ 为横轴）

3. 拟合：
   $$
   |J_{\text{cross}}|\approx C_1\,n(\log n)^k+C_2\,t
   $$

4. 运行时取：
   $$
   J_\star \approx 0.8\cdot |J_{\text{cross}}|
   $$
   并受 $J_\star^{\text{mem}}$ 截断。

------

## 5. 算法分析（正确性、复杂度；三版本均覆盖）

### 5.1 数据结构语义正确性（Count/Report/Sample 的含义）

在处理 START(q) 时，对立 Range Tree 存储的是对立活跃集（第 1 维已保证重叠）。由 §2.1–§2.2 的等价性可得：

- `Count(Q#(q))` 返回的正是
  $$
  w_e=|K_e|=|J_e|.
  $$

- `Report(Q#(q))` 枚举的正是 $K_e$，因此输出 pair 恰为 $J_e$。

- `Sample(Q#(q),t')` 若按 §2.4 的 canonical+slot 实现，则输出是 $K_e$ 上 **i.i.d. 均匀有放回**的 $t'$ 个样本。

------

### 5.2 RT‑Enumerate 的正确性

**定理 1（RT‑Enumerate 输出 i.i.d. 均匀）**

- 枚举阶段：由 §1.3 的事件唯一归属，每个 $(r_c,r_{\bar c})\in J$ 会且只会在其“较晚 START”时刻被 `Report` 输出一次，因此 `Pairs` 是 $J$ 的一一枚举（不重不漏）。
- 采样阶段：对 `Pairs[0..|J|-1]` 做独立均匀下标采样，得到 i.i.d. 均匀有放回样本。

------

### 5.3 RT‑Sampling 的正确性（全局 i.i.d. 均匀）

令 Phase1 得到的总权重：
$$
W=\sum_{e\in E} w_e = |J|.
$$
对任意 $P\in J$，由 §1.3 存在唯一事件 $e^\star$ 使 $P\in J_{e^\star}$，且 $|J_{e^\star}|=w_{e^\star}$。

对任意样本位置 $j$：

- Phase2：$\Pr(E_j=e)=w_e/W$

- Phase3：给定 $E_j=e^\star$，调用 `Sample(Q#(q),t_{e^\star})` 在 $J_{e^\star}$ 上均匀，因此
  $$
  \Pr(P\mid E_j=e^\star)=\frac{1}{w_{e^\star}}.
  $$

于是：
$$
\Pr(Z_j=P)=\frac{w_{e^\star}}{W}\cdot\frac{1}{w_{e^\star}}=\frac{1}{W}=\frac{1}{|J|}.
$$
**独立性（为什么 slot 不破坏 i.i.d.）**：Phase2 中 $\{E_j\}$ 独立生成；Phase3 对每个事件 $e$ 的 `Sample` 返回的是 i.i.d. 序列；slot 回填只是把“第 u 次样本”放回该事件对应的第 u 个 slot 位置，因此整体联合分布与逐个独立采样完全一致。

------

### 5.4 Adaptive‑RT 的正确性（切换不改变分布）

- 若未切换：`AllPairs` 完整枚举 $J$，输出等价于 RT‑Enumerate，因此正确。
- 若已切换：虽然丢弃了枚举结果，但 Phase1 始终保存了每个事件的精确 $w_e$（来自 `Count`），随后 Phase2+Phase3 与 RT‑Sampling 完全相同，因此正确。
- 切换不改变事件顺序、也不改变 Insert/Delete 时刻、也不改变 $w_e$ 的计算方式，因此不会引入任何偏差。

------

### 5.5 复杂度（时间与空间；三版本分别给出）

令 $k=2(d-1)$（常数）。

#### 5.5.1 RT‑Enumerate

- **时间**

  - 预处理：事件排序 + 二级键压缩 $O(n\log n)$
  - 扫描维护：每事件 Insert/Delete $O((\log n)^k)$，共 $O(n(\log n)^k)$
  - 每个 START 做 `Report`：$O((\log n)^k+k_e)$，且 $\sum_e k_e=|J|$
  - 数组采样：$O(t)$

  合并：
  $$
  T_{\text{Enum}}=O\bigl(n(\log n)^k+|J|+t\bigr).
  $$

- **空间**
  $$
  S_{\text{Enum}}=O\bigl(n(\log n)^{k-1}+|J|+t\bigr).
  $$

#### 5.5.2 RT‑Sampling

- **时间**

  - Phase1：$O(n(\log n)^k)$
  - Phase2：alias + slot：$O(n+t)$
  - Phase3：第二遍扫描维护 $O(n(\log n)^k)$ + 总采样 $O(n(\log n)^k+t)$

  合并：
  $$
  T_{\text{Sampling}}=O\bigl(n(\log n)^k+t\bigr).
  $$

- **空间**
  $$
  S_{\text{Sampling}}=O\bigl(n(\log n)^{k-1}+t\bigr)
  $$
  （另加 $O(n)$ 存事件与权重）。

#### 5.5.3 Adaptive‑RT

- **Case A：未切换（$|J|\le J_\star$）**
  $$
  T=O\bigl(n(\log n)^k+|J|+t\bigr),\quad
  S=O\bigl(n(\log n)^{k-1}+|J|+t\bigr).
  $$

- **Case B：已切换（$|J|>J_\star$）**
   Phase1 最多额外枚举 $O(J_\star)$ 后停止，因此：
  $$
  T=O\bigl(n(\log n)^k+J_\star+t\bigr),\quad
  S=O\bigl(n(\log n)^{k-1}+\max(J_\star,t)\bigr).
  $$
  若按 §4.2 选 $J_\star=O(n(\log n)^k+t)$，则大分支渐进与 RT‑Sampling 同阶。

------

## 你可以直接拿去用的“落地检查清单”（额外赠送，避免 implementation 被卡）

1. **事件排序**是否严格实现 END-before-START + START tie-break（类优先级 + id）？
2. **查询边界**是否严格用：
   - `hiA = lower_bound(A_i, (R_i(q),0))`
   - `loB = upper_bound(B_i, (L_i(q),n))`
3. `START(q)` 时是否先 Query/Count/Sample 再 Insert(q)？
4. `Sample` 是否真的实现了“先按子块权重抽块，再块内均匀”，并支持 batch slot？
5. Adaptive 切换时是否 **释放 AllPairs** 且不再 `Report`？