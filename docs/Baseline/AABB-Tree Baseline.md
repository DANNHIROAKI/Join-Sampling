# Plane Sweep + AABB-Tree：统一算法文档

本文件合并并统一了三个版本：**Sampling**、**Enumerate+Sampling（Baseline）**、**Adaptive+Sampling**。

---

## 1. 问题定义与分析

### 1.1 $d$ 维两类半开盒子与相交对

在 $d\ge 2$ 维欧氏空间 $\mathbb{R}^d$ 中，给定两类轴对齐半开盒子集合：
$$
R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\},\quad n=n_1+n_2.
$$
每个盒子为
$$
r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\qquad L_i(r)<R_i(r).
$$
只关心跨集合相交对：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
半开区间下，两盒子相交当且仅当逐维相交：
$$
r_c\cap r_{\bar c}\neq\varnothing
\iff
\forall i,\ \max(L_i(r_c),L_i(r_{\bar c}))<\min(R_i(r_c),R_i(r_{\bar c})).
$$

------

### 1.2 采样目标

输出 $t$ 个样本
$$
Z_1,\dots,Z_t\in J,\quad Z_j=(r_c,r_{\bar c}),
$$
要求 **i.i.d. 均匀（有放回）**：
$$
\Pr(Z_j=P)=\frac{1}{|J|}\ \ (\forall P\in J),\qquad
Z_1,\dots,Z_t\ \text{相互独立}.
$$

> 与“枚举 $J$ 再数组采样”的 Baseline 不同，这里我们仍采用“两遍扫描 + 中间 alias/slot”的结构来避免显式存下全部 $J$。
> 但由于 AABB Tree 属于工程型空间索引，其**最坏情况复杂度**不再能保证你原文那种 $O(n(\log n)^{d-1}+t)$ 的强界；本报告会给出**参数化复杂度**与**最坏退化警告**。

------

### 1.3 第 1 维 Plane Sweep 与事件分块

取第 1 维 $x_1$ 做扫描。每个盒子 $r$ 产生两个事件：

- $\texttt{START}(r)$ at $x_1=L_1(r)$
- $\texttt{END}(r)$ at $x_1=R_1(r)$

**事件排序（稳定）**：

1. 按 $x_1$ 升序；
2. 同 $x_1$：$\texttt{END}$ 在 $\texttt{START}$ 前（匹配半开边界）；
3. 同 $x_1$ 且都是 $\texttt{START}$：固定集合优先级（例如先处理 $R_c$ 的 START，再处理 $R_{\bar c}$ 的 START）以获得**唯一归属**。

扫描到位置 $x_0$ 时定义活跃集：
$$
\mathrm{Active}(R_c)=\{r_c\in R_c\mid x_0\in[L_1(r_c),R_1(r_c))\},\quad
\mathrm{Active}(R_{\bar c})\ \text{同理}.
$$
对每个 START 事件 $e$，令新盒子 $q=r(e)$。定义局部集合：

- 若 $q\in R_c$：
  $$
  K_e=\{r_{\bar c}\in \mathrm{Active}(R_{\bar c})\mid q\cap r_{\bar c}\neq\varnothing\},\quad
  J_e=\{(q,r_{\bar c})\mid r_{\bar c}\in K_e\},\quad w_e=|K_e|.
  $$

- 若 $q\in R_{\bar c}$（对称）：
  $$
  K_e=\{r_c\in \mathrm{Active}(R_c)\mid q\cap r_c\neq\varnothing\},\quad
  J_e=\{(r_c,q)\mid r_c\in K_e\},\quad w_e=|K_e|.
  $$

**命题（不交分解）**
 设 $E$ 为所有 START 事件集合，则
$$
J=\biguplus_{e\in E} J_e,\qquad |J|=\sum_{e\in E} w_e.
$$
直觉：任意相交对会在两者 START 中“更靠后”的那个 START 时刻第一次同时活跃，因此只归属到唯一一个 $J_e$。

------

### 1.4 用 AABB Tree 替换“高维模式结构”的关键点

当处理某个 START($q$) 时，对立活跃集里的盒子在 $x_1$ 维与 $q$ 必重叠（因为它们都包含 $x_0=L_1(q)$）。因此：
$$
q\cap r\neq\varnothing
\iff
\bigl(q^\star \cap r^\star\neq\varnothing\bigr),
$$
其中 $q^\star,r^\star$ 是把盒子投影到维度 $2..d$ 的 $(d-1)$ 维半开盒。

所以我们只需在维度 $2..d$ 上维护对立活跃集的动态索引，支持：

- **CountIntersect($q^\star$)**：精确返回 $|K_e|$
- **SampleIntersect($q^\star,t'$)**：在 $K_e$ 上精确生成 $t'$ 个 i.i.d. 均匀样本（有放回）

这两项由动态 AABB Tree（BVH）来完成。

------

---

## 2. 核心数据结构

令
$$
m:=d-1,\qquad
r^\star:=\prod_{k=1}^{m}[a_k(r),b_k(r)),\ \ a_k(r)=L_{k+1}(r),\ b_k(r)=R_{k+1}(r).
$$

### 2.1 事件数组与句柄

- `Events[1..2n]`：每项为 $(x_1,\mathrm{type}\in\{\texttt{START},\texttt{END}\},\mathrm{class}\in\{c,\bar c\},\mathrm{id})$
- 每个盒子有唯一 `id`；插入 AABB Tree 后得到一个 `handle`（指向叶节点或其位置），用于快速删除。

------

### 2.2 两棵动态 AABB Tree（维护对立活跃投影）

维护两棵动态 AABB Tree（BVH）：

- $\mathcal{T}^{(c)}$：存当前活跃的 $\{r_c^\star\}$
- $\mathcal{T}^{(\bar c)}$：存当前活跃的 $\{r_{\bar c}^\star\}$

每个 BVH 节点包含：

- `AABB(node)`：该节点子树所有叶条目的包围盒（仍是 $m$ 维 AABB）
- `size(node)`：子树叶子数量（用于整子树计数与子树内均匀抽样）
- 若为叶节点：保存条目 $(\mathrm{id}, r^\star)$

> 这里的 AABB Tree 可以是“动态 BVH”：插入时选择 sibling 并旋转/重构以保持较低代价；删除时按 handle 移除并做局部修复。实现细节可视为黑盒。

------

### 2.3 半开相交判定与包含判定

#### 2.3.1 半开相交（精确）

对两个 $m$ 维半开盒 $A=\prod_k[\alpha_k,\beta_k)$、$B=\prod_k[\gamma_k,\delta_k)$，定义：
$$
\texttt{Intersect}(A,B)\iff \forall k:\ \max(\alpha_k,\gamma_k)<\min(\beta_k,\delta_k).
$$
**所有叶子条目必须用该精确判定**，以避免“贴边”误判相交。

#### 2.3.2 半开包含（用于整子树剪枝）

对节点包围盒 $B$ 与查询盒 $Q$，定义：
$$
\texttt{Contained}(B,Q)\iff \forall k:\ L_k(B)\ge L_k(Q)\ \land\ R_k(B)\le R_k(Q).
$$
若 $\texttt{Contained}(B,Q)$ 成立，则该子树中任意叶盒都满足 $\subseteq B\subseteq Q$，因此**必然与 $Q$ 相交**。

------

### 2.4 CountIntersect：精确计数

递归函数 `Count(node, Q)`：

1. 若 `Intersect(AABB(node), Q)` 为假：返回 0
2. 若 `Contained(AABB(node), Q)` 为真：返回 `size(node)`
3. 若 `node` 是叶：返回 $1$（当且仅当叶条目与 $Q$ 半开相交，否则 0）
4. 否则：返回 `Count(left,Q)+Count(right,Q)`

输出即
$$
\mathrm{CountIntersect}(Q)=|\{r\in \mathrm{Active}:\ r^\star\cap Q\neq\varnothing\}|.
$$

------

### 2.5 SampleIntersect：在相交结果集上精确均匀 i.i.d. 采样

AABB Tree 默认只支持 Report，不支持“从 query 结果中均匀采样”。我们用“**Guide 计数 + 按计数加权随机下行**”实现严格均匀。

#### 2.5.1 BuildGuide：一次遍历构造命中计数结构

对查询 $Q$，构造一个“GuideNode”：

- 字段 `cnt`：该子树中与 $Q$ 相交的叶条目数
- 三种类型：
  - `FULL(node)`：若 `Contained(AABB(node),Q)`，则 `cnt=size(node)`，并标记为 FULL（无需再展开）
  - `LEAF(entry)`：叶节点且条目与 $Q$ 相交，则 `cnt=1`
  - `PARTIAL(children...)`：部分命中内部节点，递归构造子 guide，并令 `cnt=cnt_L+cnt_R`

这样相交结果集合被拆成若干互不重叠的“子块”（FULL 子树或叶命中），其叶集合不交且并起来恰好是所有命中叶。

#### 2.5.2 子树内均匀抽叶：SampleSubtreeUniform

对任意 AABB Tree 节点 `node`（不是 guide），用 `size` 做随机下降直到叶：

- 在内部节点：
  $$
  \Pr(\text{走左})=\frac{\mathrm{size(left)}}{\mathrm{size(node)}},\quad
  \Pr(\text{走右})=\frac{\mathrm{size(right)}}{\mathrm{size(node)}}.
  $$

- 直到叶返回其条目 id。

得到“在该子树所有叶上均匀”的样本。

#### 2.5.3 SampleOne：在命中集合上均匀抽 1 个

对 GuideRoot：

- 若类型 `FULL(node)`：返回 `SampleSubtreeUniform(node)`

- 若类型 `LEAF(entry)`：返回该 entry

- 若类型 `PARTIAL(children)`：在子 guide 中以概率
  $$
  \Pr(\text{选子 }i)=\frac{\mathrm{cnt}_i}{\sum_j \mathrm{cnt}_j}
  $$
  选择一个子 guide 并递归

**结论：**对任意命中叶条目 $x$，有
$$
\Pr(\texttt{SampleOne}(Q)=x)=\frac{1}{|K(Q)|},
$$
其中 $K(Q)$ 为所有与 $Q$ 相交的活跃条目集合。

#### 2.5.4 SampleIntersect(Q, t')

- 先 `GuideRoot ← BuildGuide(Q)` 得到 $W=\mathrm{cnt(GuideRoot)}=|K(Q)|$
- 重复 $t'$ 次、每次用独立随机数调用 `SampleOne(GuideRoot)`
   得到 $t'$ 个 i.i.d. 均匀样本（有放回）。

------

### 补充：ReportIntersect / IntersectReport（用于枚举）

对每棵树支持：

- `Insert(id, r^*) -> handle`
- `Delete(handle)`（或 `Delete(id)` 但最好用 handle）
- `IntersectReport(Q) -> list<id>`：返回所有叶条目 `id` 使得条目盒 `r^*` 与 `Q` 半开相交
   （内部节点先用 `AABB(node)` 做剪枝）

------

---

## 3. 算法详细流程

### 3.1 Sampling（两遍扫描，不枚举）

下面给出 **Plane Sweep + AABB Tree 的三阶段（两遍扫描）采样算法**。

### 3.1 预处理：构建事件数组

输入：$R_c,R_{\bar c}$，样本数 $t$。

1. 对每个盒子 $r$ 生成 $\texttt{START}(r)$、$\texttt{END}(r)$ 两事件
2. 按 1.3 的稳定排序规则排序
3. 给每个 START 事件分配唯一编号 $e\in\{1,\dots,|E|\}$（用于存储 $w_e$ 与 slot）

------

### 3.2 第一阶段：第一次扫描（仅计数）

初始化：
$$
\mathcal{T}^{(c)}\gets\emptyset,\quad \mathcal{T}^{(\bar c)}\gets\emptyset,\quad
W\gets 0.
$$
并准备数组 `w[e]`（START 事件权重）。

按事件顺序扫描：

- 若事件为 $\texttt{END}(r)$：从对应类的 AABB Tree 中删除 $r^\star$

- 若事件为 $\texttt{START}(q)$（其编号为 $e$）：

  - 若 $q\in R_c$：
    $$
    w_e \leftarrow \mathcal{T}^{(\bar c)}.\mathrm{CountIntersect}(q^\star),\quad W\leftarrow W+w_e,
    $$
    然后插入 $\mathcal{T}^{(c)}.\mathrm{Insert}(q^\star)$

  - 若 $q\in R_{\bar c}$ 对称处理（查询 $\mathcal{T}^{(c)}$，再插入 $\mathcal{T}^{(\bar c)}$）

扫描结束：
$$
W=\sum_{e\in E} w_e = |J|.
$$
若 $W=0$，则无相交对，直接返回空输出。

------

### 3.3 第二阶段：事件级 alias + slot 分配

在 START 事件集合 $E$ 上按权重 $w_e$ 建立分布：
$$
p_e=\frac{w_e}{W}.
$$
构造 alias 表，使得一次抽样 $O(1)$ 得到事件 $e\sim p_e$。

维护每个事件的 slot 列表 $S_e\subseteq\{1,\dots,t\}$，初始为空。

对每个样本位置 $j=1..t$：

1. 抽事件 $E_j\sim p_e$
2. 把 $j$ 追加到 $S_{E_j}$

记 $t_e:=|S_e|$，则 $\sum_e t_e=t$。

------

### 3.4 第三阶段：第二次扫描（局部批量采样 + slot 回填）

重置两棵 AABB Tree 为空：
$$
\mathcal{T}^{(c)},\mathcal{T}^{(\bar c)}\gets \emptyset,
$$
准备输出数组 `Ans[1..t]`。

再次按同样事件顺序扫描：

- $\texttt{END}(r)$：从对应类树中删除

- $\texttt{START}(q)$（编号 $e$）：

  - 令 $t_e=|S_e|$

  - 若 $t_e>0$：

    - 若 $q\in R_c$：
      $$
      \text{List}\leftarrow \mathcal{T}^{(\bar c)}.\mathrm{SampleIntersect}(q^\star,t_e)
      $$
      按 slot 顺序 $S_e=\{j_1,\dots,j_{t_e}\}$ 回填：
      $$
      \mathrm{Ans}[j_u]\leftarrow (q,\ \text{List}[u]).
      $$

    - 若 $q\in R_{\bar c}$ 对称：
      $$
      \text{List}\leftarrow \mathcal{T}^{(c)}.\mathrm{SampleIntersect}(q^\star,t_e),\quad
      \mathrm{Ans}[j_u]\leftarrow (\text{List}[u],\ q).
      $$
      （始终保持输出 pair 顺序为 $(r_c,r_{\bar c})$。）

  - 最后将 $q$ 插入本侧 AABB Tree（变为活跃）

输出：
$$
(Z_1,\dots,Z_t):=(\mathrm{Ans}[1],\dots,\mathrm{Ans}[t]).
$$

------

### 3.2 Enumerate+Sampling（Baseline：一遍扫描枚举 + 数组采样）

### 3.1 总览

算法两阶段：

1. **一次 Plane Sweep：枚举并存储全部 $J$** 到数组 `Pairs`
2. **在 `Pairs` 上做均匀 i.i.d. 采样 $t$ 次** 输出 `Ans[1..t]`

------

### 3.2 预处理

输入：$R_c,R_{\bar c}$，样本数 $t$。

1. 构造 `Events`：每个盒子 $r$ 加入 $\texttt{START}(r)$、$\texttt{END}(r)$

2. 稳定排序 `Events`（规则见 1.3）k=2,d=d(Enumerat+Sampling)

3. 初始化两棵空 AABB Tree：
   $$
   \mathcal{T}^{(c)}\leftarrow\emptyset,\quad \mathcal{T}^{(\bar c)}\leftarrow\emptyset
   $$

4. 初始化动态数组 `Pairs=[]`

------

### 3.3 一次扫描：枚举 $J$

按 `Events` 顺序处理事件：

#### 情况 1：$\texttt{END}(r)$

- 若 $r\in R_c$：$\mathcal{T}^{(c)}.\texttt{Delete}(\text{handle}[r])$
- 若 $r\in R_{\bar c}$：$\mathcal{T}^{(\bar c)}.\texttt{Delete}(\text{handle}[r])$

> 因为同坐标 END 在 START 前，半开区间下“贴边不相交”的 pair 不会被枚举出来。k=2,d=d(Enumerat+Sampling)

#### 情况 2：$\texttt{START}(q)$

不失一般性先写 $q\in R_c$（对称情况同理）：

1. 计算查询投影：
   $$
   q^\star=\prod_{i=2}^{d}[L_i(q),R_i(q))
   $$

2. 在对立活跃树上做相交报告：
   $$
   \texttt{CandIDs} \leftarrow \mathcal{T}^{(\bar c)}.\texttt{IntersectReport}(q^\star)
   $$

3. 对每个返回的候选 `id = r_{\bar c}`：

   - 用 $\texttt{Intersect}^{\star}(q^\star, r_{\bar c}^\star)$ 精确过滤（半开语义）
   - 若为真，则 `Pairs.append( (q, r_{\bar c}) )`

4. 将 $q$ 插入本侧树并记录句柄：
   $$
   \text{handle}[q]\leftarrow \mathcal{T}^{(c)}.\texttt{Insert}(q^\star)
   $$

若 $q\in R_{\bar c}$，对称处理：

- `CandIDs ← T^(c).IntersectReport(q^*)`
- 过滤后追加 `(r_c, q)`（保持输出顺序统一为 $(r_c,r_{\bar c})$）
- `handle[q] ← T^(barc).Insert(q^*)`

扫描结束后：
$$
W:=|\texttt{Pairs}|=|J|.
$$

------

### 3.4 在 `Pairs` 上均匀采样 $t$ 次（有放回）

- 若 $W=0$：返回空（或返回空数组）
- 否则对 $j=1..t$：
  - 生成 $\text{idx}_j\sim \mathrm{Unif}\{0,\dots,W-1\}$
  - 输出 $\texttt{Ans}[j]=\texttt{Pairs}[\text{idx}_j]$

------

### 3.3 Adaptive+Sampling（阈值自适应：小规模枚举，大规模两遍采样）

整体仍是 **Adaptive**：Phase1 一遍扫描（永远 Count，必要时 Report）；根据是否超过阈值 $J_\star$ 选择分支；大分支再做 Phase2+Phase3。 k=2,d=d(Adaptive+Sampling)

### 3.1 Adaptive Phase1：一次扫描（精确计数 + 可选枚举）

**输入**：$R_c,R_{\bar c}$，样本数 $t$，阈值 $J_\star$

**初始化**：

- `mode = ENUMERATE`
- `AllPairs = []`（容量上限 $J_\star$）
- 两棵树清空：$\mathcal T^{(c)}\gets\emptyset,\ \mathcal T^{(\bar c)}\gets\emptyset$
- `W = 0`（累计 $|J|$）
- 为每个 START 事件准备存储 $w_e$

**扫描事件 `Events`**：

#### (1) 处理 END 事件

- 若 `END(r)` 且 $r\in R_c$：$\mathcal T^{(c)}.\texttt{Delete}(r)$
- 若 $r\in R_{\bar c}$：$\mathcal T^{(\bar c)}.\texttt{Delete}(r)$

#### (2) 处理 START 事件 $e$，新盒子 $q=r(e)$

设 $Q=q^\star$ 为 $q$ 在维度 $2..d$ 的投影。

**Step A：先精确计数（无论是否枚举都必须做）**

- 若 $q\in R_c$：
  $$
  w_e\leftarrow \mathcal T^{(\bar c)}.\texttt{CountIntersect}(Q)
  $$

- 若 $q\in R_{\bar c}$：
  $$
  w_e\leftarrow \mathcal T^{(c)}.\texttt{CountIntersect}(Q)
  $$

更新：
$$
W\leftarrow W+w_e,
$$
并保存 $w_e$。

**Step B：是否枚举并写入 `AllPairs`**

- 若 `mode==ENUMERATE` 且 $W\le J_\star$：允许枚举
  - 若 $q\in R_c$：`List ← T^(barc).ReportIntersect(Q)`；对每个 $r_{\bar c}\in\text{List}$，追加 $(q,r_{\bar c})$
  - 若 $q\in R_{\bar c}$：`List ← T^(c).ReportIntersect(Q)`；对每个 $r_c\in\text{List}$，追加 $(r_c,q)$
- 若 `mode==ENUMERATE` 但 $W>J_\star$：触发切换
  - `mode ← COUNT_ONLY`
  - 丢弃/释放 `AllPairs`
  - 后续 START 不再调用 `ReportIntersect`
- 若 `mode==COUNT_ONLY`：跳过枚举

> **先 Count 再判断阈值** 的设计非常关键：避免“一个事件 report 到一半才超阈值”，保证 Phase1 枚举量严格 $\le J_\star$。

**Step C：插入新盒子到本侧活跃树**

- 若 $q\in R_c$：$\mathcal T^{(c)}.\texttt{Insert}(q)$
- 若 $q\in R_{\bar c}$：$\mathcal T^{(\bar c)}.\texttt{Insert}(q)$

扫描结束得到：

- 精确 $|J|=W=\sum_e w_e$
- 每个 START 的 $w_e$
- 若未切换，则 `AllPairs` 完整等于 $J$

### 3.2 分支 A：Baseline（小 $|J|$）——数组采样

若最终 `mode==ENUMERATE`（即 $W\le J_\star$ 且未丢弃），则：

- `AllPairs` 是 $J$ 的完整枚举数组
- 对 $j=1..t$：独立采样 `idx ~ UniformInt(0, |AllPairs|-1)`，输出 `AllPairs[idx]`

### 3.3 分支 B：Ours（大 $|J|$）——Phase2 + Phase3（事件级，无模式）

若发生切换（`COUNT_ONLY`），继续：

#### 3.3.1 Phase2：事件 alias + slot 分配

- 若 $W=0$：返回空
- 在 START 事件集合 $E$ 上按 $p_e=w_e/W$ 建 alias
- 对 $j=1..t$：
  - 抽事件 $E_j\sim p$
  - 把位置 $j$ 放入 slot：$S_{E_j}\gets S_{E_j}\cup\{j\}$

#### 3.3.2 Phase3：第二次扫描 + 局部批量采样 + slot 回填

- 清空两棵树（活跃集置空）
- 初始化输出 `Ans[1..t]`
- 再扫一遍 `Events`：

对 START 事件 $e$，新盒 $q=r(e)$，投影 $Q=q^\star$：

- 令 $t_e=|S_e|$
- 若 $t_e>0$：
  - 若 $q\in R_c$：`List ← T^(barc).SampleIntersect(Q, t_e)`
  - 若 $q\in R_{\bar c}$：`List ← T^(c).SampleIntersect(Q, t_e)`
- 按 slot 顺序回填（保持输出为 $(r_c,r_{\bar c})$）：
  - 若 $q\in R_c$：对 $u=1..t_e$，令 `Ans[ S_e[u] ] = (q, List[u])`
  - 若 $q\in R_{\bar c}$：`Ans[ S_e[u] ] = (List[u], q)`
- 最后把 $q$ 插入本侧树

输出 `Ans[1..t]`。

------

---

## 4. Adaptive 版本的阈值选择策略

阈值 $J_\star$ 决定“枚举到什么规模就放弃 Baseline 并切到两遍采样”。在 AABB-tree 版本里（和 R-tree 类似），**复杂度更依赖数据分布**，因此建议同时考虑内存硬约束 + 经验时间交叉点。

### 4.1 内存硬约束（必须满足）

设内存预算 `MemBudget`，给 `AllPairs` 的比例 $\rho\in(0,1)$，每个 pair 存储开销 `sizeof(Pair)`（两个 id/指针）：
$$
J_\star^{\text{mem}}
=
\left\lfloor
\frac{\text{MemBudget}\cdot\rho}{\text{sizeof(Pair)}}
\right\rfloor.
$$
必须取
$$
J_\star\le J_\star^{\text{mem}}.
$$

### 4.2 时间权衡（建议用“可测的交叉点模型”）

在大分支中，Phase1 最多额外枚举 $J_\star$ 个 pair（之后停止 `ReportIntersect`），因此“大分支的额外输出开销”大约是 $O(J_\star)$。

但 AABB-tree 的 `CountIntersect / SampleIntersect` 成本与分布相关，难以给漂亮的 $\log^{d-1}n$ 最坏界，因此更建议：

- 用 benchmark 找到 **Baseline（全枚举）** 与 **Ours 分支（不枚举，两遍采样）** 的耗时交叉点 $|J_{\text{cross}}|$

- 设 $J_\star \approx 0.5\sim 0.8\cdot |J_{\text{cross}}|$

- 再用 $J_\star^{\text{mem}}$ 截断：
  $$
  J_\star=\min\bigl(J_\star^{\text{mem}},\ 0.8|J_{\text{cross}}|\bigr).
  $$

### 4.3 工程细节建议

- **必须先 Count 再判断是否超阈值**，这样保证枚举量 $\le J_\star$，逻辑干净且可控（你原稿里也强调了这一点）。 k=2,d=d(Adaptive+Sampling)
- 如果 AABB-tree 使用“fat AABB”，只要最终叶子用半开精确判定过滤，正确性仍成立；但“full-accept（node ⊆ Q）”的优化在 fat AABB 下仍是安全的（Q 包含 fat，必包含真实盒）。

------

---

## 5. 算法分析（正确性与复杂度）

### 5.1 Sampling 版本分析

### 4.1 正确性

#### 4.1.1 分块正确：$J=\biguplus_e J_e$

由半开边界与事件排序（END-before-START + START tie-break）可保证：任意相交对在两者 START 中“更靠后”的那一刻第一次同时活跃，因此只出现在唯一的 $J_e$ 中，不重不漏。

#### 4.1.2 Phase1 计数正确：$w_e=|J_e|$

处理 START($q$) 时，对立 AABB Tree 中保存的条目恰为对立活跃集在维度 $2..d$ 的投影。又因为活跃性保证第 1 维必相交，所以
$$
|J_e|=|\{r\in \mathrm{Active}(\text{Other}) : r^\star\cap q^\star\neq\varnothing\}|
$$
而 `CountIntersect(q^\star)` 在叶子用半开精确判定、内部用安全剪枝与“全包含子树”规则返回精确命中数，因此返回值就是 $w_e$。

#### 4.1.3 局部采样正确：SampleIntersect 在 $J_e$ 上均匀 i.i.d.

固定事件 $e$ 与查询 $q^\star$。`BuildGuide` 把命中叶集合 $K_e$ 表示为若干互不交叠的块（FULL 子树或叶命中），每块权重等于其叶数。`SampleOne` 以“块权重/总权重”选块，再在块内均匀选叶，因此任意叶条目概率为 $1/|K_e|$。重复 $t_e$ 次且随机数独立，则得到 i.i.d.（有放回）。

因此对 $J_e=\{(q,r): r\in K_e\}$，局部输出的 pair 也在 $J_e$ 上均匀 i.i.d.

#### 4.1.4 全局均匀性

对任意 $P\in J$，存在唯一事件 $e^*$ 使 $P\in J_{e^*}$，且 $|J_{e^*}|=w_{e^*}$。对任意位置 $j$：
$$
\Pr(Z_j=P)=\Pr(E_j=e^*)\cdot \Pr(P\mid E_j=e^*)
=\frac{w_{e^*}}{W}\cdot \frac{1}{w_{e^*}}
=\frac{1}{|J|}.
$$

#### 4.1.5 独立性

- Phase2 中 $(E_j)_{j=1}^t$ 独立生成；
- Phase3 中每个事件 $e$ 内部用独立随机数生成 $t_e$ 个 i.i.d. 样本；
- slot 回填是确定性的索引写入，不引入额外相关性。

因此 $(Z_1,\dots,Z_t)$ 是在 $J$ 上的 i.i.d. 均匀样本序列。

------

### 4.2 复杂度（参数化 + 最坏退化警告）

AABB Tree 的复杂度强依赖数据分布与树平衡质量，下面用“访问节点数”刻画一次 query 的代价。

记：

- $h$：AABB Tree 高度（工程上通常接近 $O(\log n)$，但最坏可能退化）
- 对某次查询盒 $Q$：
  - $V(Q)$：在 `CountIntersect/BuildGuide` 中访问到的节点数（其节点 AABB 与 $Q$ 相交且未被剪掉）
  - $E(Q)$：被检查的叶条目数（叶过滤精确判定的次数）

则：

- `CountIntersect(Q)`：时间 $O(V(Q)+E(Q))$

- `BuildGuide(Q)`：时间 $O(V(Q)+E(Q))$

- `SampleIntersect(Q,t')`：先建 guide，再采样：
  $$
  O\bigl(V(Q)+E(Q)+t'\cdot h\bigr)
  $$

#### Phase1（第一次扫描，仅计数）

$$
T_1 = O(n\log n)\ +\ O(\text{Insert/Delete 总成本})\ +\ \sum_{e\in E} O\bigl(V(q_e^\star)+E(q_e^\star)\bigr).
$$

#### Phase2（alias + slot）

$$
T_2 = O(n+t).
$$

#### Phase3（第二次扫描，采样回填）

令事件 $e$ 的采样需求为 $t_e=|S_e|$，则：
$$
T_3 = O(\text{Insert/Delete 总成本})
+\sum_{e\in E: t_e>0} O\bigl(V(q_e^\star)+E(q_e^\star)+t_e\cdot h\bigr).
$$
由于 $\sum_e t_e=t$，采样下降项汇总为 $O(t\cdot h)$。

#### 合并总时间

$$
T = O(n\log n)\ +\ O(\text{总更新})\ +\ \sum_{e\in E}O(V+E)\ +\ \sum_{e:t_e>0}O(V+E)\ +\ O(t\cdot h).
$$

- **常见（剪枝有效、树较平衡）情形**：可近似认为 $h=O(\log n)$，且多数查询的 $V+E$ 也接近 $O(\log n)$ 或 $O(\log n + \text{output})$，则整体常见表现接近
  $$
  T \approx O((n+t)\log n).
  $$

- **最坏情况警告（必须写清）**：若大量 AABB 重叠导致剪枝失败，则可能出现 $V(Q)=\Theta(n)$，从而两遍扫描在理论上可能退化到 $\Theta(n^2)$ 级别。

#### 空间复杂度

- 两棵 AABB Tree：$O(n)$ 节点（常数因子与实现有关）
- 事件数组、权重、alias：$O(n)$
- slot 与输出：$O(t)$
- guide 是“按查询临时构建”，峰值 $O(\max_Q V(Q))\le O(n)$

因此峰值空间：
$$
S = O(n+t).
$$

------

### 一句话总结

> **Plane Sweep + AABB Tree（两遍扫描）**：用 $x_1$ 扫线把全体相交对按 START 事件分块；维度 $2..d$ 上用两棵动态 AABB Tree 维护对立活跃投影，Phase1 精确计数得到 $w_e$，Phase2 用事件级 alias 分配 slot，Phase3 再扫一遍并用 “BuildGuide + 按命中计数加权随机下行” 在每个局部块上精确生成 i.i.d. 均匀样本，最终得到全局 $J$ 上的 i.i.d. 均匀采样序列。
> 正确性严格；复杂度为工程/数据相关，最坏可能退化。

### 5.2 Enumerate+Sampling（Baseline）版本分析

### 4.1 正确性

#### 4.1.1 枚举不重不漏：`Pairs` 恰好等于 $J$

**不漏：** 取任意 $(r_c,r_{\bar c})\in J$。在第 1 维上区间相交，因此两者 START 事件中存在一个在事件序里更靠后，记为 $e^*$ 对应盒子 $q=r(e^*)$。处理 $e^*$ 时另一个盒子必处于对立活跃集；又因为在维度 $2..d$ 上它们投影相交，AABB Tree 的 `IntersectReport`（保守剪枝不漏真命中）会访问到该叶条目，且叶子处用 $\texttt{Intersect}^{\star}$ 过滤后会被保留，于是 pair 被加入 `Pairs`。k=2,d=d(Enumerat+Sampling)

**不重：** 同一 pair 不会在更早 START 时产生（那时更靠后的盒子尚未活跃），也不会在另一个 START 再次被添加；事件排序的 tie-break 保证“更靠后 START”唯一，从而归属唯一。k=2,d=d(Enumerat+Sampling)

因此 `Pairs` 一一枚举 $J$。

#### 4.1.2 采样均匀性与独立性

给定 `Pairs` 是 $J$ 的一一枚举数组：

- 单次采样：对任意 $P\in J$，
  $$
  \Pr(\texttt{Ans}[j]=P)=\frac{1}{W}=\frac{1}{|J|}.
  $$

- 不同 $j$ 的 $\text{idx}_j$ 独立生成，因此样本序列 i.i.d.（有放回）。

------

### 4.2 复杂度

由于 baseline 显式存储并输出全部 $J$，其复杂度必然含 $|J|$ 项（最坏 $\Theta(n^2)$）。

#### 4.2.1 时间复杂度

1. **事件排序：** $O(n\log n)$
2. **扫描阶段：**

- 每个事件一次 `Insert/Delete`：动态 AABB Tree 工程上常见为近似 $O(\log n)$（与树质量/启发式有关），但最坏可退化到 $O(n)$。

- 每个 START 一次 `IntersectReport`：设该次查询访问的 BVH 节点数为 $V_e$，叶候选数为 $C_e$，真实输出（追加到 `Pairs` 的数量）为 $k_e$，则：
  $$
  T_e = O(V_e + C_e) + O(C_e)\quad(\text{叶过滤})
  $$
  且 $k_e\le C_e$。

汇总所有 START 的真实输出：
$$
\sum_e k_e = |J|.
$$
因此扫描总时间可写成参数化形式：
$$
T_{\text{sweep}}
=
O(n\log n)
+
O\Big(\sum \text{UpdateCost}\Big)
+
O\Big(\sum_{e\in E}(V_e + C_e)\Big)
+
O(|J|).
$$
在剪枝有效、树较平衡的常见情形下，实践中常接近“输出敏感”：
$$
T_{\text{sweep}} \approx O(n\log n + |J|).
$$
最坏情形可能退化到二次级（例如大量重叠导致 $V_e=\Theta(n)$）。但 baseline 本身在 $|J|=\Theta(n^2)$ 时也不可避免是二次量级。

1. **采样阶段：** $O(t)$

综合：
$$
T_{\text{baseline(AABB)}}
=
O(n\log n)
+
O\Big(\sum \text{UpdateCost}\Big)
+
O\Big(\sum_{e}(V_e + C_e)\Big)
+
O(|J|+t).
$$

------

#### 4.2.2 空间复杂度

- 两棵 AABB Tree：$O(n)$ 级节点空间（实现常数因子取决于节点结构）
- `Events`：$O(n)$
- `Pairs`：$\Theta(|J|)$
- 输出 `Ans`：$O(t)$

因此：
$$
S_{\text{baseline(AABB)}} = O(n + |J| + t),
$$
最坏由 $|J|=\Theta(n^2)$ 主导（这是“枚举型 baseline”的必然代价）。

------

### Baseline 总结一句话（Plane Sweep + AABB Tree 版）

> **用 $x_1$ 做 plane sweep 赋予每个相交对唯一归属；用两棵动态 AABB Tree 维护对立活跃集在维度 $2..d$ 的投影，START 时做相交报告并用半开精确判定过滤，得到完整相交对数组 `Pairs`，最后在 `Pairs` 上独立均匀抽下标得到 i.i.d. 样本。**

### 5.3 Adaptive+Sampling 版本分析

### 5.1 正确性分析

#### 5.1.1 不交分解 $J=\biguplus_e J_e$ 仍成立

半开区间 + END-before-START + START tie-break 给出事件全序；任意跨集合相交对 $(r_c,r_{\bar c})$ 在第 1 维必存在“更靠后”的 START 事件 $e^*$，在处理 $e^*$ 时另一盒已活跃，因此该 pair 只会归入唯一 $J_{e^*}$，不重不漏。

#### 5.1.2 `CountIntersect` 得到的 $w_e$ 是精确的

处理 START(q) 时：

- 对立活跃集在对应 AABB-tree 中一一存储为叶子
- `CountIntersect(Q=q^\star)`：
  - 剪枝只排除一定不相交的子树（不会漏）
  - `node ⊆ Q` 的 full-accept 返回 `size(node)` 是正确的（子树所有叶都在 Q 内，必相交）
  - 叶子用半开精确判定
- 因此返回值恰为 $|K_e|$，即 $w_e=|J_e|$。

#### 5.1.3 `SampleIntersect(Q,t_e)` 的均匀性与 i.i.d.

BuildGuide 得到 `cnt(v)`（子树命中数）后，`SampleOne` 采用“按命中数加权选择子树 + 子树内递归”的方式：

- 对于非 full-accept 的内部节点：
  $$
  \Pr(\text{选到左子树})=\frac{cnt(L)}{cnt(v)}
  $$

- 到达某个命中叶时返回该叶

可证明对任意命中叶 $\ell\in K(Q)$：
$$
\Pr(\text{SampleOne 返回 }\ell)=\frac{1}{|K(Q)|}.
$$
重复 $t_e$ 次且每次使用独立随机数，得到有放回 i.i.d. 均匀样本序列。

#### 5.1.4 全局均匀性与独立性（Ours 分支）

对任意 $P\in J$，存在唯一事件 $e^*$ 使 $P\in J_{e^*}$，且 $|J_{e^*}|=w_{e^*}$。

对任意输出位置 $j$：

- Phase2：$\Pr(E_j=e)=w_e/|J|$
- Phase3：给定 $E_j=e$，局部均匀采样 $\Pr(P\mid E_j=e)=1/w_e$

故：
$$
\Pr(Z_j=P)=\frac{w_{e^*}}{|J|}\cdot\frac{1}{w_{e^*}}=\frac{1}{|J|}.
$$
独立性来自两层独立随机源：

- $(E_j)$ 独立抽取
- 每个事件内 `SampleIntersect` 产生 i.i.d. 序列
- slot 回填仅是确定性赋值，不引入相关性

因此 $(Z_1,\dots,Z_t)$ 全局 i.i.d. 均匀。

#### 5.1.5 切换不影响正确性

切换只改变是否继续 `ReportIntersect` / 是否保留 `AllPairs`，但不改变：

- 事件顺序
- Insert/Delete 时刻（活跃集演化）
- `CountIntersect` 得到的精确 $w_e$

因此：

- 未切换：等价于完整 Baseline（枚举 + 数组采样）
- 已切换：等价于从一开始就只计数，然后做 Phase2+Phase3 的 Ours 分支

两条路径都严格正确。 k=2,d=d(Adaptive+Sampling)

------

### 5.2 时间复杂度（参数化 + 最坏退化说明）

AABB-tree 与 R-tree 一样：**最坏情况可能退化**（例如大量包围盒重叠导致剪枝失效），因此这里给出“参数化”复杂度表达，并明确最坏退化风险。

记：

- $h$：AABB-tree 高度（平衡时 $h=O(\log n)$）

- 对一次查询 $Q$，令

  - $V(Q)$：遍历中访问到的节点数
  - $L(Q)$：检查到的叶子数（命中候选叶数）

- 则典型实现中：

  - `CountIntersect(Q)`：$O(V(Q)+L(Q))$

  - `ReportIntersect(Q)`：$O(V(Q)+L(Q)+k)$（$k$ 为真实输出数）

  - `SampleIntersect(Q,t')`：先 BuildGuide $O(V(Q)+L(Q))$，再采样 $t'$ 次，每次随机下降 $O(h)$：
    $$
    O(V(Q)+L(Q)+t'\cdot h).
    $$

- `Insert/Delete`：常见均摊 $O(h)$，但最坏可到 $O(n)$

#### Case A：未切换（$|J|\le J_\star$，Baseline 分支）

Phase1：

- 事件排序：$O(n\log n)$
- 扫描更新：约 $O(n\cdot h)$
- 每个 START：
  - Count：$O(V(Q_e)+L(Q_e))$
  - Report：$O(V(Q_e)+L(Q_e)+w_e)$

汇总（注意 $\sum_e w_e=|J|\le J_\star$）：
$$
T_{\text{base-branch}}
=
O(n\log n + n\cdot h)
+
\sum_{\text{START }e} O(V(Q_e)+L(Q_e))
+
O(|J|)
+
O(t).
$$

#### Case B：发生切换（$|J|>J_\star$，Ours 分支）

Phase1（最多枚举 $J_\star$ 个 pair）：
$$
T_1=
O(n\log n + n\cdot h)
+
\sum_{\text{START }e} O(V(Q_e)+L(Q_e))
+
O(J_\star).
$$
Phase2：$O(n+t)$

Phase3（第二次扫描 + 局部采样）：

- 更新：$O(n\cdot h)$
- 对有 slot 的事件 $e$：`SampleIntersect(Q_e,t_e)=O(V(Q_e)+L(Q_e)+t_e h)`
- 且 $\sum_e t_e=t$

所以：
$$
T_3
=
O(n\cdot h)
+
\sum_{e: t_e>0} O(V(Q_e)+L(Q_e))
+
O(t\cdot h).
$$
合并：
$$
T_{\text{ours-branch}}
=
O(n\log n + n\cdot h)
+
\sum_{\text{START }e} O(V(Q_e)+L(Q_e))
+
\sum_{e: t_e>0} O(V(Q_e)+L(Q_e))
+
O(t\cdot h)
+
O(J_\star).
$$

#### 最坏退化必须写清

在最坏数据（大量节点 AABB 高重叠）下，可能出现 $V(Q)=\Theta(n)$，从而两遍扫描的总成本可退化到 $\Theta(n^2)$。这也是用 AABB-tree（工程索引）替代你原来的 $\log^{d-1}n$ 最坏界结构时不可避免的代价：**分布正确性仍严格成立，但最坏复杂度不再有强保证**。

------

### 5.3 空间复杂度

公共部分：

- 两棵 AABB-tree：$O(n)$
- 事件数组与 $w_e$：$O(n)$
- alias：$O(n)$
- slot 与输出：$O(t)$

Baseline 分支额外：

- `AllPairs`：$\Theta(|J|)\le \Theta(J_\star)$

因此峰值空间：
$$
S_{\text{peak}}=O\bigl(n+\max(t,J_\star)\bigr).
$$

------

### 一句话总结（Plane Sweep + AABB-tree Adaptive）

> **Adaptive-AABBTree**：用 $x_1$ 的 plane sweep 把相交对分解到唯一 START 事件；用两棵动态 AABB-tree 维护对立活跃集在维度 $2..d$ 的投影；Phase1 永远精确 Count 得到事件权重 $w_e$，当累计 $|J|\le J_\star$ 时顺便 Report 枚举并数组采样；一旦超阈值就停止枚举并走两遍扫描：Phase2 事件 alias 分配 slot，Phase3 在每个事件用 AABB-tree 的 “BuildGuide(计数) + 按计数加权随机下降” 实现局部集合上的 **精确 i.i.d. 均匀采样**，最终得到全局 $J$ 上的 i.i.d. 均匀样本序列。