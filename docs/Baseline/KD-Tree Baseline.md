# KD‑Join Sampling Baseline 算法文档

**（Plane Sweep + 2(d−1) 维点嵌入 + 静态 KD‑tree（Active on/off） + 事件级 alias+slot）**

------

## 1. 问题定义与分析

### 1.1 输入：两类 d 维半开轴对齐盒

在 $d\ge 2$ 维欧氏空间 $\mathbb{R}^d$ 中，给定两类轴对齐半开盒子集合：
$$
R_c=\{r_{c1},\dots,r_{c n_1}\},\qquad
R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\},\qquad
n=n_1+n_2.
$$
每个盒子为半开形式：
$$
r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\qquad L_i(r)<R_i(r).
$$
一维半开区间相交当且仅当：
$$
[a,b)\cap[c,d)\neq\varnothing\iff \max(a,c)<\min(b,d).
$$

------

### 1.2 Join 结果集合与采样目标

只关心跨集合相交对：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
目标输出 $t$ 个样本 $Z_1,\dots,Z_t\in J$，要求 **i.i.d. 均匀有放回**：

- **均匀性**：$\Pr(Z_j=P)=\frac{1}{|J|}$，对任意 $P\in J$；
- **独立性**：$Z_1,\dots,Z_t$ 相互独立。

------

### 1.3 第 1 维 Plane Sweep 与事件系统

选择第 1 维 $x_1$ 扫描。每个盒子 $r$ 产生两事件：

- $\texttt{START}(r)$ at $x_1=L_1(r)$
- $\texttt{END}(r)$ at $x_1=R_1(r)$

**事件全序（稳定排序，确保“唯一归属”）：**

1. 按 $x_1$ 升序；
2. 同 $x_1$：**END 在 START 之前**（与半开边界匹配，避免贴边误相交）；
3. 同 $x_1$ 且同为 START：固定集合优先级（例如先 $R_c$ 后 $R_{\bar c}$），再按全局唯一 id 打破平局。

扫描到位置 $x_0$ 时定义活跃集：
$$
\mathrm{Active}(R_c)=\{r_c\in R_c: x_0\in[L_1(r_c),R_1(r_c))\},\quad
\mathrm{Active}(R_{\bar c})\ \text{同理}.
$$
**扫描时序约定（影响正确性）：**

- 处理 $\texttt{START}(q)$：**先查询对立活跃集，再激活 $q$**；
- 处理 $\texttt{END}(r)$：先失活 $r$。

------

### 1.4 事件诱导的块分解：$J=\biguplus_{e\in E} J_e$

令 $E$ 为所有 START 事件集合。对任一 START 事件 $e$，令新盒子 $q=r(e)$，对立活跃集记为 $\mathrm{Active}(\text{Other})$。定义：
$$
K_e := \{ r\in \mathrm{Active}(\text{Other}) : q\cap r\neq\varnothing\},\qquad w_e:=|K_e|.
$$
统一输出 pair 顺序为 $(r_c,r_{\bar c})$：

- 若 $q\in R_c$：$\displaystyle J_e=\{(q,r_{\bar c})\mid r_{\bar c}\in K_e\}$
- 若 $q\in R_{\bar c}$：$\displaystyle J_e=\{(r_c,q)\mid r_c\in K_e\}$

在上述事件全序 + 半开语义 + “START 先查后插”下成立：
$$
J=\biguplus_{e\in E} J_e,\qquad |J|=\sum_{e\in E} w_e.
$$
直观理解：任意相交对 $(r_c,r_{\bar c})$ 在事件全序中存在唯一“更靠后”的 START 事件 $e^\star$。在处理 $e^\star$ 时另一盒子已处于对立活跃集，因此该 pair **只会**在 $J_{e^\star}$ 中出现一次。

------

### 1.5 维度 2..d 的过滤：点嵌入与正交范围查询

在处理 $\texttt{START}(q)$ 时，对立活跃盒子与 $q$ 在第 1 维必然真重叠（半开 + END-before-START 保证不是仅贴边）。因此整盒相交等价于维度 $2..d$ 投影相交。

令 $m=d-1$，定义
$$
k:=2m=2(d-1).
$$
对盒子 $r$，取端点向量（只取维度 $2..d$）：
$$
\mathbf{L}(r)=(L_2(r),\dots,L_d(r))\in\mathbb{R}^m,\quad
\mathbf{R}(r)=(R_2(r),\dots,R_d(r))\in\mathbb{R}^m.
$$
定义 $k$ 维点嵌入：
$$
p(r)=(L_2(r),\dots,L_d(r),\ R_2(r),\dots,R_d(r))\in\mathbb{R}^{2m}.
$$
对查询盒子 $q$，任意 $i\ge 2$ 维上有：
$$
[L_i(r),R_i(r))\cap[L_i(q),R_i(q))\neq\varnothing
\iff
L_i(r)<R_i(q)\ \land\ R_i(r)>L_i(q).
$$
因此定义正交范围（严格不等号，带无穷边界）：
$$
\mathcal{Q}(q)=
\prod_{i=2}^{d}(-\infty,R_i(q))
\times
\prod_{i=2}^{d}(L_i(q),+\infty).
$$
**等价性：**
$$
r \text{ 与 } q \text{ 在维度 }2..d\text{ 相交}\iff p(r)\in\mathcal{Q}(q).
$$

------

## 2. 核心数据结构

KD‑tree baseline 的核心结构由四部分组成：

1. plane sweep 事件数组 `Events`
2. 严格不等号的 rank 化 + 点嵌入 $p(r)$
3. 两棵静态 KD‑tree（分别维护两侧点集） + 动态 active/cnt
4. alias + slot（事件级抽样规划）

------

### 2.1 事件数组

`Events[1..2n]` 每个元素包含：
$$
(x_1,\ \mathrm{type}\in\{\texttt{START},\texttt{END}\},\ \mathrm{class}\in\{c,\bar c\},\ \mathrm{id}).
$$
按 §1.3 的全序稳定排序。

------

### 2.2 严格不等号的可实现化：二级键 rank 化

为了将 $\mathcal{Q}(q)$ 中的严格不等号稳定映射到整数闭区间，使用二级键 $(\text{value},\text{id})$ 做坐标压缩。

为所有盒子分配全局唯一 id：
$$
\mathrm{id}(r)\in\{1,2,\dots,n\}.
$$
并设哨兵：
$$
\mathrm{LOW}=0,\qquad \mathrm{HIGH}=n+1.
$$
对每个维度 $i\in\{2,\dots,d\}$，建立两份排序数组（全局口径，包含两侧所有盒子）：

- $A_i$：排序 $\{(L_i(r),\mathrm{id}(r))\}$
- $B_i$：排序 $\{(R_i(r),\mathrm{id}(r))\}$

定义 rank 坐标（0-based）：
$$
\mathrm{rank}_{A_i}(r)=\text{pair }(L_i(r),\mathrm{id}(r))\text{ 在 }A_i\text{ 中的位置},
$$
于是点嵌入变为整数点：
$$
p(r)=(
\mathrm{rank}_{A_2}(r),\dots,\mathrm{rank}_{A_d}(r),\
\mathrm{rank}_{B_2}(r),\dots,\mathrm{rank}_{B_d}(r)
)\in\mathbb{Z}^k.
$$
对查询盒子 $q$，将 $\mathcal{Q}(q)$ 转为 rank 空间中的闭区间笛卡尔积 $Q(q)$：

- 约束 $L_i(r) 在 $A_i$ 上对应：
  $$
  U_i=\mathrm{lower\_bound}_{A_i}\bigl((R_i(q),\mathrm{LOW})\bigr)-1,\quad
  \mathrm{rank}_{A_i}\in[0,U_i].
  $$

- 约束 $R_i(r)>L_i(q)$ 在 $B_i$ 上对应：
  $$
  L'_i=\mathrm{upper\_bound}_{B_i}\bigl((L_i(q),\mathrm{HIGH})\bigr),\quad
  \mathrm{rank}_{B_i}\in[L'_i,|B_i|-1].
  $$

若存在某维 $U_i<0$ 或 $L'_i>|B_i|-1$，则 $Q(q)$ 为空，查询直接返回空。

------

### 2.3 两棵 KD‑tree：静态骨架 + 动态 active 标记

构建两棵静态平衡 KD‑tree（median split 或等价平衡策略）：

- $\mathcal{KD}^{(c)}$：包含点集 $\{p(r_c)\mid r_c\in R_c\}$
- $\mathcal{KD}^{(\bar c)}$：包含点集 $\{p(r_{\bar c})\mid r_{\bar c}\in R_{\bar c}\}$

扫描阶段不进行结构性插删，仅做动态标记：

- `active(v) ∈ {0,1}`：该结点点是否活跃
- `cnt(v)`：该结点子树内活跃点数（包含结点自身点）

每棵树维护 `node_ptr[id] -> node`，用于 $O(\log n)$ 找到对应结点并更新。

#### 2.3.1 结点字段（每棵树内部）

对每个结点 $v$：

- `pt(v) = p(r)`：该结点存储的整数点
- `id(v)`：盒子 id
- `left(v), right(v), parent(v)`
- `bbox(v)`：子树点集的包围盒（每维 min/max，静态预计算，闭区间）
- `active(v)`：动态
- `cnt(v)`：动态

（可选增强）

- `Bucket(v)`：该子树当前活跃点 id 的动态数组 + 句柄（用于子树内 $O(1)$ 均匀抽样与快速报告）。启用后总空间约 $O(n\log n)$。

#### 2.3.2 Activate / Deactivate

`Activate(id)`：

1. $v\leftarrow node\_ptr[id]$
2. 若 `active(v)=0`，置为 1，并沿 parent 向上对 `cnt(u)` 加 1（若启用 Bucket，也在沿途 Bucket 中插入并维护句柄）

`Deactivate(id)` 对称（置 0，沿途减 1）。

在平衡 KD‑tree 下树高 $h\approx O(\log n)$，因此更新通常为 $O(h)$。

------

### 2.3.3 KD‑tree 查询接口：Count / Report / Sample

对 rank 空间范围 $Q=\prod_{j=1}^{k}[L_j,U_j]$（闭区间），命中集合定义为：
$$
\mathrm{Hit}(Q)=\{r\ \text{active}: p(r)\in Q\}.
$$
KD‑tree 支持：

- `Count(Q)`：返回 $|\mathrm{Hit}(Q)|$
- `Report(Q)`：枚举并返回 $\mathrm{Hit}(Q)$ 中所有 id
- `Sample(Q,t')`：从 $\mathrm{Hit}(Q)$ 中 **i.i.d. 均匀有放回**采样 $t'$ 次

#### A. Count(Q)（精确计数）

递归 `Count(v,Q)`：

- 若 $v=\varnothing$ 或 `cnt(v)=0`：返回 0

- 若 `bbox(v)` 与 $Q$ 不相交：返回 0

- 若 `bbox(v) ⊆ Q`：返回 `cnt(v)`（整子树均命中）

- 否则：
  $$
  ans = \mathbf{1}[active(v)=1\land pt(v)\in Q] + Count(left,Q)+Count(right,Q)
  $$

#### B. Report(Q)（精确枚举）

递归 `Report(v,Q)`：

- 若 $v=\varnothing$ 或 `cnt(v)=0` 或 `bbox(v)` 与 $Q$ 不相交：返回空
- 若 `bbox(v) ⊆ Q`：
  - 若启用 Bucket：直接输出 `Bucket(v)`（全部活跃点）
  - 否则：继续下行遍历并输出子树内所有活跃点（代价与 `cnt(v)` 同阶）
- 否则：
  - 若 `active(v)=1` 且 `pt(v)∈Q`：输出 `id(v)`
  - 递归左右子树并合并输出

#### C. Sample(Q,t′)：canonical pieces + piece‑slot 批量采样（保证“有序 i.i.d.”）

**核心思路：将 $\mathrm{Hit}(Q)$ 拆成不交的 canonical pieces，并按 piece 大小加权抽取。**

**Step 1：构造 pieces（不交并分解）**
 递归 `CollectPieces(v,Q)` 生成列表 $\mathcal{P}$，其中每个 piece 要么是“全包含子树结点”要么是“单点结点”。

- 若 $v=\varnothing$ 或 `cnt(v)=0`：跳过
- 若 `bbox(v)` 与 $Q$ 不相交：跳过
- 若 `bbox(v) ⊆ Q`：加入一个 **SubtreePiece(v)**，权重 $w=\mathrm{cnt}(v)$，并停止向下展开
- 否则（部分覆盖）：
  - 若 `active(v)=1` 且 `pt(v)∈Q`：加入一个 **PointPiece(v)**，权重 $w=1$
  - 继续递归左右子树

这样得到：
$$
\mathrm{Hit}(Q)=\biguplus_{P\in\mathcal{P}} S_P,\quad w(P)=|S_P|,\quad W=\sum_{P\in\mathcal{P}} w(P)=|\mathrm{Hit}(Q)|.
$$
（若 $W=0$，直接返回空。）

**Step 2：在 pieces 上按权重抽取，并用 slot 记录每次抽取的位置**
 构建 alias（或前缀和）使得 $\Pr(P)=w(P)/W$。
 对样本位置 $u=1..t′$ 独立抽 piece $P_u$，并把位置 $u$ 放入 `Slot[P_u]`。

> 这里的 **Slot 是保证“返回序列有序 i.i.d.”的关键**：后续按 piece 批量生成样本后，再回填到对应位置，输出序列分布等价于逐次独立采样。

**Step 3：piece 内部均匀采样并回填**

- 若 piece 是 `PointPiece(v)`：对 `Slot` 中每个位置直接回填 `id(v)`（允许重复）
- 若 piece 是 `SubtreePiece(v)`：需要从结点 $v$ 子树的活跃点中均匀抽样。提供两种实现：

**实现 C‑1（无需 Bucket，按 cnt 权重随机下降，单次 $O(h)$）**
 定义 `SampleSubtreeUniform(v)`：

- 令
   $w_{self}=\mathbf{1}[active(v)=1]$
   $w_L=\mathrm{cnt}(left(v))$（若无左子则 0）
   $w_R=\mathrm{cnt}(right(v))$（若无右子则 0）
   总和 $w_{self}+w_L+w_R=\mathrm{cnt}(v)$
- 以概率 $w_{self}/cnt(v)$ 直接返回 `id(v)`；否则按 $w_L,w_R$ 比例递归到相应孩子。

对 `Slot[SubtreePiece(v)]` 中每个位置调用一次 `SampleSubtreeUniform(v)` 并回填。
 总成本约 $O(|Slot|\cdot h)$。

**实现 C‑2（启用 Bucket，子树内单次 $O(1)$）**
 `Bucket(v)` 存子树活跃点列表，则均匀采样为随机下标访问即可。
 总成本约 $O(|Slot|)$（但空间约 $O(n\log n)$）。

------

### 2.4 事件级 alias + slot

Sampling / Adaptive 的两遍扫描分支需要：

- 事件权重：每个 START 事件 $e$ 的 $w_e=|K_e|=|J_e|$
- 事件级 alias：在 $E$ 上按 $\Pr(E=e)=w_e/|J|$ 构建
- 事件 slot：对每个事件 $e$，维护 $S_e\subseteq\{1,\dots,t\}$ 表示哪些输出位置分配给该事件

------

## 3. 算法详细流程

以下三版本共享相同事件系统、rank 化、KD‑tree 构造与 Active 维护方式。

### 3.0 共同预处理

输入：$R_c,R_{\bar c}$，样本数 $t$（Adaptive 另有阈值 $J_\star$）。

1. 为每个盒子分配全局唯一 id（1..n），并保存 id→盒子对象映射；
2. 构造 `Events`：每盒子产生 START/END 各一个；
3. 按 §1.3 全序稳定排序；
4. 对每维 $i=2..d$ 构建排序数组 $A_i,B_i$，并为每个盒子计算 rank 点 $p(r)\in\mathbb{Z}^k$；
5. 离线构建两棵平衡 KD‑tree：$\mathcal{KD}^{(c)}$、$\mathcal{KD}^{(\bar c)}$；初始化所有 `active=0`，`cnt(root)=0`；
6. 实现 `BuildRange(q)`：根据 §2.2 的二分规则将 $q$ 映射为 rank 范围 $Q(q)$，若为空则返回空标记。

------

### 3.1 版本一：Sampling（两遍扫描 + 事件 alias + 局部 KD.Sample）

适用：不希望时间/空间线性依赖 $|J|$（但接受 KD‑tree 的平均性能与可能退化）。

#### Phase 1：第一次扫描（只计数，得到 $w_e$ 与 $W=|J|$）

初始化：两棵 KD‑tree 全 inactive；`W=0`；为每个 START 事件存储 `w_e`。

按 `Events` 扫描：

- `END(r)`：对所属类的 KD 进行 `Deactivate(id(r))`
- `START(q)`（事件 $e$）：
  1. $Q \leftarrow BuildRange(q)$，若空则 $w_e=0$
  2. 若 $q\in R_c$：$w_e \leftarrow \mathcal{KD}^{(\bar c)}.\texttt{Count}(Q)$
      若 $q\in R_{\bar c}$：$w_e \leftarrow \mathcal{KD}^{(c)}.\texttt{Count}(Q)$
  3. `W += w_e`
  4. `Activate(id(q))` 到本侧 KD

扫描结束：由 §1.4 得 $W=\sum_e w_e=|J|$。若 $W=0$ 返回空。

#### Phase 2：事件级 alias + slot 分配

1. 在 START 事件集合 $E$ 上按权重 $w_e/W$ 构建 alias；
2. 初始化每个事件 slot：$S_e\gets\emptyset$；
3. 对每个输出位置 $j=1..t$：
   - 抽事件 $E_j\sim (w_e/W)$
   - 加入 $S_{E_j}\leftarrow S_{E_j}\cup\{j\}$

记 $t_e=|S_e|$，则 $\sum_e t_e=t$。

#### Phase 3：第二次扫描（按事件做局部批量采样 + 回填）

清空两棵 KD‑tree 活跃状态（结构保留），初始化输出数组 `Ans[1..t]`。

按 `Events` 再扫描一遍：

- `END(r)`：对所属类 `Deactivate(id(r))`
- `START(q)`（事件 $e$）：
  - 若 $t_e=0$：仅 `Activate(id(q))`
  - 若 $t_e>0$：
    1. $Q\leftarrow BuildRange(q)$，若空则返回空列表（理论上此时 $w_e=0$，slot 也应为空）
    2. 若 $q\in R_c$：`List ← KD^(barc).Sample(Q, t_e)`
        若 $q\in R_{\bar c}$：`List ← KD^(c).Sample(Q, t_e)`
    3. 令 $S_e$ 按固定顺序排列为 $(j_1,\dots,j_{t_e})$，回填：
       - 若 $q\in R_c$：`Ans[j_u]=(q, box(List[u]))`
       - 若 $q\in R_{\bar c}$：`Ans[j_u]=(box(List[u]), q)`
  - 最后 `Activate(id(q))`

输出 `Ans`。

------

### 3.2 版本二：Enumerate+Sampling（显式枚举 $J$ + 数组采样）

适用：$|J|$ 较小，可承受 $\Theta(|J|)$ 存储；实现最简单、常数通常更小。

#### Step 1：一次扫描枚举全部 $J$

初始化：`Pairs=[]`；两棵 KD‑tree 全 inactive。

按 `Events` 扫描：

- `END(r)`：对所属类 `Deactivate(id(r))`
- `START(q)`：
  1. $Q\leftarrow BuildRange(q)$，若空则 `List=[]`
  2. 若 $q\in R_c$：`List ← KD^(barc).Report(Q)`
      若 $q\in R_{\bar c}$：`List ← KD^(c).Report(Q)`
  3. 对 `List` 中每个对立盒子 id 形成 pair：
     - 若 $q\in R_c$：`Pairs.append((q, r_bar))`
     - 若 $q\in R_{\bar c}$：`Pairs.append((r_c, q))`
  4. `Activate(id(q))`

扫描结束：`|Pairs|=|J|`。

#### Step 2：数组 i.i.d. 下标采样

- 若 `Pairs` 为空：返回空；
- 否则对 $j=1..t$：
  - `idx_j ~ Unif{0..|Pairs|-1}`
  - `Ans[j]=Pairs[idx_j]`

输出 `Ans`。

------

### 3.3 版本三：Adaptive+Sampling（阈值自适应切换）

适用：希望小 $|J|$ 时只扫一遍就采样；大 $|J|$ 时避免 $|J|$ 级枚举存储。

输入额外阈值 $J_\star$。

#### Phase 1：一次扫描（永远 Count；阈值内可选枚举）

初始化：

- `mode = ENUMERATE`
- `AllPairs = []`（仅在 ENUMERATE 时增长）
- `W=0`，并存储每个 START 的 `w_e`

两棵 KD‑tree 全 inactive。按 `Events` 扫描：

- `END(r)`：`Deactivate(id(r))`
- `START(q)`（事件 $e$）：
  1. **先计数**
     - $Q\leftarrow BuildRange(q)$，若空则 $w_e=0$
     - `w_e ← KD_other.Count(Q)`，`W += w_e`
     - 保存 `w_e`
  2. **阈值判断 + 可选枚举**
     - 若 `mode==ENUMERATE` 且 `W ≤ J_star`：
       - `List ← KD_other.Report(Q)`
       - 将 `List` 形成 pairs 追加到 `AllPairs`（保持输出顺序统一）
     - 若 `mode==ENUMERATE` 且 `W > J_star`：
       - `mode ← COUNT_ONLY`
       - 释放/丢弃 `AllPairs`
       - 后续不再 Report
  3. `Activate(id(q))`

扫描结束：得到精确 $W=|J|$。若 $W=0$ 返回空。

#### 分支 A：未切换（$W\le J_\star$）

此时 `AllPairs` 完整枚举 $J$，直接数组采样（同 §3.2 Step2）。

#### 分支 B：已切换（$W>J_\star$）

执行 Sampling 的 Phase2 + Phase3（用 Phase1 已保存的 $w_e,W$ 构建事件 alias 与 slot，再二遍扫描调用 `KD.Sample` 生成样本）。

------

## 4. Adaptive 阈值 $J_\star$ 的选择策略

### 4.1 内存硬约束（必须满足）

设可用内存预算 `MemBudget` 字节，允许 `AllPairs` 使用比例 $\rho\in(0,1)$。每个 pair 的存储开销 `sizeof(Pair)`（通常为两个 32/64-bit id）。则：
$$
J_\star^{\mathrm{mem}}=\left\lfloor\frac{\rho\cdot \mathrm{MemBudget}}{\mathrm{sizeof(Pair)}}\right\rfloor,\qquad
J_\star\le J_\star^{\mathrm{mem}}.
$$

------

### 4.2 时间权衡（可落地的经验模型）

KD‑tree 查询代价高度依赖分布与维数 $k=2(d-1)$。建议用“扫描成本 + 样本数”来定阈值：

- 更新代价：$U(n)\approx O(h)$（平衡树 $h\approx O(\log n)$）
- 范围计数/遍历代价：用 $Q_{\mathrm{KD}}(n,k)$ 表示（平均情形常写作 $O(n^{1-1/k})$ 量级；最坏可到 $O(n)$）

一次“只计数扫描”的近似代价：
$$
T_{\mathrm{sweep}}\approx 2n\cdot U(n)+n\cdot Q_{\mathrm{KD}}(n,k).
$$
发生切换时，大分支需要两遍扫描，并且 Phase1 额外尝试枚举的输出量 $\le J_\star$。因此可取：
$$
J_\star^{\mathrm{time}} = C_1\cdot T_{\mathrm{sweep}} + C_2\cdot t,\qquad
J_\star=\min(J_\star^{\mathrm{mem}}, J_\star^{\mathrm{time}}),
$$
其中 $C_1,C_2$ 由 benchmark 标定。

------

### 4.3 交叉点拟合（最推荐）

1. 固定代表性分布、维数 $d$、样本数 $t$，变化 $n$；
2. 分别运行：
   - Enumerate+Sampling（§3.2）
   - Sampling（§3.1）
3. 找到耗时交叉点 $|J_{\mathrm{cross}}|$；
4. 取 $J_\star \approx 0.8\cdot |J_{\mathrm{cross}}|$，并用 $J_\star^{\mathrm{mem}}$ 截断。

------

## 5. 算法分析（正确性与复杂度）

### 5.1 正确性

#### 引理 1（事件分块不交分解）

在 §1.3 的事件全序 + 半开语义 + “START 先查后插”下：
$$
J=\biguplus_{e\in E} J_e,\qquad |J|=\sum_{e\in E} w_e.
$$
证明要点：任意 pair 的两次 START 中存在唯一“更靠后”的 START，且在该时刻两者第一次同时活跃，因此只归属该事件一次。

#### 引理 2（点嵌入与范围查询等价）

对任意 $r,q$：
$$
r \text{ 与 } q \text{ 在维度 }2..d\text{ 相交}\iff p(r)\in\mathcal{Q}(q).
$$
rank 化构造 $Q(q)$ 保证严格不等号被精确实现，因此：
$$
p(r)\in\mathcal{Q}(q)\iff p(r)\in Q(q)\ (\text{rank 空间闭区间}).
$$

#### 引理 3（KD.Count / KD.Report 的语义正确）

KD‑tree 在只统计 active 点且使用 `bbox` 剪枝 + `bbox⊆Q` 整段接受时，`Count(Q)` 精确返回 $|\mathrm{Hit}(Q)|$，`Report(Q)` 精确枚举 $\mathrm{Hit}(Q)$。

#### 引理 4（KD.Sample(Q,t′) 输出 i.i.d. 均匀）

`CollectPieces` 给出不交并分解 $\mathrm{Hit}(Q)=\biguplus_{P} S_P$，并记录 $w(P)=|S_P|$。
 对任一命中点 $x\in S_{P^\star}$：

- 抽中 piece $P^\star$ 的概率：$w(P^\star)/W$
- 给定 $P^\star$，在 $S_{P^\star}$ 内均匀抽到 $x$ 的概率：$1/w(P^\star)$

因此单次输出命中点的概率为：
$$
\Pr(\text{output}=x)=\frac{w(P^\star)}{W}\cdot\frac1{w(P^\star)}=\frac1W=\frac1{|\mathrm{Hit}(Q)|}.
$$
piece‑slot 机制保证每个位置独立抽 piece，再独立在 piece 内采样并回填，因此输出序列为 **有序 i.i.d.**。

------

### 5.2 三个版本的正确性结论

#### 定理 A（Enumerate+Sampling）

- 枚举阶段：由引理 1，所有相交对恰好在其“更靠后 START”对应的事件被报告一次，因此 `Pairs` 一一枚举 $J$。
- 采样阶段：对 `Pairs` 做独立均匀下标采样，得到 $J$ 上 i.i.d. 均匀有放回样本。

#### 定理 B（Sampling）

对任意 $P\in J$，存在唯一事件 $e^\star$ 使 $P\in J_{e^\star}$ 且 $|J_{e^\star}|=w_{e^\star}$。
 对任一输出位置 $j$：

- $\Pr(E_j=e)=w_e/W$
- 给定 $E_j=e$，局部 `KD.Sample` 在 $J_e$ 上均匀，故 $\Pr(P\mid E_j=e^\star)=1/w_{e^\star}$

于是：
$$
\Pr(Z_j=P)=\frac{w_{e^\star}}{W}\cdot\frac1{w_{e^\star}}=\frac1W=\frac1{|J|}.
$$
独立性来自事件抽样独立与局部采样独立（slot 回填仅为位置映射）。

#### 定理 C（Adaptive）

- 未切换：等价 Enumerate+Sampling
- 已切换：Phase1 保留精确 $w_e,W$，后续等价 Sampling

------

### 5.3 复杂度（参数化表达 + 退化说明）

记：

- KD‑tree 高度 $h\approx O(\log n)$（median split 平衡）
- $Q_{\mathrm{KD}}(n,k)$：一次范围遍历访问的结点/检查成本（平均情形常写作 $O(n^{1-1/k})$ 量级；最坏 $O(n)$）
- `Bucket` 方案下，子树内采样单次 $O(1)$；否则单次 $O(h)$

#### 5.3.1 预处理

- 构事件并排序：$O(n\log n)$
- 构 $A_i,B_i$ 并计算 rank：$O((d-1)\,n\log n)$（$d$ 常数则为 $O(n\log n)$）
- 建两棵 KD‑tree：典型 $O(n\log n)$

#### 5.3.2 Enumerate+Sampling

时间：

- 扫描更新：$2n\cdot O(h)=O(nh)$
- 每个 START 做一次 `Report(Q)`：$O(Q_{\mathrm{KD}}(n,k)+k_{\text{out}}(e))$
- $\sum_e k_{\text{out}}(e)=|J|$
- 数组采样：$O(t)$

因此：
$$
T_{\mathrm{enum}} = O\bigl(n\log n + nh + n\cdot Q_{\mathrm{KD}}(n,k) + |J| + t\bigr).
$$
空间：

- KD‑tree：$O(n)$（无 Bucket）或 $O(n\log n)$（有 Bucket）
- `Pairs`：$\Theta(|J|)$
- 输出：$O(t)$

$$
S_{\mathrm{enum}}=O(\mathrm{KD\text{-}space}+|J|+t).
$$

#### 5.3.3 Sampling

Phase1（Count）：
$$
T_1 = O(nh) + n\cdot O(Q_{\mathrm{KD}}(n,k)).
$$
Phase2（事件 alias + slot）：
$$
T_2 = O(n+t),\quad S_2=O(n+t).
$$
Phase3（第二遍扫描 + 采样）：

- 更新：$O(nh)$
- 设 $M=\#\{e:t_e>0\}\le \min(n,t)$
- 每个这样的事件调用一次 `KD.Sample(Q,t_e)`：
  - 构 pieces：$O(Q_{\mathrm{KD}}(n,k))$
  - 生成样本：
    - 无 Bucket：$O(t_e\cdot h)$
    - 有 Bucket：$O(t_e)$
- 汇总 $\sum_e t_e=t$

因此：
$$
T_3 = O(nh) + M\cdot O(Q_{\mathrm{KD}}(n,k)) + O(t\cdot H),
$$
其中 $H=h$（无 Bucket）或 $H=1$（有 Bucket）。

合并：
$$
T_{\mathrm{samp}}=
O\Bigl(n\log n + nh + (n+M)\cdot Q_{\mathrm{KD}}(n,k) + t\cdot H + (n+t)\Bigr).
$$
空间：
$$
S_{\mathrm{samp}}=O(\mathrm{KD\text{-}space}+n+t).
$$

#### 5.3.4 Adaptive+Sampling

- 未切换（$|J|\le J_\star$）：与 Enumerate+Sampling 同阶，且 $|J|\le J_\star$
- 已切换（$|J|>J_\star$）：Phase1 的枚举输出量 $\le J_\star$，后续同 Sampling

$$
T_{\mathrm{adaptive,\,switch}} \approx T_{\mathrm{samp}} + O(J_\star)\ (\text{输出写入量级}),
$$

空间峰值：
$$
S_{\mathrm{adaptive,\,switch}}=O(\mathrm{KD\text{-}space}+\max(J_\star,t)+n).
$$

> **最坏情况说明**：KD‑tree 的范围遍历在高维或严重重叠分布下可能退化到 $Q_{\mathrm{KD}}(n,k)=\Theta(n)$，从而整体时间可退化到 $\Theta(n^2)$ 量级。这是 KD‑tree baseline 的典型特性。