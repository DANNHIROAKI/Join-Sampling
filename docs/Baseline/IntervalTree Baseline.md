# IntervalTree Baseline++：二维 Spatial Join 的 i.i.d. 均匀采样（Sampling / Enumerate+Sampling / Adaptive+Sampling）

## 记号与约定（强烈建议放在文档开头）

- 输入为二维平面 $\mathbb{R}^2$ 上两类轴对齐 **半开**矩形集合：
  $$
  R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
  R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\},\quad
  n=n_1+n_2.
  $$

- 每个矩形采用半开语义：
  $$
  r=[L_x(r),R_x(r))\times[L_y(r),R_y(r)),\quad
  L_x(r)<R_x(r),\ L_y(r)<R_y(r).
  $$

- 一维半开区间相交判定：
  $$
  [a,b)\cap[c,d)\neq\varnothing\iff \max(a,c)<\min(b,d).
  $$

- 统一输出顺序：所有输出 pair 都写为 $(r_c,r_{\bar c})$。当 START 来自 $R_{\bar c}$ 时需要交换顺序回填输出。

- 采样均为 **有放回**，且要求 **i.i.d.**（相互独立）。

------

# 1. 问题定义与分析

## 1.1 Join 结果集合 $J$

我们只关心跨集合相交对（固定顺序 $(r_c,r_{\bar c})$）：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
令 $W:=|J|$。

## 1.2 采样目标（核心任务）

输出 $t$ 个样本 $Z_1,\dots,Z_t\in J$，要求：

- **均匀性**：$\Pr(Z_j=P)=1/W$，对任意 $P\in J$
- **独立性**：$Z_1,\dots,Z_t$ 相互独立（i.i.d.，有放回）

目标是在尽量 **不显式枚举 $J$**（最坏 $\Theta(n^2)$）的前提下，实现上述分布。

------

## 1.3 Plane sweep（按 $x$）与事件分块：$J=\biguplus_e J_e$

### 1.3.1 事件定义

对每个矩形 $r$ 建立两个事件：

- $\texttt{START}(r)$ at $x=L_x(r)$
- $\texttt{END}(r)$ at $x=R_x(r)$

事件记录可定义为：
$$
(x,\ \text{type}\in\{\texttt{END},\texttt{START}\},\ \text{cls}\in\{c,\bar c\},\ \text{id})
$$

### 1.3.2 **事件“稳定全序”**（必须写成“全序”，避免 reviewer 质疑）

定义一个全序 $\prec$ 用于排序：

1. $x$ 小者先；
2. 同 $x$：$\texttt{END}\prec \texttt{START}$（匹配半开边界，防止“贴边误相交”）；
3. 同 $x$ 且同为 START：按固定 class 优先级（例如 $c\prec \bar c$），再按唯一 id 打破平局；
4. 同 $x$ 且同为 END：任意顺序（或也按 id 稳定化，工程上更简单）。

> 这样保证：对任意两个事件，总能比较先后；并且“同点 START”的“更靠后”是严格定义的，不会因为实现细节出现随机性。

### 1.3.3 活跃集与事件块

扫描到某事件的坐标 $x_0$ 时，定义活跃集：
$$
\mathrm{Active}(R_\star)=\{r\in R_\star\mid x_0\in[L_x(r),R_x(r))\}.
$$
对每个 START 事件 $e$，其新加入矩形为 $q=r(e)$。令 “Other” 表示对立类集合（若 $q\in R_c$，则 Other=$R_{\bar c}$；反之同理）。定义：

- 对立活跃相交对象集合：
  $$
  K_e=\{r\in\mathrm{Active}(\text{Other})\mid q\cap r\neq\varnothing\},\qquad
  w_e:=|K_e|.
  $$

- 局部相交对集合（规范化输出顺序）：

  - 若 $q\in R_c$：
    $$
    J_e=\{(q,r_{\bar c})\mid r_{\bar c}\in K_e\}
    $$

  - 若 $q\in R_{\bar c}$：
    $$
    J_e=\{(r_c,q)\mid r_c\in K_e\}
    $$

### 1.3.4 引理 1（事件分块不重不漏）

在半开语义 + END-before-START + START tie-break（全序）下：
$$
J=\biguplus_{e\in E}J_e,\qquad
|J|=\sum_{e\in E} w_e,
$$
其中 $E$ 为所有 START 事件集合。

**证明要点（建议在正文中写 4~6 行）：**

- 任意相交对 $(r_c,r_{\bar c})$ 只可能在这二者 START 之后才被发现；
- 设二者 START 事件中“更靠后”的为 $e^\star$（按全序 $\prec$ 定义），则在处理 $e^\star$ 时，另一个矩形已活跃（若不是活跃，则与 $x$ 维相交矛盾；半开+END-before-START 排除了贴边误判），因此该 pair 必被计入 $J_{e^\star}$；
- 在更早的 START 时刻，后者尚未活跃，因此不可能被计入；
- 因此每个 pair 恰被计入一次，得到不交分解与计数求和。

------

## 1.4 START 时刻：$x$ 维已保证重叠，只需处理 $y$ 区间相交

处理 $\texttt{START}(q)$ 时，$x_0=L_x(q)$。若 $r\in\mathrm{Active}(\text{Other})$，则
$$
x_0\in[L_x(r),R_x(r))\quad\text{且}\quad x_0\in[L_x(q),R_x(q))
$$
因此 **$x$ 维必重叠**。二维相交等价为 $y$ 维区间相交：
$$
q\cap r\neq\varnothing\iff
[L_y(q),R_y(q))\cap[L_y(r),R_y(r))\neq\varnothing.
$$
因此我们将几何子问题归约为：

> 动态维护对立类活跃 $y$-区间集合，对查询区间 $Q=[a,b)$ 支持
>  **Count / Report / Sample（i.i.d. 均匀）**。

------

## 1.5 A/B 分解：把“区间相交”拆成两类互斥子查询

令查询区间 $Q=[a,b)$，活跃区间 $I=[\ell,r)$。则：
$$
I\cap Q\neq\varnothing
\iff (\ell\le a<r)\ \lor\ (a<\ell<b).
$$
定义：

- **A 型（stabbing / 含点）**：$\ell\le a
- **B 型（start-in-range / 起点落入开区间）**：$a<\ell<b$

二者互斥（$\ell=a$ 只属于 A），且覆f盖所有相交情形。令
$$
A(Q)=\{I:\ell\le a<r\},\qquad
B(Q)=\{I:a<\ell<b\},
$$
则
$$
\mathrm{Overlap}(Q)=A(Q)\uplus B(Q),\qquad
|\mathrm{Overlap}(Q)|=|A(Q)|+|B(Q)|.
$$

### 引理 2（A/B 分解互斥且完备）

上述等价式成立，且 $A(Q)$ 与 $B(Q)$ 不交并覆盖 $\mathrm{Overlap}(Q)$。

------

# 2. 核心数据结构

整体结构由两部分组成：

1. **事件系统**：排序后两遍扫描（或一遍枚举）。
2. **IntervalTreeSampler**：维护活跃 $y$-区间，并支持 A/B 两类查询的 Count/Report/Sample。

------

## 2.1 事件系统（Plane Sweep）

### 2.1.1 事件数组

`Events[1..2n]` 每条记录包含：

- `x`：事件坐标
- `type`：END / START
- `cls`：c / barc
- `id`：唯一编号（全局唯一）
- `yL, yR`：$L_y(r),R_y(r)$（用于后续 IT 操作）

### 2.1.2 排序 key（建议写成一行公式，避免歧义）

例如可写为：

- `key(e) = (x(e), typeOrder(e), startClassOrder(e), id(e))`
- 其中 `typeOrder(END)=0 < typeOrder(START)=1`
- `startClassOrder` 仅在 START 时有效；END 可置 0 或忽略

这样就是严格全序，工程实现可以直接 `std::stable_sort` 或 `sort`。

------

## 2.2 IntervalTreeSampler：维护活跃 $y$ 区间的 Count / Report / Sample

对每一类分别维护一个结构：
$$
\mathcal{IT}^{(c)},\ \mathcal{IT}^{(\bar c)}.
$$
当扫描到 START(q) 时，我们对 **对立类的 $\mathcal{IT}$** 查询并采样。

### 2.2.1 对外接口（建议固定）

对一个查询区间 $Q=[a,b)$：

- `Insert(r)` / `Delete(r)`：插入/删除活跃区间 $I(r)=[L_y(r),R_y(r))$
- A 型：
  - `CountA(a)`
  - `ReportA(a)`
  - `SampleA(a, t')`：在 $A(Q)$ 上 i.i.d. 均匀有放回抽 $t'$
- B 型：
  - `CountB(a,b)`（对应 $a<\ell<b$）
  - `ReportB(a,b)`
  - `SampleB(a,b,t')`
- 便捷函数：
  - `CountOverlap(a,b)=CountA(a)+CountB(a,b)`
  - `ReportOverlap(a,b)=ReportA(a) ∪ ReportB(a,b)`（并集不重，因 A/B 互斥）
  - `SampleOverlap(a,b,t')`（可选：先按 |A|/|B| 分配，再分别 SampleA/SampleB）

> Sampling / Adaptive(大分支) 需要 Count+Sample；Enumerate/Adaptive(小分支) 需要 Report。

------

## 2.2.2 IT‑A：A 型 stabbing（$\ell\le a）

### 数据结构骨架：坐标压缩 + 段树 + 桶

1. 收集所有端点 $\{L_y(r),R_y(r)\}$（全体矩形的），排序去重得：
   $$
   Y=\{y_1<y_2<\cdots<y_m\}.
   $$

2. 建一棵平衡段树 $\mathcal{T}$，叶子对应原子区间 $[y_j,y_{j+1})$，共 $m-1$ 片叶子。

3. 对每个节点 $v$，维护桶 `BucketA(v)`，内部用

   - `vec`：致密数组（存 interval id 或指针）
   - `pos[id]`：该 id 在 vec 里的位置（用于 O(1) swap-delete）

> **工程建议**：
>
> - `pos` 最好做成“每个桶一个 hash map”会很重；更常见做法是：
>    对每个 interval 记录它插入了哪些节点，以及在各桶的位置（即 `handles` 列表），Delete 时沿 handles 逐桶 O(1) 删除。
> - 由于每条 interval 只进 $O(\log n)$ 个桶，handles 总长度 $O(n\log n)$ 可控。

### 规范覆盖（canonical cover）

对区间 $I=[\ell,r)$，令
$$
i=\mathrm{lower\_bound}(Y,\ell),\quad
j=\mathrm{lower\_bound}(Y,r),
$$
对应原子叶子范围 $[i, j)$。在段树上对该范围做“最小规范覆盖” $\mathcal{C}(I)$，得到 $O(\log n)$ 个节点，使得这些节点区间两两不交并覆盖 $[i,j)$。

- `Insert(I)`：对每个 $v\in\mathcal{C}(I)$，把 id 加入 `BucketA(v)`
- `Delete(I)`：同样从这些桶移出（用 handles 实现 O(1)）

### 关键性质：路径唯一归属（决定 Count/Report/Sample 的正确性）

对任意查询点 $a$，令 $leaf(a)$ 为包含 $a$ 的原子叶子；考虑从根到 $leaf(a)$ 的路径节点集合 $\mathrm{Path}(a)$。

**引理 3A（路径唯一归属）**
 任意活跃区间 $I$ 满足 $\ell\le a 当且仅当它在路径 $\mathrm{Path}(a)$ 的桶中 **恰出现一次**。

**证明直觉：** 规范覆盖把区间分成互不重叠的节点区间；点 $a$ 落在其中唯一一个覆盖节点上，该节点必在根到叶路径上，因此恰出现一次。

### 操作实现与复杂度

- `CountA(a)`：遍历 $\mathrm{Path}(a)$ 上所有桶大小求和，$O(\log n)$
- `ReportA(a)`：输出路径桶中的全部元素，$O(\log n + k)$
- `SampleA(a,t')`（批量 i.i.d.）：
  1. 取路径节点列表 $(v_1,\dots,v_s)$，$s=O(\log n)$
  2. 权重 $w_i=|BucketA(v_i)|$，总和 $W_A=\sum w_i$
  3. 若 $W_A=0$ 返回空
  4. 建一次 alias（或前缀和）在 $\{v_i\}$ 上按权重采样
  5. 重复 $t'$ 次：先按权重选桶，再桶内 `uniform_int(0, size-1)` 取元素
      → 总成本 $O(\log n + t')$

> 备注：若嫌 alias 常数大，因桶数 $s$ 只有 $O(\log n)$（通常 20~50），也可用前缀和 + 二分实现 $O(\log n + t'\log\log n)$；但为保持理论和代码一致，推荐 alias。

------

## 2.2.3 IT‑B：B 型 start-in-range（$a<\ell<b$）

B 型只与区间左端点 $\ell$ 有关，因此把每个活跃区间映射成一个点键 $(\ell,\text{id})$。

### 关键点 1：严格不等号必须通过二级键实现

目标条件是 **开区间** $(a,b)$，即 $a<\ell<b$。若只用 $\ell$ 排序，会在重复值与边界处出错；推荐用二级键：
$$
key(r)=(L_y(r), id(r))
$$
对所有区间的 key 排序得到秩数组，秩范围为 $[0, n)$。

- 左边界（排除 $\ell=a$）：
  $$
  L = \mathrm{upper\_bound}((a,+\infty))
  $$

- 右边界（排除 $\ell\ge b$，保留 $\ell）：
  $$
  R = \mathrm{lower\_bound}((b,-\infty))
  $$

则 $\ell\in(a,b)$ 等价于 $rank(\ell)\in[L,R)$。

### 数据结构骨架：秩域段树 + 路径桶

在秩域 $[0,n)$ 建一棵段树 $\mathcal{S}$。每个节点 $v$ 维护桶 `BucketB(v)`（同样用致密数组 + handles）。

- 插入一个点（对应某活跃 interval 的左端点秩 $p$）时：沿根到叶路径把 interval id 放入所有节点桶中。
- 删除：沿同一路径移除（用 handles 做 O(1)）。

### 查询：区间规范分解 + 不重计数

对查询区间 $[L,R)$，做段树规范分解得到 $O(\log n)$ 个两两不交的节点集合 $\mathcal{C}_{[L,R)}$。

**引理 3B（分解桶不交且完备）**
 对任意活跃点（某 interval 的 left endpoint 秩）$p\in[L,R)$，它在 $\{\mathrm{BucketB}(v)\mid v\in\mathcal{C}_{[L,R)}\}$ 中 **恰出现一次**。

**证明直觉：** 分解节点区间两两不交并覆盖 $[L,R)$，叶子 $p$ 属于其中唯一一个分解节点；而点插入沿根到叶路径，必包含该分解节点；且因为分解节点间无祖先关系，不会在两个分解节点同时出现。

### 操作实现与复杂度

- `CountB(a,b)`：
  1. 计算秩边界 $[L,R)$
  2. 求 $\mathcal{C}_{[L,R)}$
  3. 求和桶大小
      → $O(\log n)$
- `ReportB(a,b)`：枚举这些桶内容 → $O(\log n + k)$
- `SampleB(a,b,t')`（批量 i.i.d.）：
  1. 得到分解节点 $v_1,\dots,v_s$
  2. 权重 $w_i=|BucketB(v_i)|$，若总和为 0 返回空
  3. 建 alias
  4. 重复 $t'$ 次：按权重选桶 + 桶内均匀取元素
      → $O(\log n + t')$

------

## 2.2.4 IntervalTreeSampler 的总体代价（单类）

- 每个活跃区间：
  - 在 IT-A 中进入 $O(\log n)$ 个桶（规范覆盖）
  - 在 IT-B 中进入 $O(\log n)$ 个桶（根到叶路径）
- 所以单类维护空间 $O(n\log n)$，两类总 $O(n\log n)$。

------

## 2.3 Alias 与 slot（Sampling / Adaptive 大分支的离散规划）

Sampling / Adaptive 大分支需要把“全局 t 个样本”拆成“每个事件 e 的局部采样次数”。

### 2.3.1 事件级 alias

对 START 事件集合 $E$，权重为 $w_e$，总和 $W=\sum_{e\in E} w_e=|J|$。构造 alias 以实现：
$$
\Pr(E=e)=\frac{w_e}{W}.
$$
**边界必须写清：**

- 若 $W=0$：直接返回空（无相交对）
- alias 构建时可忽略 $w_e=0$ 的事件，或保留但永远不会被选中（实现上更建议过滤掉，减少 alias 大小）

### 2.3.2 A/B 级选择（固定事件 e 后）

事件 $e$ 对应查询 $Q=[a,b)$，定义
$$
w_e^A:=|A(Q)|,\qquad w_e^B:=|B(Q)|,\qquad w_e=w_e^A+w_e^B.
$$
对每个样本位置独立选择：
$$
T=\begin{cases}
A & \text{以概率 } w_e^A/w_e\\
B & \text{以概率 } w_e^B/w_e
\end{cases}
$$
**边界处理：**

- 若 $w_e^A=0$ 则必选 B；若 $w_e^B=0$ 则必选 A；
- $w_e=0$ 的事件不应被事件 alias 选中（否则除零）。

### 2.3.3 slot 结构

为每个 START 事件 e 维护两个 slot 列表：

- `S_e^A`: 需要从 A 集合采样并回填的位置集合
- `S_e^B`: 需要从 B 集合采样并回填的位置集合

> slot 的作用：把多次局部采样合并成一次 `Sample*(..., t_e^*)`，从而得到 $O(\log n + t_e)$ 而不是 $t_e\cdot O(\log n)$。

------

# 3. 算法详细流程（3 个版本）

## 3.0 公共预处理（所有版本共享）

输入：$R_c,R_{\bar c}$，样本数 $t$（Adaptive 还输入阈值 $J_\star$）。

1. 生成事件数组 `Events`（每矩形 2 个）。
2. 用 2.1 的 key 规则排序。
3. 初始化两类 IntervalTreeSampler：`ITc, ITb`（分别对应 $\mathcal{IT}^{(c)},\mathcal{IT}^{(\bar c)}$）：

- IT-A：建端点坐标压缩段树骨架（一次性）
- IT-B：建 left-endpoint 秩域段树骨架（一次性）
- 活跃集清空（桶为空）

------

## 3.1 版本一：Sampling‑IT（两遍扫描，不显式枚举）

### 3.1.1 总览

- **Phase1（第一次扫描）**：只计数，得到每个 START 事件的 $(w_e^A,w_e^B,w_e)$ 与全局 $W$
- **Phase2（离散规划）**：事件 alias + A/B 选择，生成 slot
- **Phase3（第二次扫描）**：按 slot 批量采样并回填输出

### 3.1.2 Phase1：第一次扫描（只计数）

初始化：

- `ITc.clear(), ITb.clear()`
- 对每个 START 事件 e 预留数组：`wA[e], wB[e], w[e]`
- `W = 0`

扫描 `Events`：

- 若事件为 `END(r)`：
  - 若 $r\in R_c$：`ITc.Delete(r)`
  - 若 $r\in R_{\bar c}$：`ITb.Delete(r)`
- 若事件为 `START(q)`（事件 id 为 e，$Q=[a,b)=[L_y(q),R_y(q))$）：
  - 若 $q\in R_c$：
    - `wA[e] = ITb.CountA(a)`
    - `wB[e] = ITb.CountB(a,b)`
  - 若 $q\in R_{\bar c}$`（对称）：
    - `wA[e] = ITc.CountA(a)`
    - `wB[e] = ITc.CountB(a,b)`
  - `w[e] = wA[e] + wB[e]`
  - `W += w[e]`
  - 插入本侧：
    - 若 $q\in R_c$：`ITc.Insert(q)` else `ITb.Insert(q)`

结束后若 `W==0`：返回空输出（或返回长度为 t 的空标记，视接口定义）。

### 3.1.3 Phase2：采样规划（事件 alias + A/B slot）

1. 在所有 START 事件中（建议只取 `w[e]>0` 的事件）建立 alias，使 $\Pr(E=e)=w[e]/W$。
2. 初始化 slot：对每个 e，`S_e^A = []`, `S_e^B = []`。
3. 对每个位置 $j=1..t$：
   - 抽事件 `e ~ Alias(w[e]/W)`
   - 抽类型：
     - 以概率 `wA[e]/w[e]` 选 A，否则 B（零边界按 2.3.2 处理）
   - 若 A：`S_e^A.push_back(j)`，否则 `S_e^B.push_back(j)`

### 3.1.4 Phase3：第二次扫描（局部批量采样 + slot 回填）

初始化：

- `ITc.clear(), ITb.clear()`（骨架复用，桶清空）
- `Ans[1..t]`

再次扫描 `Events`：

- `END(r)`：同 Phase1 删除
- `START(q)`（事件 e，$Q=[a,b)$）：
  - 令 `tA = |S_e^A|, tB = |S_e^B|`
  - 若 $q\in R_c$：
    - 若 `tA>0`：`ListA = ITb.SampleA(a, tA)`
    - 若 `tB>0`：`ListB = ITb.SampleB(a,b, tB)`
    - 按 slot 回填：
      - 对 `u=1..tA`: `Ans[S_e^A[u]] = (q, ListA[u])`
      - 对 `u=1..tB`: `Ans[S_e^B[u]] = (q, ListB[u])`
  - 若 $q\in R_{\bar c}$`（对称）：
    - `List*` 从 `ITc` 采样
    - 回填时注意顺序：
      - A/B 采样得到的是 $r_c$，所以写 `Ans[pos] = (r_c, q)`
  - 最后插入本侧：`IT(cls(q)).Insert(q)`

输出 `Ans`。

------

## 3.2 版本二：Enumerate+Sampling（一次扫描枚举全部 $J$ + 数组采样）

### 3.2.1 枚举阶段：一次扫描生成 `Pairs`

初始化：

- `Pairs = []`
- `ITc.clear(), ITb.clear()`

扫描 `Events`：

- `END(r)`：从对应 IT 删除
- `START(q)`（$Q=[a,b)$）：
  - 若 $q\in R_c$：
    - `ListA = ITb.ReportA(a)`
    - `ListB = ITb.ReportB(a,b)`
    - 对 `r in ListA ∪ ListB`：`Pairs.push_back((q,r))`
  - 若 $q\in R_{\bar c}$`：
    - 从 `ITc` report
    - 对 `r in ...`：`Pairs.push_back((r,q))`
  - 插入本侧 IT

结束后 `W = |Pairs|`。

### 3.2.2 采样阶段：数组上 i.i.d. 均匀采样

- 若 `W==0`：返回空
- 否则对 $j=1..t$：
  - `idx ~ UniformInt(0, W-1)`
  - `Ans[j] = Pairs[idx]`

输出 `Ans`。

> 该版本时间/空间都依赖 $|J|$；优点是实现最直接、常数小，适合 $J$ 很小的情形。

------

## 3.3 版本三：Adaptive+Sampling（阈值自适应）

目标：小 $|J|$ 用 Enumerate（一次扫描）；大 $|J|$ 用 Sampling（两遍），并保证分布严格正确。

### 3.3.1 Adaptive Phase1：一次扫描（永远 Count，必要时 Report 保存）

输入阈值 $J_\star$。

初始化：

- `mode = ENUMERATE`
- `AllPairs = []`（只在 ENUMERATE 模式维护，最多存到 $J_\star$）
- `ITc.clear(), ITb.clear()`
- `W = 0`
- 仍为每个 START 事件保存：`wA[e], wB[e], w[e]`（因为一旦切换，需要 Phase2/Phase3）

扫描 `Events`：

- `END(r)`：删除

- `START(q)`（事件 e，$Q=[a,b)$）：

  **(1) 先精确计数（无论枚举与否都做）**

  - 从对立类 IT 得到 `wA[e], wB[e]`
  - `w[e]=wA[e]+wB[e]`
  - `W += w[e]`

  **(2) 决定是否枚举**

  - 若 `mode==ENUMERATE` 且 `W <= J_star`：
    - `ListA = ReportA(...)`, `ListB = ReportB(...)`
    - 把 pair 追加进 `AllPairs`
  - 若 `mode==ENUMERATE` 且 `W > J_star`：触发切换
    - `mode = COUNT_ONLY`
    - `AllPairs.clear(); shrink_to_fit()`（释放内存）
    - 之后所有 START 不再 Report，只 Count
  - 若 `mode==COUNT_ONLY`：跳过枚举

  **(3) 插入本侧 IT**

结束后得到精确 `W=|J|`；若 `W==0`：返回空。

### 3.3.2 分支 A：未切换（$|J|\le J_\star$）

- 此时 `AllPairs` 已完整等于 $J$
- 在 `AllPairs` 上做数组均匀采样 → 返回

### 3.3.3 分支 B：已切换（$|J|>J_\star$）

- Phase1 已经保存了所有 `wA,wB,w` 与 `W`
- 执行与 Sampling‑IT 完全相同的：
  - Phase2（alias + slot）
  - Phase3（第二遍扫描 + 批量采样回填）

> 注意你草稿里强调的点非常关键：必须“先 Count 再判断切换”，否则会出现“某个 START 事件枚举一半就切换”的复杂纠错问题。Baseline++ 明确把它写成算法不变量。

------

# 4. Adaptive 版本阈值 $J_\star$ 的选择策略（更可落地的写法）

阈值本质是在权衡：

- Enumerate：一次扫描但要写出 $|J|$ 个 pair
- Sampling：两遍扫描但不写出 $|J|$

我们给三层策略（从硬约束到标定）：

## 4.1 内存硬约束（必须满足）

设：

- 总内存预算 `MemBudget` 字节
- 允许给 `AllPairs` 用的比例 $\rho\in(0,1)$
- 单个 pair 的存储开销 `sizeof(Pair)`（建议写清：两个 id + 容器开销的摊销）

则
$$
J_\star^{\text{mem}}=
\left\lfloor\frac{\rho\cdot\text{MemBudget}}{\text{sizeof(Pair)}}\right\rfloor,
\qquad
J_\star\le J_\star^{\text{mem}}.
$$

> 建议工程上取 $\rho\in[0.2,0.5]$，留足空间给段树桶、handles、slot 与系统开销。

## 4.2 时间建议（保证大分支不被“枚举到阈值”拖垮）

Adaptive 大分支（发生切换时）在 Phase1 最多额外枚举 $J_\star$ 条 pair（随后丢弃），因此：
$$
T_{\text{big}}=O(n\log n + J_\star + t).
$$
若希望大分支仍与 Sampling 的 $O(n\log n+t)$ 同阶，则选
$$
J_\star = O(n\log n + t).
$$
工程上可写成：
$$
J_\star^{\text{time}}=C_1\,n\log n + C_2\,t,
\qquad
J_\star=\min(J_\star^{\text{mem}},J_\star^{\text{time}}).
$$
其中 $C_1,C_2$ 由实现常数决定（见 4.3）。

## 4.3 交叉点标定（Benchmark 拟合，推荐写法）

在代表性数据分布下（真实数据或合成数据），对多组 $(n,t)$：

1. 跑纯 Enumerate+Sampling，记录运行时间随 $|J|$ 的变化。
2. 跑纯 Sampling‑IT（两遍），记录运行时间随 $(n,t)$ 的变化。
3. 找到两者时间相当的交叉点 $|J_{\text{cross}}|$。
4. 拟合

$$
|J_{\text{cross}}|\approx C_1\,n\log n + C_2\,t.
$$

1. 运行时取

$$
J_\star = \min\{\,0.8\cdot |J_{\text{cross}}|,\ J_\star^{\text{mem}}\}.
$$

（乘 0.8 是为了留出切换与内存碎片的安全边际。）

------

# 5. 算法分析（正确性 + 复杂度；三版本均给出）

## 5.1 关键引理汇总

### 引理 1（事件分块不重不漏）

$$
J=\biguplus_{e\in E}J_e,\qquad
|J|=\sum_{e\in E}w_e.
$$

证明见 1.3.4。

### 引理 2（A/B 分解互斥且完备）

$$
I\cap Q\neq\varnothing\iff (\ell\le a<r)\lor(a<\ell<b),
\quad
\mathrm{Overlap}(Q)=A(Q)\uplus B(Q).
$$

见 1.5。

### 引理 3（IT‑A/IT‑B 的桶分解保证“块内均匀采样”）

- IT‑A：路径桶对 $A(Q)$ 形成不交并（引理 3A）
- IT‑B：分解桶对 $B(Q)$ 形成不交并（引理 3B）

因此在 A 或 B 内 **按桶大小加权选桶 + 桶内均匀** 即得到集合上的均匀采样。

------

## 5.2 Sampling‑IT 的正确性（均匀性 + 独立性）

### 5.2.1 均匀性

固定一个输出位置 $j$，取任意 $P\in J$。由引理 1，$P$ 有唯一归属事件 $e^\star$，且 $P\in J_{e^\star}$。

算法生成 $Z_j$ 的过程是：

1. 事件级采样：$\Pr(E=e)=w_e/W$
2. 条件下 A/B：$\Pr(T=A\mid e)=w_e^A/w_e$，$\Pr(T=B\mid e)=w_e^B/w_e$
3. 在被选集合内均匀采样（由引理 3）

若 $P$ 属于 A（B 同理）：
$$
\Pr(Z_j=P)=
\Pr(E=e^\star)\cdot \Pr(T=A\mid e^\star)\cdot \Pr(P\mid e^\star,A)
=
\frac{w_{e^\star}}{W}\cdot\frac{w_{e^\star}^A}{w_{e^\star}}\cdot\frac{1}{w_{e^\star}^A}
=
\frac{1}{W}.
$$
因此对任意 $P\in J$，$\Pr(Z_j=P)=1/|J|$。

### 5.2.2 独立性（i.i.d.）

- Phase2 中每个位置 $j$ 的 $(E_j,T_j)$ 独立抽取；
- Phase3 对每个 $(e,A)$ 或 $(e,B)$ 的 `Sample*(..., t')` 必须实现为 **有放回且每次独立**；
- slot 回填只是重排输出位置，不引入随机耦合。

故 $Z_1,\dots,Z_t$ i.i.d.

------

## 5.3 Enumerate+Sampling 的正确性

### 5.3.1 枚举不重不漏

- 不漏：任意 pair $(r_c,r_{\bar c})$ 在更靠后的 START 事件时刻第一次同时活跃（引理 1），且 y 相交必落在 A 或 B 中（引理 2），因此会被 Report 输出并写入 `Pairs`。
- 不重：同一个 pair 不可能在更早 START 时出现；在同一个 START 下也不会同时落入 A 和 B（互斥）。因此不重复。

故 `Pairs` 与 $J$ 一一对应。

### 5.3.2 采样均匀且独立

数组上独立均匀下标采样天然给出 i.i.d. 均匀有放回。

------

## 5.4 Adaptive+Sampling 的正确性

- 若未切换：`AllPairs` 完整枚举 $J$，退化为 Enumerate+Sampling → 正确；
- 若已切换：Phase1 虽然丢弃枚举结果，但 **每个 START 的 $(w_e^A,w_e^B,w_e)$** 与全局 $W$ 仍是精确的；后续 Phase2+Phase3 与 Sampling‑IT 完全一致 → 正确；
- 切换不影响事件顺序与 Insert/Delete 时机，因此不会改变计数语义。

------

## 5.5 复杂度分析（时间 + 空间；三版本）

下面假设段树高度为 $O(\log n)$，并且：

- IT-A：Insert/Delete/Count/Report 的“树访问”是 $O(\log n)$
- IT-A/IT-B 的批量 Sample 为 $O(\log n + t')$

### 5.5.1 Sampling‑IT

- **时间**

  - 排序：$O(n\log n)$
  - Phase1：每事件 Insert/Delete/Count 为 $O(\log n)$ → $O(n\log n)$
  - Phase2：alias + t 次规划 → $O(n+t)$
  - Phase3：第二遍扫描 $O(n\log n)$ + 总采样 $O(t)$
     总计：

  $$
  T_{\text{Sampling}}=O(n\log n+t).
  $$

- **空间**

  - 两类 IT：$O(n\log n)$
  - 事件/权重数组：$O(n)$
  - slot + 输出：$O(t)$
     总计：

  $$
  S_{\text{Sampling}}=O(n\log n+t).
  $$

### 5.5.2 Enumerate+Sampling

- **时间**

  - 排序 + 扫描维护：$O(n\log n)$
  - Report 输出必写出 $|J|$ 个 pair：$\Theta(|J|)$
  - 数组采样：$O(t)$

  $$
  T_{\text{Enum}}=O(n\log n+|J|+t).
  $$

- **空间**

  - IT：$O(n\log n)$
  - `Pairs`：$\Theta(|J|)$
  - 输出：$O(t)$

  $$
  S_{\text{Enum}}=O(n\log n+|J|+t).
  $$

### 5.5.3 Adaptive+Sampling

- 若未切换（$W=|J|\le J_\star$）：
  $$
  T_{\text{Adaptive-A}}=O(n\log n+|J|+t),\quad
  S_{\text{Adaptive-A}}=O(n\log n+|J|+t).
  $$

- 若已切换（$W=|J|>J_\star$）：

  - Phase1 最多额外枚举 $J_\star$，其余只 Count：
    $$
    O(n\log n+J_\star)
    $$

  - Phase2+Phase3：$O(n\log n+t)$

  $$
  T_{\text{Adaptive-B}}=O(n\log n+J_\star+t),
  \quad
  S_{\text{Adaptive-B}}=O(n\log n+\max(J_\star,t)).
  $$

------

## 5.6 实现注意事项（建议作为“实现细节/坑点”小节保留）

1. **半开语义一致性**：
   - 相交判定用严格不等号；
   - 事件排序 END-before-START，否则贴边会被错误计入。
2. **B 型严格不等号 $a<\ell<b$**：
   - 必须用二级键 $(L_y,\text{id})$ 做 `upper_bound / lower_bound` 映射秩域，避免重复值与边界 bug。
3. **桶实现（O(1) 桶内抽样 + O(1) 删除）**：
   - 每桶用 `vector<int> vec` + swap-delete；
   - Delete 必须依赖 handles（每 interval 记录它在哪些节点桶、位置是多少），否则会退化到桶内线性查找。
4. **随机性与 i.i.d.**：
   - `SampleA/SampleB` 必须“每次独立、有放回”；
   - 不要把批量采样误实现成“不放回抽样”。
5. **零权重边界**：
   - `W==0` 直接返回空；
   - alias 仅对 `w[e]>0` 的事件建表；
   - A/B 选择要处理 `wA==0` 或 `wB==0`。