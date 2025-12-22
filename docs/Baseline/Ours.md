# Spatial Join 的 i.i.d. 均匀抽样（Sampling / Enumerat+Sampling / Adaptive+Sampling）

## 记号与约定

- 维度：$d\ge 2$，假设 $d$ 为常数（因此 $2^{d-1}$ 视为常数因子）。 

- 两类盒子集合：
  $$
  R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
  R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\},\quad
  n=n_1+n_2.
  $$
  
- 每个盒子是轴对齐 **半开**盒：
  $$
  r=\prod_{i=1}^d [L_i(r),R_i(r)),\quad L_i(r)<R_i(r).
  $$
  半开区间用于消除“贴边是否相交”的歧义，并与事件排序规则（同坐标 END 在 START 前）配套。 

- 统一输出 pair 的顺序：任何输出均为 $(r_c,r_{\bar c})$。当 START 来自 $R_{\bar c}$ 时需交换顺序写入输出。 

------

## 1. 问题定义与分析

### 1.1 相交 join 结果集合

我们只关心跨集合相交对（spatial join result）：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
两盒相交当且仅当每一维区间相交（半开区间下等价于严格不等式）：
$$
[L_i(r_c),R_i(r_c))\cap[L_i(r_{\bar c}),R_i(r_{\bar c}))\neq\varnothing
\iff
\max(L_i(\cdot))<\min(R_i(\cdot)).
$$


### 1.2 采样目标：在 $J$ 上 i.i.d. 均匀有放回

目标输出 $t$ 个样本：
$$
Z_1,\dots,Z_t\in J,
$$
要求 **i.i.d. 均匀有放回**：
$$
\Pr(Z_j=P)=\frac{1}{|J|}\ \ (\forall P\in J),\qquad Z_1,\dots,Z_t\ \text{相互独立}.
$$
并且我们的主版本（Sampling / Adaptive 大分支）希望**最坏情况**复杂度不依赖 $|J|$（因为 $|J|$ 可能是 $\Theta(n^2)$）。 

### 1.3 扫描维与事件系统：用第 1 维做 plane sweep

选择第 1 维 $x_1$ 做扫描。对每个盒子 $r$ 构造事件：

- $\texttt{START}(r)$ 在 $x_1=L_1(r)$
- $\texttt{END}(r)$ 在 $x_1=R_1(r)$

事件排序（稳定）：

1. 按 $x_1$ 升序；
2. 同 $x_1$：**END 在 START 前**（与半开区间一致，避免边界误计）；
3. 同 $x_1$ 且同为 START：固定规则（例如先处理 $R_c$ 的 START 再处理 $R_{\bar c}$），保证“唯一归属”；
4. 同 $x_1$ 且同为 END：相对顺序无关。

扫描到某个位置时，定义活跃集（隐式存于数据结构中）：
$$
\mathrm{Active}(R_c)=\{r_c\in R_c\mid x_0\in[L_1(r_c),R_1(r_c))\},\quad
\mathrm{Active}(R_{\bar c})\ \text{同理}.
$$


### 1.4 维度 $2..d$ 的 A/B 模式分解（关键）

对任意维度 $i\ge 2$，记
$$
D_i(r):=L_i(r),\quad U_i(r):=R_i(r).
$$
当处理某个 START 事件的新盒子 $q$ 与另一侧活跃盒子 $r$ 时，在维度 $i\ge2$ 的相交条件可拆成两种**互斥且覆盖**的情况：

- **A 型（lower 被刺穿）**
  $$
  D_i(r)\le D_i(q)<U_i(r)
  $$

- **B 型（对方 lower 落入我的开区间）**
  $$
  D_i(q)<D_i(r)<U_i(q)
  $$

因此每个相交对在维度 $2..d$ 上对应唯一模式向量：
$$
g=(g_2,\dots,g_d)\in\{A,B\}^{d-1}.
$$

### 1.5 事件分块：把全局 $J$ 分解成不交并（关键命题）

设 $E$ 为所有 START 事件集合。对任一 START 事件 $e$，令新盒子为 $r(e)$。若 $r(e)\in R_c$，则“另一侧”为 $R_{\bar c}$；反之亦然。

对每个模式 $g\in\{A,B\}^{d-1}$ 定义：
$$
K_e^{(g)}=\Bigl\{ r\in \mathrm{Active}(\text{Other})\ \Big|\ 
\forall i=2..d,\ 
\begin{cases}
g_i=A \Rightarrow D_i(r)\le D_i(r(e))<U_i(r),\\
g_i=B \Rightarrow D_i(r(e))<D_i(r)<U_i(r(e))
\end{cases}
\Bigr\}.
$$
并令
$$
w_e^{(g)}:=|K_e^{(g)}|,\qquad w_e:=\sum_g w_e^{(g)}.
$$
再定义局部 join 块（注意输出顺序统一为 $(r_c,r_{\bar c})$）：

- 若 $r(e)\in R_c$：
  $$
  J_e^{(g)}=\{(r(e),r_{\bar c})\mid r_{\bar c}\in K_e^{(g)}\}
  $$

- 若 $r(e)\in R_{\bar c}$：
  $$
  J_e^{(g)}=\{(r_c,r(e))\mid r_c\in K_e^{(g)}\}
  $$

令 $J_e=\biguplus_g J_e^{(g)}$。

**命题（不交分解）**
 在“半开区间 + END-before-START + START 固定 tie-break”下有：
$$
J=\biguplus_{e\in E} J_e,\qquad |J|=\sum_{e\in E} w_e.
$$
直观解释：每个相交对会在两者 START 中“更靠后”的那个 START 时第一次同时活跃，因此只被分配给唯一 $J_e$。

### 1.6 三个版本的定位

- **Sampling（Ours 三阶段，不枚举 $|J|$**）：最坏情况时间/空间不依赖 $|J|$。 
- **Enumerat+Sampling（Baseline，枚举全部 $J$ 再采样）**：实现直观，但时间/空间依赖 $|J|$。 
- **Adaptive+Sampling（自适应阈值）**：当 $|J|$ 小时退化为 Enumerate；当 $|J|$ 大时切换为 Sampling（并确保切换不破坏正确性）。 

------

## 2. 核心数据结构

本节描述两类核心组件：

1. 扫描事件系统（第 1 维）
2. 维度 $2..d$ 上的**高维模式数据结构** $\mathcal{D}_{m,g}$，其中 $m=d-1$

### 2.1 事件系统（plane sweep on $x_1$）

维护数组 `Events[1..2n]`，元素含：
$$
(x_1,\ \text{type}\in\{\texttt{START},\texttt{END}\},\ \text{class}\in\{c,\bar c\},\ \text{id}).
$$
按 1.3 的规则稳定排序。扫描中不显式维护 Active 列表，而是通过高维结构隐式维护“当前活跃盒子集合”。

### 2.2 高维模式结构 $\mathcal{D}_{m,g}$：接口与目标集合

令 $m=d-1$。把维度 $2..d$ 重新编号为 $1..m$（纯为书写方便）。对每个模式
$$
g\in\{A,B\}^{m}
$$
我们希望维护一个动态结构，支持对活跃盒子的如下操作：

- `Insert(r)` / `Delete(r)`：把盒子加入/移出活跃集
- `Count(q*)`：返回 $|S_g(q^*)|$
- `Sample(q*, t')`：从 $S_g(q^*)$ 中 **i.i.d. 均匀**采样 $t'$ 次
- `Report(q*)`（仅 Enumerate / Adaptive 的枚举分支需要）：枚举 $S_g(q^*)$ 全部元素

其中 $q^*$ 表示查询盒子在维度 $2..d$ 的投影；目标集合定义为：
$$
S_g(q)=\Bigl\{ r\in S_{\text{act}} :
\forall k=1..m,\ 
\begin{cases}
g_k=A \Rightarrow a_k(r)\le a_k(q)<b_k(r),\\
g_k=B \Rightarrow a_k(q)<a_k(r)<b_k(q)
\end{cases}\Bigr\}.
$$
这与前述 $K_e^{(g)}$ 完全一致（只不过把“Other 侧活跃集”抽象成 $S_{\text{act}}$）。

### 2.3 一维基础结构（作为递归构造的基石）

高维结构按维度递归嵌套实现，底层需要两种一维动态结构：

#### 2.3.1 一维 stabbing 结构（模式 A）

模型：动态维护活跃区间集合 $S_{\text{act}}\subseteq\{[\ell_i,r_i)\}$。对查询点 $q$，目标：
$$
J(q)=\{[\ell,r)\in S_{\text{act}}:\ \ell\le q<r\}.
$$
典型实现（骨架 + 活跃前缀桶）：

- 用所有端点建立平衡段树（叶为原子区间）。
- 每个区间 $I$ 在段树上有唯一的**最小规范覆盖** $\mathcal{C}(I)$，大小 $O(\log N)$。
- 对每个节点维护桶（dense array + handle），桶前缀表示活跃元素，支持：
  - 插入/删除：沿规范覆盖的 $O(\log N)$ 个桶，桶内 swap-delete，均摊 $O(1)$
  - 计数：沿查询点路径求和 $O(\log N)$
  - 采样：把路径上的桶大小当权重，在桶集合上做 alias，再桶内均匀取样
  - 批量采样 $t'$：用“slot/批量分配”避免每次递归，做到 $O(\log N + t')$

该结构是模式 A（$\ell\le q <r$）的一维实现。 

#### 2.3.2 一维 range 结构（模式 B）

模型：动态维护活跃点（用秩区分重复值）。对查询区间 $(\ell,r)$ 或 $[\ell,r)$，目标：
$$
J(q)=\{x\in X_{\text{act}}:\ x\in q\}.
$$
典型实现：

- 按 $(x_i, i)$ 排序得到秩数组（用二级键避免严格不等号问题）。
- 在秩域 $[0,N)$ 上建段树。
- 插入/删除更新根到叶路径 $O(\log N)$。
- 查询区间映射为秩区间 $[L,R)$，做规范分解得到 $O(\log N)$ 个节点桶：
  - 计数求和 $O(\log N)$
  - 采样同样用桶权重 alias，并可用 slot 做批量采样 $O(\log N+t')$

该结构对应模式 B（$\ell < x < r$）的一维实现。 

### 2.4 从 1 维到 $m$ 维：$\mathcal{D}_{m,g}$ 的递归构造要点

对 $m>1$，设尾模式 $g'=(g_2,\dots,g_m)$，核心思想是：

- 在第 1 维根据 $g_1$ 选择外层结构：
  - 若 $g_1=A$：用 **stabbing** 外层（区间刺点），节点对应规范覆盖/路径；
  - 若 $g_1=B$：用 **range** 外层（点落区间），节点对应秩段树的规范分解；
- 外层每个节点都挂一个 $(m-1)$ 维子结构 $\mathcal{D}_{m-1,g'}(\cdot)$，维护剩余维度上的投影；
- 查询时，外层把目标集合拆成若干个互不相交的子块（路径桶或规范分解桶）；对子块分别在子结构中 Count/Sample；再用“权重组合”实现全局 Count/Sample；批量采样用 slot 技巧把多次 Sample 合并为一次递归调用，从而得到 $+t'$ 而不是 $\times t'$。 

> **定理（高维模式数据结构复杂度）**
>  对常数 $m$ 与任意模式 $g\in\{A,B\}^m$，存在结构 $\mathcal{D}_{m,g}$ 支持：
>
> - Insert/Delete：$O((\log N)^m)$
> - Count：$O((\log N)^m)$
> - Sample $t'$ 次：$O((\log N)^m+t')$
> - Report（若需要）：$O((\log N)^m+k)$
> - 空间：$O(N(\log N)^m)$
>    其中 $k=|S_g(q)|$。

### 2.5 Alias 与 slot 容器（Sampling / Adaptive 大分支使用）

- **事件级 alias**：在 START 事件集合 $E$ 上按权重 $w_e$ 构造分布
  $$
  \Pr(E=e)=w_e/W,\quad W=\sum_{e\in E} w_e=|J|.
  $$
  单次采样 $O(1)$，构建 $O(|E|)=O(n)$。

- **模式级选择**：对固定事件 $e$，在常数大小的模式集合 $\{A,B\}^{d-1}$ 上按 $w_e^{(g)}$ 选择 $g$。因为模式数是常数，可直接线性选择或 alias，均视作 $O(1)$。

- **slot 列表**：对每个 $(e,g)$ 维护
  $$
  S_e^{(g)}\subseteq\{1,\dots,t\},
  $$
  表示哪些样本位置分配给该 $(e,g)$。第三阶段在处理 $e$ 时一次性采样 $|S_e^{(g)}|$ 个并回填这些位置。

------

## 3. 算法详细流程

这一章给出三个版本的**可执行级流程**。为避免重复，先给“公共预处理”，然后分别给三种主流程。

### 3.0 公共预处理（所有版本共享）

输入：$R_c, R_{\bar c}$，样本数 $t$（以及 Adaptive 的阈值 $J_\star$）。

1. 生成事件数组 `Events`（每盒 2 个事件），按规则排序。
2. 对两类集合分别、对每个模式 $g\in\{A,B\}^{d-1}$，构建 $\mathcal{D}^{(c)}_g$、$\mathcal{D}^{(\bar c)}_g$ 的“骨架”，并清空活跃集。预处理时间/空间为 $O(n(\log n)^{d-1})$。

------

### 3.1 版本一：Sampling（Ours 三阶段，不枚举 $|J|$）

该版本目标是最坏情形时间/空间不依赖 $|J|$。 

#### Phase 1：第一次扫描（只计数，得到所有 $w_e^{(g)}$、$w_e$、$W$）

初始化：`W = 0`。对事件按序扫描：

- `END(r)`：从其所属类别的所有模式结构中 `Delete(r)`

- `START(r)`（记为事件 $e$，新盒 $q=r(e)$）：

  1. 令 `Other` 为另一类集合；对每个模式 $g$：
     $$
     w_e^{(g)} \leftarrow \mathcal{D}^{(\text{Other})}_g.\texttt{Count}(q^*)
     $$
     并令 $w_e=\sum_g w_e^{(g)}$。保存 $\{w_e^{(g)}\},w_e$。

  2. `W += w_e`

  3. 将 $q$ 插入本侧所有模式结构：`Insert(q)`

扫描结束得到 $W=|J|$。若 $W=0$，直接返回空结果。 

#### Phase 2：采样规划（事件 alias + 模式分配 + slot 归档）

1. 在 START 事件集合 $E$ 上按权重 $w_e/W$ 构建 alias。
2. 对每个样本位置 $j=1..t$：
   - 抽事件：$E_j\sim w_e/W$
   - 条件下抽模式：$\Pr(G_j=g\mid E_j=e)=w_e^{(g)}/w_e$
   - 将 $j$ 加入 slot：$S_{E_j}^{(G_j)}\leftarrow S_{E_j}^{(G_j)}\cup\{j\}$

阶段二时间 $O(n+t)$，空间 $O(t)$。 

#### Phase 3：第二次扫描（按 slot 批量采样并回填 Ans）

1. 清空所有 $\mathcal{D}$ 的活跃集（骨架保留），准备输出数组 `Ans[1..t]`。

2. 按同一事件顺序再扫描一遍：

   - `END(r)`：同 Phase1 删除

   - `START(q=r(e))`：对每个模式 $g$：

     - 设 $t_e^{(g)}=|S_e^{(g)}|$

     - 若 $t_e^{(g)}>0$，调用
       $$
       \text{List}\leftarrow \mathcal{D}^{(\text{Other})}_g.\texttt{Sample}(q^*,\ t_e^{(g)})
       $$
       得到 $t_e^{(g)}$ 个 i.i.d. 均匀的“另一侧盒子”。

     - 将这些样本按 slot 写入 `Ans`（保持 pair 顺序为 $(r_c,r_{\bar c})$）：

       - 若 $q\in R_c$，回填 $(q, r)$
       - 若 $q\in R_{\bar c}$，回填 $(r, q)$

   - 最后 `Insert(q)` 到本侧所有模式结构

输出 `Ans[1..t]`。 

------

### 3.2 版本二：Enumerat+Sampling（先枚举 $J$ 再数组采样）

该版本是直接 baseline：先构造 `Pairs[0..|J|-1]`，再均匀取随机下标。采样天然 i.i.d.；关键是枚举要不重不漏。 

#### 枚举阶段：一次扫描 + Report 输出全部相交对

初始化：`Pairs = []`。扫描事件：

- `END(r)`：从本类所有模式结构中 `Delete(r)`

- `START(q=r(e))`：

  1. 令 `Other` 为另一类集合；对每个模式 $g$：
     $$
     \text{List}\leftarrow \mathcal{D}^{(\text{Other})}_g.\texttt{Report}(q^*)
     $$
     对 `List` 中每个盒子 $r$，向 `Pairs` 追加一条：

     - 若 $q\in R_c$：追加 $(q,r)$
     - 若 $q\in R_{\bar c}$：追加 $(r,q)$

  2. `Insert(q)` 到本侧所有模式结构

扫描结束后 `W = |Pairs| = |J|`。 

#### 采样阶段：数组上独立均匀下标采样（有放回）

- 若 `W==0`：返回空
- 否则对 $j=1..t$：
  - `idx ~ Unif{0,...,W-1}`
  - `Ans[j] = Pairs[idx]`

输出 `Ans`。 

------

### 3.3 版本三：Adaptive+Sampling（阈值 $J_\star$ 自适应切换）

Adaptive 用一次扫描同时做“精确计数”，并在 $|J|$ 小时顺便枚举保存；一旦发现累计相交对数会超过阈值 $J_\star$，立刻停止枚举并丢弃已枚举结果，转而只计数；最后走 Ours 的 Phase2+Phase3。 

#### Adaptive Phase1：单次扫描（始终计数 + 可选枚举）

初始化：

- `mode = ENUMERATE`
- `AllPairs = []`（最多允许增长到 $J_\star$，用于小分支）
- `W = 0`
- 所有模式结构活跃集清空（骨架已建好）

扫描事件：

- `END(r)`：同前，`Delete(r)`

- `START(q=r(e))`：

  1. **先计数（无论枚不枚举都做）**
      对每个模式 $g$：
     $$
     w_e^{(g)}\leftarrow \mathcal{D}^{(\text{Other})}_g.\texttt{Count}(q^*)
     $$
     令 $w_e=\sum_g w_e^{(g)}$，保存它们，并更新 `W += w_e`。

  2. **决定是否枚举**

     - 若 `mode==ENUMERATE` 且 `W <= J_star`：允许枚举该事件的全部对并写入 `AllPairs`：
        对每个模式 $g$：`Report` 并按 $(r_c,r_{\bar c})$ 追加到 `AllPairs`
     - 若 `mode==ENUMERATE` 但 `W > J_star`：
       - 立刻 `mode = COUNT_ONLY`
       - 释放/丢弃 `AllPairs`（后续不保留任何枚举结果）
       - 从此只 Count，不再 Report
     - 若 `mode==COUNT_ONLY`：跳过枚举

  3. `Insert(q)` 到本侧结构

扫描结束后得到精确 $W=|J|$ 及所有 $\{w_e^{(g)}\}$。 

#### 分支判定

- 若最终 `mode==ENUMERATE`（等价于 $|J|=W\le J_\star$）：
   直接在 `AllPairs` 上做数组采样（与 Enumerat+Sampling 相同）。 
- 若 `mode==COUNT_ONLY`（等价于 $|J|>J_\star$）：
   执行 Ours 的 Phase2 + Phase3（注意：Phase1 已经做完且计数已保存）。

------

## 4. Adaptive 版本阈值 $J_\star$ 的选择策略

阈值的核心作用：

- 小 $|J|$：枚举 + 数组采样更快（省掉第二遍扫描）
- 大 $|J|$：必须避免 `AllPairs` 爆内存/爆时间，及时切换到 Sampling

下面给出三层策略（从硬约束到软优化）。 

### 4.1 内存硬约束（必须满足）

设可用内存预算为 `MemBudget` 字节，允许给 `AllPairs` 使用比例 $\rho\in(0,1)$。每个 pair 的存储开销记为 `sizeof(Pair)`（例如两条 box-id/指针）。则：
$$
J_\star^{\text{mem}}=
\left\lfloor\frac{\text{MemBudget}\cdot \rho}{\text{sizeof(Pair)}}\right\rfloor.
$$
必须保证：
$$
J_\star \le J_\star^{\text{mem}}.
$$
否则会有 OOM 风险。 

### 4.2 保持“大分支”渐进最坏界的时间建议

在切换发生的情况下，Adaptive 在 Phase1 最多额外枚举 $O(J_\star)$ 个 pair（因为一旦 `W>J_star` 就停止 Report）。因此大分支总时间为：
$$
T_{\text{ours-branch}}=O\bigl(n(\log n)^{d-1}+J_\star+t\bigr).
$$
若希望大分支保持与 Sampling 相同的渐进：
$$
O\bigl(n(\log n)^{d-1}+t\bigr),
$$
则应选择
$$
J_\star=O\bigl(n(\log n)^{d-1}+t\bigr).
$$
实践写法常用：
$$
J_\star^{\text{time}}=C_1\,n(\log n)^{d-1}+C_2\,t,\qquad
J_\star=\min(J_\star^{\text{mem}},\ J_\star^{\text{time}}).
$$


### 4.3 工程标定（交叉点拟合，推荐）

做一组基准测试（不同数据分布/尺度/维度/样本数）：

1. 分别跑 “纯 Enumerat+Sampling” 与 “纯 Sampling”。

2. 找到两者耗时相当的交叉点 $J_{\text{cross}}$（以 $|J|$ 为横轴）。

3. 拟合
   $$
   |J_{\text{cross}}|\approx C_1\,n(\log n)^{d-1}+C_2\,t.
   $$

4. 运行时取略小于交叉点的 $J_\star$，并再受 $J_\star^{\text{mem}}$ 截断。 

------

## 5. 算法分析（正确性、复杂度等；三个版本都包含）

### 5.1 正确性基础：两层“不交并”分解

#### 5.1.1 模式分解的互斥与覆盖

对任意维 $i\ge2$，两区间相交当且仅当 A 或 B 之一成立，且二者互斥。由此每个相交对对应唯一 $g\in\{A,B\}^{d-1}$，因此对固定事件 $e$：
$$
J_e=\biguplus_g J_e^{(g)}.
$$

#### 5.1.2 事件分解 $J=\biguplus_e J_e$

由半开区间与事件排序（END-before-START + START tie-break），任一相交对会在两者 START 中“较晚”出现的那次第一次同时活跃，从而只会被计入唯一的 $J_e$，既不漏也不重：
$$
J=\biguplus_{e\in E}J_e,\qquad |J|=\sum_{e\in E}w_e.
$$

------

### 5.2 Sampling（Ours 三阶段）的正确性：全局 i.i.d. 均匀

设输出第 $j$ 个样本为 $Z_j=\text{Ans}[j]$。取任意 $P\in J$。由“事件不交分解 + 模式唯一性”，存在唯一 $(e^*,g^*)$ 使得 $P\in J_{e^*}^{(g^*)}$。

算法对位置 $j$ 的生成可拆成三段概率：

1. 抽事件：$\Pr(E_j=e)=w_e/W$
2. 条件下抽模式：$\Pr(G_j=g\mid E_j=e)=w_e^{(g)}/w_e$
3. 在 $K_e^{(g)}$ 上局部均匀采样：$\Pr(P\mid e,g)=1/w_e^{(g)}$

于是：
$$
\Pr(Z_j=P)=\frac{w_{e^*}}{W}\cdot\frac{w_{e^*}^{(g^*)}}{w_{e^*}}\cdot\frac{1}{w_{e^*}^{(g^*)}}
=\frac{1}{W}=\frac{1}{|J|}.
$$
因此每个位置单独均匀。

独立性来自两点：

- Phase2 中 $(E_j,G_j)$ 是对每个 $j$ 独立生成；
- Phase3 中对每个 $(e,g)$ 批量采样得到的序列是 i.i.d.，且不同 $(e,g)$ 之间独立；每个位置 $j$ 只读取其唯一对应 $(E_j,G_j)$ 的一个样本。

故 $Z_1,\dots,Z_t$ 为 i.i.d. 均匀有放回样本。 

------

### 5.3 Enumerat+Sampling 的正确性：i.i.d. 均匀

1. 枚举阶段的正确性（不重不漏）：

- 不漏：任意相交对会在“较晚 START”时刻被发现；当处理该 START 时另一侧盒子必活跃且在某个唯一模式 $g$ 下满足约束，因此会被 `Report` 输出并加入 `Pairs`。
- 不重：更早的 START 不可能看到该 pair；同一 START 下不同模式块互斥，不会重复输出。

因此 `Pairs` 是 $J$ 的一次一一枚举。 

1. 采样阶段：对 `Pairs[0..|J|-1]` 独立均匀抽下标，自然得到 i.i.d. 均匀有放回样本。 

------

### 5.4 Adaptive+Sampling 的正确性：切换不改变分布

Adaptive 的关键点：**Phase1 始终做精确 Count 并保存 $\{w_e^{(g)}\}$**；切换只影响“是否将 `Report` 的枚举结果保存”。因此：

- 若未切换：`AllPairs` 完整枚举 $J$，等价于 Enumerat+Sampling → 正确；
- 若已切换：丢弃所有枚举结果，但保留精确计数，随后执行与 Sampling 相同的 Phase2+Phase3 → 正确；

切换不会影响事件顺序、Insert/Delete 时刻、以及每个事件的计数值，因此不会引入任何分布偏差。

------

### 5.5 复杂度分析（时间与空间；三个版本分别给出）

以下令 $m=d-1$。

#### 5.5.1 Sampling（Ours 三阶段）

- **时间**

  - 预处理建骨架：$O(n(\log n)^m)$

  - Phase1 扫描：每事件插删/Count 为 $O((\log n)^m)$，总 $O(n(\log n)^m)$

  - Phase2 规划：$O(n+t)$

  - Phase3 扫描 + 批量采样：
    $$
    \sum_{e,g} O((\log n)^m+t_e^{(g)})=O(n(\log n)^m+t)
    $$

  合并：
  $$
  T_{\text{Sampling}}=O\bigl(n(\log n)^{d-1}+t\bigr).
  $$
  
- **空间**

  - 所有模式、两类集合的结构：$O(n(\log n)^m)$（模式数常数）
  - slot + 输出：$O(t)$
     故：

  $$
  S_{\text{Sampling}}=O\bigl(n(\log n)^{d-1}+t\bigr).
  $$

  

#### 5.5.2 Enumerat+Sampling

- **时间**

  - 预处理与结构维护：$O(n(\log n)^m)$
  - 枚举输出：$|J|$ 条 pair 必须写出，代价 $O(|J|)$
  - 数组采样：$O(t)$
     故：

  $$
  T_{\text{Enum}}=O\bigl(n(\log n)^{d-1}+|J|+t\bigr).
  $$

  

- **空间**

  - 数据结构：$O(n(\log n)^m)$
  - 显式存 `Pairs`：$\Theta(|J|)$
  - 输出：$O(t)$
     故：

  $$
  S_{\text{Enum}}=O\bigl(n(\log n)^{d-1}+|J|+t\bigr).
  $$

  

#### 5.5.3 Adaptive+Sampling

分两种情况： 

- **Case A：未切换（$|J|\le J_\star$）**

  - Phase1：计数 + 枚举，总时间 $O(n(\log n)^{d-1}+|J|)$
  - 数组采样：$O(t)$

  $$
  T_{\text{Adaptive-A}}=O\bigl(n(\log n)^{d-1}+|J|+t\bigr),
  \quad
  S_{\text{Adaptive-A}}=O\bigl(n(\log n)^{d-1}+|J|+t\bigr).
  $$

  

- **Case B：发生切换（$|J|>J_\star$）**

  - Phase1：最多枚举 $O(J_\star)$ 条（之后只 Count），时间
    $$
    O\bigl(n(\log n)^{d-1}+J_\star\bigr)
    $$

  - Phase2 + Phase3：与 Sampling 相同，$O(n(\log n)^{d-1}+t)$
     合并：

  $$
  T_{\text{Adaptive-B}}=O\bigl(n(\log n)^{d-1}+J_\star+t\bigr).
  $$

  空间峰值：
  $$
  S_{\text{Adaptive-B}}=O\bigl(n(\log n)^{d-1}+\max(J_\star,t)\bigr).
  $$
  若按建议取 $J_\star=O(n(\log n)^{d-1}+t)$，则大分支仍保持与 Sampling 同阶的渐进上界。