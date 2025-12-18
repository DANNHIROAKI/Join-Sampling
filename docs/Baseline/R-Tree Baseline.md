## 1. 问题定义与分析

### 1.1 输入：两类半开轴对齐盒子

在 $d\ge 2$ 维欧氏空间 $\mathbb{R}^d$ 中，给定两类轴对齐半开盒子集合：
$$
R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\},\quad n=n_1+n_2.
$$
每个盒子
$$
r=\prod_{i=1}^d [L_i(r),R_i(r)),\quad L_i(r)<R_i(r).
$$
我们只关心**跨集合相交对**：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
半开区间下，两盒子相交当且仅当每一维投影区间都相交：
$$
[a,b)\cap[c,d)\neq\varnothing \iff \max(a,c)<\min(b,d).
$$

------

### 1.2 输出：对 $J$ 的 i.i.d. 均匀有放回采样

目标输出 $t$ 个样本
$$
Z_1,\dots,Z_t\in J
$$
要求 **i.i.d. 均匀有放回**：
$$
\Pr(Z_j=P)=\frac{1}{|J|}\quad(\forall P\in J),\qquad Z_1,\dots,Z_t\ \text{相互独立}.
$$

------

### 1.3 第 1 维 Plane Sweep：事件系统与“唯一归属”的关键性质

选第 1 维 $x_1$ 做扫描。每个盒子 $r$ 产生两个事件：

- $\texttt{START}(r)$ at $x_1=L_1(r)$
- $\texttt{END}(r)$ at $x_1=R_1(r)$

**稳定排序规则（非常关键，匹配半开边界）：**

1. 按 $x_1$ 升序；
2. 同一 $x_1$：$\texttt{END}$ 在 $\texttt{START}$ 之前（避免“贴边误相交”）；
3. 同一 $x_1$ 且都是 $\texttt{START}$：用固定 tie-break（例如先处理 $R_c$ 的 START，再处理 $R_{\bar c}$），保证“唯一归属”。

扫描到位置 $x_0$ 时定义活跃集：
$$
\mathrm{Active}(R_c)=\{r_c\in R_c\mid x_0\in[L_1(r_c),R_1(r_c))\},\quad
\mathrm{Active}(R_{\bar c})\ \text{同理}.
$$
**关键观察 1（维度分离）：**
 当处理某个 $\texttt{START}(q)$ 时，当前 $x_0=L_1(q)$。对任一对立活跃盒 $r\in \mathrm{Active}(\text{Other})$，由于 END 已先于 START 处理且区间半开，$r$ 与 $q$ 在第 1 维必然“真重叠”（不是仅贴边）。因此判断整盒相交只需检查维度 $2..d$ 是否相交。

------

### 1.4 事件诱导的分块：$J=\biguplus_e J_e$

令 $E$ 为所有 START 事件集合。对每个 START 事件 $e$，设其对应的新盒子为 $q=r(e)$，定义对立活跃集的相交对象集合：
$$
K_e := \{ r\in \mathrm{Active}(\text{Other}) \mid q\cap r\neq\varnothing\},\qquad w_e:=|K_e|.
$$
为了统一输出 pair 顺序为 $(r_c,r_{\bar c})$，定义：

- 若 $q\in R_c$：$J_e=\{(q,r_{\bar c})\mid r_{\bar c}\in K_e\}$
- 若 $q\in R_{\bar c}$：$J_e=\{(r_c,q)\mid r_c\in K_e\}$

在上述事件排序下成立：
$$
J=\biguplus_{e\in E} J_e,\qquad |J|=\sum_{e\in E} w_e.
$$
直觉：每个相交对会在两者 START 事件中“更靠后”的那一个时刻首次同时活跃，从而恰好归属到唯一的 $J_e$。

------

## 2. 核心数据结构

R-tree 在这里承担的是：动态维护“对立活跃集”在维度 $2..d$ 上的投影，并支持 **Report / Count / Sample** 三类查询接口。

### 2.1 事件数组与盒子句柄

- `Events[1..2n]`：事件元素 $(x_1,\text{type},\text{class},\text{box-id})$
- `type ∈ {START, END}`, `class ∈ {c, \bar c}`
- 每个盒子有唯一 `box-id`，用于插入/删除时定位（最好还能保存 `handle` 指向叶条目位置以加速删除）

------

### 2.2 只存维度 $2..d$ 的投影：$r^\star$

因为第 1 维由 sweep 保证重叠，我们在 R-tree 里存盒子在维度 $2..d$ 的 $(d-1)$ 维投影：
$$
r^\star := \prod_{i=2}^d [L_i(r),R_i(r)).
$$
维护两棵动态 R-tree：

- $\mathcal{RT}^{(c)}$：存 $\mathrm{Active}(R_c)$ 的投影 $\{r_c^\star\}$
- $\mathcal{RT}^{(\bar c)}$：存 $\mathrm{Active}(R_{\bar c})$ 的投影 $\{r_{\bar c}^\star\}$

每个叶条目：$(\text{box-id}, r^\star)$。

------

### 2.3 半开语义的精确相交判定（叶子过滤必需）

对两个 $(d-1)$ 维半开盒
 $A=\prod_{k=2}^d[\alpha_k,\beta_k)$、$B=\prod_{k=2}^d[\gamma_k,\delta_k)$，定义精确判定：
$$
\texttt{Intersect}^\star(A,B)\iff \forall k=2..d:\ \max(\alpha_k,\gamma_k)<\min(\beta_k,\delta_k).
$$
注意：R-tree 内部节点用 MBR 做剪枝时，可以用“保守相交”（避免漏报），但**输出/计数/采样必须在叶子用 $\texttt{Intersect}^\star$ 过滤**，否则会把“贴边”当作相交。

------

### 2.4 R-tree 需要支持的接口

下面按三个版本的需求列出接口（从弱到强）：

#### 2.4.1 动态维护：Insert / Delete

- `Insert(r)`：插入条目 $(id(r), r^\star)$
- `Delete(r)`：删除条目（最好利用 handle；否则需 FindLeaf 搜索）

#### 2.4.2 枚举：ReportIntersect(Q)

`ReportIntersect(Q)`：返回所有活跃条目 $r$ 使得 $r^\star\cap Q\neq\varnothing$。

实现：标准 R-tree 相交搜索（MBR 剪枝） + 叶子用半开精确判定过滤。用于 **Enumerate+Sampling** 与 **Adaptive 的枚举分支**。

#### 2.4.3 计数：CountIntersect(Q)

`CountIntersect(Q)`：精确返回
$$
|\{r\in \mathrm{Active}: r^\star\cap Q\neq\varnothing\}|.
$$
典型递归规则（精确、带剪枝）：

- 若 `MBR(node)` 与 $Q$ 不相交：贡献 0
- 若 `MBR(node) ⊆ Q`：直接贡献 `size(node)`（整子树必命中）
- 若叶：逐条用 $\texttt{Intersect}^\star$ 检查并计数
- 否则递归子节点求和

其中 `size(node)` 是子树叶条目数，需随 Insert/Delete 维护。

#### 2.4.4 采样：SampleIntersect(Q, t′)

`SampleIntersect(Q, t′)`：从命中集合
$$
K(Q)=\{r\in \mathrm{Active}: r^\star\cap Q\neq\varnothing\}
$$
中**均匀有放回**抽取 $t′$ 个 box-id（i.i.d）。

统一写法：**“前沿/引导结构（Guide/Frontier）+ 按计数加权选择 + 子树内均匀下降”**。两份文档里分别以 “BuildGuide” 或 “Frontier 分解”描述，本质一致：先构造一个覆盖命中集合的、不重叠的组件分解，每个组件带权重，然后做加权抽取，组件内再均匀抽。

**一种可直接落地的批量采样实现（推荐）：**

- **Step A：构造前沿组件集合 $\mathcal{F}$**
   DFS 遍历 R-tree：

  - 若 `MBR(node)` 与 $Q$ 不相交：跳过
  - 若 `MBR(node) ⊆ Q`：把该 node 作为一个 **FullyAccepted 组件**加入 $\mathcal{F}$，权重 $w=\text{size(node)}$，停止向下展开
  - 若叶：扫描叶条目，收集命中条目列表 `HitList` 作为 **Leaf 组件**，权重 $w=|\text{HitList}|$
  - 若内部：递归孩子

  得到
  $$
  K(Q)=\biguplus_{C\in\mathcal{F}} S_C,\qquad w_C=|S_C|,\qquad W=\sum_{C}w_C=|K(Q)|.
  $$

- **Step B：组件级别抽样（alias/前缀和）并分配 slot**
   需要 $t′$ 次样本：独立抽组件 $C$ ，满足 $\Pr(C)=w_C/W$，并统计每个组件被抽到次数 $t_C$。

- **Step C：组件内部均匀抽样**

  - 若 $C$ 是 FullyAccepted node：执行 `SampleSubtreeUniform(node)`，每次在子树内按 `size(child)` 比例随机下降到叶，再在叶内均匀选条目（单次 $O(h)$，$h$ 为树高）
  - 若 $C$ 是 Leaf 组件：直接在 `HitList` 中均匀选一个（单次 $O(1)$）

  重复生成 $t_C$ 次，汇总得到 $t′$ 个 i.i.d 样本。

------

## 3. 算法详细流程（3 个版本）

下面所有版本共享相同的事件系统与两棵 R-tree 维护方式，只是“何时 Report / Count / Sample，以及是否需要第二遍扫描”不同。

------

### 3.1 版本 A：Enumerate+Sampling（显式枚举 $J$ + 数组采样）

**适用场景：** $|J|$ 小（否则时间/空间都被 $|J|$ 主导）。

#### 3.1.1 输入输出

- 输入：$R_c,R_{\bar c}$，样本数 $t$
- 输出：$t$ 个 i.i.d. 均匀样本 $(r_c,r_{\bar c})\in J$

#### 3.1.2 算法步骤

**预处理：**

1. 构造 `Events`：每盒子产生 START/END；
2. 按 §1.3 排序；
3. 初始化两棵空 R-tree：$\mathcal{RT}^{(c)},\mathcal{RT}^{(\bar c)}$
4. 初始化动态数组 `Pairs=[]`

**一次扫描（枚举所有相交对）：**

对每个事件按序处理：

- 若 `END(r)`：
  - 若 $r\in R_c$：$\mathcal{RT}^{(c)}.\texttt{Delete}(r)$
  - 若 $r\in R_{\bar c}$：$\mathcal{RT}^{(\bar c)}.\texttt{Delete}(r)$
- 若 `START(q)`：
  - 若 $q\in R_c$：
    1. `Cand ← RT^(barc).ReportIntersect(q^*)`
    2. 对每个候选 $r_{\bar c}\in$Cand：若 $\texttt{Intersect}^\star(r_{\bar c}^\star,q^\star)$ 为真，则 `Pairs.append((q,r_barc))`
    3. `RT^(c).Insert(q)`
  - 若 $q\in R_{\bar c}$（对称）：
    1. `Cand ← RT^(c).ReportIntersect(q^*)`
    2. 过滤后追加 `(r_c, q)`（保持 pair 顺序统一）
    3. `RT^(barc).Insert(q)`

扫描结束：$W=|\texttt{Pairs}|=|J|$。

**数组采样：**

- 若 $W=0$：返回空
- 否则对 $j=1..t$：生成 $\mathrm{idx}_j\sim\mathrm{Unif}\{0,\dots,W-1\}$，输出 `Pairs[idx_j]`。

------

### 3.2 版本 B：Sampling（不枚举 $J$，两遍扫描 + 事件级 alias + 局部均匀采样）

**适用场景：** $|J|$ 大或不可枚举时，仍要得到 $J$ 上全局 i.i.d. 均匀样本。

#### 3.2.1 Phase 1：第一次扫描（只计数）

目标：计算每个 START 事件的权重 $w_e=|J_e|$，并累计
$$
W=\sum_{e\in E} w_e=|J|.
$$
初始化：$\mathcal{RT}^{(c)},\mathcal{RT}^{(\bar c)}\gets\emptyset$，`W=0`，并为每个 START 事件存储 `w_e`。

扫描规则：

- `END(r)`：与版本 A 相同删除
- `START(q)`：
  - 若 $q\in R_c$：
     $w_e \leftarrow \mathcal{RT}^{(\bar c)}.\texttt{CountIntersect}(q^\star)$，`W += w_e`，再 `RT^(c).Insert(q)`
  - 若 $q\in R_{\bar c}$`：对称

若最终 $W=0$，直接返回空。

#### 3.2.2 Phase 2：事件级 alias + slot 分配

在 START 事件集合 $E$ 上按
$$
p_e=\frac{w_e}{W}
$$
构造 alias 表。

为每个事件维护 slot 列表 $S_e\subseteq\{1,\dots,t\}$。对每个样本位置 $j$：

1. 抽事件 $E_j\sim p_e$
2. 把位置 $j$ 放入 $S_{E_j}$

令 $t_e:=|S_e|$，有 $\sum_e t_e=t$。

#### 3.2.3 Phase 3：第二次扫描（按事件做局部批量采样 + 回填）

重置两棵 R-tree 为空，初始化 `Ans[1..t]`。

第二遍扫描：

- `END(r)`：删除同上
- `START(q)`（事件 $e$）：
  - 若 $t_e=|S_e|=0$：只插入，不采样
  - 若 $t_e>0$：
    - 若 $q\in R_c$：
       `List ← RT^(barc).SampleIntersect(q^*, t_e)` 得到 $t_e$ 个 i.i.d. 的 $r_{\bar c}$，按 slot 顺序回填 `Ans[j_u]=(q, List[u])`
    - 若 $q\in R_{\bar c}$`：对称回填为 `(List[u], q)`
  - 最后插入 $q$ 到本侧 R-tree

输出 `Ans`。

------

### 3.3 版本 C：Adaptive+Sampling（阈值自适应：小 $J$ 枚举，大 $J$ 两遍采样）

**适用场景：** 你希望“小输出时低常数、只扫一遍；大输出时避免枚举爆炸”。

#### 3.3.1 输入输出

- 输入：$R_c,R_{\bar c}$，样本数 $t$，阈值 $J_\star$
- 输出：$t$ 个 i.i.d. 均匀样本

#### 3.3.2 Phase 1：一次扫描（永远 Count；在阈值内可选 Report 枚举）

初始化：

- `mode = ENUMERATE`
- `AllPairs = []`（最多保存到 $J_\star$）
- `W=0`（累计 $|J|$）
- 两棵 R-tree 清空
- 记录每个 START 的权重 `w_e`

扫描每个事件：

- `END(r)`：删除同前

- `START(q)`（事件 $e$）：

  **Step A：先精确计数（无论是否枚举都要）**

  - 若 $q\in R_c$：$w_e \leftarrow \mathcal{RT}^{(\bar c)}.\texttt{CountIntersect}(q^\star)$
  - 否则对称
     更新 `W += w_e` 并存储 `w_e`。

  **Step B：阈值判断 + 可选枚举**

  - 若 `mode==ENUMERATE` 且 $W\le J_\star$：允许枚举本事件贡献
    - 若 $q\in R_c$：`List ← RT^(barc).ReportIntersect(q^*)`，把每个命中追加到 `AllPairs` 形成 $(q,r_{\bar c})$
    - 若 $q\in R_{\bar c}$`：对称追加 $(r_c,q)$
  - 若 `mode==ENUMERATE` 但 $W>J_\star$：**触发切换**
    - `mode ← COUNT_ONLY`
    - 丢弃/释放 `AllPairs`
    - 后续不再调用 `ReportIntersect`（只 Count）
  - 若 `mode==COUNT_ONLY`：跳过枚举

  **Step C：插入 $q$**

  - 插入到本侧 R-tree

扫描结束得到精确 $W=|J|$。

> 重要细节：**“先 Count 再判断阈值”** 的写法可以保证不会出现“一个 START 事件枚举到一半才超阈值”的尴尬，同时 Phase1 枚举输出量 $\le J_\star$。

#### 3.3.3 分支选择与后续流程

- 若 $W=0$：返回空
- 若最终 `mode==ENUMERATE`（等价于 $W\le J_\star$ 且未切换）：
  - `AllPairs` 已完整枚举 $J$
  - 直接做数组均匀采样（同版本 A），**无需第二遍扫描**
- 若最终 `mode==COUNT_ONLY`（发生切换，$W>J_\star$）：
  - 执行版本 B 的 Phase2（事件 alias + slot）与 Phase3（第二遍扫描 + SampleIntersect 回填）



------

## 4. Adaptive 版本阈值 $J_\star$ 的选择策略

阈值的目标是：

- **不爆内存**（枚举数组 AllPairs 可控）
- **不浪费时间**（大 $J$ 时不要无意义枚举太多）

### 4.1 内存硬约束（必须满足）

给定内存预算 `MemBudget`，给 `AllPairs` 的比例 $\rho\in(0,1)$，每个 pair 存储开销 `sizeof(Pair)`（两个 id/指针），则：
$$
J_\star^{\text{mem}}
=
\left\lfloor
\frac{\text{MemBudget}\cdot\rho}{\text{sizeof(Pair)}}
\right\rfloor.
$$
必须取 $J_\star \le J_\star^{\text{mem}}$。

------

### 4.2 时间权衡（建议的“成本模型”）

在 **切换分支**（$W>J_\star$）中，Phase1 最多枚举 $J_\star$ 个 pair，形成额外开销 $O(J_\star)$（之后就只 Count），其余开销与“纯 Sampling”同阶。

因此一个实用原则是让 $J_\star$ 不超过“两遍扫描几何查询的量级”。文档给出的经验写法是：
$$
J_\star^{\text{time}} \approx C_1\cdot(\text{一次扫描的平均几何访问成本规模}) + C_2\cdot t,
\qquad
J_\star=\min(J_\star^{\text{mem}},\ J_\star^{\text{time}}).
$$
其中常数 $C_1,C_2$ 用基准测试标定。

------

### 4.3 工程标定：交叉点拟合（最推荐）

建议做 benchmark：

1. 固定 $d$，多组 $n,t$，多种分布（均匀/聚簇/长条等）
2. 分别跑：
   - “纯 Baseline 枚举+采样”（版本 A）
   - “纯 Sampling 两遍采样”（版本 B）
3. 找耗时相当的交叉点 $|J_{\text{cross}}|$
4. 设置 $J_\star \approx 0.8\cdot |J_{\text{cross}}|$，并再用 $J_\star^{\text{mem}}$ 截断

这样 $J_\star$ 会自动贴近你实现与数据分布下的真实最优分界。

------

## 5. 算法分析（正确性 + 复杂度，三版本都包含）

> 本节把三版本共用的正确性“骨架”先说明，再分别落到每个版本的结论与复杂度表达。

### 5.1 正确性分析

#### 5.1.1 共同基础：START 事件处只需检查维度 $2..d$

处理 $\texttt{START}(q)$ 时，对立活跃盒子与 $q$ 在第 1 维必然重叠（END-before-START + 半开语义），因此整盒相交等价于投影 $q^\star$ 与 $r^\star$ 在维度 $2..d$ 相交。

#### 5.1.2 共同基础：事件分块不交分解 $J=\biguplus_e J_e$

对任意相交对 $(r_c,r_{\bar c})\in J$，在事件全序下存在唯一“更靠后”的 START 事件 $e^\star$。处理 $e^\star$ 时另一个盒子必在对立活跃集中，因此该 pair 只会归属到唯一的 $J_{e^\star}$。从而：
$$
J=\biguplus_{e\in E} J_e,\qquad |J|=\sum_e w_e.
$$

#### 5.1.3 R-tree 查询语义正确性（Report/Count/Sample）

- `ReportIntersect`：R-tree 搜索保证不漏（无 false negative）；叶子用半开精确判定过滤后输出集合恰为命中集合。
- `CountIntersect`：用“不相交剪枝 + MBR⊆Q 整子树计数 + 叶子精确判定”，计数精确等于命中数。
- `SampleIntersect`：通过“前沿分解成不交组件 + 按组件大小加权选组件 + 组件内均匀抽取”，保证每个命中条目概率为 $1/|K(Q)|$，重复独立抽取得到 i.i.d.。

------

### 5.2 版本 A：Enumerate+Sampling 的正确性

**结论：** 输出是 $J$ 上 i.i.d. 均匀有放回样本。

- 枚举正确（不重不漏）：任意 pair 只会在“更靠后 START”时刻被报告并追加一次，因此 `Pairs` 精确枚举 $J$。
- 采样正确：对 `Pairs[0..|J|-1]` 做独立均匀下标采样，显然得到 i.i.d. 均匀分布。

------

### 5.3 版本 B：Sampling 的正确性

**结论：** 输出是 $J$ 上 i.i.d. 均匀有放回样本。

关键链条：

1. Phase1 的 `CountIntersect` 精确得到每个事件权重 $w_e=|J_e|$，且 $W=\sum_e w_e=|J|$。

2. Phase2 抽事件 $E_j$ 的分布为 $\Pr(E_j=e)=w_e/W$。

3. Phase3 在给定事件 $e$ 下，`SampleIntersect(q^*, t_e)` 在 $K_e$ 上均匀有放回抽样，因此对任意 $P\in J_e$，
   $$
   \Pr(P\mid E_j=e)=\frac{1}{w_e}.
   $$

4. 对任意 $P\in J$，存在唯一 $e^\star$ 使 $P\in J_{e^\star}$，于是
   $$
   \Pr(Z_j=P)=\Pr(E_j=e^\star)\cdot\Pr(P\mid E_j=e^\star)=\frac{w_{e^\star}}{W}\cdot\frac{1}{w_{e^\star}}=\frac{1}{|J|}.
   $$

5. 独立性：事件抽样对不同 $j$ 独立；局部采样按“有放回 + 独立随机数”生成；slot 回填不引入相关性，因此全局 i.i.d.。

------

### 5.4 版本 C：Adaptive+Sampling 的正确性

**结论：** 无论是否触发切换，输出都是 $J$ 上 i.i.d. 均匀样本。

- 未切换：`AllPairs` 完整枚举 $J$，再数组均匀采样 → 正确。
- 已切换：丢弃枚举结果不影响 Phase1 得到的精确 $w_e$ 与 $W$，后续完全等价于版本 B 的 Phase2+Phase3 → 正确。
- “切换不破坏正确性”的原因：切换只改变是否保存 `Report` 的输出，不改变事件顺序、Insert/Delete 时刻与 `CountIntersect` 的精确权重。

------

## 5.5 复杂度分析（R-tree 参数化 + 最坏退化提醒）

R-tree 的时间高度依赖数据分布与 MBR 重叠程度，因此文档采用“访问节点数/候选数”刻画一次查询代价，并明确指出最坏可能退化到二次。

### 5.5.1 预处理公共开销

- 构造并排序事件：$O(n\log n)$
- 两棵 R-tree 初始化：$O(1)$（空树）
- 扫描过程中每个事件一次 Insert/Delete（常见实现平均 $O(\log_B n)$，最坏可能 $O(n)$）

------

### 5.5.2 版本 A：Enumerate+Sampling 的时间/空间

**时间（参数化表达）：**

- 排序：$O(n\log n)$

- 扫描：

  - 更新：$\sum \text{UpdateCost}$

  - 每个 START 一次报告：令该次访问节点数为 $V_e$，候选数为 $C_e$，输出真命中数为 $k_e$（且 $\sum_e k_e=|J|$），则总计
    $$
    O\Big(\sum_e (V_e+C_e)\Big)+O(|J|).
    $$

- 数组采样：$O(t)$

综合可写为：
$$
T_{\text{Enum}}
=
O(n\log n)
+
O\Big(\sum_e (V_e+C_e)\Big)
+
O(|J|+t),
$$
最坏可能退化到 $\Theta(n^2)$ 级别。

**空间：**

- 两棵 R-tree + 事件：$O(n)$
- `Pairs`：$\Theta(|J|)$
- 输出：$O(t)$

$$
S_{\text{Enum}}=O(n+|J|+t).
$$



------

### 5.5.3 版本 B：Sampling 的时间/空间

引入记号：对查询盒 $q^\star$：

- $V(q^\star)$：访问节点数
- $E(q^\star)$：检查的叶条目数
- 树高 $h\approx O(\log_B n)$

则：

- `CountIntersect(q^\star)`：$O(V(q^\star)+E(q^\star))$

- `SampleIntersect(q^\star,t′)`：先构建 guide/frontier，再采样：
  $$
  O(V(q^\star)+E(q^\star)+t′\cdot h).
  $$



**总时间：**

- Phase1：
  $$
  T_1 = O(n\log n) + O(\text{updates}) + \sum_{e\in E} O(V(q_e^\star)+E(q_e^\star)).
  $$

- Phase2：$T_2=O(n+t)$

- Phase3：
  $$
  T_3 = O(\text{updates})+\sum_{e:t_e>0} O(V(q_e^\star)+E(q_e^\star)+t_e\cdot h).
  $$

且 $\sum_e t_e=t$，因此采样下降项汇总为 $O(t\cdot h)$。

**空间：**

- 两棵 R-tree：$O(n)$
- 事件权重 + alias：$O(n)$
- slot + 输出：$O(t)$

$$
S_{\text{Sampling}}=O(n+t).
$$



**最坏情况警告：** 若大量 MBR 重叠导致剪枝失败，可能出现 $V(q)=\Theta(n)$，从而总体退化到 $\Theta(n^2)$ 级。

------

### 5.5.4 版本 C：Adaptive+Sampling 的时间/空间

分两种分支：

#### Case A：未切换（$|J|\le J_\star$）

等价于版本 A，但保证 $|J|\le J_\star$：
$$
T_{\text{Adaptive,base}}
=
O(n\log n)+O(\text{sweep queries})+O(|J|)+O(t),
\quad
S_{\text{Adaptive,base}}=O(n+|J|+t).
$$


#### Case B：发生切换（$|J|>J_\star$）

- Phase1：始终 Count；Report 总输出量 $\le J_\star$，并会丢弃
- 后续执行版本 B 的 Phase2+Phase3

因此可写为：
$$
T_{\text{Adaptive,switch}}
\approx
T_{\text{Sampling}}
+
O(J_\star),
$$
