# 1. 问题定义与分析

## 1.1 输入与相交对集合

在二维平面 \(\mathbb{R}^2\) 中给定两类轴对齐半开矩形集合：

\[
R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\},\quad
n=n_1+n_2.
\]

每个矩形（半开语义）：
\[
r=[L_x(r),R_x(r))\times[L_y(r),R_y(r)),\quad L_x(r)<R_x(r),\ L_y(r)<R_y(r).
\]

我们只关心跨集合相交对（按固定顺序输出 \((r_c,r_{\bar c})\)）：
\[
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
\]

半开区间下，一维区间相交判定为：
\[
[a,b)\cap[c,d)\neq\varnothing\iff \max(a,c)<\min(b,d).
\]

---

## 1.2 采样目标

输出 \(t\) 个样本 \(Z_1,\dots,Z_t\in J\)，要求：

- **均匀性**：\(\Pr(Z_j=P)=1/|J|\)，对任意 \(P\in J\)
- **独立性**：\(Z_1,\dots,Z_t\) 相互独立（i.i.d.，有放回）

目标是在尽量不显式枚举 \(J\)（最坏可达 \(\Theta(n^2)\)）的前提下实现上述分布。

---

## 1.3 Plane Sweep（按 x）与事件分块：\(J=\biguplus_e J_e\)

取 x 轴做 plane sweep。每个矩形 \(r\) 产生两个事件：

- \(\texttt{START}(r)\) at \(x=L_x(r)\)
- \(\texttt{END}(r)\) at \(x=R_x(r)\)

**事件稳定排序（非常关键）**：

1. 按 \(x\) 升序；
2. 同一 \(x\)：**END 在 START 前**（匹配半开边界，避免“贴边误相交”）；
3. 同一 \(x\) 且都是 START：固定集合优先级（例如先处理 \(R_c\) 再处理 \(R_{\bar c}\)），再用 id 打破平局——确保每个相交对“唯一归属”。

扫描到位置 \(x_0\) 时，活跃集：
\[
\mathrm{Active}(R_\star)=\{r\in R_\star\mid x_0\in[L_x(r),R_x(r))\}.
\]

对每个 START 事件 \(e\)（新加入矩形 \(q=r(e)\)），定义“对立活跃相交对象集合”：
\[
K_e=\{r\in\mathrm{Active}(\text{Other})\mid q\cap r\neq\varnothing\},\qquad w_e=|K_e|.
\]

并把局部相交对集规范化为输出顺序 \((r_c,r_{\bar c})\)：

- 若 \(q\in R_c\)，则 \(J_e=\{(q,r_{\bar c})\mid r_{\bar c}\in K_e\}\)
- 若 \(q\in R_{\bar c}\)，则 \(J_e=\{(r_c,q)\mid r_c\in K_e\}\)

**关键事实（唯一归属 / 不交分解）**：
\[
J=\biguplus_{e\in E}J_e,\qquad |J|=\sum_{e\in E} w_e,
\]
其中 \(E\) 为所有 START 事件集合。

直觉：任意相交对 \((r_c,r_{\bar c})\) 会在两者 START 事件中“更靠后”的那个 START 时刻第一次同时活跃，因此只会被记入一次。

---

## 1.4 START 时刻：x 维已保证重叠，只需处理 y 维区间相交

处理 \(\texttt{START}(q)\) 时，扫描位置 \(x_0=L_x(q)\)。对任意 \(r\in\mathrm{Active}(\text{Other})\)，由活跃定义可知 \(x_0\in[L_x(r),R_x(r))\)，同时 \(x_0\in[L_x(q),R_x(q))\) 恒成立，故 **x 维必重叠**。

于是二维相交等价于 y 区间相交：
\[
q\cap r\neq\varnothing\iff [L_y(q),R_y(q))\cap[L_y(r),R_y(r))\neq\varnothing.
\]

因此在 Interval Tree 版本中，所有几何子问题都退化为：

> 动态维护“活跃 y 区间集合”，对查询区间 \([a,b)\) 进行 **Count / Report / Sample**。

---

## 1.5 A/B 分解：把“区间-区间相交”拆成两个易实现子查询

令查询区间 \(Q=[a,b)\)，任一活跃区间 \(I=[\ell,r)\)。则
\[
I\cap Q\neq\varnothing
\iff
(\ell\le a<r)\ \lor\ (a<\ell<b).
\]

定义：

- **A 型（stabbing / 含点）**：\(\ell\le a<r\)（区间包含点 \(a\)）
- **B 型（start-in-range / 起点落入）**：\(a<\ell<b\)（区间起点落在开区间 \((a,b)\)）

二者互斥（\(\ell=a\) 只属于 A），且覆盖所有相交情形。

记
\[
A(Q)=\{I\in\mathcal{I}_{\text{act}}: \ell\le a<r\},\qquad
B(Q)=\{I\in\mathcal{I}_{\text{act}}: a<\ell<b\}.
\]
则
\[
\mathrm{Overlap}(Q)=A(Q)\uplus B(Q),\qquad
|\mathrm{Overlap}(Q)|=|A(Q)|+|B(Q)|.
\]

这一步的意义是：

- A 型是“点刺查询”（stabbing），很适合 segment tree/interval tree；
- B 型是“起点范围查询”（range on left endpoints），也可用 segment tree + 桶实现；
- 两者都能自然实现 **计数** 与 **均匀采样**，并且可以批量采样到 \(O(\log n + t')\)。

---

# 2. 核心数据结构

Interval Tree 统一版本建议把几何部分封装为一个对象（对每个类各维护一份），并提供统一接口。

## 2.1 事件系统（Plane Sweep）

- 事件数组 `Events[1..2n]`：每条记录
  \[
  (x,\ \text{type}\in\{\texttt{START},\texttt{END}\},\ \text{class}\in\{c,\bar c\},\ \text{id}).
  \]
- 稳定排序：`x` 升序；同 `x`：END 在 START 前；同为 START：按 class 优先级再按 id。

扫描过程中，活跃集不显式存成列表，而是由下述 IntervalTree 结构“隐式表示”。

---

## 2.2 IntervalTreeSampler：维护活跃 y 区间的 Count / Report / Sample

对每一类 \(R_c\) 与 \(R_{\bar c}\)，维护一个结构 \(\mathcal{IT}^{(c)}\)、\(\mathcal{IT}^{(\bar c)}\)，内部由两部分组成：

- **IT-A（stabbing）**：支持 \(A\) 型（包含点 \(a\)）
- **IT-B（start-in-range）**：支持 \(B\) 型（起点落入 \((a,b)\)）

### 2.2.1 对外接口（建议固定下来）

令查询区间 \(I(q)=[a,b)\)：

- `Insert(r)` / `Delete(r)`：插入/删除活跃区间 \(I(r)=[L_y(r),R_y(r))\)
- `CountA(a)`：返回 \(|A(I(q))|\)
- `CountB(a,b)`：返回 \(|B(I(q))|\)
- `CountOverlap(a,b)`：返回 \(|A|+|B|\)
- `ReportA(a)` / `ReportB(a,b)` / `ReportOverlap(a,b)`：枚举对应集合（只在枚举分支使用）
- `SampleA(a,t')` / `SampleB(a,b,t')`：在对应集合上 **i.i.d. 均匀有放回**抽取 \(t'\) 个对象
- `SampleOverlap(a,b,t')`：可先按 \(|A|,|B|\) 将 \(t'\) 个样本位置分配到 A/B，再分别调用 SampleA / SampleB

> 说明：Sampling/Adaptive 的两遍采样阶段需要 `Count` 与 `Sample`；Enumerate+Sampling 与 Adaptive 的小分支需要 `Report`。

---

### 2.2.2 IT-A：点刺（stabbing）结构

**目标**：动态维护活跃区间集合 \(\mathcal{I}_{\text{act}}\)，对点 \(a\) 支持

- `CountA(a) = |{I=[\ell,r)\in\mathcal{I}_{\text{act}}: \ell\le a<r}|`
- `ReportA(a)`：枚举这些区间
- `SampleA(a,t')`：从这些区间 i.i.d. 均匀抽样 \(t'\) 次

**一种推荐实现（可复现、且支持批量采样）**：

1. 收集所有 y 端点（所有 \(L_y\) 与 \(R_y\)），排序去重：\(y_1<y_2<\cdots<y_m\)。
2. 以叶子单元区间 \([y_j,y_{j+1})\) 为叶构造平衡段树 \(\mathcal{T}\)。
3. 对任意区间 \(I\)，在树上求其 **最小规范覆盖**（canonical cover）\(\mathcal{C}(I)\)，大小 \(O(\log n)\)。
4. 每个结点 \(v\) 维护一个“活跃桶” `BucketA(v)`（致密数组 + 句柄表，支持 O(1) 均匀抽取 + O(1) 删除）。
5. `Insert/Delete(I)`：将 \(I\) 加入/移出 \(\mathcal{C}(I)\) 覆盖到的各桶。

**关键性质（路径唯一归属）**：对任意点 \(a\)，所有包含 \(a\) 的区间在根到叶路径上的桶形成**不交并**（每个区间恰好落在路径上的一个桶里），因此：

- `CountA(a)`：访问根到叶路径 \(O(\log n)\) 个桶求和
- `ReportA(a)`：枚举这些桶中所有元素，总成本 \(O(\log n + k)\)
- `SampleA(a,t')`：
  - 先取路径上的桶集合 \(\{v_1,\dots,v_s\}\)，令权重 \(w_i=|\text{BucketA}(v_i)|\)
  - 建一次 alias（或前缀和），然后做 \(t'\) 次：按权重选桶，再桶内均匀抽一个元素
  - 总成本 \(O(\log n + t')\)

---

### 2.2.3 IT-B：起点范围（start-in-range）结构

**目标**：动态维护活跃区间的左端点 \(\ell=L_y(r)\) 作为“点”，对开区间 \((a,b)\) 支持

- `CountB(a,b) = |{I\in\mathcal{I}_{\text{act}}: a<\ell<b}|`
- `ReportB(a,b)`：枚举这些区间
- `SampleB(a,b,t')`：从这些区间 i.i.d. 均匀抽样 \(t'\) 次

**一种推荐实现（与上面 IT-A 风格统一，且可批量采样）**：

1. 把所有可能的左端点键 \((L_y(r), \text{id}(r))\) 排序得到秩数组（用 id 解决重复值）。
2. 在秩区间 \([1..N]\) 上建立平衡线段树 \(\mathcal{S}\)。
3. 每个活跃点（即某个活跃区间的左端点）插入到 **根到叶路径上的所有结点桶** `BucketB(v)`。
   - 插入/删除：沿路径更新 \(O(\log n)\) 个桶。
4. 对查询开区间 \((a,b)\)：
   - 用二级键把严格不等号转为秩边界：
     - 左边界用 `upper_bound((a, +∞))`，得到第一个 \(>a\) 的秩；
     - 右边界用 `lower_bound((b, -∞))`，得到第一个 \(\ge b\) 的秩；
     - 得到秩区间 \([L,R)\)，对应 \((a,b)\)。
   - 将 \([L,R)\) 做规范分解，得到 \(O(\log n)\) 个两两不交的结点区间集合 \(\mathcal{C}_{[L,R)}\)。

**关键性质（分解桶不交）**：由于每个点被插入到其根到叶路径所有桶中，而 \([L,R)\) 的规范分解结点集合两两不交且不包含祖先-后代关系，因此每个点（在 \([L,R)\) 内）恰好出现在其中一个结点桶里，从而：

- `CountB(a,b)`：对分解桶求和，\(O(\log n)\)
- `ReportB(a,b)`：枚举分解桶中所有元素，\(O(\log n + k)\)
- `SampleB(a,b,t')`：
  - 对分解桶按桶大小做一次 alias
  - 做 \(t'\) 次：按权重抽桶，桶内均匀抽元素
  - 总成本 \(O(\log n + t')\)

---

## 2.3 Alias 与 slot（事件级 / A-B 级）

Sampling 与 Adaptive 的两遍分支需要离散采样规划：

- **事件级 alias**：在 START 事件集合上按权重 \(w_e\) 建 alias，使得
  \(\Pr(E=e)=w_e/|J|\)，单次 \(O(1)\)
- **A/B 级选择**：固定事件 \(e\) 后，按 \((w_e^{(A)},w_e^{(B)})\) 选类型
- **slot**：记录每个事件（以及 A/B）被分配到的输出位置
  - `S_e^{(A)}`, `S_e^{(B)}` 记录应由该事件该类型生成的样本位置集合

---

# 3. 算法详细流程

本节给出三种版本的流程，并尽量复用同一套符号与数据结构。

为清晰起见，假设我们对每个矩形都分配了唯一 `id`，并且在扫描时总能根据 `id` 取到其 \(x\)-端点与 \(y\)-端点。

---

## 3.1 版本一：Sampling（两遍扫描，不显式枚举）

### 3.1.1 输入输出

- 输入：\(R_c, R_{\bar c}\)，样本数 \(t\)
- 输出：`Ans[1..t]`，每个元素为一个相交对 \((r_c,r_{\bar c})\in J\)

### 3.1.2 总览

- **Phase1（第一次扫描）**：只做计数，得到每个 START 的 \(w_e^{(A)},w_e^{(B)},w_e\) 与全局 \(W=|J|\)
- **Phase2（离散规划）**：事件 alias + A/B 选择，产生 slot：`S_e^A, S_e^B`
- **Phase3（第二次扫描）**：在每个事件 e 的局部集合上批量均匀采样，并按 slot 回填输出

### 3.1.3 Phase1：第一次扫描（只计数）

初始化：

- 构建并排序 `Events`
- 清空两类的 IntervalTreeSampler：\(\mathcal{IT}^{(c)},\mathcal{IT}^{(\bar c)}\leftarrow\emptyset\)
- 为每个 START 事件预留存储：`wA[e], wB[e], w[e]`
- 全局累计 \(W\leftarrow 0\)

扫描 `Events`：

- `END(r)`：
  - 若 \(r\in R_c\)：\(\mathcal{IT}^{(c)}.Delete(r)\)
  - 若 \(r\in R_{\bar c}\)：\(\mathcal{IT}^{(\bar c)}.Delete(r)\)

- `START(q)`（事件 \(e\)，且 \(I(q)=[a,b)\)）：
  - 若 \(q\in R_c\)：
    - `wA[e] = IT^(barc).CountA(a)`
    - `wB[e] = IT^(barc).CountB(a,b)`
  - 若 \(q\in R_{\bar c}\)`：对称从 `IT^(c)` 计数
  - `w[e] = wA[e] + wB[e]`
  - \(W\leftarrow W+w[e]\)
  - 将 \(q\) 插入本侧：\(\mathcal{IT}^{(class(q))}.Insert(q)\)

扫描结束后：\(W=|J|\)。若 \(W=0\)，直接返回空输出。

### 3.1.4 Phase2：采样规划（事件 alias + A/B slot）

1. 在所有 START 事件集合上按权重 \(w[e]/W\) 建 alias。
2. 初始化所有 slot：`S_e^A,S_e^B ← ∅`。
3. 对每个样本位置 \(j=1..t\)：
   - 抽事件：\(E_j\sim w[e]/W\)
   - 条件下抽类型：
     \[
     T_j=\begin{cases}
     A & \text{以概率 } wA[E_j]/w[E_j],\\
     B & \text{以概率 } wB[E_j]/w[E_j].
     \end{cases}
     \]
   - 把 \(j\) 加入 `S_{E_j}^{T_j}`。

### 3.1.5 Phase3：第二次扫描（局部批量采样 + slot 回填）

初始化：

- 清空两类 IntervalTreeSampler 的活跃集（骨架可复用）
- `Ans[1..t]`

再次扫描 `Events`：

- `END(r)`：同 Phase1 删除

- `START(q)`（事件 \(e\)，\(I(q)=[a,b)\)）：
  - 令 \(t_e^A=|S_e^A|\)，\(t_e^B=|S_e^B|\)
  - 若 \(q\in R_c\)：
    - 若 \(t_e^A>0\)：`ListA ← IT^(barc).SampleA(a, t_e^A)`
    - 若 \(t_e^B>0\)：`ListB ← IT^(barc).SampleB(a,b, t_e^B)`
    - 按 slot 回填：
      - 对 `S_e^A = {j_1,...,j_{t_e^A}}`：`Ans[j_u] = (q, ListA[u])`
      - 对 `S_e^B = {j_1,...,j_{t_e^B}}`：`Ans[j_u] = (q, ListB[u])`
  - 若 \(q\in R_{\bar c}\)`：对称从 `IT^(c)` 采样，但输出顺序保持 \((r_c,r_{\bar c})\)：
    - `Ans[j_u] = (List[ u ], q)`
  - 最后插入 \(q\) 到本侧活跃结构

输出 `Ans`。

---

## 3.2 版本二：Enumerate+Sampling（一次扫描显式枚举 + 数组采样）

该版本是经典 baseline：

1. 一次扫描显式枚举所有相交对 \(J\) 存入数组 `Pairs`；
2. 对 `Pairs` 做独立均匀下标采样得到输出。

### 3.2.1 预处理

- 构建并排序 `Events`
- 初始化两类的 IntervalTree（只需支持 `Insert/Delete/ReportOverlap`）
- `Pairs ← []`

### 3.2.2 一次 Plane Sweep：枚举并存储全部 \(J\)

扫描 `Events`：

- `END(r)`：从所属类的 IntervalTree 中删除

- `START(q)`：
  - 若 \(q\in R_c\)：
    - `Cand ← IT^(barc).ReportOverlap(I(q))` （报告所有 y 区间与 \([a,b)\) 相交的对立活跃对象）
    - 对每个 `r` in `Cand`：
      - 二维情形：直接 `Pairs.append((q,r))`
      - 若把 IntervalTree 仅作“第 2 维候选过滤”（\(d>2\)）：再检查剩余维度相交后再 append
    - `IT^(c).Insert(q)`
  - 若 \(q\in R_{\bar c}\)`：对称，回填保持输出顺序：`Pairs.append((r,q))`

扫描结束得到 \(W=|Pairs|=|J|\)。

### 3.2.3 在 `Pairs` 上做 \(t\) 次 i.i.d. 均匀采样

- 若 \(W=0\)：返回空。
- 否则对 \(j=1..t\)：
  - `idx ~ UniformInt(0, W-1)`
  - `Ans[j] = Pairs[idx]`

输出 `Ans`。

---

## 3.3 版本三：Adaptive+Sampling（阈值自适应：小 J 枚举，大 J 两遍采样）

该版本把上面两种方案融合：

- 若 \(|J|\) 很小：走 Enumerate+Sampling，一遍扫描即可完成（常数因子好）；
- 若 \(|J|\) 很大：走 Sampling 的两遍方案（不显式存 \(J\)，可扩展）。

关键在于一个阈值 \(J_\star\)。

### 3.3.1 Phase1：一次扫描（永远计数，必要时枚举）

输入：\(R_c,R_{\bar c}\)、样本数 \(t\)、阈值 \(J_\star\)。

初始化：

- `mode = ENUMERATE`
- `AllPairs = []`（最多存到 \(J_\star\)）
- 清空 \(\mathcal{IT}^{(c)},\mathcal{IT}^{(\bar c)}\)
- \(W\leftarrow 0\)
- 仍为每个 START 事件记录 \(w_e^{(A)},w_e^{(B)},w_e\)

扫描 `Events`：

- `END(r)`：从所属类 \(\mathcal{IT}\) 删除

- `START(q)`（事件 e，\(I(q)=[a,b)\)）：

  1) **先精确计数（无论是否枚举都做）**

  - `wA[e] = IT^(Other).CountA(a)`
  - `wB[e] = IT^(Other).CountB(a,b)`
  - `w[e] = wA[e] + wB[e]`
  - \(W\leftarrow W+w[e]\)

  2) **阈值判断 + 可选枚举（只在 mode==ENUMERATE 时）**

  - 若 `mode==ENUMERATE` 且 \(W\le J_\star\)：
    - `ListA ← IT^(Other).ReportA(a)`
    - `ListB ← IT^(Other).ReportB(a,b)`
    - 对 `ListA ∪ ListB` 中每个对立对象 `r`：
      - 若 \(q\in R_c\)：`AllPairs.append((q,r))`
      - 若 \(q\in R_{\bar c}\)`：`AllPairs.append((r,q))`

  - 若 `mode==ENUMERATE` 且 \(W>J_\star\)：触发切换
    - `mode ← COUNT_ONLY`
    - 释放/丢弃 `AllPairs`（不再占用内存）
    - 从此刻起后续 START 不再 `Report`，只 `Count`

  > 注意：必须“先 Count 再判断”，这样不会出现“一个 START 事件只枚举一半就切换”的复杂情况，同时保证 Phase1 枚举量严格 \(\le J_\star\)。

  3) 插入 \(q\) 到本侧 \(\mathcal{IT}\)。

Phase1 结束：得到精确 \(W=|J|\)。若 \(W=0\)，返回空。

### 3.3.2 分支 A：未切换（\(|J|\le J_\star\)）→ 直接数组采样

若最终 `mode==ENUMERATE`，则 `AllPairs` 已完整枚举 \(J\)。

对 `AllPairs` 做 \(t\) 次独立均匀下标采样，返回 `Ans`。

### 3.3.3 分支 B：已切换（\(|J|>J_\star\)）→ 走 Sampling 的 Phase2 + Phase3

若最终 `mode==COUNT_ONLY`，则 Phase1 已经为每个 START 事件保存了 \(wA,wB,w\) 以及 \(W\)。

- **Phase2**：事件 alias + A/B slot（同 Sampling 版本）
- **Phase3**：第二次扫描 + 局部批量采样 + slot 回填（同 Sampling 版本）

输出 `Ans`。

---

# 4. Adaptive 版本阈值的选择策略

阈值 \(J_\star\) 的作用是：

- 在 \(|J|\) 小时尽量走 baseline（一次扫描、低常数）；
- 在 \(|J|\) 大时尽快停止枚举，转入两遍采样（可扩展、避免爆内存）。

推荐从三方面确定 \(J_\star\)：

## 4.1 内存硬约束（必须满足）

若每个 pair 只存两个 32/64-bit id，记 `sizeof(Pair)` 字节。给 `AllPairs` 分配内存预算 `MemBudget·ρ`（\(0<\rho<1\)）：
\[
J_\star^{\text{mem}}=
\left\lfloor\frac{\rho\cdot\text{MemBudget}}{\text{sizeof(Pair)}}\right\rfloor,
\qquad
J_\star\le J_\star^{\text{mem}}.
\]

这是防止 OOM 的硬上界。

## 4.2 时间建议（保证“大分支”不被枚举拖垮）

在 Interval Tree 实现中：

- 两遍采样（Sampling）主开销约 \(O(n\log n + t)\)
- Adaptive 的大分支在 Phase1 最多额外枚举 \(J_\star\) 个 pair（然后丢弃），因此
\[
T_{\text{big}}=O(n\log n + J_\star + t).
\]

若希望大分支仍与 Sampling 同阶 \(O(n\log n+t)\)，可选
\[
J_\star=O(n\log n+t).
\]
工程上常用：
\[
J_\star^{\text{time}} = C_1\,n\log n + C_2\,t,
\qquad
J_\star=\min(J_\star^{\text{mem}},J_\star^{\text{time}}).
\]

## 4.3 交叉点标定（benchmark 拟合）

如果你有真实数据分布（或生成数据），可以标定“baseline 与两遍采样的交叉点”：

1. 在多组 \((n,t)\) 与代表性数据分布下，分别跑：
   - 纯 baseline（一直枚举）
   - 纯两遍采样（从一开始就只计数、无枚举）
2. 找到耗时相当的输出规模 \(|J_{\text{cross}}|\)
3. 拟合 \(|J_{\text{cross}}|\approx C_1 n\log n + C_2 t\)
4. 运行时取略小于交叉点的 \(J_\star\)，并再受 \(J_\star^{\text{mem}}\) 截断

---

# 5. 算法的分析（正确性、复杂度等）

本节对三个版本分别给出：

- 关键引理
- 正确性（均匀性 + 独立性）
- 复杂度（时间、空间）

---

## 5.1 关键引理

### 引理 1（事件分块不重不漏）

在半开语义 + END-before-START + START tie-break 的事件全序下：
\[
J=\biguplus_{e\in E} J_e,\qquad |J|=\sum_{e\in E} w_e.
\]

理由：任意相交对会在两者 START 中“更靠后”的那个 START 事件时刻第一次同时活跃，因此只会被记入一次。

### 引理 2（A/B 分解互斥且完备）

对查询区间 \(Q=[a,b)\) 与任意活跃区间 \(I=[\ell,r)\)：
\[
I\cap Q\neq\varnothing\iff (\ell\le a<r)\lor(a<\ell<b).
\]
因此 \(\mathrm{Overlap}(Q)=A(Q)\uplus B(Q)\)，且 \(|\mathrm{Overlap}(Q)|=|A|+|B|\)。

### 引理 3（局部混合采样是均匀的）

固定事件 \(e\) 与查询 \(Q_e=[a,b)\)，令 \(w_e^A=|A(Q_e)|\)，\(w_e^B=|B(Q_e)|\)，\(w_e=w_e^A+w_e^B\)。

若一次局部采样按以下方式生成：

1. 以概率 \(w_e^A/w_e\) 选 A，否则选 B；
2. 在选中的集合内做均匀采样；

则任意 \(r\in\mathrm{Overlap}(Q_e)\) 被选中的概率为 \(1/w_e\)，即在 \(J_e\) 内均匀。

---

## 5.2 版本一：Sampling（两遍扫描）

### 5.2.1 正确性

- Phase1 精确得到每个事件权重 \(w_e\)
- Phase2 以 \(\Pr(E=e)=w_e/|J|\) 抽事件
- Phase3 在给定事件 \(e\) 后，根据引理 3 在 \(J_e\) 内均匀抽样

因此对任意 \(P\in J\)，设其唯一归属事件为 \(e^\star\)，则
\[
\Pr(Z_j=P)=\Pr(E_j=e^\star)\cdot\Pr(P\mid E_j=e^\star)
=\frac{w_{e^\star}}{|J|}\cdot\frac{1}{w_{e^\star}}=\frac{1}{|J|}.
\]

独立性来自：Phase2 中每个位置 \(j\) 的 \((E_j,T_j)\) 独立生成；Phase3 中每个 `SampleA/SampleB` 输出 i.i.d.；slot 回填不引入耦合。

### 5.2.2 时间复杂度

- 事件排序：\(O(n\log n)\)
- Phase1：每事件插删 \(O(\log n)\)，每 START 两次计数 \(O(\log n)\) → \(O(n\log n)\)
- Phase2：alias + 规划 \(O(n+t)\)
- Phase3：第二遍扫描插删 \(O(n\log n)\)，且所有局部采样汇总为 \(O(n\log n + t)\)

总时间：
\[
T=O(n\log n+t).
\]

### 5.2.3 空间复杂度

- 两类 IntervalTreeSampler：每个活跃对象被存到 \(O(\log n)\) 个桶 → \(O(n\log n)\)
- 事件/权重：\(O(n)\)
- slot 与输出：\(O(t)\)

总空间：
\[
S=O(n\log n+t).
\]

---

## 5.3 版本二：Enumerate+Sampling（枚举 + 数组采样）

### 5.3.1 正确性

- 枚举阶段：由引理 1，且每个 START 事件仅枚举对立活跃集中的相交对象，故 `Pairs` 不重不漏且恰为 \(J\)
- 采样阶段：对 `Pairs` 做独立均匀下标采样，显然得到 \(J\) 上 i.i.d. 均匀样本

### 5.3.2 时间复杂度（二维）

- 排序：\(O(n\log n)\)
- 扫描：每 START 的 report 成本 \(O(\log n + k_e)\)，总输出 \(\sum_e k_e=|J|\)
- 数组采样：\(O(t)\)

因此：
\[
T=O(n\log n+|J|+t).
\]

（若把 Interval Tree 仅作候选过滤用于 \(d>2\)，则把 \(|J|\) 换成候选总量 \(\sum_e \mathrm{cand}_e\)，最坏仍可达 \(\Theta(n^2)\)。）

### 5.3.3 空间复杂度

- 活跃索引结构：典型 \(O(n)\) 或 \(O(n\log n)\)（取决于实现）
- 显式存储 `Pairs`：\(\Theta(|J|)\)
- 输出：\(O(t)\)

总空间：
\[
S=O(n+|J|+t)\quad(\text{或 }O(n\log n+|J|+t)).
\]

---

## 5.4 版本三：Adaptive+Sampling

### 5.4.1 正确性

Adaptive 的正确性来自“两条分支都正确 + 切换不改变计数语义”。

- 若未切换（\(|J|\le J_\star\)）：Phase1 完整枚举得到 `AllPairs=J`，数组采样正确
- 若已切换（\(|J|>J_\star\)）：Phase1 始终通过 `CountA/CountB` 得到精确 \(w_e\)，丢弃枚举结果不影响权重；后续 Phase2+Phase3 与 Sampling 完全一致，因此正确

切换只影响“是否保留枚举输出、是否继续调用 Report”，不影响事件顺序、Insert/Delete 时刻、以及每个 START 的精确计数 \(w_e\)。

### 5.4.2 时间复杂度

设阈值为 \(J_\star\)。

- 未切换：
  \[
  T_{\text{base}}=O(n\log n+|J|+t).
  \]
- 已切换：
  - Phase1：计数 \(O(n\log n)\) + 最多枚举 \(J_\star\)
  - Phase2+Phase3：\(O(n\log n+t)\)

  因此
  \[
  T_{\text{big}}=O(n\log n+J_\star+t).
  \]

若取 \(J_\star=O(n\log n+t)\)，则大分支仍为 \(O(n\log n+t)\)。

### 5.4.3 空间复杂度

- 几何结构与权重：\(O(n\log n)\)
- 未切换：还需 `AllPairs` 的 \(|J|\) 存储
- 已切换：`AllPairs` 峰值最多 \(J_\star\) 且会释放，另有 slot/output \(O(t)\)

因此峰值空间：
\[
S_{\text{base}}=O(n\log n+|J|+t),\qquad
S_{\text{big}}=O(n\log n+\max(J_\star,t)).
\]

---

## 5.5 实现注意事项（容易踩坑但非常关键）

1. **半开语义的一致性**：
   - 相交判定必须用严格不等号 \(\max(L_1,L_2)<\min(R_1,R_2)\)。
   - 事件排序必须 END-before-START，否则会把“贴边”误认为相交。

2. **B 型严格不等号** \(a<\ell<b\)：
   - 强烈建议用二级键 \((y,\text{id})\) 做秩边界（left 用 upper_bound，right 用 lower_bound）。

3. **桶的“致密数组 + 句柄”**：
   - 为支持 O(1) 桶内均匀抽样与 O(1) 删除，建议每个桶维护 vector + position map（swap-delete）。

4. **随机性与独立性**：
   - `SampleA/SampleB` 每次抽样应视为独立调用（或使用独立随机源），确保理论 i.i.d. 与实现一致。

---

## 一句话总结

> **Sampling‑IT**：两遍扫描，Phase1 精确计数得到每个 START 的权重与 A/B 分解，Phase2 事件 alias + A/B slot，Phase3 二次扫描做局部均匀采样并回填，整体 \(O(n\log n+t)\) 且输出 i.i.d. 均匀。
>
> **Enumerate‑IT**：一遍扫描枚举所有相交对形成数组，再做 \(t\) 次独立均匀下标采样，二维下 \(O(n\log n+|J|+t)\)。
>
> **Adaptive‑IT**：用阈值 \(J_\star\) 在小 \(|J|\) 自动走枚举 baseline，大 \(|J|\) 自动切换两遍采样，保持分布严格正确并兼顾常数与可扩展性。