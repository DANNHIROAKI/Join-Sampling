# 纯 Range Tree Join Sampling 统一算法文档（RT‑Sampling / RT‑Enumerate / Adaptive‑RT）

## 记号与约定

- 维度：$d\ge 2$，并假设 $d$ 为常数（因此 $k=2(d-1)$ 也是常数）。

- 两类半开轴对齐盒子集合：
  $$
  R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
  R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\},\quad n=n_1+n_2.
  $$
  每个盒子
  $$
  r=\prod_{i=1}^d [L_i(r),R_i(r)),\qquad L_i(r)<R_i(r).
  $$

- 输出 pair 统一写成 $(r_c,r_{\bar c})$。当 START 来自 $R_{\bar c}$ 时需交换顺序写入输出。

- 半开区间相交判定：
  $$
  [a,b)\cap[c,d)\neq\varnothing \iff \max(a,c)<\min(b,d).
  $$

------

## 1. 问题定义与分析

### 1.1 Join 结果集合

我们只关心跨集合相交对：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
目标是从 $J$ 中输出 $t$ 个样本 $Z_1,\dots,Z_t\in J$，要求 **i.i.d. 均匀有放回**：
$$
\Pr(Z_j=P)=\frac{1}{|J|}\ (\forall P\in J),\qquad Z_1,\dots,Z_t\ \text{相互独立}.
$$

### 1.2 第 1 维 plane sweep 与事件系统（保证“唯一归属”）

选择第 1 维 $x_1$ 做扫描。每个盒子 $r$ 产生两事件：

- $\texttt{START}(r)$ at $x_1=L_1(r)$
- $\texttt{END}(r)$ at $x_1=R_1(r)$

稳定排序规则（必须与半开边界匹配）：

1. 按 $x_1$ 升序；
2. 同 $x_1$：**END 在 START 前**（避免贴边误相交）；
3. 同 $x_1$ 且同为 START：固定集合优先级（如先 $R_c$ 后 $R_{\bar c}$），再用 id 打破平局，以便证明“唯一归属”。

扫描位置为 $x_0$ 时定义活跃集：
$$
\mathrm{Active}(R_c)=\{r_c\in R_c:\ x_0\in[L_1(r_c),R_1(r_c))\},\quad
\mathrm{Active}(R_{\bar c})\ \text{同理}.
$$
**关键观察（START 时第 1 维必然重叠）**：处理 $\texttt{START}(q)$ 时 $x_0=L_1(q)$，对任何对立活跃盒子 $r\in\mathrm{Active}(\text{Other})$ 有 $x_0\in[L_1(r),R_1(r))$。由于“END 在 START 前”且半开，因此不会把仅贴边的盒子误认为活跃。于是：在 START 时，只要能判断维度 $2..d$ 相交，就能推出整盒相交。

### 1.3 START 事件诱导的局部块分解：$J=\biguplus_e J_e$

令 $E$ 为所有 START 事件集合。对一个 START 事件 $e$，设新盒子 $q=r(e)$，对立活跃集为 $\mathrm{Active}(\text{Other})$。定义：
$$
K_e=\{r\in \mathrm{Active}(\text{Other})\mid q\cap r\neq\varnothing\},\qquad w_e:=|K_e|.
$$
并定义局部相交对集合（输出统一为 $(r_c,r_{\bar c})$）：

- 若 $q\in R_c$：$J_e=\{(q,r_{\bar c})\mid r_{\bar c}\in K_e\}$
- 若 $q\in R_{\bar c}$：$J_e=\{(r_c,q)\mid r_c\in K_e\}$

**命题（事件分块不交分解）**：
$$
J=\biguplus_{e\in E} J_e,\qquad |J|=\sum_{e\in E}w_e.
$$
证明要点：任意相交对的两个 START 中必有“更靠后”的那个（同坐标由 tie‑break 决定）；该时刻另一盒子必在对立活跃集中，因此该对恰在该 START 被发现一次且仅一次。

> 这条分解是 Sampling/Adaptive 版本“按事件权重抽样”的核心。

------

## 2. 核心数据结构

RT 版本的本质是：**用一个高维正交范围查询（维数 $k=2(d-1)$）一次性表达“维度 $2..d$ 相交”**，然后用动态 $k$-维 Range Tree 在活跃集上支持 `Count/Report/Sample`。

### 2.1 从盒子相交到 $k=2(d-1)$ 维正交范围查询

令 $m=d-1$，$k=2m=2(d-1)$。对任意盒子 $r$，取其在维度 $2..d$ 的端点向量：
$$
\mathbf{L}(r)=(L_2(r),\dots,L_d(r))\in\mathbb{R}^m,\quad
\mathbf{R}(r)=(R_2(r),\dots,R_d(r))\in\mathbb{R}^m.
$$
定义嵌入点：
$$
p(r)=(\mathbf{L}(r),\mathbf{R}(r))=
(L_2(r),\dots,L_d(r),R_2(r),\dots,R_d(r))\in\mathbb{R}^{2m}.
$$
对查询盒子 $q$，维度 $i\ge 2$ 的相交条件等价于：
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
**引理（等价性）**：
$$
r\text{ 与 }q\text{ 在维度 }2..d\text{ 相交}\iff p(r)\in \mathcal{Q}(q).
$$

> 注意：在 plane sweep 的 START 时刻，活跃性保证第 1 维也已重叠，所以“维度 $2..d$ 相交”即可推出“整盒相交”。

### 2.2 严格不等号的可实现化：二级键 + rank 边界

因为 $\mathcal{Q}(q)$ 里存在严格不等号（如 $L_i(r)、$R_i(r)>L_i(q)$），建议用“二级键坐标压缩”把它变成稳定、无歧义的 rank 区间：

- 给每个盒子唯一 id：$\mathrm{id}(r)\in\{1,\dots,n\}$；

- 每个坐标值用二级键 $(\text{value},\mathrm{id})$ 表示并做字典序排序压缩为 rank；

- 例如 $L_i(r) 可实现为
  $$
  (L_i(r),\mathrm{id}(r)) < (R_i(q),0),
  $$
  对应上界 `lower_bound((R_i(q),0))`；

- $R_i(r)>L_i(q)$ 可实现为
  $$
  (R_i(r),\mathrm{id}(r)) > (L_i(q),n),
  $$
  对应下界 `upper_bound((L_i(q),n))`。

这样 $\mathcal{Q}(q)$ 在 rank 空间中可当作“闭区间笛卡尔积”喂给 Range Tree。

### 2.3 两棵动态 $k$-维 Range Tree：维护对立活跃集

维护两棵动态 $k$-维 Range Tree：

- $\mathcal{RT}^{(c)}$：维护当前 $\mathrm{Active}(R_c)$ 的点集 $\{p(r_c)\}$
- $\mathcal{RT}^{(\bar c)}$：维护当前 $\mathrm{Active}(R_{\bar c})$ 的点集 $\{p(r_{\bar c})\}$

每个点携带 `box-id` 以便从点回溯盒子对象。

#### 2.3.1 需要的接口

对任一 $\mathcal{RT}$ 支持：

- `Insert(r)` / `Delete(r)`
- `Count(q)`：返回 $\#\{r:\ p(r)\in\mathcal{Q}(q)\}$
- `Report(q)`：枚举所有 $r$ 使 $p(r)\in\mathcal{Q}(q)$
- `Sample(q,t')`：从集合 $\{r:\ p(r)\in\mathcal{Q}(q)\}$ 中 **i.i.d. 均匀有放回**采样 $t'$ 次

对应的标准复杂度假设（常数维 $k$）：
$$
\texttt{Insert/Delete/Count}=O((\log n)^k),\quad
\texttt{Report}=O((\log n)^k+k_{\text{out}}),\quad
\text{Space}=O(n(\log n)^{k-1}).
$$

### 2.4 如何在 Range Tree 上实现 `Sample(q,t')`（核心补强）

Range Tree 原生只提供范围计数/报告。为了实现“均匀采样”，采用**规范分解（canonical decomposition）+ 子块加权混合**，并用 **slot/批量分配** 达到 $O((\log n)^k+t')$ 而不是 $t'(\log n)^k$。这部分是 Sampling/Adaptive 版本正确性与复杂度的关键。

#### 2.4.1 单次采样：按子块大小加权

对一个 $k$-维正交范围 $\mathcal{Q}$，Range Tree 在最外层坐标会生成一个规范分解结点集合
$$
\mathcal{C}=\{v_1,\dots,v_s\},\quad s=O(\log n),
$$
它们对应的子树点集两两不交，且覆盖所有满足外层坐标约束的点。对每个 $v\in\mathcal{C}$，其关联结构是一个 $(k-1)$ 维 Range Tree $\mathcal{RT}_v$。令 $\mathcal{Q}_{\text{tail}}$ 是去掉外层坐标后的剩余范围约束，定义：
$$
w_v := \texttt{Count}_{\mathcal{RT}_v}(\mathcal{Q}_{\text{tail}}),\qquad W:=\sum_{v\in\mathcal{C}}w_v.
$$
**单次 `Sample`**：

1. 以概率 $w_v/W$ 选择一个结点 $v\in\mathcal{C}$；
2. 递归在 $\mathcal{RT}_v$ 上对 $\mathcal{Q}_{\text{tail}}$ 做均匀采样，得到一个点。

因为每个候选点属于且只属于一个子块，所以全局输出概率为 $1/W$，从而单次均匀成立。

#### 2.4.2 批量采样：slot 聚合，复杂度 $O((\log n)^k+t')$

若直接重复单次过程 $t'$ 次，会得到 $t'(\log n)^k$。为了得到 $+t'$，做法是：

1. 先计算所有 $w_v$，并在 $\mathcal{C}$ 上按权重建 alias（或前缀和）；
2. 独立抽 $t'$ 次结点并统计每个结点被抽到次数 $t_v$（slot 数）；
3. 对每个 $v$ 只递归调用一次：

$$
\mathcal{RT}_v.\texttt{Sample}(\mathcal{Q}_{\text{tail}},t_v),
$$

把返回的 $t_v$ 个点写回对应 slot。

这样总成本为：
$$
O((\log n)^k)+\sum_{v\in\mathcal{C}}O(t_v)=O((\log n)^k+t').
$$
slot 只是把“独立试验按结点分桶再批量执行”的实现技巧，因此输出仍是 i.i.d. 均匀有放回。

------

## 3. 算法详细流程（三个版本）

本章给出三个版本：**RT‑Sampling**、**RT‑Enumerate**、**Adaptive‑RT**。它们共享同一套预处理与事件系统。

### 3.0 公共预处理

输入：$R_c,R_{\bar c}$，样本数 $t$（Adaptive 还要 $J_\star$）。

1. 构造事件数组 `Events[1..2n]`：每盒产生 START/END；
2. 按 §1.2 的稳定规则排序；
3. 对维度 $2..d$ 的端点做二级键坐标压缩（§2.2），以便实现严格不等号；
4. 初始化两棵空 Range Tree：$\mathcal{RT}^{(c)},\mathcal{RT}^{(\bar c)}$。

------

### 3.1 版本一：RT‑Enumerate（Enumerat+Sampling）

**思路**：一次 sweep 显式枚举全部 $J$ 到数组 `Pairs`，再对数组做 $t$ 次独立均匀下标采样。 

#### 3.1.1 枚举阶段：一次 plane sweep + Report

初始化 `Pairs=[]`。

依序扫描 `Events`：

- 若 `END(r)`：
  - 若 $r\in R_c$：$\mathcal{RT}^{(c)}.\texttt{Delete}(r)$
  - 若 $r\in R_{\bar c}$：$\mathcal{RT}^{(\bar c)}.\texttt{Delete}(r)$
- 若 `START(q)`：
  - 若 $q\in R_c$：
    1. `List ← RT^(barc).Report(q)`（即所有 $r_{\bar c}$ 使 $p(r_{\bar c})\in\mathcal{Q}(q)$）
    2. 对每个 $r_{\bar c}\in$List：`Pairs.append((q,r_{\bar c}))`
    3. $\mathcal{RT}^{(c)}.\texttt{Insert}(q)$
  - 若 $q\in R_{\bar c}$（对称）：
    1. `List ← RT^(c).Report(q)`
    2. 对每个 $r_c\in$List：`Pairs.append((r_c,q))`
    3. $\mathcal{RT}^{(\bar c)}.\texttt{Insert}(q)$

扫描结束后 $W=|Pairs|=|J|$。 

#### 3.1.2 采样阶段：数组 i.i.d. 下标采样

- 若 $W=0$：返回空；
- 否则对 $j=1..t$：
   $\text{idx}\sim\mathrm{Unif}\{0,\dots,W-1\}$，输出 `Ans[j]=Pairs[idx]`。



------

### 3.2 版本二：RT‑Sampling（Sampling，不显式枚举 $J$）

**思路**：利用事件分块 $J=\biguplus_e J_e$，先在 Phase1 只算每个 START 的权重 $w_e=|J_e|$ 并得到 $W=|J|$，再在事件层做 alias 抽样并用 slot 汇总，最后 Phase3 第二遍扫描时局部调用 `RT.Sample` 真正生成样本。 

#### Phase 1：第一次扫描（只计数）

初始化 `W=0`，为每个 START 事件准备记录 `w_e`。

扫描 `Events`：

- `END(r)`：按所属集合在对应 $\mathcal{RT}$ 中删除。
- `START(q)`（事件 $e$）：
  - 若 $q\in R_c$：$w_e \leftarrow \mathcal{RT}^{(\bar c)}.\texttt{Count}(q)$
  - 若 $q\in R_{\bar c}$：$w_e \leftarrow \mathcal{RT}^{(c)}.\texttt{Count}(q)$
     更新 `W += w_e`，并把 $q$ 插入本侧 $\mathcal{RT}$。

结束后由 §1.3 命题得 $W=\sum_e w_e=|J|$。若 $W=0$ 直接返回空。 

#### Phase 2：事件级 alias + slot 分配

1. 在 START 事件集合 $E$ 上按权重 $w_e/W$ 构建 alias（可只保留 $w_e>0$ 的事件）。
2. 初始化每个事件 slot：$S_e\gets\emptyset$。
3. 对每个样本位置 $j=1..t$：抽事件 $E_j\sim w_e/W$，将 $j$ 放入 $S_{E_j}$。

记 $t_e:=|S_e|$，则 $\sum_e t_e=t$。 

#### Phase 3：第二次扫描（按 slot 批量采样并回填）

1. 清空两棵 $\mathcal{RT}$（活跃集置空，骨架可复用）。
2. 初始化输出数组 `Ans[1..t]`。
3. 再扫描一遍 `Events`：

- `END(r)`：删除
- `START(q)`（事件 $e$）：
  - 若 $t_e=0$：仅执行 Insert(q)
  - 若 $t_e>0$：
    - 若 $q\in R_c$：`List ← RT^(barc).Sample(q, t_e)` 得到 $t_e$ 个 i.i.d. 均匀的 $r_{\bar c}$，按 slot 回填 `Ans[j_u]=(q, List[u])`
    - 若 $q\in R_{\bar c}$：从 $\mathcal{RT}^{(c)}$ sample，回填 `Ans[j_u]=(List[u], q)`
  - 最后插入 $q$ 到本侧 $\mathcal{RT}$

输出 `Ans`。 

------

### 3.3 版本三：Adaptive‑RT（Adaptive+Sampling）

Adaptive‑RT 目标：当 $|J|$ 很小，用“一遍扫描枚举 + 数组采样”常数更好；当 $|J|$ 很大，用 RT‑Sampling 避免 $|J|$ 依赖。切换由阈值 $J_\star$ 控制。 

#### Phase 1：一次扫描（永远 Count，必要时 Enumerate）

初始化：

- `mode = ENUMERATE`
- `AllPairs = []`（最多允许增长到 $J_\star$）
- `W=0`
- Range Tree 活跃集清空

扫描 `Events`：

- `END(r)`：删除

- `START(q)`（事件 $e$）：

  1. **先计数**（无论是否枚举都要做）：
     $$
     w_e\leftarrow
     \begin{cases}
     \mathcal{RT}^{(\bar c)}.\texttt{Count}(q), & q\in R_c\\
     \mathcal{RT}^{(c)}.\texttt{Count}(q), & q\in R_{\bar c}
     \end{cases}
     $$
     更新 `W += w_e` 并保存 $w_e$。

  2. **阈值判断 + 可选枚举**：

     - 若 `mode==ENUMERATE` 且 `W ≤ J_star`：允许枚举本事件产生的全部 pairs：
        用 `Report` 得到对立活跃盒子并 append 到 `AllPairs`（保持输出顺序 $(r_c,r_{\bar c})$）。
     - 若 `mode==ENUMERATE` 且 `W > J_star`：发生切换：
        `mode ← COUNT_ONLY`，并丢弃/释放 `AllPairs`，后续不再 Report，只 Count。
     - 若 `mode==COUNT_ONLY`：跳过枚举。

  3. `Insert(q)` 到本侧树

扫描结束，得到精确 $|J|=W=\sum_e w_e$。若 $W=0$ 返回空。 

#### 分支 A（未切换，$|J|\le J_\star$）：Baseline 数组采样

此时 `AllPairs` 已完整枚举 $J$，对其做 $t$ 次独立均匀下标采样即可。

#### 分支 B（已切换，$|J|>J_\star$）：执行 RT‑Sampling 的 Phase2+Phase3

- Phase2：在事件层按 $w_e/W$ 建 alias 并分配 slot；
- Phase3：第二遍扫描，按 slot 调用 `RT.Sample` 回填输出。

注意：这条路径等价于“从一开始就不枚举，只计数”的 RT‑Sampling，只是多付出 $O(J_\star)$ 的有限枚举尝试成本。

------

## 4. Adaptive 版本阈值 $J_\star$ 的选择策略

阈值决定“枚举到多大就放弃、转采样”。推荐按 **内存硬约束 + 时间建议** 两层确定。 

### 4.1 内存硬约束（必须满足）

设可用内存预算 `MemBudget`（字节），允许给 `AllPairs` 使用比例 $\rho\in(0,1)$。每个 pair 的存储开销约为 `sizeof(Pair)`（例如两个 32/64-bit id）。则：
$$
J_\star^{\text{mem}}=\left\lfloor\frac{\rho\cdot \text{MemBudget}}{\text{sizeof(Pair)}}\right\rfloor,\qquad
J_\star\le J_\star^{\text{mem}}.
$$


### 4.2 时间权衡（理论建议）

当发生切换时，Adaptive‑RT 的 Phase1 最多额外枚举 $O(J_\star)$ 个 pair（之后不再 Report），因此大分支时间为：
$$
T_{\text{ours-branch}}=O\bigl(n(\log n)^k+J_\star+t\bigr).
$$
若希望大分支仍与 RT‑Sampling 同阶：
$$
O\bigl(n(\log n)^k+t\bigr),
$$
则应取：
$$
J_\star=O\bigl(n(\log n)^k+t\bigr),
\quad k=2(d-1).
$$
工程上可写成：
$$
J_\star^{\text{time}}=C_1\,n(\log n)^k+C_2\,t,\qquad
J_\star=\min(J_\star^{\text{mem}},J_\star^{\text{time}}).
$$


### 4.3 工程标定（推荐做法）

在多组数据/参数上分别运行“纯枚举版（RT‑Enumerate）”与“纯采样版（RT‑Sampling）”，找到耗时交叉点 $|J_{\text{cross}}|$，用
$$
|J_{\text{cross}}|\approx C_1\,n(\log n)^k+C_2\,t
$$
拟合常数，然后运行时用略小于交叉点的值作为 $J_\star$，并受 $J_\star^{\text{mem}}$ 截断。 

------

## 5. 算法分析（正确性与复杂度；三版本均包含）

### 5.1 数据结构层正确性：`Count/Report/Sample` 的语义

在处理 START(q) 时，对立 Range Tree 存储的是对立活跃集。由 §1.2 “START 时第 1 维必重叠”与 §2.1 引理（维度 $2..d$ 相交 $\Leftrightarrow p(r)\in\mathcal{Q}(q)$），可得：

- `Count(q)` 返回的恰是 $w_e=|K_e|=|J_e|$；
- `Report(q)` 枚举的恰是 $K_e$（因此枚举 pairs 恰是 $J_e$）；
- `Sample(q,t')` 按 §2.4 的 canonical+slot 实现后，输出是 $K_e$ 上 **i.i.d. 均匀有放回**的 $t'$ 个样本。

### 5.2 RT‑Enumerate 的正确性

**(1) 不重不漏枚举 $J$**
 任取 $(r_c,r_{\bar c})\in J$。令两者 START 事件中在事件全序里“更靠后”的那个为 $e^\star$，对应 $q=r(e^\star)$。在处理 $e^\star$ 时，另一个盒子必在对立活跃集中，且维度 $2..d$ 相交 $\Rightarrow p(\cdot)\in\mathcal{Q}(q)$，因此会被 `Report(q)` 报告并写入 `Pairs`。它不可能在更早 START 时出现，因为那时 $q$ 尚未插入活跃结构。故每对恰好一次。 

**(2) 数组采样 i.i.d. 均匀**
 `Pairs` 是 $J$ 的一一枚举，对下标独立均匀采样即得到 $J$ 上 i.i.d. 均匀有放回样本。 

### 5.3 RT‑Sampling 的正确性（全局 i.i.d. 均匀）

利用 §1.3 的不交分解 $J=\biguplus_e J_e$。任取 $P\in J$，存在唯一事件 $e^\star$ 使 $P\in J_{e^\star}$，且 $|J_{e^\star}|=w_{e^\star}$。对任意位置 $j$：

- Phase2：$\Pr(E_j=e)=w_e/W$
- Phase3：给定 $E_j=e$，`Sample(q,t_e)` 在 $J_e$ 上均匀，故 $\Pr(P\mid E_j=e^\star)=1/w_{e^\star}$

于是：
$$
\Pr(Z_j=P)=\frac{w_{e^\star}}{W}\cdot\frac{1}{w_{e^\star}}=\frac{1}{W}=\frac{1}{|J|}.
$$
独立性来自：$\{E_j\}$ 独立生成；每个事件内 `Sample` 输出 i.i.d.；slot 回填只是重排映射，因此整体输出 i.i.d.。 

### 5.4 Adaptive‑RT 的正确性（切换不改变分布）

- 若未切换：`AllPairs` 完整枚举 $J$，输出等价于 RT‑Enumerate → 正确；
- 若已切换：丢弃枚举结果，但仍保留所有 $w_e$ 的精确计数（来自 `Count`），随后执行 RT‑Sampling 的 Phase2+Phase3 → 正确；
- 切换仅影响“是否保留 `Report` 的输出”，不改变事件顺序、Insert/Delete 时刻、也不改变 $w_e$ 的计算方式，因此不会引入偏差。

------

### 5.5 复杂度（时间与空间；三版本）

令 $k=2(d-1)$（常数）。

#### 5.5.1 RT‑Enumerate

- **时间**

  - 事件排序与压缩：$O(n\log n)$（相对主项可忽略）
  - 扫描：每事件 Insert/Delete $O((\log n)^k)$；每 START 做一次 `Report`：$O((\log n)^k+k_e)$
  - $\sum_e k_e=|J|$
  - 采样：$O(t)$
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
  - Phase2：alias+slot $O(n+t)$
  - Phase3：更新 $O(n(\log n)^k)$ + 批量采样总计 $O(n(\log n)^k+t)$
     合并：

  $$
  T_{\text{Sampling}}=O\bigl(n(\log n)^k+t\bigr).
  $$

- **空间**
  $$
  S_{\text{Sampling}}=O\bigl(n(\log n)^{k-1}+t\bigr)
  $$
  （外加 $O(n)$ 存事件与权重）。
  

#### 5.5.3 Adaptive‑RT

分两种情况： 

- **Case A：未切换（$|J|\le J_\star$）**
  $$
  T=O\bigl(n(\log n)^k+|J|+t\bigr),\quad
  S=O\bigl(n(\log n)^{k-1}+|J|+t\bigr).
  $$

- **Case B：已切换（$|J|>J_\star$）**
   Phase1 最多枚举 $O(J_\star)$ 后停止，因此：
  $$
  T=O\bigl(n(\log n)^k+J_\star+t\bigr),\quad
  S=O\bigl(n(\log n)^{k-1}+\max(J_\star,t)\bigr).
  $$
  若按 §4.2 取 $J_\star=O(n(\log n)^k+t)$，则大分支渐进与 RT‑Sampling 同阶。