# AGR‑BoxJoin Baseline（升级版 v2）

**网格上界 μ + Alias 加权 + Rejection，输出严格 i.i.d. 均匀有放回的 box‑intersection join samples**

> Baseline 定位：
>  这是一个**理论干净、实现相对简单、非常适合对照**的 i.i.d. uniform join sampling baseline。它的核心是构造一个“易采样的候选超集”并用 rejection 把分布纠正回真实 join 结果。
>  代价是：**性能强依赖上界 μ 的紧致程度（接受率）**，并且网格参数选择不当会显著退化（这也是 baseline 的价值所在：清晰展示“上界松→接受率低→采样慢”的风险曲线）。

------

## 1. 问题定义与分析

### 1.1 输入对象：两组 d 维轴对齐半开盒（MBR）

在维度 $d\ge 2$ 的空间中给定两组盒子集合：
$$
A=\{r_1,\dots,r_{n_1}\},\qquad B=\{s_1,\dots,s_{n_2}\}.
$$
每个盒子为半开超矩形：
$$
x\in r \iff x\in \prod_{i=1}^d [L_i(r),R_i(r)),\quad L_i(r)<R_i(r).
$$
**半开相交判定**（推荐作为唯一标准判定函数）：
$$
\texttt{Intersect}(r,s)\equiv \bigwedge_{i=1}^d\left(\max(L_i(r),L_i(s))<\min(R_i(r),R_i(s))\right).
$$
半开与严格不等号的组合可以消除“贴边算不算相交”的歧义。

------

### 1.2 Join 结果集合与采样目标

我们关心 cross‑set intersection join：
$$
J=\{(r,s)\mid r\in A,\ s\in B,\ \texttt{Intersect}(r,s)=\texttt{true}\}.
$$
给定样本数 $t$，目标输出
$$
Z_1,\dots,Z_t\in J
$$
并满足**严格 i.i.d. 均匀有放回**：
$$
\Pr(Z_j=P)=\frac1{|J|}\ \ (\forall P\in J),\qquad Z_1,\dots,Z_t\ \text{相互独立}.
$$

> 工程边界：若 $|J|=0$，则“在空集上均匀采样”数学上无定义。我们规定：**返回空**并显式标记 join 为空（本报告在算法部分给出“安全终止”的默认做法）。

------

### 1.3 为什么“先均匀选 r，再均匀选其匹配 s”会偏

定义度数：
$$
\deg(r)=|\{s\in B\mid \texttt{Intersect}(r,s)\}|.
$$
若先均匀选 $r$，再从其匹配集合均匀选 $s$，则 pair 的概率 $\propto 1/n_1\cdot 1/\deg(r)$，会**显著偏向小度数对象**，无法做到对所有 pair 均匀。
 因此要获得 uniform join samples，必须使用“按匹配数或其上界”进行加权（weighted sampling）——这也是 Amagata’25 及经典 join sampling 文献反复强调的原因。

------

## 2. 引用文章内容与我们继承/改造的框架

这一节专门回答你要求的“引用文章标题 + 引用内容详细描述”。

### 2.1 被引用论文信息与核心结论

**论文**：Daichi Amagata, *Random Sampling over Spatial Range Joins*（arXiv 版本 2025‑08‑20；文中说明已被 ICDE 2025 接收）。

论文讨论的是 **spatial range join**：给定点集 $R,S$，每个 $r\in R$ 关联一个窗口 $w(r)$，join 结果为 $\{(r,s)\mid s\in w(r)\}$。论文明确指出：为了后续分析/学习任务有效，join samples 需要 **uniform 且 independent**。

并且论文在形式化定义中（Definition 2）将任务描述为：返回 $t$ 个来自 join 结果集合 $J$ 的样本，“each of which is picked uniformly at random”，并且文中强调 uniform 与 independent 的重要性，且默认“无损一般性”假设 $|J|\ge 1$。

------

### 2.2 Amagata’25 的 baseline 框架：Alias 加权 + 上界 μ + Rejection

论文在 baseline 部分设计了两类 KDS‑based baseline，其中第二个（KDS‑rejection）特别关键：它用**易算上界 $\mu(r)\ge |S(w(r))|$** 替代精确计数，并用 rejection 纠正分布；它仍可保证 uniform & independent。

其典型结构是：

1. 为每个 $r$ 计算上界 $\mu(r)$；
2. 用 Walker alias 让 $r$ 以概率 $\mu(r)/\sum_{r'}\mu(r')$ 被采到；
3. 再在对应范围内采样 $s$ 并进行 accept/reject，从而保证每个 join pair 被输出的概率一致。

论文还解释了 grid mapping 如何让 $\mu(r)$ 以常数/低代价得到，并指出：上界若不够紧，会导致接受率低、采样迭代次数增加，从而变慢。

------

### 2.3 我们的“box‑intersection join”如何继承该正确性结构

Amagata 论文是 “点‑窗口” 的 range join；我们这里是 “盒‑盒相交” join。
 但我们仍可以继承同一条理论主线：

> 构造一个候选超集 $U\supseteq J$，并让算法对 $U$ 上任意候选 pair 的 **提案概率为常数**；
>  然后用谓词 `Intersect` 做 rejection（只接受属于 $J$ 的 pair），则条件化后输出在 $J$ 上均匀；由于每次迭代独立，最终样本 i.i.d.。

这就是本 baseline（AGR‑BoxJoin）的设计核心：

- “网格 + 上界 μ” 用于构造候选超集并实现加权抽样；
- “rejection” 用于把分布纠正回真实 join。

------

## 3. 核心数据结构（可实现级描述）

下面固定记号：

- 查询侧：$A$（我们按 $r\in A$ 来驱动采样/枚举）
- 被索引侧：$B$（我们对 $B$ 建网格索引并从中抽 $s$）

> 工程建议：可以交换 $A/B$（见 3.8 与 5.4），以提升上界紧致度与接受率。

------

### 3.1 Box 结构与数值约定

每个 box 记录：

- `id`：唯一编号（建议 int64）
- `L[1..d]`, `R[1..d]`：端点

数值建议：

- 尽量用 **int64**（例如将原浮点坐标统一乘以 $10^k$ 并四舍五入成整数），避免浮点比较导致边界不一致；
- 所有区间均按半开 $[L,R)$ 语义处理，与 `Intersect` 的严格不等式一致。

------

### 3.2 代表点 $\mathrm{rep}(s)$

我们为每个 $s\in B$ 选一个代表点（用于网格映射）。默认：
$$
\mathrm{rep}(s) := (L_1(s),\dots,L_d(s))\quad\text{(lower‑left / minimum corner)}.
$$
为什么选 lower‑left？

- 每个 $s$ 的 lower‑left 是确定且易算；
- 用它可以构造一个**简单且可证明覆盖的扩张区域** $E(r)$，从而保证不漏真 join（见 6.1 的覆盖引理）。

> 可扩展讨论：也可以用中心点、upper‑right 等，但需要同步修改覆盖证明与扩张方式。Baseline 主体固定用 lower‑left，避免引入额外复杂性。

------

### 3.3 稀疏哈希网格 CellMap（对 $B$ 建索引）

为每个维度设定网格边长 $g_i>0$（各维可不同），定义 cell key：
$$
k_i(x)=\left\lfloor \frac{x_i}{g_i}\right\rfloor,\qquad 
\mathbf{k}(x)=(k_1(x),\dots,k_d(x)).
$$
cell 的几何区域定义为半开：
$$
\text{cell}(\mathbf{k})=\prod_{i=1}^d [k_i g_i,\ (k_i+1)g_i).
$$
数据结构：

- `CellMap[key] -> vector<sid>`：存放所有满足 $\mathbf{k}(\mathrm{rep}(s))=\text{key}$ 的 $s\in B$
- `CellSize[key] = |CellMap[key]|`

**负坐标注意**：必须使用数学意义的 floor（向 $-\infty$ 取整）。不同语言对负数整除/取模规则不同，需要显式实现。否则 key 会错位，进而影响覆盖性与性能（虽然“分布均匀性”的证明形式仍成立，但实际会漏 join 或误判）。

------

### 3.4 上界 μ 的构造：最大长度 $M_i$ 与扩张区域 $E(r)$

#### (1) 预处理：被索引侧每维最大长度

$$
M_i := \max_{s\in B}\big(R_i(s)-L_i(s)\big),\quad i=1..d.
$$

#### (2) 对每个 $r\in A$ 的扩张候选区域

$$
E(r):=\prod_{i=1}^d[\,L_i(r)-M_i,\ R_i(r)\,).
$$

这与草稿一致，但升级版会在分析部分给出更严格的“覆盖引理”证明与半开边界说明。

------

### 3.5 候选 cell 集合 $C(r)$（务必用“cell 相交”定义）

为了让“提案概率为常数”的证明与实现完全一致，我们在升级版里**统一**采用如下定义：
$$
C(r):=\big\{\mathbf{k}\ \big|\ \mathbf{k}\in \text{KeysNonEmpty},\ \text{cell}(\mathbf{k})\cap E(r)\neq\varnothing\big\}.
$$
也就是说：

- $C(r)$ 是 **非空 cell** 的集合（来自 CellMap 的 key 集合）；
- 保留条件是 cell 的区域与 $E(r)$ 相交（而不是只看某个点是否落入）。

这种定义有两大好处：

1. 与代码逻辑一致（我们实际是“选 cell 再选 cell 内的 s”）；
2. 直接导出候选支持集 $U_{\text{prop}}$ 的“提案概率常数”性质（见 6.2）。

**计算 key 范围（每一维）**：
 令 $E_i(r)=[a_i,b_i)=[L_i(r)-M_i,\ R_i(r))$。
 则一维上与 $[a_i,b_i)$ 相交的 cell index 范围为：
$$
k_i^{\min}=\left\lfloor\frac{a_i}{g_i}\right\rfloor,\qquad
k_i^{\max}=\left\lceil\frac{b_i}{g_i}\right\rceil-1.
$$

------

### 3.6 上界与权重

对每个 $r$：
$$
\mu(r)=\sum_{\mathbf{k}\in C(r)} |CellMap[\mathbf{k}]|,\qquad 
\mu(r,\mathbf{k})=|CellMap[\mathbf{k}]|.
$$
直观解释：$\mu(r)$ 等于“所有候选 cell 内对象数之和”，即候选超集里 $r$ 可配到的候选 $s$ 的数量（注意：候选超集可能包含大量 false positives，这是 baseline 的性能瓶颈来源）。

------

### 3.7 Alias 结构（Walker alias）

为了常数时间 weighted sampling，构建两级 alias：

- `AliasR`：在 $A'=\{r\in A\mid \mu(r)>0\}$ 上建 alias，权重为 $\mu(r)$
- `AliasCell[r]`：对每个 $r\in A'$，在 $C(r)$ 上建 alias，权重为 $|CellMap[\mathbf{k}]|$

若 $\mu(r)=0$，则 $r$ 不可能产生任何候选 pair，直接剔除（既省空间又避免无意义采样）。

> 论文中明确使用 Walker alias 来实现按权重的离散抽样，并强调其 $O(1)$ 抽样代价与 $O(n)$ 构建代价。

------

### 3.8（强烈推荐写进报告）计算 $C(r)$ 的两种实现口径

这是你们草稿中容易被审稿人质疑、也是我之前给你们扣分的点：**“我们只保留非空 cell”不等于“我们不需要扫描空 cell”。**

升级版明确给两种口径，读者可按工程实现选择其一：

**实现 A：朴素枚举 key 超矩形（简单但可能爆炸）**

- 枚举所有组合 $\prod_i[k_i^{\min},k_i^{\max}]$，对每个 key 查询 CellMap 是否存在。
- 复杂度与该超矩形体积成正比：可能非常大（尤其是高维/小 cell）。
   适合：低维（2D/3D）且参数调得很粗的情况。

**实现 B：稀疏 cell slab 索引（推荐，仍保持 baseline 简洁）**
 增加一个辅助结构：

- `SlabIndex[k1] -> vector<key>`：按第一维的 $k_1$ 将所有非空 cell key 分桶存起来（存全 key 向量）。
- 对每个 $r$：只遍历 $k_1\in[k_1^{\min},k_1^{\max}]$ 的这些桶，把 key 再按其它维检查是否落在 $[k_i^{\min},k_i^{\max}]$ 内。

这样构造 $C(r)$ 的时间是：
$$
O\!\left(\sum_{k_1=k_1^{\min}}^{k_1^{\max}} |\text{SlabIndex}[k_1]|\right),
$$
即只与“这些 slab 里实际存在的非空 cell 数”相关，避免扫描大量空 cell。

> 这一改动非常“baseline 友好”：仍是 hash‑grid 思路，没有引入复杂的 range tree；但能显著减少空 cell 扫描开销，从而让 baseline 更像一个可落地系统组件（也是“90+”的关键补强点）。

------

## 4. 算法详细流程（三版本）

> 三版本完全对应你们草稿：Sampling（AGR‑S）、Enumerate+Sampling（AGR‑ES）、Adaptive+Sampling（AGR‑AS）。
>  升级点在于：**支持集定义统一**、**空 join 安全终止策略明确**、**伪代码与工程细节更完整**。

------

### 4.0 通用预处理 `Preprocess(A,B,g)`

输入：$A,B$，网格参数 $g_1..g_d$
 输出：$M_i$、`CellMap`、每个 $r$ 的 $C(r),\mu(r)$、alias 结构

步骤：

1. **计算 $M_i$**：扫描 $B$，得到每维最大长度 $M_i$。
2. **建立 CellMap**：对每个 $s\in B$，计算 $\mathbf{k}(\mathrm{rep}(s))$，插入 `CellMap[key]`。
3. **构造 KeysNonEmpty 与可选 SlabIndex**：
   - `KeysNonEmpty = keys(CellMap)`
   - 可选：构建 `SlabIndex[k1]`（推荐）。
4. 对每个 $r\in A$**构造候选 cell 集合 $C(r)$**：
   - 计算 $E(r)$ 与每维 key 范围 $k_i^{\min},k_i^{\max}$；
   - 用“实现 A 或 B”生成 $C(r)\subseteq KeysNonEmpty$；
   - 计算 $\mu(r)=\sum_{k\in C(r)}|CellMap[k]|$。
5. 对每个 $\mu(r)>0$ 的 $r$，构建 `AliasCell[r]`；
6. 在所有 $\mu(r)>0$ 的 $r$ 上建 `AliasR`（权重 $\mu(r)$）；
7. 计算 $\mathrm{MuSum}=\sum_{r\in A}\mu(r)$。

**必要条件快速判空**：若 `MuSum==0`，则必有 $|J|=0$。
 （注意：`MuSum>0` 并不能推出 $|J|>0$，只是“可能非空”。）

------

### 4.1 版本一：Sampling（AGR‑S）

**用途定位**：当你确信 $|J|>0$（或容许 join 为空时返回空前先跑 Adaptive），AGR‑S 提供最简 i.i.d. uniform 采样。

**流程**：

1. `Preprocess(A,B,g)`
2. 若 `MuSum==0`：返回空
3. 循环直到输出 $t$ 个样本：
   - $r\leftarrow$ `AliasSample(AliasR)`
   - $\mathbf{k}\leftarrow$ `AliasSample(AliasCell[r])`
   - $s\leftarrow$ 在 `CellMap[k]` 中均匀选一个元素
   - 若 `Intersect(r,s)`：accept 并输出 $(r,s)$，否则 reject 继续

伪代码：

```
AGR-S(A,B,t):
  Preprocess(A,B)
  if MuSum == 0: return empty

  Ans = []
  while |Ans| < t:
      r = AliasSample(AliasR)
      k = AliasSample(AliasCell[r])
      s = UniformPick(CellMap[k])
      if Intersect(r,s):
          Ans.append((r,s))
  return Ans
```

**安全终止声明（必须写进报告）**：

- 若 $|J|>0$，则接受率 $p_{\text{acc}}>0$，算法以概率 1 终止且期望迭代有限（见 6.2）。
- 若 $|J|=0$ 但 `MuSum>0`，则 $p_{\text{acc}}=0$，AGR‑S 将无限循环。
   因此：**默认推荐在系统里用 AGR‑AS（Adaptive）作为 baseline 主入口**；AGR‑S 作为其“大 join 分支”的子过程使用（此时已知 join 非空）。

------

### 4.2 版本二：Enumerate+Sampling（AGR‑ES）

**用途定位**：当 join 结果小/可存储时，AGR‑ES 直接枚举 $J$ 并在数组上均匀采样，是最直观的正确 baseline。

**枚举阶段**（构造 `AllPairs`）：
 对每个 $r\in A$：

1. 计算 $C(r)$（同预处理逻辑）
2. 对每个 $\mathbf{k}\in C(r)$：遍历 `CellMap[k]` 中每个 $s$
3. 若 `Intersect(r,s)`：把 pair push 到 `AllPairs`

由于每个 $s$ 只属于一个 cell（按 rep 映射），对固定 $r$ 不会在不同 cell 里重复扫描同一个 $s$，因此不会产生重复 pair。

**采样阶段**（数组均匀）：

- 若 `W = |AllPairs| == 0`：返回空
- 否则对 $j=1..t$：独立均匀抽 `idx_j ∈ [0, W-1]`，返回 `AllPairs[idx_j]`

伪代码：

```
AGR-ES(A,B,t):
  Build CellMap, M_i (alias 可不建)
  AllPairs = []
  for r in A:
      build E(r), enumerate candidate cells C(r)
      for k in C(r):
          for s in CellMap[k]:
              if Intersect(r,s):
                  AllPairs.push((r,s))

  if |AllPairs| == 0: return empty
  Ans = []
  for j=1..t:
      idx = UniformInt(0, |AllPairs|-1)
      Ans.append(AllPairs[idx])
  return Ans
```

------

### 4.3 版本三：Adaptive+Sampling（AGR‑AS，默认推荐）

**用途定位**：这是 baseline 的“可落地版本”，解决两件事：

- join 为空：能安全终止并返回空；
- join 很大：避免存储爆炸，自动切换到 AGR‑S。

**输入额外参数**：阈值 $J_\star$。

**流程**：

1. 先 `Preprocess(A,B,g)`（因为若转 Sampling 需要 alias）
2. 初始化 `AllPairs=[]`
3. 枚举真实 join pair：
   - 每发现一个 `Intersect(r,s)` 为真，就 push 到 `AllPairs`
   - 若 `AllPairs.size > J_star`：
     - 立刻 `AllPairs.clear()`（避免前缀偏差争议）
     - 直接调用 `AGR-S`（复用已建 alias），输出 $t$ 个样本并返回
4. 若枚举结束仍未超过阈值：
   - 若 `AllPairs` 为空：返回空
   - 否则在 `AllPairs` 上数组均匀采样输出 $t$

伪代码：

```
AGR-AS(A,B,t,J_star):
  Preprocess(A,B,g)

  AllPairs = []
  for r in A:
      enumerate C(r)
      for k in C(r):
          for s in CellMap[k]:
              if Intersect(r,s):
                  AllPairs.push((r,s))
                  if |AllPairs| > J_star:
                      AllPairs.clear()
                      return AGR-S-UsingExistingAlias(t)

  if |AllPairs| == 0: return empty
  return ArrayUniformSampling(AllPairs, t)
```

**关键点（写论文一定要写）**：

- 一旦切换，完全丢弃 AllPairs 前缀，确保输出分布与“纯 AGR‑S”一致，不受枚举路径影响。
- 一旦切换发生，说明至少发现 $J_\star+1$ 个真实 pair，因此 $|J|>0$，从而 AGR‑S 的终止性得到保证。

------

## 5. Adaptive 阈值 $J_\star$ 的选择策略（更论文/工程化）

这一节给一个“可写进论文、可复现实验”的阈值策略：**硬约束 + 软约束 + 标定**，并附带诊断工具。

### 5.1 内存硬约束（必须）

设给 `AllPairs` 的内存预算为 `MemBudgetPairs` 字节；每个 pair 仅存两个 int64 id，则理论 16B（实际 vector 有额外开销，可按 24B/32B 估算更保守）。
$$
J_\star^{\text{mem}}=\left\lfloor \frac{\text{MemBudgetPairs}}{\text{BytesPerPair}} \right\rfloor,\qquad
J_\star\le J_\star^{\text{mem}}.
$$
论文写法建议：明确你们实验机器内存、`BytesPerPair` 取值和 `MemBudgetPairs` 比例（比如 0.1×RAM）。

------

### 5.2 时间控制（避免“枚举到阈值后丢弃”浪费太大）

在大 join 情况下，AGR‑AS 会先枚举到 $J_\star$ 再丢弃，这段浪费大致 $\propto J_\star$。一个常用工程约束是：
$$
J_\star^{\text{time}}=\alpha\cdot t,
$$
$\alpha$ 是小常数（如 1/2/5/10），表示你最多愿意为 $t$ 个样本付出多少“无效枚举”级别的额外工作。

最终取：
$$
J_\star=\min(J_\star^{\text{mem}},\ J_\star^{\text{time}}).
$$

------

### 5.3 基于 benchmark 的交叉点标定（推荐写进论文）

为了让阈值选择不显得“拍脑袋”，推荐在实验 setup 里加入：

1. 固定 $t$，在代表性数据集上跑 AGR‑ES 与 AGR‑S
2. 找到它们运行时间相近的 join size 区间，记为 $|J_{\text{cross}}|$
3. 令 $J_\star\approx |J_{\text{cross}}|$（同时受 $J_\star^{\text{mem}}$ 约束）

这能让审稿人接受“阈值来自经验标定”而不是人为挑选。

------

### 5.4 Pilot 接受率诊断（用于调参，不影响正确性）

定义接受率：
$$
p_{\text{acc}}=\Pr(\text{accept})=\frac{|J|}{\mathrm{MuSum}}.
$$
在不知 $|J|$ 的情况下，你可以运行一个与最终采样**随机源隔离**的 pilot（例如做 10k 次 proposal），用 `accept/try` 估计 $p_{\text{acc}}$，用于：

- 调整网格尺度 $g_i$（太小→$C(r)$ 多、预处理/枚举重；太大→cell 太大、μ 松、接受率低）
- 决定是否交换 $A/B$（选择更“紧”的一侧做被索引侧 $B$，通常能减小 $M_i$ 从而减小 $E(r)$）

这些操作不改变正确性，只提升效率。

------

## 6. 算法分析（正确性 + 复杂度；三版本都包含）

### 6.1 覆盖性引理：不漏真 join（关键）

**引理 1（覆盖性）**
 令 $\mathrm{rep}(s)=L(s)$，且 $M_i=\max_{s\in B}(R_i(s)-L_i(s))$。
 若 $\texttt{Intersect}(r,s)$ 为真，则
$$
\mathrm{rep}(s)\in E(r)=\prod_{i=1}^d [L_i(r)-M_i,\ R_i(r)).
$$
**证明（逐维）**：固定维度 $i$。相交意味着
$$
\max(L_i(r),L_i(s))<\min(R_i(r),R_i(s)).
$$
特别地，有 $L_i(s)。另一方面，相交也意味着 $L_i(r)，即
$$
L_i(s)>L_i(r)-(R_i(s)-L_i(s))\ge L_i(r)-M_i.
$$
合并两式得到
$$
L_i(r)-M_i < L_i(s) < R_i(r),
$$
从而 $L_i(s)\in [L_i(r)-M_i,\ R_i(r))$ 成立（左端严格大于，但落在半开区间内）。对所有维度同时成立，故 $\mathrm{rep}(s)\in E(r)$。□

**推论 1**：若 $(r,s)\in J$，则 $\text{cell}(\mathbf{k}(\mathrm{rep}(s)))\cap E(r)\neq\varnothing$，因此 $\mathbf{k}(\mathrm{rep}(s))\in C(r)$。
 这确保 join 的每个真 pair 都在我们的候选支持集中（下一节将用到）。

------

### 6.2 AGR‑S 的正确性：uniform & i.i.d.（关键修正版）

你们草稿中最容易被挑刺的点是候选集合 $U$ 的定义。升级版统一使用**与实现一致**的支持集：

- 定义候选支持集（算法可能提出的所有 pair）：

$$
U_{\text{prop}}=\{(r,s)\mid r\in A,\ s\in B,\ \mathbf{k}(\mathrm{rep}(s))\in C(r)\}.
$$

显然
$$
|U_{\text{prop}}|=\sum_{r\in A}\mu(r)=\mathrm{MuSum}.
$$
并由 6.1 推论得 $J\subseteq U_{\text{prop}}$。

------

#### 引理 2（提案概率为常数）

对任意 $(r,s)\in U_{\text{prop}}$，设 $\mathbf{k}=\mathbf{k}(\mathrm{rep}(s))$，一次 proposal 提出该 pair 的概率为：
$$
\Pr(\text{propose }(r,s))
=
\underbrace{\frac{\mu(r)}{\mathrm{MuSum}}}_{\text{选 }r}
\cdot
\underbrace{\frac{|CellMap[\mathbf{k}]|}{\mu(r)}}_{\text{选 cell}}
\cdot
\underbrace{\frac{1}{|CellMap[\mathbf{k}]|}}_{\text{cell 内均匀选 }s}
=\frac{1}{\mathrm{MuSum}}.
$$
与具体 $(r,s)$ 无关。□

------

#### 定理 1（条件于 accept，输出在 $J$ 上均匀）

接受谓词为 `Intersect(r,s)`，因此 accept 等价于 “pair 属于 $J$”。
 一次 proposal 的接受概率：
$$
p_{\text{acc}}=\Pr(\text{accept})=\sum_{(r,s)\in J}\Pr(\text{propose }(r,s))
=|J|\cdot\frac{1}{\mathrm{MuSum}}=\frac{|J|}{\mathrm{MuSum}}.
$$
于是对任意 $P\in J$：
$$
\Pr(\text{output }P\mid \text{accept})
=
\frac{\Pr(\text{propose }P)}{\Pr(\text{accept})}
=
\frac{1/\mathrm{MuSum}}{|J|/\mathrm{MuSum}}
=\frac1{|J|}.
$$
因此每个输出样本在 $J$ 上均匀。□

------

#### 定理 2（独立性：i.i.d.）

每轮 proposal 使用独立随机源生成 $(r,s)$，accept/reject 仅依赖该轮 $(r,s)$。
 标准 rejection sampling 结论：被接受的样本序列是 i.i.d.，分布为“proposal 分布在 accept 条件下的条件分布”，即上面的均匀分布。□

------

#### 终止性说明

- 若 $|J|>0$，则 $p_{\text{acc}}>0$，期望 proposal 次数为 $t/p_{\text{acc}}=t\cdot \mathrm{MuSum}/|J|$。
- 若 $|J|=0$，则 $p_{\text{acc}}=0$，AGR‑S 不终止。
   因此实验/系统默认采用 AGR‑AS 来保证安全终止。

------

### 6.3 AGR‑ES 的正确性：枚举正确 + 数组采样 i.i.d.

**命题 1（枚举不漏）**：任意 $(r,s)\in J$，由覆盖引理，$\mathbf{k}(\mathrm{rep}(s))\in C(r)$，因此在枚举 $r$ 时一定会扫描到 $s$ 所在 cell 并做精确相交判定输出。

**命题 2（枚举不重）**：每个 $s$ 只属于一个 cell（按 rep 定义），对固定 $r$ 每个 cell 只扫描一次，因此不会重复输出同一 $(r,s)$。

因此 `AllPairs` 恰好枚举 $J$。随后对数组下标独立均匀采样，自然得到 i.i.d. uniform with replacement。□

------

### 6.4 AGR‑AS 的正确性：切换不引入偏差 + 终止安全

AGR‑AS 有两分支：

- 未切换：等价 AGR‑ES → i.i.d. uniform
- 切换：丢弃 `AllPairs` 前缀，完全执行 AGR‑S → i.i.d. uniform

切换条件只依赖枚举计数与阈值，不会对 AGR‑S 的 proposal 分布施加任何条件（因为我们丢弃前缀），因此切换不引入偏差。
 并且切换发生时已发现至少 $J_\star+1$ 个真 pair，保证 $|J|>0$，从而 AGR‑S 必然终止。□

------

### 6.5 复杂度分析（三版本分别给出）

#### 6.5.1 记号

- $n_1=|A|,\ n_2=|B|,\ n=n_1+n_2$
- 非空 cell 数：$m_{\text{cell}}=|\text{KeysNonEmpty}|\le n_2$
- $C_r=|C(r)|$
- $\mu(r)=\sum_{k\in C(r)}|CellMap[k]|$
- $\mathrm{MuSum}=\sum_{r\in A}\mu(r)=|U_{\text{prop}}|$
- $|J|$：真实 join size
- $p_{\text{acc}}=|J|/\mathrm{MuSum}$（当 $|J|>0$）

------

#### 6.5.2 预处理 `Preprocess`

- 计算 $M_i$：$O(n_2\cdot d)$
- 建 CellMap：$O(n_2)$
- 构造所有 $C(r)$、计算 $\mu(r)$：取决于你用的实现口径：

**实现 A（朴素枚举 key 超矩形）**
 令
$$
V(r)=\prod_{i=1}^d \bigl(k_i^{\max}-k_i^{\min}+1\bigr)
$$
则总时间 $O(\sum_{r\in A} V(r))$，可能非常大。

**实现 B（SlabIndex 稀疏过滤，推荐）**
 总时间
$$
O\!\left(\sum_{r\in A}\sum_{k_1=k_1^{\min}}^{k_1^{\max}} |\text{SlabIndex}[k_1]|\right)
$$
≈ 遍历“落在该 slab 范围内的非空 cell 数”。

- alias 构建：`AliasR` 为 $O(n_1)$，所有 `AliasCell[r]` 总计 $O(\sum_r C_r)$

空间：

- CellMap：$O(n_2)$
- 存储所有 $C(r)$+alias：$O(n_1+\sum_r C_r)$

------

#### 6.5.3 AGR‑S（Sampling）期望时间

一次 proposal 成本：

- alias 抽样 $O(1)$
- cell 内均匀取样 $O(1)$
- `Intersect` 判定 $O(d)$

所以每次 proposal 为 $O(d)$。期望 proposal 次数为 $t/p_{\text{acc}}=t\cdot \mathrm{MuSum}/|J|$。因此：
$$
T_{\text{AGR-S}}=
T_{\text{prep}}
+
O\!\left(t\cdot d\cdot\frac{\mathrm{MuSum}}{|J|}\right).
$$
空间（除索引外）：$O(t)$。

------

#### 6.5.4 AGR‑ES（Enumerate+Sampling）时间与空间

枚举扫描候选对总数为 $\sum_r \mu(r)=\mathrm{MuSum}$，每次检查 $O(d)$：
$$
T_{\text{AGR-ES}}=
T_{\text{prep (或部分 prep)}}
+
O(d\cdot \mathrm{MuSum})
+
O(t).
$$
空间：
$$
S_{\text{AGR-ES}}=O(n_2+n_1+\sum_r C_r+|J|).
$$

------

#### 6.5.5 AGR‑AS（Adaptive）两情形

- 若未切换（$|J|\le J_\star$）：
  $$
  T= T_{\text{prep}} + O(d\cdot \mathrm{MuSum}) + O(t),
  \quad
  S=O(n_2+n_1+\sum_r C_r+|J|).
  $$
  （$\mathrm{MuSum}$ 仍是扫描候选的代价上界。）

- 若切换（$|J|>J_\star$）：

  - 枚举阶段最多存到 $J_\star$ 就停止存储，随后进入 AGR‑S

  $$
  T=
  T_{\text{prep}}
  +
  O(d\cdot N_{\text{scan}}(J_\star))
  +
  O\!\left(t\cdot d\cdot\frac{\mathrm{MuSum}}{|J|}\right),
  $$

  其中 $N_{\text{scan}}(J_\star)$ 为“为找到 $J_\star$ 个真相交 pair 所需扫描的候选对数”（依赖数据分布，最坏可到 $\mathrm{MuSum}$）。
  $$
  S=O(n_2+n_1+\sum_r C_r+\max(J_\star,t)).
  $$