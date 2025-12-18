# AGR‑BoxJoin Baseline（Amagata’25 风格：网格上界 μ + Alias 加权 + Rejection）

> Baseline 目的：提供一个“**严格 i.i.d. 均匀有放回**”的 join sampling 方案，思想源自 Amagata’25 的“上界 μ + rejection 修正仍能保证 uniform & independent”的框架。
>  这个 baseline **非常适合作对照**：实现简单、理论干净，但性能强依赖上界紧致程度（接受率）。

------

## 1. 问题定义与分析 + 文献引用

### 1.1 Spatial Join Sampling（与你们问题对齐）

在 $d\ge2$ 维欧氏空间 $\mathbb{R}^d$ 中给定两组轴对齐半开盒子（MBR / hyper-rectangle）：
$$
R_c=\{r_{c1},\dots,r_{c n_1}\},\qquad
R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\}.
$$
每个盒子：
$$
r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\quad L_i(r)<R_i(r).
$$
空间相交连接（intersection join / filter step）输出集合：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
相交判定（半开）等价于对所有维度 $i$：
$$
\max(L_i(r_c),L_i(r_{\bar c}))<\min(R_i(r_c),R_i(r_{\bar c})).
$$
采样目标：输出 $t$ 个样本对
$$
Z_1,\dots,Z_t\in J
$$
并满足 **i.i.d. 均匀有放回**：
$$
\Pr(Z_j=P)=\frac1{|J|}\ \ (\forall P\in J),\qquad Z_1,\dots,Z_t \text{ 相互独立}.
$$
这正是你们文稿设定的标准目标。

> 备注：若 $|J|=0$，则“在空集合上均匀采样”在数学上无定义；工程上一般返回空并标记 join 为空。后文会说明如何在 baseline 内做到**安全终止**。

------

### 1.2 为什么不能“先均匀选 r，再均匀选其匹配 s”

令
$$
\deg(r)=|\{s\in R_{\bar c}\mid r\cap s\neq\varnothing\}|.
$$
若你先均匀选 $r$，再在其匹配集合中均匀选 $s$，则 pair 的输出概率会与 $\deg(r)$ 强相关，无法保证对所有 pair 均匀（会偏向大度数对象）。因此必须引入**按匹配数（或其可用上界）加权**的机制，这也是 Amagata’25 明确强调的核心原因之一。

------

### 1.3 文献：Amagata’25 可被我们继承的关键点

**论文**：Daichi Amagata, *Random Sampling over Spatial Range Joins* (2025).

你上面已经抓到了它对我们最关键的三点（这里系统化写成“可继承条目”）：

1. **采样目标必须 uniform（并按问题设定独立）**。文中定义强调返回 $t$ 个来自 $J$ 的样本，每个样本“picked uniformly at random”。
2. **需要 weighted sampling**（例如 Walker alias）来处理不同对象匹配数的差异。
3. **用易算上界 $\mu(r)\ge \deg(r)$ + rejection sampling 仍能保证 uniform & independent**，这是 KDS-rejection 一类方法的理论核心。
4. **主框架**可以抽象为：
    “按 $\mu(r)$ 选 $r$ → 按 $\mu(r,c)$ 选 cell → 从 cell 内抽 $s$ → 若满足谓词则 accept，否则 reject”。

------

### 1.4 从 “range join sampling” 迁移到 “box intersection join sampling”

Amagata’25 的思路本质是：构造一个**易采样**的候选超集 $U\supseteq J$，并让在 $U$ 上的“提案概率”对每个候选对为常数，再通过“相交谓词”做 rejection 把分布修正回 $J$ 上的均匀分布。

对盒子相交 join，我们用“**代表点 + 扩张候选区域**”实现 $U$ 的构造（下一节详述），并保持相同的正确性结构。

------

## 2. 核心数据结构

下面固定记号：

- 查询侧 $A:=R_c$
- 被索引侧 $B:=R_{\bar c}$（建立网格索引的一侧）
   （实现时可交换 $A/B$，见 §4.4 的工程建议。）

------

### 2.1 Box 结构

每个盒子存：

- `id`：唯一编号（建议 int64）
- `L[1..d]`、`R[1..d]`：端点（建议统一缩放成 int64；浮点会引入边界歧义）

半开相交判定（建议作为唯一判定函数）：
$$
\texttt{Intersect}(r,s)=\bigwedge_{i=1}^d \left(\max(L_i(r),L_i(s)) < \min(R_i(r),R_i(s))\right).
$$

------

### 2.2 代表点 rep(s)

默认选择被索引侧盒子的 **lower-left（最小端点向量）**：
$$
\mathrm{rep}(s):=(L_1(s),\dots,L_d(s)).
$$

> 你在原稿里也提到可换成 center，但需要同步修改上界构造；这可以作为“扩展讨论”，不影响 baseline 主体正确性。

------

### 2.3 稀疏哈希网格 CellMap

为每个维度设定 cell 边长 $g_i>0$（可不同维度不同）。

cell key 定义为整数向量 $\mathbf{k}=(k_1,\dots,k_d)$，其中
$$
k_i=\left\lfloor \frac{\mathrm{rep}(s)_i}{g_i}\right\rfloor.
$$
数据结构：

- `CellMap[ key ] -> vector<sid>`：该 cell 内所有 $s\in B$ 的 id 列表
- `CellSize[key] = |CellMap[key]|`（或用列表长度即时取）

> **注意负坐标**：必须使用数学意义的 floor。不同语言对负数整除的行为不同，需要显式处理，否则会产生系统性偏差（虽然不影响均匀性证明的形式，但会影响正确性/覆盖性）。

------

### 2.4 上界相关：候选区域 E(r)、候选 cell 集合 C(r)、上界 μ(r)

**(1) 预处理：被索引侧每维最大长度**
$$
M_i := \max_{s\in B}\big(R_i(s)-L_i(s)\big),\quad i=1..d.
$$
你前面推导用到这一步（并用于保证“不漏真解”）。

**(2) 候选区域（保证覆盖）**

对每个 $r\in A$ 定义：
$$
E(r):=\prod_{i=1}^d [\,L_i(r)-M_i,\ R_i(r)\,).
$$
（这就是你原稿里的关键构造。）

**(3) 候选 cell 集合**
$$
C(r):=\{\mathbf{k}\mid \text{cell}(\mathbf{k})\cap E(r)\neq\varnothing\}.
$$
实现上先算每一维的 cell index 范围：
$$
k_i^{\min}=\left\lfloor\frac{L_i(r)-M_i}{g_i}\right\rfloor,\qquad
k_i^{\max}=\left\lceil\frac{R_i(r)}{g_i}\right\rceil-1.
$$
然后枚举所有组合 $\prod_i[k_i^{\min},k_i^{\max}]$，并仅保留 `CellMap` 中存在的非空 cell key（稀疏网格）。

记 $C_r:=|C(r)|$，它决定预处理与 per‑r alias 的开销。

**(4) 上界与权重**
$$
\mu(r):=\sum_{\mathbf{k}\in C(r)} |CellMap[\mathbf{k}]|,\qquad
\mu(r,\mathbf{k}) := |CellMap[\mathbf{k}]|.
$$
并由覆盖性质推出 $\deg(r)\le \mu(r)$。

------

### 2.5 Alias 结构（Walker alias）

为了常数时间 weighted sampling（与 Amagata’25 叙述一致）：

- 全局 alias：`AliasR`，元素为 $r\in A$，权重 $\mu(r)$
- 每个 $r$ 的 alias：`AliasCell[r]`，元素为 $\mathbf{k}\in C(r)$，权重 $\mu(r,\mathbf{k})=|CellMap[\mathbf{k}]|$
- 若 $\mu(r)=0$，则该 $r$ 从 `AliasR` 中剔除（它不可能产生样本）。

------

### 2.6（Enumerate/Adaptive 用）AllPairs 缓冲区

- `AllPairs`: 动态数组，存储真实 join pair `(rid, sid)`
- Adaptive 版本最多存 $J_\star$ 个 pair

------

## 3. 算法详细流程（Sampling / Enumerate+Sampling / Adaptive+Sampling）

先给一个**通用预处理**，之后给三版本流程。为保证可复现，所有遍历顺序都固定（按 id、按维度序等）。

------

### 3.0 通用预处理 Preprocess(A,B)

**输入**：$A,B$，维度 $d$，网格参数 ${g_i}$
 **输出**：`CellMap`、$M_i$、对每个 $r$ 的 $C(r),\mu(r)$，alias 结构

步骤：

1. 计算 $M_i=\max_{s\in B}(R_i(s)-L_i(s))$。
2. 建 `CellMap`：对每个 $s\in B$，插入到 `CellMap[key(rep(s))]`。
3. 对每个 $r\in A$：

- 构造 $E(r)$，枚举非空候选 cells 得到 $C(r)$；
- 计算 $\mu(r)=\sum_{\mathbf{k}\in C(r)}|CellMap[\mathbf{k}]|$；
- 若 $\mu(r)>0$：构建 `AliasCell[r]`。

1. 在所有 $\mu(r)>0$ 的 $r$ 上按权重 $\mu(r)$ 建 `AliasR`。
2. 计算 $\sum_{r \in A} \mu(r)=\sum_r\mu(r)$（可选但推荐：用于空 join 的“必要条件判断”）。

> 若 $\sum_{r \in A} \mu(r)=0$，则 $|J|=0$（因为对任意 $r$，候选区域里都没有任何 $s$，不可能相交）。

------

### 3.1 版本一：Sampling（AGR‑S）

> 核心：**两级 alias + cell 内均匀 + 相交谓词 rejection**，完全对应 Amagata’25 的主框架。

#### 流程

1. `Preprocess(A,B)`
2. 若 `MuSum==0`：返回空（join 必为空）
3. 重复直到输出 $t$ 个样本：

- $r \leftarrow \textsf{WeightedSample}(\texttt{AliasR})$（权重 $\mu(r)$）
- $\mathbf{k}\leftarrow \textsf{WeightedSample}(\texttt{AliasCell}[r])$（权重 cell size）
- $s \leftarrow$ 从 `CellMap[k]` 中均匀取一个元素
- 若 `Intersect(r,s)`：accept，输出 $(r,s)$；否则 reject，继续

这一步与原始 baseline 描述一致。

#### 伪代码（Sampling）

```
AGR-Sampling(A,B,t):
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

#### 重要备注：如何保证“安全终止”

- 若 $|J|>0$：accept 概率 $p_{\text{acc}}=|J|/\sum_{r \in A} \mu(r)>0$，算法几乎必然终止，且期望迭代次数有限。
- 若 $|J|=0$ 但 $\sum_{r \in A} \mu(r)>0$：accept 概率为 0，纯 rejection 会无限循环。
   **因此实际工程建议默认使用 Adaptive 版本（§3.3）**：它会先尝试枚举确认是否为空/是否很小，超过阈值则再转 Sampling；此时若发生转 Sampling，必然已经找到了至少 $J_\star+1$ 个真实 join pair，保证 $|J|>0$。

------

### 3.2 版本二：Enumerate+Sampling（AGR‑ES）

> 核心：用同一套 “候选区域 E(r) + 网格 cell 扫描 + 精确相交判定” **显式枚举全部 $J$** 并存入数组，然后做数组下标的 i.i.d. 均匀采样。

#### 枚举阶段（生成 AllPairs）

对每个 $r\in A$：

1. 计算 $E(r)$、枚举候选 cell $C(r)$
2. 对每个候选 cell $\mathbf{k}\in C(r)$：遍历 `CellMap[k]` 内所有 $s$
3. 若 `Intersect(r,s)`：把 pair push 到 `AllPairs`

由于每个 $s$ 只属于一个 cell，且对固定 $r$ 每个 cell 只扫一次，因此不会出现“同一 (r,s) 在不同 cell 重复输出”的情况。

#### 采样阶段（数组采样）

令 $W=|AllPairs|=|J|$：

- 若 $W=0$：返回空

- 否则对 $j=1..t$：
  $$
  \text{idx}_j\sim \mathrm{Unif}\{0,\dots,W-1\},\qquad \texttt{Ans}[j]=AllPairs[\text{idx}_j]
  $$

得到 i.i.d. 均匀有放回样本。

#### 伪代码（Enumerate+Sampling）

```
AGR-EnumerateSampling(A,B,t):
  Build CellMap and M_i (alias 可不建)

  AllPairs = []
  for r in A (fixed order):
     build E(r), enumerate candidate cells C(r)
     for k in C(r):
        for s in CellMap[k]:
           if Intersect(r,s):
              AllPairs.push((r,s))

  if |AllPairs| == 0: return empty
  Ans[1..t]
  for j=1..t:
      idx = UniformInt(0, |AllPairs|-1)
      Ans[j] = AllPairs[idx]
  return Ans
```

------

### 3.3 版本三：Adaptive+Sampling（AGR‑AS）

> 核心：引入阈值 $J_\star$
>
> - 若 $|J|\le J_\star$：直接枚举并数组采样（快且稳定）
> - 若 $|J|>J_\star$：一旦发现超过阈值，立刻丢弃 `AllPairs`，转 AGR‑S（避免存储爆炸）

#### 流程

1. 先 `Preprocess(A,B)`（要建 alias，可能要转 Sampling）
2. 初始化 `AllPairs=[]`
3. 枚举真实 join pair，并计数：

- 发现一个真相交 pair：push 到 `AllPairs`
- 若 `AllPairs.size > J_*`：
  - 释放/清空 `AllPairs`（避免任何“前缀偏差”风险）
  - 直接运行 AGR‑S 输出 $t$ 个样本并返回

1. 若枚举结束仍未超过阈值：

- 若 `AllPairs` 为空：返回空
- 否则在 `AllPairs` 上做数组均匀采样输出 $t$

#### 伪代码（Adaptive）

```
AGR-Adaptive(A,B,t,J_star):
  Preprocess(A,B)    // 建 CellMap, M_i, mu, AliasR, AliasCell

  AllPairs = []
  for r in A:
     for k in C(r):
        for s in CellMap[k]:
           if Intersect(r,s):
              AllPairs.push((r,s))
              if |AllPairs| > J_star:
                 AllPairs.clear()
                 return AGR-Sampling-UsingExistingAlias(A,B,t)  // 直接采样

  // 未切换
  if |AllPairs| == 0: return empty
  return ArrayUniformSampling(AllPairs, t)
```

> 关键点：**一旦切换，完全不用 AllPairs 前缀结果**，以保证“分布完全等同于 AGR‑S”而不受“枚举到阈值的过程”影响。

------

## 4. Adaptive 版本阈值 $J_\star$ 的选择策略

你的问题要求“阈值选择策略”要写得像论文一样可落地，这里给一个**三层策略**：先保内存，再控时间，再可选用 pilot 标定常数。

------

### 4.1 内存硬约束（必须）

设可给 `AllPairs` 的内存预算为 `MemBudgetPairs` 字节；每个 pair 存储成本约 `sizeof(Pair)`（比如两个 int64 是 16B，外加 vector 开销）。
$$
J_\star^{\text{mem}}=\left\lfloor \frac{\text{MemBudgetPairs}}{\text{sizeof(Pair)}} \right\rfloor.
$$
必须保证
$$
J_\star \le J_\star^{\text{mem}}.
$$

------

### 4.2 时间控制（避免大 join 时“枚举到阈值后丢弃”太浪费）

当 $|J|$ 很大时，Adaptive 会先枚举到 $J_\star$ 再丢弃，然后再做 sampling 输出 $t$ 个样本。

一个非常实用的经验约束是：
$$
J_\star^{\text{time}} = \alpha \cdot t,
$$
其中 $\alpha$ 是小常数（例如 1、2、5、10），表示你最多愿意在大 join 情况下多做 $\alpha t$ 量级的无效枚举。

直觉（与 Amagata’25 的 acceptance 结构兼容）：Sampling 的期望 proposal 次数随 $t$ 线性增长，而“先枚举再丢弃”的浪费 roughly 与 $J_\star$ 线性增长，所以限制 $J_\star$ 与 $t$ 同阶可以避免切换浪费压过采样主体。

最终建议：
$$
J_\star=\min(J_\star^{\text{mem}},\ J_\star^{\text{time}}).
$$

------

### 4.3 可选：基于 benchmark 的交叉点标定（更论文/工程化）

你可以在实验设置里给出一个“阈值默认值”的来源：

1. 在代表性数据集上固定 $t$
2. 分别测 AGR‑ES（全枚举）与 AGR‑S（纯 sampling）的耗时
3. 观测两者耗时相当的 $|J|$ 区间，记为 $|J_{\text{cross}}|$
4. 令 $J_\star^{\text{time}}\approx |J_{\text{cross}}|$（并受内存上界约束）

------

### 4.4 可选：pilot 估计接受率，用于诊断/调参

Amagata’25 风格 rejection 的核心指标是接受率：
$$
p_{\text{acc}}=\Pr(\text{accept})=\frac{|J|}{\sum_r\mu(r)}.
$$
并且每个样本的期望迭代次数是
$$
\mathbb{E}[\#\text{iters per sample}] = \frac{\sum_r \mu(r)}{|J|}. 
$$


虽然 $|J|$未知，但你可以做一个与最终采样**随机源隔离**的 pilot（不把 pilot 产生的样本用于最终输出），用 accept/try 的比例估计 $p_{\text{acc}}$，来决定是否需要调整 $g_i$、是否交换 $A/B$、是否做更紧的上界（例如局部上界）。这不改变正确性，只影响效率。

------

## 5. 算法分析（三版本：正确性 + 复杂度）

### 5.1 正确性：关键引理（覆盖性 + 上界）

**引理 1（覆盖性：不漏真解）**
 对被索引侧 $B$ 定义 $M_i=\max_{s\in B}(R_i(s)-L_i(s))$，并令 $\mathrm{rep}(s)=L(s)$。若 $r\cap s\neq\varnothing$，则
$$
\mathrm{rep}(s)\in E(r)=\prod_{i=1}^d[\,L_i(r)-M_i,\ R_i(r)\,).
$$
证明就是你原稿里那段逐维推导。

**推论**：对每个 $r$，真实匹配数
$$
\deg(r)\le \mu(r)=\sum_{\mathbf{k}\in C(r)}|CellMap[\mathbf{k}]|.
$$


------

### 5.2 Sampling（AGR‑S）正确性：uniform & i.i.d.

定义候选集合
$$
U=\{(r,s)\mid r\in A,\ s\in B,\ \mathrm{rep}(s)\in E(r)\}.
$$
由引理 1 得 $J\subseteq U$。

**引理 2（提案概率为常数）**
 对任意 $(r,s)\in U$，设 $s$ 所在 cell 为 $\mathbf{k}$，一次迭代提出该 pair 的概率：
$$
\Pr(\text{propose }(r,s))
=
\underbrace{\frac{\mu(r)}{\sum_{r'}\mu(r')}}_{\text{选 }r}
\cdot
\underbrace{\frac{\mu(r,\mathbf{k})}{\mu(r)}}_{\text{选 cell}}
\cdot
\underbrace{\frac{1}{|CellMap[\mathbf{k}]|}}_{\text{cell 内均匀选 }s}
=
\frac{1}{\sum_{r'}\mu(r')},
$$
与具体 $(r,s)$ 无关。

**定理 1（条件于 accept，输出在 $J$ 上均匀）**
 一次迭代条件于 accept 时：
$$
\Pr(\text{output }(r,s)\mid \text{accept})=\frac1{|J|}.
$$
证明结构与 Amagata’25 版本一致：accept 概率为 $|J|/\sum\mu$，做条件化即得 $1/|J|$。

**定理 2（独立性：i.i.d.）**
 每次尝试使用独立随机过程产生提案，accept/reject 仅依赖本次提案，重复直到收集到 $t$ 个 accept，属于标准 rejection sampling 流程，因此 $t$ 个输出相互独立同分布。

> 唯一需要额外声明的是：若 $|J|=0$，accept 概率为 0，算法不会终止；因此工程上应默认用 Adaptive 或先做“空 join 判定”。

------

### 5.3 Enumerate+Sampling（AGR‑ES）正确性

- 枚举阶段输出的集合等于 $J$：因为对任意相交 pair $(r,s)$，由引理 1 有 $\mathrm{rep}(s)\in E(r)$，因此 $s$ 一定会出现在 $r$ 扫描的某个候选 cell 的列表中，并会被精确相交判定输出。
- 数组采样阶段：独立均匀抽数组下标，得到 i.i.d. 均匀有放回样本。

------

### 5.4 Adaptive（AGR‑AS）正确性

Adaptive 只有两种分支：

- 未切换：等价 AGR‑ES → i.i.d. uniform
- 切换：丢弃枚举前缀，完全按 AGR‑S 重新采样 → i.i.d. uniform

分支选择由数据与阈值决定（不依赖采样随机性），因此整体输出仍严格 i.i.d. uniform。

------

### 5.5 复杂度分析（三版本都给）

记号：

- $n_1=|A|,\ n_2=|B|,\ n=n_1+n_2$
- $C_r=|C(r)|$（候选非空 cell 数）
- $\mu(r)=\sum_{\mathbf{k}\in C(r)}|CellMap[\mathbf{k}]|$
- $\sum_{r \in A} \mu(r)=\sum_{r\in A}\mu(r)$
- $|J|$ 为真实 join size
- $p_{\text{acc}}=|J|/\sum_{r \in A} \mu(r)$（若 $|J|>0$）

#### 5.5.1 预处理（共同）

- 计算 $M_i$：$O(n_2\cdot d)$
- 建 CellMap：$O(n_2)$
- 枚举 $C(r)$ 并求 $\mu(r)$：$O(\sum_{r\in A} C_r)$（按稀疏哈希查 cell 是否非空）
- alias 构建：全局 $O(n_1)$，每个 $r$ 为 $O(C_r)$，总计 $O(\sum_r C_r)$

空间：

- CellMap：$O(n_2)$
- 存 $C(r)$/alias：$O(n_1+\sum_r C_r)$

这些与您原始 baseline 复杂度写法一致。

------

#### 5.5.2 Sampling（AGR‑S）期望时间

一次 proposal 代价：alias 抽样 $O(1)$ + 相交测试 $O(d)$，所以 $O(d)$。

每个样本的期望迭代次数为 $\sum_{r \in A} \mu(r)/|J|$，因此输出 $t$ 个样本的期望时间：
$$
T_{\text{AGR-S}} 
= T_{\text{prep}} + O\!\left(t\cdot d\cdot \frac{\sum_{r \in A} \mu(r)}{|J|}\right).
$$
这正是你原稿（并对齐 Amagata’25）里写的形式。

空间（不含索引本体）：$O(t)$，索引与 alias 为 $O(n_2+n_1+\sum_r C_r)$。

------

#### 5.5.3 Enumerate+Sampling（AGR‑ES）时间与空间

枚举阶段需要检查全部候选对数量 $\sum_{r \in A} \mu(r)$，每次检查 $O(d)$：
$$
T_{\text{enum}}=T_{\text{prep}}+O(d\cdot \sum_{r \in A} \mu(r)).
$$
存储全部输出对：
$$
S_{\text{AGR-ES}}=O(n_2+n_1+\sum_r C_r+|J|).
$$
数组采样额外 $O(t)$。

------

#### 5.5.4 Adaptive（AGR‑AS）两情形

- 若未切换（$|J|\le J_\star$）：

$$
T=T_{\text{prep}}+O(d\cdot \sum_{r \in A} \mu(r)+t),\quad
S=O(n_2+n_1+\sum_r C_r+|J|)\le O(\dots+J_\star).
$$

- 若切换（$|J|>J_\star$）：

先枚举直到发现 $J_\star+1$ 个真实 join pair（此时至少保证 $|J|>0$），再运行 AGR‑S：
$$
T
=
T_{\text{prep}}
+ O(d\cdot N_{\text{scan}}(J_\star))
+ O\!\left(t\cdot d\cdot \frac{\sum_{r \in A} \mu(r)}{|J|}\right),
$$
其中 $N_{\text{scan}}(J_\star)$ 为“为了找到 $J_\star$ 个真相交 pair，需要扫描的候选对数量”（依赖数据分布；最坏可到 $\sum_{r \in A} \mu(r)$）。

空间峰值：
$$
S=O(n_2+n_1+\sum_r C_r+\max(J_\star,t)).
$$

------

### 5.6（写实验讨论时很好用）优点与局限

**优点（baseline 非常干净）**

- 结构极简：hash‑grid + alias + 常数采样 + 相交判定。
- 输出严格满足 i.i.d. uniform（理论证明短且直观）。
- 空间不依赖 $|J|$（Sampling/Adaptive 的大分支），对比“枚举型 baseline”很有说服力。

**局限（也是它作为 baseline 的意义）**

- 速度强依赖上界紧致程度：上界松 → 接受率低 → 迭代次数爆炸，这是  讨论的核心风险点之一。
- 若 $M_i$ 很大导致 $E(r)$ 过大，$\mu(r)$ 可能接近 $n_2$，退化为“几乎全局拒绝采样”。
- 高维下候选 cell 枚举可能指数膨胀（$d$ 很大时不友好）。