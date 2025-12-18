# TS23‑TLSOP Spatial Join Sampling Baseline（Two‑layer SOP, Upgraded）

> **定位**：这是一个“强枚举型、工程可复现”的 baseline。
>  Join 生成部分复现 Tsitsigkos’23 的核心结论：**two‑layer SOP 的 A/B/C/D 分类 + 仅评估 9 个 mini‑join，即可在不做 reference‑point duplicate elimination 的情况下枚举无重复 join 结果**。
>  在此之上，我们补齐三种采样版本，均实现目标：**对 join 结果集合 $J$ 的 i.i.d. 均匀有放回采样**。Tsitsigkos’23

------

## 1. 问题定义与分析（含引用文章内容）

### 1.1 任务定义（与 Spatial Join Sampling 对齐）

给定两组空间对象（仅考虑其 MBR，等价于过滤步 filter step）：

- $R=\{r_1,\dots,r_{n_R}\}$
- $S=\{s_1,\dots,s_{n_S}\}$

每个对象用二维轴对齐矩形（MBR）表示（建议统一用**半开区间**，避免贴边歧义）：
$$
r=[x_l(r),x_u(r))\times[y_l(r),y_u(r)),\quad
s=[x_l(s),x_u(s))\times[y_l(s),y_u(s)).
$$
过滤步的相交连接结果：
$$
J=\{(r,s)\mid r\in R,\ s\in S,\ r\cap s\neq \varnothing\}. 
$$
两矩形相交（半开）等价于：
$$
\max(x_l(r),x_l(s))<\min(x_u(r),x_u(s)) \ \land\ \max(y_l(r),y_l(s))<\min(y_u(r),y_u(s)).
$$

### 1.2 采样目标（Sampling）

输出 $t$ 个样本对 $Z_1,\dots,Z_t\in J$，要求：

- **均匀性**：$\Pr(Z_j=P)=1/|J|$（任意 $P\in J$）
- **独立性**：$Z_1,\dots,Z_t$ i.i.d.，允许重复抽到同一 pair（with replacement）Tsitsigkos’23

### 1.3 Baseline 论文来源与“可直接继承”的关键点（Tsitsigkos’23）

**论文信息（写进论文/实验设置）**
 *Two-layer Space-oriented Partitioning for Non-point Data* — Dimitrios Tsitsigkos, Panagiotis Bouros, Konstantinos Lampropoulos, Nikos Mamoulis, Manolis Terrovitis（Tsitsigkos’23）。Tsitsigkos’23

**本文 baseline 继承点（按论文贡献拆解）**

1. **Two‑layer SOP（主分区 + 次分区）**：在 SOP（例如规则网格）基础上，每个 tile 内把对象按“begin 值是否在 tile 起点之后”分为 A/B/C/D 四类并分别存储。
2. **16 个 mini‑join 中 7 个只产 duplicates**：因此 tile 内只需评估剩余 **9 个 mini‑join** 即得到**无重复 join 结果**，无需 reference‑point 去重。
3. **mini‑join 枚举器**：提供 plane‑sweep（排序 + forward scan）的实现框架（Algorithm 2），并给 reduced/batch 优化（Algorithm 3/4）。Tsitsigkos’23
4. **可推广到更高维**：次分区可扩展到 $m$ 维（需 $2^m$ 类）。Tsitsigkos’23

> 我们的补强点：把上述 join 枚举“工程化”为一个**确定性 JoinStream**，然后在 JoinStream 之上给出三种采样版本（Sampling / Enumerate+Sampling / Adaptive+Sampling），并补齐 correctness/complexity 与阈值策略。Tsitsigkos’23

------

## 2. 核心数据结构

### 2.1 主分区：规则网格 tiles（SOP + multi‑assignment）

设空间域（或数据包围盒）为：
$$
\mathcal{D}=[X_{\min},X_{\max})\times[Y_{\min},Y_{\max}).
$$
选网格分辨率 $(n_x,n_y)$（列数/行数），tile：
$$
T_{i,j}=[X_i,X_{i+1})\times[Y_j,Y_{j+1}),\quad 0\le i<n_x,\ 0\le j<n_y.
$$
**multi‑assignment**：对象分配到所有与其 MBR 相交的 tile（SOP 多重分配）。

> 工程建议：用“两遍分配”（先统计每 tile/每类容量，再一次性分配数组）可以避免大量 push_back 扩容；但 baseline 里简单 vector 也可。

### 2.2 次分区：tile 内四类 A/B/C/D（Two‑layer SOP）

对 tile $T$，记起点为 $(T.x_l,T.y_l)$。对被分配到 $T$ 的矩形 $r$，按 $x_l(r),y_l(r)$ 与起点比较分类：

- **A 类**：$x_l(r)\ge T.x_l$ 且 $y_l(r)\ge T.y_l$
- **B 类**：$x_l(r)\ge T.x_l$ 且 $y_l(r)< T.y_l$
- **C 类**：$x_l(r)< T.x_l$ 且 $y_l(r)\ge T.y_l$
- **D 类**：$x_l(r)< T.x_l$ 且 $y_l(r)< T.y_l$

于是每个 tile 存 8 个列表：
$$
R_T^A,R_T^B,R_T^C,R_T^D;\quad S_T^A,S_T^B,S_T^C,S_T^D
$$

### 2.3 关键：9 个 mini‑join 组合（无重复枚举）

每个 tile $T$ 只需评估下列 9 个 mini‑join（其他 7 个跳过，因为只产 duplicates）：

1. $R_T^A \Join S_T^A$
2. $R_T^A \Join S_T^B$
3. $R_T^A \Join S_T^C$
4. $R_T^A \Join S_T^D$
5. $R_T^B \Join S_T^A$
6. $R_T^B \Join S_T^C$
7. $R_T^C \Join S_T^A$
8. $R_T^C \Join S_T^B$
9. $R_T^D \Join S_T^A$

### 2.4 mini‑join 枚举器：Plane‑sweep（Algorithm 2 思路）

每个 mini‑join 使用 plane‑sweep（排序 + forward scan）枚举相交对；并可选 reduced/batch 优化减少排序与比较（Algorithm 3/4）。Tsitsigkos’23

**为保证 Sampling 两遍一致性**，需要固定：

- 排序键（例如按 $x_l$）
- tie‑break（例如按对象 id）
- 遍历顺序（从小到大）

### 2.5 Sampling 相关辅助结构（按版本）

- **Sampling（2‑pass Index）**：`Idx[1..t]`（随机索引及其样本位置）、`Ans[1..t]`
- **Enumerate+Sampling**：`AllPairs[0..|J|-1]`
- **Adaptive**：阈值 $J_\star$、`mode ∈ {ENUMERATE, COUNT_ONLY}`、`AllPairs`（最多临时存到 $J_\star$）

------

## 3. 算法详细流程（Sampling / Enumerate+Sampling / Adaptive+Sampling）

### 3.1 JoinStream：把论文 join 枚举“工程化”为确定性输出流

目标：构造生成器
$$
\texttt{JOIN\_STREAM()} \to (r,s)\in J \quad(\text{每次 yield 一个 pair})
$$
并保证每个 pair 恰好一次且顺序可复现。

#### 3.1.1 建索引 `BuildIndex(R,S)`

1. 建 $n_x\times n_y$ tiles
2. 对每个 $r\in R$：计算覆盖 tile 集合；对每个 tile 按 A/B/C/D 规则插入到 $R_T^*$
3. 对每个 $s\in S$：同理插入到 $S_T^*$
4. （推荐）对每个 tile、每个 class 列表做一次排序（如按 $x_l$，tie by id），后续复用

> 复用排序很关键：例如 $R_T^A$ 会参与 4 个 mini‑join，如果不复用会重复排序 4 次。

#### 3.1.2 mini‑join 内枚举：`PlaneSweepMiniJoin(A,B)`

一个 baseline 级的 deterministic 枚举器（足够写论文、可复现）：

- 输入：两个已排序列表 $A,B$（按 $x_l$，tie by id）
- 输出：枚举所有 $(a,b)$ 满足矩形相交

实现上可以用典型 sweepline 思路（维护一个“仍可能与当前 a 相交”的 B 候选容器），并对候选做 y 检查。

> 论文里给的是更具体的 plane‑sweep 伪代码（Algorithm 2）以及 reduced/batch 优化（Algorithm 3/4）；你们实现 baseline 时可以先上 Algorithm 2 的简单版本，再逐步加优化。Tsitsigkos’23

#### 3.1.3 JoinStream 的全局顺序（必须固定）

- tile 顺序：按 tile id（如 row-major）
- tile 内 mini‑join 顺序：固定为这 9 个组合
- mini‑join 内输出顺序：由排序键 + tie‑break 决定

这样两次 `JOIN_STREAM()` 的输出顺序完全一致，才能做 index sampling。Tsitsigkos’23

------

### 3.2 版本一：Sampling（Two‑Pass Index Sampling）

**名字**：TS23‑TLSOP‑2PassIndex（对应你要求的 Sampling 版本）

**核心思想**：

- Pass1：只计数 $N=|J|$（不存 pair）
- 生成 $t$ 个独立均匀随机下标（有放回）
- Pass2：再次枚举 JoinStream，命中下标时输出对应 pair

该思路在 Tsitsigkos’23 baseline 文档中已有描述，我们这里补齐“确定性要求 + 伪代码 + 正确性/复杂度”。

#### 3.2.1 流程

1. `Index = BuildIndex(R,S)`
2. **Pass‑1 Count**：`N=0; for pair in JOIN_STREAM(Index): N++`
3. 生成 i.i.d. 下标 $K_1,\dots,K_t \sim \mathrm{Unif}\{1..N\}$，并把 $(K_j,j)$ 按 $K_j$ 排序（便于第二遍单次扫描回填）Tsitsigkos’23
4. **Pass‑2 Select**：再次扫描 JoinStream，用计数器 idx=1..N 编号流元素，命中则写入 `Ans[j]`
5. 输出 `Ans`

#### 3.2.2 伪代码

```
TS23_TLSOP_2PassIndex(R, S, t):

1) Index = BuildIndex(R, S, n_x, n_y)

2) // Pass-1: count
   N = 0
   for pair in JOIN_STREAM(Index):
       N++

   if N == 0: return empty

3) // draw indices with replacement
   for j = 1..t:
       K[j] ~ UniformInt(1, N)
   sort (K[j], j) by K ascending

4) // Pass-2: select by index
   idx = 0
   p = 1
   for pair in JOIN_STREAM(Index):
       idx++
       while p <= t and K_sorted[p].value == idx:
           Ans[K_sorted[p].pos] = pair
           p++
       if p > t: break

5) return Ans
```

------

### 3.3 版本二：Enumerate+Sampling（显式枚举后数组采样）

**名字**：TS23‑TLSOP‑EnumerateThenSample（对应你要求的 Enumerate+Sampling 版本）

**核心思想**：
 只做**一遍** JoinStream，把所有 pairs 存进 `AllPairs`（大小 $N=|J|$），然后对数组做 i.i.d. 均匀采样。

#### 3.3.1 流程

1. `Index = BuildIndex(R,S)`
2. 枚举 JoinStream：`AllPairs.append(pair)`
3. 采样：对每个 j，取 `idx ~ Unif{0..N-1}`，输出 `AllPairs[idx]`

#### 3.3.2 伪代码

```
TS23_TLSOP_EnumerateThenSample(R, S, t):

1) Index = BuildIndex(R, S, n_x, n_y)

2) AllPairs = []
   for pair in JOIN_STREAM(Index):
       AllPairs.append(pair)

   N = |AllPairs|
   if N == 0: return empty

3) for j = 1..t:
       idx ~ UniformInt(0, N-1)
       Ans[j] = AllPairs[idx]

4) return Ans
```

> 这个版本通常是最快的 baseline（只扫一遍 join），但空间 $O(|J|)$ 可能爆。

------

### 3.4 版本三：Adaptive+Sampling（阈值自适应）

**名字**：TS23‑TLSOP‑Adaptive（对应你要求的 Adaptive+Sampling 版本）

**目标**：

- 若 $|J|\le J_\star$：退化为 Enumerate+Sampling（**1 pass**，最快）
- 若 $|J|>J_\star$：退化为 Sampling（**2 pass**，但不存 $J$，避免爆内存）

#### 3.4.1 关键切换逻辑（推荐实现）

第一遍扫描 JoinStream 时：

- 永远做 `N++`
- 若 `mode==ENUMERATE` 且 `N <= J_star`：把 pair push 到 `AllPairs`
- 一旦 `N > J_star`：
  - `mode = COUNT_ONLY`
  - 立刻释放 `AllPairs`
  - 之后只计数，不再存 pair（仍然会枚举 pair 来计数——baseline 本质上仍依赖输出规模）

扫描结束得到精确 $N=|J|$。

- 若未切换：直接数组采样
- 若已切换：生成随机索引并做第二遍索引定位（同 Sampling）

#### 3.4.2 伪代码

```
TS23_TLSOP_Adaptive(R, S, t, J_star):

1) Index = BuildIndex(R, S, n_x, n_y)

2) mode = ENUMERATE
   AllPairs = []
   N = 0

   for pair in JOIN_STREAM(Index):
       N++
       if mode == ENUMERATE:
           if N <= J_star:
               AllPairs.append(pair)
           else:
               mode = COUNT_ONLY
               AllPairs.clear()   // release memory

3) if N == 0: return empty

4) if mode == ENUMERATE:
       // N <= J_star, AllPairs is full J
       for j = 1..t:
           idx ~ UniformInt(0, N-1)
           Ans[j] = AllPairs[idx]
       return Ans

5) // mode == COUNT_ONLY: do 2-pass index sampling
   for j = 1..t:
       K[j] ~ UniformInt(1, N)
   sort (K[j], j) by K ascending

   idx = 0
   p = 1
   for pair in JOIN_STREAM(Index):
       idx++
       while p <= t and K_sorted[p].value == idx:
           Ans[K_sorted[p].pos] = pair
           p++
       if p > t: break

6) return Ans
```

------

## 4. Adaptive 版本阈值 $J_\star$ 的选择策略

你需要一个能写进论文的策略：**硬约束（内存） + 软约束（时间交叉点/经验）**。

### 4.1 内存硬约束（必须满足）

设给 `AllPairs` 的可用内存预算为 `MemBudgetPairs` 字节。每个 pair 的存储开销近似：

- 两个 32-bit id：8 字节
- 两个 64-bit id / 指针：16 字节
   再乘容器与对齐开销系数 $c_{\text{over}}\in[1.2,2]$。

则：
$$
J_\star^{\text{mem}}
=
\left\lfloor
\frac{\text{MemBudgetPairs}}{\text{PairBytes}\cdot c_{\text{over}}}
\right\rfloor.
$$
推荐保守取：
$$
J_\star = \left\lfloor 0.7\cdot J_\star^{\text{mem}}\right\rfloor
$$
给索引结构、排序缓冲、线程等留余量。

### 4.2 时间软约束（建议做一次离线标定）

对给定硬件与实现：

- Enumerate+Sampling：1 次 JoinStream + 写 `AllPairs`
- Sampling：2 次 JoinStream + 排序 $t$ 个索引

如果你们发现写 `AllPairs` 的内存带宽成本很高（尤其是大对象/GC 语言），则应把 $J_\star$ 设小一些，让大规模时倾向两遍不物化。

一个稳妥写法：

1. 选 2–3 个代表性数据集做 pilot
2. 找到 Enumerate 与 2‑pass 耗时交叉点对应的 $|J|$
3. 把交叉点作为 $J_\star^{\text{time}}$，最终：

$$
J_\star = \min(J_\star^{\text{mem}},\ J_\star^{\text{time}}).
$$

> 如果你们不想引入 $J_\star^{\text{time}}$，只用内存上界也完全合理（论文里写“由内存预算决定”即可）。

------

## 5. 算法分析（正确性 + 复杂度；三版本都包含）

### 5.1 正确性

#### 引理 1：JoinStream 不重不漏（每个 pair 恰好一次）

依据 Tsitsigkos’23 的 two‑layer SOP 结论：tile 内 16 个 class‑to‑class mini‑join 中有 7 个只会产生 duplicates，故只评估剩余 9 个 mini‑join 即可获得**无重复的 join 枚举**，无需 reference‑point 去重。
 因此 `JOIN_STREAM()` 的输出集合与 $J$ 一致，且每个 pair 只出现一次。Tsitsigkos’23

> 这一步对 sampling 至关重要：若 join 枚举有重复，则采样会对重复出现的 pair 过度采样导致偏差。Tsitsigkos’23

#### 引理 2：Sampling（2‑pass Index）输出 i.i.d. 均匀样本

把 JoinStream 的输出顺序视为对 $J$ 的一个固定编号：
$$
J=\{P_1,\dots,P_N\},\quad N=|J|.
$$
Sampling 版本生成独立均匀索引 $K_j\sim\mathrm{Unif}\{1..N\}$，并输出 $Z_j=P_{K_j}$。
 则对任意 $P_i\in J$：$\Pr(Z_j=P_i)=1/N$，且因 $K_j$ 独立，故 $Z_j$ i.i.d.。

#### 引理 3：Enumerate+Sampling 输出 i.i.d. 均匀样本

若 `AllPairs` 恰好包含每个 pair 一次，则对数组下标做独立均匀采样同样得到 i.i.d. 均匀样本（同上）。

#### 定理 1：Adaptive 的切换不影响正确性

Adaptive 只有两种终态：

- 未切换：等价于 Enumerate+Sampling
- 切换：等价于 Sampling（2‑pass Index）

切换仅影响“是否保留 AllPairs”，不改变 JoinStream，也不改变最终 $N=|J|$ 的精确计数与随机索引分布，因此输出仍为 $J$ 上 i.i.d. 均匀有放回样本。

------

### 5.2 复杂度

#### 5.2.1 记号

- $n=n_R+n_S$
- replication factor：$\rho_R,\rho_S$
- two‑layer 索引规模：$n' = n_R\rho_R+n_S\rho_S$
- 9 个 mini‑join 集合：$\mathcal{M}$
- tile $T$ 内各类规模：$|R_T^X|,|S_T^Y|$

#### 5.2.2 建索引（two‑layer grid）

时间：
$$
O(n_R\rho_R+n_S\rho_S)
$$
空间：
$$
O(n_R\rho_R+n_S\rho_S
$$

#### 5.2.3 一次 JoinStream 遍历代价 

论文/实现上可用“排序 + forward scan + 输出”的组合写成：
$$
T_{\text{stream}} = \sum_T \sum_{(X,Y)\in\mathcal{M}} \Big(\text{sort}(|R_T^X|,|S_T^Y|)+\text{sweep}(|R_T^X|,|S_T^Y|)\Big).
$$

#### 5.2.4 Sampling（2‑pass Index）

时间：
$$
T_{\text{Sampling}} = 2\cdot T_{\text{stream}} + O(t\log t)
$$
空间：

$$
S_{\text{Sampling}} = O(n' + t) 
$$

#### 5.2.5 Enumerate+Sampling 时间： $ T_{\text{Enum+Sampling}} = T_{\text{stream}} + O(|J|) + O(t)$

空间：
$$
S_{\text{Enum+Sampling}} = O(n' + |J| + t).
$$

#### 5.2.6 Adaptive+Sampling

- 若 $|J|\le J_\star$：等价 Enumerate+Sampling
  $$
  T = T_{\text{stream}} + O(|J|) + O(t),\quad
  S = O(n' + |J| + t)\le O(n' + J_\star + t)
  $$

- 若 $|J|>J_\star$：等价 2‑pass + 额外最多白存 $J_\star$ 个 pair
  $$
  T = 2\cdot T_{\text{stream}} + O(J_\star) + O(t\log t),\quad
  S_{\max} = O(n' + J_\star + t)
  $$

------

## 额外补强：复现实验必须写清的“确定性规则”

为了保证 Sampling/Adaptive 的两遍流一致（否则索引定位会错），建议在实现说明中显式写：

1. tile 顺序固定（tile id 升序）
2. mini‑join 顺序固定（9 个组合固定列表）Tsitsigkos’23
3. mini‑join 内排序稳定（按 $x_l$，tie by object id）
4. plane‑sweep 输出顺序完全由上述规则决定（不可用不确定迭代容器）Tsitsigkos’23