# 1. 问题定义与分析

## 1.1 输入：两类半开轴对齐盒子

在 $d\ge 2$ 维欧氏空间 $\mathbb{R}^d$ 中给定两类轴对齐**半开**盒集合：
$$
R_c=\{r_{c1},\dots,r_{c n_1}\},\quad
R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar c n_2}\},\quad n=n_1+n_2.
$$
每个盒子
$$
r=\prod_{i=1}^d [L_i(r),R_i(r)),\qquad L_i(r)<R_i(r).
$$
我们只关心**跨集合相交对**（spatial join）：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$

### 半开相交判定（维度独立）

半开区间下，一维相交：
$$
[a,b)\cap[c,d)\neq\varnothing \iff \max(a,c)<\min(b,d).
$$
因此两盒相交当且仅当每一维投影都相交。

------

## 1.2 输出：对 $J$ 的 i.i.d. 均匀有放回采样

输出 $t$ 个样本
$$
Z_1,\dots,Z_t\in J
$$
要求 **i.i.d. 均匀有放回**：
$$
\Pr(Z_j=P)=\frac{1}{|J|}\ \ (\forall P\in J),\qquad Z_1,\dots,Z_t\ \text{相互独立}.
$$

------

## 1.3 第 1 维 Plane Sweep：事件系统与严格“全序”（关键）

选择第 1 维 $x_1$ 扫描。每个盒子 $r$ 产生两个事件：

- $\texttt{START}(r)$ at $x_1=L_1(r)$
- $\texttt{END}(r)$ at $x_1=R_1(r)$

### 事件全序 $\prec$（必须写成**确定的总排序**）

对任意事件 $e$，定义键值：
$$
\mathrm{key}(e)=\big(x_1(e),\ \mathrm{typeRank}(e),\ \mathrm{classRank}(e),\ \mathrm{id}(e)\big)
$$
并用字典序比较。

推荐设置：

- $\mathrm{typeRank}(\texttt{END})=0,\ \mathrm{typeRank}(\texttt{START})=1$
   **同坐标 END 在 START 前**（与半开一致，避免“贴边误相交”）

- $\mathrm{classRank}$ 只在 “同坐标且都是 START” 时生效：例如
  $$
  \mathrm{classRank}(c)=0,\quad \mathrm{classRank}(\bar c)=1
  $$
  即**先处理 $R_c$ 的 START，再处理 $R_{\bar c}$** 的 START

- $\mathrm{id}(e)$（box-id）用于把排序变成严格全序（避免“仍然并列”的实现不确定性）

> 你原稿里写“固定 tie-break”是对的；这里把它升级成**可直接实现且可证明**的“全序定义”。

### 活跃集定义（与实现一致）

扫描按事件顺序逐个处理。定义在处理到某一时刻（某个事件之前）的活跃集：
$$
\mathrm{Active}(R_c)=\{r_c\in R_c \mid \texttt{START}(r_c)\ \text{已处理且}\ \texttt{END}(r_c)\ \text{未处理}\}
$$
对 $R_{\bar c}$ 同理。

等价地，在坐标意义上：当前处理 $\texttt{START}(q)$ 时，$x_0=L_1(q)$，对任意 $r\in\mathrm{Active}(\text{Other})$ 必有
$$
x_0\in [L_1(r),R_1(r))
$$
且不会出现 $x_0=R_1(r)$（因为 END 已先处理并删除）。

------

## 1.4 关键引理：START 处“维度分离”（只需检查维度 2..d）

**引理 1（维度分离）**
 处理某个 $\texttt{START}(q)$ 时，令 $x_0=L_1(q)$。对任意对立活跃盒 $r\in\mathrm{Active}(\text{Other})$，$q$ 与 $r$ 在第 1 维必然“真重叠”（非贴边）：
$$
[L_1(q),R_1(q))\cap[L_1(r),R_1(r))\neq\varnothing.
$$
**证明要点：**
 因为 $r$ 活跃，说明其 END 还未发生，因此 $x_0；而 $q$ 是 START，必有 $x_0。于是
$$
x_0<\min(R_1(q),R_1(r)),
$$
从而交集非空。∎

因此，在处理 $\texttt{START}(q)$ 时要判断 $q\cap r\neq\varnothing$，只需检查投影到维度 $2..d$ 的相交。

------

## 1.5 事件诱导分块：$J=\biguplus_e J_e$（唯一归属的严格版本）

对每个 START 事件 $e\in E$，设新盒子 $q=r(e)$，定义对立活跃集中的命中集合：
$$
K_e:=\{\, r\in \mathrm{Active}(\text{Other})\mid q^\star\cap r^\star\neq\varnothing \,\},
\qquad
w_e:=|K_e|.
$$
其中投影
$$
r^\star:=\prod_{i=2}^d [L_i(r),R_i(r)).
$$
为了统一输出顺序为 $(r_c,r_{\bar c})$，定义：

- 若 $q\in R_c$：$J_e=\{(q,r_{\bar c})\mid r_{\bar c}\in K_e\}$
- 若 $q\in R_{\bar c}$：$J_e=\{(r_c,q)\mid r_c\in K_e\}$

**命题 1（不交并分解 / 唯一归属）**
 在上述事件全序 $\prec$ 下：
$$
J=\biguplus_{e\in E}J_e,\qquad |J|=\sum_{e\in E}w_e.
$$
**证明要点（更严谨版）：**
 取任意相交对 $P=(r_c,r_{\bar c})\in J$。比较 $\texttt{START}(r_c)$ 与 $\texttt{START}(r_{\bar c})$ 在全序下谁更靠后，记为 $e^\star$。

- 在处理 $e^\star$ 时，另一个盒子必已 START 且其 END 未到（否则二者不可能相交），因此在对立活跃集中，且投影相交，因此 $P\in J_{e^\star}$（不漏）。
- 在处理较早的 START 时，较晚的那个盒子尚未入活跃集，因此不可能输出 $P$（不重）。
   故 $J$ 被唯一划分为 $\{J_e\}$ 的不交并。∎

------

# 2. 核心数据结构

Baseline 的核心是：**两棵动态 R-tree 维护对立活跃集在维度 $2..d$ 的投影**，并支持：

- `ReportIntersect(Q)`：枚举命中
- `CountIntersect(Q)`：精确计数
- `SampleIntersect(Q,t′)`：从命中集合 i.i.d. 均匀采样

> 这三者是三种版本算法（枚举 / 两遍采样 / 自适应）的最小支撑接口。

------

## 2.1 事件数组与索引（建议做到“可复现实现”）

维护排序后的事件数组：

- `Events[1..2n]`：每项为
  $$
  (x_1,\ type\in\{\texttt{START},\texttt{END}\},\ class\in\{c,\bar c\},\ box\_id)
  $$

- 排序严格使用 $\mathrm{key}(e)$（见 §1.3）

为了支持 Sampling / Adaptive 的 Phase2（alias）与 Phase3（第二遍回填），建议为每个 START 事件分配一个连续编号：

- `start_id(e) ∈ {1..n}`（按 Events 顺序第一次看到 START 就 ++）

并维护数组：

- `w[start_id]`：该 START 的权重 $w_e$
- `slots[start_id]`：该 START 被分到的输出位置列表（或计数 + 位置列表）

------

## 2.2 R-tree 只存维度 $2..d$ 投影：$r^\star$

如 §1.4，扫描保证第 1 维重叠，所以 R-tree 只存：
$$
r^\star=\prod_{i=2}^d [L_i(r),R_i(r)).
$$
维护两棵动态树：

- $\mathcal{RT}^{(c)}$：存 $\mathrm{Active}(R_c)$ 的投影条目 $(id,r^\star)$
- $\mathcal{RT}^{(\bar c)}$：存 $\mathrm{Active}(R_{\bar c})$

并在 sweep 中对 START/END 做 Insert/Delete。

------

## 2.3 半开语义：精确相交判定与剪枝判定（必须区分）

给定两个 $(d-1)$ 维半开盒
$$
A=\prod_{k=2}^d[\alpha_k,\beta_k),\qquad
B=\prod_{k=2}^d[\gamma_k,\delta_k),
$$
定义精确相交：
$$
\texttt{Intersect}^\star(A,B)\iff \forall k=2..d:\ \max(\alpha_k,\gamma_k)<\min(\beta_k,\delta_k).
$$

### 内部节点剪枝的安全性

- 用 `MBR(node)` 与 $Q$ 的相交来剪枝：若 `MBR(node)` 与 $Q$ 不相交，则子树无命中（安全剪枝，无漏报）。
- 用 `MBR(node) ⊆ Q` 来整子树接受：则该子树所有条目都与 $Q$ 相交（对半开集合也成立）。

> **关键要求：** 无论是 Report / Count / Sample，最终是否命中必须在叶层用 $\texttt{Intersect}^\star$ 过滤，避免“贴边”被误认为相交。

------

## 2.4 R-tree 结构增广：必须维护 `size(node)`（为采样服务）

为了让 `SampleIntersect` 能在 FullyAccepted 子树内均匀抽取，需要每个节点维护：

- `size(node)` = 该节点子树中**叶条目数**（即 active 条目数）

并在以下操作中保持一致：

- Insert：沿插入路径 `size += 1`，若分裂则重算两个新节点的 size 并继续向上修正
- Delete：沿删除路径 `size -= 1`，若触发 condense / reinsertion，同样修正 size

> 若 `size(node)` 不准确，采样分布会偏（这是 baseline 里最容易“隐形出错”的地方）。

------

## 2.5 Delete 的工程选择：handle vs FindLeaf（给出两种可落地方案）

### 方案 A：维护 handle（更快，但实现更复杂）

- `handle[id] = (leaf_node_ptr, entry_index)`
- 删除时 $O(1)$ 定位到叶并移除
- 但需要在**节点分裂、合并、重插入**时正确更新 handle

适合你们自己实现 R-tree 且愿意仔细处理指针有效性。

### 方案 B：不维护 handle（更稳健）

- `Delete(id, r*)` 使用标准 `FindLeaf` 搜索定位叶条目
- 平均代价接近一次查询（通常 $O(\log_B n)$），实现简单可靠

如果 baseline 目标是“稳 + 好复现”，推荐方案 B；如果目标是极致速度，用方案 A。

------

## 2.6 三个接口的“实现级”定义与伪代码

下面用统一函数：

- `Overlap*(A,B)`：半开精确相交（同 $\texttt{Intersect}^\star$）

- `Contained*(A,B)`：半开包含
  $$
  A\subseteq B \iff \forall k:\ L_k(A)\ge L_k(B)\ \wedge\ R_k(A)\le R_k(B)
  $$

### 2.6.1 ReportIntersect(Q)

返回命中集合：
$$
K(Q)=\{r\in\mathrm{Active}: r^\star\cap Q\neq\varnothing\}.
$$

```
ReportIntersect(node, Q):
    if not Overlap*(MBR(node), Q): return []
    if node is leaf:
        ans = []
        for entry in node.entries:
            if Overlap*(entry.mbr, Q):   # exact half-open
                ans.append(entry.id)
        return ans
    else:
        ans = []
        for child in node.children:
            ans.extend( ReportIntersect(child, Q) )
        return ans
```

### 2.6.2 CountIntersect(Q)（精确计数 + 关键剪枝）

```
CountIntersect(node, Q):
    if not Overlap*(MBR(node), Q): return 0
    if Contained*(MBR(node), Q): return size(node)
    if node is leaf:
        cnt = 0
        for entry in node.entries:
            if Overlap*(entry.mbr, Q):
                cnt += 1
        return cnt
    else:
        s = 0
        for child in node.children:
            s += CountIntersect(child, Q)
        return s
```

### 2.6.3 SampleIntersect(Q, t′)：Frontier 不交分解 + 加权抽组件 + 组件内均匀

核心思想：把 $K(Q)$ 分解成互不重叠的组件集合 $\mathcal{F}$：

- FullyAccepted 组件：某个节点 `node` 满足 `MBR(node) ⊆ Q`，对应集合 $S_C$ 是该子树所有条目，权重 $w_C=size(node)$
- LeafList 组件：某个叶节点里命中的 entry 列表，权重是列表长度

#### Step A：BuildFrontier(Q)

```
BuildFrontier(node, Q, F):
    if not Overlap*(MBR(node), Q): return
    if Contained*(MBR(node), Q):
        F.append( Component(type="FULL", node=node, w=size(node)) )
        return
    if node is leaf:
        hit = []
        for entry in node.entries:
            if Overlap*(entry.mbr, Q):
                hit.append(entry.id)
        if len(hit) > 0:
            F.append( Component(type="LEAF", list=hit, w=len(hit)) )
        return
    else:
        for child in node.children:
            BuildFrontier(child, Q, F)
```

**性质（必须在文档里写清楚）：**
$$
K(Q)=\biguplus_{C\in\mathcal{F}}S_C,\qquad |S_C|=w_C,\qquad \sum_C w_C=|K(Q)|.
$$
不交性来自：FullyAccepted 节点处停止向下；不同节点子树不交；不同叶的条目不交。

#### Step B：组件级抽样（alias + slot）

对 $\mathcal{F}$ 按权重 $w_C$ 构建 alias（或前缀和），独立抽 $t′$ 次组件，统计每个组件被抽中次数 `t_C`。

#### Step C：组件内均匀抽样

- 若组件是 `LEAF(list)`：直接在 `list` 中均匀取 1 个 id
- 若组件是 `FULL(node)`：在该子树内按 `size(child)` 比例随机下沉到叶，再在叶中均匀挑 entry

```
SampleSubtreeUniform(node):
    while node is not leaf:
        choose child i with prob size(child_i) / size(node)
        node = child_i
    return UniformChoice(node.entries).id
```

最终把每个组件贡献的 `t_C` 次样本合并输出，得到 **i.i.d. 均匀**序列（见 §5.1.4 证明）。

------

# 3. 算法详细流程（3 个版本）

三个版本共享同一事件系统与两棵 R-tree，仅在“是否枚举 / 是否需要第二遍 / 是否自适应切换”不同。

下面统一记：

- `RTc` = $\mathcal{RT}^{(c)}$，`RTb` = $\mathcal{RT}^{(\bar c)}$
- `OtherRT(q)`：若 $q\in R_c$ 则 OtherRT=RTb，否则 OtherRT=RTc
- `InsertSelf(q)`：插入到本侧树

------

## 3.1 版本 A：Enumerate+Sampling（显式枚举 $J$ + 数组采样）

**适用：** $|J|$ 小，可存下所有 pair。

### 3.1.1 步骤

**预处理：**

1. 构造 `Events`，按 §1.3 全序排序
2. `RTc, RTb` 初始化为空
3. `Pairs = []`

**一次扫描枚举：**

```
for e in Events:
    if e.type == END:
        if e.class == c:    RTc.Delete(e.box)
        else:               RTb.Delete(e.box)

    else: # START(q)
        q = e.box
        Q = q.star  # 投影到维度 2..d

        if q.class == c:
            cand = RTb.ReportIntersect(Q)
            for id in cand:
                # 这里 cand 已经是叶层 half-open 精确过滤后的命中
                Pairs.append( (q.id, id) )
            RTc.Insert(q)
        else:
            cand = RTc.ReportIntersect(Q)
            for id in cand:
                Pairs.append( (id, q.id) )
            RTb.Insert(q)
```

扫描结束：`W = len(Pairs) = |J|`。

**数组采样：**
 若 `W==0` 返回空；否则对 $j=1..t$：
$$
idx_j\sim \mathrm{Unif}\{0,\dots,W-1\},\quad Z_j=Pairs[idx_j].
$$

------

## 3.2 版本 B：Sampling（两遍扫描 + 事件 alias + 局部均匀采样）

**适用：** $|J|$ 大、不可枚举，但需要全局 i.i.d. 均匀样本。

### Phase 1：第一次扫描（只计数 $w_e$，得到 $W=|J|$）

初始化：`RTc, RTb` 为空，`W=0`，`start_id=0`。

```
for e in Events:
    if e.type == END:
        delete like version A
    else: # START
        start_id += 1
        q = e.box; Q = q.star
        if q.class == c:
            w[start_id] = RTb.CountIntersect(Q)
            W += w[start_id]
            RTc.Insert(q)
        else:
            w[start_id] = RTc.CountIntersect(Q)
            W += w[start_id]
            RTb.Insert(q)
```

若 `W==0`：返回空。

> 工程细节：构建 alias 时建议只保留 `w>0` 的 start 事件，避免 0 权重事件进入离散分布。

### Phase 2：事件级 alias + slot 分配

构建分布：
$$
p_e = \frac{w_e}{W}.
$$
对每个样本位置 $j\in\{1..t\}$：

1. 抽 `sid ~ p`（START 编号）
2. `slots[sid].append(j)`
    并记 `t_sid = len(slots[sid])`

### Phase 3：第二次扫描（按事件批量 SampleIntersect 回填）

重置：`RTc, RTb` 为空，`start_id=0`，输出数组 `Ans[1..t]`。

```
for e in Events:
    if e.type == END:
        delete
    else: # START
        start_id += 1
        q = e.box; Q = q.star
        need = len(slots[start_id])
        if need > 0:
            if q.class == c:
                list = RTb.SampleIntersect(Q, need)  # 返回 need 个 iid 样本 id
                for u in 1..need:
                    j = slots[start_id][u]
                    Ans[j] = (q.id, list[u])
            else:
                list = RTc.SampleIntersect(Q, need)
                for u in 1..need:
                    j = slots[start_id][u]
                    Ans[j] = (list[u], q.id)
        # 无论 need 是否为 0，都要插入 q
        if q.class == c: RTc.Insert(q)
        else:            RTb.Insert(q)
```

------

## 3.3 版本 C：Adaptive+Sampling（阈值 $J_\star$ 自适应切换）

**目标：** 小 $J$ 时“一遍扫 + 枚举存下来”，大 $J$ 时“避免爆炸 + 两遍采样”，且保证正确性不受影响。

### Adaptive Phase 1：一次扫描（始终 Count；在阈值内可选 Report）

初始化：

- `mode = ENUMERATE`
- `AllPairs = []`（最多存到 $J_\star$）
- `W=0`
- `RTc, RTb` 为空
- `w[start_id]` 仍需完整记录（为后续 alias）

```
start_id = 0
for e in Events:
    if e.type == END:
        delete
    else:
        start_id += 1
        q = e.box; Q = q.star

        # Step A: 永远先 Count，得到精确 w 和 W
        if q.class == c:
            w[start_id] = RTb.CountIntersect(Q)
        else:
            w[start_id] = RTc.CountIntersect(Q)
        W += w[start_id]

        # Step B: 阈值内可枚举，否则切换
        if mode == ENUMERATE:
            if W <= J_star:
                if q.class == c:
                    ids = RTb.ReportIntersect(Q)
                    for id in ids: AllPairs.append( (q.id, id) )
                else:
                    ids = RTc.ReportIntersect(Q)
                    for id in ids: AllPairs.append( (id, q.id) )
            else:
                mode = COUNT_ONLY
                AllPairs.clear()   # 丢弃已枚举结果，避免继续占内存

        # Step C: 插入
        if q.class == c: RTc.Insert(q)
        else:            RTb.Insert(q)
```

### Phase 1 后分支

- 若 `W==0`：返回空
- 若 `mode==ENUMERATE`（等价于 $W\le J_\star$）：
   `AllPairs` 就是完整 $J$，直接数组采样（同版本 A）
- 若 `mode==COUNT_ONLY`（等价于 $W>J_\star$）：
   执行版本 B 的 Phase2 + Phase3（注意：`w[]` 已经有了）

------

# 4. Adaptive 阈值 $J_\star$ 的选择策略（更可落地版）

阈值目标：

- **内存不爆**：`AllPairs` 不 OOM
- **时间不浪费**：大 $J$ 时尽快切换，别“枚举一堆又丢弃”
- **贴近真实最优**：随实现与数据分布自适应

------

## 4.1 内存硬约束（必须满足）

设可用于 `AllPairs` 的内存预算为 `MemBudgetPairs`（字节），每个 pair 存储开销为 `sizeof(Pair)`（两个 32-bit id 也至少 8B，但在 C++/Java 容器中往往 16B~24B 甚至更高）。则：
$$
J_\star^{\text{mem}}=
\left\lfloor\frac{\text{MemBudgetPairs}}{\text{sizeof(Pair)}}\right\rfloor.
$$
必须取
$$
J_\star \le J_\star^{\text{mem}}.
$$
**工程建议：** 如果用 `vector<pair<int,int>>`，请把 `sizeof(Pair)` 用实际编译环境测出来（含对齐）。

------

## 4.2 大分支时间上界约束（保证切换后不“拖累两遍采样”）

切换分支下，Phase1 最多额外枚举 $O(J_\star)$ 个 pair（超过阈值立刻停止 Report）。因此
$$
T_{\text{Adaptive,switch}} \approx T_{\text{Sampling}} + O(J_\star).
$$
实用建议是让 $J_\star$ 不超过“两遍扫描的主要几何成本量级”，例如：
$$
J_\star^{\text{time}} = C_1\cdot(\text{一次 sweep 的平均查询成本}) + C_2\cdot t.
$$
最终
$$
J_\star = \min(J_\star^{\text{mem}},\ J_\star^{\text{time}}).
$$

------

## 4.3 交叉点标定（推荐：可复现实验写法）

做 benchmark（不同数据分布 / 不同 $n,t$）：

1. 运行纯 Enumerate+Sampling（版本 A），记录耗时与峰值内存
2. 运行纯 Sampling（版本 B），记录耗时与峰值内存
3. 找到“耗时相近”的交叉点 $|J_{\text{cross}}|$
4. 设置

$$
J_\star \approx 0.8\cdot |J_{\text{cross}}|
$$

并再受 $J_\star^{\text{mem}}$ 截断

这样 $J_\star$ 会贴合你实现与数据分布，不会拍脑袋。

------

# 5. 算法分析（正确性 + 复杂度；三个版本都包含）

## 5.1 正确性

### 5.1.1 引理：START 时维度分离成立

已在 §1.4 给出（引理 1）。因此在 START(q) 时：
$$
q\cap r\neq\varnothing\ \Longleftrightarrow\ q^\star\cap r^\star\neq\varnothing
\quad(\forall r\in\mathrm{Active}(\text{Other})).
$$

------

### 5.1.2 引理：事件分块 $J=\biguplus_e J_e$

已在 §1.5 给出（命题 1）。关键依赖：
 **半开 + END-before-START + START 全序 tie-break**。

------

### 5.1.3 引理：R-tree 的 Report / Count 精确性

- `ReportIntersect(Q)`：内部节点仅用 `MBR` 做剪枝，不会漏掉任何真实命中；叶层用半开精确相交过滤，因此输出集合恰为 $K(Q)$。
- `CountIntersect(Q)`：
  - 不相交剪枝贡献 0 正确
  - `MBR ⊆ Q` 时整子树必命中，贡献 `size(node)` 正确
  - 其余递归到叶并精确过滤，计数正确

------

### 5.1.4 引理：SampleIntersect(Q,t′) 的 i.i.d. 均匀性

记命中集合 $K(Q)$。Frontier 构造得到：
$$
K(Q)=\biguplus_{C\in\mathcal{F}}S_C,\quad |S_C|=w_C,\quad \sum_C w_C=|K(Q)|.
$$
算法对单个样本的生成过程等价于：

1. 抽组件 $C$ ：$\Pr(C)=w_C/|K(Q)|$
2. 在 $S_C$ 内均匀抽一个元素：$\Pr(x\mid C)=1/w_C$

于是对任意 $x\in K(Q)$（属于唯一组件 $C(x)$）：
$$
\Pr(\text{输出 }x)=\frac{w_{C(x)}}{|K(Q)|}\cdot\frac{1}{w_{C(x)}}=\frac{1}{|K(Q)|}.
$$
重复抽样使用独立随机数，因此得到 i.i.d. 均匀有放回。∎

------

## 5.2 版本 A：Enumerate+Sampling 的正确性

- 枚举阶段：由命题 1，每个相交对只在“更靠后 START”处被加入一次，因此 `Pairs` 精确等于 $J$
- 采样阶段：对 `Pairs[0..|J|-1]` 独立均匀抽下标，显然得到 $J$ 上 i.i.d. 均匀

------

## 5.3 版本 B：Sampling 的正确性

设 Phase1 得到每个 START 权重 $w_e=|J_e|$，并有 $W=\sum_e w_e=|J|$。

对任意 $P\in J$，存在唯一 $e^\star$ 使 $P\in J_{e^\star}$。Phase2 选中事件的概率：
$$
\Pr(E_j=e^\star)=\frac{w_{e^\star}}{W}.
$$
Phase3 在给定事件 $e^\star$ 下，对 $J_{e^\star}$ 内均匀：
$$
\Pr(P\mid E_j=e^\star)=\frac{1}{w_{e^\star}}.
$$
因此
$$
\Pr(Z_j=P)=\frac{w_{e^\star}}{W}\cdot\frac{1}{w_{e^\star}}=\frac{1}{|J|}.
$$
独立性来自：事件抽样对不同 $j$ 独立；每次局部采样也是有放回独立；slot 回填只是重排不引入相关性。∎

------

## 5.4 版本 C：Adaptive+Sampling 的正确性

- 未切换：`AllPairs` 完整枚举 $J$，等价版本 A
- 已切换：Phase1 虽丢弃枚举结果，但**保留精确 $w_e$ 与 $W$**，后续完全等价版本 B
   切换不改变事件顺序、插删时刻与计数，因此不改变分布。∎

------

## 5.5 复杂度分析（参数化 + 最坏退化）

R-tree 性能依赖 MBR 重叠。我们用参数化刻画一次查询代价：

对一次查询 $Q$：

- $V(Q)$：访问的节点数（含内部与叶）
- $E(Q)$：检查的叶条目数（叶内逐条精确过滤）
- 树高 $h\approx O(\log_B n)$，$B$ 为扇出

则：

- `ReportIntersect(Q)`：$O(V(Q)+E(Q)+k)$，$k=|K(Q)|$

- `CountIntersect(Q)`：$O(V(Q)+E(Q))$

- `SampleIntersect(Q,t′)`：

  - BuildFrontier：$O(V(Q)+E(Q))$
  - 子树均匀下沉：最朴素实现 $O(t′\cdot h)$
     总计：

  $$
  O(V(Q)+E(Q)+t′h).
  $$

------

### 5.5.1 版本 A：Enumerate+Sampling

- 排序：$O(n\log n)$

- sweep 更新：$O(\sum \text{Insert/Delete})$（通常 $\approx O(n\log_B n)$，最坏可退化）

- 每个 START 一次 Report：总代价
  $$
  O\Big(\sum_{e\in E}(V_e+E_e)\Big)+O(|J|).
  $$

- 数组采样：$O(t)$

空间：
$$
S=O(n+|J|+t).
$$

------

### 5.5.2 版本 B：Sampling

- Phase1（Count）：
  $$
  T_1=O(n\log n)+O(\text{updates})+\sum_{e\in E}O(V_e+E_e).
  $$

- Phase2：alias + slots：$O(n+t)$

- Phase3：对所有 $t_e>0$ 的事件执行 SampleIntersect：
  $$
  T_3=O(\text{updates})+\sum_{e:t_e>0}O(V_e+E_e+t_e h)
  $$

且 $\sum_e t_e=t$，故总下沉项为 $O(t h)$。

空间：
$$
S=O(n+t)
$$
（外加 `w[]` 与 alias 的 $O(n)$）。

------

### 5.5.3 版本 C：Adaptive+Sampling

- 若未切换（$|J|\le J_\star$）：等价版本 A，但 $|J|$ 被阈值上界限制

- 若切换（$|J|>J_\star$）：
  $$
  T_{\text{switch}} \approx T_{\text{Sampling}} + O(J_\star),\qquad
  S_{\text{switch}}=O(n+\max(t,J_\star)).
  $$

------

### 5.5.4 最坏情况警告（务必保留）

在大量 MBR 重叠、剪枝几乎失效时，可能出现 $V(Q)=\Theta(n)$，导致整体退化到 $\Theta(n^2)$ 级别。这是 R-tree baseline 的固有弱点，应在实验分析中明确说明。

------

## 你可以直接用的“增强点清单”（确保 baseline 真能到 90+）

1. **事件排序写成严格全序**（含 box-id tie-break）
2. **明确活跃集语义与 Insert/Delete 时刻**（和实现一致）
3. **CountIntersect 必须有 `MBR ⊆ Q` 快捷分支 + `size(node)`**
4. **SampleIntersect 用 Frontier 不交分解**，并在 FullyAccepted 子树用 size 比例下沉保证均匀
5. **0 权重 START 事件从 alias 中剔除**，避免概率/实现细节坑
6. **Delete 的 handle/FindLeaf 两方案都写清楚**（给读者可复现路径）
7. 正确性证明至少写到：
   - 维度分离引理
   - 唯一归属命题
   - SampleIntersect 的均匀性证明
   - Sampling 两遍的全局均匀性推导
