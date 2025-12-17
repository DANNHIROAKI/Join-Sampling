# 04_extending_to_high_dim.md

当前仓库默认只显式实例化/编译 **Dim=2**（二维矩形）。但大部分模块是按 `template<int Dim>` 写的，因此扩展到高维主要是“补齐显式实例化 + 工厂 + 算法数据结构”。

下面给出一个**工程可落地**的扩展路线（按优先级排序）。

---

## 0) 先明确：两类“维度扩展”

1) **工程维度**：让程序在 `--dim=3/4/...` 时能编译 & 运行（不一定性能最优）。  
2) **算法维度**：让 Ours 的复杂度界与抽样正确性在高维也成立（这需要高维模式数据结构）。

建议先做 (1)，再做 (2)。

---

## 1) 工程维度：让 Dim>2 可以跑起来

### 1.1 数据结构层（基本已经支持）

- `geometry/point.h`, `geometry/box.h` 已是模板
- `io/dataset.h`, `binary_io.h`, `csv_io.h` 已是模板
- 合成数据生成器也已按 `Dim` 写接口（`GenerateDataset<Dim>`）

需要注意：

- 二进制格式 header 中有 `dim` 字段，但读/写 API 仍是编译期 Dim：  
  `ReadRelationBinary<Dim>(...)` 会检查 header 的 `dim == Dim`，否则报错。
- 因此要支持 runtime 的 `--dim`，你需要**在二进制读入前就决定 Dim**。

### 1.2 工厂层：为每个 Dim 增加 baseline_factory

当前只有：

- `src/baselines/baseline_factory_2d.{h,cpp}`
- `src/index/index_factory_2d.cpp`
- `src/data/synthetic/generator_factory_2d.cpp`

扩展方式（推荐）：

- 为每个维度增加一套工厂：
  - `baseline_factory_3d.*`
  - `index_factory_3d.cpp`
  - `generator_factory_3d.cpp`（如果需要）
- 或者写一个 `baseline_factory.h`，内部 `switch(dim)` 分发到具体工厂（便于 apps 调用）

### 1.3 显式实例化（减少模板膨胀）

当前目录里很多 `instantiations_2d.cpp`：

- `src/geometry/instantiations_2d.cpp`
- `src/index/instantiations_2d.cpp`
- `src/join/instantiations_2d.cpp`
- ...

扩展方式：

- 为 Dim=3 新增对应文件 `instantiations_3d.cpp`
- 内容是 `template class ...<3>;` 或 `template bool Foo<3>(...);`

> 如果你不显式实例化，纯模板也能编译，但会导致编译时间大幅增加（尤其 baseline 多、模板层级深时）。

### 1.4 Apps：去掉 Dim=2 的硬限制

目前 `apps/*.cpp` 里通常有：

```cpp
if (cfg.dataset.dim != 2) { ... return; }
```

扩展到多维后应改成：

- `switch(cfg.dataset.dim)`：
  - case 2：用 `BaselineFactory2D` + `Dataset<2>`
  - case 3：用 `BaselineFactory3D` + `Dataset<3>`
  - ...
- 对于暂不支持的 baseline / variant：在工厂里返回清晰错误信息（“not implemented for Dim=3”）

---

## 2) 算法维度：Ours 在高维的路线

### 2.1 从 2D 的 A/B 到 dD 的 2^(d-1) 模式

在 Dim=2 时，相交在非扫描维（y）可拆成两种互斥情况（A/B）。

在 Dim=d 时，选择第 1 维 sweep 后，剩余 \(d-1\) 维每一维都有 A/B 两种情况，因此模式数为：

\[
|\{A,B\}^{d-1}| = 2^{d-1}.
\]

实现上需要：

- 对每个 START 事件：
  - 计算所有模式的计数 \(w_e^{(g)}\)
  - Phase3 时按模式进行局部均匀采样

### 2.2 高维模式数据结构（核心工作量）

2D 版本用 BIT + buckets 就能做：

- `Count`：点刺/区间范围计数
- `Sample`：按 rank 选 bucket，再 bucket 内均匀取

高维需要一个支持：

- 动态 Insert/Delete（active set）
- Count / Sample / Report（可选）

的递归结构（例如“外层 segment tree + 内层子结构”的套路）。

如果你们论文里已经有“高维模式数据结构”设计（例如按 A/B 递归构造 stabbing/range 结构），工程上可以按以下步骤落地：

1. 先实现 Dim=3 的版本（最小非平凡高维）
2. 再模板化到 Dim=4/5（常数维）

建议 API 维持与 2D 相同：

```cpp
Count(q_projection)
Sample(q_projection, k)
Report(q_projection)   // optional
Insert(r)
Delete(r)
```

---

## 3) Baseline 在高维的优先级建议

如果目标是尽快做出“能跑的高维实验”，建议优先扩展：

1) `r_tree`：天然支持高维 MBR（很多实现就是 dD）  
2) `aabb`：AABB-tree 也天然支持 dD  
3) `kd_tree` / `range_tree`：若 embedding 已定义  
4) `interval_tree`：更偏 2D/1D 的结构，扩展需要重新设计“扫描轴与 interval 轴”的组合

---

## 4) 测试与验证（强烈建议）

- 新增 `tests/test_geometry_dim3.cpp`：随机盒子的相交判定
- 新增 `sjs_verify --dim=3` 小规模 oracle 对比（n=1e3 等）
- 对 2D 与 3D 都做采样质量检验（至少 χ²/KS 抽查）

---

## 5) 最小可行里程碑（MVP）

- ✅ Dim=3 的数据 I/O：能读写 `.bin`，能跑 `sjs_run --dim=3`  
- ✅ 至少 1 个 baseline（例如 `r_tree/enum_sampling`）能在 3D 跑通  
- ✅ `sjs_verify` 能在 3D 做 oracle 对比  
- ✅ 有一组 “Dim=2 vs Dim=3” 的 scale 实验曲线（哪怕 baseline 少）

做到以上，论文里就可以把“高维扩展”作为一个坚实的工程点来写。
