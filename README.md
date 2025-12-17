# SpatialJoinSamplingCpp

一个**全 C++**的实验框架，用于评估“在空间 Join 结果集合 \(J\) 上进行 **i.i.d. 均匀（有放回）抽样**”的多种方法（Ours + 多个 Baseline），并提供可复现的 **SIGMOD 风格实验流水线**（单跑 / 扫参 / 校验 / 性能剖析 / 数据生成与格式转换）。

> 当前默认只编译 **二维（Dim=2）** 版本，但所有核心数据结构与接口均以模板设计，方便后续扩展到高维与真实数据集。

---

## 目录结构

```
SpatialJoinSamplingCpp/
├── CMakeLists.txt
├── README.md
│
├── include/                   # 头文件（大多 header-only）
│   └── sjs/...
├── src/                       # .cpp 编译单元（工厂/显式实例化/封装）
│   └── ...
│
├── apps/                      # 全部可执行程序入口（无 Python）
│   ├── sjs_gen_dataset.cpp
│   ├── sjs_run.cpp
│   ├── sjs_sweep.cpp
│   ├── sjs_verify.cpp
│   ├── sjs_profile.cpp
│   └── sjs_convert_realdata.cpp
│
├── tests/                     # 单元测试/冒烟测试（无第三方框架）
├── config/                    # 默认配置与 sweep 配置样例
└── docs/                      # 项目文档（实验协议/数据格式/扩展指南）
```

---

## 1) 编译

### 依赖

- CMake >= 3.16
- C++17 编译器（GCC/Clang/MSVC）
- `Threads`（CMake 会自动找）

### 构建命令

```bash
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
```

编译后生成的可执行文件在 `build/` 目录下（例如 `build/sjs_run`）。

---

## 2) 快速开始（合成数据）

### 2.1 生成合成数据（并导出内部二进制格式）

```bash
./sjs_gen_dataset \
  --dataset_source=synthetic \
  --gen=stripe_ctrl_alpha \
  --dataset=demo \
  --n_r=100000 --n_s=100000 \
  --alpha=1e-6 --gen_seed=1 \
  --out_dir=data/synthetic/demo \
  --write_csv=0
```

生成输出：
- `data/synthetic/demo/demo_R.bin`
- `data/synthetic/demo/demo_S.bin`
- `data/synthetic/demo/demo_gen_report.json`（轻量 JSON）

> 二进制格式说明见 `docs/01_data_format.md`。

### 2.2 单次运行：指定方法 + 版本 + 抽样数 t

```bash
./sjs_run \
  --dataset_source=binary \
  --path_r=data/synthetic/demo/demo_R.bin \
  --path_s=data/synthetic/demo/demo_S.bin \
  --dataset=demo \
  --method=ours \
  --variant=sampling \
  --t=100000 \
  --seed=1 \
  --repeats=3 \
  --out_dir=results/raw/demo
```

输出 CSV（每次 repeat 一行），并带 phase breakdown / 估计误差等字段。

### 2.3 扫参：读取 `config/sweeps/*.json` 批量跑

```bash
./sjs_sweep --config=config/sweeps/sweep_alpha.json
```

典型输出：
- `.../sweep_raw.csv`（每个 repeat 一行）
- `.../sweep_summary.csv`（按 key 聚合统计）

### 2.4 小规模正确性验证 + 采样质量检验

```bash
./sjs_verify \
  --dataset_source=synthetic \
  --gen=stripe_ctrl_alpha \
  --dataset=verify_small \
  --n_r=2000 --n_s=2000 \
  --alpha=1e-4 --gen_seed=1 \
  --method=ours --variant=sampling \
  --t=20000 --seed=1 \
  --out_dir=results/verify
```

### 2.5 剖析：输出 phase breakdown / 计数器

```bash
./sjs_profile \
  --dataset_source=synthetic --gen=stripe_ctrl_alpha --dataset=prof \
  --n_r=200000 --n_s=200000 --alpha=1e-6 --gen_seed=1 \
  --method=ours --variant=adaptive --t=100000 --seed=1 \
  --out_dir=results/profile
```

---

## 3) 方法（method）与三种版本（variant）

### method

项目内置方法名（`--method=`）包含：

- `ours`（本文方法）
- `interval_tree`
- `aabb`
- `kd_tree`
- `r_tree`
- `range_tree`
- `pbsm`
- `tlsop`
- `sirs`
- `rejection`
- `tsunami`（若已注册到工厂/CLI）

### variant

每个 baseline 必须实现三种版本（统一 Runner）：

- `sampling`：纯 Sampling（不显式枚举全部 join）
- `enum_sampling`：Enumerate + Sampling（全枚举 join，再数组采样）
- `adaptive`：Adaptive（阈值触发：小 join 用枚举，大 join 用 sampling）

详细解释见 `docs/03_baselines.md` 与 `docs/02_experiment_protocol.md`。

---

## 4) 复现性与输出约定

- **随机种子分离**：合成数据使用 `--gen_seed`，采样与算法随机性使用 `--seed`。
- **半开区间**：默认几何是 half-open box `[lo, hi)`；扫描线事件 tie-break 对应 half-open 语义。
- **输出字段**：`results/raw/*.csv` 里包含：
  - `count_value / count_ci_low / count_ci_high`
  - `wall_ms`
  - `used_enumeration / enum_truncated / adaptive_branch`
  - `phases_json`（各阶段耗时，JSON-lite）

---

## 5) 文档导航

- `docs/00_overview.md`：项目总览（算法/代码结构/工作流）
- `docs/01_data_format.md`：二进制/CSV 数据格式与字段
- `docs/02_experiment_protocol.md`：实验协议（公平性、指标、sweep 输出）
- `docs/03_baselines.md`：所有 baseline 的实现口径与可调参数
- `docs/04_extending_to_high_dim.md`：扩展到高维的工程路线
- `docs/05_adding_real_datasets.md`：接入真实数据的步骤（OSM/TIGER 等）

---

## 6) 运行测试

```bash
cd build
ctest --output-on-failure
```

---

## 7) FAQ

### Q: 为什么不直接用 JSON 库 / CLI 库？
为了让实验环境更“可移植、可控、少依赖”。当前解析的是 **最小可用 JSON 子集**（用于 sweep 文件），并保留 `extra` 字段以方便扩展。

### Q: 现在只支持二维吗？
是的，默认显式实例化的是 Dim=2。但数据结构/接口都是模板化的，文档中给了扩展到高维的路线（见 `docs/04_extending_to_high_dim.md`）。

