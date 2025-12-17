# 01_data_format.md

本项目支持三类数据源（`--dataset_source=`）：

- `synthetic`：运行时生成（便于 sweep）
- `binary`：项目自定义二进制格式（推荐用于大规模实验）
- `csv`：简单 CSV/TSV（便于调试/可视化）

---

## 1) 内存中的数据结构

核心容器在 `include/sjs/io/dataset.h`：

- `Relation<Dim>`：一组 `Box<Dim>` + 可选稳定 `ids`
- `Dataset<Dim>`：一对关系 `R` / `S` + 元信息（name、half_open）

### 1.1 Box 表示

`Box<Dim>` 由两个点组成：

- `lo`：下界点
- `hi`：上界点

默认 **half-open**：每维 `[lo[d], hi[d])`。

要求：

- `lo[d] < hi[d]`（proper box；否则为空/非法）
- 数据读入/生成时建议执行 `Validate(require_proper=true)`

---

## 2) 自定义二进制格式（.bin）

二进制读写在 `include/sjs/io/binary_io.h`。

### 2.1 设计目标

- 解析快、体积小（相比 CSV）
- 版本化头部，支持向前兼容
- 保存维度、标志位（是否 half-open，是否带 ids）
- 支持 float32 / float64 存储

### 2.2 文件布局（Little-Endian）

```
FileHeader (fixed)
name_len: u32
name bytes [name_len]
records[count]:
  lo[dim] (scalar)
  hi[dim] (scalar)
  (optional) id (u32)   if flags&kHasIds
```

`FileHeader` 关键字段：

- `magic = "SJSBOX\0\0"`
- `version = 1`
- `dim`：维度
- `scalar_bits`：32 或 64
- `flags`：bitset
  - `kHalfOpen`：几何语义是否 half-open
  - `kHasIds`：是否包含 id
- `count`：记录数
- `endian`：端序 marker（用于检测错误端序）

> 目前默认写出 `float64` 且 `kHasIds=true`。

### 2.3 Relation vs Dataset 的存储

- 一个 `.bin` 文件对应一个 `Relation`（即 R 或 S）
- `Dataset` 通常用两个文件表示：`*_R.bin` 与 `*_S.bin`

推荐命名：

```
<dataset>_R.bin
<dataset>_S.bin
```

---

## 3) CSV/TSV 格式（调试用）

CSV 读写在 `include/sjs/io/csv_io.h`。

### 3.1 推荐列格式

默认支持（带 header）：

- Dim=2：
  - `x_lo, y_lo, x_hi, y_hi, id`
- 若无 `id` 列，读入时可自动生成连续 id

分隔符：

- CSV：`,`
- TSV：`\t`（可用 `--csv_sep=tab` 或 `--csv_sep=\t`）

### 3.2 注意事项（非常重要）

- CSV 是文本格式：精度与解析速度受限，不建议在大规模 benchmark 中直接使用
- 若数据来自外部系统，请尽量先转成 `.bin`

---

## 4) 合成数据输出规范

`sjs_gen_dataset` 默认把合成数据写到：

- `data/synthetic/<dataset>/<dataset>_R.bin`
- `data/synthetic/<dataset>/<dataset>_S.bin`
- （可选）`*_R.csv` / `*_S.csv`
- `*_gen_report.json`：生成器参数与派生统计（轻量 JSON-lite，便于记录）

---

## 5) 真实数据接入建议

真实数据常见来源：

- OSM 建筑/道路等矢量边界
- TIGER/Line（美国行政区等）
- 自己的轨迹/地块矩形化数据

建议统一流程：

1. 读入原始数据（shp/pbf/geojson/…）
2. 转成内部 `Relation<Dim>`（半开、lo<hi、坐标归一或投影）
3. 写出 `.bin`（`WriteRelationBinary`）
4. 实验阶段直接 `--dataset_source=binary`

详细步骤见 `docs/05_adding_real_datasets.md`。
