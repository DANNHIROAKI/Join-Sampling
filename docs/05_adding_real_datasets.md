# 05_adding_real_datasets.md

目前你们的实验主要基于合成数据；但为了让论文更“SIGMOD 级别可信”，后续很大概率需要：

- 至少 1~3 个真实空间数据集（OSM / TIGER / 城市 POI / 建筑物等）
- 对比在真实数据上的速度/内存/误差
- 与合成数据的趋势做 cross-check

本仓库已经预留了接入真实数据的工程位置（stub），你只需要按本文档补全即可。

---

## 1) 推荐策略：统一转成内部 .bin

真实数据的源格式五花八门（pbf/shp/geojson/wkt/csv），而 baseline 的输入统一是：

- `Dataset<Dim>`（两关系 R/S 的 MBR 集合）

因此最稳定的策略是：

> **真实数据只做一次离线转换**：原始格式 → MBR（矩形）→ 内部 `.bin`  
> 后续所有实验都用 `--dataset_source=binary` 直接读 `.bin`，避免每次解析真实格式影响 benchmark。

---

## 2) 需要补全的代码入口

### 2.1 realdata stub

- `include/sjs/io/realdata_stub.h`：声明接口
- `src/io/realdata_stub.cpp`：当前是占位实现

你需要实现：

- `realdata::LoadRelation(Source src, path, Relation<2>* out)`
- （可选）`realdata::LoadDatasetPair(...)`

并支持至少一个真实格式：
- OSM PBF（需要 libosmium）
- 或 Shapefile（需要 GDAL/OGR）
- 或 GeoJSON / WKT CSV（可用轻量解析器）

### 2.2 转换工具

- `apps/sjs_convert_realdata.cpp`：预留的 CLI 工具  
  目标：把真实数据转成内部 `.bin`。

推荐输出结构：

```
data/real/<dataset_name>/
  R.bin
  S.bin
  meta.json
```

---

## 3) 真实数据“几何到矩形”的口径

真实数据往往不是矩形，而是点/线/多边形。

统一规则：

- 如果输入是几何体（polygon/line），先计算其 MBR（min/max x,y）
- 得到矩形后，转换为 half-open：
  - `hi` 使用 `nextafter(hi, +inf)` 或加一个极小 epsilon（工程上更简单：保持 closed 表示也可，但需要保证整个系统一致）
  - 本项目默认 half-open；建议最终写入 `.bin` 时标记 `kHalfOpen`

**边界处理很关键**：否则审稿人很容易追问“贴边是否算相交”。

---

## 4) 坐标系统与归一化

真实数据可能在：

- 经纬度（WGS84）
- Web Mercator
- 本地投影（UTM 等）
- 甚至是整数网格

建议统一：

1) 转到同一坐标系（至少 R 和 S 必须一致）
2) 可选：归一化到 \([0,1]^2\) 以便不同数据集可比较（不是必须，但利于控制数值范围）

归一化做法：

- 先算全局 bounding box（覆盖 R 和 S）
- 对每个坐标做线性缩放到 \([0,1]\)

---

## 5) 去噪与过滤（强烈建议）

真实数据常见问题：

- 空几何 / 无效几何 → 生成空 MBR
- 极小/极大尺度差异 → baseline 的行为会很极端
- 数据量过大 → 你可能需要抽样或裁剪

建议在转换阶段支持：

- `--drop_empty=1`
- `--limit=<N>`：最多保留 N 个对象
- `--clip=<xmin,ymin,xmax,ymax>`：只保留与窗口相交的对象（城市区域）
- `--shuffle=1`：打乱输入顺序（避免数据源的顺序偏差）

这些参数可以作为 `sjs_convert_realdata` 的 CLI 选项。

---

## 6) 如何选择真实数据集（论文角度）

建议至少覆盖两类：

1) **相对稀疏**：例如道路 vs POI（\(|J|\) 较小）
2) **相对稠密**：例如建筑 vs 地块（\(|J|\) 较大）

这样能展示：

- Enum+Sampling 在小 join 的优势
- Sampling/Adaptive 在大 join 的必要性

---

## 7) 最小落地计划（两周内可完成）

- 第 1 步：先支持一种最容易解析的真实格式（例如 GeoJSON 或 WKT CSV）
- 第 2 步：用 `sjs_convert_realdata` 转出 `.bin`
- 第 3 步：用 `sjs_run --dataset_source=binary` 跑 2~3 个 baseline
- 第 4 步：补一张 “synthetic vs real” 的趋势对比图（scale 或 alpha 类似）

做到以上，你的论文实验部分就能“更像 SIGMOD”。
