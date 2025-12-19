
# SpatialJoinSamplingCpp

A lightweight **C++17** experiment framework for **counting** and drawing **i.i.d. uniform samples (with replacement)** from the result of a **2D spatial join** between two relations of **axis-aligned rectangles**.

This project is designed for SIGMOD-style evaluations: multiple baselines, repeatable runs, sweep harness, and simple output files (CSV/TSV + custom binary).

## Key points

- Geometry uses **half-open boxes**: `[lo, hi)` on each axis.
- Current build focuses on **Dim = 2** (most code is templated for future extension).
- No heavy dependencies (no external JSON library; `sjs_sweep` parses a small JSON subset).

## Build

```bash
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
````

Executables (built by default):

* `sjs_gen_dataset` – generate synthetic datasets and export to binary/CSV
* `sjs_run` – run a single experiment (with repeats)
* `sjs_sweep` – run a parameter sweep from a JSON config
* `sjs_verify` – brute-force oracle checks + sampling-quality tests (small datasets)
* `sjs_profile` – profiling / microbenchmarks
* `sjs_convert_realdata` – conversion entry point (real-data loaders are stubs)

## Quick start

### 1) Generate a synthetic dataset (optional)

```bash
./sjs_gen_dataset \
  --dataset_source=synthetic --gen=stripe --dataset=demo \
  --n_r=100000 --n_s=100000 --alpha=1e-6 --gen_seed=1 \
  --out_dir=data/synthetic --write_csv=1
```

Outputs (default names):

* `data/synthetic/demo_R.bin`
* `data/synthetic/demo_S.bin`
* (optional) `demo_R.csv`, `demo_S.csv`

### 2) Run a single experiment

Synthetic (generate on the fly):

```bash
./sjs_run \
  --dataset_source=synthetic --gen=stripe --dataset=demo \
  --n_r=100000 --n_s=100000 --alpha=1e-6 --gen_seed=1 \
  --method=ours --variant=sampling --t=100000 --seed=1 --repeats=3 \
  --out_dir=results/raw/demo
```

Binary input:

```bash
./sjs_run \
  --dataset_source=binary --dataset=demo \
  --path_r=data/synthetic/demo_R.bin --path_s=data/synthetic/demo_S.bin \
  --method=ours --variant=sampling --t=100000 --seed=1 --repeats=3 \
  --out_dir=results/raw/demo_bin
```

Help:

```bash
./sjs_run --help
```

### 3) Run a sweep

```bash
./sjs_sweep --config=config/sweeps/sweep_alpha.json
```

### 4) Verify correctness / sampling quality (small datasets)

```bash
./sjs_verify \
  --dataset_source=synthetic --gen=stripe --dataset=verify \
  --n_r=2000 --n_s=2000 --alpha=1e-3 --gen_seed=1 \
  --method=ours --variant=sampling --t=20000 --seed=1 --repeats=1
```

## Methods and variants

### Variants

* `sampling` – 2-pass sampling pipeline (target: exact i.i.d. uniform samples)
* `enum_sampling` – enumerate join stream then rank-sample
* `adaptive` – enumerate when join is small, otherwise fallback to sampling

### Methods (Dim=2)

* `ours`
* `aabb`
* `interval_tree`
* `kd_tree`
* `r_tree`
* `range_tree`
* `pbsm`
* `tlsop`
* `sirs`
* `rejection`
* `tsunami`

Run `./sjs_run --help` to see the supported `method/variant` combinations.

## Datasets

### Synthetic generators

Pass via `--gen=...`:

* `stripe` (alias of `stripe_ctrl_alpha`)
* `uniform`
* `clustered`
* `hetero_sizes`

### Binary format (recommended)

Custom, versioned binary format (magic `SJSBOX`). Used for fast, deterministic IO.

### CSV input

A simple reader is provided (intended for debugging/small data). Separator can be set via `--csv_sep=tab` or `--csv_sep=,`.

## Outputs

### `sjs_run`

* `<out_dir>/run.csv` (one row per repeat)
* Optional samples: `<out_dir>/samples/*.tsv` (enable with `--write_samples=1`)

### `sjs_sweep`

* `<out_dir>/sweep_raw.csv` (one row per run/repeat)
* `<out_dir>/sweep_summary.csv` (aggregated per configuration)

## Testing

```bash
cd build
ctest --output-on-failure
```

## Known limitations

* Current build targets **2D** (`--dim=2`).
* `sjs_convert_realdata` is wired, but real-data importers are **stubs** (see `include/sjs/io/realdata_stub.h`).
