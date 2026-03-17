#!/usr/bin/env python3
"""Plot main result figures for the v2 light alpha100 sweep."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.ticker import FuncFormatter


OURS_COLOR = "#0B6E4F"
KD_COLOR = "#C44E52"
DIM_COLORS = {
    "d2": "#0B6E4F",
    "d3": "#1F77B4",
    "d4": "#E17C05",
    "d5": "#7A3E9D",
}
METHOD_LABELS = {
    "ours": "SJS",
    "kd_tree": "KDS Baseline",
}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    default_summary = root / "out" / "alpha100_light" / "all_dims_summary.csv"
    default_manifest = root / "out" / "alpha100_light" / "manifest.tsv"
    default_out = root / "out" / "alpha100_light" / "figures"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=default_summary)
    parser.add_argument("--manifest", type=Path, default=default_manifest)
    parser.add_argument("--out_dir", type=Path, default=default_out)
    return parser.parse_args()


def load_dim_map(manifest_path: Path) -> dict[int, str]:
    dim_map: dict[int, str] = {}
    with manifest_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            dim_map[int(row["n"])] = f"d{int(row['dim'])}"
    if not dim_map:
        raise ValueError(f"No rows found in manifest: {manifest_path}")
    return dim_map


def load_rows(summary_path: Path, dim_map: dict[int, str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with summary_path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            n_r = int(row["n_r"])
            if n_r not in dim_map:
                raise ValueError(f"Cannot map n_r={n_r} to a dimension from manifest")
            rows.append(
                {
                    "dim": dim_map[n_r],
                    "n_r": n_r,
                    "n_s": int(row["n_s"]),
                    "t": int(row["t"]),
                    "method": row["method"],
                    "wall_mean_ms": float(row["wall_mean_ms"]),
                    "wall_stdev_ms": float(row["wall_stdev_ms"]),
                    "count_mean_key": row["count_mean"],
                    "ok_rate": float(row["ok_rate"]),
                }
            )
    if not rows:
        raise ValueError(f"No rows found in summary CSV: {summary_path}")
    return rows


def geometric_mean(values: list[float]) -> float:
    if not values:
        raise ValueError("geometric_mean requires at least one value")
    if any(v <= 0 for v in values):
        raise ValueError("geometric_mean requires positive values")
    return math.exp(sum(math.log(v) for v in values) / len(values))


def format_speedup(value: float) -> str:
    if value >= 1:
        return f"{value:.2f}x"
    if value >= 0.1:
        return f"{value:.2f}x"
    if value >= 0.01:
        return f"{value:.3f}x"
    return f"{value:.4f}x"


def build_aggregates(
    rows: list[dict[str, object]],
) -> tuple[
    dict[str, dict[str, dict[int, list[float]]]],
    dict[str, dict[int, dict[str, float]]],
    list[dict[str, object]],
]:
    runtime: dict[str, dict[str, dict[int, list[float]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )
    pair_buckets: dict[tuple[str, int, str], dict[str, float]] = defaultdict(dict)

    for row in rows:
        dim = str(row["dim"])
        method = str(row["method"])
        t_value = int(row["t"])
        wall = float(row["wall_mean_ms"])
        runtime[dim][method][t_value].append(wall)

        pair_key = (dim, t_value, str(row["count_mean_key"]))
        pair_buckets[pair_key][method] = wall

    speedup_pairs: list[dict[str, object]] = []
    for (dim, t_value, count_key), pair in sorted(pair_buckets.items()):
        if "ours" not in pair or "kd_tree" not in pair:
            continue
        speedup = pair["kd_tree"] / pair["ours"]
        speedup_pairs.append(
            {
                "dim": dim,
                "t": t_value,
                "count_key": count_key,
                "ours_wall_ms": pair["ours"],
                "kd_wall_ms": pair["kd_tree"],
                "speedup": speedup,
            }
        )

    speedups: dict[str, dict[int, dict[str, float]]] = defaultdict(dict)
    grouped_speedups: dict[tuple[str, int], list[float]] = defaultdict(list)
    for pair in speedup_pairs:
        grouped_speedups[(str(pair["dim"]), int(pair["t"]))].append(float(pair["speedup"]))

    for (dim, t_value), values in sorted(grouped_speedups.items()):
        speedups[dim][t_value] = {
            "geomean": geometric_mean(values),
            "min": min(values),
            "max": max(values),
        }

    return runtime, speedups, speedup_pairs


def comma_formatter() -> FuncFormatter:
    return FuncFormatter(lambda x, _: f"{int(x):,}" if x >= 1 else f"{x:g}")


def save_figure(fig: plt.Figure, base_path: Path) -> None:
    fig.savefig(base_path.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(base_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_runtime_by_dim(
    out_dir: Path,
    runtime: dict[str, dict[str, dict[int, list[float]]]],
    dim_to_n: dict[str, int],
) -> None:
    dims = sorted(runtime.keys(), key=lambda d: int(d[1:]))
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2), sharex=True, sharey=True)
    axes = axes.ravel()

    legend_handles = None
    legend_labels = None

    for ax, dim in zip(axes, dims):
        for method, color in [("ours", OURS_COLOR), ("kd_tree", KD_COLOR)]:
            t_values = sorted(runtime[dim][method].keys())
            y_mean = [mean(runtime[dim][method][t]) for t in t_values]
            y_low = [min(runtime[dim][method][t]) for t in t_values]
            y_high = [max(runtime[dim][method][t]) for t in t_values]

            ax.fill_between(t_values, y_low, y_high, color=color, alpha=0.12)
            line = ax.plot(
                t_values,
                y_mean,
                marker="o",
                markersize=6,
                linewidth=2.4,
                color=color,
                label=METHOD_LABELS[method],
            )

        ax.set_title(f"{dim.upper()}  (n = {dim_to_n[dim]:,})", fontsize=12, weight="bold")
        ax.set_yscale("log")
        ax.grid(True, which="major", linestyle="--", alpha=0.32)
        ax.grid(True, which="minor", linestyle=":", alpha=0.12)
        ax.xaxis.set_major_formatter(comma_formatter())
        if legend_handles is None:
            legend_handles, legend_labels = ax.get_legend_handles_labels()

    for ax in axes[:]:
        ax.set_xlabel("sample size t")
        ax.set_ylabel("wall time (ms)")

    fig.suptitle(
        "v2 main result: wall time vs sample size",
        fontsize=16,
        weight="bold",
        y=0.98,
    )
    fig.text(
        0.5,
        0.01,
        "Solid lines show the mean over generator seeds; shaded bands show min-max across seeds.",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    if legend_handles and legend_labels:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            ncol=2,
            frameon=False,
            bbox_to_anchor=(0.5, 0.94),
        )
    fig.tight_layout(rect=(0, 0.04, 1, 0.92))
    save_figure(fig, out_dir / "main_runtime_by_dim")


def plot_speedup_vs_t(
    out_dir: Path,
    speedups: dict[str, dict[int, dict[str, float]]],
) -> None:
    dims = sorted(speedups.keys(), key=lambda d: int(d[1:]))
    fig, ax = plt.subplots(figsize=(9.5, 5.8))

    for dim in dims:
        t_values = sorted(speedups[dim].keys())
        y_values = [speedups[dim][t]["geomean"] for t in t_values]
        y_low = [speedups[dim][t]["min"] for t in t_values]
        y_high = [speedups[dim][t]["max"] for t in t_values]
        color = DIM_COLORS.get(dim, "#333333")

        ax.fill_between(t_values, y_low, y_high, color=color, alpha=0.10)
        ax.plot(
            t_values,
            y_values,
            marker="o",
            markersize=6,
            linewidth=2.5,
            color=color,
            label=dim.upper(),
        )

    ax.axhline(1.0, color="#444444", linewidth=1.2, linestyle="--")
    ax.set_xscale("linear")
    ax.set_yscale("log")
    ax.set_xlabel("sample size t")
    ax.set_ylabel("speedup = KDS / SJS")
    ax.set_title("v2 main result: relative speedup by dimension", fontsize=15, weight="bold")
    ax.xaxis.set_major_formatter(comma_formatter())
    ax.grid(True, which="major", linestyle="--", alpha=0.32)
    ax.grid(True, which="minor", linestyle=":", alpha=0.12)
    ax.legend(frameon=False, ncol=4, loc="upper center")
    fig.text(
        0.5,
        0.01,
        "Values above 1 mean SJS is faster. Bands show the min-max spread across generator seeds.",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save_figure(fig, out_dir / "main_speedup_vs_t")


def plot_dim_geomean_bar(
    out_dir: Path,
    speedups: dict[str, dict[int, dict[str, float]]],
) -> None:
    dims = sorted(speedups.keys(), key=lambda d: int(d[1:]))
    values = [geometric_mean([speedups[dim][t]["geomean"] for t in sorted(speedups[dim])]) for dim in dims]
    colors = [DIM_COLORS.get(dim, "#333333") for dim in dims]

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    bars = ax.bar([dim.upper() for dim in dims], values, color=colors, width=0.64)
    ax.axhline(1.0, color="#444444", linewidth=1.2, linestyle="--")
    ax.set_yscale("log")
    ax.set_ylabel("geometric mean speedup (KDS / SJS)")
    ax.set_title("v2 main result: overall speedup by dimension", fontsize=15, weight="bold")
    ax.grid(True, which="major", axis="y", linestyle="--", alpha=0.32)
    ax.grid(True, which="minor", axis="y", linestyle=":", alpha=0.12)

    for bar, value in zip(bars, values):
        label = format_speedup(value)
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value * (1.12 if value >= 1 else 1.28),
            label,
            ha="center",
            va="bottom",
            fontsize=10,
            color="#222222",
        )

    fig.text(
        0.5,
        0.01,
        "Each bar aggregates all sample sizes and both generator seeds for that dimension.",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save_figure(fig, out_dir / "main_dim_geomean_speedup")


def plot_speedup_heatmap(
    out_dir: Path,
    speedups: dict[str, dict[int, dict[str, float]]],
) -> None:
    dims = sorted(speedups.keys(), key=lambda d: int(d[1:]))
    t_values = sorted({t for dim in dims for t in speedups[dim]})
    matrix = []
    for dim in dims:
        matrix.append([math.log10(speedups[dim][t]["geomean"]) for t in t_values])

    vmax = max(abs(v) for row in matrix for v in row)
    cmap = LinearSegmentedColormap.from_list(
        "speedup_balance",
        ["#B83B3B", "#F6E7C1", "#0E7490"],
        N=256,
    )
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(8.6, 4.7))
    image = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(range(len(t_values)))
    ax.set_xticklabels([f"{t:,}" for t in t_values])
    ax.set_yticks(range(len(dims)))
    ax.set_yticklabels([dim.upper() for dim in dims])
    ax.set_xlabel("sample size t")
    ax.set_ylabel("dimension")
    ax.set_title("v2 main result: speedup heatmap", fontsize=15, weight="bold")

    for row_idx, dim in enumerate(dims):
        for col_idx, t_value in enumerate(t_values):
            speedup = speedups[dim][t_value]["geomean"]
            text_color = "white" if abs(matrix[row_idx][col_idx]) > vmax * 0.55 else "#222222"
            ax.text(
                col_idx,
                row_idx,
                format_speedup(speedup),
                ha="center",
                va="center",
                fontsize=10,
                color=text_color,
                weight="bold",
            )

    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label("log10(KDS / SJS)")
    fig.text(
        0.5,
        0.01,
        "Blue cells favor SJS; red cells favor the KDS baseline.",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save_figure(fig, out_dir / "main_speedup_heatmap")


def write_plot_data(
    out_dir: Path,
    runtime: dict[str, dict[str, dict[int, list[float]]]],
    speedups: dict[str, dict[int, dict[str, float]]],
    dim_to_n: dict[str, int],
) -> None:
    out_path = out_dir / "main_plot_data.csv"
    fieldnames = [
        "dim",
        "n_r",
        "t",
        "ours_mean_ms",
        "ours_min_ms",
        "ours_max_ms",
        "kd_mean_ms",
        "kd_min_ms",
        "kd_max_ms",
        "speedup_geomean",
        "speedup_min",
        "speedup_max",
    ]

    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for dim in sorted(runtime.keys(), key=lambda d: int(d[1:])):
            for t_value in sorted(runtime[dim]["ours"].keys()):
                ours_values = runtime[dim]["ours"][t_value]
                kd_values = runtime[dim]["kd_tree"][t_value]
                writer.writerow(
                    {
                        "dim": dim,
                        "n_r": dim_to_n[dim],
                        "t": t_value,
                        "ours_mean_ms": f"{mean(ours_values):.6f}",
                        "ours_min_ms": f"{min(ours_values):.6f}",
                        "ours_max_ms": f"{max(ours_values):.6f}",
                        "kd_mean_ms": f"{mean(kd_values):.6f}",
                        "kd_min_ms": f"{min(kd_values):.6f}",
                        "kd_max_ms": f"{max(kd_values):.6f}",
                        "speedup_geomean": f"{speedups[dim][t_value]['geomean']:.6f}",
                        "speedup_min": f"{speedups[dim][t_value]['min']:.6f}",
                        "speedup_max": f"{speedups[dim][t_value]['max']:.6f}",
                    }
                )


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    dim_map_by_n = load_dim_map(args.manifest)
    rows = load_rows(args.summary, dim_map_by_n)
    runtime, speedups, _ = build_aggregates(rows)

    dim_to_n = {dim: n_r for n_r, dim in dim_map_by_n.items()}
    plot_runtime_by_dim(args.out_dir, runtime, dim_to_n)
    plot_speedup_vs_t(args.out_dir, speedups)
    plot_dim_geomean_bar(args.out_dir, speedups)
    plot_speedup_heatmap(args.out_dir, speedups)
    write_plot_data(args.out_dir, runtime, speedups, dim_to_n)

    print(f"[plot_main_results] wrote figures to: {args.out_dir}")


if __name__ == "__main__":
    main()
