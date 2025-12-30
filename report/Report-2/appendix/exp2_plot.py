#!/usr/bin/env python3
"""exp2_plot.py

Paper-friendly plotting for EXP-2 (Runtime vs t).

This is a drop-in replacement of the EXP-2 plotter with two extra goals:
  1) Paper figures should use *consistent* repeats across points.
  2) The script should *not silently hide issues*: it writes a drop report if any
     (variant,method,t) points are excluded due to insufficient repeats.

Compared with the baseline plotter:
  - default t0=1000
  - new arg: --min_repeats (paper plots only include points with n>=min_repeats)
  - new arg: --paper_with_sample_phase (paper mode also writes sample-phase figure)
  - error_mode=auto uses --min_repeats (if set) to decide p95 vs stdev.

Inputs (under --out_dir):
  sweep_raw.csv
  sweep_summary.csv (optional)

Outputs (written into --out_dir):
  exp2_paper_runtime_vs_t.pdf/.png
  exp2_paper_delta_vs_t.pdf/.png
  exp2_paper_sample_phase_vs_t.pdf/.png  (optional)
  exp2_ns_per_sample.csv
  EXP2_README.txt
  exp2_dropped_points.csv  (if any)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple


# ----------------------------
# Small helpers
# ----------------------------
def to_int(x: str, default: int = 0) -> int:
    try:
        return int(float(x))
    except Exception:
        return default


def to_float(x: str, default: float = float("nan")) -> float:
    try:
        return float(x)
    except Exception:
        return default


def percentile_linear(vals: Sequence[float], p: float) -> float:
    """Linear interpolation percentile. p in [0,1]."""
    if not vals:
        return float("nan")
    xs = sorted(vals)
    if len(xs) == 1:
        return xs[0]
    pos = p * (len(xs) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return xs[lo]
    w = pos - lo
    return xs[lo] * (1 - w) + xs[hi] * w


def safe_stdev(vals: Sequence[float]) -> float:
    if len(vals) < 2:
        return 0.0
    return statistics.stdev(vals)


def ensure_file(path: str) -> None:
    if not os.path.isfile(path):
        raise FileNotFoundError(path)


def split_csv_list(s: str) -> List[str]:
    s = (s or "").strip()
    if not s:
        return []
    parts = []
    for p in s.split(","):
        p = p.strip()
        if p:
            parts.append(p)
    return parts


# ----------------------------
# CLI
# ----------------------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", required=True, help="EXP-2 output directory containing sweep_raw.csv")
    ap.add_argument("--t0", type=int, default=1000, help="Baseline t0 for Δruntime plot (default: 1000)")
    ap.add_argument(
        "--error",
        choices=["auto", "p95", "stdev"],
        default="auto",
        help="Error bar type (default: auto => p95 if repeats>=10 else stdev)",
    )
    ap.add_argument(
        "--warmup_reps",
        type=int,
        default=0,
        help="Exclude rep < warmup_reps from all plotted stats (default: 0).",
    )
    ap.add_argument(
        "--min_repeats",
        type=int,
        default=0,
        help="If >0, paper plots only include points with n>=min_repeats (after warmup filtering).",
    )
    ap.add_argument(
        "--mode",
        choices=["paper", "full"],
        default="paper",
        help="paper: concise figures; full: also output per-variant diagnostics.",
    )
    ap.add_argument(
        "--topk",
        type=int,
        default=6,
        help="In paper mode, plot only the top-K fastest methods (ranked on sampling @ max t), plus always_include.",
    )
    ap.add_argument(
        "--always_include",
        default="ours",
        help="Comma list of methods that must always be included in paper plots (default: ours)",
    )
    ap.add_argument(
        "--paper_methods",
        default="",
        help="If non-empty, use exactly this comma-list of methods for paper plots (overrides topk ranking).",
    )
    ap.add_argument(
        "--exclude_methods",
        default="",
        help="Comma list of methods to exclude from plots (e.g., rejection).",
    )
    ap.add_argument(
        "--paper_with_sample_phase",
        type=int,
        default=1,
        help="If 1, also output exp2_paper_sample_phase_vs_t.(pdf|png) when available (default: 1).",
    )
    return ap.parse_args()


# ----------------------------
# Data model
# ----------------------------
@dataclass(frozen=True)
class Agg:
    variant: str
    method: str
    t: int
    n: int
    median: float
    p95: float
    stdev: float


# ----------------------------
# Read sweep_raw.csv
# ----------------------------
def read_raw_rows(raw_path: str, warmup_reps: int) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(raw_path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            ok = to_int(row.get("ok", "0"), 0)
            if ok != 1:
                continue
            rep = to_int(row.get("rep", "0"), 0)
            if rep < warmup_reps:
                continue

            t = to_int(row.get("t", "0"), 0)
            if t <= 0:
                continue
            wall = to_float(row.get("wall_ms", "nan"))
            if math.isnan(wall) or wall <= 0:
                continue

            rows.append(row)
    return rows


def choose_error_mode(aggs: List[Agg], requested: str, min_repeats: int) -> str:
    if requested != "auto":
        return requested
    # If the runner tells us the effective repeats (min_repeats), prefer that.
    if min_repeats and min_repeats >= 10:
        return "p95"
    # Otherwise fall back to observed min-n.
    min_n = min((a.n for a in aggs), default=0)
    return "p95" if min_n >= 10 else "stdev"


def aggregate_metric(rows: List[Dict[str, str]], key: str) -> List[Agg]:
    grouped: Dict[Tuple[str, str, int], List[float]] = {}
    for row in rows:
        variant = row.get("variant", "")
        method = row.get("method", "")
        t = to_int(row.get("t", "0"), 0)
        v = to_float(row.get(key, "nan"))
        if t <= 0 or math.isnan(v) or v <= 0:
            continue
        grouped.setdefault((variant, method, t), []).append(v)

    out: List[Agg] = []
    for (variant, method, t), vals in grouped.items():
        out.append(
            Agg(
                variant=variant,
                method=method,
                t=t,
                n=len(vals),
                median=statistics.median(vals),
                p95=percentile_linear(vals, 0.95),
                stdev=safe_stdev(vals),
            )
        )
    return out


# ----------------------------
# Robust sample-phase extraction
# ----------------------------
_PHASE_KEYS_PRIORITY: Tuple[str, ...] = (
    "run_sample_ms",
    "sample_ms",
    "phase3_sample_ms",
    "phase3_ms",
    # Adaptive/fallback variants
    "run_fallback_sample_ms",
    "fallback_sample_ms",
    "fallback_sampling_ms",
    "adaptive_fallback_sampling_ms",
    "adaptive_sampling_ms",
    "phase3_fallback_sample_ms",
    "phase3_fallback_ms",
    "sample_phase_ms",
)


def _extract_sample_ms(phases: object) -> Optional[float]:
    if not isinstance(phases, dict):
        return None

    for k in _PHASE_KEYS_PRIORITY:
        if k in phases:
            try:
                v = float(phases[k])
            except Exception:
                continue
            if v > 0:
                return v

    cand: List[float] = []
    for k, v in phases.items():
        if not isinstance(k, str):
            continue
        lk = k.lower()
        if not lk.endswith("_ms"):
            continue
        if "sample" not in lk:
            continue
        try:
            fv = float(v)
        except Exception:
            continue
        if fv > 0:
            cand.append(fv)

    return max(cand) if cand else None


def aggregate_sample_phase(rows: List[Dict[str, str]]) -> List[Agg]:
    grouped: Dict[Tuple[str, str, int], List[float]] = {}
    for row in rows:
        variant = row.get("variant", "")
        method = row.get("method", "")
        t = to_int(row.get("t", "0"), 0)
        if t <= 0:
            continue

        phases_str = row.get("phases_json", "")
        if not phases_str:
            continue
        try:
            phases = json.loads(phases_str)
        except Exception:
            continue

        sample_ms = _extract_sample_ms(phases)
        if sample_ms is None:
            continue
        grouped.setdefault((variant, method, t), []).append(sample_ms)

    out: List[Agg] = []
    for (variant, method, t), vals in grouped.items():
        out.append(
            Agg(
                variant=variant,
                method=method,
                t=t,
                n=len(vals),
                median=statistics.median(vals),
                p95=percentile_linear(vals, 0.95),
                stdev=safe_stdev(vals),
            )
        )
    return out


# ----------------------------
# Δruntime (per-rep baseline)
# ----------------------------
def compute_delta_rows(rows: List[Dict[str, str]], t0: int, require_t0: bool) -> List[Agg]:
    by_vm: Dict[Tuple[str, str], List[Tuple[int, int, float]]] = {}
    for row in rows:
        variant = row.get("variant", "")
        method = row.get("method", "")
        t = to_int(row.get("t", "0"), 0)
        rep = to_int(row.get("rep", "0"), 0)
        wall = to_float(row.get("wall_ms", "nan"))
        if t <= 0 or math.isnan(wall) or wall <= 0:
            continue
        by_vm.setdefault((variant, method), []).append((t, rep, wall))

    grouped: Dict[Tuple[str, str, int], List[float]] = {}

    for (variant, method), pts in by_vm.items():
        ts = sorted({t for (t, _rep, _wall) in pts})
        if not ts:
            continue
        if require_t0 and (t0 not in ts):
            continue
        t0_used = t0 if (t0 in ts) else ts[0]

        baseline: Dict[int, float] = {}
        for (t, rep, wall) in pts:
            if t == t0_used:
                baseline[rep] = wall

        for (t, rep, wall) in pts:
            b = baseline.get(rep)
            if b is None:
                continue
            grouped.setdefault((variant, method, t), []).append(wall - b)

    out: List[Agg] = []
    for (variant, method, t), deltas in grouped.items():
        out.append(
            Agg(
                variant=variant,
                method=method,
                t=t,
                n=len(deltas),
                median=statistics.median(deltas),
                p95=percentile_linear(deltas, 0.95),
                stdev=safe_stdev(deltas),
            )
        )
    return out


# ----------------------------
# Filtering + reporting
# ----------------------------
def filter_by_min_repeats(aggs: List[Agg], min_repeats: int) -> Tuple[List[Agg], List[Tuple[str, str, int, int]]]:
    if min_repeats <= 0:
        return aggs, []
    kept: List[Agg] = []
    dropped: List[Tuple[str, str, int, int]] = []
    for a in aggs:
        if a.n >= min_repeats:
            kept.append(a)
        else:
            dropped.append((a.variant, a.method, a.t, a.n))
    return kept, dropped


def write_drop_report(out_dir: str, dropped: List[Tuple[str, str, int, int]], min_repeats: int) -> None:
    if not dropped:
        return
    path = os.path.join(out_dir, "exp2_dropped_points.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variant", "method", "t", "n_ok_after_warmup", "min_repeats_required"])
        for (variant, method, t, n) in sorted(dropped):
            w.writerow([variant, method, t, n, min_repeats])
    print(f"[EXP2][PLOT] wrote {path} (dropped {len(dropped)} points)")


# ----------------------------
# ns/sample table
# ----------------------------
def regression_slope(xs: List[float], ys: List[float]) -> float:
    if len(xs) < 2:
        return float("nan")
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    var = sum((x - mx) ** 2 for x in xs)
    if var <= 0:
        return float("nan")
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return cov / var


def write_ns_per_sample(out_dir: str, phase_aggs: List[Agg]) -> str:
    path = os.path.join(out_dir, "exp2_ns_per_sample.csv")
    grouped: Dict[Tuple[str, str], List[Tuple[int, float]]] = {}
    for a in phase_aggs:
        grouped.setdefault((a.variant, a.method), []).append((a.t, a.median))

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "variant",
                "method",
                "n_points",
                "t_min",
                "t_max",
                "sample_ms_min",
                "sample_ms_max",
                "ns_per_sample_endpoints",
                "ns_per_sample_regression",
                "ns_per_sample_tail3",
            ]
        )
        for (variant, method), pts in sorted(grouped.items()):
            pts = sorted(pts, key=lambda x: x[0])
            if len(pts) < 2:
                continue
            t_min, y_min = pts[0]
            t_max, y_max = pts[-1]

            slope_end_ms = (y_max - y_min) / (t_max - t_min) if t_max > t_min else float("nan")
            ns_end = slope_end_ms * 1e6

            xs = [float(t) for t, _y in pts]
            ys = [float(y) for _t, y in pts]
            slope_all = regression_slope(xs, ys)
            ns_reg = slope_all * 1e6

            if len(pts) >= 3:
                tail = pts[-3:]
                xs_t = [float(t) for t, _y in tail]
                ys_t = [float(y) for _t, y in tail]
                slope_tail = regression_slope(xs_t, ys_t)
                ns_tail = slope_tail * 1e6
            else:
                ns_tail = float("nan")

            w.writerow([variant, method, len(pts), t_min, t_max, y_min, y_max, ns_end, ns_reg, ns_tail])

    return path


# ----------------------------
# Plotting
# ----------------------------
def _yerr(a: Agg, mode: str) -> Tuple[float, float]:
    if mode == "stdev":
        return (0.0, max(0.0, a.stdev))
    return (0.0, max(0.0, a.p95 - a.median))


def _collect_series(aggs: List[Agg]) -> Dict[Tuple[str, str], List[Agg]]:
    out: Dict[Tuple[str, str], List[Agg]] = {}
    for a in aggs:
        out.setdefault((a.variant, a.method), []).append(a)
    for k in out:
        out[k] = sorted(out[k], key=lambda x: x.t)
    return out


def _variant_sort_key(v: str) -> int:
    order = {"sampling": 0, "enum_sampling": 1, "adaptive": 2}
    return order.get(v, 99)


def select_paper_methods(
    wall_aggs: List[Agg],
    variants: List[str],
    topk: int,
    always_include: List[str],
    paper_methods: List[str],
    exclude: List[str],
) -> List[str]:
    all_methods = sorted({a.method for a in wall_aggs})
    exclude_set = set(exclude)

    if paper_methods:
        methods = [m for m in paper_methods if (m in all_methods and m not in exclude_set)]
        for m in always_include:
            if m in all_methods and m not in exclude_set and m not in methods:
                methods.append(m)
        return methods

    rank_variant = "sampling" if ("sampling" in variants) else (variants[0] if variants else "")
    cand = [a for a in wall_aggs if a.variant == rank_variant]
    if not cand:
        cand = wall_aggs[:]

    t_max = max((a.t for a in cand), default=None)
    rank_pts = [a for a in cand if (t_max is not None and a.t == t_max)]
    if not rank_pts:
        rank_pts = cand

    rank_pts = [a for a in rank_pts if a.method not in exclude_set]
    rank_pts.sort(key=lambda x: x.median)

    methods: List[str] = []
    for a in rank_pts:
        if a.method not in methods:
            methods.append(a.method)
        if len(methods) >= max(1, topk):
            break

    for m in always_include:
        if m in all_methods and m not in exclude_set and m not in methods:
            methods.append(m)

    return methods


def plot_paper_panels(
    out_dir: str,
    aggs: List[Agg],
    methods: List[str],
    error_mode: str,
    filename: str,
    ylabel: str,
    yscale: str,
    title: str,
    *,
    add_hline0: bool = False,
    symlog_linthresh: float = 1.0,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[EXP2][PLOT] matplotlib not available; skip plots. ({e})")
        return

    variants = sorted({a.variant for a in aggs}, key=_variant_sort_key)
    series = _collect_series(aggs)

    fig, axes = plt.subplots(1, len(variants), figsize=(4.4 * max(1, len(variants)), 3.2), squeeze=False)
    axes_list = axes[0]

    for ax, variant in zip(axes_list, variants):
        for method in methods:
            pts = series.get((variant, method))
            if not pts:
                continue
            ts = [p.t for p in pts]
            ys = [p.median for p in pts]
            low = []
            upp = []
            for p in pts:
                l, u = _yerr(p, error_mode)
                low.append(l)
                upp.append(u)

            ax.errorbar(ts, ys, yerr=[low, upp], marker="o", linewidth=1, markersize=3, label=method)

        ax.set_xscale("log")
        if yscale == "symlog":
            ax.set_yscale("symlog", linthresh=symlog_linthresh)
        else:
            ax.set_yscale(yscale)

        if add_hline0:
            ax.axhline(0.0, linewidth=0.8)

        ax.set_xlabel("t (samples)")
        ax.set_title(variant)

    axes_list[0].set_ylabel(ylabel)
    handles, labels = axes_list[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=min(len(labels), 6), fontsize=8, frameon=False)
        fig.subplots_adjust(bottom=0.28)

    fig.suptitle(title, y=0.98)
    fig.tight_layout()

    pdf = os.path.join(out_dir, f"{filename}.pdf")
    png = os.path.join(out_dir, f"{filename}.png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=220)
    plt.close(fig)


def _plot_per_variant(out_dir: str, aggs: List[Agg], filename_prefix: str, ylabel: str, yscale: str, error_mode: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return

    variants = sorted({a.variant for a in aggs}, key=_variant_sort_key)
    series = _collect_series(aggs)

    for variant in variants:
        methods = sorted({a.method for a in aggs if a.variant == variant})
        plt.figure(figsize=(5.2, 3.6))
        for method in methods:
            pts = series.get((variant, method))
            if not pts:
                continue
            ts = [p.t for p in pts]
            ys = [p.median for p in pts]
            low = []
            upp = []
            for p in pts:
                l, u = _yerr(p, error_mode)
                low.append(l)
                upp.append(u)
            plt.errorbar(ts, ys, yerr=[low, upp], marker="o", linewidth=1, markersize=3, label=method)

        plt.xscale("log")
        if yscale == "symlog":
            plt.yscale("symlog", linthresh=1.0)
            plt.axhline(0.0, linewidth=0.8)
        else:
            plt.yscale(yscale)

        plt.xlabel("t (samples)")
        plt.ylabel(ylabel)
        plt.title(f"EXP-2 {filename_prefix} ({variant})")
        plt.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"exp2_{filename_prefix}_{variant}.pdf"))
        plt.savefig(os.path.join(out_dir, f"exp2_{filename_prefix}_{variant}.png"), dpi=200)
        plt.close()


def write_readme(
    out_dir: str,
    args: argparse.Namespace,
    variants: List[str],
    methods_paper: List[str],
    error_mode: str,
    n_raw_rows: int,
    n_wall_points: int,
    dropped_points: int,
) -> None:
    path = os.path.join(out_dir, "EXP2_README.txt")
    with open(path, "w") as f:
        f.write("EXP-2 (Runtime vs t) - plotting artifact\n\n")
        f.write(f"mode={args.mode}\n")
        f.write(f"warmup_reps={args.warmup_reps}\n")
        f.write(f"min_repeats={args.min_repeats}\n")
        f.write(f"t0={args.t0}\n")
        f.write(f"error_mode={error_mode} (requested={args.error})\n")
        f.write(f"paper_with_sample_phase={args.paper_with_sample_phase}\n")
        f.write(f"variants={','.join(variants)}\n\n")

        f.write(f"raw_ok_rows_used(after warmup)={n_raw_rows}\n")
        f.write(f"wall_points_aggregated={n_wall_points}\n")
        f.write(f"points_dropped_by_min_repeats={dropped_points}\n\n")

        f.write("Paper plots:\n")
        f.write("  exp2_paper_runtime_vs_t.(pdf|png)\n")
        f.write("  exp2_paper_delta_vs_t.(pdf|png)\n")
        if args.paper_with_sample_phase:
            f.write("  exp2_paper_sample_phase_vs_t.(pdf|png)  [if phases_json available]\n")
        f.write("\nTables:\n")
        f.write("  exp2_ns_per_sample.csv\n\n")

        f.write("Paper method set:\n")
        f.write("  " + ", ".join(methods_paper) + "\n\n")

        if args.mode == "full":
            f.write("Full mode also writes per-variant diagnostics:\n")
            f.write("  exp2_runtime_vs_t_<variant>.(pdf|png)\n")
            f.write("  exp2_delta_runtime_vs_t_<variant>.(pdf|png)\n")
            f.write("  exp2_sample_phase_vs_t_<variant>.(pdf|png)\n\n")

        if dropped_points > 0:
            f.write("Note:\n")
            f.write("  Some points were dropped from paper plots due to insufficient successful repeats.\n")
            f.write("  See exp2_dropped_points.csv.\n")

    print(f"[EXP2][PLOT] wrote {path}")


def main() -> None:
    args = parse_args()
    out_dir = os.path.abspath(args.out_dir)
    raw_path = os.path.join(out_dir, "sweep_raw.csv")
    ensure_file(raw_path)

    raw_rows = read_raw_rows(raw_path, warmup_reps=max(0, args.warmup_reps))
    if not raw_rows:
        print("[EXP2][PLOT] No successful rows in sweep_raw.csv after warmup filtering; nothing to plot.")
        return

    wall_aggs_all = aggregate_metric(raw_rows, "wall_ms")
    phase_aggs_all = aggregate_sample_phase(raw_rows)

    # Δruntime uses wall_ms only
    delta_aggs_paper_all = compute_delta_rows(raw_rows, t0=args.t0, require_t0=True)
    delta_aggs_full_all = compute_delta_rows(raw_rows, t0=args.t0, require_t0=False)

    # Apply min_repeats filtering for PAPER plots (and for slope table to keep it reliable).
    wall_aggs, dropped_wall = filter_by_min_repeats(wall_aggs_all, args.min_repeats)
    delta_aggs_paper, dropped_delta = filter_by_min_repeats(delta_aggs_paper_all, args.min_repeats)
    phase_aggs, dropped_phase = filter_by_min_repeats(phase_aggs_all, args.min_repeats)

    dropped = dropped_wall + dropped_delta + dropped_phase
    write_drop_report(out_dir, dropped, args.min_repeats)

    error_mode = choose_error_mode(wall_aggs if wall_aggs else wall_aggs_all, args.error, args.min_repeats)

    variants = sorted({a.variant for a in wall_aggs_all}, key=_variant_sort_key)

    always_include = split_csv_list(args.always_include)
    paper_methods_arg = split_csv_list(args.paper_methods)
    exclude = split_csv_list(args.exclude_methods)

    # Method selection should be based on the filtered (reliable) wall aggs if available.
    wall_for_rank = wall_aggs if wall_aggs else wall_aggs_all
    methods_paper = select_paper_methods(
        wall_aggs=wall_for_rank,
        variants=variants,
        topk=max(1, args.topk),
        always_include=always_include,
        paper_methods=paper_methods_arg,
        exclude=exclude,
    )

    print(f"[EXP2][PLOT] mode={args.mode} warmup_reps={args.warmup_reps} t0={args.t0} min_repeats={args.min_repeats}")
    print(f"[EXP2][PLOT] error_mode={error_mode} (requested={args.error})")
    print(f"[EXP2][PLOT] paper_methods={methods_paper}")

    # Always write ns/sample table (use filtered phase points if min_repeats is set).
    ns_path = write_ns_per_sample(out_dir, phase_aggs if phase_aggs else phase_aggs_all)
    print(f"[EXP2][PLOT] wrote {ns_path}")

    # Paper plots (concise)
    plot_paper_panels(
        out_dir=out_dir,
        aggs=wall_aggs if wall_aggs else wall_aggs_all,
        methods=methods_paper,
        error_mode=error_mode,
        filename="exp2_paper_runtime_vs_t",
        ylabel="wall runtime (median ms)",
        yscale="log",
        title="EXP-2 Runtime vs t",
    )

    plot_paper_panels(
        out_dir=out_dir,
        aggs=delta_aggs_paper if delta_aggs_paper else delta_aggs_paper_all,
        methods=methods_paper,
        error_mode=error_mode,
        filename="exp2_paper_delta_vs_t",
        ylabel=f"Δ wall runtime (median ms)  [baseline t0={args.t0}]",
        yscale="symlog",
        title="EXP-2 ΔRuntime vs t",
        add_hline0=True,
        symlog_linthresh=1.0,
    )

    if args.paper_with_sample_phase and (phase_aggs or phase_aggs_all):
        plot_paper_panels(
            out_dir=out_dir,
            aggs=phase_aggs if phase_aggs else phase_aggs_all,
            methods=methods_paper,
            error_mode=error_mode,
            filename="exp2_paper_sample_phase_vs_t",
            ylabel="sample-phase time (median ms)",
            yscale="log",
            title="EXP-2 Sample-phase vs t",
        )

    # Full diagnostics if requested
    if args.mode == "full":
        _plot_per_variant(out_dir, wall_aggs_all, "runtime_vs_t", "wall runtime (median ms)", "log", error_mode)
        _plot_per_variant(out_dir, delta_aggs_full_all, "delta_runtime_vs_t", f"Δ wall runtime (median ms) [t0={args.t0}]", "symlog", error_mode)
        if phase_aggs_all:
            _plot_per_variant(out_dir, phase_aggs_all, "sample_phase_vs_t", "sample-phase time (median ms)", "log", error_mode)

    write_readme(
        out_dir=out_dir,
        args=args,
        variants=variants,
        methods_paper=methods_paper,
        error_mode=error_mode,
        n_raw_rows=len(raw_rows),
        n_wall_points=len(wall_aggs_all),
        dropped_points=len(dropped),
    )

    print(f"[EXP2][PLOT] DONE. Out dir: {out_dir}")


if __name__ == "__main__":
    main()
