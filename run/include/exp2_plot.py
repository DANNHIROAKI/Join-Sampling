#!/usr/bin/env python3
"""exp2_plot.py

EXP-2 plotting (Runtime vs t) with more "paper-friendly" defaults.

Compared to the earlier version, this script addresses the two most common
EXP-2 presentation pitfalls:

  1) repeats 太少却画 p95：
     - supports --error=auto (default)
     - auto picks p95 if repeats>=10, otherwise stdev

  2) Δruntime 噪声大、甚至出现负值：
     - Δruntime is computed from sweep_raw.csv **per rep** by subtracting the
       matched baseline run at t0 (same method, variant, rep).
     - then we aggregate deltas with median/p95/stdev.

Inputs (under --out_dir):
  sweep_summary.csv
  sweep_raw.csv

Outputs (written into --out_dir):
  exp2_runtime_vs_t_<variant>.pdf/.png
  exp2_delta_runtime_vs_t_<variant>.pdf/.png
  exp2_sample_phase_vs_t_<variant>.pdf/.png
  exp2_ns_per_sample.csv
  EXP2_README.txt
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple


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


def percentile_linear(vals: List[float], p: float) -> float:
    """Linear interpolation percentile. p in [0,1]."""
    if not vals:
        return float("nan")
    vals = sorted(vals)
    if len(vals) == 1:
        return vals[0]
    pos = p * (len(vals) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    w = pos - lo
    return vals[lo] * (1 - w) + vals[hi] * w


def safe_stdev(vals: List[float]) -> float:
    if len(vals) < 2:
        return 0.0
    return statistics.stdev(vals)


def ensure_file(path: str) -> None:
    if not os.path.isfile(path):
        raise FileNotFoundError(path)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", required=True, help="EXP-2 output directory containing sweep_raw.csv and sweep_summary.csv")
    ap.add_argument("--t0", type=int, default=1000, help="Baseline t0 for Δruntime plot (default: 1000)")
    ap.add_argument(
        "--error",
        choices=["auto", "p95", "stdev"],
        default="auto",
        help="Error bar type (default: auto => p95 if repeats>=10 else stdev)",
    )
    return ap.parse_args()


def read_summary(summary_path: str) -> List[Dict[str, str]]:
    """Read sweep_summary.csv and keep only fully successful points."""
    rows: List[Dict[str, str]] = []
    with open(summary_path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            ok_rate = to_float(row.get("ok_rate", "0"), 0.0)
            if ok_rate < 1.0:
                continue
            wall_med = to_float(row.get("wall_median_ms", "-1"), -1.0)
            if wall_med <= 0:
                continue
            rows.append(row)
    return rows


@dataclass
class DeltaAgg:
    variant: str
    method: str
    t: int
    n: int
    median_ms: float
    p95_ms: float
    stdev_ms: float
    t0_used: int


def read_raw_wall(raw_path: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(raw_path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            ok = to_int(row.get("ok", "0"), 0)
            if ok != 1:
                continue
            wall = to_float(row.get("wall_ms", "-1"), -1.0)
            if wall <= 0:
                continue
            rows.append(row)
    return rows


def choose_error_mode(summary_rows: List[Dict[str, str]], requested: str) -> str:
    if requested != "auto":
        return requested
    # If repeats is missing, default to stdev.
    reps: List[int] = []
    for r in summary_rows:
        reps.append(to_int(r.get("repeats", "0"), 0))
    min_rep = min(reps) if reps else 0
    # Use p95 only if it is meaningful.
    if min_rep >= 10:
        # also require p95 column to be present and valid
        any_valid = False
        for r in summary_rows:
            p95 = to_float(r.get("wall_p95_ms", "nan"))
            if not math.isnan(p95) and p95 > 0:
                any_valid = True
                break
        if any_valid:
            return "p95"
    return "stdev"


def get_yerr_summary(row: Dict[str, str], mode: str) -> Tuple[float, float]:
    """Upper-only errorbar: return (lower, upper)."""
    if mode == "stdev":
        err = max(0.0, to_float(row.get("wall_stdev_ms", "0"), 0.0))
        return (0.0, err)
    # p95
    med = to_float(row.get("wall_median_ms", "nan"))
    p95 = to_float(row.get("wall_p95_ms", "nan"))
    if math.isnan(med) or math.isnan(p95):
        return (0.0, 0.0)
    return (0.0, max(0.0, p95 - med))


def compute_delta_from_raw(raw_rows: List[Dict[str, str]], t0: int) -> List[DeltaAgg]:
    """Compute Δwall per (variant,method,t) by matching baseline at t0 per rep."""
    # Key baseline: (variant, method, rep) -> wall_ms at t0_used
    # We may need to fall back to min t if requested t0 not present.

    # Collect all points per (variant,method)
    by_vm: Dict[Tuple[str, str], List[Tuple[int, int, float]]] = {}
    for r in raw_rows:
        variant = r.get("variant", "")
        method = r.get("method", "")
        t = to_int(r.get("t", "0"), 0)
        rep = to_int(r.get("rep", "0"), 0)
        wall = to_float(r.get("wall_ms", "nan"))
        if t <= 0 or math.isnan(wall):
            continue
        by_vm.setdefault((variant, method), []).append((t, rep, wall))

    out: List[DeltaAgg] = []

    for (variant, method), pts in sorted(by_vm.items()):
        # Determine which t0 to use for this (variant,method)
        ts = sorted({t for (t, _rep, _wall) in pts})
        if not ts:
            continue
        t0_used = t0 if t0 in ts else ts[0]

        baseline: Dict[int, float] = {}
        for (t, rep, wall) in pts:
            if t == t0_used:
                baseline[rep] = wall

        # For each t, collect deltas for reps that have a baseline
        deltas_by_t: Dict[int, List[float]] = {}
        for (t, rep, wall) in pts:
            b = baseline.get(rep)
            if b is None:
                continue
            deltas_by_t.setdefault(t, []).append(wall - b)

        for t, deltas in sorted(deltas_by_t.items()):
            if not deltas:
                continue
            med = statistics.median(deltas)
            p95 = percentile_linear(deltas, 0.95)
            sd = safe_stdev(deltas)
            out.append(
                DeltaAgg(
                    variant=variant,
                    method=method,
                    t=t,
                    n=len(deltas),
                    median_ms=med,
                    p95_ms=p95,
                    stdev_ms=sd,
                    t0_used=t0_used,
                )
            )
    return out


@dataclass
class PhaseAgg:
    variant: str
    method: str
    t: int
    n: int
    median_ms: float
    p95_ms: float
    stdev_ms: float


def read_sample_phase(raw_path: str) -> List[PhaseAgg]:
    """Extract sample-phase time from phases_json in sweep_raw.csv."""
    samples: Dict[Tuple[str, str, int], List[float]] = {}
    phase_keys = ("run_sample_ms", "sample_ms", "phase3_ms", "phase3_sample_ms")

    with open(raw_path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            ok = to_int(row.get("ok", "0"), 0)
            if ok != 1:
                continue

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

            sample_ms: Optional[float] = None
            for k in phase_keys:
                if k in phases:
                    try:
                        sample_ms = float(phases[k])
                    except Exception:
                        sample_ms = None
                    break
            if sample_ms is None or sample_ms <= 0:
                continue

            samples.setdefault((variant, method, t), []).append(sample_ms)

    out: List[PhaseAgg] = []
    for (variant, method, t), vals in samples.items():
        if not vals:
            continue
        out.append(
            PhaseAgg(
                variant=variant,
                method=method,
                t=t,
                n=len(vals),
                median_ms=statistics.median(vals),
                p95_ms=percentile_linear(vals, 0.95),
                stdev_ms=safe_stdev(vals),
            )
        )
    return out


def get_yerr_delta(row: DeltaAgg, mode: str) -> Tuple[float, float]:
    if mode == "stdev":
        return (0.0, max(0.0, row.stdev_ms))
    # p95
    return (0.0, max(0.0, row.p95_ms - row.median_ms))


def get_yerr_phase(row: PhaseAgg, mode: str) -> Tuple[float, float]:
    if mode == "stdev":
        return (0.0, max(0.0, row.stdev_ms))
    return (0.0, max(0.0, row.p95_ms - row.median_ms))


def regression_slope(xs: List[float], ys: List[float]) -> float:
    """Return slope in (y units)/(x units)."""
    if len(xs) < 2:
        return float("nan")
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    var = sum((x - mx) ** 2 for x in xs)
    if var <= 0:
        return float("nan")
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return cov / var


def main() -> None:
    args = parse_args()
    out_dir = os.path.abspath(args.out_dir)
    summary_path = os.path.join(out_dir, "sweep_summary.csv")
    raw_path = os.path.join(out_dir, "sweep_raw.csv")

    ensure_file(summary_path)
    ensure_file(raw_path)

    # Matplotlib optional
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[EXP2][PLOT] matplotlib not available; will only write CSV/README. ({e})")
        plt = None  # type: ignore

    summary_rows = read_summary(summary_path)
    if not summary_rows:
        print("[EXP2][PLOT] No successful rows in sweep_summary.csv; nothing to plot.")
        return

    error_mode = choose_error_mode(summary_rows, args.error)
    min_repeats = min(to_int(r.get("repeats", "0"), 0) for r in summary_rows)
    print(f"[EXP2][PLOT] error_mode={error_mode} (requested={args.error}, min_repeats={min_repeats})")

    raw_rows = read_raw_wall(raw_path)
    delta_agg = compute_delta_from_raw(raw_rows, args.t0)
    phase_agg = read_sample_phase(raw_path)

    variants = sorted({row["variant"] for row in summary_rows})

    # ----------------------------
    # Plot 1: wall runtime vs t (from summary)
    # ----------------------------
    if plt is not None:
        for variant in variants:
            sub = [row for row in summary_rows if row["variant"] == variant]
            methods = sorted({row["method"] for row in sub})

            plt.figure()
            for method in methods:
                g = [row for row in sub if row["method"] == method]
                g.sort(key=lambda x: to_int(x["t"]))
                ts = [to_int(x["t"]) for x in g]
                ys = [to_float(x["wall_median_ms"]) for x in g]
                low: List[float] = []
                upp: List[float] = []
                for x in g:
                    l, u = get_yerr_summary(x, error_mode)
                    low.append(l)
                    upp.append(u)

                plt.errorbar(ts, ys, yerr=[low, upp], marker="o", label=method, linewidth=1)

            plt.xscale("log")
            plt.yscale("log")
            plt.xlabel("t (samples)")
            plt.ylabel("runtime (median ms)")
            plt.title(f"EXP-2 Runtime vs t ({variant})")
            plt.legend(fontsize=7)
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"exp2_runtime_vs_t_{variant}.pdf"))
            plt.savefig(os.path.join(out_dir, f"exp2_runtime_vs_t_{variant}.png"), dpi=200)
            plt.close()

    # ----------------------------
    # Plot 2: Δ wall runtime vs t (from raw, per-rep delta)
    # ----------------------------
    if plt is not None and delta_agg:
        for variant in variants:
            sub = [row for row in delta_agg if row.variant == variant]
            if not sub:
                continue
            methods = sorted({row.method for row in sub})

            # Determine t0_used per (variant,method) (should be consistent across t)
            t0_map: Dict[str, int] = {}
            for r in sub:
                t0_map.setdefault(r.method, r.t0_used)

            plt.figure()
            for method in methods:
                g = [row for row in sub if row.method == method]
                g.sort(key=lambda x: x.t)
                ts = [x.t for x in g]
                ys = [x.median_ms for x in g]
                low: List[float] = []
                upp: List[float] = []
                for x in g:
                    l, u = get_yerr_delta(x, error_mode)
                    low.append(l)
                    upp.append(u)

                plt.errorbar(ts, ys, yerr=[low, upp], marker="o", label=method, linewidth=1)

            plt.xscale("log")
            plt.yscale("linear")
            # Prefer the requested baseline t0; if missing, we used per-method min-t.
            # For the y-label, show the mode baseline if consistent.
            t0_useds = sorted(set(t0_map.values()))
            if len(t0_useds) == 1:
                t0_label = t0_useds[0]
            else:
                t0_label = args.t0

            plt.xlabel("t (samples)")
            plt.ylabel(f"Δruntime (median ms)  [per-rep baseline @ t={t0_label}]")
            plt.title(f"EXP-2 ΔRuntime vs t ({variant})")
            plt.axhline(0.0, linewidth=0.8)
            plt.legend(fontsize=7)
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"exp2_delta_runtime_vs_t_{variant}.pdf"))
            plt.savefig(os.path.join(out_dir, f"exp2_delta_runtime_vs_t_{variant}.png"), dpi=200)
            plt.close()

    if plt is not None and not delta_agg:
        print("[EXP2][PLOT] Δruntime: no usable raw wall_ms data; skip delta plots.")

    # ----------------------------
    # Plot 3: sample-phase vs t (from raw phases_json)
    # ----------------------------
    if plt is not None and phase_agg:
        phase_variants = sorted({row.variant for row in phase_agg})
        for variant in phase_variants:
            sub = [row for row in phase_agg if row.variant == variant]
            methods = sorted({row.method for row in sub})

            plt.figure()
            for method in methods:
                g = [row for row in sub if row.method == method]
                g.sort(key=lambda x: x.t)
                ts = [x.t for x in g]
                ys = [x.median_ms for x in g]
                low: List[float] = []
                upp: List[float] = []
                for x in g:
                    l, u = get_yerr_phase(x, error_mode)
                    low.append(l)
                    upp.append(u)

                plt.errorbar(ts, ys, yerr=[low, upp], marker="o", label=method, linewidth=1)

            plt.xscale("log")
            plt.yscale("log")
            plt.xlabel("t (samples)")
            plt.ylabel("sample-phase time (median ms)")
            plt.title(f"EXP-2 Sample-phase vs t ({variant})")
            plt.legend(fontsize=7)
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"exp2_sample_phase_vs_t_{variant}.pdf"))
            plt.savefig(os.path.join(out_dir, f"exp2_sample_phase_vs_t_{variant}.png"), dpi=200)
            plt.close()

    if plt is not None and not phase_agg:
        print("[EXP2][PLOT] phases_json sample-phase not found; skip sample-phase plots.")

    # ----------------------------
    # Report: ns/sample slopes from sample-phase medians
    # ----------------------------
    slopes_path = os.path.join(out_dir, "exp2_ns_per_sample.csv")
    grouped: Dict[Tuple[str, str], List[Tuple[int, float]]] = {}
    for row in phase_agg:
        grouped.setdefault((row.variant, row.method), []).append((row.t, row.median_ms))

    with open(slopes_path, "w", newline="") as f:
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
            # endpoints
            slope_ms = (y_max - y_min) / (t_max - t_min) if t_max > t_min else float("nan")
            ns_end = slope_ms * 1e6

            # regression over all points
            xs = [float(t) for t, _y in pts]
            ys = [float(y) for _t, y in pts]
            slope_all = regression_slope(xs, ys)
            ns_reg = slope_all * 1e6

            # tail-3 regression (more representative for large t)
            if len(pts) >= 3:
                tail = pts[-3:]
                xs_t = [float(t) for t, _y in tail]
                ys_t = [float(y) for _t, y in tail]
                slope_tail = regression_slope(xs_t, ys_t)
                ns_tail = slope_tail * 1e6
            else:
                ns_tail = float("nan")

            w.writerow([variant, method, len(pts), t_min, t_max, y_min, y_max, ns_end, ns_reg, ns_tail])

    print(f"[EXP2][PLOT] wrote {slopes_path}")

    # ----------------------------
    # README artifact
    # ----------------------------
    readme_path = os.path.join(out_dir, "EXP2_README.txt")
    with open(readme_path, "w") as f:
        f.write("EXP-2 (Runtime vs t) - plots\n\n")
        f.write("Inputs:\n")
        f.write("  sweep_raw.csv      : one row per run/repeat (contains phases_json)\n")
        f.write("  sweep_summary.csv  : aggregated stats per (method,variant,t)\n\n")
        f.write("Outputs:\n")
        f.write("  exp2_runtime_vs_t_<variant>.pdf/.png          : wall median vs t (log-x, log-y, upper-only error)\n")
        f.write("  exp2_delta_runtime_vs_t_<variant>.pdf/.png    : Δ wall median vs t (log-x, linear-y)\n")
        f.write("  exp2_sample_phase_vs_t_<variant>.pdf/.png     : sample-phase (median run_sample_ms) vs t (log-x, log-y)\n")
        f.write("  exp2_ns_per_sample.csv                        : ns/sample (endpoints + regression)\n\n")
        f.write("Error bars:\n")
        f.write(f"  requested={args.error}\n")
        f.write(f"  chosen={error_mode}\n")
        f.write(f"  min_repeats_over_plotted_points={min_repeats}\n\n")
        f.write("Δruntime definition:\n")
        f.write("  For each (variant,method,rep), subtract the matched baseline run at t0 (or min t if t0 missing).\n")
        f.write("  This reduces noise compared to subtracting aggregated medians.\n\n")
        f.write("Notes:\n")
        f.write("  - Points with ok_rate < 1.0 are filtered out.\n")
        f.write("  - Some method/variant combos may be intentionally skipped by sjs_sweep (e.g., rejection+sampling).\n")

        # Quick warning about negative delta medians
        neg = [d for d in delta_agg if d.median_ms < -1e-9]
        if neg:
            f.write("\nWarnings:\n")
            f.write(f"  - {len(neg)} Δruntime points have negative median (likely measurement noise).\n")

    print(f"[EXP2][PLOT] wrote {readme_path}")
    print(f"[EXP2][PLOT] DONE. Out dir: {out_dir}")


if __name__ == "__main__":
    main()
