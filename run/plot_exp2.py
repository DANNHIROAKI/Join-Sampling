#!/usr/bin/env python3
"""
run/plot_exp2.py

Enhanced EXP-2 plotting:
  - runtime vs t (wall median; log-x log-y)
  - Δruntime vs t (baseline-subtracted; log-x linear-y)
  - sample-phase vs t (median run_sample_ms from phases_json; log-x log-y)
  - ns/sample slope report (CSV)

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
from typing import Dict, List, Tuple, Any


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


def ensure_file(path: str) -> None:
    if not os.path.isfile(path):
        raise FileNotFoundError(path)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", required=True, help="EXP-2 output directory containing sweep_raw.csv and sweep_summary.csv")
    ap.add_argument("--t0", type=int, default=1000, help="Baseline t0 for Δruntime plot (default: 1000)")
    ap.add_argument("--error", choices=["p95", "stdev"], default="p95", help="Error bar type (default: p95)")
    return ap.parse_args()


def read_summary(summary_path: str) -> List[Dict[str, str]]:
    """
    Read sweep_summary.csv and keep only fully successful points.
    We also drop sentinel/invalid runtime values to keep log-y safe.
    """
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


def read_sample_phase(raw_path: str) -> List[Dict[str, Any]]:
    """
    Parse sweep_raw.csv and extract sample-phase time from phases_json.
    Aggregate per (variant, method, t) with median and p95.
    """
    samples: Dict[Tuple[str, str, int], List[float]] = {}

    # common keys we may see in phases_json across versions
    phase_keys = ("run_sample_ms", "sample_ms", "phase3_ms", "phase3_sample_ms")

    with open(raw_path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            ok = to_int(row.get("ok", "0"))
            if ok != 1:
                continue

            variant = row.get("variant", "")
            method = row.get("method", "")
            t = to_int(row.get("t", "0"))

            phases_str = row.get("phases_json", "")
            if not phases_str:
                continue

            try:
                phases = json.loads(phases_str)
            except Exception:
                continue

            sample_ms = None
            for k in phase_keys:
                if k in phases:
                    sample_ms = phases[k]
                    break
            if sample_ms is None:
                continue

            try:
                v = float(sample_ms)
            except Exception:
                continue

            # ignore non-positive to keep log-y safe
            if v <= 0:
                continue

            samples.setdefault((variant, method, t), []).append(v)

    agg: List[Dict[str, Any]] = []
    for (variant, method, t), vals in samples.items():
        if not vals:
            continue
        med = statistics.median(vals)
        p95 = percentile_linear(vals, 0.95)
        agg.append(
            {
                "variant": variant,
                "method": method,
                "t": t,
                "sample_median_ms": med,
                "sample_p95_ms": p95,
                "n": len(vals),
            }
        )
    return agg


def get_yerr(row: Dict[str, str], mode: str) -> Tuple[float, float]:
    """
    Return (lower, upper) errors for upper-only errorbar.
    For log-y safety, we default lower=0.
    """
    if mode == "stdev":
        err = max(0.0, to_float(row.get("wall_stdev_ms", "0"), 0.0))
        return (0.0, err)
    # p95
    med = to_float(row.get("wall_median_ms", "nan"))
    p95 = to_float(row.get("wall_p95_ms", "nan"))
    if math.isnan(med) or math.isnan(p95):
        return (0.0, 0.0)
    return (0.0, max(0.0, p95 - med))


def get_phase_yerr(row: Dict[str, Any], mode: str) -> Tuple[float, float]:
    if mode == "stdev":
        # raw-phase stdev isn't stored in our agg; use 0
        return (0.0, 0.0)
    # p95
    med = float(row["sample_median_ms"])
    p95 = float(row["sample_p95_ms"])
    return (0.0, max(0.0, p95 - med))


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
        print(f"[EXP2][PLOT] matplotlib not available; skip plotting. ({e})")
        return

    summary_rows = read_summary(summary_path)
    if not summary_rows:
        print("[EXP2][PLOT] No successful rows in sweep_summary.csv; nothing to plot.")
        return

    phase_agg = read_sample_phase(raw_path)

    variants = sorted({row["variant"] for row in summary_rows})

    # ----------------------------
    # Plot 1: wall runtime vs t
    # ----------------------------
    for variant in variants:
        sub = [row for row in summary_rows if row["variant"] == variant]
        methods = sorted({row["method"] for row in sub})

        plt.figure()
        for method in methods:
            g = [row for row in sub if row["method"] == method]
            g.sort(key=lambda x: to_int(x["t"]))
            ts = [to_int(x["t"]) for x in g]
            ys = [to_float(x["wall_median_ms"]) for x in g]
            low = []
            upp = []
            for x in g:
                l, u = get_yerr(x, args.error)
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
    # Plot 2: Δ wall runtime vs t
    # ----------------------------
    for variant in variants:
        sub = [row for row in summary_rows if row["variant"] == variant]
        methods = sorted({row["method"] for row in sub})

        plt.figure()
        for method in methods:
            g = [row for row in sub if row["method"] == method]
            g.sort(key=lambda x: to_int(x["t"]))
            ts = [to_int(x["t"]) for x in g]
            ys = [to_float(x["wall_median_ms"]) for x in g]
            low = []
            upp = []
            for x in g:
                l, u = get_yerr(x, args.error)
                low.append(l)
                upp.append(u)

            # baseline selection
            if args.t0 in ts:
                y0 = ys[ts.index(args.t0)]
                t0_used = args.t0
            else:
                y0 = ys[0]
                t0_used = ts[0]

            dys = [y - y0 for y in ys]

            plt.errorbar(ts, dys, yerr=[low, upp], marker="o", label=method, linewidth=1)

        plt.xscale("log")
        plt.yscale("linear")  # delta can be 0
        plt.xlabel("t (samples)")
        plt.ylabel(f"Δruntime (median ms)  [baseline @ t={t0_used}]")
        plt.title(f"EXP-2 ΔRuntime vs t ({variant})")
        plt.axhline(0.0, linewidth=0.8)
        plt.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"exp2_delta_runtime_vs_t_{variant}.pdf"))
        plt.savefig(os.path.join(out_dir, f"exp2_delta_runtime_vs_t_{variant}.png"), dpi=200)
        plt.close()

    # ----------------------------
    # Plot 3: sample-phase vs t (from raw phases_json)
    # ----------------------------
    if phase_agg:
        phase_variants = sorted({row["variant"] for row in phase_agg})
        for variant in phase_variants:
            sub = [row for row in phase_agg if row["variant"] == variant]
            methods = sorted({row["method"] for row in sub})

            plt.figure()
            for method in methods:
                g = [row for row in sub if row["method"] == method]
                g.sort(key=lambda x: int(x["t"]))
                ts = [int(x["t"]) for x in g]
                ys = [float(x["sample_median_ms"]) for x in g]
                low = []
                upp = []
                for x in g:
                    l, u = get_phase_yerr(x, args.error)
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
    else:
        print("[EXP2][PLOT] phases_json sample-phase not found; skip sample-phase plots.")

    # ----------------------------
    # Report: ns/sample slope from sample-phase medians (endpoints)
    # ----------------------------
    slopes_path = os.path.join(out_dir, "exp2_ns_per_sample.csv")
    if phase_agg:
        grouped: Dict[Tuple[str, str], List[Tuple[int, float]]] = {}
        for row in phase_agg:
            key = (row["variant"], row["method"])
            grouped.setdefault(key, []).append((int(row["t"]), float(row["sample_median_ms"])))

        with open(slopes_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(
                [
                    "variant",
                    "method",
                    "t_min",
                    "t_max",
                    "sample_ms_min",
                    "sample_ms_max",
                    "ns_per_sample_endpoints",
                ]
            )
            for (variant, method), pts in sorted(grouped.items()):
                pts = sorted(pts, key=lambda x: x[0])
                if len(pts) < 2:
                    continue
                t0, y0 = pts[0]
                t1, y1 = pts[-1]
                if t1 <= t0:
                    continue
                slope_ms = (y1 - y0) / (t1 - t0)
                ns_per_sample = slope_ms * 1e6
                w.writerow([variant, method, t0, t1, y0, y1, ns_per_sample])

        print(f"[EXP2][PLOT] wrote {slopes_path}")
    else:
        # still write an empty file to avoid surprises in pipelines
        with open(slopes_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["variant", "method", "t_min", "t_max", "sample_ms_min", "sample_ms_max", "ns_per_sample_endpoints"])
        print(f"[EXP2][PLOT] wrote empty {slopes_path} (no phase data)")

    # ----------------------------
    # README artifact
    # ----------------------------
    readme_path = os.path.join(out_dir, "EXP2_README.txt")
    with open(readme_path, "w") as f:
        f.write("EXP-2 (Runtime vs t) - enhanced plots\n\n")
        f.write("Inputs:\n")
        f.write("  sweep_raw.csv      : one row per run/repeat (contains phases_json)\n")
        f.write("  sweep_summary.csv  : aggregated stats per (method,variant,t)\n\n")
        f.write("Outputs:\n")
        f.write("  exp2_runtime_vs_t_<variant>.pdf/.png          : wall median vs t (log-x, log-y, upper-only error)\n")
        f.write("  exp2_delta_runtime_vs_t_<variant>.pdf/.png    : Δ wall median vs t (log-x, linear-y)\n")
        f.write("  exp2_sample_phase_vs_t_<variant>.pdf/.png     : sample-phase (median run_sample_ms) vs t (log-x, log-y)\n")
        f.write("  exp2_ns_per_sample.csv                        : slope-based ns/sample from sample-phase medians\n\n")
        f.write("Notes:\n")
        f.write("  - Points with ok_rate < 1.0 are filtered out (e.g., skipped configurations).\n")
        f.write("  - Error bars: upper-only (p95-median) by default; avoids log-y issues.\n")

    print(f"[EXP2][PLOT] wrote {readme_path}")
    print(f"[EXP2][PLOT] DONE. Out dir: {out_dir}")


if __name__ == "__main__":
    main()
