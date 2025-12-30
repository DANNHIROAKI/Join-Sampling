#!/usr/bin/env python3
"""Plot EXP-3 (RQ3) results: runtime vs alpha + adaptive branch ratio.

Inputs:
  - sweep_summary.csv (aggregated: alpha/method/variant -> wall_median_ms, wall_p95_ms, ok_rate)
  - adaptive_branch_ratio.csv (derived from sweep_raw.csv)

Per EXP-3.md:
  - x axis: alpha, use symlog (alpha includes 0)
  - y axis: wall_median_ms, use log-y
  - MUST filter points where ok_rate==0 or wall_median_ms<=0 (log-y safety)
  - Also output a CSV listing skipped points (so we don't silently hide them).

The script gracefully skips plotting if matplotlib is not available.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict
from typing import DefaultDict, Dict, List, Optional, Tuple


def try_float(x: str) -> Optional[float]:
    try:
        v = float(x)
        if math.isfinite(v):
            return v
        return None
    except Exception:
        return None


def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot EXP-3 results")
    ap.add_argument("--summary", required=True, help="Path to sweep_summary.csv")
    ap.add_argument("--branch", required=True, help="Path to adaptive_branch_ratio.csv")
    ap.add_argument("--out_dir", required=True, help="Output directory for plots")
    ap.add_argument(
        "--include_enum_sampling",
        action="store_true",
        help="Also include enum_sampling curves in per-method plots",
    )
    ap.add_argument(
        "--linthresh",
        type=float,
        default=0.03,
        help="symlog linthresh for alpha axis (default: 0.03)",
    )
    args = ap.parse_args()

    summary_path = args.summary
    branch_path = args.branch
    out_dir = args.out_dir
    ensure_dir(out_dir)

    # Try to import matplotlib; if unavailable, just exit gracefully.
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as e:
        print(f"[EXP3] matplotlib not available; skip plotting. ({e})")
        return 0

    if not os.path.exists(summary_path):
        raise SystemExit(f"[EXP3][ERROR] Missing sweep_summary.csv: {summary_path}")
    if not os.path.exists(branch_path):
        raise SystemExit(f"[EXP3][ERROR] Missing adaptive_branch_ratio.csv: {branch_path}")

    variants = ["sampling", "adaptive"]
    if args.include_enum_sampling:
        variants.append("enum_sampling")

    # ---------- Read summary into series ----------
    # series[(method, variant)] -> list[(alpha, wall_median_ms, wall_p95_ms, ok_rate)]
    series: DefaultDict[Tuple[str, str], List[Tuple[float, float, float, float]]] = defaultdict(list)
    skipped: List[Dict[str, str]] = []

    with open(summary_path, "r", newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        required = {"alpha", "method", "variant", "wall_median_ms", "wall_p95_ms", "ok_rate"}
        missing = required - set(rd.fieldnames or [])
        if missing:
            raise SystemExit(f"[EXP3][ERROR] sweep_summary.csv missing columns: {sorted(missing)}")

        for row in rd:
            method = (row.get("method") or "").strip()
            variant = (row.get("variant") or "").strip()
            if not method or not variant:
                continue

            a = try_float(row.get("alpha", ""))
            med = try_float(row.get("wall_median_ms", ""))
            p95 = try_float(row.get("wall_p95_ms", ""))
            okr = try_float(row.get("ok_rate", ""))

            if a is None or med is None or p95 is None or okr is None:
                continue

            # Record skipped points for log-y safety.
            reason = None
            if okr <= 0:
                reason = "ok_rate==0"
            elif med <= 0:
                reason = "wall_median_ms<=0"

            if reason is not None:
                skipped.append(
                    {
                        "method": method,
                        "variant": variant,
                        "alpha": str(a),
                        "ok_rate": str(okr),
                        "wall_median_ms": str(med),
                        "wall_p95_ms": str(p95),
                        "reason": reason,
                    }
                )
                # Still keep it out of series (cannot be plotted on log-y)
                continue

            series[(method, variant)].append((a, med, p95, okr))

    # Sort points by alpha
    for k in list(series.keys()):
        series[k].sort(key=lambda x: x[0])

    methods = sorted({m for (m, _v) in series.keys()})

    # Save skipped points CSV (do not silently hide failures/skips)
    if skipped:
        skipped_csv = os.path.join(out_dir, "plot_skipped_points.csv")
        with open(skipped_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(
                f,
                fieldnames=[
                    "method",
                    "variant",
                    "alpha",
                    "ok_rate",
                    "wall_median_ms",
                    "wall_p95_ms",
                    "reason",
                ],
            )
            w.writeheader()
            w.writerows(skipped)
        print("[EXP3] Wrote:", skipped_csv)

    # ---------- Plot 1: runtime vs alpha (per method) ----------
    for m in methods:
        plt.figure()
        any_line = False
        for v in variants:
            pts = series.get((m, v))
            if not pts:
                continue
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            plt.plot(xs, ys, marker="o", label=v)
            any_line = True

        if not any_line:
            plt.close()
            continue

        plt.xscale("symlog", linthresh=args.linthresh)
        plt.yscale("log")
        plt.xlabel("alpha = |J|/(n_r+n_s)")
        plt.ylabel("wall_median_ms (log)")
        plt.title(f"EXP-3 Runtime vs alpha — {m}")
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.legend()
        out = os.path.join(out_dir, f"runtime_vs_alpha_{m}.png")
        plt.tight_layout()
        plt.savefig(out, dpi=200)
        plt.close()
        print("[EXP3] Wrote:", out)

    # ---------- Plot 1b: all methods (sampling) ----------
    plt.figure(figsize=(12, 7))
    any_line = False
    for (m, v), pts in sorted(series.items()):
        if v != "sampling":
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        if not xs:
            continue
        plt.plot(xs, ys, marker="o", label=m)
        any_line = True

    if any_line:
        plt.xscale("symlog", linthresh=args.linthresh)
        plt.yscale("log")
        plt.xlabel("alpha = |J|/(n_r+n_s)")
        plt.ylabel("wall_median_ms (log)")
        plt.title("EXP-3 Runtime vs alpha — sampling (all methods)")
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.legend(ncol=2, fontsize=8)
        out = os.path.join(out_dir, "runtime_vs_alpha_all_methods_sampling.png")
        plt.tight_layout()
        plt.savefig(out, dpi=200)
        plt.close()
        print("[EXP3] Wrote:", out)
    else:
        plt.close()

    # ---------- Plot 2: adaptive branch ratio (enumerate fraction) ----------
    branch: DefaultDict[str, List[Tuple[float, float]]] = defaultdict(list)  # method -> [(alpha, enumerate_frac)]

    with open(branch_path, "r", newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        required = {"method", "alpha", "enumerate_frac"}
        missing = required - set(rd.fieldnames or [])
        if missing:
            raise SystemExit(
                f"[EXP3][ERROR] adaptive_branch_ratio.csv missing columns: {sorted(missing)}"
            )

        for row in rd:
            method = (row.get("method") or "").strip()
            a = try_float(row.get("alpha", ""))
            enumf = try_float(row.get("enumerate_frac", ""))
            if not method or a is None or enumf is None:
                continue
            branch[method].append((a, enumf))

    plt.figure(figsize=(12, 6))
    any_line = False
    for m, pts in sorted(branch.items()):
        pts.sort(key=lambda x: x[0])
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        if not xs:
            continue
        plt.plot(xs, ys, marker="o", label=m)
        any_line = True

    if any_line:
        plt.xscale("symlog", linthresh=args.linthresh)
        plt.ylim(-0.05, 1.05)
        plt.xlabel("alpha = |J|/(n_r+n_s)")
        plt.ylabel("adaptive enumerate fraction")
        plt.title("EXP-3 Adaptive branch ratio vs alpha (enumerate fraction)")
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.legend(ncol=2, fontsize=8)
        out = os.path.join(out_dir, "adaptive_branch_ratio_enumerate_frac.png")
        plt.tight_layout()
        plt.savefig(out, dpi=200)
        plt.close()
        print("[EXP3] Wrote:", out)
    else:
        plt.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
