#!/usr/bin/env python3
"""EXP-5 plotting helper (runtime vs t).

This file exists to keep run/run_exp5.sh free of embedded Python.

Inputs:
  - sweep_summary.csv produced by sjs_sweep

Output:
  - runtime_t.png (log-x log-y), one curve per (method, variant)

Filtering policy (aligned with EXP-5.md):
  - keep ok_rate == 1.0 (fully successful points)
  - drop rows with note starting with 'SKIPPED'
  - drop non-positive t or wall_median_ms (log-scale safety)

If matplotlib is unavailable, the script prints a warning and exits 0.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Any


def _read_rows(summary_csv: str) -> List[Dict[str, str]]:
    # Prefer pandas if available, but keep a pure-csv fallback.
    try:
        import pandas as pd  # type: ignore

        df = pd.read_csv(summary_csv)
        return df.to_dict(orient="records")  # type: ignore[return-value]
    except Exception:
        rows: List[Dict[str, str]] = []
        with open(summary_csv, "r", newline="", encoding="utf-8") as f:
            rd = csv.DictReader(f)
            for r in rd:
                rows.append({k: (v if v is not None else "") for k, v in r.items()})
        return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary_csv", required=True, help="Path to sweep_summary.csv")
    ap.add_argument("--out_png", required=True, help="Output PNG path")
    ap.add_argument("--title", default="EXP-5 runtime vs t", help="Plot title")
    args = ap.parse_args()

    summary_csv = args.summary_csv
    out_png = args.out_png

    if not os.path.isfile(summary_csv):
        print(f"[exp5_plot][ERROR] summary_csv not found: {summary_csv}", file=sys.stderr)
        return 2

    # Matplotlib is optional; skip gracefully if missing.
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as e:
        print(f"[exp5_plot][WARN] matplotlib not available; skip plot. ({e})")
        return 0

    rows = _read_rows(summary_csv)

    # Required columns (best effort; if missing, we still try).
    needed = {"method", "variant", "t", "wall_median_ms"}
    missing = needed - set(rows[0].keys() if rows else [])
    if rows and missing:
        print(f"[exp5_plot][WARN] missing columns in summary: {sorted(missing)}")

    def as_float(x: Any) -> float:
        try:
            return float(x)
        except Exception:
            return float("nan")

    def as_str(x: Any) -> str:
        return "" if x is None else str(x)

    filtered: List[Dict[str, Any]] = []
    for r in rows:
        ok_rate = as_float(r.get("ok_rate", "1"))
        note = as_str(r.get("note", "")).strip()
        t = as_float(r.get("t", "nan"))
        wall = as_float(r.get("wall_median_ms", "nan"))

        # Filter rules (aligned with EXP-5.md recommendations)
        if ok_rate != ok_rate:
            ok_rate = 1.0
        if ok_rate < 1.0:
            continue
        if note.startswith("SKIPPED"):
            continue
        if not (t > 0.0 and wall > 0.0):
            continue

        filtered.append(
            {
                "method": as_str(r.get("method", "")),
                "variant": as_str(r.get("variant", "")),
                "t": t,
                "wall_median_ms": wall,
            }
        )

    if not filtered:
        print(f"[exp5_plot][WARN] empty dataframe after filtering; no plot written: {out_png}")
        return 0

    series: Dict[Tuple[str, str], List[Tuple[float, float]]] = defaultdict(list)
    for r in filtered:
        key = (r["method"], r["variant"])
        series[key].append((r["t"], r["wall_median_ms"]))

    for key in list(series.keys()):
        series[key].sort(key=lambda p: p[0])

    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)

    plt.figure()
    for (method, variant), pts in sorted(series.items()):
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        plt.plot(xs, ys, marker="o", label=f"{method}-{variant}")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("t (#samples)")
    plt.ylabel("median wall time (ms)")
    plt.title(args.title)
    plt.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"[exp5_plot] Wrote: {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
