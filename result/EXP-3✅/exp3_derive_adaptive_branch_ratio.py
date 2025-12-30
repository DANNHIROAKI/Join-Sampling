#!/usr/bin/env python3
"""Derive EXP-3 adaptive branch ratios from sweep_raw.csv.

Per EXP-3.md (RQ3):
  - Only consider successful runs: variant == 'adaptive' AND ok == 1
  - enumerate_frac = count(adaptive_branch == 'enumerate_all') / runs
  - fallback_frac  = 1 - enumerate_frac (we also expose the breakdown)

The output CSV is grouped by:
  (dataset, generator, n_r, n_s, method, alpha)

Example:
  python3 run/include/exp3_derive_adaptive_branch_ratio.py \
    --raw run/temp/exp3/sweep_raw.csv \
    --out run/temp/exp3/derived/adaptive_branch_ratio.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict
from typing import Dict, Tuple


def parse_ok(v: str) -> bool:
    v = (v or "").strip().lower()
    return v in ("1", "true", "yes", "y")


def alpha_key(a: str) -> float:
    try:
        return float(a)
    except Exception:
        return math.inf


def main() -> int:
    ap = argparse.ArgumentParser(description="Derive adaptive branch ratios for EXP-3")
    ap.add_argument("--raw", required=True, help="Path to sweep_raw.csv")
    ap.add_argument("--out", required=True, help="Output CSV path")
    args = ap.parse_args()

    raw_path = args.raw
    out_path = args.out

    if not os.path.exists(raw_path):
        raise SystemExit(f"[EXP3][ERROR] sweep_raw.csv not found: {raw_path}")

    # Group adaptive runs by (dataset,generator,n_r,n_s,method,alpha)
    Key = Tuple[str, str, str, str, str, str]
    grp: Dict[Key, Dict[str, int]] = defaultdict(
        lambda: {
            "total": 0,
            "enumerate_all": 0,
            "fallback_sampling": 0,
            "fallback_sampling_no_pilot": 0,
            "other": 0,
        }
    )

    with open(raw_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {
            "variant",
            "method",
            "alpha",
            "adaptive_branch",
            "ok",
            "dataset",
            "generator",
            "n_r",
            "n_s",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"[EXP3][ERROR] sweep_raw.csv missing columns: {sorted(missing)}")

        for row in reader:
            if (row.get("variant") or "") != "adaptive":
                continue
            if not parse_ok(row.get("ok", "")):
                continue

            key: Key = (
                row.get("dataset", ""),
                row.get("generator", ""),
                row.get("n_r", ""),
                row.get("n_s", ""),
                row.get("method", ""),
                row.get("alpha", ""),
            )

            g = grp[key]
            g["total"] += 1

            b = (row.get("adaptive_branch") or "").strip()
            if b == "enumerate_all":
                g["enumerate_all"] += 1
            elif b == "fallback_sampling":
                g["fallback_sampling"] += 1
            elif b == "fallback_sampling_no_pilot":
                g["fallback_sampling_no_pilot"] += 1
            else:
                g["other"] += 1

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)

    fieldnames = [
        "dataset",
        "generator",
        "n_r",
        "n_s",
        "method",
        "alpha",
        "runs",
        "enumerate_all",
        "fallback_sampling",
        "fallback_sampling_no_pilot",
        "other",
        "enumerate_frac",
        "fallback_frac",
        "fallback_no_pilot_frac",
        "other_frac",
        "fallback_total_frac",
    ]

    rows = []
    for (dataset, generator, n_r, n_s, method, alpha), g in grp.items():
        total = g["total"] or 1
        enum_f = g["enumerate_all"] / total
        fb_f = g["fallback_sampling"] / total
        fb_np_f = g["fallback_sampling_no_pilot"] / total
        other_f = g["other"] / total
        rows.append(
            {
                "dataset": dataset,
                "generator": generator,
                "n_r": n_r,
                "n_s": n_s,
                "method": method,
                "alpha": alpha,
                "runs": g["total"],
                "enumerate_all": g["enumerate_all"],
                "fallback_sampling": g["fallback_sampling"],
                "fallback_sampling_no_pilot": g["fallback_sampling_no_pilot"],
                "other": g["other"],
                "enumerate_frac": enum_f,
                "fallback_frac": fb_f,
                "fallback_no_pilot_frac": fb_np_f,
                "other_frac": other_f,
                "fallback_total_frac": fb_f + fb_np_f,
            }
        )

    rows.sort(key=lambda r: (r["method"], alpha_key(r["alpha"])))

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print("[EXP3] Wrote adaptive branch ratio:", out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
