#!/usr/bin/env python3
"""EXP-6 plotting helper (phase diagram + ratio heatmap + optional branch heatmap).

This file exists to keep run/run_exp6.sh free of embedded Python.

Inputs:
  - sweep_summary.csv: produced by sjs_sweep
  - sweep_raw.csv    : produced by sjs_sweep (optional; enables adaptive_branch heatmap)

Outputs (written to --out_dir):
  - exp6_phase_winner_<method>.csv
  - exp6_ratio_<method>.csv
  - exp6_adaptive_enumerate_frac_<method>.csv (if raw csv contains adaptive_branch)
  - exp6_phase_<method>.png (if matplotlib is available)
  - exp6_ratio_<method>.png (if matplotlib is available)
  - exp6_adaptive_branch_<method>.png (if matplotlib is available and raw supports it)

Plotting policy (aligned with EXP-6.md):
  - Only use fully successful points (ok_rate == 1) for phase/ratio.
  - Treat SKIPPED/unsupported points as NA.
  - NA must stay NA (do not silently drop); CSV matrices record NA as blank/NaN and
    winner as -1.

If matplotlib is unavailable, the script writes CSV matrices only and exits 0.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Tuple


def _as_float(x: Any, default: float = math.nan) -> float:
    try:
        if x is None:
            return default
        if isinstance(x, float):
            return x
        s = str(x).strip()
        if s == "":
            return default
        return float(s)
    except Exception:
        return default


def _as_int(x: Any, default: Optional[int] = None) -> Optional[int]:
    try:
        if x is None:
            return default
        s = str(x).strip()
        if s == "":
            return default
        return int(float(s))
    except Exception:
        return default


def _fmt_alpha(a: float) -> str:
    # Avoid "0.30000000000004" style strings.
    try:
        return f"{a:g}"
    except Exception:
        return str(a)


def _read_rows(path: str) -> List[Dict[str, Any]]:
    """Read a CSV file into a list-of-dicts.

    Prefers pandas if available, but keeps a pure-csv fallback.
    """

    try:
        import pandas as pd  # type: ignore

        df = pd.read_csv(path)
        return df.to_dict(orient="records")  # type: ignore[return-value]
    except Exception:
        rows: List[Dict[str, Any]] = []
        with open(path, "r", newline="", encoding="utf-8") as f:
            rd = csv.DictReader(f)
            for r in rd:
                rows.append({k: v for k, v in r.items()})
        return rows


def _write_matrix_csv(path: str, row_labels: List[str], col_labels: List[str], matrix: List[List[Any]]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["alpha\\t"] + col_labels)
        for i, a in enumerate(row_labels):
            w.writerow([a] + list(matrix[i]))


def _unique_sorted(vals: Iterable[Any], key=None) -> List[Any]:
    return sorted(set(vals), key=key)


def _extract_summary_times(
    summary_rows: List[Dict[str, Any]],
    method: str,
    variants: List[str],
) -> Tuple[List[float], List[int], Dict[Tuple[float, int, str], float]]:
    """Return (alphas, ts, time_map).

    time_map[(alpha, t, variant)] = wall_median_ms (only if ok_rate==1 and wall>0 and not SKIPPED)
    otherwise missing entries imply NA.
    """

    sub = [r for r in summary_rows if str(r.get("method", "")) == method]
    if not sub:
        return [], [], {}

    alphas: List[float] = []
    ts: List[int] = []
    time_map: Dict[Tuple[float, int, str], float] = {}

    for r in sub:
        a = _as_float(r.get("alpha"))
        t = _as_int(r.get("t"))
        v = str(r.get("variant", ""))
        ok_rate = _as_float(r.get("ok_rate"), default=1.0)
        wall = _as_float(r.get("wall_median_ms"))
        note = str(r.get("note", "") or "").strip()

        if math.isnan(a) or t is None:
            continue
        alphas.append(a)
        ts.append(t)

        # Only accept fully successful points, as required by EXP-6.md.
        # Also drop non-positive wall times (log-scale safety + guards for SKIPPED rows).
        if ok_rate < 1.0:
            continue
        if wall <= 0 or math.isnan(wall):
            continue
        note_up = note.upper()
        # Treat explicit skips/unsupported as NA, and also treat truncation/cap hits
        # as NA (EXP-6.md requires excluding truncated enum points from oracle(min)).
        if note_up.startswith("SKIPPED") or note_up.startswith("UNSUPPORTED") or ("TRUNC" in note_up) or ("CAP" in note_up):
            continue
        if v not in variants:
            continue

        time_map[(a, t, v)] = wall

    alphas_u = _unique_sorted(alphas, key=float)
    ts_u = _unique_sorted(ts, key=int)
    return alphas_u, ts_u, time_map


def _compute_phase_and_ratio(
    alphas: List[float],
    ts: List[int],
    time_map: Dict[Tuple[float, int, str], float],
    variants: List[str],
) -> Tuple[List[List[int]], List[List[float]]]:
    idx = {v: i for i, v in enumerate(variants)}

    winner: List[List[int]] = [[-1 for _ in ts] for _ in alphas]
    ratio: List[List[float]] = [[math.nan for _ in ts] for _ in alphas]

    for i, a in enumerate(alphas):
        for j, t in enumerate(ts):
            vals = [math.nan] * len(variants)
            for v in variants:
                if (a, t, v) in time_map:
                    vals[idx[v]] = time_map[(a, t, v)]

            # Winner among finite vals
            best_k = -1
            best_val = float("inf")
            for k, x in enumerate(vals):
                if x is None or math.isnan(x):
                    continue
                if x < best_val:
                    best_val = x
                    best_k = k
            winner[i][j] = best_k

            # Ratio: adaptive / min(sampling, enum_sampling) among available baselines
            try:
                adaptive = vals[idx["adaptive"]]
            except Exception:
                adaptive = math.nan

            denom = float("inf")
            for base_v in ("sampling", "enum_sampling"):
                if base_v not in idx:
                    continue
                x = vals[idx[base_v]]
                if x is None or math.isnan(x):
                    continue
                denom = min(denom, x)

            if adaptive is not None and not math.isnan(adaptive) and denom != float("inf") and denom > 0:
                ratio[i][j] = adaptive / denom
            else:
                ratio[i][j] = math.nan

    return winner, ratio


def _try_import_matplotlib():
    try:
        import matplotlib.pyplot as plt  # type: ignore

        return plt
    except Exception as e:
        print(f"[exp6_plot][WARN] matplotlib not available; write CSV only. ({e})")
        return None


def _plot_imshow(
    plt,
    matrix: List[List[float]],
    alphas: List[float],
    ts: List[int],
    title: str,
    out_png: str,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> None:
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)

    # Use numpy if available (better NaN handling), but keep a fallback.
    try:
        import numpy as np  # type: ignore

        arr = np.array(matrix, dtype=float)
        plt.figure()
        if vmin is None and vmax is None:
            im = plt.imshow(arr, aspect="auto", interpolation="nearest")
        else:
            im = plt.imshow(arr, aspect="auto", interpolation="nearest", vmin=vmin, vmax=vmax)
        plt.yticks(range(len(alphas)), [_fmt_alpha(a) for a in alphas])
        plt.xticks(range(len(ts)), [str(t) for t in ts], rotation=45, ha="right")
        plt.xlabel("t")
        plt.ylabel("alpha")
        plt.title(title)
        plt.colorbar(im)
        plt.tight_layout()
        plt.savefig(out_png, dpi=200)
        plt.close()
        return
    except Exception:
        pass

    # Pure-python fallback
    plt.figure()
    if vmin is None and vmax is None:
        im = plt.imshow(matrix, aspect="auto", interpolation="nearest")
    else:
        im = plt.imshow(matrix, aspect="auto", interpolation="nearest", vmin=vmin, vmax=vmax)
    plt.yticks(range(len(alphas)), [_fmt_alpha(a) for a in alphas])
    plt.xticks(range(len(ts)), [str(t) for t in ts], rotation=45, ha="right")
    plt.xlabel("t")
    plt.ylabel("alpha")
    plt.title(title)
    plt.colorbar(im)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def _compute_adaptive_enumerate_fraction(
    raw_rows: List[Dict[str, Any]],
    method: str,
    alphas: List[float],
    ts: List[int],
) -> Optional[List[List[float]]]:
    """Return enumerate_frac matrix for adaptive runs (ok==1), or None if unsupported."""

    # Required columns for branch heatmap.
    required = {"method", "variant", "alpha", "t", "ok", "adaptive_branch"}
    if not raw_rows:
        return None
    if not required.issubset(set(raw_rows[0].keys())):
        return None

    # Count per (alpha,t): ok_total, enum_ok
    ok_total: Dict[Tuple[float, int], int] = defaultdict(int)
    enum_ok: Dict[Tuple[float, int], int] = defaultdict(int)

    for r in raw_rows:
        if str(r.get("method", "")) != method:
            continue
        if str(r.get("variant", "")) != "adaptive":
            continue
        a = _as_float(r.get("alpha"))
        t = _as_int(r.get("t"))
        ok = _as_int(r.get("ok"), default=0) or 0
        if math.isnan(a) or t is None:
            continue
        if ok != 1:
            continue

        ok_total[(a, t)] += 1

        branch = str(r.get("adaptive_branch", "") or "").strip()
        if branch.startswith("enumerate"):
            enum_ok[(a, t)] += 1

    mat: List[List[float]] = [[math.nan for _ in ts] for _ in alphas]
    for i, a in enumerate(alphas):
        for j, t in enumerate(ts):
            tot = ok_total.get((a, t), 0)
            if tot <= 0:
                mat[i][j] = math.nan
            else:
                mat[i][j] = enum_ok.get((a, t), 0) / tot
    return mat


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary_csv", required=True, help="Path to sweep_summary.csv")
    ap.add_argument("--raw_csv", default="", help="Path to sweep_raw.csv (optional; enables branch heatmap)")
    ap.add_argument("--out_dir", required=True, help="Directory to write figures and CSV matrices")
    args = ap.parse_args()

    summary_csv = args.summary_csv
    raw_csv = args.raw_csv
    out_dir = args.out_dir

    if not os.path.isfile(summary_csv):
        print(f"[exp6_plot][ERROR] summary_csv not found: {summary_csv}", file=sys.stderr)
        return 2

    os.makedirs(out_dir, exist_ok=True)

    summary_rows = _read_rows(summary_csv)
    if not summary_rows:
        print(f"[exp6_plot][ERROR] empty summary: {summary_csv}", file=sys.stderr)
        return 3

    # Determine methods present in summary
    methods = _unique_sorted([str(r.get("method", "")) for r in summary_rows if str(r.get("method", ""))])
    if not methods:
        print(f"[exp6_plot][ERROR] no method column found or empty methods in: {summary_csv}", file=sys.stderr)
        return 4

    # EXP-6 fixed comparison set
    variants = ["sampling", "enum_sampling", "adaptive"]

    # Optional raw rows
    raw_rows: List[Dict[str, Any]] = []
    if raw_csv and os.path.isfile(raw_csv):
        try:
            raw_rows = _read_rows(raw_csv)
        except Exception as e:
            print(f"[exp6_plot][WARN] failed to read raw_csv ({raw_csv}): {e}")
            raw_rows = []

    plt = _try_import_matplotlib()

    for method in methods:
        alphas, ts, time_map = _extract_summary_times(summary_rows, method, variants)
        if not alphas or not ts:
            print(f"[exp6_plot][WARN] no usable (alpha,t) grid for method={method}; skip")
            continue

        winner, ratio = _compute_phase_and_ratio(alphas, ts, time_map, variants)

        # CSV matrices (always)
        winner_csv = os.path.join(out_dir, f"exp6_phase_winner_{method}.csv")
        ratio_csv = os.path.join(out_dir, f"exp6_ratio_{method}.csv")
        _write_matrix_csv(winner_csv, [_fmt_alpha(a) for a in alphas], [str(t) for t in ts], winner)
        _write_matrix_csv(ratio_csv, [_fmt_alpha(a) for a in alphas], [str(t) for t in ts], ratio)

        # Optional branch heatmap
        branch_mat = _compute_adaptive_enumerate_fraction(raw_rows, method, alphas, ts)
        branch_csv = ""
        if branch_mat is not None:
            branch_csv = os.path.join(out_dir, f"exp6_adaptive_enumerate_frac_{method}.csv")
            _write_matrix_csv(branch_csv, [_fmt_alpha(a) for a in alphas], [str(t) for t in ts], branch_mat)

        if plt is None:
            msg = f"[exp6_plot] Wrote CSV matrices for method={method}: {winner_csv}, {ratio_csv}"
            if branch_csv:
                msg += f", {branch_csv}"
            print(msg)
            continue

        # Winner plot: convert -1 -> NaN for display
        winner_float: List[List[float]] = []
        for row in winner:
            winner_float.append([float(x) if x >= 0 else math.nan for x in row])

        phase_png = os.path.join(out_dir, f"exp6_phase_{method}.png")
        _plot_imshow(
            plt,
            winner_float,
            alphas,
            ts,
            title=(
                "EXP-6 Phase diagram (winner: 0 sampling / 1 enum_sampling / 2 adaptive) "
                f"— method={method}"
            ),
            out_png=phase_png,
            vmin=-0.5,
            vmax=2.5,
        )

        ratio_png = os.path.join(out_dir, f"exp6_ratio_{method}.png")
        _plot_imshow(
            plt,
            ratio,
            alphas,
            ts,
            title=(
                "EXP-6 ratio = T(adaptive)/min(T(sampling),T(enum_sampling)) "
                f"— method={method}"
            ),
            out_png=ratio_png,
        )

        if branch_mat is not None:
            branch_png = os.path.join(out_dir, f"exp6_adaptive_branch_{method}.png")
            _plot_imshow(
                plt,
                branch_mat,
                alphas,
                ts,
                title=f"EXP-6 adaptive enumerate fraction (ok==1) — method={method}",
                out_png=branch_png,
                vmin=0.0,
                vmax=1.0,
            )

        print(f"[exp6_plot] Wrote PNGs/CSVs for method={method} into: {out_dir}")

    print("[exp6_plot] done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
