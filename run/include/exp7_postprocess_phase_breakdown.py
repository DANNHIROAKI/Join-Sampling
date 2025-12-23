#!/usr/bin/env python3
"""Post-process EXP-7 sweep_raw.csv into a phase breakdown table + plots.

Inputs:
  - One or two sweep_raw.csv files (sampling+adaptive; optionally enum_sampling)

Outputs:
  - merged raw CSV (for convenience)
  - breakdown CSV (median over repeats) with stage times and stage percentages
  - breakdown Markdown table (quick paste into README/paper)
  - stacked-bar figures per alpha (requires matplotlib; optional)

The breakdown follows EXP-7.md's mapping rules:
  sampling:
    Build = run_build_ms
    Count = run_count_ms
    Sample = run_sample_ms
    Enumerate = 0

  enum_sampling:
    Build = run_build_ms
    Enumerate = run_enum_prepare_ms + run_enum_pass1_count_ms
    Sample = run_draw_ranks_ms + run_enum_pass2_select_ms
    Count = 0

  adaptive:
    Build = run_build_ms
    Enumerate(pilot) = run_pilot_enum_prepare_ms + run_pilot_enum_scan_ms
    Count(fallback) = run_fallback_count_ms
    Sample = run_fallback_sample_ms + run_small_join_sample_from_list_ms

For percentages, total = sum(phases_json[run_*_ms]) by default.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


def _fnum(x: Any, default: float = 0.0) -> float:
    try:
        return float(x)
    except Exception:
        return float(default)


def _parse_phases(s: str) -> Dict[str, Any]:
    if not s:
        return {}
    try:
        obj = json.loads(s)
        return obj if isinstance(obj, dict) else {}
    except Exception:
        return {}


def _stage_breakdown(variant: str, phases: Dict[str, Any]) -> Tuple[float, float, float, float, float]:
    """Return (build_ms, count_ms, enum_ms, sample_ms, total_ms_from_phases)."""

    g = lambda k: _fnum(phases.get(k, 0.0), 0.0)

    build = g("run_build_ms")

    if variant == "sampling":
        count = g("run_count_ms")
        enum = 0.0
        sample = g("run_sample_ms")

    elif variant == "enum_sampling":
        # enum_sampling runner phases:
        # run_enum_prepare_ms, run_enum_pass1_count_ms, run_draw_ranks_ms, run_enum_pass2_select_ms
        count = 0.0
        enum = g("run_enum_prepare_ms") + g("run_enum_pass1_count_ms")
        sample = g("run_draw_ranks_ms") + g("run_enum_pass2_select_ms")

    elif variant == "adaptive":
        # adaptive runner phases:
        # run_pilot_enum_prepare_ms, run_pilot_enum_scan_ms,
        # run_small_join_sample_from_list_ms (if small-join),
        # run_fallback_count_ms, run_fallback_sample_ms (if fallback)
        enum = g("run_pilot_enum_prepare_ms") + g("run_pilot_enum_scan_ms")
        count = g("run_fallback_count_ms")
        sample = g("run_fallback_sample_ms") + g("run_small_join_sample_from_list_ms")

    else:
        # best-effort for unknown variants
        count = 0.0
        enum = 0.0
        sample = 0.0

    # Total = sum of all recorded run_*_ms phases (preferred; sums to 100%).
    total = 0.0
    for k, v in phases.items():
        if isinstance(k, str) and k.startswith("run_") and k.endswith("_ms"):
            total += _fnum(v, 0.0)

    if total <= 0.0:
        total = build + count + enum + sample

    return build, count, enum, sample, total


def _median(vals: Sequence[float]) -> float:
    cleaned = [float(v) for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    return statistics.median(cleaned) if cleaned else 0.0


def _safe_float_key(x: Any) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description="EXP-7 postprocess: phase breakdown")
    ap.add_argument("--raw_sa", required=True, help="Path to sampling_adaptive/sweep_raw.csv")
    ap.add_argument("--raw_enum", default="", help="Optional path to enum_sparse/sweep_raw.csv")

    ap.add_argument("--out_breakdown", required=True, help="Output breakdown CSV path")
    ap.add_argument("--out_md", required=True, help="Output breakdown Markdown path")
    ap.add_argument("--fig_dir", required=True, help="Directory for plots (optional; requires matplotlib)")
    ap.add_argument("--merged_raw", required=True, help="Output merged raw CSV path")

    ap.add_argument(
        "--exclude_enum_truncated",
        type=int,
        default=1,
        help="Exclude rows where enum_truncated==1 from the breakdown (default: 1)",
    )

    args = ap.parse_args()

    raw_files: List[str] = []
    raw_sa = Path(args.raw_sa)
    if raw_sa.exists():
        raw_files.append(str(raw_sa))
    else:
        raise SystemExit(f"[EXP7][ERROR] Missing required raw file: {raw_sa}")

    raw_enum = Path(args.raw_enum) if args.raw_enum else None
    if raw_enum and raw_enum.exists():
        raw_files.append(str(raw_enum))

    out_breakdown = Path(args.out_breakdown)
    out_md = Path(args.out_md)
    fig_dir = Path(args.fig_dir)
    merged_raw = Path(args.merged_raw)

    out_breakdown.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    merged_raw.parent.mkdir(parents=True, exist_ok=True)

    # ----------------------
    # Read & merge rows
    # ----------------------
    rows: List[Dict[str, Any]] = []
    for path in raw_files:
        with open(path, "r", newline="", encoding="utf-8") as f:
            rdr = csv.DictReader(f)
            for r in rdr:
                r = dict(r)
                r["_raw_file"] = path
                rows.append(r)

    if not rows:
        raise SystemExit("[EXP7][ERROR] No rows found in input sweep_raw.csv files")

    # Write merged raw (union of fieldnames for robustness).
    all_fields: List[str] = sorted({k for r in rows for k in r.keys()})
    with open(merged_raw, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=all_fields)
        w.writeheader()
        w.writerows(rows)

    # Filter: ok==1; optionally exclude enum_truncated==1.
    def is_ok(r: Dict[str, Any]) -> bool:
        ok = str(r.get("ok", "0")).strip()
        if ok != "1":
            return False
        if args.exclude_enum_truncated:
            if str(r.get("enum_truncated", "0")).strip() == "1":
                return False
        return True

    ok_rows = [r for r in rows if is_ok(r)]

    # Compute stage times.
    for r in ok_rows:
        phases = _parse_phases(str(r.get("phases_json", "")))
        variant = str(r.get("variant", ""))

        b, c, e, s, tot = _stage_breakdown(variant, phases)
        r["build_ms"] = b
        r["count_ms"] = c
        r["enum_ms"] = e
        r["sample_ms"] = s
        r["run_total_ms_from_phases"] = tot
        r["wall_ms"] = _fnum(r.get("wall_ms", 0.0), 0.0)

    # Group by (dataset,generator,alpha,n_r,n_s,method,variant,t)
    gcols = ["dataset", "generator", "alpha", "n_r", "n_s", "method", "variant", "t"]
    groups: Dict[Tuple[str, ...], List[Dict[str, Any]]] = {}
    for r in ok_rows:
        key = tuple(str(r.get(k, "")) for k in gcols)
        groups.setdefault(key, []).append(r)

    records: List[Dict[str, Any]] = []
    for key, rs in groups.items():
        rec: Dict[str, Any] = {k: v for k, v in zip(gcols, key)}
        rec["repeats_ok"] = len(rs)

        for m in ["build_ms", "count_ms", "enum_ms", "sample_ms", "run_total_ms_from_phases", "wall_ms"]:
            rec[m] = _median([_fnum(r.get(m, 0.0), 0.0) for r in rs])

        total = rec["run_total_ms_from_phases"]
        if total <= 0.0:
            total = rec["build_ms"] + rec["count_ms"] + rec["enum_ms"] + rec["sample_ms"]

        rec["build_pct"] = (rec["build_ms"] / total * 100.0) if total > 0 else 0.0
        rec["count_pct"] = (rec["count_ms"] / total * 100.0) if total > 0 else 0.0
        rec["enum_pct"] = (rec["enum_ms"] / total * 100.0) if total > 0 else 0.0
        rec["sample_pct"] = (rec["sample_ms"] / total * 100.0) if total > 0 else 0.0

        records.append(rec)

    # Sort for readability: alpha, variant, method
    records.sort(
        key=lambda r: (
            _safe_float_key(r.get("alpha", "nan")),
            r.get("variant", ""),
            r.get("method", ""),
        )
    )

    # Write breakdown CSV
    fieldnames = gcols + [
        "repeats_ok",
        "build_ms",
        "count_ms",
        "enum_ms",
        "sample_ms",
        "run_total_ms_from_phases",
        "wall_ms",
        "build_pct",
        "count_pct",
        "enum_pct",
        "sample_pct",
    ]

    with open(out_breakdown, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(records)

    # Write markdown table
    lines: List[str] = []
    lines.append(
        "| alpha | method | variant | Build% | Count% | Enum% | Sample% | total_ms(phases) | wall_ms |"
    )
    lines.append("|---:|---|---|---:|---:|---:|---:|---:|---:|")
    for r in records:
        lines.append(
            f"| {r['alpha']} | {r['method']} | {r['variant']} | "
            f"{r['build_pct']:.1f} | {r['count_pct']:.1f} | {r['enum_pct']:.1f} | {r['sample_pct']:.1f} | "
            f"{r['run_total_ms_from_phases']:.1f} | {r['wall_ms']:.1f} |"
        )

    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[EXP7] Wrote merged raw   : {merged_raw}")
    print(f"[EXP7] Wrote breakdown CSV: {out_breakdown}")
    print(f"[EXP7] Wrote breakdown MD : {out_md}")

    # ----------------------
    # Optional plotting
    # ----------------------
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as e:
        print(f"[EXP7] matplotlib not available -> skipping plots ({e})")
        return 0

    # Bucket by alpha
    by_alpha: Dict[str, List[Dict[str, Any]]] = {}
    for r in records:
        by_alpha.setdefault(str(r.get("alpha", "")), []).append(r)

    def _sanitize_alpha(a: str) -> str:
        # Keep filenames reasonably portable.
        return a.replace("/", "_").replace(" ", "").replace("+", "plus")

    for alpha, rs in sorted(by_alpha.items(), key=lambda kv: _safe_float_key(kv[0])):
        rs = sorted(rs, key=lambda x: (x.get("variant", ""), x.get("method", "")))

        labels = [f"{x['method']}\n({x['variant']})" for x in rs]
        build = [float(x["build_pct"]) for x in rs]
        count = [float(x["count_pct"]) for x in rs]
        enum = [float(x["enum_pct"]) for x in rs]
        sample = [float(x["sample_pct"]) for x in rs]

        width = max(10.0, 0.7 * len(labels))
        plt.figure(figsize=(width, 5.2))
        x = list(range(len(labels)))

        plt.bar(x, build, label="Build")
        bottom = build[:]
        plt.bar(x, count, bottom=bottom, label="Count")
        bottom = [bottom[i] + count[i] for i in range(len(bottom))]
        plt.bar(x, enum, bottom=bottom, label="Enumerate")
        bottom = [bottom[i] + enum[i] for i in range(len(bottom))]
        plt.bar(x, sample, bottom=bottom, label="Sample")

        plt.xticks(x, labels, rotation=0)
        plt.ylim(0, 100)
        plt.ylabel("Percent of total (sum of run_* phases)")
        plt.title(f"EXP-7 Phase Breakdown (alpha={alpha})")
        plt.legend()
        plt.tight_layout()

        out_png = fig_dir / f"exp7_phase_breakdown_alpha_{_sanitize_alpha(alpha)}.png"
        plt.savefig(out_png, dpi=200)
        plt.close()
        print(f"[EXP7] Wrote plot: {out_png}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
