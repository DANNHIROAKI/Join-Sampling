#!/usr/bin/env python3
"""EXP-1 log parser: sjs_verify logs -> CSV tables.

This script is intentionally standalone (stdlib-only) so it can run in any
artifact/CI environment.

Inputs
------
  --log_dir : directory containing per-(method__variant).log files produced by sjs_verify
  --out_dir : output directory for CSV tables (and where FAILURES.txt is expected)

Outputs (written to --out_dir)
-----------------------------
  exp1_quality_raw.csv
  exp1_quality_summary.csv

The parsing logic matches the output format used by the repository's sjs_verify.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import statistics
from typing import Any, Dict, List, Optional


def _median(values: List[Any]) -> Any:
    vals: List[float] = []
    for v in values:
        if v is None or v == "":
            continue
        try:
            fv = float(v)
        except Exception:
            continue
        if fv != fv:  # NaN
            continue
        vals.append(fv)
    return statistics.median(vals) if vals else ""


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--log_dir", required=True, help="Directory containing *.log files")
    ap.add_argument("--out_dir", required=True, help="Directory to write CSV tables")
    ap.add_argument(
        "--fail_file",
        default=None,
        help="Optional FAILURES.txt path (defaults to <out_dir>/FAILURES.txt)",
    )
    args = ap.parse_args(argv)

    log_dir = args.log_dir
    out_dir = args.out_dir
    fail_path = args.fail_file or os.path.join(out_dir, "FAILURES.txt")

    logs = sorted(glob.glob(os.path.join(log_dir, "*.log")))

    # Regex patterns that match sjs_verify output
    re_run = re.compile(r"^---- run rep=(\d+)\s+seed=(\d+)\s+----$")
    re_method = re.compile(r"^method=([^\s]+)\s+variant=([^\s]+)\s+t=(\d+)\s*$")
    re_count = re.compile(
        r"^count=([0-9.eE+-]+)\s+\((exact|est)\)\s+oracle=([0-9.eE+-]+)\s+rel_err=([0-9.eE+-]+)\s*$"
    )
    re_samples = re.compile(r"^samples=(\d+)\s*$")
    re_failed = re.compile(r"^FAILED:\s*(.*)$")
    re_quality_skipped = re.compile(r"^quality:\s*skipped\b")
    re_missing = re.compile(r"^\s+missing_in_universe=([0-9.eE+-]+)\s*$")
    re_chi2 = re.compile(r"^\s+chi2_stat=.*\s+p_value=([0-9.eE+-]+)\s*$")
    re_ac1 = re.compile(r"^\s+autocorr_hash_lag1=([0-9.eE+-]+)\s*$")
    re_ks = re.compile(r"^\s+ks_hash_uniform01\s+D=.*\s+p=([0-9.eE+-]+)\s*$")

    rows: List[Dict[str, Any]] = []

    def push_row(cur: Optional[Dict[str, Any]]) -> None:
        if not cur:
            return
        cur.setdefault("status", "OK")
        rows.append(cur)

    # Parse each log file
    for path in logs:
        base = os.path.basename(path)

        # filename: <method>__<variant>.log
        m_guess, v_guess = "", ""
        if "__" in base:
            m_guess = base.split("__", 1)[0]
            v_guess = base.split("__", 1)[1].rsplit(".", 1)[0]

        cur: Optional[Dict[str, Any]] = None
        saw_any_run = False

        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.rstrip("\n")

                m = re_run.match(line)
                if m:
                    saw_any_run = True
                    # flush previous run record
                    push_row(cur)
                    cur = {
                        "method": m_guess,
                        "variant": v_guess,
                        "rep": int(m.group(1)),
                        "seed": int(m.group(2)),
                        "log_file": base,
                        "status": "OK",
                    }
                    continue

                if cur is None:
                    continue  # ignore prologue before first run

                m = re_failed.match(line)
                if m:
                    cur["status"] = "FAILED"
                    cur["error"] = m.group(1).strip()
                    continue

                m = re_method.match(line)
                if m:
                    cur["method"] = m.group(1)
                    cur["variant"] = m.group(2)
                    cur["t"] = int(m.group(3))
                    continue

                m = re_count.match(line)
                if m:
                    cur["count"] = float(m.group(1))
                    cur["count_kind"] = m.group(2)
                    cur["oracle"] = float(m.group(3))
                    cur["rel_err"] = float(m.group(4))
                    continue

                m = re_samples.match(line)
                if m:
                    cur["samples"] = int(m.group(1))
                    continue

                if re_quality_skipped.match(line):
                    # Oracle universe likely truncated; sjs_verify marks it explicitly.
                    cur["status"] = "QUALITY_SKIPPED"
                    continue

                m = re_missing.match(line)
                if m:
                    cur["missing_in_universe"] = float(m.group(1))
                    continue

                m = re_chi2.match(line)
                if m:
                    cur["chi2_p"] = float(m.group(1))
                    continue

                m = re_ac1.match(line)
                if m:
                    cur["autocorr_lag1"] = float(m.group(1))
                    continue

                m = re_ks.match(line)
                if m:
                    cur["ks_p"] = float(m.group(1))
                    continue

        push_row(cur)

        # If a log has no "---- run ..." blocks at all, keep a placeholder row.
        if not saw_any_run:
            rows.append(
                {
                    "method": m_guess,
                    "variant": v_guess,
                    "rep": "",
                    "seed": "",
                    "log_file": base,
                    "status": "NO_RUN_BLOCKS",
                }
            )

    # Incorporate top-level failures (non-zero exit codes) so they appear in the CSV too.
    if os.path.exists(fail_path):
        with open(fail_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                m = re.search(
                    r"rc=(\d+)\s+method=([^\s]+)\s+variant=([^\s]+)\s+log=(.+)$",
                    line,
                )
                if not m:
                    continue
                rc = int(m.group(1))
                method = m.group(2)
                variant = m.group(3)
                log = os.path.basename(m.group(4))
                rows.append(
                    {
                        "method": method,
                        "variant": variant,
                        "rep": "",
                        "seed": "",
                        "log_file": log,
                        "status": "FAILED_CREATE_OR_ARGS",
                        "exit_code": rc,
                    }
                )

    os.makedirs(out_dir, exist_ok=True)

    # Write raw table
    raw_csv = os.path.join(out_dir, "exp1_quality_raw.csv")
    all_fields = sorted({k for r in rows for k in r.keys()})
    with open(raw_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=all_fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # Summary (median by method×variant)
    summary_csv = os.path.join(out_dir, "exp1_quality_summary.csv")
    keys = sorted(set((r.get("method", ""), r.get("variant", "")) for r in rows))

    with open(summary_csv, "w", newline="", encoding="utf-8") as f:
        fields = [
            "method",
            "variant",
            "n_rows",
            "n_ok",
            "n_failed",
            "n_failed_create_or_args",
            "n_quality_skipped",
            "n_no_run_blocks",
            "t_median",
            "oracle_median",
            "rel_err_median",
            "missing_median",
            "chi2_p_median",
            "ks_p_median",
            "autocorr_median",
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()

        for method, variant in keys:
            grp = [
                r for r in rows if r.get("method") == method and r.get("variant") == variant
            ]
            n_rows = len(grp)
            n_ok = sum(1 for r in grp if r.get("status") == "OK")
            n_failed = sum(1 for r in grp if r.get("status") == "FAILED")
            n_fc = sum(1 for r in grp if r.get("status") == "FAILED_CREATE_OR_ARGS")
            n_qs = sum(1 for r in grp if r.get("status") == "QUALITY_SKIPPED")
            n_nrb = sum(1 for r in grp if r.get("status") == "NO_RUN_BLOCKS")

            w.writerow(
                {
                    "method": method,
                    "variant": variant,
                    "n_rows": n_rows,
                    "n_ok": n_ok,
                    "n_failed": n_failed,
                    "n_failed_create_or_args": n_fc,
                    "n_quality_skipped": n_qs,
                    "n_no_run_blocks": n_nrb,
                    "t_median": _median([r.get("t") for r in grp]),
                    "oracle_median": _median([r.get("oracle") for r in grp]),
                    "rel_err_median": _median([r.get("rel_err") for r in grp]),
                    "missing_median": _median([r.get("missing_in_universe") for r in grp]),
                    "chi2_p_median": _median([r.get("chi2_p") for r in grp]),
                    "ks_p_median": _median([r.get("ks_p") for r in grp]),
                    "autocorr_median": _median([r.get("autocorr_lag1") for r in grp]),
                }
            )

    print("[exp1_parse_verify_logs] Wrote:", raw_csv)
    print("[exp1_parse_verify_logs] Wrote:", summary_csv)

    # Quick sanity: missing should be 0 when quality is computed.
    bad_missing = [
        r
        for r in rows
        if r.get("status") == "OK" and float(r.get("missing_in_universe", 0.0) or 0.0) != 0.0
    ]
    if bad_missing:
        print(
            "[exp1_parse_verify_logs] WARNING: missing_in_universe != 0 found in",
            len(bad_missing),
            "OK runs (correctness failure).",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
