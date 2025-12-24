#!/usr/bin/env python3
# ==============================================================================
# run/include/csv_check.py
#
# Small CSV helper for bash-driven tests.
# Purpose: avoid fragile "cut -d, ..." parsing when CSV has quoted JSON columns.
#
# Typical usage:
#   python3 run/include/csv_check.py path/to/run.csv --nonempty
#   python3 run/include/csv_check.py path/to/run.csv --last-row-col-contains adaptive_branch fallback
#   python3 run/include/csv_check.py path/to/run.csv --any-row-col-int-eq enum_truncated 1
#   python3 run/include/csv_check.py path/to/run.csv --any-row-col-int-eq ok 0
#
# Exit code:
#   0 = check passed
#   2 = check failed (assertion)
# ==============================================================================
from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Any, Dict, List, Optional, Tuple


def eprint(*args: Any) -> None:
    print(*args, file=sys.stderr)


def die(msg: str, code: int = 2) -> None:
    eprint(msg)
    raise SystemExit(code)


def read_csv(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    if not os.path.exists(path):
        die(f"[csv_check] File not found: {path}")
    try:
        with open(path, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fieldnames = list(reader.fieldnames or [])
            rows = list(reader)
            return fieldnames, rows
    except UnicodeDecodeError:
        # Best-effort fallback for non-utf8
        with open(path, "r", newline="", encoding="latin-1") as f:
            reader = csv.DictReader(f)
            fieldnames = list(reader.fieldnames or [])
            rows = list(reader)
            return fieldnames, rows
    except Exception as ex:
        die(f"[csv_check] Failed to read CSV: {path}\n  error: {ex}")


def coerce_int(x: str) -> Optional[int]:
    s = (x or "").strip()
    if s == "":
        return None
    # Accept "1", "1.0", "true", "false"
    if s.lower() in ("true", "t", "yes", "y"):
        return 1
    if s.lower() in ("false", "f", "no", "n"):
        return 0
    try:
        if "." in s:
            return int(float(s))
        return int(s)
    except Exception:
        return None


def main() -> int:
    ap = argparse.ArgumentParser(description="CSV check helper (robust to quoted JSON columns).")
    ap.add_argument("csv_path", help="Path to a CSV file.")
    ap.add_argument("--nonempty", action="store_true", help="Assert CSV has at least one data row.")
    ap.add_argument(
        "--require-columns",
        nargs="+",
        default=None,
        metavar="COL",
        help="Assert these columns exist in header.",
    )
    ap.add_argument(
        "--last-row-col-contains",
        nargs=2,
        metavar=("COL", "SUBSTR"),
        help="Assert last row's COL contains SUBSTR (case-insensitive).",
    )
    ap.add_argument(
        "--last-row-col-eq",
        nargs=2,
        metavar=("COL", "VALUE"),
        help="Assert last row's COL equals VALUE (string compare after strip).",
    )
    ap.add_argument(
        "--any-row-col-int-eq",
        nargs=2,
        metavar=("COL", "INT"),
        help="Assert any row has int(COL) == INT (supports true/false).",
    )
    ap.add_argument(
        "--all-rows-col-int-eq",
        nargs=2,
        metavar=("COL", "INT"),
        help="Assert all rows have int(COL) == INT (supports true/false).",
    )
    args = ap.parse_args()

    cols, rows = read_csv(args.csv_path)

    if args.nonempty and len(rows) == 0:
        die(f"[csv_check] CSV is empty (no rows): {args.csv_path}")

    if args.require_columns:
        missing = [c for c in args.require_columns if c not in cols]
        if missing:
            die(f"[csv_check] Missing required columns in {args.csv_path}: {missing}\n  have: {cols}")

    if args.last_row_col_contains is not None:
        col, substr = args.last_row_col_contains
        if col not in cols:
            die(f"[csv_check] Column not found: {col}\n  file: {args.csv_path}\n  have: {cols}")
        if len(rows) == 0:
            die(f"[csv_check] No rows in CSV for last-row check: {args.csv_path}")
        v = (rows[-1].get(col, "") or "").strip()
        if substr.lower() not in v.lower():
            die(
                f"[csv_check] last-row check failed:\n"
                f"  file : {args.csv_path}\n"
                f"  col  : {col}\n"
                f"  want : contains '{substr}'\n"
                f"  got  : '{v}'"
            )

    if args.last_row_col_eq is not None:
        col, want = args.last_row_col_eq
        if col not in cols:
            die(f"[csv_check] Column not found: {col}\n  file: {args.csv_path}\n  have: {cols}")
        if len(rows) == 0:
            die(f"[csv_check] No rows in CSV for last-row check: {args.csv_path}")
        v = (rows[-1].get(col, "") or "").strip()
        if v != want.strip():
            die(
                f"[csv_check] last-row equality failed:\n"
                f"  file : {args.csv_path}\n"
                f"  col  : {col}\n"
                f"  want : '{want.strip()}'\n"
                f"  got  : '{v}'"
            )

    if args.any_row_col_int_eq is not None:
        col, want_s = args.any_row_col_int_eq
        want = coerce_int(want_s)
        if want is None:
            die(f"[csv_check] Invalid INT for --any-row-col-int-eq: {want_s}")
        if col not in cols:
            die(f"[csv_check] Column not found: {col}\n  file: {args.csv_path}\n  have: {cols}")
        ok = False
        for r in rows:
            got = coerce_int(r.get(col, ""))
            if got is not None and got == want:
                ok = True
                break
        if not ok:
            die(
                f"[csv_check] any-row int equality failed:\n"
                f"  file : {args.csv_path}\n"
                f"  col  : {col}\n"
                f"  want : {want}\n"
                f"  note : no row matched"
            )

    if args.all_rows_col_int_eq is not None:
        col, want_s = args.all_rows_col_int_eq
        want = coerce_int(want_s)
        if want is None:
            die(f"[csv_check] Invalid INT for --all-rows-col-int-eq: {want_s}")
        if col not in cols:
            die(f"[csv_check] Column not found: {col}\n  file: {args.csv_path}\n  have: {cols}")
        for idx, r in enumerate(rows):
            got = coerce_int(r.get(col, ""))
            if got is None or got != want:
                die(
                    f"[csv_check] all-rows int equality failed:\n"
                    f"  file : {args.csv_path}\n"
                    f"  row  : {idx}\n"
                    f"  col  : {col}\n"
                    f"  want : {want}\n"
                    f"  got  : {r.get(col, '')!r}"
                )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
