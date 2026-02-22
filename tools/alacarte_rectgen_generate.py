#!/usr/bin/env python3
"""tools/alacarte_rectgen_generate.py

Generate controllable output-density synthetic axis-aligned hyper-rectangles
using the `alacarte-rectgen` package, then export to the SJS-HighDims binary
format (SJSBOX v1).

This is meant to be the *single source of truth* for synthetic dataset
generation in this repo when you choose:

  --gen=alacarte_rectgen

It matches the fixed parameterization described in your spec:

  - universe=None          (default [0,1)^d)
  - volume_dist="normal"
  - volume_cv=0.25
  - shape_sigma=0.5
  - tune_tol_rel=0.01
  - dtype=float32

It also performs a generation quality audit:
  - records info["alpha_expected_est"] and info["coverage"]
  - estimates realized alpha via pair sampling

Outputs
-------
Writes two binary relation files (R and S) in SJSBOX v1 format.
Optionally writes CSV files for debugging.
Always writes a JSON report (if --report_path is provided).

Binary format note
------------------
We intentionally do *not* write explicit ids into the binary files
(to keep the format simple and streaming-friendly). The C++ reader can
re-generate sequential ids on load via:
  BinaryReadOptions.generate_ids_if_missing = true.
"""

from __future__ import annotations

import argparse
import json
import os
import struct
import sys
import time
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple


def _eprint(*args: Any) -> None:
    print(*args, file=sys.stderr)


# ----------------------------
# SJSBOX v1 binary writer
# ----------------------------


@dataclass(frozen=True)
class SJSBoxHeader:
    magic: bytes = b"SJSBOX\x00\x00"
    version: int = 1
    endian_marker: int = 0x0102030405060708


FLAG_HALF_OPEN = 1 << 0
FLAG_HAS_IDS = 1 << 1


def _ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def write_sjsbox_relation(
    path: str,
    lower,
    upper,
    *,
    name: str = "",
    scalar_bits: int = 32,
    half_open: bool = True,
    write_ids: bool = False,
    chunk_rows: int = 1_000_000,
) -> None:
    """Write one relation in SJSBOX v1 format.

    Parameters
    ----------
    lower, upper
        NumPy arrays shaped (n, d). They will be converted to float32/float64.
    write_ids
        If True, writes an explicit u32 id per record.
        NOTE: Writing ids requires per-record interleaving (lo/hi/id), which
        makes streaming slower. This repo generally doesn't need explicit ids
        in files.
    """

    try:
        import numpy as np
    except Exception as ex:
        raise RuntimeError("numpy is required to write SJSBOX files") from ex

    lower = np.asarray(lower)
    upper = np.asarray(upper)
    if lower.ndim != 2 or upper.ndim != 2:
        raise ValueError("lower/upper must be 2D arrays")
    if lower.shape != upper.shape:
        raise ValueError(f"lower/upper shape mismatch: {lower.shape} vs {upper.shape}")

    n, d = lower.shape
    if d <= 0:
        raise ValueError("dimension d must be > 0")

    if scalar_bits == 32:
        dt = np.dtype("<f4")
    elif scalar_bits == 64:
        dt = np.dtype("<f8")
    else:
        raise ValueError("scalar_bits must be 32 or 64")

    # Ensure contiguous little-endian floats.
    lower = np.ascontiguousarray(lower, dtype=dt)
    upper = np.ascontiguousarray(upper, dtype=dt)

    flags = 0
    if half_open:
        flags |= FLAG_HALF_OPEN
    if write_ids:
        flags |= FLAG_HAS_IDS

    header = SJSBoxHeader()
    hdr_bytes = struct.pack(
        "<8sIIIIQQQQQQ",
        header.magic,
        header.version,
        int(d),
        int(scalar_bits),
        int(flags),
        int(n),
        header.endian_marker,
        0,
        0,
        0,
        0,
    )
    name_bytes = name.encode("utf-8")
    if len(name_bytes) > 2**32 - 1:
        raise ValueError("relation name too long")

    _ensure_parent_dir(path)
    with open(path, "wb") as f:
        f.write(hdr_bytes)
        f.write(struct.pack("<I", len(name_bytes)))
        if name_bytes:
            f.write(name_bytes)

        # Record layout per row: [lo_0..lo_{d-1}, hi_0..hi_{d-1}] (+ optional id)
        # For performance and memory, write in chunks.
        if chunk_rows <= 0:
            chunk_rows = 1_000_000

        if write_ids:
            # Interleaving ids with floats is necessary. We support it, but it is
            # intentionally not the default.
            id_dt = np.dtype("<u4")
            for i in range(0, n, chunk_rows):
                j = min(n, i + chunk_rows)
                m = j - i
                inter = np.empty((m, 2 * d), dtype=dt)
                inter[:, :d] = lower[i:j]
                inter[:, d:] = upper[i:j]
                ids = np.arange(i, j, dtype=id_dt)

                # Write record-by-record (still chunked to keep peak memory bounded).
                # This is slower than the no-id path, but correct.
                for r in range(m):
                    f.write(inter[r].tobytes(order="C"))
                    f.write(ids[r].tobytes(order="C"))
        else:
            # Fast path: write only floats (no ids), fully vectorized.
            for i in range(0, n, chunk_rows):
                j = min(n, i + chunk_rows)
                m = j - i
                inter = np.empty((m, 2 * d), dtype=dt)
                inter[:, :d] = lower[i:j]
                inter[:, d:] = upper[i:j]
                f.write(inter.tobytes(order="C"))


def write_relation_csv(
    path: str,
    lower,
    upper,
    *,
    sep: str = ",",
    include_header: bool = True,
    chunk_rows: int = 200_000,
) -> None:
    """Write a relation as a simple CSV/TSV file.

    Format matches sjs/io/csv_io.h expectations:
      id, lo0..lo(d-1), hi0..hi(d-1)
    """

    try:
        import numpy as np
    except Exception as ex:
        raise RuntimeError("numpy is required to write CSV") from ex

    if sep == "tab" or sep == "\\t":
        sep = "\t"
    if len(sep) != 1:
        raise ValueError("sep must be a single character (or 'tab')")
    sep_ch = sep

    lower = np.asarray(lower)
    upper = np.asarray(upper)
    if lower.shape != upper.shape or lower.ndim != 2:
        raise ValueError("lower/upper must be 2D arrays with the same shape")

    n, d = lower.shape

    _ensure_parent_dir(path)
    with open(path, "w", encoding="utf-8") as f:
        if include_header:
            cols = ["id"] + [f"lo{i}" for i in range(d)] + [f"hi{i}" for i in range(d)]
            f.write(sep_ch.join(cols) + "\n")

        if chunk_rows <= 0:
            chunk_rows = 200_000

        # Stream rows in chunks to keep memory bounded.
        for i in range(0, n, chunk_rows):
            j = min(n, i + chunk_rows)
            ids = np.arange(i, j, dtype=np.int64)
            # Build chunk as strings.
            # This is for debugging only; do not use for giant datasets.
            for r, idx in enumerate(ids):
                row = [str(int(idx))]
                row += [repr(float(x)) for x in lower[i + r]]
                row += [repr(float(x)) for x in upper[i + r]]
                f.write(sep_ch.join(row) + "\n")


# ----------------------------
# RectGen driver
# ----------------------------


def _safe_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate synthetic box datasets using alacarte-rectgen and export to SJSBOX binary.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Scale/shape controls.
    p.add_argument("--nR", type=int, required=True, help="|R|")
    p.add_argument("--nS", type=int, required=True, help="|S|")
    p.add_argument("--d", type=int, required=True, help="dimension d")
    p.add_argument("--alpha_out", type=float, required=True, help="target alpha_out = |J|/(|R|+|S|)")
    p.add_argument("--seed", type=int, default=1, help="generator seed")

    # Output.
    p.add_argument("--out_r", type=str, required=True, help="output path for R (SJSBOX binary)")
    p.add_argument("--out_s", type=str, required=True, help="output path for S (SJSBOX binary)")
    p.add_argument("--dataset_name", type=str, default="synthetic", help="semantic dataset name for reports")
    p.add_argument("--report_path", type=str, default="", help="write JSON audit report to this path")
    p.add_argument("--print_report", action="store_true", help="print JSON audit report to stdout")

    # Optional CSV export (debug only).
    p.add_argument("--write_csv", action="store_true", help="also write CSV relations (debug; can be huge)")
    p.add_argument("--csv_r", type=str, default="", help="override CSV R path")
    p.add_argument("--csv_s", type=str, default="", help="override CSV S path")
    p.add_argument("--csv_sep", type=str, default=",", help="CSV separator: ',' or 'tab' or '\\t'")

    # Audit.
    p.add_argument(
        "--audit_pairs",
        type=int,
        default=2_000_000,
        help="num_pairs for pair-sampling audit (0 disables audit)",
    )
    p.add_argument("--audit_seed", type=int, default=0, help="seed for pair-sampling audit")

    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)

    if args.nR <= 0 or args.nS <= 0:
        _eprint("[rectgen][FATAL] nR and nS must be > 0")
        return 2
    if args.d <= 0:
        _eprint("[rectgen][FATAL] d must be > 0")
        return 2
    if not (args.alpha_out >= 0.0):
        _eprint("[rectgen][FATAL] alpha_out must be >= 0")
        return 2

    try:
        import numpy as np
    except Exception:
        _eprint("[rectgen][FATAL] numpy is required")
        return 2

    try:
        import alacarte_rectgen as ar
    except Exception as ex:
        _eprint("[rectgen][FATAL] Failed to import alacarte_rectgen.")
        _eprint("Install with: pip install alacarte-rectgen")
        _eprint(f"Import error: {ex}")
        return 2

    # ----------------------------
    # Fixed configuration (as specified)
    # ----------------------------
    universe = None
    volume_dist = "normal"
    volume_cv = 0.25
    shape_sigma = 0.5
    tune_tol_rel = 0.01
    dtype = np.float32

    t0 = time.time()
    R, S, info = ar.make_rectangles_R_S(
        nR=int(args.nR),
        nS=int(args.nS),
        alpha_out=float(args.alpha_out),
        d=int(args.d),
        universe=universe,
        volume_dist=volume_dist,
        volume_cv=float(volume_cv),
        shape_sigma=float(shape_sigma),
        seed=int(args.seed),
        tune_tol_rel=float(tune_tol_rel),
        dtype=dtype,
    )
    gen_sec = time.time() - t0

    # Audit realized alpha via pair sampling (optional).
    alpha_hat = None
    p_hat = None
    audit_sec = 0.0
    if int(args.audit_pairs) > 0:
        t1 = time.time()
        alpha_hat, p_hat = ar.estimate_alpha_by_pair_sampling(
            R,
            S,
            num_pairs=int(args.audit_pairs),
            seed=int(args.audit_seed),
        )
        audit_sec = time.time() - t1

    # Export to SJSBOX binary.
    t2 = time.time()
    write_sjsbox_relation(
        args.out_r,
        R.lower,
        R.upper,
        name="R",
        scalar_bits=32,
        half_open=True,
        write_ids=False,
    )
    write_sjsbox_relation(
        args.out_s,
        S.lower,
        S.upper,
        name="S",
        scalar_bits=32,
        half_open=True,
        write_ids=False,
    )
    io_sec = time.time() - t2

    # Optional CSV (debug only).
    csv_sec = 0.0
    csv_r_path = ""
    csv_s_path = ""
    if args.write_csv:
        csv_r_path = args.csv_r or (os.path.splitext(args.out_r)[0] + ".csv")
        csv_s_path = args.csv_s or (os.path.splitext(args.out_s)[0] + ".csv")
        t3 = time.time()
        write_relation_csv(csv_r_path, R.lower, R.upper, sep=args.csv_sep, include_header=True)
        write_relation_csv(csv_s_path, S.lower, S.upper, sep=args.csv_sep, include_header=True)
        csv_sec = time.time() - t3

    # Prepare report.
    alpha_target = float(args.alpha_out)
    eps_alpha = None
    if alpha_hat is not None and alpha_target != 0.0:
        eps_alpha = abs(float(alpha_hat) - alpha_target) / abs(alpha_target)

    report: Dict[str, Any] = {
        "generator": "alacarte_rectgen",
        "dataset": str(args.dataset_name),
        "dim": int(args.d),
        "n_r": int(args.nR),
        "n_s": int(args.nS),
        "alpha_target": alpha_target,
        "alpha_expected_est": _safe_float(info.get("alpha_expected_est")),
        "coverage": _safe_float(info.get("coverage")),
        "pair_intersection_prob_est": _safe_float(info.get("pair_intersection_prob_est")),
        "audit_num_pairs": int(args.audit_pairs),
        "audit_seed": int(args.audit_seed),
        "alpha_hat_est": _safe_float(alpha_hat),
        "p_hat_est": _safe_float(p_hat),
        "epsilon_alpha": _safe_float(eps_alpha),
        "fixed_params": {
            "universe": None,
            "volume_dist": volume_dist,
            "volume_cv": float(volume_cv),
            "shape_sigma": float(shape_sigma),
            "tune_tol_rel": float(tune_tol_rel),
            "dtype": "float32",
        },
        "paths": {
            "out_r": os.path.abspath(args.out_r),
            "out_s": os.path.abspath(args.out_s),
            "csv_r": os.path.abspath(csv_r_path) if csv_r_path else "",
            "csv_s": os.path.abspath(csv_s_path) if csv_s_path else "",
        },
        "timing_sec": {
            "generation": float(gen_sec),
            "audit": float(audit_sec),
            "binary_io": float(io_sec),
            "csv_io": float(csv_sec),
        },
    }

    if args.report_path:
        _ensure_parent_dir(args.report_path)
        with open(args.report_path, "w", encoding="utf-8") as f:
            json.dump(report, f, ensure_ascii=False, indent=2)
            f.write("\n")

    if args.print_report:
        print(json.dumps(report, ensure_ascii=False))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
