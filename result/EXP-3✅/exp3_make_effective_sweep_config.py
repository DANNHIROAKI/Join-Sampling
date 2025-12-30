#!/usr/bin/env python3
"""Create an "effective" EXP-3 sweep config JSON.

EXP-3 (RQ3) sweeps the density parameter alpha = |J|/(n_r+n_s) using the
synthetic stripe generator.

This helper:
  1) Loads an input sweep JSON (typically config/sweeps/sweep_alpha.json)
  2) Pins base.output.out_dir to a specified directory
  3) Applies optional overrides from environment variables so that the caller
     can tweak parameters without editing JSON.

Recognized environment variables (all optional):
  EXP3_NR, EXP3_NS       : base.dataset.synthetic.n_r / n_s
  EXP3_T                 : base.run.t
  EXP3_REPEATS           : base.run.repeats
  EXP3_JSTAR             : base.run.j_star
  EXP3_ENUM_CAP          : base.run.enum_cap
  EXP3_THREADS           : base.sys.threads
  EXP3_ALPHA_LIST        : comma-separated list -> sweep.alpha (floats)
  EXP3_METHODS           : comma-separated list -> sweep.method
  EXP3_VARIANTS          : comma-separated list -> sweep.variant

Example:
  python3 run/include/exp3_make_effective_sweep_config.py \
    --in config/sweeps/sweep_alpha.json \
    --out run/temp/exp3/sweep_exp3_effective.json \
    --out_dir run/temp/exp3
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Any, Callable, Dict, List, Optional


def get_env(name: str) -> Optional[str]:
    v = os.environ.get(name)
    if v is None:
        return None
    v = v.strip()
    return v if v else None


def parse_csv_list(s: str, cast: Callable[[str], Any]) -> List[Any]:
    out: List[Any] = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        out.append(cast(part))
    return out


def ensure_dict(spec: Dict[str, Any], key: str) -> Dict[str, Any]:
    val = spec.get(key)
    if not isinstance(val, dict):
        val = {}
        spec[key] = val
    return val


def main() -> int:
    ap = argparse.ArgumentParser(description="Build an effective sweep config for EXP-3")
    ap.add_argument("--in", dest="in_path", required=True, help="Input sweep JSON")
    ap.add_argument("--out", dest="out_path", required=True, help="Output effective JSON")
    ap.add_argument(
        "--out_dir",
        required=True,
        help="Pinned output directory (writes base.output.out_dir = this)",
    )
    args = ap.parse_args()

    with open(args.in_path, "r", encoding="utf-8") as f:
        spec = json.load(f)

    if not isinstance(spec, dict):
        raise SystemExit(f"[EXP3][ERROR] Input JSON must be an object: {args.in_path}")

    base = ensure_dict(spec, "base")
    output = ensure_dict(base, "output")
    output["out_dir"] = args.out_dir

    # Optional overrides (no JSON editing needed)
    nr = get_env("EXP3_NR")
    ns = get_env("EXP3_NS")
    t = get_env("EXP3_T")
    rep = get_env("EXP3_REPEATS")
    j_star = get_env("EXP3_JSTAR")
    enum_cap = get_env("EXP3_ENUM_CAP")
    threads = get_env("EXP3_THREADS")

    alpha_list = get_env("EXP3_ALPHA_LIST")
    methods = get_env("EXP3_METHODS")
    variants = get_env("EXP3_VARIANTS")

    # Dataset sizes (synthetic)
    if nr is not None or ns is not None:
        dataset = ensure_dict(base, "dataset")
        synthetic = ensure_dict(dataset, "synthetic")
        if nr is not None:
            synthetic["n_r"] = int(nr)
        if ns is not None:
            synthetic["n_s"] = int(ns)

    # Run parameters
    if any(x is not None for x in (t, rep, j_star, enum_cap)):
        run = ensure_dict(base, "run")
        if t is not None:
            run["t"] = int(t)
        if rep is not None:
            run["repeats"] = int(rep)
        if j_star is not None:
            run["j_star"] = int(j_star)
        if enum_cap is not None:
            run["enum_cap"] = int(enum_cap)

    # System threads (fairness)
    if threads is not None:
        syscfg = ensure_dict(base, "sys")
        syscfg["threads"] = int(threads)

    # Sweep lists
    sweep = ensure_dict(spec, "sweep")
    if alpha_list is not None:
        sweep["alpha"] = parse_csv_list(alpha_list, float)
    if methods is not None:
        sweep["method"] = parse_csv_list(methods, str)
    if variants is not None:
        sweep["variant"] = parse_csv_list(variants, str)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_path)), exist_ok=True)
    with open(args.out_path, "w", encoding="utf-8") as f:
        json.dump(spec, f, indent=2, ensure_ascii=False)

    print("[EXP3] Wrote effective sweep config:", args.out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
