#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# EXP-3: Runtime vs alpha (RQ3)
#
# Goal:
#   Sweep density alpha = |J|/(n_r+n_s) from 0 -> 300 using the synthetic
#   stripe_ctrl_alpha generator, and measure end-to-end runtime for each method
#   (sampling + adaptive). Also derive the adaptive branch ratio (enumerate vs
#   fallback_sampling) from sweep_raw.csv.
#
# Outputs (by default):
#   results/sweeps/exp3_alpha_<timestamp>/
#     sweep_raw.csv
#     sweep_summary.csv
#     derived/adaptive_branch_ratio.csv
#     plots/*.png  (if matplotlib is available)
#
# This script is self-contained: it builds the project (Release), runs the sweep,
# and produces derived tables/plots in one shot.
#
# Usage:
#   bash run/run_exp3.sh
#
# Optional overrides (environment variables):
#   EXP3_CONFIG            : path to a sweep JSON (default: config/sweeps/sweep_alpha.json)
#   EXP3_OUT_DIR           : output directory (default: results/sweeps/exp3_alpha_<ts>)
#   EXP3_BUILD_DIR         : build directory (default: build_exp3)
#   EXP3_JOBS              : build parallelism (default: nproc)
#
#   # Override key experiment params WITHOUT editing JSON:
#   EXP3_NR, EXP3_NS       : set n_r / n_s
#   EXP3_T                 : set sample size t
#   EXP3_REPEATS           : set repeats
#   EXP3_JSTAR             : set adaptive threshold j_star
#   EXP3_ENUM_CAP          : set enum_cap (0 = no cap)
#   EXP3_THREADS           : set sys.threads (recommended 1 for fairness)
#   EXP3_ALPHA_LIST        : comma-separated alpha list, e.g. "0,1,2,3,5,10,30,100,300"
#   EXP3_METHODS           : comma-separated method list
#   EXP3_VARIANTS          : comma-separated variant list (default: from config)
#
#   EXP3_RUN_TESTS=1       : run ctest after build (default: 0)
#   EXP3_PLOT=0            : skip plotting even if matplotlib exists (default: 1)
# ------------------------------------------------------------------------------

set -euo pipefail

# -----------------------------
# Resolve repo root
# -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# -----------------------------
# Dependency checks
# -----------------------------
for cmd in cmake python3; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "[EXP3][ERROR] Required command not found: ${cmd}" >&2
    exit 1
  fi
done

# -----------------------------
# Config & output locations
# -----------------------------
DEFAULT_CONFIG="${ROOT_DIR}/config/sweeps/sweep_alpha.json"
CONFIG_PATH="${EXP3_CONFIG:-${DEFAULT_CONFIG}}"
if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "[EXP3][ERROR] Sweep config not found: ${CONFIG_PATH}" >&2
  exit 2
fi

TS="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${EXP3_OUT_DIR:-${ROOT_DIR}/results/sweeps/exp3_alpha_${TS}}"
BUILD_DIR="${EXP3_BUILD_DIR:-${ROOT_DIR}/build_exp3}"
JOBS="${EXP3_JOBS:-$(command -v nproc >/dev/null 2>&1 && nproc || echo 4)}"
RUN_TESTS="${EXP3_RUN_TESTS:-0}"
DO_PLOT="${EXP3_PLOT:-1}"

mkdir -p "${OUT_DIR}"

# -----------------------------
# Repro / fairness knobs
# -----------------------------
# Keep everything single-threaded unless you explicitly override EXP3_THREADS.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# -----------------------------
# Build (Release)
# -----------------------------
echo "[EXP3] Building project (Release) ..."
rm -rf "${BUILD_DIR}"
cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" -j "${JOBS}"

if [[ "${RUN_TESTS}" == "1" ]]; then
  echo "[EXP3] Running tests (ctest) ..."
  (cd "${BUILD_DIR}" && ctest --output-on-failure)
fi

# -----------------------------
# Locate binaries
# -----------------------------
SJS_SWEEP=""
if [[ -x "${BUILD_DIR}/sjs_sweep" ]]; then
  SJS_SWEEP="${BUILD_DIR}/sjs_sweep"
elif [[ -x "${ROOT_DIR}/sjs_sweep" ]]; then
  SJS_SWEEP="${ROOT_DIR}/sjs_sweep"
else
  SJS_SWEEP="$(find "${BUILD_DIR}" "${ROOT_DIR}" -maxdepth 3 -type f -name sjs_sweep -perm -111 2>/dev/null | head -n 1 || true)"
fi

if [[ -z "${SJS_SWEEP}" ]]; then
  echo "[EXP3][ERROR] Cannot find executable: sjs_sweep" >&2
  exit 3
fi

echo "[EXP3] Using sjs_sweep: ${SJS_SWEEP}"

# -----------------------------
# Prepare an "effective" sweep JSON (copy + optional overrides)
# -----------------------------
EFFECTIVE_CONFIG="${OUT_DIR}/sweep_exp3_effective.json"

# Export vars needed by the embedded Python snippets
export CONFIG_PATH
export EFFECTIVE_CONFIG
export OUT_DIR

python3 - <<'PY'
import json, os, sys

in_path  = os.environ["CONFIG_PATH"]
out_path = os.environ["EFFECTIVE_CONFIG"]
out_dir  = os.environ["OUT_DIR"]

def get_env(name):
    v = os.environ.get(name)
    return v if v is not None and v != "" else None

def parse_csv_list(s, cast=str):
    # "a,b,c" -> [cast(a), cast(b), cast(c)]
    return [cast(x.strip()) for x in s.split(",") if x.strip() != ""]

with open(in_path, "r", encoding="utf-8") as f:
    spec = json.load(f)

# Always pin output dir to this run.
spec.setdefault("base", {}).setdefault("output", {})["out_dir"] = out_dir

# --- Optional overrides (no JSON editing needed) ---
nr = get_env("EXP3_NR")
ns = get_env("EXP3_NS")
t  = get_env("EXP3_T")
rep = get_env("EXP3_REPEATS")
j_star = get_env("EXP3_JSTAR")
enum_cap = get_env("EXP3_ENUM_CAP")
threads = get_env("EXP3_THREADS")

alpha_list = get_env("EXP3_ALPHA_LIST")
methods    = get_env("EXP3_METHODS")
variants   = get_env("EXP3_VARIANTS")

# Dataset sizes (synthetic)
if nr is not None:
    spec["base"]["dataset"]["synthetic"]["n_r"] = int(nr)
if ns is not None:
    spec["base"]["dataset"]["synthetic"]["n_s"] = int(ns)

# Run parameters
if t is not None:
    spec["base"]["run"]["t"] = int(t)
if rep is not None:
    spec["base"]["run"]["repeats"] = int(rep)
if j_star is not None:
    spec["base"]["run"]["j_star"] = int(j_star)
if enum_cap is not None:
    spec["base"]["run"]["enum_cap"] = int(enum_cap)

# System threads (fairness)
if threads is not None:
    spec["base"].setdefault("sys", {})["threads"] = int(threads)

# Sweep lists
if alpha_list is not None:
    # allow floats; store as numbers in JSON
    spec.setdefault("sweep", {})["alpha"] = parse_csv_list(alpha_list, float)
if methods is not None:
    spec.setdefault("sweep", {})["method"] = parse_csv_list(methods, str)
if variants is not None:
    spec.setdefault("sweep", {})["variant"] = parse_csv_list(variants, str)

with open(out_path, "w", encoding="utf-8") as f:
    json.dump(spec, f, indent=2, ensure_ascii=False)

print("[EXP3] Wrote effective sweep config:", out_path)
PY


# -----------------------------
# Run the sweep
# -----------------------------
echo "[EXP3] Running sweep (alpha: 0 -> 300) ..."
LOG_PATH="${OUT_DIR}/exp3_sweep.log"
"${SJS_SWEEP}" --config="${EFFECTIVE_CONFIG}" 2>&1 | tee "${LOG_PATH}"

RAW_CSV="${OUT_DIR}/sweep_raw.csv"
SUMMARY_CSV="${OUT_DIR}/sweep_summary.csv"

if [[ ! -f "${RAW_CSV}" || ! -f "${SUMMARY_CSV}" ]]; then
  echo "[EXP3][ERROR] Expected outputs not found under: ${OUT_DIR}" >&2
  echo "  Missing: ${RAW_CSV} or ${SUMMARY_CSV}" >&2
  exit 4
fi

echo "[EXP3] Sweep done."
echo "[EXP3] Raw     : ${RAW_CSV}"
echo "[EXP3] Summary : ${SUMMARY_CSV}"

# -----------------------------
# Derive: adaptive branch ratio
# -----------------------------
DERIVED_DIR="${OUT_DIR}/derived"
PLOTS_DIR="${OUT_DIR}/plots"
mkdir -p "${DERIVED_DIR}" "${PLOTS_DIR}"

BRANCH_CSV="${DERIVED_DIR}/adaptive_branch_ratio.csv"
export RAW_CSV BRANCH_CSV

python3 - <<'PY'
import csv, os, math
from collections import defaultdict

raw_path = os.environ["RAW_CSV"]
out_path = os.environ["BRANCH_CSV"]

# Group adaptive runs by (dataset,generator,n_r,n_s,method,alpha)
grp = defaultdict(lambda: {"total":0, "enumerate_all":0, "fallback_sampling":0, "fallback_sampling_no_pilot":0, "other":0})

with open(raw_path, "r", newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    required = {"variant","method","alpha","adaptive_branch","ok","dataset","generator","n_r","n_s"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"[EXP3][ERROR] sweep_raw.csv missing columns: {sorted(missing)}")

    for row in reader:
        if row["variant"] != "adaptive":
            continue
        key = (
            row["dataset"],
            row["generator"],
            row["n_r"],
            row["n_s"],
            row["method"],
            row["alpha"],
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

# Write derived CSV
fieldnames = [
    "dataset","generator","n_r","n_s","method","alpha",
    "runs",
    "enumerate_all","fallback_sampling","fallback_sampling_no_pilot","other",
    "enumerate_frac","fallback_frac","fallback_no_pilot_frac","other_frac",
]
rows = []
for (dataset,generator,n_r,n_s,method,alpha), g in grp.items():
    total = g["total"] or 1
    rows.append({
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
        "enumerate_frac": g["enumerate_all"]/total,
        "fallback_frac": g["fallback_sampling"]/total,
        "fallback_no_pilot_frac": g["fallback_sampling_no_pilot"]/total,
        "other_frac": g["other"]/total,
    })

# Sort by method then numeric alpha
def alpha_key(a):
    try:
        return float(a)
    except:
        return math.inf

rows.sort(key=lambda r: (r["method"], alpha_key(r["alpha"])))

with open(out_path, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

print("[EXP3] Wrote adaptive branch ratio:", out_path)
PY

# -----------------------------
# Optional: plotting (matplotlib)
# -----------------------------
if [[ "${DO_PLOT}" == "1" ]]; then
  export SUMMARY_CSV PLOTS_DIR BRANCH_CSV
  python3 - <<'PY'
import os, csv, math
from collections import defaultdict

summary_path = os.environ["SUMMARY_CSV"]
branch_path  = os.environ["BRANCH_CSV"]
plots_dir    = os.environ["PLOTS_DIR"]

# Try to import matplotlib; if unavailable, just exit gracefully.
try:
    import matplotlib.pyplot as plt
except Exception as e:
    print("[EXP3] matplotlib not available; skip plotting. (", e, ")")
    raise SystemExit(0)

# ---------- Plot 1: runtime vs alpha (per method, sampling vs adaptive) ----------
# Read summary into groups: (method,variant) -> [(alpha, wall_median_ms, wall_p95_ms)]
series = defaultdict(list)
with open(summary_path, "r", newline="", encoding="utf-8") as f:
    rd = csv.DictReader(f)
    required = {"alpha","method","variant","wall_median_ms","wall_p95_ms","ok_rate"}
    missing = required - set(rd.fieldnames or [])
    if missing:
        raise SystemExit(f"[EXP3][ERROR] sweep_summary.csv missing columns: {sorted(missing)}")
    for row in rd:
        # Keep even if ok_rate<1; plot what we have
        try:
            a = float(row["alpha"])
            med = float(row["wall_median_ms"])
            p95 = float(row["wall_p95_ms"])
        except:
            continue
        series[(row["method"], row["variant"])].append((a, med, p95, float(row.get("ok_rate","1") or "1")))

# Sort points
for k in list(series.keys()):
    series[k].sort(key=lambda x: x[0])

methods = sorted({m for (m, v) in series.keys()})
for m in methods:
    plt.figure()
    for v in ["sampling","adaptive","enum_sampling"]:
        pts = series.get((m,v))
        if not pts:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        plt.plot(xs, ys, marker="o", label=v)
    plt.xscale("symlog", linthresh=0.03)
    plt.yscale("log")
    plt.xlabel("alpha = |J|/(n_r+n_s)")
    plt.ylabel("wall_median_ms (log)")
    plt.title(f"EXP-3 Runtime vs alpha — {m}")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    out = os.path.join(plots_dir, f"runtime_vs_alpha_{m}.png")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    print("[EXP3] Wrote plot:", out)

# Also make an "all methods" plot (may be crowded).
plt.figure(figsize=(12,7))
for (m,v), pts in series.items():
    if v != "sampling":
        continue
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    plt.plot(xs, ys, marker="o", label=m)
plt.xscale("symlog", linthresh=0.03)
plt.yscale("log")
plt.xlabel("alpha = |J|/(n_r+n_s)")
plt.ylabel("wall_median_ms (log)")
plt.title("EXP-3 Runtime vs alpha — sampling (all methods)")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend(ncol=2, fontsize=8)
out = os.path.join(plots_dir, "runtime_vs_alpha_all_methods_sampling.png")
plt.tight_layout()
plt.savefig(out, dpi=200)
plt.close()
print("[EXP3] Wrote plot:", out)

# ---------- Plot 2: adaptive branch ratio vs alpha (enumerate_frac) ----------
branch = defaultdict(list)  # method -> [(alpha, enumerate_frac, fallback_frac)]
with open(branch_path, "r", newline="", encoding="utf-8") as f:
    rd = csv.DictReader(f)
    required = {"method","alpha","enumerate_frac","fallback_frac","fallback_no_pilot_frac"}
    missing = required - set(rd.fieldnames or [])
    if missing:
        raise SystemExit(f"[EXP3][ERROR] adaptive_branch_ratio.csv missing columns: {sorted(missing)}")
    for row in rd:
        try:
            a = float(row["alpha"])
            enumf = float(row["enumerate_frac"])
            fallf = float(row["fallback_frac"]) + float(row["fallback_no_pilot_frac"])
        except:
            continue
        branch[row["method"]].append((a, enumf, fallf))

plt.figure(figsize=(12,6))
for m, pts in sorted(branch.items()):
    pts.sort(key=lambda x: x[0])
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]  # enumerate fraction
    plt.plot(xs, ys, marker="o", label=m)
plt.xscale("symlog", linthresh=0.03)
plt.xlabel("alpha = |J|/(n_r+n_s)")
plt.ylabel("adaptive enumerate fraction")
plt.title("EXP-3 Adaptive branch ratio vs alpha (enumerate fraction)")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend(ncol=2, fontsize=8)
out = os.path.join(plots_dir, "adaptive_branch_ratio_enumerate_frac.png")
plt.tight_layout()
plt.savefig(out, dpi=200)
plt.close()
print("[EXP3] Wrote plot:", out)
PY
fi

echo "[EXP3] Done."
echo "----------------------------------------"
echo "EXP-3 output directory:"
echo "  ${OUT_DIR}"
echo "Key files:"
echo "  ${RAW_CSV}"
echo "  ${SUMMARY_CSV}"
echo "  ${BRANCH_CSV}"
echo "  ${LOG_PATH}"
echo "----------------------------------------"
