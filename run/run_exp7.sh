#!/usr/bin/env bash
# run/run_exp7.sh
#
# EXP-7: Phase breakdown & explanatory profiling
# ----------------------------------------------
# One-shot script that (1) builds the project, (2) runs EXP-7 sweeps on 2–3
# representative density points, and (3) post-processes sweep_raw.csv into
# a phase breakdown table + stacked-bar figures.
#
# Outputs (default):
#   results/sweeps/exp7_profile/sampling_adaptive/sweep_raw.csv
#   results/sweeps/exp7_profile/sampling_adaptive/sweep_summary.csv
#   results/sweeps/exp7_profile/enum_sparse/sweep_raw.csv            (optional)
#   results/sweeps/exp7_profile/enum_sparse/sweep_summary.csv        (optional)
#   results/sweeps/exp7_profile/exp7_breakdown_median.csv
#   results/sweeps/exp7_profile/figs/exp7_phase_breakdown_alpha_*.png
#
# Usage:
#   bash run/run_exp7.sh
#
# Optional environment overrides:
#   BUILD_DIR=build            (default: <repo>/build)
#   BUILD_TYPE=Release         (default: Release)
#   OUT_BASE=results/...       (default: <repo>/results/sweeps/exp7_profile)
#
#   NR=100000 NS=100000        (relation sizes)
#   T=100000                   (samples per run)
#   REPEATS=3                  (repeats per configuration)
#   GEN_SEED=1 SEED=1          (dataset/run seeds)
#   THREADS=1                  (fairness: single-thread by default)
#
#   ALPHAS="0.1 3 30"          (representative sparse/medium/dense points)
#   METHODS="ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection tsunami"
#
#   J_STAR=1000000             (adaptive threshold)
#   ENUM_CAP=0                 (0 => no cap; only used by enum-based runners)
#
#   INCLUDE_ENUM_SPARSE=1      (run enum_sampling at SPARSE_ALPHA)
#   SPARSE_ALPHA=0.1
#
# Notes:
# - This script uses the built-in sweep harness (sjs_sweep) and expects a small JSON subset.
# - Post-processing is done via python3 stdlib; plotting requires matplotlib (optional).
set -euo pipefail
IFS=$'\n\t'

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"
}

find_bin() {
  # Usage: find_bin <build_dir> <name>
  local bdir="$1"
  local name="$2"
  local cand
  for cand in \
    "${bdir}/${name}" \
    "${bdir}/Release/${name}" \
    "${bdir}/Debug/${name}" \
    "${bdir}/bin/${name}" \
    "${bdir}/src/${name}" \
  ; do
    if [[ -x "${cand}" ]]; then
      echo "${cand}"
      return 0
    fi
  done
  return 1
}

# ---------------------------
# Locate repo root
# ---------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ---------------------------
# User-tunable parameters
# ---------------------------
BUILD_DIR="${BUILD_DIR:-${REPO_ROOT}/build}"
BUILD_TYPE="${BUILD_TYPE:-Release}"
OUT_BASE="${OUT_BASE:-${REPO_ROOT}/results/sweeps/exp7_profile}"

NR="${NR:-100000}"
NS="${NS:-100000}"
T="${T:-100000}"
REPEATS="${REPEATS:-3}"
GEN_SEED="${GEN_SEED:-1}"
SEED="${SEED:-1}"
THREADS="${THREADS:-1}"

ALPHAS="${ALPHAS:-0.1 3 30}"
METHODS="${METHODS:-ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection tsunami}"

J_STAR="${J_STAR:-1000000}"
ENUM_CAP="${ENUM_CAP:-0}"

INCLUDE_ENUM_SPARSE="${INCLUDE_ENUM_SPARSE:-1}"
SPARSE_ALPHA="${SPARSE_ALPHA:-0.1}"

# ---------------------------
# Dependency checks
# ---------------------------
need_cmd cmake
need_cmd python3

log "Repo root   : ${REPO_ROOT}"
log "Build dir   : ${BUILD_DIR} (${BUILD_TYPE})"
log "Out base    : ${OUT_BASE}"
log "N (R,S)     : ${NR}, ${NS}"
log "t, repeats  : ${T}, ${REPEATS}"
log "alphas      : ${ALPHAS}"
log "methods     : ${METHODS}"
log "threads     : ${THREADS}"
log "enum sparse : INCLUDE_ENUM_SPARSE=${INCLUDE_ENUM_SPARSE} SPARSE_ALPHA=${SPARSE_ALPHA}"

# ---------------------------
# Build (Release)
# ---------------------------
log "== Build =="
mkdir -p "${BUILD_DIR}"
cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
cmake --build "${BUILD_DIR}" -j

SJS_SWEEP="$(find_bin "${BUILD_DIR}" sjs_sweep || true)"
[[ -n "${SJS_SWEEP}" ]] || die "Could not locate sjs_sweep binary under ${BUILD_DIR}"
log "Using sjs_sweep: ${SJS_SWEEP}"

# ---------------------------
# Generate sweep configs
# ---------------------------
GEN_DIR="${REPO_ROOT}/run/_generated_exp7"
mkdir -p "${GEN_DIR}"

OUT_SA="${OUT_BASE}/sampling_adaptive"
OUT_ENUM="${OUT_BASE}/enum_sparse"

CFG_SA="${GEN_DIR}/exp7_sampling_adaptive.json"
CFG_ENUM="${GEN_DIR}/exp7_enum_sparse.json"

log "== Generate sweep config: sampling+adaptive =="
REPO_ROOT="${REPO_ROOT}" OUT_DIR="${OUT_SA}" CFG_PATH="${CFG_SA}" \
NR="${NR}" NS="${NS}" T="${T}" REPEATS="${REPEATS}" GEN_SEED="${GEN_SEED}" SEED="${SEED}" \
THREADS="${THREADS}" J_STAR="${J_STAR}" ENUM_CAP="${ENUM_CAP}" \
ALPHAS="${ALPHAS}" METHODS="${METHODS}" VARIANTS_JSON='["sampling","adaptive"]' \
python3 - <<'PY'
import json, os
from pathlib import Path

def env_int(k: str) -> int:
    return int(os.environ[k])

def env_float_list(k: str):
    return [float(x) for x in os.environ[k].split()]

def env_str_list(k: str):
    return [x for x in os.environ[k].split() if x]

repo = Path(os.environ["REPO_ROOT"])
out_dir = os.environ["OUT_DIR"]
cfg_path = Path(os.environ["CFG_PATH"])

nr = env_int("NR")
ns = env_int("NS")
t = env_int("T")
repeats = env_int("REPEATS")
gen_seed = env_int("GEN_SEED")
seed = env_int("SEED")
threads = env_int("THREADS")
j_star = env_int("J_STAR")
enum_cap = env_int("ENUM_CAP")

alphas = env_float_list("ALPHAS")
methods = env_str_list("METHODS")
variants = json.loads(os.environ["VARIANTS_JSON"])

cfg = {
  "base": {
    "dataset": {
      "source": "synthetic",
      "name": "exp7_stripe",
      "dim": 2,
      "synthetic": {
        "generator": "stripe_ctrl_alpha",
        "n_r": nr,
        "n_s": ns,
        "alpha": alphas[0] if alphas else 0.1,
        "seed": gen_seed,
        "params": {
          "control_axis": 1,
          "core_lo": 0.45,
          "core_hi": 0.55,
          "gap_factor": 0.1,
          "delta_factor": 0.25,
          "shuffle_strips": True,
          "shuffle_r": False,
          "swap_sides": False
        }
      }
    },
    "run": {
      "method": "ours",
      "variant": "sampling",
      "t": t,
      "seed": seed,
      "repeats": repeats,
      "j_star": j_star,
      "enum_cap": enum_cap,
      "write_samples": False,
      "verify": False,
      "extra": {
        "pbsm_scheme": "stripes",
        "pbsm_k": 0,
        "pbsm_part_axis": 0,
        "pbsm_sweep": "orthogonal",

        "tlsop_nx": 128,
        "tlsop_ny": 128,
        "tlsop_reuse_sort": True,

        "sirs_outer": "",
        "sirs_leaf_size": 32,

        "rej_index": "S",
        "rej_rep_center": False,
        "rej_count_draws": 50000,
        "rej_max_factor": 10000
      }
    },
    "output": {
      "out_dir": out_dir,
      "run_tag": "exp7_sampling_adaptive"
    },
    "logging": {
      "level": "info",
      "with_timestamp": True,
      "with_thread_id": False
    },
    "sys": { "threads": threads }
  },
  "sweep": {
    "alpha": alphas,
    "method": methods,
    "variant": variants
  },
  "files": {
    "raw": "sweep_raw.csv",
    "summary": "sweep_summary.csv"
  },
  "meta": {
    "note": "EXP-7 (phase breakdown): sampling+adaptive at representative alpha points."
  }
}

cfg_path.parent.mkdir(parents=True, exist_ok=True)
cfg_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
print(f"Wrote {cfg_path}")
PY

log "== Run sweep: sampling+adaptive =="
mkdir -p "${OUT_SA}"
"${SJS_SWEEP}" --config="${CFG_SA}"

RAW_SA="${OUT_SA}/sweep_raw.csv"
[[ -f "${RAW_SA}" ]] || die "Expected raw file not found: ${RAW_SA}"
log "Wrote: ${RAW_SA}"

# Optional: enum_sampling only at sparse alpha (to avoid enumeration blow-ups at dense points)
RAW_ENUM=""
if [[ "${INCLUDE_ENUM_SPARSE}" == "1" ]]; then
  log "== Generate sweep config: enum_sampling (sparse only) =="
  REPO_ROOT="${REPO_ROOT}" OUT_DIR="${OUT_ENUM}" CFG_PATH="${CFG_ENUM}" \
  NR="${NR}" NS="${NS}" T="${T}" REPEATS="${REPEATS}" GEN_SEED="${GEN_SEED}" SEED="${SEED}" \
  THREADS="${THREADS}" J_STAR="${J_STAR}" ENUM_CAP="${ENUM_CAP}" \
  ALPHAS="${SPARSE_ALPHA}" METHODS="${METHODS}" VARIANTS_JSON='["enum_sampling"]' \
  python3 - <<'PY'
import json, os
from pathlib import Path

def env_int(k: str) -> int:
    return int(os.environ[k])

def env_float_list(k: str):
    return [float(x) for x in os.environ[k].split()]

def env_str_list(k: str):
    return [x for x in os.environ[k].split() if x]

out_dir = os.environ["OUT_DIR"]
cfg_path = Path(os.environ["CFG_PATH"])

nr = env_int("NR")
ns = env_int("NS")
t = env_int("T")
repeats = env_int("REPEATS")
gen_seed = env_int("GEN_SEED")
seed = env_int("SEED")
threads = env_int("THREADS")
j_star = env_int("J_STAR")
enum_cap = env_int("ENUM_CAP")

alphas = env_float_list("ALPHAS")  # should be a single sparse alpha
methods = env_str_list("METHODS")
variants = json.loads(os.environ["VARIANTS_JSON"])

cfg = {
  "base": {
    "dataset": {
      "source": "synthetic",
      "name": "exp7_stripe",
      "dim": 2,
      "synthetic": {
        "generator": "stripe_ctrl_alpha",
        "n_r": nr,
        "n_s": ns,
        "alpha": alphas[0] if alphas else 0.1,
        "seed": gen_seed,
        "params": {
          "control_axis": 1,
          "core_lo": 0.45,
          "core_hi": 0.55,
          "gap_factor": 0.1,
          "delta_factor": 0.25,
          "shuffle_strips": True,
          "shuffle_r": False,
          "swap_sides": False
        }
      }
    },
    "run": {
      "method": "ours",
      "variant": "enum_sampling",
      "t": t,
      "seed": seed,
      "repeats": repeats,
      "j_star": j_star,
      "enum_cap": enum_cap,
      "write_samples": False,
      "verify": False,
      "extra": {
        "pbsm_scheme": "stripes",
        "pbsm_k": 0,
        "pbsm_part_axis": 0,
        "pbsm_sweep": "orthogonal",

        "tlsop_nx": 128,
        "tlsop_ny": 128,
        "tlsop_reuse_sort": True,

        "sirs_outer": "",
        "sirs_leaf_size": 32,

        "rej_index": "S",
        "rej_rep_center": False,
        "rej_count_draws": 50000,
        "rej_max_factor": 10000
      }
    },
    "output": {
      "out_dir": out_dir,
      "run_tag": "exp7_enum_sparse"
    },
    "logging": {
      "level": "info",
      "with_timestamp": True,
      "with_thread_id": False
    },
    "sys": { "threads": threads }
  },
  "sweep": {
    "alpha": alphas,
    "method": methods,
    "variant": variants
  },
  "files": {
    "raw": "sweep_raw.csv",
    "summary": "sweep_summary.csv"
  },
  "meta": {
    "note": "EXP-7 (phase breakdown): enum_sampling only at sparse alpha to keep runtime bounded."
  }
}

cfg_path.parent.mkdir(parents=True, exist_ok=True)
cfg_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
print(f"Wrote {cfg_path}")
PY

  log "== Run sweep: enum_sampling (sparse only) =="
  mkdir -p "${OUT_ENUM}"
  "${SJS_SWEEP}" --config="${CFG_ENUM}"

  RAW_ENUM="${OUT_ENUM}/sweep_raw.csv"
  [[ -f "${RAW_ENUM}" ]] || die "Expected raw file not found: ${RAW_ENUM}"
  log "Wrote: ${RAW_ENUM}"
else
  log "Skipping enum_sampling sweep (INCLUDE_ENUM_SPARSE=${INCLUDE_ENUM_SPARSE})"
fi

# ---------------------------
# Post-process: phase breakdown
# ---------------------------
log "== Post-process: phase breakdown table + plots =="

OUT_BREAKDOWN="${OUT_BASE}/exp7_breakdown_median.csv"
OUT_MD="${OUT_BASE}/exp7_breakdown_median.md"
FIG_DIR="${OUT_BASE}/figs"
MERGED_RAW="${OUT_BASE}/exp7_merged_sweep_raw.csv"

RAW_SA="${RAW_SA}" RAW_ENUM="${RAW_ENUM}" OUT_BREAKDOWN="${OUT_BREAKDOWN}" OUT_MD="${OUT_MD}" \
FIG_DIR="${FIG_DIR}" MERGED_RAW="${MERGED_RAW}" \
python3 - <<'PY'
import csv, json, os, math, statistics
from pathlib import Path

raw_files = []
raw_sa = os.environ.get("RAW_SA", "")
raw_enum = os.environ.get("RAW_ENUM", "")

if raw_sa and Path(raw_sa).exists():
    raw_files.append(raw_sa)
if raw_enum and Path(raw_enum).exists():
    raw_files.append(raw_enum)

if not raw_files:
    raise SystemExit("No sweep_raw.csv files found to post-process.")

out_breakdown = Path(os.environ["OUT_BREAKDOWN"])
out_md = Path(os.environ["OUT_MD"])
fig_dir = Path(os.environ["FIG_DIR"])
merged_raw = Path(os.environ["MERGED_RAW"])

out_breakdown.parent.mkdir(parents=True, exist_ok=True)
fig_dir.mkdir(parents=True, exist_ok=True)

def parse_phases(s: str):
    if not s:
        return {}
    try:
        return json.loads(s)
    except Exception:
        return {}

def fnum(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return float(default)

def stage_breakdown(variant: str, ph: dict):
    g = lambda k: fnum(ph.get(k, 0.0), 0.0)

    build = g("run_build_ms")

    if variant == "sampling":
        count = g("run_count_ms")
        enum = 0.0
        sample = g("run_sample_ms")

    elif variant == "enum_sampling":
        # enum_sampling runner phases:
        # run_enum_prepare, run_enum_pass1_count, run_draw_ranks, run_enum_pass2_select
        count = 0.0
        enum = g("run_enum_prepare_ms") + g("run_enum_pass1_count_ms")
        sample = g("run_draw_ranks_ms") + g("run_enum_pass2_select_ms")

    elif variant == "adaptive":
        # adaptive runner phases:
        # run_pilot_enum_prepare, run_pilot_enum_scan,
        # run_small_join_sample_from_list (if small-join),
        # run_fallback_count, run_fallback_sample (if fallback)
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
    for k, v in ph.items():
        if isinstance(k, str) and k.startswith("run_") and k.endswith("_ms"):
            total += fnum(v, 0.0)
    if total <= 0:
        total = build + count + enum + sample

    return build, count, enum, sample, total

# Read & merge rows
rows = []
for path in raw_files:
    with open(path, "r", newline="", encoding="utf-8") as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            r["_raw_file"] = path
            rows.append(r)

# Write merged raw for convenience
if rows:
    merged_raw.parent.mkdir(parents=True, exist_ok=True)
    with open(merged_raw, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

ok_rows = [r for r in rows if str(r.get("ok", "0")).strip() == "1"]

# Compute stage times
for r in ok_rows:
    ph = parse_phases(r.get("phases_json", ""))
    variant = str(r.get("variant", ""))
    b, c, e, s, tot = stage_breakdown(variant, ph)
    r["build_ms"] = b
    r["count_ms"] = c
    r["enum_ms"] = e
    r["sample_ms"] = s
    r["run_total_ms_from_phases"] = tot
    # wall_ms exists in sweep_raw; keep it as reference
    r["wall_ms"] = fnum(r.get("wall_ms", 0.0), 0.0)

# Group by config key (median over repeats)
gcols = ["dataset", "generator", "alpha", "n_r", "n_s", "method", "variant", "t"]
groups = {}
for r in ok_rows:
    key = tuple(r.get(k, "") for k in gcols)
    groups.setdefault(key, []).append(r)

def med(vals):
    vals = [float(v) for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    return statistics.median(vals) if vals else 0.0

records = []
for key, rs in groups.items():
    rec = {k: v for k, v in zip(gcols, key)}
    rec["repeats_ok"] = len(rs)
    for m in ["build_ms", "count_ms", "enum_ms", "sample_ms", "run_total_ms_from_phases", "wall_ms"]:
        rec[m] = med([r[m] for r in rs])
    total = rec["run_total_ms_from_phases"] if rec["run_total_ms_from_phases"] > 0 else (rec["build_ms"] + rec["count_ms"] + rec["enum_ms"] + rec["sample_ms"])
    rec["build_pct"] = (rec["build_ms"] / total * 100.0) if total > 0 else 0.0
    rec["count_pct"] = (rec["count_ms"] / total * 100.0) if total > 0 else 0.0
    rec["enum_pct"] = (rec["enum_ms"] / total * 100.0) if total > 0 else 0.0
    rec["sample_pct"] = (rec["sample_ms"] / total * 100.0) if total > 0 else 0.0
    records.append(rec)

# Sort for readability
def ffloat(x):
    try:
        return float(x)
    except Exception:
        return float("nan")

records.sort(key=lambda r: (ffloat(r.get("alpha", "nan")), r.get("variant",""), r.get("method","")))

# Write breakdown CSV
out_breakdown.parent.mkdir(parents=True, exist_ok=True)
with open(out_breakdown, "w", newline="", encoding="utf-8") as f:
    fieldnames = gcols + [
        "repeats_ok",
        "build_ms", "count_ms", "enum_ms", "sample_ms",
        "run_total_ms_from_phases", "wall_ms",
        "build_pct", "count_pct", "enum_pct", "sample_pct"
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(records)

# Write a compact markdown table
lines = []
lines.append("| alpha | method | variant | Build% | Count% | Enum% | Sample% | total_ms(phases) | wall_ms |")
lines.append("|---:|---|---|---:|---:|---:|---:|---:|---:|")
for r in records:
    lines.append(
        f"| {r['alpha']} | {r['method']} | {r['variant']} | "
        f"{r['build_pct']:.1f} | {r['count_pct']:.1f} | {r['enum_pct']:.1f} | {r['sample_pct']:.1f} | "
        f"{r['run_total_ms_from_phases']:.1f} | {r['wall_ms']:.1f} |"
    )
out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

print(f"Wrote breakdown CSV: {out_breakdown}")
print(f"Wrote breakdown MD : {out_md}")
print(f"Wrote merged raw   : {merged_raw}")

# Plot stacked bars (optional; matplotlib not required for the table)
try:
    import matplotlib.pyplot as plt
except Exception as e:
    print(f"matplotlib not available -> skipping plots ({e})")
    raise SystemExit(0)

# Bucket by alpha
by_alpha = {}
for r in records:
    by_alpha.setdefault(r["alpha"], []).append(r)

for alpha, rs in by_alpha.items():
    # order bars: variant then method
    rs = sorted(rs, key=lambda x: (x.get("variant",""), x.get("method","")))
    labels = [f"{x['method']}\n({x['variant']})" for x in rs]
    build = [x["build_pct"] for x in rs]
    count = [x["count_pct"] for x in rs]
    enum = [x["enum_pct"] for x in rs]
    sample = [x["sample_pct"] for x in rs]

    # dynamic width
    w = max(10.0, 0.7 * len(labels))
    plt.figure(figsize=(w, 5.2))
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

    out_png = fig_dir / f"exp7_phase_breakdown_alpha_{alpha}.png"
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"Wrote plot: {out_png}")
PY

log "== Done =="
echo
echo "EXP-7 outputs:"
echo "  - Sampling+Adaptive raw   : ${RAW_SA}"
if [[ -n "${RAW_ENUM}" ]]; then
  echo "  - EnumSampling raw (sparse): ${RAW_ENUM}"
fi
echo "  - Merged raw              : ${MERGED_RAW}"
echo "  - Breakdown CSV           : ${OUT_BREAKDOWN}"
echo "  - Breakdown MD            : ${OUT_MD}"
echo "  - Figures (if matplotlib) : ${FIG_DIR}/exp7_phase_breakdown_alpha_*.png"
echo
log "Tip: open the breakdown table and plots to explain who wins on which alpha and why."
