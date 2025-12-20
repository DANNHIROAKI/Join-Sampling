#!/usr/bin/env bash
# run/run_exp5.sh
#
# EXP-5 (RQ5): Distribution robustness.
# - Switch synthetic distribution to: uniform / clustered / hetero_sizes
# - Repeat the EXP-2 style sweep: runtime vs t (sampling + adaptive)
#
# Output:
#   results/sweeps/exp5/<TAG>/{uniform,clustered,hetero}_t/sweep_{raw,summary}.csv
#   results/sweeps/exp5/<TAG>/{uniform,clustered,hetero}_t/runtime_t.png
#
# Usage:
#   bash run/run_exp5.sh
#
# Optional env overrides:
#   EXP5_TAG=exp5_mylabel
#   EXP5_BUILD_TYPE=Release|RelWithDebInfo
#   EXP5_NR=100000
#   EXP5_NS=100000
#   EXP5_REPEATS=3
#   EXP5_THREADS=1
#   EXP5_T_LIST="1000 3000 10000 30000 100000 300000 1000000"
#   EXP5_METHODS="ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection"
#   EXP5_VARIANTS="sampling adaptive"
#   EXP5_RUN_TESTS=0|1
#
set -euo pipefail

# --------------------------
# Repo root + basic paths
# --------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

BUILD_DIR="${ROOT_DIR}/build"
BUILD_TYPE="${EXP5_BUILD_TYPE:-Release}"

TAG="${EXP5_TAG:-exp5_$(date +%Y%m%d_%H%M%S)}"
OUT_BASE_REL="results/sweeps/exp5/${TAG}"
OUT_BASE_ABS="${ROOT_DIR}/${OUT_BASE_REL}"

mkdir -p "${OUT_BASE_ABS}"

# --------------------------
# Experiment knobs (defaults match config/sweeps/sweep_t.json where applicable)
# --------------------------
NR="${EXP5_NR:-100000}"
NS="${EXP5_NS:-100000}"
REPEATS="${EXP5_REPEATS:-3}"
THREADS="${EXP5_THREADS:-1}"

# Default t list (EXP-2 style)
read -r -a T_LIST <<< "${EXP5_T_LIST:-1000 3000 10000 30000 100000 300000 1000000}"

# Default method list (Dim=2)
read -r -a METHODS <<< "${EXP5_METHODS:-ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection}"

# For EXP-5 we focus on the sampling pipelines (sampling + adaptive).
read -r -a VARIANTS <<< "${EXP5_VARIANTS:-sampling adaptive}"

RUN_TESTS="${EXP5_RUN_TESTS:-0}"

# Base seeds
GEN_SEED="${EXP5_GEN_SEED:-1}"
RUN_SEED="${EXP5_SEED:-1}"

# Adaptive threshold (used by variant=adaptive)
J_STAR="${EXP5_J_STAR:-1000000}"

# Enumeration cap (kept at 0 since we do not run enum_sampling here)
ENUM_CAP="${EXP5_ENUM_CAP:-0}"

# Extra baseline parameters (copied from config/sweeps/sweep_t.json)
PDSM_SCHEME="${EXP5_PBSM_SCHEME:-stripes}"
PDSM_K="${EXP5_PBSM_K:-0}"
PDSM_PART_AXIS="${EXP5_PBSM_PART_AXIS:-0}"
PDSM_SWEEP="${EXP5_PBSM_SWEEP:-orthogonal}"

TLSOP_NX="${EXP5_TLSOP_NX:-128}"
TLSOP_NY="${EXP5_TLSOP_NY:-128}"
TLSOP_REUSE_SORT="${EXP5_TLSOP_REUSE_SORT:-true}"

SIRS_OUTER="${EXP5_SIRS_OUTER:-}"
SIRS_LEAF_SIZE="${EXP5_SIRS_LEAF_SIZE:-32}"

REJ_INDEX="${EXP5_REJ_INDEX:-S}"
REJ_REP_CENTER="${EXP5_REJ_REP_CENTER:-false}"
REJ_COUNT_DRAWS="${EXP5_REJ_COUNT_DRAWS:-50000}"
REJ_MAX_FACTOR="${EXP5_REJ_MAX_FACTOR:-10000}"

# --------------------------
# Helpers
# --------------------------
log() { echo "[run_exp5] $*"; }

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: required command not found: $1" >&2
    exit 1
  fi
}

find_exe() {
  local name="$1"
  local candidates=(
    "${BUILD_DIR}/${name}"
    "${BUILD_DIR}/apps/${name}"
    "${BUILD_DIR}/bin/${name}"
    "${BUILD_DIR}/${BUILD_TYPE}/${name}"
  )
  for p in "${candidates[@]}"; do
    if [[ -x "$p" ]]; then
      echo "$p"
      return 0
    fi
  done
  # Fallback: search a bit (kept shallow so it doesn't crawl forever).
  local p
  p="$(find "${BUILD_DIR}" -maxdepth 3 -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  if [[ -n "${p}" && -x "${p}" ]]; then
    echo "$p"
    return 0
  fi
  return 1
}

json_array_u64() {
  # shellcheck disable=SC2206
  local arr=($*)
  local out="["
  local first=1
  for x in "${arr[@]}"; do
    if [[ "${first}" -eq 1 ]]; then
      first=0
    else
      out+=", "
    fi
    out+="${x}"
  done
  out+="]"
  echo "${out}"
}

json_array_str() {
  # shellcheck disable=SC2206
  local arr=($*)
  local out="["
  local first=1
  for x in "${arr[@]}"; do
    if [[ "${first}" -eq 1 ]]; then
      first=0
    else
      out+=", "
    fi
    out+="\"${x}\""
  done
  out+="]"
  echo "${out}"
}

write_sweep_config() {
  local dist_name="$1"      # e.g., uniform
  local generator="$2"      # e.g., uniform
  local params_json="$3"    # JSON object (string) for dataset.synthetic.params
  local out_rel="$4"        # relative out_dir for this sweep
  local cfg_path="$5"       # where to write the JSON config

  local t_json methods_json variants_json
  t_json="$(json_array_u64 "${T_LIST[*]}")"
  methods_json="$(json_array_str "${METHODS[*]}")"
  variants_json="$(json_array_str "${VARIANTS[*]}")"

  cat > "${cfg_path}" <<EOF
{
  "base": {
    "dataset": {
      "source": "synthetic",
      "name": "${dist_name}_t",
      "dim": 2,
      "synthetic": {
        "generator": "${generator}",
        "n_r": ${NR},
        "n_s": ${NS},
        "alpha": 1.0,
        "seed": ${GEN_SEED},
        "params": ${params_json}
      }
    },
    "run": {
      "method": "ours",
      "variant": "sampling",
      "t": 10000,
      "seed": ${RUN_SEED},
      "repeats": ${REPEATS},
      "j_star": ${J_STAR},
      "enum_cap": ${ENUM_CAP},
      "write_samples": false,
      "verify": false,
      "extra": {
        "pbsm_scheme": "${PDSM_SCHEME}",
        "pbsm_k": ${PDSM_K},
        "pbsm_part_axis": ${PDSM_PART_AXIS},
        "pbsm_sweep": "${PDSM_SWEEP}",
        "tlsop_nx": ${TLSOP_NX},
        "tlsop_ny": ${TLSOP_NY},
        "tlsop_reuse_sort": ${TLSOP_REUSE_SORT},
        "sirs_outer": "${SIRS_OUTER}",
        "sirs_leaf_size": ${SIRS_LEAF_SIZE},
        "rej_index": "${REJ_INDEX}",
        "rej_rep_center": ${REJ_REP_CENTER},
        "rej_count_draws": ${REJ_COUNT_DRAWS},
        "rej_max_factor": ${REJ_MAX_FACTOR}
      }
    },
    "output": { "out_dir": "${out_rel}", "run_tag": "${TAG}_${dist_name}_t" },
    "logging": { "level": "info", "with_timestamp": true, "with_thread_id": false },
    "sys": { "threads": ${THREADS} }
  },

  "sweep": {
    "t": ${t_json},
    "method": ${methods_json},
    "variant": ${variants_json}
  },

  "files": { "raw": "sweep_raw.csv", "summary": "sweep_summary.csv" },
  "meta": { "note": "EXP-5 (RQ5): distribution robustness; runtime vs t; dist=${dist_name}" }
}
EOF
}

plot_runtime_t() {
  local summary_csv="$1"
  local out_png="$2"

  python3 - "$summary_csv" "$out_png" <<'PY'
import sys
import pandas as pd
import matplotlib.pyplot as plt

summary_csv = sys.argv[1]
out_png = sys.argv[2]

df = pd.read_csv(summary_csv)

# Keep only successful rows
if "ok_rate" in df.columns:
    df = df[df["ok_rate"] > 0.0]

# Some rows are placeholders with notes ("SKIPPED: ...")
if "note" in df.columns:
    df = df[~df["note"].fillna("").str.startswith("SKIPPED")]

if df.empty:
    print("WARN: empty dataframe after filtering; no plot written:", out_png)
    raise SystemExit(0)

# Each curve = method + variant
for (method, variant), g in df.groupby(["method", "variant"], dropna=False):
    g = g.sort_values("t")
    plt.plot(g["t"], g["wall_median_ms"], label=f"{method}-{variant}")

plt.xscale("log")
plt.yscale("log")
plt.xlabel("t (#samples)")
plt.ylabel("median wall time (ms)")
plt.title("EXP-5 runtime vs t")
plt.legend(fontsize=7)
plt.tight_layout()
plt.savefig(out_png, dpi=200)
print("Wrote:", out_png)
PY
}


# --------------------------
# Main
# --------------------------
need_cmd cmake
need_cmd python3

log "Repo root: ${ROOT_DIR}"
log "Build dir: ${BUILD_DIR} (type=${BUILD_TYPE})"
log "Output base: ${OUT_BASE_ABS}"
log "NR=${NR}, NS=${NS}, repeats=${REPEATS}, threads=${THREADS}"
log "t_list: ${T_LIST[*]}"
log "methods: ${METHODS[*]}"
log "variants: ${VARIANTS[*]}"
log "tag: ${TAG}"

cd "${ROOT_DIR}"

# Build (from scratch)
log "Configuring..."
cmake -S . -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
log "Building..."
cmake --build "${BUILD_DIR}" -j

# Optional unit tests
if [[ "${RUN_TESTS}" == "1" ]]; then
  log "Running tests (ctest)..."
  (cd "${BUILD_DIR}" && ctest --output-on-failure)
fi

# Locate the sweep executable
SJS_SWEEP="$(find_exe sjs_sweep || true)"
if [[ -z "${SJS_SWEEP}" ]]; then
  echo "ERROR: cannot find sjs_sweep in ${BUILD_DIR}" >&2
  exit 2
fi
log "Using sjs_sweep: ${SJS_SWEEP}"

# Prepare sweep configs (written into the output folder for reproducibility)
CFG_DIR_ABS="${OUT_BASE_ABS}/configs"
mkdir -p "${CFG_DIR_ABS}"

# Distribution params (kept explicit for reproducibility)
UNIFORM_PARAMS='{"w_min":0.005, "w_max":0.02, "same_size_all_dims":false, "shuffle_r":false, "shuffle_s":false}'
CLUSTERED_PARAMS='{"clusters":10, "sigma":0.05, "share_clusters":true, "w_min":0.003, "w_max":0.02, "shuffle_r":false, "shuffle_s":false}'
HETERO_PARAMS='{"anisotropic":true, "p_large":0.1, "w_small_min":0.002, "w_small_max":0.01, "w_large_min":0.1, "w_large_max":0.3, "shuffle_r":false, "shuffle_s":false}'

# Write configs
CFG_UNIFORM="${CFG_DIR_ABS}/exp5_t_uniform.json"
CFG_CLUSTERED="${CFG_DIR_ABS}/exp5_t_clustered.json"
CFG_HETERO="${CFG_DIR_ABS}/exp5_t_hetero.json"

write_sweep_config "uniform"    "uniform"      "${UNIFORM_PARAMS}"   "${OUT_BASE_REL}/uniform_t"   "${CFG_UNIFORM}"
write_sweep_config "clustered"  "clustered"    "${CLUSTERED_PARAMS}" "${OUT_BASE_REL}/clustered_t" "${CFG_CLUSTERED}"
write_sweep_config "hetero"     "hetero_sizes" "${HETERO_PARAMS}"    "${OUT_BASE_REL}/hetero_t"    "${CFG_HETERO}"

# Run sweeps
log "Running sweep: uniform ..."
mkdir -p "${OUT_BASE_ABS}/uniform_t"
"${SJS_SWEEP}" --config="${CFG_UNIFORM}" 2>&1 | tee "${OUT_BASE_ABS}/uniform_t/sweep.log"

log "Running sweep: clustered ..."
mkdir -p "${OUT_BASE_ABS}/clustered_t"
"${SJS_SWEEP}" --config="${CFG_CLUSTERED}" 2>&1 | tee "${OUT_BASE_ABS}/clustered_t/sweep.log"

log "Running sweep: hetero_sizes ..."
mkdir -p "${OUT_BASE_ABS}/hetero_t"
"${SJS_SWEEP}" --config="${CFG_HETERO}" 2>&1 | tee "${OUT_BASE_ABS}/hetero_t/sweep.log"

# Plot runtime vs t (one figure per distribution)
log "Plotting runtime-vs-t figures ..."
plot_runtime_t "${OUT_BASE_ABS}/uniform_t/sweep_summary.csv"   "${OUT_BASE_ABS}/uniform_t/runtime_t.png"
plot_runtime_t "${OUT_BASE_ABS}/clustered_t/sweep_summary.csv" "${OUT_BASE_ABS}/clustered_t/runtime_t.png"
plot_runtime_t "${OUT_BASE_ABS}/hetero_t/sweep_summary.csv"    "${OUT_BASE_ABS}/hetero_t/runtime_t.png"

log "DONE."
log "Results:"
log "  ${OUT_BASE_ABS}/uniform_t/sweep_summary.csv"
log "  ${OUT_BASE_ABS}/clustered_t/sweep_summary.csv"
log "  ${OUT_BASE_ABS}/hetero_t/sweep_summary.csv"
log "Figures:"
log "  ${OUT_BASE_ABS}/uniform_t/runtime_t.png"
log "  ${OUT_BASE_ABS}/clustered_t/runtime_t.png"
log "  ${OUT_BASE_ABS}/hetero_t/runtime_t.png"
