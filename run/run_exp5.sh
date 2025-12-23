#!/usr/bin/env bash
# run/run_exp5.sh
#
# EXP-5 (RQ5): Robustness to distribution shift (runtime vs t).
#
# Aligned with EXP-5.md:
#   - Fixed NR=NS (default 100k), sweep t: 1k -> 1M
#   - Compare variants: sampling + adaptive (no enum_sampling; enum_cap=0)
#   - Distributions: uniform / clustered / hetero_sizes (2D)
#
# Repo-wide runner conventions (your 5 requirements):
#   1) Logic/params match EXP-5.md
#   2) Build dir is always under: <repo_root>/build/<build_type>/
#   3) Final results are written to: <repo_root>/results/raw/exp5 (overwrite)
#   4) No embedded Python in this bash: any Python lives under <repo_root>/run/include/
#   5) All generated files (configs/logs/csv/png/...) are first written to:
#        <repo_root>/run/temp/exp5
#      then copied to results/raw/exp5 on success.
#
# Usage:
#   bash run/run_exp5.sh
#
# Common overrides (env vars):
#   EXP5_BUILD_TYPE=Release|Debug|RelWithDebInfo|MinSizeRel
#   EXP5_CLEAN_BUILD=0|1
#   EXP5_RUN_TESTS=0|1
#
#   EXP5_NR=100000
#   EXP5_NS=100000
#   EXP5_REPEATS=3
#   EXP5_THREADS=1
#   EXP5_T_LIST="1000 3000 10000 30000 100000 300000 1000000"
#   EXP5_METHODS="ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection"
#   EXP5_VARIANTS="sampling adaptive"
#
#   EXP5_GEN_SEED=1
#   EXP5_SEED=1
#   EXP5_J_STAR=1000000
#   EXP5_ENUM_CAP=0
#
set -euo pipefail

# --------------------------
# helpers
# --------------------------
log() { echo "[run_exp5] $*"; }
die() { echo "[run_exp5][FATAL] $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"
}

nproc_safe() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
  elif command -v getconf >/dev/null 2>&1; then
    getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4
  else
    echo 4
  fi
}

lower() {
  # lowercase a string (ASCII)
  echo "$1" | tr '[:upper:]' '[:lower:]'
}

build_subdir_from_type() {
  local t="$1"
  case "$t" in
    Release) echo "release";;
    Debug) echo "debug";;
    RelWithDebInfo) echo "relwithdebinfo";;
    MinSizeRel) echo "minsizerel";;
    *) lower "$t";;
  esac
}

resolve_exe() {
  # Find an executable produced by CMake, regardless of output directory layout.
  local name="$1"
  local root="$2"

  local candidates=(
    "$root/$name"
    "$root/apps/$name"
    "$root/bin/$name"
    "$root/src/apps/$name"
  )

  for p in "${candidates[@]}"; do
    if [[ -x "$p" ]]; then
      echo "$p"
      return 0
    fi
  done

  local found
  found="$(find "$root" -maxdepth 4 -type f -name "$name" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "$found" && -x "$found" ]] || return 1
  echo "$found"
}

json_array_u64() {
  local out="["
  local first=1
  for x in "$@"; do
    if [[ $first -eq 1 ]]; then first=0; else out+=", "; fi
    out+="$x"
  done
  out+="]"
  echo "$out"
}

json_array_str() {
  local out="["
  local first=1
  for x in "$@"; do
    if [[ $first -eq 1 ]]; then first=0; else out+=", "; fi
    out+="\"$x\""
  done
  out+="]"
  echo "$out"
}

write_sweep_config() {
  # Args:
  #   dist_name, generator, params_json, out_rel, cfg_path, run_tag
  local dist_name="$1"
  local generator="$2"
  local params_json="$3"
  local out_rel="$4"
  local cfg_path="$5"
  local run_tag="$6"

  local t_json methods_json variants_json
  t_json="$(json_array_u64 "${T_LIST[@]}")"
  methods_json="$(json_array_str "${METHODS[@]}")"
  variants_json="$(json_array_str "${VARIANTS[@]}")"

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
        "pbsm_scheme": "${PBSM_SCHEME}",
        "pbsm_k": ${PBSM_K},
        "pbsm_part_axis": ${PBSM_PART_AXIS},
        "pbsm_sweep": "${PBSM_SWEEP}",
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
    "output": { "out_dir": "${out_rel}", "run_tag": "${run_tag}" },
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

copy_if_missing() {
  local src="$1"
  local dst="$2"
  if [[ ! -f "$dst" ]]; then
    cp -f "$src" "$dst"
  fi
}

# --------------------------
# locate repo root
# --------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# --------------------------
# required tools
# --------------------------
need_cmd cmake
need_cmd python3

# --------------------------
# build + output dirs (convention)
# --------------------------
BUILD_TYPE="${EXP5_BUILD_TYPE:-Release}"
BUILD_SUBDIR="$(build_subdir_from_type "$BUILD_TYPE")"
BUILD_DIR="${ROOT_DIR}/build/${BUILD_SUBDIR}"

CLEAN_BUILD="${EXP5_CLEAN_BUILD:-0}"
RUN_TESTS="${EXP5_RUN_TESTS:-0}"
JOBS="${EXP5_JOBS:-$(nproc_safe)}"

TEMP_ROOT="${ROOT_DIR}/run/temp/exp5"
RESULT_ROOT="${ROOT_DIR}/results/raw/exp5"
INCLUDE_DIR="${ROOT_DIR}/run/include"

# temp layout
CFG_DIR="${TEMP_ROOT}/configs"
META_DIR="${TEMP_ROOT}/meta"
LOG_DIR="${TEMP_ROOT}/logs"

DIST_UNIFORM_DIR="${TEMP_ROOT}/uniform_t"
DIST_CLUSTERED_DIR="${TEMP_ROOT}/clustered_t"
DIST_HETERO_DIR="${TEMP_ROOT}/hetero_t"

# --------------------------
# experiment knobs (defaults = EXP-5.md)
# --------------------------
NR="${EXP5_NR:-100000}"
NS="${EXP5_NS:-100000}"
REPEATS="${EXP5_REPEATS:-3}"
THREADS="${EXP5_THREADS:-1}"

read -r -a T_LIST <<< "${EXP5_T_LIST:-1000 3000 10000 30000 100000 300000 1000000}"
read -r -a METHODS <<< "${EXP5_METHODS:-ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection}"
read -r -a VARIANTS <<< "${EXP5_VARIANTS:-sampling adaptive}"

GEN_SEED="${EXP5_GEN_SEED:-1}"
RUN_SEED="${EXP5_SEED:-1}"

J_STAR="${EXP5_J_STAR:-1000000}"
ENUM_CAP="${EXP5_ENUM_CAP:-0}"

# Baseline parameters (kept consistent with sweep_t defaults)
PBSM_SCHEME="${EXP5_PBSM_SCHEME:-stripes}"
PBSM_K="${EXP5_PBSM_K:-0}"
PBSM_PART_AXIS="${EXP5_PBSM_PART_AXIS:-0}"
PBSM_SWEEP="${EXP5_PBSM_SWEEP:-orthogonal}"

TLSOP_NX="${EXP5_TLSOP_NX:-128}"
TLSOP_NY="${EXP5_TLSOP_NY:-128}"
TLSOP_REUSE_SORT="${EXP5_TLSOP_REUSE_SORT:-true}"

SIRS_OUTER="${EXP5_SIRS_OUTER:-}"
SIRS_LEAF_SIZE="${EXP5_SIRS_LEAF_SIZE:-32}"

REJ_INDEX="${EXP5_REJ_INDEX:-S}"
REJ_REP_CENTER="${EXP5_REJ_REP_CENTER:-false}"
REJ_COUNT_DRAWS="${EXP5_REJ_COUNT_DRAWS:-50000}"
REJ_MAX_FACTOR="${EXP5_REJ_MAX_FACTOR:-10000}"

# run tag (metadata only; output directory is fixed by convention)
TAG="${EXP5_TAG:-exp5}"
RUN_ID="$(date +%Y%m%d_%H%M%S)"

# --------------------------
# fairness / reproducibility
# --------------------------
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${THREADS}}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-${THREADS}}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-${THREADS}}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-${THREADS}}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-${THREADS}}"

# --------------------------
# clean temp dir (always), prep folders
# --------------------------
rm -rf "${TEMP_ROOT}"
mkdir -p "${CFG_DIR}" "${META_DIR}" "${LOG_DIR}"
mkdir -p "${DIST_UNIFORM_DIR}" "${DIST_CLUSTERED_DIR}" "${DIST_HETERO_DIR}"

# --------------------------
# record manifest + sysinfo (in temp)
# --------------------------
{
  echo "EXP=5"
  echo "RUN_ID=${RUN_ID}"
  echo "TAG=${TAG}"
  echo "ROOT_DIR=${ROOT_DIR}"
  echo "BUILD_TYPE=${BUILD_TYPE}"
  echo "BUILD_DIR=${BUILD_DIR}"
  echo "TEMP_ROOT=${TEMP_ROOT}"
  echo "RESULT_ROOT=${RESULT_ROOT}"
  echo ""
  echo "NR=${NR}"
  echo "NS=${NS}"
  echo "REPEATS=${REPEATS}"
  echo "THREADS=${THREADS}"
  echo "T_LIST=${T_LIST[*]}"
  echo "METHODS=${METHODS[*]}"
  echo "VARIANTS=${VARIANTS[*]}"
  echo "GEN_SEED=${GEN_SEED}"
  echo "RUN_SEED=${RUN_SEED}"
  echo "J_STAR=${J_STAR}"
  echo "ENUM_CAP=${ENUM_CAP}"
  echo ""
  echo "PBSM_SCHEME=${PBSM_SCHEME}"
  echo "PBSM_K=${PBSM_K}"
  echo "PBSM_PART_AXIS=${PBSM_PART_AXIS}"
  echo "PBSM_SWEEP=${PBSM_SWEEP}"
  echo "TLSOP_NX=${TLSOP_NX}"
  echo "TLSOP_NY=${TLSOP_NY}"
  echo "TLSOP_REUSE_SORT=${TLSOP_REUSE_SORT}"
  echo "SIRS_OUTER=${SIRS_OUTER}"
  echo "SIRS_LEAF_SIZE=${SIRS_LEAF_SIZE}"
  echo "REJ_INDEX=${REJ_INDEX}"
  echo "REJ_REP_CENTER=${REJ_REP_CENTER}"
  echo "REJ_COUNT_DRAWS=${REJ_COUNT_DRAWS}"
  echo "REJ_MAX_FACTOR=${REJ_MAX_FACTOR}"
  echo ""
  echo "ENV: OMP_NUM_THREADS=${OMP_NUM_THREADS}"
  echo "ENV: MKL_NUM_THREADS=${MKL_NUM_THREADS}"
  echo "ENV: OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS}"
  echo "ENV: NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS}"
} > "${META_DIR}/manifest.txt"

{
  echo "date: $(date -Is)"
  echo "pwd:  ${ROOT_DIR}"
  uname -a || true
  echo
  cmake --version || true
  echo
  if command -v c++ >/dev/null 2>&1; then c++ --version || true; fi
  if command -v g++ >/dev/null 2>&1; then g++ --version || true; fi
  if command -v clang++ >/dev/null 2>&1; then clang++ --version || true; fi
  echo
  if command -v git >/dev/null 2>&1 && [[ -d "${ROOT_DIR}/.git" ]]; then
    echo "git_head: $(git -C "${ROOT_DIR}" rev-parse HEAD || true)"
    echo "git_status:"
    git -C "${ROOT_DIR}" status --porcelain || true
  fi
} > "${META_DIR}/sysinfo.txt"

# Save the exact runner used (for artifact completeness)
SCRIPT_PATH="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"
cp -f "${SCRIPT_PATH}" "${META_DIR}/run_exp5.sh"

# --------------------------
# build
# --------------------------
log "Repo root   : ${ROOT_DIR}"
log "Build type  : ${BUILD_TYPE}"
log "Build dir   : ${BUILD_DIR}"
log "Temp root   : ${TEMP_ROOT}"
log "Result root : ${RESULT_ROOT}"

if [[ "${CLEAN_BUILD}" == "1" ]]; then
  log "CLEAN_BUILD=1 -> rm -rf ${BUILD_DIR}"
  rm -rf "${BUILD_DIR}"
fi

mkdir -p "${BUILD_DIR}"
log "Configuring (CMake) ..."
cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" 2>&1 | tee "${LOG_DIR}/cmake_configure.log"
log "Building (-j ${JOBS}) ..."
cmake --build "${BUILD_DIR}" -j "${JOBS}" 2>&1 | tee "${LOG_DIR}/cmake_build.log"

if [[ "${RUN_TESTS}" == "1" ]]; then
  need_cmd ctest
  log "Running tests (ctest) ..."
  (cd "${BUILD_DIR}" && ctest --output-on-failure) 2>&1 | tee "${LOG_DIR}/ctest.log"
fi

SJS_SWEEP="$(resolve_exe sjs_sweep "${BUILD_DIR}" || true)"
[[ -n "${SJS_SWEEP}" ]] || die "Cannot find executable: sjs_sweep under ${BUILD_DIR}"
log "Using sjs_sweep: ${SJS_SWEEP}"

# --------------------------
# write configs (in temp)
# --------------------------
# Distribution params (must match EXP-5.md defaults)
UNIFORM_PARAMS='{"w_min":0.005, "w_max":0.02, "same_size_all_dims":false, "shuffle_r":false, "shuffle_s":false}'
CLUSTERED_PARAMS='{"clusters":10, "sigma":0.05, "share_clusters":true, "w_min":0.003, "w_max":0.02, "shuffle_r":false, "shuffle_s":false}'
HETERO_PARAMS='{"anisotropic":true, "p_large":0.1, "w_small_min":0.002, "w_small_max":0.01, "w_large_min":0.1, "w_large_max":0.3, "shuffle_r":false, "shuffle_s":false}'

CFG_UNIFORM="${CFG_DIR}/exp5_t_uniform.json"
CFG_CLUSTERED="${CFG_DIR}/exp5_t_clustered.json"
CFG_HETERO="${CFG_DIR}/exp5_t_hetero.json"

# out_dir is a repo-relative path (so configs remain portable).
OUT_UNIFORM_REL="run/temp/exp5/uniform_t"
OUT_CLUSTERED_REL="run/temp/exp5/clustered_t"
OUT_HETERO_REL="run/temp/exp5/hetero_t"

write_sweep_config "uniform"   "uniform"      "${UNIFORM_PARAMS}"   "${OUT_UNIFORM_REL}"   "${CFG_UNIFORM}"   "${TAG}_uniform_t"
write_sweep_config "clustered" "clustered"    "${CLUSTERED_PARAMS}" "${OUT_CLUSTERED_REL}" "${CFG_CLUSTERED}" "${TAG}_clustered_t"
write_sweep_config "hetero"    "hetero_sizes" "${HETERO_PARAMS}"    "${OUT_HETERO_REL}"    "${CFG_HETERO}"    "${TAG}_hetero_t"

# --------------------------
# run sweeps (output goes into run/temp/exp5/<dist>_t)
# --------------------------
cd "${ROOT_DIR}"

log "Running sweep: uniform ..."
"${SJS_SWEEP}" --config="${CFG_UNIFORM}" 2>&1 | tee "${DIST_UNIFORM_DIR}/sweep.log"
copy_if_missing "${CFG_UNIFORM}" "${DIST_UNIFORM_DIR}/sweep_config.json"

log "Running sweep: clustered ..."
"${SJS_SWEEP}" --config="${CFG_CLUSTERED}" 2>&1 | tee "${DIST_CLUSTERED_DIR}/sweep.log"
copy_if_missing "${CFG_CLUSTERED}" "${DIST_CLUSTERED_DIR}/sweep_config.json"

log "Running sweep: hetero_sizes ..."
"${SJS_SWEEP}" --config="${CFG_HETERO}" 2>&1 | tee "${DIST_HETERO_DIR}/sweep.log"
copy_if_missing "${CFG_HETERO}" "${DIST_HETERO_DIR}/sweep_config.json"

# Basic output checks (fail fast if sweeps didn't produce expected files)
for d in "${DIST_UNIFORM_DIR}" "${DIST_CLUSTERED_DIR}" "${DIST_HETERO_DIR}"; do
  [[ -f "${d}/sweep_raw.csv" ]] || die "Missing: ${d}/sweep_raw.csv"
  [[ -f "${d}/sweep_summary.csv" ]] || die "Missing: ${d}/sweep_summary.csv"
done

# --------------------------
# plots (python code is external, under run/include)
# --------------------------
PLOT_PY="${INCLUDE_DIR}/exp5_plot_runtime_t.py"
[[ -f "${PLOT_PY}" ]] || die "Missing plot script: ${PLOT_PY}"

log "Plotting runtime-vs-t figures ..."
python3 "${PLOT_PY}" --summary_csv "${DIST_UNIFORM_DIR}/sweep_summary.csv"   --out_png "${DIST_UNIFORM_DIR}/runtime_t.png"   --title "EXP-5 runtime vs t (Uniform)"
python3 "${PLOT_PY}" --summary_csv "${DIST_CLUSTERED_DIR}/sweep_summary.csv" --out_png "${DIST_CLUSTERED_DIR}/runtime_t.png" --title "EXP-5 runtime vs t (Clustered)"
python3 "${PLOT_PY}" --summary_csv "${DIST_HETERO_DIR}/sweep_summary.csv"    --out_png "${DIST_HETERO_DIR}/runtime_t.png"    --title "EXP-5 runtime vs t (Hetero-sizes)"

# --------------------------
# sync to results/raw/exp5 (overwrite on success)
# --------------------------
log "Syncing results -> ${RESULT_ROOT} (overwrite)"
rm -rf "${RESULT_ROOT}"
mkdir -p "${RESULT_ROOT}"
cp -a "${TEMP_ROOT}/." "${RESULT_ROOT}/"

log "DONE."
log "Results written to: ${RESULT_ROOT}"
log "Key files:"
log "  ${RESULT_ROOT}/uniform_t/sweep_summary.csv"
log "  ${RESULT_ROOT}/clustered_t/sweep_summary.csv"
log "  ${RESULT_ROOT}/hetero_t/sweep_summary.csv"
