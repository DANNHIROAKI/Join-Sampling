#!/usr/bin/env bash
# run/run_exp6.sh
#
# EXP-6 (RQ6): Adaptivity effectiveness on an (alpha, t) grid.
#
# Aligned with EXP-6.md:
#   - Data: synthetic "stripe" (alias of stripe_ctrl_alpha), where
#       alpha = |J| / (n_r + n_s)  and  |J| = round(alpha * (n_r+n_s))
#   - Sweep: a grid over (alpha, t)
#   - Compare variants (fixed): sampling vs enum_sampling vs adaptive
#   - Metric: end-to-end wall time (summary uses median/p95 + ok_rate)
#   - Artifacts: sweep_raw.csv, sweep_summary.csv, phase/ratio heatmaps
#
# Repo-wide runner conventions:
#   1) Experiment logic/params align with EXP-6.md (RQ6 intent)
#   2) Build output under: <repo_root>/build/<build_type>/
#   3) Final results under: <repo_root>/results/raw/exp6  (OVERWRITE each run) [default]
#   4) No embedded python in this bash: python lives under <repo_root>/run/include/
#   5) All generated configs/logs/etc are written to: <repo_root>/run/temp/exp6
#      and copied to results/raw/exp6 on success.
#
# Usage:
#   bash run/run_exp6.sh
#
# Common overrides (env vars):
#   EXP6_BUILD_TYPE=Release|Debug|RelWithDebInfo|MinSizeRel
#   EXP6_CLEAN_BUILD=0|1
#   EXP6_RUN_TESTS=0|1
#   EXP6_BUILD_JOBS=8
#
#   # Dataset size (NOTE: EXP6_N sets each side, i.e., n_r=n_s=EXP6_N; total=2*EXP6_N)
#   EXP6_N=100000
#   EXP6_N_R=100000 EXP6_N_S=100000
#   EXP6_GEN_SEED=1
#
#   # Grid sweep
#   EXP6_ALPHAS="0.01,0.03,0.1,0.3,1,3,10,30"
#   EXP6_TS="1000,3000,10000,30000,100000,300000,1000000"
#
#   # Which methods to include (comma-separated)
#   EXP6_METHODS="ours,kd_tree,r_tree"
#
#   # Run controls
#   EXP6_REPEATS=3
#   EXP6_RUN_SEED=1
#   EXP6_THREADS=1
#   EXP6_J_STAR=1000000
#   EXP6_ENUM_CAP=0
#
#   # Output dirs (optional)
#   EXP6_OUT_DIR="results/raw/exp6"       # relative-to-repo-root or absolute
#   EXP6_TEMP_DIR="run/temp/exp6"         # relative-to-repo-root or absolute
#
set -euo pipefail
IFS=$'\n\t'

trap 'echo -e "[EXP6][FATAL] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

# --------------------------
# helpers
# --------------------------
log() { echo "[EXP6] $*"; }
die() { echo "[EXP6][FATAL] $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"
}

nproc_safe() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
  elif command -v getconf >/dev/null 2>&1; then
    getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4
  elif [[ "$(uname -s)" == "Darwin" ]] && command -v sysctl >/dev/null 2>&1; then
    sysctl -n hw.ncpu 2>/dev/null || echo 4
  else
    echo 4
  fi
}

lower() {
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

trim_spaces() {
  echo "$1" | tr -d '[:space:]'
}

json_bool() {
  case "$(trim_spaces "$1")" in
    1|true|TRUE|True|yes|YES|y|Y) echo "true";;
    0|false|FALSE|False|no|NO|n|N) echo "false";;
    *) echo "false";;
  esac
}

json_array_numbers_from_csv() {
  local s
  s="$(trim_spaces "$1")"
  IFS=',' read -r -a arr <<< "$s"
  local out="["
  local first=1
  for v in "${arr[@]}"; do
    [[ -z "$v" ]] && continue
    if [[ $first -eq 0 ]]; then out+=", "; fi
    out+="$v"
    first=0
  done
  out+="]"
  echo "$out"
}

json_array_strings_from_csv() {
  local s
  s="$(trim_spaces "$1")"
  IFS=',' read -r -a arr <<< "$s"
  local out="["
  local first=1
  for v in "${arr[@]}"; do
    [[ -z "$v" ]] && continue
    if [[ $first -eq 0 ]]; then out+=", "; fi
    out+="\"$v\""
    first=0
  done
  out+="]"
  echo "$out"
}

resolve_exe() {
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

abs_path_under_root() {
  local p="$1"
  local root="$2"
  if [[ "$p" == /* ]]; then
    echo "$p"
  else
    echo "${root}/${p}"
  fi
}

iso_timestamp() {
  # Portable ISO timestamp for Linux/macOS
  date +"%Y-%m-%dT%H:%M:%S%z"
}

# --------------------------
# locate repo root
# --------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# --------------------------
# deps
# --------------------------
need_cmd cmake
need_cmd python3
need_cmd tee
need_cmd find
need_cmd awk

# --------------------------
# parameters (override via env)
# --------------------------

# Build
EXP6_BUILD_TYPE="${EXP6_BUILD_TYPE:-Release}"
EXP6_BUILD_JOBS="${EXP6_BUILD_JOBS:-$(nproc_safe)}"
EXP6_CLEAN_BUILD="${EXP6_CLEAN_BUILD:-0}"
EXP6_RUN_TESTS="${EXP6_RUN_TESTS:-0}"

BUILD_SUBDIR="$(build_subdir_from_type "${EXP6_BUILD_TYPE}")"
BUILD_DIR="${EXP6_BUILD_DIR:-${ROOT_DIR}/build/${BUILD_SUBDIR}}"

# Dataset (synthetic stripe)
EXP6_N="${EXP6_N:-}"
EXP6_N_R="${EXP6_N_R:-100000}"
EXP6_N_S="${EXP6_N_S:-100000}"
if [[ -n "${EXP6_N}" ]]; then
  EXP6_N_R="${EXP6_N}"
  EXP6_N_S="${EXP6_N}"
fi
EXP6_GEN_SEED="${EXP6_GEN_SEED:-1}"

# Generator params (stripe_ctrl_alpha)
EXP6_CONTROL_AXIS="${EXP6_CONTROL_AXIS:-1}"
EXP6_CORE_LO="${EXP6_CORE_LO:-0.45}"
EXP6_CORE_HI="${EXP6_CORE_HI:-0.55}"
EXP6_GAP_FACTOR="${EXP6_GAP_FACTOR:-0.10}"
EXP6_DELTA_FACTOR="${EXP6_DELTA_FACTOR:-0.25}"
EXP6_SHUFFLE_STRIPS="${EXP6_SHUFFLE_STRIPS:-true}"
EXP6_SHUFFLE_R="${EXP6_SHUFFLE_R:-false}"
EXP6_SWAP_SIDES="${EXP6_SWAP_SIDES:-false}"

# Sweep grid
EXP6_ALPHAS="${EXP6_ALPHAS:-0.01,0.03,0.1,0.3,1,3,10,30}"
EXP6_TS="${EXP6_TS:-1000,3000,10000,30000,100000,300000,1000000}"

# Methods to run
EXP6_METHODS="${EXP6_METHODS:-ours}"

# Variants are fixed for EXP-6; allow override only for debugging.
EXP6_VARIANTS="${EXP6_VARIANTS:-sampling,enum_sampling,adaptive}"

# Run control
EXP6_REPEATS="${EXP6_REPEATS:-3}"
EXP6_RUN_SEED="${EXP6_RUN_SEED:-1}"
EXP6_THREADS="${EXP6_THREADS:-1}"

# Adaptive / enumeration knobs
EXP6_J_STAR="${EXP6_J_STAR:-1000000}"
EXP6_ENUM_CAP="${EXP6_ENUM_CAP:-0}"

if [[ "${EXP6_ENUM_CAP}" != "0" ]]; then
  log "[WARN] EXP6_ENUM_CAP=${EXP6_ENUM_CAP} != 0. This can truncate enum_sampling and may invalidate oracle(min) unless you exclude truncated points in plotting."
fi

# --------------------------
# directories (default matches repo convention; allow override via env)
# --------------------------
EXP6_TEMP_DIR="${EXP6_TEMP_DIR:-run/temp/exp6}"
EXP6_OUT_DIR="${EXP6_OUT_DIR:-results/raw/exp6}"

TEMP_ROOT="$(abs_path_under_root "${EXP6_TEMP_DIR}" "${ROOT_DIR}")"
FINAL_OUT="$(abs_path_under_root "${EXP6_OUT_DIR}" "${ROOT_DIR}")"

META_DIR="${TEMP_ROOT}/meta"
LOG_DIR="${TEMP_ROOT}/logs"
FIGS_DIR="${TEMP_ROOT}/figs"

PLOT_HELPER="${ROOT_DIR}/run/include/exp6_plot.py"

# Clean temp every run (keeps runs deterministic + avoids stale artifacts)
rm -rf "${TEMP_ROOT}"
mkdir -p "${META_DIR}" "${LOG_DIR}" "${FIGS_DIR}"

# Ensure results/raw exists
mkdir -p "${ROOT_DIR}/results/raw"

# --------------------------
# fairness / reproducibility
# --------------------------
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${EXP6_THREADS}}"
export MKL_NUM_THREADS="${MKO_NUM_THREADS:-${EXP6_THREADS}}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-${EXP6_THREADS}}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-${EXP6_THREADS}}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-${EXP6_THREADS}}"

# --------------------------
# manifest + sysinfo
# --------------------------
{
  echo "timestamp=$(iso_timestamp)"
  echo "root_dir=${ROOT_DIR}"
  echo "build_type=${EXP6_BUILD_TYPE}"
  echo "build_dir=${BUILD_DIR}"
  echo "methods=${EXP6_METHODS}"
  echo "variants=${EXP6_VARIANTS}"
  echo "alphas=${EXP6_ALPHAS}"
  echo "t_list=${EXP6_TS}"
  echo "n_r=${EXP6_N_R}"
  echo "n_s=${EXP6_N_S}"
  echo "gen_seed=${EXP6_GEN_SEED}"
  echo "run_seed=${EXP6_RUN_SEED}"
  echo "repeats=${EXP6_REPEATS}"
  echo "threads=${EXP6_THREADS}"
  echo "j_star=${EXP6_J_STAR}"
  echo "enum_cap=${EXP6_ENUM_CAP}"
  echo "stripe.control_axis=${EXP6_CONTROL_AXIS}"
  echo "stripe.core_lo=${EXP6_CORE_LO}"
  echo "stripe.core_hi=${EXP6_CORE_HI}"
  echo "stripe.gap_factor=${EXP6_GAP_FACTOR}"
  echo "stripe.delta_factor=${EXP6_DELTA_FACTOR}"
  echo "stripe.shuffle_strips=${EXP6_SHUFFLE_STRIPS}"
  echo "stripe.shuffle_r=${EXP6_SHUFFLE_R}"
  echo "stripe.swap_sides=${EXP6_SWAP_SIDES}"
  if command -v git >/dev/null 2>&1 && git -C "${ROOT_DIR}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "git_sha=$(git -C "${ROOT_DIR}" rev-parse HEAD)"
    echo "git_dirty=$(git -C "${ROOT_DIR}" status --porcelain | wc -l | awk '{print $1}')"
  else
    echo "git_sha=unknown"
    echo "git_dirty=unknown"
  fi
} > "${META_DIR}/manifest.txt"

{
  echo "date: $(iso_timestamp)"
  uname -a || true
  echo
  cmake --version || true
  echo
  if command -v g++ >/dev/null 2>&1; then g++ --version || true; fi
  if command -v clang++ >/dev/null 2>&1; then clang++ --version || true; fi
} > "${META_DIR}/sysinfo.txt"

# --------------------------
# build
# --------------------------
if [[ "${EXP6_CLEAN_BUILD}" == "1" ]]; then
  log "Cleaning build dir: ${BUILD_DIR}"
  rm -rf "${BUILD_DIR}"
fi

log "Configuring (CMake, ${EXP6_BUILD_TYPE}) -> ${BUILD_DIR}"
cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE="${EXP6_BUILD_TYPE}" \
  2>&1 | tee "${LOG_DIR}/cmake_configure.log"

log "Building (jobs=${EXP6_BUILD_JOBS})"
cmake --build "${BUILD_DIR}" -j "${EXP6_BUILD_JOBS}" \
  2>&1 | tee "${LOG_DIR}/cmake_build.log"

if [[ "${EXP6_RUN_TESTS}" == "1" ]]; then
  log "Running tests (ctest)"
  (cd "${BUILD_DIR}" && ctest --output-on-failure) 2>&1 | tee "${LOG_DIR}/ctest.log"
fi

SJS_SWEEP="$(resolve_exe sjs_sweep "${BUILD_DIR}" || true)"
[[ -n "${SJS_SWEEP}" ]] || die "Could not find executable 'sjs_sweep' under: ${BUILD_DIR}"
log "Using sjs_sweep: ${SJS_SWEEP}"

# --------------------------
# Sanity check (FIXED): ensure sjs_sweep supports JSON config (--config)
# NOTE: do NOT use `help | grep -q` under pipefail; it can false-fail due to SIGPIPE.
# --------------------------
HELP_TXT="$("${SJS_SWEEP}" --help 2>&1 || true)"
if [[ "${HELP_TXT}" != *"--config"* ]]; then
  die "The resolved sjs_sweep does not appear to mention --config in its help output. If you built the legacy src/apps sweep binary, rebuild the root apps version (apps/sjs_sweep.cpp) or adjust the runner."
fi

# --------------------------
# generate sweep config (JSON) -> TEMP_ROOT
# --------------------------
ALPHAS_JSON="$(json_array_numbers_from_csv "${EXP6_ALPHAS}")"
TS_JSON="$(json_array_numbers_from_csv "${EXP6_TS}")"
METHODS_JSON="$(json_array_strings_from_csv "${EXP6_METHODS}")"
VARIANTS_JSON="$(json_array_strings_from_csv "${EXP6_VARIANTS}")"

SHUFFLE_STRIPS_JSON="$(json_bool "${EXP6_SHUFFLE_STRIPS}")"
SHUFFLE_R_JSON="$(json_bool "${EXP6_SHUFFLE_R}")"
SWAP_SIDES_JSON="$(json_bool "${EXP6_SWAP_SIDES}")"

CONFIG_PATH="${TEMP_ROOT}/exp6_alpha_t.json"

cat > "${CONFIG_PATH}" <<EOF
{
  "base": {
    "dataset": {
      "source": "synthetic",
      "name": "exp6_stripe",
      "dim": 2,
      "synthetic": {
        "generator": "stripe",
        "n_r": ${EXP6_N_R},
        "n_s": ${EXP6_N_S},
        "alpha": 0.1,
        "seed": ${EXP6_GEN_SEED},
        "params": {
          "control_axis": ${EXP6_CONTROL_AXIS},
          "core_lo": ${EXP6_CORE_LO},
          "core_hi": ${EXP6_CORE_HI},
          "gap_factor": ${EXP6_GAP_FACTOR},
          "delta_factor": ${EXP6_DELTA_FACTOR},
          "shuffle_strips": ${SHUFFLE_STRIPS_JSON},
          "shuffle_r": ${SHUFFLE_R_JSON},
          "swap_sides": ${SWAP_SIDES_JSON}
        }
      }
    },

    "run": {
      "method": "ours",
      "variant": "sampling",
      "t": 100000,
      "seed": ${EXP6_RUN_SEED},
      "repeats": ${EXP6_REPEATS},
      "enum_cap": ${EXP6_ENUM_CAP},
      "j_star": ${EXP6_J_STAR},
      "write_samples": false,
      "verify": false,
      "extra": {}
    },

    "output": {
      "out_dir": "${TEMP_ROOT}"
    },

    "logging": {
      "level": "info"
    },

    "sys": {
      "threads": ${EXP6_THREADS}
    }
  },

  "sweep": {
    "alpha": ${ALPHAS_JSON},
    "t": ${TS_JSON},
    "method": ${METHODS_JSON},
    "variant": ${VARIANTS_JSON}
  },

  "files": {
    "raw": "sweep_raw.csv",
    "summary": "sweep_summary.csv"
  }
}
EOF

log "Wrote sweep config: ${CONFIG_PATH}"

# --------------------------
# run sweep
# --------------------------
SWEEP_LOG="${LOG_DIR}/sjs_sweep.log"
log "Running sweep"
log "Log: ${SWEEP_LOG}"

(
  cd "${ROOT_DIR}"
  "${SJS_SWEEP}" --config="${CONFIG_PATH}"
) 2>&1 | tee "${SWEEP_LOG}"

RAW_CSV="${TEMP_ROOT}/sweep_raw.csv"
SUMMARY_CSV="${TEMP_ROOT}/sweep_summary.csv"
[[ -f "${RAW_CSV}" ]] || die "Missing expected output: ${RAW_CSV}"
[[ -f "${SUMMARY_CSV}" ]] || die "Missing expected output: ${SUMMARY_CSV}"

log "Sweep done."
log "  raw:     ${RAW_CSV}"
log "  summary: ${SUMMARY_CSV}"

# --------------------------
# plot (phase + ratio + adaptive_branch heatmap)
# --------------------------
[[ -f "${PLOT_HELPER}" ]] || die "Missing plot helper: ${PLOT_HELPER} (expected under run/include/)"

log "Plotting via: ${PLOT_HELPER}"
python3 "${PLOT_HELPER}" \
  --summary_csv "${SUMMARY_CSV}" \
  --raw_csv "${RAW_CSV}" \
  --out_dir "${FIGS_DIR}" \
  2>&1 | tee "${LOG_DIR}/plot.log"

# --------------------------
# finalize -> FINAL_OUT (overwrite)
# --------------------------
log "Copying artifacts to final results dir (overwrite): ${FINAL_OUT}"
rm -rf "${FINAL_OUT}"
mkdir -p "${FINAL_OUT}"
cp -a "${TEMP_ROOT}/." "${FINAL_OUT}/"

log "DONE ✅"
log "Final results: ${FINAL_OUT}"
log "Key files:"
log "  ${FINAL_OUT}/sweep_raw.csv"
log "  ${FINAL_OUT}/sweep_summary.csv"
log "  ${FINAL_OUT}/figs/exp6_phase_<method>.png (if matplotlib)"
log "  ${FINAL_OUT}/figs/exp6_ratio_<method>.png (if matplotlib)"
log "  ${FINAL_OUT}/figs/exp6_adaptive_branch_<method>.png (if raw supports adaptive_branch)"
