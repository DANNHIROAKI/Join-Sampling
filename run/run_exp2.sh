#!/usr/bin/env bash
# run/run_exp2.sh
#
# EXP-2: Runtime vs t (RQ2)
# - Fixed dataset + alpha, sweep t from 1k -> 1M
# - Runs ALL methods x variants listed in config/sweeps/sweep_t.json
# - Writes results to a dedicated output directory under results/
#
# This version calls run/plot_exp2.py to generate:
#   1) runtime vs t (wall median, log-x log-y, upper-only p95)
#   2) Δruntime vs t (baseline-subtracted, log-x linear-y)
#   3) sample-phase vs t (median run_sample_ms from phases_json, log-x log-y)
#   4) ns/sample slopes (CSV)
#
# Usage:
#   bash run/run_exp2.sh
#
# Optional overrides:
#   bash run/run_exp2.sh --out_dir results/exp2_runtime_vs_t
#   bash run/run_exp2.sh --config config/sweeps/sweep_t.json
#   bash run/run_exp2.sh --clean        # remove build dir first
#   bash run/run_exp2.sh --no-plot      # skip plotting step
#   bash run/run_exp2.sh --no-build     # skip build step

set -euo pipefail

# ----------------------------
# Helpers
# ----------------------------
die() { echo "[EXP2][FATAL] $*" >&2; exit 1; }
log() { echo "[EXP2] $*"; }

usage() {
  cat <<'EOF'
EXP-2 runner (Runtime vs t)

Usage:
  bash run/run_exp2.sh [options]

Options:
  --config <path>     Sweep JSON (default: config/sweeps/sweep_t.json)
  --build_dir <path>  Build directory (default: build_exp2_release)
  --out_dir <path>    Output directory (default: results/exp2/runtime_vs_t/<timestamp>)
  --clean             Remove build_dir before building
  --no-build          Skip build step (assumes sjs_sweep already built)
  --no-plot           Skip plotting step
  -h, --help          Show help
EOF
}

# Determine CPU parallelism (portable-ish)
nproc_safe() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
  elif [[ "$(uname -s)" == "Darwin" ]] && command -v sysctl >/dev/null 2>&1; then
    sysctl -n hw.ncpu
  else
    echo 4
  fi
}

find_exe() {
  local build_dir="$1"
  local name="$2"
  local p
  p="$(find "${build_dir}" -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "${p}" ]] || die "Cannot find executable '${name}' under ${build_dir}. Did the build succeed?"
  echo "${p}"
}

# ----------------------------
# Paths (relative to this script)
# ----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DEFAULT_CONFIG="${REPO_ROOT}/config/sweeps/sweep_t.json"
DEFAULT_BUILD_DIR="${REPO_ROOT}/build_exp2_release"
TS="$(date +%Y%m%d_%H%M%S)"
DEFAULT_OUT_DIR="${REPO_ROOT}/results/exp2/runtime_vs_t/${TS}"

CONFIG="${DEFAULT_CONFIG}"
BUILD_DIR="${DEFAULT_BUILD_DIR}"
OUT_DIR="${DEFAULT_OUT_DIR}"

DO_CLEAN=0
DO_BUILD=1
DO_PLOT=1

# ----------------------------
# Parse args
# ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)    CONFIG="$2"; shift 2;;
    --build_dir) BUILD_DIR="$2"; shift 2;;
    --out_dir)   OUT_DIR="$2"; shift 2;;
    --clean)     DO_CLEAN=1; shift;;
    --no-build)  DO_BUILD=0; shift;;
    --no-plot)   DO_PLOT=0; shift;;
    -h|--help)   usage; exit 0;;
    *) die "Unknown argument: $1 (try --help)";;
  esac
done

# Make paths absolute for robustness (if python3 exists)
if command -v python3 >/dev/null 2>&1; then
  CONFIG="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${CONFIG}")"
  BUILD_DIR="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${BUILD_DIR}")"
  OUT_DIR="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${OUT_DIR}")"
  REPO_ROOT="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${REPO_ROOT}")"
fi

[[ -f "${CONFIG}" ]] || die "Sweep config not found: ${CONFIG}"
[[ -f "${REPO_ROOT}/CMakeLists.txt" ]] || die "CMakeLists.txt not found at repo root: ${REPO_ROOT}"

log "Repo root : ${REPO_ROOT}"
log "Config    : ${CONFIG}"
log "Build dir : ${BUILD_DIR}"
log "Out dir   : ${OUT_DIR}"

# Force single-thread behavior for common libs (even if unused)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# ----------------------------
# Build (Release)
# ----------------------------
if [[ "${DO_BUILD}" -eq 1 ]]; then
  if ! command -v cmake >/dev/null 2>&1; then
    die "cmake not found. Please install cmake first."
  fi

  if [[ "${DO_CLEAN}" -eq 1 ]]; then
    log "Cleaning build dir: ${BUILD_DIR}"
    rm -rf "${BUILD_DIR}"
  fi

  mkdir -p "${BUILD_DIR}"
  log "Configuring (Release)..."
  cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release

  log "Building..."
  cmake --build "${BUILD_DIR}" -j "$(nproc_safe)"
else
  log "Skipping build step (--no-build)."
fi

# Locate sjs_sweep
SJS_SWEEP="$(find_exe "${BUILD_DIR}" "sjs_sweep")"
log "Using sjs_sweep: ${SJS_SWEEP}"

# ----------------------------
# Run EXP-2 sweep
# ----------------------------
mkdir -p "${OUT_DIR}"
mkdir -p "${OUT_DIR}/logs"

# Save environment info (helps reproducibility)
{
  echo "timestamp=${TS}"
  echo "repo_root=${REPO_ROOT}"
  echo "config=${CONFIG}"
  echo "build_dir=${BUILD_DIR}"
  echo "out_dir=${OUT_DIR}"
  echo
  echo "uname:"
  uname -a || true
  echo
  echo "compiler:"
  (c++ --version || g++ --version || clang++ --version) 2>/dev/null || true
  echo
  echo "cmake:"
  cmake --version 2>/dev/null || true
  echo
  if command -v git >/dev/null 2>&1 && [[ -d "${REPO_ROOT}/.git" ]]; then
    echo "git:"
    (cd "${REPO_ROOT}" && git rev-parse HEAD && git status --porcelain) || true
  fi
} > "${OUT_DIR}/logs/env.txt"

log "Running sweep (methods x variants x t x repeats)..."
log "Command: ${SJS_SWEEP} --config ${CONFIG} --out_dir ${OUT_DIR}"

# Run from repo root so relative paths in JSON work as expected.
pushd "${REPO_ROOT}" >/dev/null
set +e
"${SJS_SWEEP}" --config "${CONFIG}" --out_dir "${OUT_DIR}" 2>&1 | tee "${OUT_DIR}/logs/sjs_sweep.log"
rc="${PIPESTATUS[0]}"
set -e
popd >/dev/null

if [[ "${rc}" -ne 0 ]]; then
  die "sjs_sweep failed with exit code ${rc}. See ${OUT_DIR}/logs/sjs_sweep.log"
fi

# Basic output checks
[[ -f "${OUT_DIR}/sweep_raw.csv" ]] || die "Missing expected output: ${OUT_DIR}/sweep_raw.csv"
[[ -f "${OUT_DIR}/sweep_summary.csv" ]] || die "Missing expected output: ${OUT_DIR}/sweep_summary.csv"

log "Sweep finished OK."
log "Raw     : ${OUT_DIR}/sweep_raw.csv"
log "Summary : ${OUT_DIR}/sweep_summary.csv"

# ----------------------------
# (Optional) Plot (enhanced)
# ----------------------------
if [[ "${DO_PLOT}" -eq 1 ]]; then
  if command -v python3 >/dev/null 2>&1; then
    PLOT_SCRIPT="${REPO_ROOT}/run/plot_exp2.py"
    [[ -f "${PLOT_SCRIPT}" ]] || die "Missing plot script: ${PLOT_SCRIPT} (did you create run/plot_exp2.py?)"

    log "Generating enhanced plots via: ${PLOT_SCRIPT}"
    # If matplotlib is missing, plot_exp2.py will skip gracefully.
    python3 "${PLOT_SCRIPT}" --out_dir "${OUT_DIR}" --t0 1000 --error p95
  else
    log "python3 not found; skipping plot step."
  fi
else
  log "Plot step disabled (--no-plot)."
fi

log "DONE. Results in: ${OUT_DIR}"
