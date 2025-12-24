#!/usr/bin/env bash
# run/run_exp2.sh
#
# EXP-2：Runtime vs t（对应 RQ2）
#
# 本脚本需与 Docs/Experiment/EXP-2.md 的表述/原理/思想严格对齐：
#   - 固定同一数据集 R,S 与固定密度参数 alpha
#   - 只 sweep t（典型：1k → 1M；具体列表由 config/sweeps/sweep_t.json 决定）
#   - 对每个 method × variant（sampling / enum_sampling / adaptive）运行 sweep
#   - 主指标：端到端 wall-clock 时间（summary 中按 method×variant×t 聚合 median/p95/stdev 等）
#   - 公平性：默认单线程（--threads=1）并默认关闭 write_samples（--write_samples=0）避免 I/O 噪声
#
# 统一目录规范（对齐你的全局要求）：
#   1) Build 统一放在 <repo_root>/build/<type>/ 下（默认 Release -> build/release）
#   2) 运行产生的日志/中间文件/CSV/图，先写到 <repo_root>/run/temp/exp2/
#   3) 成功后将 <repo_root>/run/temp/exp2/ 覆盖同步到 <repo_root>/results/raw/exp2/
#   4) Bash 内不包含任何“内嵌 Python”（不使用 python -c / heredoc）
#   5) 本脚本运行产生的 json/csv/png/pdf/log 等“必要文件”均先落到 run/temp/exp2
#
# 用法：
#   bash run/run_exp2.sh
#
# 可选参数：
#   --config <path>         Sweep JSON（默认：config/sweeps/sweep_t.json）
#   --build_type <type>     Release|Debug|RelWithDebInfo|MinSizeRel（默认：Release）
#   --clean                 删除 build/<type> 后重建
#   --no-build              跳过构建（假设 sjs_sweep 已存在于 build/<type>）
#   --no-plot               跳过绘图（仍生成 sweep_raw.csv / sweep_summary.csv）
#   --threads <int>         覆盖 sweep.json 的 sys.threads（默认：1）
#   --write_samples <0|1>   覆盖 sweep.json 的 run.write_samples（默认：0）
#
# 注意：
#   - 本实验的最终结果目录固定为 results/raw/exp2（成功后覆盖旧结果），不提供 --out_dir。
#   - 若要调整数据/alpha/t 列表/方法集合/重复次数等，请编辑 sweep_t.json。

set -euo pipefail
IFS=$'\n\t'

trap 'echo -e "[EXP2][FATAL] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2; exit 1' ERR

# ----------------------------
# Helpers
# ----------------------------
log() { echo -e "[EXP2] $*"; }
die() { echo -e "[EXP2][FATAL] $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"
}

usage() {
  cat <<'EOF'
EXP-2 runner (Runtime vs t)

Usage:
  bash run/run_exp2.sh [options]

Options:
  --config <path>         Sweep JSON (default: config/sweeps/sweep_t.json)
  --build_type <type>     Release|Debug|RelWithDebInfo|MinSizeRel (default: Release)
  --clean                 Remove build/<type> before building
  --no-build              Skip build step
  --no-plot               Skip plotting step
  --threads <int>         Override sys.threads (default: 1)
  --write_samples <0|1>   Override run.write_samples (default: 0)
  -h, --help              Show help

Outputs (fixed):
  - Temp (always overwritten):  <repo_root>/run/temp/exp2/
  - Final (overwritten on success): <repo_root>/results/raw/exp2/
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

build_subdir_from_type() {
  case "$1" in
    Release) echo "release";;
    Debug) echo "debug";;
    RelWithDebInfo) echo "relwithdebinfo";;
    MinSizeRel) echo "minsizerel";;
    *) die "Unknown --build_type '$1' (expected Release|Debug|RelWithDebInfo|MinSizeRel)";;
  esac
}

find_exe() {
  local build_dir="$1"
  local name="$2"

  # Fast paths.
  if [[ -x "${build_dir}/${name}" ]]; then
    echo "${build_dir}/${name}"
    return
  fi
  if [[ -x "${build_dir}/apps/${name}" ]]; then
    echo "${build_dir}/apps/${name}"
    return
  fi

  # Fallback: search.
  local p
  p="$(find "${build_dir}" -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "${p}" ]] || die "Cannot find executable '${name}' under ${build_dir}. Did the build succeed?"
  echo "${p}"
}

# ----------------------------
# Resolve paths (relative to this script)
# ----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ----------------------------
# Defaults (aligned with EXP-2.md)
# ----------------------------
CONFIG="config/sweeps/sweep_t.json"
BUILD_TYPE="Release"
THREADS=1
WRITE_SAMPLES=0

DO_CLEAN=0
DO_BUILD=1
DO_PLOT=1

# ----------------------------
# Parse args
# ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)        CONFIG="$2"; shift 2;;
    --build_type)    BUILD_TYPE="$2"; shift 2;;
    --threads)       THREADS="$2"; shift 2;;
    --write_samples) WRITE_SAMPLES="$2"; shift 2;;
    --clean)         DO_CLEAN=1; shift;;
    --no-build)      DO_BUILD=0; shift;;
    --no-plot)       DO_PLOT=0; shift;;
    -h|--help)       usage; exit 0;;
    *) die "Unknown argument: $1 (try --help)";;
  esac
done

# Resolve CONFIG relative to repo root (no embedded Python).
if [[ "${CONFIG}" != /* ]]; then
  CONFIG="${REPO_ROOT}/${CONFIG}"
fi

# Fixed directories (per your global requirements).
TEMP_DIR="${REPO_ROOT}/run/temp/exp2"
RESULT_DIR="${REPO_ROOT}/results/raw/exp2"

BUILD_SUBDIR="$(build_subdir_from_type "${BUILD_TYPE}")"
BUILD_DIR="${REPO_ROOT}/build/${BUILD_SUBDIR}"

# ----------------------------
# Preflight
# ----------------------------
need_cmd cmake
need_cmd tee
need_cmd find

[[ -f "${CONFIG}" ]] || die "Sweep config not found: ${CONFIG}"
[[ -f "${REPO_ROOT}/CMakeLists.txt" ]] || die "CMakeLists.txt not found at repo root: ${REPO_ROOT}"

if ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || [[ "${THREADS}" -le 0 ]]; then
  die "--threads must be a positive integer (got: '${THREADS}')"
fi
if [[ "${WRITE_SAMPLES}" != "0" && "${WRITE_SAMPLES}" != "1" ]]; then
  die "--write_samples must be 0 or 1 (got: '${WRITE_SAMPLES}')"
fi

# Clean temp (always overwrite temp).
rm -rf "${TEMP_DIR}"
mkdir -p "${TEMP_DIR}/logs"

TS="$(date +%Y%m%d_%H%M%S)"

log "Repo root      : ${REPO_ROOT}"
log "Sweep config   : ${CONFIG}"
log "Build type     : ${BUILD_TYPE}"
log "Build dir      : ${BUILD_DIR}"
log "Threads        : ${THREADS}"
log "write_samples  : ${WRITE_SAMPLES}"
log "Temp dir       : ${TEMP_DIR}"
log "Final results  : ${RESULT_DIR} (will be overwritten on success)"

# Force single-thread behavior for common libs (fairness + reproducibility).
export OMP_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export VECLIB_MAXIMUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"

# Keep an immutable copy of the sweep config for provenance.
cp -f "${CONFIG}" "${TEMP_DIR}/sweep_t_used.json"

# Save environment + manifest.
{
  echo "EXP-2 manifest"
  echo "timestamp=${TS}"
  echo "repo_root=${REPO_ROOT}"
  echo "config=${CONFIG}"
  echo "build_type=${BUILD_TYPE}"
  echo "build_dir=${BUILD_DIR}"
  echo "threads=${THREADS}"
  echo "write_samples=${WRITE_SAMPLES}"
  echo "temp_dir=${TEMP_DIR}"
  echo "result_dir=${RESULT_DIR}"
} > "${TEMP_DIR}/MANIFEST.txt"

{
  echo "timestamp=${TS}"
  echo "repo_root=${REPO_ROOT}"
  echo "config=${CONFIG}"
  echo "build_dir=${BUILD_DIR}"
  echo "temp_dir=${TEMP_DIR}"
  echo "threads=${THREADS}"
  echo "write_samples=${WRITE_SAMPLES}"
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
} > "${TEMP_DIR}/logs/env.txt"

# ----------------------------
# Build (Release by default)
# ----------------------------
if [[ "${DO_BUILD}" -eq 1 ]]; then
  if [[ "${DO_CLEAN}" -eq 1 ]]; then
    log "Cleaning build dir: ${BUILD_DIR}"
    rm -rf "${BUILD_DIR}"
  fi

  mkdir -p "${BUILD_DIR}"

  log "Configuring (${BUILD_TYPE})..."
  cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
    -DSJS_BUILD_ROOT_APPS=ON \
    -DSJS_BUILD_TESTS=OFF \
    2>&1 | tee "${TEMP_DIR}/logs/cmake_configure.log"

  log "Building..."
  cmake --build "${BUILD_DIR}" -j "$(nproc_safe)" \
    2>&1 | tee "${TEMP_DIR}/logs/cmake_build.log"
else
  log "Skipping build step (--no-build)."
fi

SJS_SWEEP="$(find_exe "${BUILD_DIR}" "sjs_sweep")"
log "Using sjs_sweep: ${SJS_SWEEP}"

# ----------------------------
# Run EXP-2 sweep
# ----------------------------
OUT_DIR="${TEMP_DIR}"  # sweep outputs (raw/summary + plots) land directly in temp
RAW_FILE="sweep_raw.csv"
SUMMARY_FILE="sweep_summary.csv"

log "Running sweep (methods × variants × t × repeats)..."
log "Command: ${SJS_SWEEP} --config ${CONFIG} --out_dir ${OUT_DIR} --threads=${THREADS} --write_samples=${WRITE_SAMPLES}"

# Run from repo root so relative paths in JSON resolve as expected.
pushd "${REPO_ROOT}" >/dev/null
set +e
"${SJS_SWEEP}" \
  --config="${CONFIG}" \
  --out_dir="${OUT_DIR}" \
  --raw_file="${RAW_FILE}" \
  --summary_file="${SUMMARY_FILE}" \
  --threads="${THREADS}" \
  --write_samples="${WRITE_SAMPLES}" \
  2>&1 | tee "${TEMP_DIR}/logs/sjs_sweep.log"
rc="${PIPESTATUS[0]}"
set -e
popd >/dev/null

if [[ "${rc}" -ne 0 ]]; then
  die "sjs_sweep failed with exit code ${rc}. See ${TEMP_DIR}/logs/sjs_sweep.log"
fi

# Basic output checks (must exist per EXP-2.md)
[[ -f "${OUT_DIR}/${RAW_FILE}" ]] || die "Missing expected output: ${OUT_DIR}/${RAW_FILE}"
[[ -f "${OUT_DIR}/${SUMMARY_FILE}" ]] || die "Missing expected output: ${OUT_DIR}/${SUMMARY_FILE}"

log "Sweep finished OK."
log "Raw     : ${OUT_DIR}/${RAW_FILE}"
log "Summary : ${OUT_DIR}/${SUMMARY_FILE}"

# ----------------------------
# (Optional) Plot
# ----------------------------
if [[ "${DO_PLOT}" -eq 1 ]]; then
  if command -v python3 >/dev/null 2>&1; then
    PLOT_SCRIPT="${REPO_ROOT}/run/include/exp2_plot.py"
    [[ -f "${PLOT_SCRIPT}" ]] || die "Missing plot script: ${PLOT_SCRIPT}"

    log "Generating plots via: ${PLOT_SCRIPT}"
    # Note: plots are written into OUT_DIR (= run/temp/exp2)
    python3 "${PLOT_SCRIPT}" --out_dir "${OUT_DIR}" --t0 1000 --error p95 \
      2>&1 | tee "${TEMP_DIR}/logs/plot_exp2.log"
  else
    log "python3 not found; skipping plot step."
  fi
else
  log "Plot step disabled (--no-plot)."
fi

# ----------------------------
# Sync temp -> results/raw/exp2 (overwrite)
# ----------------------------
log "Syncing results to: ${RESULT_DIR} (overwrite)"
rm -rf "${RESULT_DIR}"
mkdir -p "${RESULT_DIR}"
cp -a "${TEMP_DIR}/." "${RESULT_DIR}/"

log "DONE."
log "Temp dir     : ${TEMP_DIR}"
log "Final results: ${RESULT_DIR}"
