#!/usr/bin/env bash
# run/run_exp3.sh
# ------------------------------------------------------------------------------
# EXP-3: Runtime vs 密度参数 α (RQ3)
#
# 目标（严格对齐 EXP-3.md）：
#   - 扫描 alpha = |J|/(n_r+n_s)（synthetic stripe/可控密度条带构造）
#   - 记录端到端 wall time（median + p95）与 ok_rate
#   - 从 sweep_raw.csv 导出 adaptive 分支比例（enumerate_all vs fallback）
#   - （可选）生成 runtime-vs-alpha 与 branch-ratio 图（symlog-x, log-y）
#
# 目录规范（对齐你提出的 5 条要求）：
#   - build:   <repo_root>/build/<build_type_subdir>
#   - temp:    <repo_root>/run/temp/exp3        （本脚本所有产物先写这里）
#   - results: <repo_root>/results/raw/exp3     （本脚本结束时覆盖同步到这里）
#
# 依赖：cmake, python3 （以及可选 matplotlib 用于画图）
#
# 用法：
#   bash run/run_exp3.sh
#
# 可选环境变量（保持与旧脚本一致 + 增补 build 规范）：
#   EXP3_CONFIG            : sweep JSON 路径 (默认: config/sweeps/sweep_alpha.json)
#   EXP3_BUILD_TYPE        : Release|Debug|RelWithDebInfo|MinSizeRel (默认: Release)
#   EXP3_CLEAN_BUILD=1     : 清理 build/<type> 后重建 (默认: 0)
#   EXP3_JOBS              : 编译并行度 (默认: nproc)
#   EXP3_RUN_TESTS=1       : 编译后跑 ctest (默认: 0)
#   EXP3_PLOT=0            : 跳过画图 (默认: 1)
#   EXP3_PLOT_ENUM=1       : 画图时包含 enum_sampling 曲线 (默认: 0)
#
#   # 不改 JSON 的参数覆盖（由 run/include/exp3_make_effective_sweep_config.py 读取）：
#   EXP3_NR, EXP3_NS
#   EXP3_T
#   EXP3_REPEATS
#   EXP3_JSTAR
#   EXP3_ENUM_CAP
#   EXP3_THREADS           : sys.threads（默认强制 1，保证公平；可显式设为 >1）
#   EXP3_ALPHA_LIST        : 逗号分隔 alpha 列表（float）
#   EXP3_METHODS           : 逗号分隔 method 列表
#   EXP3_VARIANTS          : 逗号分隔 variant 列表
# ------------------------------------------------------------------------------

set -euo pipefail
IFS=$'\n\t'

trap 'echo "[EXP3][FATAL] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

log()  { echo "[EXP3] $*"; }
die()  { echo "[EXP3][ERROR] $*" >&2; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

nproc_safe() {
  if [[ -n "${EXP3_JOBS:-}" ]]; then
    echo "${EXP3_JOBS}"; return
  fi
  if command -v nproc >/dev/null 2>&1; then
    nproc; return
  fi
  if [[ "$(uname -s)" == "Darwin" ]] && command -v sysctl >/dev/null 2>&1; then
    sysctl -n hw.ncpu 2>/dev/null || echo 4
    return
  fi
  echo 4
}

build_subdir_from_type() {
  case "$1" in
    Release) echo "release";;
    Debug) echo "debug";;
    RelWithDebInfo) echo "relwithdebinfo";;
    MinSizeRel) echo "minsizerel";;
    *) die "Unsupported EXP3_BUILD_TYPE='$1' (use Release|Debug|RelWithDebInfo|MinSizeRel)";;
  esac
}

find_exe() {
  local build_dir="$1"
  local name="$2"
  local p=""

  # Common direct locations
  if [[ -x "${build_dir}/${name}" ]]; then
    echo "${build_dir}/${name}"; return
  fi
  if [[ -x "${build_dir}/apps/${name}" ]]; then
    echo "${build_dir}/apps/${name}"; return
  fi

  p="$(find "${build_dir}" -maxdepth 4 -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "${p}" ]] || die "Cannot find executable '${name}' under ${build_dir}. Did the build succeed?"
  echo "${p}"
}

# -----------------------------
# Resolve repo root
# -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# -----------------------------
# Paths (fixed by repo policy)
# -----------------------------
BUILD_ROOT="${ROOT_DIR}/build"
TEMP_DIR="${ROOT_DIR}/run/temp/exp3"
RESULT_DIR="${ROOT_DIR}/results/raw/exp3"
INCLUDE_DIR="${ROOT_DIR}/run/include"

# -----------------------------
# Inputs / switches
# -----------------------------
DEFAULT_CONFIG="${ROOT_DIR}/config/sweeps/sweep_alpha.json"
CONFIG_PATH="${EXP3_CONFIG:-${DEFAULT_CONFIG}}"
if [[ "${CONFIG_PATH}" != /* ]]; then
  CONFIG_PATH="${ROOT_DIR}/${CONFIG_PATH}"
fi
[[ -f "${CONFIG_PATH}" ]] || die "Sweep config not found: ${CONFIG_PATH}"

BUILD_TYPE="${EXP3_BUILD_TYPE:-Release}"
BUILD_SUBDIR="$(build_subdir_from_type "${BUILD_TYPE}")"
BUILD_DIR="${BUILD_ROOT}/${BUILD_SUBDIR}"

DO_CLEAN_BUILD="${EXP3_CLEAN_BUILD:-0}"
RUN_TESTS="${EXP3_RUN_TESTS:-0}"
DO_PLOT="${EXP3_PLOT:-1}"
PLOT_ENUM="${EXP3_PLOT_ENUM:-0}"

JOBS="$(nproc_safe)"
TS="$(date +%Y%m%d_%H%M%S)"

# Default fairness: single-thread unless user explicitly overrides.
export EXP3_THREADS="${EXP3_THREADS:-1}"

# -----------------------------
# Dependency checks
# -----------------------------
need_cmd cmake
need_cmd python3
need_cmd tee
need_cmd find

# -----------------------------
# Init temp dir (all produced files go here first)
# -----------------------------
rm -rf "${TEMP_DIR}"
mkdir -p "${TEMP_DIR}/logs" "${TEMP_DIR}/derived" "${TEMP_DIR}/plots" "${TEMP_DIR}/meta"

log "Repo root   : ${ROOT_DIR}"
log "Config      : ${CONFIG_PATH}"
log "Build dir   : ${BUILD_DIR} (type=${BUILD_TYPE})"
log "Temp dir    : ${TEMP_DIR}"
log "Result dir  : ${RESULT_DIR} (will be overwritten on success)"
log "Jobs        : ${JOBS}"
log "sys.threads : ${EXP3_THREADS} (default=1 for fairness; override by setting EXP3_THREADS)"

# -----------------------------
# Repro / fairness knobs
# -----------------------------
# Cap common BLAS/OpenMP thread pools for fairness.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# -----------------------------
# Manifest (what we ran)
# -----------------------------
{
  echo "exp=exp3"
  echo "timestamp=${TS}"
  echo "repo_root=${ROOT_DIR}"
  echo "config=${CONFIG_PATH}"
  echo "build_type=${BUILD_TYPE}"
  echo "build_dir=${BUILD_DIR}"
  echo "jobs=${JOBS}"
  echo "sys_threads=${EXP3_THREADS}"
  echo "run_tests=${RUN_TESTS}"
  echo "plot=${DO_PLOT}"
  echo "plot_enum_sampling=${PLOT_ENUM}"
  echo ""
  echo "uname:"; uname -a || true
  echo ""
  echo "cmake:"; cmake --version 2>/dev/null || true
  echo ""
  if command -v git >/dev/null 2>&1 && [[ -d "${ROOT_DIR}/.git" ]]; then
    echo "git:";
    (cd "${ROOT_DIR}" && git rev-parse HEAD && git status --porcelain) || true
  fi
  echo ""
  echo "EXP3 overrides (env, optional):"
  echo "  EXP3_NR=${EXP3_NR:-}"
  echo "  EXP3_NS=${EXP3_NS:-}"
  echo "  EXP3_T=${EXP3_T:-}"
  echo "  EXP3_REPEATS=${EXP3_REPEATS:-}"
  echo "  EXP3_JSTAR=${EXP3_JSTAR:-}"
  echo "  EXP3_ENUM_CAP=${EXP3_ENUM_CAP:-}"
  echo "  EXP3_ALPHA_LIST=${EXP3_ALPHA_LIST:-}"
  echo "  EXP3_METHODS=${EXP3_METHODS:-}"
  echo "  EXP3_VARIANTS=${EXP3_VARIANTS:-}"
} > "${TEMP_DIR}/meta/manifest.txt"

# -----------------------------
# Build
# -----------------------------
log "Building project (${BUILD_TYPE}) ..."
mkdir -p "${BUILD_ROOT}"

if [[ "${DO_CLEAN_BUILD}" == "1" ]]; then
  log "EXP3_CLEAN_BUILD=1 -> rm -rf ${BUILD_DIR}"
  rm -rf "${BUILD_DIR}"
fi

# Configure
cmake_args=(
  -S "${ROOT_DIR}"
  -B "${BUILD_DIR}"
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
  -DSJS_BUILD_ROOT_APPS=ON
)
if [[ "${RUN_TESTS}" == "1" ]]; then
  cmake_args+=( -DSJS_BUILD_TESTS=ON )
fi

cmake "${cmake_args[@]}" 2>&1 | tee "${TEMP_DIR}/logs/cmake_configure.log"
cmake --build "${BUILD_DIR}" -j "${JOBS}" 2>&1 | tee "${TEMP_DIR}/logs/cmake_build.log"

if [[ "${RUN_TESTS}" == "1" ]]; then
  need_cmd ctest
  log "Running tests (ctest) ..."
  (cd "${BUILD_DIR}" && ctest --output-on-failure) 2>&1 | tee "${TEMP_DIR}/logs/ctest.log"
fi

SJS_SWEEP="$(find_exe "${BUILD_DIR}" "sjs_sweep")"
log "Using sjs_sweep: ${SJS_SWEEP}"

# -----------------------------
# Prepare effective sweep config (JSON) under run/temp/exp3
# -----------------------------
EFFECTIVE_CONFIG="${TEMP_DIR}/sweep_exp3_effective.json"
[[ -f "${INCLUDE_DIR}/exp3_make_effective_sweep_config.py" ]] || die "Missing: ${INCLUDE_DIR}/exp3_make_effective_sweep_config.py"

python3 "${INCLUDE_DIR}/exp3_make_effective_sweep_config.py" \
  --in "${CONFIG_PATH}" \
  --out "${EFFECTIVE_CONFIG}" \
  --out_dir "${TEMP_DIR}" \
  2>&1 | tee "${TEMP_DIR}/logs/make_effective_config.log"

# -----------------------------
# Run sweep
# -----------------------------
log "Running sweep (alpha scan) ..."
SWEEP_LOG="${TEMP_DIR}/logs/exp3_sweep.log"

pushd "${ROOT_DIR}" >/dev/null
set +e
"${SJS_SWEEP}" --config="${EFFECTIVE_CONFIG}" 2>&1 | tee "${SWEEP_LOG}"
rc="${PIPESTATUS[0]}"
set -e
popd >/dev/null

if [[ "${rc}" -ne 0 ]]; then
  die "sjs_sweep failed with exit code ${rc}. See ${SWEEP_LOG}"
fi

RAW_CSV="${TEMP_DIR}/sweep_raw.csv"
SUMMARY_CSV="${TEMP_DIR}/sweep_summary.csv"

# Some sweep implementations may nest outputs; fall back to a quick search.
if [[ ! -f "${RAW_CSV}" ]]; then
  RAW_CSV_FOUND="$(find "${TEMP_DIR}" -maxdepth 3 -name sweep_raw.csv | head -n 1 || true)"
  [[ -n "${RAW_CSV_FOUND}" ]] && RAW_CSV="${RAW_CSV_FOUND}"
fi
if [[ ! -f "${SUMMARY_CSV}" ]]; then
  SUMMARY_CSV_FOUND="$(find "${TEMP_DIR}" -maxdepth 3 -name sweep_summary.csv | head -n 1 || true)"
  [[ -n "${SUMMARY_CSV_FOUND}" ]] && SUMMARY_CSV="${SUMMARY_CSV_FOUND}"
fi

[[ -f "${RAW_CSV}" ]] || die "Missing expected output: sweep_raw.csv under ${TEMP_DIR}"
[[ -f "${SUMMARY_CSV}" ]] || die "Missing expected output: sweep_summary.csv under ${TEMP_DIR}"

log "Sweep done."
log "  Raw     : ${RAW_CSV}"
log "  Summary : ${SUMMARY_CSV}"

# -----------------------------
# Derive: adaptive branch ratio (from raw)
# -----------------------------
BRANCH_CSV="${TEMP_DIR}/derived/adaptive_branch_ratio.csv"
[[ -f "${INCLUDE_DIR}/exp3_derive_adaptive_branch_ratio.py" ]] || die "Missing: ${INCLUDE_DIR}/exp3_derive_adaptive_branch_ratio.py"

python3 "${INCLUDE_DIR}/exp3_derive_adaptive_branch_ratio.py" \
  --raw "${RAW_CSV}" \
  --out "${BRANCH_CSV}" \
  2>&1 | tee "${TEMP_DIR}/logs/derive_branch_ratio.log"

# -----------------------------
# Optional: plotting (matplotlib)
# -----------------------------
if [[ "${DO_PLOT}" == "1" ]]; then
  [[ -f "${INCLUDE_DIR}/exp3_plot.py" ]] || die "Missing: ${INCLUDE_DIR}/exp3_plot.py"

  plot_args=(
    --summary "${SUMMARY_CSV}"
    --branch "${BRANCH_CSV}"
    --out_dir "${TEMP_DIR}/plots"
  )
  if [[ "${PLOT_ENUM}" == "1" ]]; then
    plot_args+=( --include_enum_sampling )
  fi

  python3 "${INCLUDE_DIR}/exp3_plot.py" "${plot_args[@]}" \
    2>&1 | tee "${TEMP_DIR}/logs/plot.log"
else
  log "EXP3_PLOT=0 -> skip plotting"
fi

# -----------------------------
# Sync to results/raw/exp3 (overwrite old results)
# -----------------------------
log "Syncing results to ${RESULT_DIR} (overwrite) ..."
rm -rf "${RESULT_DIR}"
mkdir -p "${RESULT_DIR}"
cp -a "${TEMP_DIR}/." "${RESULT_DIR}/"

log "DONE ✅"
log "----------------------------------------"
log "EXP-3 final results directory:"
log "  ${RESULT_DIR}"
log "Key files:"
log "  ${RESULT_DIR}/sweep_raw.csv"
log "  ${RESULT_DIR}/sweep_summary.csv"
log "  ${RESULT_DIR}/derived/adaptive_branch_ratio.csv"
log "  ${RESULT_DIR}/meta/manifest.txt"
log "----------------------------------------"
