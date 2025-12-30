#!/usr/bin/env bash
# run/run_exp3.sh
# ------------------------------------------------------------------------------
# EXP-3: Runtime vs 密度参数 α (RQ3)
#
# ✅ 这份脚本只做“实验叙事增强”，不做任何不公平的比较：
#   - 不改 C++ 实现、不改计时口径、不删 baseline
#   - 仅增强 alpha 点位覆盖（尤其是 alpha* 附近）与输出结构（更 paper-friendly）
#
# 目标（严格对齐 EXP-3.md）：
#   - 扫描 alpha = |J|/(n_r+n_s)（synthetic stripe/可控密度条带构造）
#   - 记录端到端 wall time（median + p95）与 ok_rate
#   - 从 sweep_raw.csv 导出 adaptive 分支比例（enumerate_all vs fallback）
#   - （可选）生成 runtime-vs-alpha 与 branch-ratio 图（symlog-x, log-y）
#
# 目录规范：
#   - build:   <repo_root>/build/<build_type_subdir>
#   - temp:    <repo_root>/run/temp/exp3
#   - results: <repo_root>/results/raw/exp3
#
# 用法：
#   bash run/run_exp3.sh
#
# 可选环境变量（保持旧脚本一致 + 增补）：
#   EXP3_CONFIG            : sweep JSON 路径 (默认: config/sweeps/sweep_alpha.json)
#   EXP3_BUILD_TYPE        : Release|Debug|RelWithDebInfo|MinSizeRel (默认: Release)
#   EXP3_CLEAN_BUILD=1     : 清理 build/<type> 后重建 (默认: 0)
#   EXP3_JOBS              : 编译并行度 (默认: nproc)
#   EXP3_RUN_TESTS=1       : 编译后跑 ctest (默认: 0)
#
#   EXP3_PLOT=0            : 跳过画图 (默认: 1)
#   EXP3_PLOT_ENUM=1       : 额外生成包含 enum_sampling 的图（默认: 0）
#
#   # 不改 JSON 的参数覆盖（由 run/include/exp3_make_effective_sweep_config.py 读取）：
#   EXP3_NR, EXP3_NS
#   EXP3_T
#   EXP3_REPEATS
#   EXP3_JSTAR
#   EXP3_ENUM_CAP
#   EXP3_THREADS           : sys.threads（默认=1；可显式设为 >1）
#   EXP3_ALPHA_LIST        : 逗号分隔 alpha 列表（float）
#   EXP3_METHODS           : 逗号分隔 method 列表
#   EXP3_VARIANTS          : 逗号分隔 variant 列表
#
#   # 新增（仅影响脚本行为，不影响算法公平性）：
#   EXP3_AUTO_ALPHA_LIST=1 : 若未显式指定 EXP3_ALPHA_LIST，则自动设置更“paper-friendly”的默认列表
#   EXP3_KEEP_HISTORY=1    : 结果不覆盖写入 results/raw/exp3/<timestamp>/ 并更新 latest 软链
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
# Paths
# -----------------------------
BUILD_ROOT="${ROOT_DIR}/build"
TEMP_DIR="${ROOT_DIR}/run/temp/exp3"
RESULT_ROOT="${ROOT_DIR}/results/raw/exp3"
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

AUTO_ALPHA="${EXP3_AUTO_ALPHA_LIST:-1}"
KEEP_HISTORY="${EXP3_KEEP_HISTORY:-0}"

JOBS="$(nproc_safe)"
TS="$(date +%Y%m%d_%H%M%S)"

# Default fairness: single-thread unless user explicitly overrides.
export EXP3_THREADS="${EXP3_THREADS:-1}"

# -----------------------------
# If user did NOT specify alpha list, use a more paper-friendly default
# (covers alpha* vicinity; EXP-3.md also recommends adding 4/5/6).
# -----------------------------
if [[ -z "${EXP3_ALPHA_LIST:-}" && "${AUTO_ALPHA}" == "1" ]]; then
  # If user provided NR/NS/JSTAR, we can compute a rough alpha* for logging & optional densify.
  # Otherwise we fall back to the recommended list that includes 4/5/6.
  ALPHA_STAR=""
  if [[ -n "${EXP3_NR:-}" && -n "${EXP3_NS:-}" && -n "${EXP3_JSTAR:-}" ]]; then
    export _EXP3_NR="${EXP3_NR}"
    export _EXP3_NS="${EXP3_NS}"
    export _EXP3_JSTAR="${EXP3_JSTAR}"
    ALPHA_STAR="$(python3 - <<'PY'
import os, math
nr=float(os.environ["_EXP3_NR"])
ns=float(os.environ["_EXP3_NS"])
j=float(os.environ["_EXP3_JSTAR"])
if nr+ns <= 0:
    print("")
else:
    print(j/(nr+ns))
PY
)"
  fi

  if [[ -n "${ALPHA_STAR}" ]]; then
    # Build: base + dense around floor(alpha*) ±2 + tail
    export _EXP3_ALPHA_STAR="${ALPHA_STAR}"
    export EXP3_ALPHA_LIST="$(python3 - <<'PY'
import os, math
a=float(os.environ["_EXP3_ALPHA_STAR"])
base=[0,0.03,0.1,0.3,1,3]
s=max(0, int(math.floor(a))-2)
e=max(0, int(math.floor(a))+2)
dense=[float(i) for i in range(s,e+1)]
tail=[10,30,100,300]
seen=set()
out=[]
for v in base+dense+tail:
    if v not in seen:
        seen.add(v); out.append(v)
# nice formatting
def fmt(x):
    if abs(x-round(x))<1e-12:
        return str(int(round(x)))
    return ("%g" % x)
print(",".join(fmt(x) for x in out))
PY
)"
    log "AUTO alpha list enabled. Computed alpha*≈${ALPHA_STAR}; using EXP3_ALPHA_LIST=${EXP3_ALPHA_LIST}"
  else
    export EXP3_ALPHA_LIST="0,0.03,0.1,0.3,1,3,4,5,6,10,30,100,300"
    log "AUTO alpha list enabled. Using EXP3_ALPHA_LIST=${EXP3_ALPHA_LIST}"
  fi
elif [[ -n "${EXP3_ALPHA_LIST:-}" ]]; then
  log "Using user-provided EXP3_ALPHA_LIST=${EXP3_ALPHA_LIST}"
else
  log "AUTO alpha list disabled (EXP3_AUTO_ALPHA_LIST=${AUTO_ALPHA}). Using JSON alpha list."
fi

# -----------------------------
# Dependency checks
# -----------------------------
need_cmd cmake
need_cmd python3
need_cmd tee
need_cmd find

# -----------------------------
# Init temp dir
# -----------------------------
rm -rf "${TEMP_DIR}"
mkdir -p "${TEMP_DIR}/logs" "${TEMP_DIR}/derived" "${TEMP_DIR}/plots/main" "${TEMP_DIR}/plots/with_enum" "${TEMP_DIR}/meta"

log "Repo root   : ${ROOT_DIR}"
log "Config      : ${CONFIG_PATH}"
log "Build dir   : ${BUILD_DIR} (type=${BUILD_TYPE})"
log "Temp dir    : ${TEMP_DIR}"
log "Result root : ${RESULT_ROOT}"
log "Jobs        : ${JOBS}"
log "sys.threads : ${EXP3_THREADS} (default=1 for fairness; override by setting EXP3_THREADS)"
log "Keep history: ${KEEP_HISTORY}"
log "Plot        : ${DO_PLOT} (extra enum plots: ${PLOT_ENUM})"

# -----------------------------
# Repro / fairness knobs
# -----------------------------
# Cap common BLAS/OpenMP thread pools. If you run multi-thread (EXP3_THREADS>1),
# we follow that cap to avoid "asked for >1 but got forced to 1".
THREAD_CAP="${EXP3_THREADS}"
if ! [[ "${THREAD_CAP}" =~ ^[0-9]+$ ]]; then THREAD_CAP=1; fi
export OMP_NUM_THREADS="${THREAD_CAP}"
export OPENBLAS_NUM_THREADS="${THREAD_CAP}"
export MKL_NUM_THREADS="${THREAD_CAP}"
export VECLIB_MAXIMUM_THREADS="${THREAD_CAP}"
export NUMEXPR_NUM_THREADS="${THREAD_CAP}"

# -----------------------------
# Manifest
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
  echo "thread_cap=${THREAD_CAP}"
  echo "run_tests=${RUN_TESTS}"
  echo "plot=${DO_PLOT}"
  echo "plot_enum_extra=${PLOT_ENUM}"
  echo "keep_history=${KEEP_HISTORY}"
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
# Prepare effective sweep config
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
# Derive: adaptive branch ratio
# -----------------------------
BRANCH_CSV="${TEMP_DIR}/derived/adaptive_branch_ratio.csv"
[[ -f "${INCLUDE_DIR}/exp3_derive_adaptive_branch_ratio.py" ]] || die "Missing: ${INCLUDE_DIR}/exp3_derive_adaptive_branch_ratio.py"

python3 "${INCLUDE_DIR}/exp3_derive_adaptive_branch_ratio.py" \
  --raw "${RAW_CSV}" \
  --out "${BRANCH_CSV}" \
  2>&1 | tee "${TEMP_DIR}/logs/derive_branch_ratio.log"

# -----------------------------
# Plotting
# -----------------------------
if [[ "${DO_PLOT}" == "1" ]]; then
  [[ -f "${INCLUDE_DIR}/exp3_plot.py" ]] || die "Missing: ${INCLUDE_DIR}/exp3_plot.py"

  # Main plots (recommended for paper): default DOES NOT include enum_sampling
  log "Plotting (main: no enum_sampling) ..."
  python3 "${INCLUDE_DIR}/exp3_plot.py" \
    --summary "${SUMMARY_CSV}" \
    --branch "${BRANCH_CSV}" \
    --out_dir "${TEMP_DIR}/plots/main" \
    2>&1 | tee "${TEMP_DIR}/logs/plot_main.log"

  # Optional appendix plots: include enum_sampling (if requested)
  if [[ "${PLOT_ENUM}" == "1" ]]; then
    log "Plotting (appendix: include enum_sampling) ..."
    python3 "${INCLUDE_DIR}/exp3_plot.py" \
      --summary "${SUMMARY_CSV}" \
      --branch "${BRANCH_CSV}" \
      --out_dir "${TEMP_DIR}/plots/with_enum" \
      --include_enum_sampling \
      2>&1 | tee "${TEMP_DIR}/logs/plot_with_enum.log"
  fi
else
  log "EXP3_PLOT=0 -> skip plotting"
fi

# -----------------------------
# Sync results
# -----------------------------
if [[ "${KEEP_HISTORY}" == "1" ]]; then
  FINAL_DIR="${RESULT_ROOT}/${TS}"
  log "Syncing results to ${FINAL_DIR} (keep history) ..."
  mkdir -p "${FINAL_DIR}"
  cp -a "${TEMP_DIR}/." "${FINAL_DIR}/"
  # update latest link (best effort)
  ln -sfn "${FINAL_DIR}" "${RESULT_ROOT}/latest" 2>/dev/null || true
else
  FINAL_DIR="${RESULT_ROOT}"
  log "Syncing results to ${FINAL_DIR} (overwrite) ..."
  rm -rf "${FINAL_DIR}"
  mkdir -p "${FINAL_DIR}"
  cp -a "${TEMP_DIR}/." "${FINAL_DIR}/"
fi

log "DONE ✅"
log "----------------------------------------"
log "EXP-3 final results directory:"
log "  ${FINAL_DIR}"
log "Key files:"
log "  ${FINAL_DIR}/sweep_raw.csv"
log "  ${FINAL_DIR}/sweep_summary.csv"
log "  ${FINAL_DIR}/derived/adaptive_branch_ratio.csv"
log "  ${FINAL_DIR}/meta/manifest.txt"
log "Plots:"
log "  ${FINAL_DIR}/plots/main"
if [[ "${PLOT_ENUM}" == "1" ]]; then
  log "  ${FINAL_DIR}/plots/with_enum"
fi
log "----------------------------------------"
