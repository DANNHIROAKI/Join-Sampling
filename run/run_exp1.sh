#!/usr/bin/env bash
# run/run_exp1.sh
#
# EXP-1：Correctness & Sample Quality（对应 RQ1）
#
# 本脚本需与 Docs/Experiment/EXP-1.md 的表述/原理/思想严格对齐：
#   - 在“小规模、可算 oracle”的数据上，用 sjs_verify 做：
#       method × variant × seed(repeats)
#   - 输出并汇总采样质量指标：
#       (1) correctness: missing_in_universe（必须为 0）
#       (2) uniformity:  χ² / KS p-value（统计意义下不拒绝均匀）
#       (3) independence: autocorrelation（sanity check，应接近 0）
#
# 统一目录规范（对齐你的全局要求）：
#   1) Build 统一放在 <repo_root>/build/ 下
#   2) 运行过程中产生的日志/中间文件先写到 <repo_root>/run/temp/exp1/
#   3) 最终实验结果同步到 <repo_root>/results/raw/exp1/（每次覆盖旧结果）
#   4) Bash 内嵌的 Python 已剥离到 <repo_root>/run/include/
#      - 本实验使用：run/include/exp1_parse_verify_logs.py
#
# 用法：
#   bash run/run_exp1.sh
#
# 可选环境变量（覆盖默认参数）：
#   BUILD_TYPE=Release|Debug|RelWithDebInfo|MinSizeRel
#   CLEAN_BUILD=0|1            (默认 0：不清理 build/<type>)
#   RUN_TESTS=0|1              (默认 1：运行 ctest)
#   JOBS=8
#   THREADS=1                  (传给 sjs_verify 的 --threads)
#
#   # 数据设置（EXP-1 推荐小规模，可算 oracle）
#   GEN=stripe                 (SYN-CTRL 条带构造；stripe 是 stripe_ctrl_alpha 的 alias)
#   NR=2000
#   NS=2000
#   ALPHA=1                    (alpha = |J|/(NR+NS))
#   GEN_SEED=1
#
#   # 采样设置
#   T=100000
#   SEED0=1
#   REPEATS=5                  (实际 seed=SEED0..SEED0+REPEATS-1)
#
#   # oracle / universe 预算（用于 correctness + χ²）
#   ORACLE_MAX_CHECKS=50000000
#   ORACLE_COLLECT_LIMIT=1000000
#   ORACLE_CAP=0
#
#   # methods/variants 子集（可选：空格分隔）
#   METHODS="ours aabb interval_tree"
#   VARIANTS="sampling enum_sampling"
#
set -euo pipefail
IFS=$'\n\t'

trap 'echo -e "[EXP1][FATAL] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2; exit 1' ERR

# --------------------------
# helpers
# --------------------------
log()  { echo -e "[EXP1] $*"; }
warn() { echo -e "[EXP1][WARN] $*" >&2; }
die()  { echo -e "[EXP1][FATAL] $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"
}

build_subdir_for_type() {
  local bt="$1"
  case "${bt}" in
    Release|release) echo "release" ;;
    Debug|debug) echo "debug" ;;
    RelWithDebInfo|relwithdebinfo) echo "relwithdebinfo" ;;
    MinSizeRel|minsizerel) echo "minsizerel" ;;
    *) die "Unsupported BUILD_TYPE='${bt}' (expected Release|Debug|RelWithDebInfo|MinSizeRel)" ;;
  esac
}

find_exe() {
  local build_dir="$1"
  local name="$2"
  local c
  for c in \
    "${build_dir}/${name}" \
    "${build_dir}/apps/${name}" \
    "${build_dir}/src/apps/${name}" \
    ; do
    if [[ -x "${c}" ]]; then
      echo "${c}"
      return 0
    fi
  done
  return 1
}

# Resolve repo root (script lives in <root>/run/run_exp1.sh)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXP_ID="exp1"

# --------------------------
# Parameters (override via env)
# --------------------------
BUILD_TYPE="${BUILD_TYPE:-Release}"
CLEAN_BUILD="${CLEAN_BUILD:-0}"
RUN_TESTS="${RUN_TESTS:-1}"

JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)}"
THREADS="${THREADS:-1}"

GEN="${GEN:-stripe}"
NR="${NR:-2000}"
NS="${NS:-2000}"
ALPHA="${ALPHA:-1}"
GEN_SEED="${GEN_SEED:-1}"

T="${T:-100000}"
SEED0="${SEED0:-1}"
REPEATS="${REPEATS:-5}"

ORACLE_MAX_CHECKS="${ORACLE_MAX_CHECKS:-50000000}"
ORACLE_COLLECT_LIMIT="${ORACLE_COLLECT_LIMIT:-1000000}"
ORACLE_CAP="${ORACLE_CAP:-0}"

# Methods / variants in this Dim=2 build
METHODS_DEFAULT=(ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection tsunami)
VARIANTS_DEFAULT=(sampling enum_sampling adaptive)

# Allow overriding METHODS / VARIANTS as space-separated strings
if [[ -n "${METHODS:-}" ]]; then
  read -r -a METHODS_ARR <<< "${METHODS}"
else
  METHODS_ARR=("${METHODS_DEFAULT[@]}")
fi

if [[ -n "${VARIANTS:-}" ]]; then
  read -r -a VARIANTS_ARR <<< "${VARIANTS}"
else
  VARIANTS_ARR=("${VARIANTS_DEFAULT[@]}")
fi

# --------------------------
# Paths (统一：build / run/temp / results/raw)
# --------------------------
BUILD_SUBDIR="$(build_subdir_for_type "${BUILD_TYPE}")"
BUILD_DIR="${ROOT}/build/${BUILD_SUBDIR}"

TEMP_ROOT="${ROOT}/run/temp/${EXP_ID}"
TEMP_LOG_DIR="${TEMP_ROOT}/logs"

RESULT_ROOT="${ROOT}/results/raw/${EXP_ID}"

RUN_INCLUDE_DIR="${ROOT}/run/include"
PY_PARSER="${RUN_INCLUDE_DIR}/exp1_parse_verify_logs.py"

DATASET_NAME="exp1_${GEN}_nr${NR}_ns${NS}_a${ALPHA}_g${GEN_SEED}"

# --------------------------
# preflight
# --------------------------
need_cmd cmake
need_cmd python3
need_cmd tee
if [[ "${RUN_TESTS}" == "1" ]]; then
  need_cmd ctest
fi
[[ -f "${PY_PARSER}" ]] || die "Missing python parser: ${PY_PARSER} (required by rule #4: no embedded Python in bash)"

# Thread caps (fairness + reproducibility)
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"

log "Repo root:        ${ROOT}"
log "Build dir:        ${BUILD_DIR} (BUILD_TYPE=${BUILD_TYPE}, CLEAN_BUILD=${CLEAN_BUILD})"
log "Temp dir:         ${TEMP_ROOT}"
log "Results dir:      ${RESULT_ROOT} (will overwrite)"
log "Parser:           ${PY_PARSER}"
log "Params: JOBS=${JOBS} THREADS=${THREADS} RUN_TESTS=${RUN_TESTS}"
log "Params: GEN=${GEN} NR=${NR} NS=${NS} ALPHA=${ALPHA} GEN_SEED=${GEN_SEED}"
log "Params: T=${T} SEED0=${SEED0} REPEATS=${REPEATS}"
log "Params: ORACLE_MAX_CHECKS=${ORACLE_MAX_CHECKS} ORACLE_COLLECT_LIMIT=${ORACLE_COLLECT_LIMIT} ORACLE_CAP=${ORACLE_CAP}"
log "Params: METHODS=${METHODS_ARR[*]}"
log "Params: VARIANTS=${VARIANTS_ARR[*]}"

# Reset temp/results (覆盖旧结果)
rm -rf "${TEMP_ROOT}"
mkdir -p "${TEMP_LOG_DIR}"

rm -rf "${RESULT_ROOT}"
mkdir -p "${RESULT_ROOT}"

# Guard: ensure oracle is feasible (EXP-1 必须小规模)
if [[ "${NR}" =~ ^[0-9]+$ && "${NS}" =~ ^[0-9]+$ && "${ORACLE_MAX_CHECKS}" =~ ^[0-9]+$ ]]; then
  ORACLE_CHECKS=$(( NR * NS ))
  if (( ORACLE_CHECKS > ORACLE_MAX_CHECKS )); then
    die "NR*NS=${ORACLE_CHECKS} exceeds ORACLE_MAX_CHECKS=${ORACLE_MAX_CHECKS}. Reduce NR/NS or raise ORACLE_MAX_CHECKS."
  fi
else
  log "Note: NR/NS/ORACLE_MAX_CHECKS not all integers; skipping NR*NS<=ORACLE_MAX_CHECKS guard."
fi

# Record environment + resolved config (先写 temp，最后同步到 results/raw/exp1)
{
  echo "EXP-1 manifest"
  echo "-------------"
  echo "date: $(date -Is || date)"
  echo "uname: $(uname -a || true)"
  echo "cmake: $(cmake --version | head -n1 || true)"
  if command -v git >/dev/null 2>&1 && [[ -d "${ROOT}/.git" ]]; then
    echo "git_head: $(git -C "${ROOT}" rev-parse --short HEAD || true)"
    echo "git_status:"
    git -C "${ROOT}" status --porcelain || true
  fi
  echo ""
  echo "paths:"
  echo "  repo_root=${ROOT}"
  echo "  build_dir=${BUILD_DIR}"
  echo "  temp_dir=${TEMP_ROOT}"
  echo "  results_dir=${RESULT_ROOT}"
  echo ""
  echo "dataset_source=synthetic"
  echo "dataset_label=${DATASET_NAME}"
  echo "gen=${GEN}"
  echo "n_r=${NR}"
  echo "n_s=${NS}"
  echo "alpha=${ALPHA}   # alpha=|J|/(n_r+n_s)"
  echo "gen_seed=${GEN_SEED}"
  echo ""
  echo "t=${T}"
  echo "seed0=${SEED0}"
  echo "repeats=${REPEATS}"
  echo ""
  echo "oracle_max_checks=${ORACLE_MAX_CHECKS}"
  echo "oracle_collect_limit=${ORACLE_COLLECT_LIMIT}"
  echo "oracle_cap=${ORACLE_CAP}"
  echo ""
  echo "methods=${METHODS_ARR[*]}"
  echo "variants=${VARIANTS_ARR[*]}"
  echo ""
  echo "thread_caps:"
  echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS}"
  echo "  OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS}"
  echo "  MKL_NUM_THREADS=${MKL_NUM_THREADS}"
  echo "  NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS}"
  echo "  VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS}"
  echo "  THREADS(flag)=${THREADS}"
} > "${TEMP_ROOT}/MANIFEST.txt"

# --------------------------
# Build
# --------------------------
log "Step 1/3: Configure + build"
if [[ "${CLEAN_BUILD}" == "1" ]]; then
  log "CLEAN_BUILD=1 -> rm -rf ${BUILD_DIR}"
  rm -rf "${BUILD_DIR}"
fi

cmake -S "${ROOT}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DSJS_BUILD_ROOT_APPS=ON \
  -DSJS_BUILD_TESTS=ON \
  2>&1 | tee "${TEMP_LOG_DIR}/cmake_configure.log"

cmake --build "${BUILD_DIR}" -j "${JOBS}" 2>&1 | tee "${TEMP_LOG_DIR}/cmake_build.log"

if [[ "${RUN_TESTS}" == "1" ]]; then
  log "Step 2/3: Run tests (ctest)"
  ctest --test-dir "${BUILD_DIR}" --output-on-failure 2>&1 | tee "${TEMP_LOG_DIR}/ctest.log"
else
  log "Step 2/3: RUN_TESTS=0 -> skip ctest"
fi

SJS_VERIFY="$(find_exe "${BUILD_DIR}" "sjs_verify" || true)"
[[ -n "${SJS_VERIFY}" && -x "${SJS_VERIFY}" ]] || die "sjs_verify not found under build dir: ${BUILD_DIR}"

# --------------------------
# Run EXP-1 (method × variant × seed)
# --------------------------
log "Step 3/3: Run sjs_verify for method×variant with repeats (seeds)"

FAIL_FILE="${TEMP_ROOT}/FAILURES.txt"
: > "${FAIL_FILE}"

total_jobs=0
failed_jobs=0

for m in "${METHODS_ARR[@]}"; do
  for v in "${VARIANTS_ARR[@]}"; do
    total_jobs=$((total_jobs + 1))
    log_file="${TEMP_LOG_DIR}/${m}__${v}.log"

    log "(${total_jobs}) method=${m} variant=${v} -> ${log_file}"

    # Keep going even if one combo fails (failures are recorded & surfaced).
    if ! "${SJS_VERIFY}" \
      --dataset_source=synthetic --gen="${GEN}" --dataset="${DATASET_NAME}" \
      --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
      --method="${m}" --variant="${v}" --t="${T}" --seed="${SEED0}" --repeats="${REPEATS}" \
      --threads="${THREADS}" \
      --oracle_max_checks="${ORACLE_MAX_CHECKS}" --oracle_collect_limit="${ORACLE_COLLECT_LIMIT}" --oracle_cap="${ORACLE_CAP}" \
      > "${log_file}" 2>&1
    then
      rc=$?
      failed_jobs=$((failed_jobs + 1))
      echo "[EXP1] FAIL rc=${rc} method=${m} variant=${v} log=${log_file}" | tee -a "${FAIL_FILE}" >&2
      continue
    fi
  done
done

# Parse logs -> CSV tables (no embedded Python; call run/include/*.py)
log "Parsing logs -> CSV (written to temp)"
python3 "${PY_PARSER}" \
  --log_dir "${TEMP_LOG_DIR}" \
  --out_dir "${TEMP_ROOT}" \
  --fail_file "${FAIL_FILE}" \
  2>&1 | tee "${TEMP_LOG_DIR}/parse_logs.log"

# --------------------------
# Sync: temp → results/raw/exp1 (覆盖旧结果)
# --------------------------
log "Syncing artifacts to ${RESULT_ROOT} (overwrite)"

# results/raw/exp1 当前已是空目录；直接整体复制即可
cp -a "${TEMP_ROOT}/." "${RESULT_ROOT}/"

log "EXP-1 DONE"
log "Artifacts:"
log "  - Results: ${RESULT_ROOT}"
log "  - Logs:    ${RESULT_ROOT}/logs"
log "  - Tables:  ${RESULT_ROOT}/exp1_quality_raw.csv"
log "             ${RESULT_ROOT}/exp1_quality_summary.csv"
log "  - Manifest:${RESULT_ROOT}/MANIFEST.txt"
log "  - Temp:    ${TEMP_ROOT}"

if [[ "${failed_jobs}" -gt 0 ]]; then
  warn "${failed_jobs}/${total_jobs} (method×variant) jobs failed. See ${RESULT_ROOT}/FAILURES.txt"
  exit 2
fi

exit 0
