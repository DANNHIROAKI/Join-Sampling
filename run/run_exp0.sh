#!/usr/bin/env bash
# run/run_exp0.sh
#
# EXP-0：基础 sanity（单元级 + 端到端 smoke）
#
# 本脚本需与 Docs/Experiment/EXP-0.md 的流程严格对齐：
#   Step 1/5: Clean build + compile
#   Step 2/5: ctest
#   Step 3/5: 生成小数据集到磁盘（验证生成器 + IO）
#   Step 4/5: sjs_run 端到端 smoke（synthetic & binary）
#   Step 5/5: sjs_verify（oracle + sample-quality 最低限度校验）
#
# 统一目录规范（对齐你的全局要求）：
#   - Build 统一放在 <repo_root>/build/ 下
#   - 运行过程中产生的日志/中间文件先写到 <repo_root>/run/temp/exp0/
#   - 最终实验结果同步到 <repo_root>/results/raw/exp0/（每次覆盖旧结果）
#   - 生成的数据集（bin/csv/report）同步到 <repo_root>/data/synthetic/exp0/（每次覆盖旧数据）
#
# 用法：
#   bash run/run_exp0.sh
#
# 可选环境变量：
#   BUILD_TYPE=Release|Debug
#   CLEAN_BUILD=0|1            (默认 1：清理并重建 build/<type>)
#   JOBS=8
#   THREADS=1
#   NR=2000  NS=2000  ALPHA=1e-3  GEN_SEED=1
#   T=20000  SEED=1
#   ORACLE_MAX_CHECKS=50000000
#   RUN_ASAN=0|1
#
#   # 对 (est) 的 count 估计是否当成 gate
#   CHECK_EST_COUNT=0|1        (默认 0：只 warning，不 fail)
#   EST_REL_ERR_MAX=0.5        (CHECK_EST_COUNT=1 时使用，默认 0.5)
#   EST_REL_ERR_WARN=0.5       (CHECK_EST_COUNT=0 时使用，默认 0.5)
#
set -euo pipefail
IFS=$'\n\t'

trap 'echo -e "[EXP0][FATAL] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2; exit 1' ERR

# --------------------------
# helpers
# --------------------------
log()  { echo -e "[EXP0] $*"; }
warn() { echo -e "[EXP0][WARN] $*" >&2; }
die()  { echo -e "[EXP0][FATAL] $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"
}

# Extract key=value-like fields from a summary line such as:
# "count=4 (exact)  oracle=4  rel_err=0"
extract_field() {
  local line="$1"
  local key="$2"
  # Split on spaces, '=', '(', ')' -> tokens: count 4 exact oracle 4 rel_err 0
  echo "${line}" | awk -v k="${key}" -F'[ =()]+' '{
    for (i=1; i<=NF; ++i) if ($i==k) { print $(i+1); exit }
  }'
}

# Extract count kind: "exact" or "est" from:
# "count=5.18 (est) oracle=4 rel_err=..."
extract_count_kind() {
  local line="$1"
  echo "${line}" | awk -F'[ =()]+' '{
    for (i=1; i<=NF; ++i) if ($i=="count") { print $(i+2); exit }
  }'
}

float_is_zero() {
  local x="$1"
  awk -v v="${x}" 'BEGIN { exit ((v+0)==0 ? 0 : 1) }'
}

float_leq() {
  local x="$1"
  local y="$2"
  awk -v a="${x}" -v b="${y}" 'BEGIN { exit ((a+0) <= (b+0) ? 0 : 1) }'
}

float_eq() {
  local x="$1"
  local y="$2"
  awk -v a="${x}" -v b="${y}" 'BEGIN { exit ((a+0) == (b+0) ? 0 : 1) }'
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

# Resolve repo root (script lives in <root>/run/run_exp0.sh)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXP_ID="exp0"
DATASET_NAME="exp0"   # used for generated file names

# --------------------------
# Parameters (override via env)
# --------------------------
BUILD_TYPE="${BUILD_TYPE:-Release}"
CLEAN_BUILD="${CLEAN_BUILD:-1}"

JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)}"
THREADS="${THREADS:-1}"

NR="${NR:-2000}"
NS="${NS:-2000}"
ALPHA="${ALPHA:-1e-3}"
GEN_SEED="${GEN_SEED:-1}"

T="${T:-20000}"
SEED="${SEED:-1}"

ORACLE_MAX_CHECKS="${ORACLE_MAX_CHECKS:-50000000}"
RUN_ASAN="${RUN_ASAN:-0}"

CHECK_EST_COUNT="${CHECK_EST_COUNT:-0}"
EST_REL_ERR_MAX="${EST_REL_ERR_MAX:-0.5}"
EST_REL_ERR_WARN="${EST_REL_ERR_WARN:-0.5}"

# --------------------------
# Paths (统一：build / run/temp / results/raw)
# --------------------------
BUILD_SUBDIR="$(build_subdir_for_type "${BUILD_TYPE}")"
BUILD_DIR="${ROOT}/build/${BUILD_SUBDIR}"
ASAN_BUILD_DIR="${ROOT}/build/asan"

TEMP_ROOT="${ROOT}/run/temp/${EXP_ID}"
TEMP_LOG_DIR="${TEMP_ROOT}/logs"
TEMP_DATA_DIR="${TEMP_ROOT}/data"

RESULT_ROOT="${ROOT}/results/raw/${EXP_ID}"

DATA_DIR="${ROOT}/data/synthetic/${EXP_ID}"

# --------------------------
# preflight
# --------------------------
need_cmd cmake
need_cmd ctest
need_cmd grep
need_cmd awk
need_cmd tee

log "Repo root:        ${ROOT}"
log "Build dir:        ${BUILD_DIR} (BUILD_TYPE=${BUILD_TYPE}, CLEAN_BUILD=${CLEAN_BUILD})"
log "Temp dir:         ${TEMP_ROOT}"
log "Results dir:      ${RESULT_ROOT} (will overwrite)"
log "Data dir:         ${DATA_DIR} (will overwrite)"
log "Params: JOBS=${JOBS} THREADS=${THREADS}"
log "Params: NR=${NR} NS=${NS} ALPHA=${ALPHA} GEN_SEED=${GEN_SEED} T=${T} SEED=${SEED}"
log "Params: ORACLE_MAX_CHECKS=${ORACLE_MAX_CHECKS} RUN_ASAN=${RUN_ASAN}"
log "Params: CHECK_EST_COUNT=${CHECK_EST_COUNT} EST_REL_ERR_MAX=${EST_REL_ERR_MAX} EST_REL_ERR_WARN=${EST_REL_ERR_WARN}"

# Reset temp/results/data (覆盖旧结果)
rm -rf "${TEMP_ROOT}"
mkdir -p "${TEMP_LOG_DIR}" "${TEMP_DATA_DIR}"

rm -rf "${RESULT_ROOT}"
mkdir -p "${RESULT_ROOT}"

rm -rf "${DATA_DIR}"
mkdir -p "${DATA_DIR}"

# sanity: oracle checks budget (only if NR/NS are integers)
if [[ "${NR}" =~ ^[0-9]+$ && "${NS}" =~ ^[0-9]+$ && "${ORACLE_MAX_CHECKS}" =~ ^[0-9]+$ ]]; then
  ORACLE_CHECKS=$(( NR * NS ))
  if (( ORACLE_CHECKS > ORACLE_MAX_CHECKS )); then
    die "NR*NS=${ORACLE_CHECKS} exceeds ORACLE_MAX_CHECKS=${ORACLE_MAX_CHECKS}. Reduce NR/NS or raise ORACLE_MAX_CHECKS."
  fi
else
  log "Note: NR/NS/ORACLE_MAX_CHECKS not all integers; skipping NR*NS<=ORACLE_MAX_CHECKS guard."
fi

# record environment (先写 temp，最后同步到 results/raw/exp0/logs/)
{
  echo "date: $(date -Is || date)"
  echo "uname: $(uname -a || true)"
  echo "cmake: $(cmake --version | head -n1 || true)"
  echo "ctest: $(ctest --version | head -n1 || true)"
  if command -v git >/dev/null 2>&1 && [[ -d "${ROOT}/.git" ]]; then
    echo "git_head: $(git -C "${ROOT}" rev-parse --short HEAD || true)"
    echo "git_status:"
    git -C "${ROOT}" status --porcelain || true
  fi
  echo "resolved_params:"
  echo "  BUILD_TYPE=${BUILD_TYPE}"
  echo "  CLEAN_BUILD=${CLEAN_BUILD}"
  echo "  JOBS=${JOBS}"
  echo "  THREADS=${THREADS}"
  echo "  NR=${NR}"
  echo "  NS=${NS}"
  echo "  ALPHA=${ALPHA}"
  echo "  GEN_SEED=${GEN_SEED}"
  echo "  T=${T}"
  echo "  SEED=${SEED}"
  echo "  ORACLE_MAX_CHECKS=${ORACLE_MAX_CHECKS}"
  echo "  RUN_ASAN=${RUN_ASAN}"
  echo "  CHECK_EST_COUNT=${CHECK_EST_COUNT}"
  echo "  EST_REL_ERR_MAX=${EST_REL_ERR_MAX}"
  echo "  EST_REL_ERR_WARN=${EST_REL_ERR_WARN}"
} > "${TEMP_LOG_DIR}/env.txt"

# --------------------------
# Step 1/5) clean build + compile
# --------------------------
log "Step 1/5: Clean build + compile"

if [[ "${CLEAN_BUILD}" == "1" ]]; then
  rm -rf "${BUILD_DIR}"
fi
mkdir -p "${BUILD_DIR}"

pushd "${BUILD_DIR}" >/dev/null

cmake "${ROOT}" \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DSJS_BUILD_ROOT_APPS=ON \
  -DSJS_BUILD_TESTS=ON \
  2>&1 | tee "${TEMP_LOG_DIR}/cmake_configure.log"

cmake --build . -j "${JOBS}" 2>&1 | tee "${TEMP_LOG_DIR}/cmake_build.log"

popd >/dev/null

# executables (位置可能因 CMake 配置而不同，做一层探测)
SJS_GEN="$(find_exe "${BUILD_DIR}" "sjs_gen_dataset" || true)"
SJS_RUN="$(find_exe "${BUILD_DIR}" "sjs_run" || true)"
SJS_VERIFY="$(find_exe "${BUILD_DIR}" "sjs_verify" || true)"

[[ -n "${SJS_GEN}" && -x "${SJS_GEN}" ]]       || die "Missing executable: sjs_gen_dataset under ${BUILD_DIR}"
[[ -n "${SJS_RUN}" && -x "${SJS_RUN}" ]]       || die "Missing executable: sjs_run under ${BUILD_DIR}"
[[ -n "${SJS_VERIFY}" && -x "${SJS_VERIFY}" ]] || die "Missing executable: sjs_verify under ${BUILD_DIR}"

# --------------------------
# Step 2/5) unit-level sanity: ctest
# --------------------------
log "Step 2/5: Run unit tests (ctest)"
pushd "${BUILD_DIR}" >/dev/null
ctest --output-on-failure 2>&1 | tee "${TEMP_LOG_DIR}/ctest.log"
popd >/dev/null

# --------------------------
# Step 3/5) end-to-end: generate small dataset to disk (binary + csv)
# --------------------------
log "Step 3/5: Generate a small dataset (binary + csv) via sjs_gen_dataset"
# Note: generator name "stripe" is an alias of stripe_ctrl_alpha.
"${SJS_GEN}" \
  --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}" \
  --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
  --out_dir="${TEMP_DATA_DIR}" --write_csv=1 \
  --threads="${THREADS}" \
  2>&1 | tee "${TEMP_LOG_DIR}/gen_dataset.log"

R_BIN_TEMP="${TEMP_DATA_DIR}/${DATASET_NAME}_R.bin"
S_BIN_TEMP="${TEMP_DATA_DIR}/${DATASET_NAME}_S.bin"
[[ -f "${R_BIN_TEMP}" ]] || die "Expected output not found: ${R_BIN_TEMP}"
[[ -f "${S_BIN_TEMP}" ]] || die "Expected output not found: ${S_BIN_TEMP}"

# 生成器报告（若存在则检查；不同版本可能不写 report，不强制 fail）
GEN_REPORT_TEMP="${TEMP_DATA_DIR}/${DATASET_NAME}_gen_report.json"
if [[ ! -f "${GEN_REPORT_TEMP}" ]]; then
  warn "Generator report not found: ${GEN_REPORT_TEMP} (will continue; check sjs_gen_dataset implementation if needed)."
fi

# 同步到 canonical 数据目录：data/synthetic/exp0（覆盖旧数据）
cp -a "${TEMP_DATA_DIR}/." "${DATA_DIR}/."

R_BIN="${DATA_DIR}/${DATASET_NAME}_R.bin"
S_BIN="${DATA_DIR}/${DATASET_NAME}_S.bin"
[[ -f "${R_BIN}" ]] || die "Expected output not found after sync: ${R_BIN}"
[[ -f "${S_BIN}" ]] || die "Expected output not found after sync: ${S_BIN}"

# --------------------------
# Step 4/5) end-to-end: sjs_run smoke (synthetic + binary)
# --------------------------
log "Step 4/5: End-to-end smoke via sjs_run"

run_one() {
  local tag="$1"; shift
  local log_file="${TEMP_LOG_DIR}/${tag}.log"
  log "  - Running: ${tag}"
  "${SJS_RUN}" "$@" 2>&1 | tee "${log_file}"
}

# (A) Synthetic on-the-fly, run all 3 variants for a representative method (ours)
for VARIANT in sampling enum_sampling adaptive; do
  OUT_DIR="${TEMP_ROOT}/run_ours_${VARIANT}_synthetic"
  mkdir -p "${OUT_DIR}"
  run_one "run_ours_${VARIANT}_synthetic" \
    --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}_synthetic" \
    --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
    --method=ours --variant="${VARIANT}" --t="${T}" --seed="${SEED}" --repeats=1 \
    --out_dir="${OUT_DIR}" --write_samples=1 \
    --threads="${THREADS}"

  [[ -f "${OUT_DIR}/run.csv" ]] || die "Missing run.csv for ${OUT_DIR}"
done

# (B) Binary input path (tests IO), just sampling variant (ours)
OUT_DIR="${TEMP_ROOT}/run_ours_sampling_binary"
mkdir -p "${OUT_DIR}"
run_one "run_ours_sampling_binary" \
  --dataset_source=binary --dataset="${DATASET_NAME}_binary" \
  --path_r="${R_BIN}" --path_s="${S_BIN}" \
  --method=ours --variant=sampling --t="${T}" --seed="${SEED}" --repeats=1 \
  --out_dir="${OUT_DIR}" --write_samples=1 \
  --threads="${THREADS}"
[[ -f "${OUT_DIR}/run.csv" ]] || die "Missing run.csv for ${OUT_DIR}"

# (C) Smoke r_tree sampling once (regression guard; should never crash)
OUT_DIR="${TEMP_ROOT}/run_r_tree_sampling_synthetic"
mkdir -p "${OUT_DIR}"
run_one "run_r_tree_sampling_synthetic" \
  --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}_rtree_smoke" \
  --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
  --method=r_tree --variant=sampling --t="${T}" --seed="${SEED}" --repeats=1 \
  --out_dir="${OUT_DIR}" --write_samples=1 \
  --threads="${THREADS}"
[[ -f "${OUT_DIR}/run.csv" ]] || die "Missing run.csv for ${OUT_DIR}"

# --------------------------
# Step 5/5) end-to-end: sjs_verify (oracle + quality); enforce correctness
# --------------------------
log "Step 5/5: Correctness & sample-quality sanity via sjs_verify"

verify_one() {
  local method="$1"
  local variant="$2"

  local tag="verify_${method}_${variant}"
  local log_file="${TEMP_LOG_DIR}/${tag}.log"

  log "  - Verifying: method=${method} variant=${variant}"
  "${SJS_VERIFY}" \
    --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}_verify" \
    --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
    --method="${method}" --variant="${variant}" --t="${T}" --seed="${SEED}" --repeats=1 \
    --threads="${THREADS}" \
    --oracle_max_checks="${ORACLE_MAX_CHECKS}" \
    2>&1 | tee "${log_file}"

  # Require "missing_in_universe=0"
  local missing
  missing="$(grep -E "missing_in_universe=" "${log_file}" | tail -n1 | awk -F= '{print $2}' | tr -d '[:space:]' || true)"
  if [[ -z "${missing}" ]]; then
    die "No 'missing_in_universe=' line found in ${log_file}. (Was universe collection skipped?)"
  fi
  if [[ "${missing}" != "0" ]]; then
    die "Correctness FAILED for ${method}/${variant}: missing_in_universe=${missing} (see ${log_file})"
  fi

  # Parse summary line: count / (exact|est) / oracle / rel_err
  local sum_line count oracle rel_err count_kind
  sum_line="$(grep -E "count=.*oracle=" "${log_file}" | tail -n1 || true)"
  if [[ -z "${sum_line}" ]]; then
    die "No summary line containing 'count=' and 'oracle=' found in ${log_file}."
  fi

  count="$(extract_field "${sum_line}" "count")"
  oracle="$(extract_field "${sum_line}" "oracle")"
  rel_err="$(extract_field "${sum_line}" "rel_err")"
  count_kind="$(extract_count_kind "${sum_line}")"

  if [[ -z "${count}" || -z "${oracle}" || -z "${rel_err}" || -z "${count_kind}" ]]; then
    die "Failed to parse count/oracle/rel_err/count_kind from: '${sum_line}' (see ${log_file})"
  fi

  if [[ "${count_kind}" == "exact" ]]; then
    # Exact methods must match oracle.
    if ! float_eq "${count}" "${oracle}"; then
      die "Correctness FAILED for ${method}/${variant}: count=${count} != oracle=${oracle} (exact) (see ${log_file})"
    fi
    if ! float_is_zero "${rel_err}"; then
      die "Correctness FAILED for ${method}/${variant}: rel_err=${rel_err} != 0 (exact) (see ${log_file})"
    fi
  else
    # Estimated methods (e.g., rejection) do not guarantee count==oracle.
    # Default: warn only; optionally gate with CHECK_EST_COUNT=1.
    if [[ "${CHECK_EST_COUNT}" == "1" ]]; then
      if ! float_leq "${rel_err}" "${EST_REL_ERR_MAX}"; then
        die "EST count check FAILED for ${method}/${variant}: rel_err=${rel_err} > ${EST_REL_ERR_MAX} (see ${log_file})"
      fi
    else
      if ! float_leq "${rel_err}" "${EST_REL_ERR_WARN}"; then
        warn "EST count warning for ${method}/${variant}: rel_err=${rel_err} > ${EST_REL_ERR_WARN} (count_kind=${count_kind})."
      else
        log "    NOTE: ${method}/${variant} reports estimated count (count_kind=${count_kind}), rel_err=${rel_err} (not gated)."
      fi
    fi
  fi
}

# Verify all registered methods in sampling mode (core claim)
METHODS=(ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection tsunami)
for m in "${METHODS[@]}"; do
  verify_one "${m}" "sampling"
done

# Optional: verify binary path as well (only if sjs_verify supports --path_r/--path_s)
if "${SJS_VERIFY}" --help 2>&1 | grep -q -- "--path_r"; then
  log "  - Verifying binary dataset path (ours/sampling)"
  tag="verify_ours_sampling_binary"
  log_file="${TEMP_LOG_DIR}/${tag}.log"
  "${SJS_VERIFY}" \
    --dataset_source=binary --dataset="${DATASET_NAME}_binary_verify" \
    --path_r="${R_BIN}" --path_s="${S_BIN}" \
    --method=ours --variant=sampling --t="${T}" --seed="${SEED}" --repeats=1 \
    --threads="${THREADS}" \
    --oracle_max_checks="${ORACLE_MAX_CHECKS}" \
    2>&1 | tee "${log_file}"

  missing="$(grep -E "missing_in_universe=" "${log_file}" | tail -n1 | awk -F= '{print $2}' | tr -d '[:space:]' || true)"
  [[ "${missing}" == "0" ]] || die "Binary verify FAILED: missing_in_universe=${missing} (see ${log_file})"
else
  log "  - Skipping binary verify (sjs_verify --help does not show --path_r)."
fi

# --------------------------
# Optional) ASAN/UBSAN build + ctest + r_tree verify
# --------------------------
if [[ "${RUN_ASAN}" == "1" ]]; then
  log "Optional: RUN_ASAN=1 → ASAN/UBSAN build + tests + r_tree verify"

  rm -rf "${ASAN_BUILD_DIR}"
  mkdir -p "${ASAN_BUILD_DIR}"
  pushd "${ASAN_BUILD_DIR}" >/dev/null

  cmake "${ROOT}" -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DSJS_BUILD_ROOT_APPS=ON -DSJS_BUILD_TESTS=ON \
    -DCMAKE_CXX_FLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer" \
    -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined" \
    2>&1 | tee "${TEMP_LOG_DIR}/cmake_asan_configure.log"

  cmake --build . -j "${JOBS}" 2>&1 | tee "${TEMP_LOG_DIR}/cmake_asan_build.log"

  export ASAN_OPTIONS="halt_on_error=1:abort_on_error=1:detect_leaks=1:strict_string_checks=1:symbolize=1"
  export UBSAN_OPTIONS="halt_on_error=1:abort_on_error=1:print_stacktrace=1"

  ctest --output-on-failure 2>&1 | tee "${TEMP_LOG_DIR}/ctest_asan.log"

  # r_tree verify under ASAN (the historical crash site)
  ASAN_SJS_VERIFY="$(find_exe "${ASAN_BUILD_DIR}" "sjs_verify" || true)"
  [[ -n "${ASAN_SJS_VERIFY}" && -x "${ASAN_SJS_VERIFY}" ]] || die "Missing asan executable: sjs_verify under ${ASAN_BUILD_DIR}"

  "${ASAN_SJS_VERIFY}" \
    --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}_verify_asan" \
    --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
    --method=r_tree --variant=sampling --t="${T}" --seed="${SEED}" --repeats=1 \
    --threads="${THREADS}" \
    --oracle_max_checks="${ORACLE_MAX_CHECKS}" \
    2>&1 | tee "${TEMP_LOG_DIR}/verify_asan_r_tree_sampling.log"

  popd >/dev/null
fi

# --------------------------
# Sync: temp → results/raw/exp0 (覆盖旧结果)
# --------------------------
log "Syncing artifacts to ${RESULT_ROOT} (overwrite)"

mkdir -p "${RESULT_ROOT}/logs"
cp -a "${TEMP_LOG_DIR}/." "${RESULT_ROOT}/logs/"

# Copy run output directories (run_*)
shopt -s nullglob
for d in "${TEMP_ROOT}"/run_*; do
  if [[ -d "${d}" ]]; then
    cp -a "${d}" "${RESULT_ROOT}/"
  fi
done
shopt -u nullglob

log "EXP-0 PASSED ✅"
log "Artifacts:"
log "  - Logs:    ${RESULT_ROOT}/logs"
log "  - Results: ${RESULT_ROOT}"
log "  - Data:    ${DATA_DIR}"
log "  - Temp:    ${TEMP_ROOT}"
