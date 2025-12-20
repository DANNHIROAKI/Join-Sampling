#!/usr/bin/env bash
# run/run_exp0.sh
#
# EXP-0：基础 sanity（单元级 + 端到端 smoke）
#  1) Clean build + compile（Release/Debug）
#  2) ctest（几何/事件扫描/Join oracle/采样质量/baseline smoke/rtree split 回归）
#  3) 生成小数据集到磁盘（binary + csv）
#  4) 端到端 smoke：sjs_run
#       - ours: sampling / enum_sampling / adaptive（synthetic）
#       - r_tree: sampling（synthetic）  ← 之前的崩溃回归重点
#       - ours: sampling（binary IO）
#       - r_tree: sampling（binary IO）
#  5) 端到端 correctness：sjs_verify（oracle + quality）
#       - 所有 methods（sampling）必须满足：
#           missing_in_universe=0 且 count==oracle
#  6) （可选，默认开启）ASAN/UBSAN pass：
#       - 再建一个 sanitizer build，跑 ctest
#       - 以及跑 verify（默认 ours + r_tree；可选全方法）
#
# 用法（仓库根目录或任意目录均可）：
#   bash run/run_exp0.sh
#
# 可选环境变量（不设则用默认值）：
#   BUILD_TYPE=Release|Debug|RelWithDebInfo
#   JOBS=8
#   THREADS=1
#   NR=2000  NS=2000  ALPHA=1e-3  GEN_SEED=1
#   T=20000  SEED=1  REPEATS=1
#   ORACLE_MAX_CHECKS=50000000
#
#   RUN_ASAN=1|0              # 默认 1（建议在 EXP-0 打开）
#   ASAN_JOBS=4               # sanitizer build 的并行度（默认同 JOBS）
#   ASAN_VERIFY_ALL=1|0       # 默认 0；设 1 会在 ASAN 下 verify 全方法
#
set -Eeuo pipefail
IFS=$'\n\t'

# --------------------------
# helpers
# --------------------------
log()  { echo -e "[EXP0] $*"; }
die()  { echo -e "[EXP0][FATAL] $*" >&2; exit 1; }

on_err() {
  local ec=$?
  echo -e "[EXP0][FATAL] Error (exit=${ec}) at line ${BASH_LINENO[0]}: ${BASH_COMMAND}" >&2
  exit "${ec}"
}
trap on_err ERR

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"
}

# Resolve repo root (script lives in <root>/run/run_exp0.sh)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --------------------------
# parameters (override via env)
# --------------------------
BUILD_TYPE="${BUILD_TYPE:-Release}"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)}"
THREADS="${THREADS:-1}"

NR="${NR:-2000}"
NS="${NS:-2000}"
ALPHA="${ALPHA:-1e-3}"
GEN_SEED="${GEN_SEED:-1}"

T="${T:-20000}"
SEED="${SEED:-1}"
REPEATS="${REPEATS:-1}"

ORACLE_MAX_CHECKS="${ORACLE_MAX_CHECKS:-50000000}"

RUN_ASAN="${RUN_ASAN:-1}"
ASAN_JOBS="${ASAN_JOBS:-${JOBS}}"
ASAN_VERIFY_ALL="${ASAN_VERIFY_ALL:-0}"

# --------------------------
# paths
# --------------------------
BUILD_DIR="${ROOT}/build_exp0"
OUT_ROOT="${ROOT}/results/raw/exp0"
DATA_DIR="${ROOT}/data/synthetic/exp0"
LOG_DIR="${OUT_ROOT}/logs"

DATASET_NAME="exp0"  # used for generated file names

mkdir -p "${OUT_ROOT}" "${DATA_DIR}" "${LOG_DIR}"

# --------------------------
# preflight
# --------------------------
need_cmd cmake
need_cmd ctest
need_cmd grep
need_cmd awk
need_cmd tee

log "Repo root: ${ROOT}"
log "Build dir: ${BUILD_DIR}"
log "Output dir: ${OUT_ROOT}"
log "Data dir: ${DATA_DIR}"
log "Params: BUILD_TYPE=${BUILD_TYPE} JOBS=${JOBS} THREADS=${THREADS}"
log "Params: NR=${NR} NS=${NS} ALPHA=${ALPHA} GEN_SEED=${GEN_SEED}"
log "Params: T=${T} SEED=${SEED} REPEATS=${REPEATS}"
log "Params: ORACLE_MAX_CHECKS=${ORACLE_MAX_CHECKS}"
log "Params: RUN_ASAN=${RUN_ASAN} ASAN_JOBS=${ASAN_JOBS} ASAN_VERIFY_ALL=${ASAN_VERIFY_ALL}"

# record environment
{
  echo "date: $(date -Is || date)"
  echo "uname: $(uname -a || true)"
  echo "c++: $(c++ --version 2>/dev/null | head -n1 || true)"
  echo "cmake: $(cmake --version | head -n1 || true)"
  echo "ctest: $(ctest --version | head -n1 || true)"
  if command -v git >/dev/null 2>&1 && git -C "${ROOT}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "git_head: $(git -C "${ROOT}" rev-parse HEAD)"
    echo "git_dirty_files: $(git -C "${ROOT}" status --porcelain | wc -l | tr -d ' ')"
  fi
} > "${LOG_DIR}/env.txt"

# --------------------------
# shared helpers for run/verify
# --------------------------
run_one() {
  local sjs_run_bin="$1"; shift
  local log_dir="$1"; shift
  local tag="$1"; shift

  local log_file="${log_dir}/${tag}.log"
  log "  - Running: ${tag}"
  "${sjs_run_bin}" "$@" 2>&1 | tee "${log_file}"
}

verify_one() {
  local sjs_verify_bin="$1"; shift
  local log_dir="$1"; shift
  local method="$1"; shift
  local variant="$1"; shift

  local tag="verify_${method}_${variant}"
  local log_file="${log_dir}/${tag}.log"

  log "  - Verifying: method=${method} variant=${variant}"
  "${sjs_verify_bin}" \
    --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}_verify" \
    --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
    --method="${method}" --variant="${variant}" --t="${T}" --seed="${SEED}" --repeats="${REPEATS}" \
    --oracle_max_checks="${ORACLE_MAX_CHECKS}" \
    --threads="${THREADS}" \
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

  # Require count == oracle (EXP-0 strict correctness gate)
  # Example line:
  #   count=4 (exact)  oracle=4  rel_err=0
  local count_line count_val oracle_val
  count_line="$(grep -E "count=[0-9]+.*oracle=[0-9]+" "${log_file}" | tail -n1 || true)"
  if [[ -z "${count_line}" ]]; then
    die "No 'count=... oracle=...' line found in ${log_file}."
  fi
  count_val="$(echo "${count_line}"  | sed -n 's/.*count=\([0-9]\+\).*/\1/p' | head -n1 || true)"
  oracle_val="$(echo "${count_line}" | sed -n 's/.*oracle=\([0-9]\+\).*/\1/p' | head -n1 || true)"
  if [[ -z "${count_val}" || -z "${oracle_val}" ]]; then
    die "Failed to parse count/oracle from line: ${count_line}"
  fi
  if [[ "${count_val}" != "${oracle_val}" ]]; then
    die "Correctness FAILED for ${method}/${variant}: count=${count_val} != oracle=${oracle_val} (see ${log_file})"
  fi
}

# --------------------------
# 1) clean build + compile
# --------------------------
log "Step 1/6: Clean build + compile"
rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"
pushd "${BUILD_DIR}" >/dev/null

cmake "${ROOT}" \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DSJS_BUILD_ROOT_APPS=ON \
  -DSJS_BUILD_TESTS=ON \
  2>&1 | tee "${LOG_DIR}/cmake_configure.log"

cmake --build . -j "${JOBS}" 2>&1 | tee "${LOG_DIR}/cmake_build.log"

popd >/dev/null

# executables
SJS_GEN="${BUILD_DIR}/sjs_gen_dataset"
SJS_RUN="${BUILD_DIR}/sjs_run"
SJS_VERIFY="${BUILD_DIR}/sjs_verify"

[[ -x "${SJS_GEN}" ]]    || die "Missing executable: ${SJS_GEN}"
[[ -x "${SJS_RUN}" ]]    || die "Missing executable: ${SJS_RUN}"
[[ -x "${SJS_VERIFY}" ]] || die "Missing executable: ${SJS_VERIFY}"

# --------------------------
# 2) unit-level sanity: ctest
# --------------------------
log "Step 2/6: Run unit tests (ctest)"
pushd "${BUILD_DIR}" >/dev/null
ctest --output-on-failure 2>&1 | tee "${LOG_DIR}/ctest.log"
popd >/dev/null

# --------------------------
# 3) generate dataset to disk (binary + csv)
# --------------------------
log "Step 3/6: Generate a small dataset (binary + csv) via sjs_gen_dataset"
# Note: generator name "stripe" is an alias of stripe_ctrl_alpha.
"${SJS_GEN}" \
  --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}" \
  --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
  --out_dir="${DATA_DIR}" --write_csv=1 \
  --threads="${THREADS}" \
  2>&1 | tee "${LOG_DIR}/gen_dataset.log"

R_BIN="${DATA_DIR}/${DATASET_NAME}_R.bin"
S_BIN="${DATA_DIR}/${DATASET_NAME}_S.bin"
[[ -f "${R_BIN}" ]] || die "Expected output not found: ${R_BIN}"
[[ -f "${S_BIN}" ]] || die "Expected output not found: ${S_BIN}"

# --------------------------
# 4) end-to-end smoke via sjs_run (synthetic + binary)
# --------------------------
log "Step 4/6: End-to-end smoke via sjs_run"

# (A) Synthetic on-the-fly: ours (all 3 variants)
for VARIANT in sampling enum_sampling adaptive; do
  OUT_DIR="${OUT_ROOT}/run_ours_${VARIANT}_synthetic"
  mkdir -p "${OUT_DIR}"
  run_one "${SJS_RUN}" "${LOG_DIR}" "run_ours_${VARIANT}_synthetic" \
    --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}_synthetic" \
    --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
    --method=ours --variant="${VARIANT}" --t="${T}" --seed="${SEED}" --repeats="${REPEATS}" \
    --out_dir="${OUT_DIR}" --write_samples=1 \
    --threads="${THREADS}"
  [[ -f "${OUT_DIR}/run.csv" ]] || die "Missing run.csv for ${OUT_DIR}"
done

# (B) Synthetic on-the-fly: r_tree sampling (regression focus)
OUT_DIR="${OUT_ROOT}/run_rtree_sampling_synthetic"
mkdir -p "${OUT_DIR}"
run_one "${SJS_RUN}" "${LOG_DIR}" "run_rtree_sampling_synthetic" \
  --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}_synthetic" \
  --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
  --method=r_tree --variant=sampling --t="${T}" --seed="${SEED}" --repeats="${REPEATS}" \
  --out_dir="${OUT_DIR}" --write_samples=1 \
  --threads="${THREADS}"
[[ -f "${OUT_DIR}/run.csv" ]] || die "Missing run.csv for ${OUT_DIR}"

# (C) Binary input path (tests IO): ours sampling
OUT_DIR="${OUT_ROOT}/run_ours_sampling_binary"
mkdir -p "${OUT_DIR}"
run_one "${SJS_RUN}" "${LOG_DIR}" "run_ours_sampling_binary" \
  --dataset_source=binary --dataset="${DATASET_NAME}_binary" \
  --path_r="${R_BIN}" --path_s="${S_BIN}" \
  --method=ours --variant=sampling --t="${T}" --seed="${SEED}" --repeats="${REPEATS}" \
  --out_dir="${OUT_DIR}" --write_samples=1 \
  --threads="${THREADS}"
[[ -f "${OUT_DIR}/run.csv" ]] || die "Missing run.csv for ${OUT_DIR}"

# (D) Binary input path (tests IO): r_tree sampling
OUT_DIR="${OUT_ROOT}/run_rtree_sampling_binary"
mkdir -p "${OUT_DIR}"
run_one "${SJS_RUN}" "${LOG_DIR}" "run_rtree_sampling_binary" \
  --dataset_source=binary --dataset="${DATASET_NAME}_binary" \
  --path_r="${R_BIN}" --path_s="${S_BIN}" \
  --method=r_tree --variant=sampling --t="${T}" --seed="${SEED}" --repeats="${REPEATS}" \
  --out_dir="${OUT_DIR}" --write_samples=1 \
  --threads="${THREADS}"
[[ -f "${OUT_DIR}/run.csv" ]] || die "Missing run.csv for ${OUT_DIR}"

# --------------------------
# 5) correctness & sample-quality sanity via sjs_verify
# --------------------------
log "Step 5/6: Correctness & sample-quality sanity via sjs_verify"

METHODS=(ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection tsunami)
for m in "${METHODS[@]}"; do
  verify_one "${SJS_VERIFY}" "${LOG_DIR}" "${m}" "sampling"
done

# --------------------------
# 6) optional: ASAN/UBSAN build + tests (+ verify)
# --------------------------
if [[ "${RUN_ASAN}" == "1" ]]; then
  log "Step 6/6: ASAN/UBSAN pass (enabled)"

  BUILD_DIR_ASAN="${ROOT}/build_exp0_asan"
  LOG_DIR_ASAN="${LOG_DIR}/asan"
  mkdir -p "${LOG_DIR_ASAN}"

  rm -rf "${BUILD_DIR_ASAN}"
  mkdir -p "${BUILD_DIR_ASAN}"
  pushd "${BUILD_DIR_ASAN}" >/dev/null

  cmake "${ROOT}" -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DSJS_BUILD_ROOT_APPS=ON -DSJS_BUILD_TESTS=ON \
    -DCMAKE_CXX_FLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer" \
    -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined" \
    2>&1 | tee "${LOG_DIR_ASAN}/cmake_configure_asan.log"

  cmake --build . -j "${ASAN_JOBS}" 2>&1 | tee "${LOG_DIR_ASAN}/cmake_build_asan.log"

  export ASAN_OPTIONS="${ASAN_OPTIONS:-halt_on_error=1:abort_on_error=1:detect_leaks=1:strict_string_checks=1:symbolize=1}"
  export UBSAN_OPTIONS="${UBSAN_OPTIONS:-halt_on_error=1:abort_on_error=1:print_stacktrace=1}"

  ctest --output-on-failure 2>&1 | tee "${LOG_DIR_ASAN}/ctest_asan.log"

  SJS_VERIFY_ASAN="${BUILD_DIR_ASAN}/sjs_verify"
  [[ -x "${SJS_VERIFY_ASAN}" ]] || die "Missing executable: ${SJS_VERIFY_ASAN}"

  if [[ "${ASAN_VERIFY_ALL}" == "1" ]]; then
    log "  - ASAN verify: ALL methods"
    for m in "${METHODS[@]}"; do
      verify_one "${SJS_VERIFY_ASAN}" "${LOG_DIR_ASAN}" "${m}" "sampling"
    done
  else
    log "  - ASAN verify: ours + r_tree (set ASAN_VERIFY_ALL=1 for all methods)"
    verify_one "${SJS_VERIFY_ASAN}" "${LOG_DIR_ASAN}" "ours" "sampling"
    verify_one "${SJS_VERIFY_ASAN}" "${LOG_DIR_ASAN}" "r_tree" "sampling"
  fi

  popd >/dev/null
else
  log "Step 6/6: ASAN/UBSAN pass (skipped, RUN_ASAN=0)"
fi

log "EXP-0 PASSED ✅"
log "Artifacts:"
log "  - Logs:    ${LOG_DIR}"
log "  - Results: ${OUT_ROOT}"
log "  - Data:    ${DATA_DIR}"
if [[ "${RUN_ASAN}" == "1" ]]; then
  log "  - ASAN build: ${ROOT}/build_exp0_asan"
  log "  - ASAN logs:  ${LOG_DIR}/asan"
fi
