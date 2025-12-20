#!/usr/bin/env bash
# run/run_exp0.sh
#
# EXP-0：基础 sanity（单元级 + 端到端 smoke）
# - 编译（Release/Debug）
# - 跑 ctest（几何/事件扫描/Join oracle/采样质量/所有 baselines 的 smoke）
# - 端到端：生成小数据集（写 bin/csv）→ sjs_run（synthetic & binary）
# - 端到端：sjs_verify（oracle + 质量指标），并强制检查 missing_in_universe=0
#
# 用法（仓库根目录或任意目录均可）：
#   bash run/run_exp0.sh
#
# 可选环境变量（不设则用默认值）：
#   BUILD_TYPE=Release|Debug
#   JOBS=8
#   THREADS=1
#   NR=2000  NS=2000  ALPHA=1e-3  GEN_SEED=1
#   T=20000  SEED=1
#
set -euo pipefail
IFS=$'\n\t'

# --------------------------
# helpers
# --------------------------
log()  { echo -e "[EXP0] $*"; }
die()  { echo -e "[EXP0][FATAL] $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"
}

# Resolve repo root (script lives in <root>/run/run_exp0.sh)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Parameters (override via env)
BUILD_TYPE="${BUILD_TYPE:-Release}"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)}"
THREADS="${THREADS:-1}"

NR="${NR:-2000}"
NS="${NS:-2000}"
ALPHA="${ALPHA:-1e-3}"
GEN_SEED="${GEN_SEED:-1}"

T="${T:-20000}"
SEED="${SEED:-1}"

# Paths
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

log "Repo root: ${ROOT}"
log "Build dir: ${BUILD_DIR}"
log "Output dir: ${OUT_ROOT}"
log "Data dir: ${DATA_DIR}"

log "Params: BUILD_TYPE=${BUILD_TYPE} JOBS=${JOBS} THREADS=${THREADS}"
log "Params: NR=${NR} NS=${NS} ALPHA=${ALPHA} GEN_SEED=${GEN_SEED} T=${T} SEED=${SEED}"

# record environment
{
  echo "date: $(date -Is || date)"
  echo "uname: $(uname -a || true)"
  echo "cmake: $(cmake --version | head -n1 || true)"
  echo "ctest: $(ctest --version | head -n1 || true)"
} > "${LOG_DIR}/env.txt"

# --------------------------
# 1) clean build + compile
# --------------------------
log "Step 1/5: Clean build + compile"
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
log "Step 2/5: Run unit tests (ctest)"
pushd "${BUILD_DIR}" >/dev/null
ctest --output-on-failure 2>&1 | tee "${LOG_DIR}/ctest.log"
popd >/dev/null

# --------------------------
# 3) end-to-end: generate small dataset to disk (binary + csv)
# --------------------------
log "Step 3/5: Generate a small dataset (binary + csv) via sjs_gen_dataset"
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
# 4) end-to-end: sjs_run smoke (synthetic + binary)
# --------------------------
log "Step 4/5: End-to-end smoke via sjs_run"

run_one() {
  local tag="$1"; shift
  local log_file="${LOG_DIR}/${tag}.log"
  log "  - Running: ${tag}"
  "${SJS_RUN}" "$@" 2>&1 | tee "${log_file}"
}

# (A) Synthetic on-the-fly, run all 3 variants for a representative method (ours)
for VARIANT in sampling enum_sampling adaptive; do
  OUT_DIR="${OUT_ROOT}/run_ours_${VARIANT}_synthetic"
  mkdir -p "${OUT_DIR}"
  run_one "run_ours_${VARIANT}_synthetic" \
    --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}_synthetic" \
    --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
    --method=ours --variant="${VARIANT}" --t="${T}" --seed="${SEED}" --repeats=1 \
    --out_dir="${OUT_DIR}" --write_samples=1 \
    --threads="${THREADS}"

  [[ -f "${OUT_DIR}/run.csv" ]] || die "Missing run.csv for ${OUT_DIR}"
done

# (B) Binary input path (tests IO), just sampling variant
OUT_DIR="${OUT_ROOT}/run_ours_sampling_binary"
mkdir -p "${OUT_DIR}"
run_one "run_ours_sampling_binary" \
  --dataset_source=binary --dataset="${DATASET_NAME}_binary" \
  --path_r="${R_BIN}" --path_s="${S_BIN}" \
  --method=ours --variant=sampling --t="${T}" --seed="${SEED}" --repeats=1 \
  --out_dir="${OUT_DIR}" --write_samples=1 \
  --threads="${THREADS}"

[[ -f "${OUT_DIR}/run.csv" ]] || die "Missing run.csv for ${OUT_DIR}"

# --------------------------
# 5) end-to-end: sjs_verify (oracle + quality); enforce missing_in_universe=0
# --------------------------
log "Step 5/5: Correctness & sample-quality sanity via sjs_verify"

verify_one() {
  local method="$1"
  local variant="$2"

  local tag="verify_${method}_${variant}"
  local log_file="${LOG_DIR}/${tag}.log"

  log "  - Verifying: method=${method} variant=${variant}"
  "${SJS_VERIFY}" \
    --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}_verify" \
    --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
    --method="${method}" --variant="${variant}" --t="${T}" --seed="${SEED}" --repeats=1 \
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
}

# Verify all registered methods in sampling mode (core claim)
METHODS=(ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection tsunami)
for m in "${METHODS[@]}"; do
  verify_one "${m}" "sampling"
done

log "EXP-0 PASSED ✅"
log "Artifacts:"
log "  - Logs:    ${LOG_DIR}"
log "  - Results: ${OUT_ROOT}"
log "  - Data:    ${DATA_DIR}"
