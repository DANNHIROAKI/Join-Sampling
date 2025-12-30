#!/usr/bin/env bash
# run/run_exp7.sh
#
# EXP-7：阶段分解与解释性 Profiling（Phase Breakdown）
# ----------------------------------------------------
# 修正点（重要）：
#   (A) 默认不再用 --exclude_enum_truncated 1 全局过滤，
#       因为 adaptive/fallback 的 pilot 截断（enum_truncated==1）是“预期行为”，
#       EXP-7 需要展示它来解释 dense regime 的机制。
#       -> 新增 EXCLUDE_ENUM_TRUNCATED（默认 0），用于控制后处理过滤。
#
#   (B) 默认 ALPHAS 改为 "0.1 3 5 30"（与 EXP-7.md 的“推荐增强”一致，
#       在阈值附近加点更能解释 adaptive 拐点）。
#
# 其余目录/规范保持与你现有脚本一致。
#
# 用法：
#   bash run/run_exp7.sh
#
# 常用覆盖参数（环境变量）：
#   BUILD_TYPE=Release|Debug|RelWithDebInfo|MinSizeRel
#   CLEAN_BUILD=0|1
#   NR=100000  NS=100000
#   T=100000   REPEATS=3
#   GEN_SEED=1 SEED=1
#   THREADS=1
#   ALPHAS="0.1 3 5 30"
#   METHODS="ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection tsunami"
#   J_STAR=1000000
#   ENUM_CAP=0
#   INCLUDE_ENUM_SPARSE=1  SPARSE_ALPHA=0.1
#   EXCLUDE_ENUM_TRUNCATED=0   # ✅关键：默认 0，避免误删 adaptive fallback
#
set -euo pipefail
IFS=$'\n\t'

trap 'echo -e "[EXP7][FATAL] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2; exit 1' ERR

# --------------------------
# helpers
# --------------------------
log() { echo -e "[EXP7] $*"; }
die() { echo -e "[EXP7][FATAL] $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"
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

find_exe() {
  local build_dir="$1"
  local name="$2"

  # Fast paths
  if [[ -x "${build_dir}/${name}" ]]; then
    echo "${build_dir}/${name}"; return 0
  fi
  if [[ -x "${build_dir}/apps/${name}" ]]; then
    echo "${build_dir}/apps/${name}"; return 0
  fi

  # Fallback search
  local p
  p="$(find "${build_dir}" -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "${p}" ]] || return 1
  echo "${p}"
}

# --------------------------
# resolve paths
# --------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# --------------------------
# parameters (aligned with EXP-7.md; with fixes)
# --------------------------
BUILD_TYPE="${EXP7_BUILD_TYPE:-${BUILD_TYPE:-Release}}"
CLEAN_BUILD="${EXP7_CLEAN_BUILD:-${CLEAN_BUILD:-0}}"

NR="${NR:-100000}"
NS="${NS:-100000}"
T="${T:-100000}"
REPEATS="${REPEATS:-3}"
GEN_SEED="${GEN_SEED:-1}"
SEED="${SEED:-1}"
THREADS="${THREADS:-1}"

# ✅改动：默认加入阈值附近点 5（更容易解释 adaptive 拐点）
ALPHAS="${ALPHAS:-0.1 3 5 30}"

METHODS="${METHODS:-ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection tsunami}"

J_STAR="${J_STAR:-1000000}"
ENUM_CAP="${ENUM_CAP:-0}"

INCLUDE_ENUM_SPARSE="${INCLUDE_ENUM_SPARSE:-1}"
SPARSE_ALPHA="${SPARSE_ALPHA:-0.1}"

# ✅关键改动：默认不排除 enum_truncated（避免误删 adaptive/fallback 的 pilot 截断）
EXCLUDE_ENUM_TRUNCATED="${EXP7_EXCLUDE_ENUM_TRUNCATED:-${EXCLUDE_ENUM_TRUNCATED:-0}}"

# --------------------------
# fixed directories (per your global requirements)
# --------------------------
BUILD_SUBDIR="$(build_subdir_from_type "${BUILD_TYPE}")"
BUILD_DIR="${REPO_ROOT}/build/${BUILD_SUBDIR}"

TEMP_DIR="${REPO_ROOT}/run/temp/exp7"
RESULT_DIR="${REPO_ROOT}/results/raw/exp7"

PY_MAKE_CFG="${REPO_ROOT}/run/include/exp7_make_sweep_config.py"
PY_POST="${REPO_ROOT}/run/include/exp7_postprocess_phase_breakdown.py"

# --------------------------
# preflight
# --------------------------
need_cmd cmake
need_cmd python3
need_cmd tee
need_cmd find

[[ -f "${REPO_ROOT}/CMakeLists.txt" ]] || die "CMakeLists.txt not found at repo root: ${REPO_ROOT}"
[[ -f "${PY_MAKE_CFG}" ]] || die "Missing helper: ${PY_MAKE_CFG}"
[[ -f "${PY_POST}" ]] || die "Missing helper: ${PY_POST}"

# numeric sanity checks
for v in NR NS T REPEATS GEN_SEED SEED THREADS J_STAR ENUM_CAP EXCLUDE_ENUM_TRUNCATED; do
  val="${!v}"
  if ! [[ "${val}" =~ ^[0-9]+$ ]]; then
    die "${v} must be an integer (got: '${val}')"
  fi
  if [[ "${val}" -lt 0 ]]; then
    die "${v} must be >= 0 (got: '${val}')"
  fi
done
if [[ "${THREADS}" -le 0 ]]; then
  die "THREADS must be >= 1 (got: '${THREADS}')"
fi
if [[ "${EXCLUDE_ENUM_TRUNCATED}" -ne 0 && "${EXCLUDE_ENUM_TRUNCATED}" -ne 1 ]]; then
  die "EXCLUDE_ENUM_TRUNCATED must be 0 or 1 (got: '${EXCLUDE_ENUM_TRUNCATED}')"
fi

# clean temp (always overwrite temp)
rm -rf "${TEMP_DIR}"
mkdir -p "${TEMP_DIR}/logs" "${TEMP_DIR}/configs"

TS="$(date +%Y%m%d_%H%M%S)"

log "Repo root    : ${REPO_ROOT}"
log "Build type   : ${BUILD_TYPE}"
log "Build dir    : ${BUILD_DIR}"
log "Temp dir     : ${TEMP_DIR}"
log "Final results: ${RESULT_DIR} (will be overwritten on success)"
log "N(R,S)       : ${NR}, ${NS}"
log "t, repeats   : ${T}, ${REPEATS}"
log "alphas       : ${ALPHAS}"
log "methods      : ${METHODS}"
log "threads      : ${THREADS}"
log "J_STAR       : ${J_STAR}"
log "ENUM_CAP     : ${ENUM_CAP}"
log "enum sparse  : INCLUDE_ENUM_SPARSE=${INCLUDE_ENUM_SPARSE} SPARSE_ALPHA=${SPARSE_ALPHA}"
log "postprocess  : EXCLUDE_ENUM_TRUNCATED=${EXCLUDE_ENUM_TRUNCATED} (0 keeps adaptive fallback pilot; 1 may hide it)"

if [[ "${INCLUDE_ENUM_SPARSE}" == "1" ]]; then
  if ! echo " ${ALPHAS} " | grep -q " ${SPARSE_ALPHA} "; then
    log "[WARN] SPARSE_ALPHA=${SPARSE_ALPHA} is not in ALPHAS='${ALPHAS}'. enum_sparse will run at an extra alpha point."
  fi
fi

# Thread caps (fairness + reproducibility)
export OMP_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"

# Save manifest + env
{
  echo "EXP-7 manifest"
  echo "timestamp=${TS}"
  echo "repo_root=${REPO_ROOT}"
  echo "build_type=${BUILD_TYPE}"
  echo "build_dir=${BUILD_DIR}"
  echo "temp_dir=${TEMP_DIR}"
  echo "result_dir=${RESULT_DIR}"
  echo
  echo "NR=${NR}"
  echo "NS=${NS}"
  echo "T=${T}"
  echo "REPEATS=${REPEATS}"
  echo "GEN_SEED=${GEN_SEED}"
  echo "SEED=${SEED}"
  echo "THREADS=${THREADS}"
  echo
  echo "ALPHAS=${ALPHAS}"
  echo "METHODS=${METHODS}"
  echo "J_STAR=${J_STAR}"
  echo "ENUM_CAP=${ENUM_CAP}"
  echo "INCLUDE_ENUM_SPARSE=${INCLUDE_ENUM_SPARSE}"
  echo "SPARSE_ALPHA=${SPARSE_ALPHA}"
  echo
  echo "EXCLUDE_ENUM_TRUNCATED=${EXCLUDE_ENUM_TRUNCATED}"
} > "${TEMP_DIR}/MANIFEST.txt"

{
  echo "timestamp=${TS}"
  uname -a || true
  echo
  echo "cmake:"; cmake --version 2>/dev/null || true
  echo
  echo "compiler:"; (c++ --version || g++ --version || clang++ --version) 2>/dev/null || true
  echo
  if command -v git >/dev/null 2>&1 && [[ -d "${REPO_ROOT}/.git" ]]; then
    echo "git:";
    (cd "${REPO_ROOT}" && git rev-parse HEAD && git status --porcelain) || true
  fi
} > "${TEMP_DIR}/logs/env.txt"

# --------------------------
# build
# --------------------------
log "== Build =="
if [[ "${CLEAN_BUILD}" == "1" ]]; then
  log "Cleaning build dir: ${BUILD_DIR}"
  rm -rf "${BUILD_DIR}"
fi
mkdir -p "${BUILD_DIR}"

cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DSJS_BUILD_ROOT_APPS=ON \
  -DSJS_BUILD_TESTS=OFF \
  2>&1 | tee "${TEMP_DIR}/logs/cmake_configure.log"

cmake --build "${BUILD_DIR}" -j "$(nproc_safe)" \
  2>&1 | tee "${TEMP_DIR}/logs/cmake_build.log"

SJS_SWEEP="$(find_exe "${BUILD_DIR}" sjs_sweep || true)"
[[ -n "${SJS_SWEEP}" ]] || die "Could not locate sjs_sweep under ${BUILD_DIR}"
log "Using sjs_sweep: ${SJS_SWEEP}"

# --------------------------
# generate sweep configs (JSON in run/temp)
# --------------------------
OUT_SA="${TEMP_DIR}/sampling_adaptive"
OUT_ENUM="${TEMP_DIR}/enum_sparse"

CFG_SA="${TEMP_DIR}/configs/exp7_sampling_adaptive.json"
CFG_ENUM="${TEMP_DIR}/configs/exp7_enum_sparse.json"

log "== Generate sweep config: sampling + adaptive (all ALPHAS) =="
python3 "${PY_MAKE_CFG}" \
  --out "${CFG_SA}" \
  --out_dir "${OUT_SA}" \
  --run_tag "exp7_sampling_adaptive" \
  --note "EXP-7 (phase breakdown): sampling+adaptive on representative alphas." \
  --nr "${NR}" --ns "${NS}" --t "${T}" --repeats "${REPEATS}" \
  --gen_seed "${GEN_SEED}" --seed "${SEED}" \
  --threads "${THREADS}" --j_star "${J_STAR}" --enum_cap "${ENUM_CAP}" \
  --alphas "${ALPHAS}" \
  --methods "${METHODS}" \
  --variants "sampling adaptive" \
  2>&1 | tee "${TEMP_DIR}/logs/make_cfg_sampling_adaptive.log"

log "== Run sweep: sampling + adaptive =="
mkdir -p "${OUT_SA}"
cp -f "${CFG_SA}" "${OUT_SA}/sweep_config.json"

pushd "${REPO_ROOT}" >/dev/null
set +e
"${SJS_SWEEP}" --config="${CFG_SA}" 2>&1 | tee "${TEMP_DIR}/logs/sjs_sweep_sampling_adaptive.log"
rc="${PIPESTATUS[0]}"
set -e
popd >/dev/null
[[ "${rc}" -eq 0 ]] || die "sjs_sweep (sampling+adaptive) failed with exit code ${rc}"

RAW_SA="${OUT_SA}/sweep_raw.csv"
SUM_SA="${OUT_SA}/sweep_summary.csv"
[[ -f "${RAW_SA}" ]] || die "Missing expected output: ${RAW_SA}"
[[ -f "${SUM_SA}" ]] || die "Missing expected output: ${SUM_SA}"
log "Wrote: ${RAW_SA}"
log "Wrote: ${SUM_SA}"

RAW_ENUM=""  # optional
if [[ "${INCLUDE_ENUM_SPARSE}" == "1" ]]; then
  log "== Generate sweep config: enum_sampling (sparse only) =="
  python3 "${PY_MAKE_CFG}" \
    --out "${CFG_ENUM}" \
    --out_dir "${OUT_ENUM}" \
    --run_tag "exp7_enum_sparse" \
    --note "EXP-7 (phase breakdown): enum_sampling only at sparse alpha to keep runtime bounded." \
    --nr "${NR}" --ns "${NS}" --t "${T}" --repeats "${REPEATS}" \
    --gen_seed "${GEN_SEED}" --seed "${SEED}" \
    --threads "${THREADS}" --j_star "${J_STAR}" --enum_cap "${ENUM_CAP}" \
    --alphas "${SPARSE_ALPHA}" \
    --methods "${METHODS}" \
    --variants "enum_sampling" \
    2>&1 | tee "${TEMP_DIR}/logs/make_cfg_enum_sparse.log"

  log "== Run sweep: enum_sampling (sparse only) =="
  mkdir -p "${OUT_ENUM}"
  cp -f "${CFG_ENUM}" "${OUT_ENUM}/sweep_config.json"

  pushd "${REPO_ROOT}" >/dev/null
  set +e
  "${SJS_SWEEP}" --config="${CFG_ENUM}" 2>&1 | tee "${TEMP_DIR}/logs/sjs_sweep_enum_sparse.log"
  rc="${PIPESTATUS[0]}"
  set -e
  popd >/dev/null
  [[ "${rc}" -eq 0 ]] || die "sjs_sweep (enum_sparse) failed with exit code ${rc}"

  RAW_ENUM="${OUT_ENUM}/sweep_raw.csv"
  SUM_ENUM="${OUT_ENUM}/sweep_summary.csv"
  [[ -f "${RAW_ENUM}" ]] || die "Missing expected output: ${RAW_ENUM}"
  [[ -f "${SUM_ENUM}" ]] || die "Missing expected output: ${SUM_ENUM}"
  log "Wrote: ${RAW_ENUM}"
  log "Wrote: ${SUM_ENUM}"
else
  log "Skipping enum_sampling sweep (INCLUDE_ENUM_SPARSE=${INCLUDE_ENUM_SPARSE})"
fi

# --------------------------
# postprocess
# --------------------------
log "== Post-process: phase breakdown table + plots =="
OUT_BREAKDOWN="${TEMP_DIR}/exp7_breakdown_median.csv"
OUT_MD="${TEMP_DIR}/exp7_breakdown_median.md"
FIG_DIR="${TEMP_DIR}/figs"
MERGED_RAW="${TEMP_DIR}/exp7_merged_sweep_raw.csv"

python3 "${PY_POST}" \
  --raw_sa "${RAW_SA}" \
  --raw_enum "${RAW_ENUM}" \
  --out_breakdown "${OUT_BREAKDOWN}" \
  --out_md "${OUT_MD}" \
  --fig_dir "${FIG_DIR}" \
  --merged_raw "${MERGED_RAW}" \
  --exclude_enum_truncated "${EXCLUDE_ENUM_TRUNCATED}" \
  2>&1 | tee "${TEMP_DIR}/logs/postprocess.log"

[[ -f "${OUT_BREAKDOWN}" ]] || die "Missing expected output: ${OUT_BREAKDOWN}"
[[ -f "${OUT_MD}" ]] || die "Missing expected output: ${OUT_MD}"
[[ -f "${MERGED_RAW}" ]] || die "Missing expected output: ${MERGED_RAW}"

# --------------------------
# sync temp -> results/raw/exp7 (overwrite)
# --------------------------
log "== Sync results to ${RESULT_DIR} (overwrite) =="
rm -rf "${RESULT_DIR}"
mkdir -p "${RESULT_DIR}"
cp -a "${TEMP_DIR}/." "${RESULT_DIR}/"

log "DONE ✅"
log "Temp dir      : ${TEMP_DIR}"
log "Final results : ${RESULT_DIR}"

cat <<EOF

EXP-7 outputs (final):
  - sampling_adaptive raw      : ${RESULT_DIR}/sampling_adaptive/sweep_raw.csv
  - sampling_adaptive summary  : ${RESULT_DIR}/sampling_adaptive/sweep_summary.csv
EOF

if [[ -n "${RAW_ENUM}" ]]; then
  cat <<EOF
  - enum_sparse raw            : ${RESULT_DIR}/enum_sparse/sweep_raw.csv
  - enum_sparse summary        : ${RESULT_DIR}/enum_sparse/sweep_summary.csv
EOF
fi

cat <<EOF
  - merged raw                 : ${RESULT_DIR}/exp7_merged_sweep_raw.csv
  - breakdown CSV              : ${RESULT_DIR}/exp7_breakdown_median.csv
  - breakdown MD               : ${RESULT_DIR}/exp7_breakdown_median.md
  - figures (if matplotlib)    : ${RESULT_DIR}/figs/exp7_phase_breakdown_alpha_*.png
EOF
