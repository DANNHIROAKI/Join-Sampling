#!/usr/bin/env bash
# run/run_exp2.sh (optimized & aligned)
#
# EXP-2: Runtime vs t (RQ2)
# 目标：在“固定数据 + 固定 α”下，仅扫描 t，测端到端 wall-clock，统计 median/p95。
# 语义约束：二维 half-open 盒；从 J 中均匀 i.i.d. 有放回采样；不计数据生成/IO 等非主体时间。
# 设计依据：见 EXP-2.md（固定项/扫描项/三种模式/计时口径/统计方式/输出结构/三张图）。
#
# 统一规范（满足你的 5 条硬要求）：
#  1) 与 EXP-2.md 表述/原理/思想严格对齐（仅扫 t；single-thread；端到端 wall；用 sweep_t.json）。
#  2) Build 目录统一在仓库根的 build/ 下：默认 ${ROOT}/build（Release）。
#  3) 实验结果统一落在 ${ROOT}/results/raw/exp2；每次运行覆盖旧结果。
#  4) Bash 中不嵌 Python；绘图脚本放在 run/include/plot_exp2.py（若不存在则回退 run/plot_exp2.py）。
#  5) 运行过程产出的 json/csv/log 等“必要文件”全部先落在 run/temp/exp2，再整体同步到 results/raw/exp2。
#
# 用法：
#   bash run/run_exp2.sh
# 可选参数：
#   --config <path>   指定 sweep 配置（默认 config/sweeps/sweep_t.json，只扫 t）
#   --no-build        跳过构建（假设可执行已在 build/ 下）
#   --no-plot         跳过绘图
#   --clean           先清空 run/temp/exp2 和 results/raw/exp2
#   --build-type {Release|Debug|RelWithDebInfo|MinSizeRel}  (默认 Release)
#   -h, --help        显示帮助
#
set -euo pipefail
IFS=$'\n\t'

# ----------------------------
# helpers
# ----------------------------
die() { echo "[EXP2][FATAL] $*" >&2; exit 1; }
log() { echo "[EXP2] $*"; }

nproc_safe() {
  if command -v nproc >/dev/null 2>&1; then nproc
  elif [[ "$(uname -s)" == "Darwin" ]] && command -v sysctl >/dev/null 2>&1; then sysctl -n hw.ncpu
  else echo 4; fi
}

abs_path() {
  python3 - "$1" <<'PY' || true
import os,sys; print(os.path.abspath(sys.argv[1]))
PY
}

find_exe() {
  local build_dir="$1"; local name="$2"
  local cands=(
    "${build_dir}/${name}"
    "${build_dir}/apps/${name}"
  )
  for p in "${cands[@]}"; do
    [[ -x "$p" ]] && { echo "$p"; return; }
  done
  # fallback: search
  local p="$(command -v "${build_dir}/${name}" 2>/dev/null || true)"
  [[ -n "$p" && -x "$p" ]] && { echo "$p"; return; }
  p="$(find "${build_dir}" -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "$p" ]] && { echo "$p"; return; }
  die "Cannot find executable '${name}' under ${build_dir}."
}

usage() {
  cat <<'EOF'
EXP-2 runner (Runtime vs t, aligned)

Usage:
  bash run/run_exp2.sh [options]

Options:
  --config <path>     Sweep JSON (default: config/sweeps/sweep_t.json)
  --no-build          Skip build step
  --no-plot           Skip plotting
  --clean             Remove temp/results dirs before run
  --build-type <t>    Release|Debug|RelWithDebInfo|MinSizeRel (default: Release)
  -h, --help          Show this help
EOF
}

# ----------------------------
# resolve paths
# ----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

CONFIG="${ROOT}/config/sweeps/sweep_t.json"
DO_BUILD=1
DO_PLOT=1
DO_CLEAN=0
BUILD_TYPE="Release"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2;;
    --no-build) DO_BUILD=0; shift;;
    --no-plot) DO_PLOT=0; shift;;
    --clean) DO_CLEAN=1; shift;;
    --build-type) BUILD_TYPE="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1 (try --help)";;
  esac
done

BUILD_DIR="${ROOT}/build"
TEMP_DIR="${ROOT}/run/temp/exp2"
RESULT_DIR="${ROOT}/results/raw/exp2"
LOG_DIR="${TEMP_DIR}/logs"
META_DIR="${TEMP_DIR}/meta"

mkdir -p "${TEMP_DIR}" "${RESULT_DIR}" "${LOG_DIR}" "${META_DIR}"

# paths to python plotting (prefer include/, fallback to legacy run/)
PLOT_SCRIPT="${ROOT}/run/include/plot_exp2.py"
[[ -f "${PLOT_SCRIPT}" ]] || PLOT_SCRIPT="${ROOT}/run/plot_exp2.py"

# absolute paths (if python3 available)
if command -v python3 >/dev/null 2>&1; then
  CONFIG="$(abs_path "${CONFIG}")"
fi
[[ -f "${CONFIG}" ]] || die "Sweep config not found: ${CONFIG}"

# ----------------------------
# single-thread fairness
# ----------------------------
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# ----------------------------
# clean temp/results (overwrite semantics)
# ----------------------------
if [[ "${DO_CLEAN}" -eq 1 ]]; then
  log "Cleaning ${TEMP_DIR} and ${RESULT_DIR}"
  rm -rf "${TEMP_DIR}" "${RESULT_DIR}"
  mkdir -p "${TEMP_DIR}" "${RESULT_DIR}" "${LOG_DIR}" "${META_DIR}"
fi

# always start from empty temp to avoid stale files
rm -rf "${TEMP_DIR}"
mkdir -p "${TEMP_DIR}" "${LOG_DIR}" "${META_DIR}"

# ----------------------------
# build under ROOT/build
# ----------------------------
if [[ "${DO_BUILD}" -eq 1 ]]; then
  command -v cmake >/dev/null 2>&1 || die "cmake not found"
  log "Configuring CMake (type=${BUILD_TYPE}) at ${BUILD_DIR}"
  cmake -S "${ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
    2>&1 | tee "${LOG_DIR}/cmake_configure.log"
  log "Building ..."
  cmake --build "${BUILD_DIR}" -j "$(nproc_safe)" \
    2>&1 | tee "${LOG_DIR}/cmake_build.log"
else
  log "Skip build (--no-build)"
fi

SJS_SWEEP="$(find_exe "${BUILD_DIR}" "sjs_sweep")"
log "Using sjs_sweep: ${SJS_SWEEP}"

# ----------------------------
# record environment (manifest)
# ----------------------------
{
  echo "exp=EXP-2 (Runtime vs t)"
  echo "timestamp=$(date -Is || date)"
  echo "repo_root=${ROOT}"
  echo "build_dir=${BUILD_DIR}"
  echo "build_type=${BUILD_TYPE}"
  echo "config=${CONFIG}"
  echo "threads=1 (OMP/BLAS capped)"
  echo
  echo "uname:"; uname -a || true; echo
  echo "compiler:"; (c++ --version || g++ --version || clang++ --version) 2>/dev/null || true; echo
  echo "cmake:"; cmake --version 2>/dev/null || true; echo
  if command -v git >/dev/null 2>&1 && [[ -d "${ROOT}/.git" ]]; then
    echo "git_rev=$(git -C "${ROOT}" rev-parse HEAD || true)"
    echo "git_status:"; git -C "${ROOT}" status --porcelain || true
  fi
} > "${META_DIR}/manifest.txt"

# also keep a commands log
COMMANDS_LOG="${TEMP_DIR}/commands.log"
: > "${COMMANDS_LOG}"

run_and_log() {
  local desc="$1"; shift
  echo "[CMD] ${desc}: $*" | tee -a "${COMMANDS_LOG}"
  "$@"
}

# ----------------------------
# run sweep (ALL methods×variants×t×repeats from config)
# ----------------------------
log "Running sweep ..."
pushd "${ROOT}" >/dev/null
set +e
run_and_log "sjs_sweep" "${SJS_SWEEP}" --config "${CONFIG}" --out_dir "${TEMP_DIR}" \
  2>&1 | tee "${LOG_DIR}/sjs_sweep.log"
rc="${PIPESTATUS[0]}"
set -e
popd >/dev/null

[[ "${rc}" -eq 0 ]] || die "sjs_sweep failed (rc=${rc}), see ${LOG_DIR}/sjs_sweep.log"

# required outputs
[[ -f "${TEMP_DIR}/sweep_raw.csv" ]] || die "Missing ${TEMP_DIR}/sweep_raw.csv"
[[ -f "${TEMP_DIR}/sweep_summary.csv" ]] || die "Missing ${TEMP_DIR}/sweep_summary.csv"

log "Sweep OK. raw+summary CSV in ${TEMP_DIR}"

# ----------------------------
# plotting (optional; no embedded python)
# ----------------------------
if [[ "${DO_PLOT}" -eq 1 ]]; then
  if command -v python3 >/dev/null 2>&1 && [[ -f "${PLOT_SCRIPT}" ]]; then
    log "Plotting via: ${PLOT_SCRIPT}"
    # 建议参数：t0=1000，误差条用 p95（EXP-2.md 推荐三张图 + ns/sample 估计）
    python3 "${PLOT_SCRIPT}" --out_dir "${TEMP_DIR}" --t0 1000 --error p95 \
      2>&1 | tee "${LOG_DIR}/plot_exp2.log" || \
      log "Plot script returned non-zero; continue without plots."
  else
    log "Skip plotting (python3 or plot script not found)."
  fi
else
  log "Plot step disabled (--no-plot)."
fi

# ----------------------------
# sync to results/raw/exp2 (overwrite old)
# ----------------------------
log "Syncing to ${RESULT_DIR} (overwrite old)"
rm -rf "${RESULT_DIR}"
mkdir -p "${RESULT_DIR}"
# 保留 temp 目录内部的层级与文件名
rsync -a --delete "${TEMP_DIR}/" "${RESULT_DIR}/"

log "DONE. Results in: ${RESULT_DIR}"
