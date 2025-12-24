#!/usr/bin/env bash
# run/run_exp2.sh
#
# EXP-2：Runtime vs t（对应 RQ2）
#
# 目标：生成“可以直接放进论文”的 EXP-2 结果（稳定、可复现、图少且清晰）。
#
# 这份脚本在原始版本的基础上，重点补强了 3 个会导致“结果不够理想”的真实因素：
#
#  (A) 冷启动/热缓存偏差 → Δruntime 出现大量负值 / 曲线噪声大
#      - 新增 warmup repeats：默认多跑 1 次（rep=0）作为 warmup
#      - 绘图时默认排除 warmup（--warmup_reps）
#
#  (B) Δruntime 图被极端 outlier 压扁
#      - 新绘图脚本默认用 symlog-y（能同时显示小趋势与大 outlier）
#      - paper 模式默认只画 2 张图（runtime & Δruntime），每张图 3 个 panel
#
#  (C) 图太多、不“论文友好”
#      - 默认 --plot_mode=paper，只输出 2 张“主图级别”的 PDF/PNG
#      - 同时输出 ns/sample 的 CSV（用于论文文字/表格）
#
# 用法：
#   bash run/run_exp2.sh
#
# 常用可选参数：
#   --config <path>         Sweep JSON（默认：config/sweeps/sweep_t.json）
#   --build_type <type>     Release|Debug|RelWithDebInfo|MinSizeRel（默认：Release）
#   --clean                 删除 build/<type> 后重建
#   --no-build              跳过构建
#   --no-plot               跳过绘图
#   --threads <int>         覆盖 base.sys.threads（默认：1）
#   --write_samples <0|1>   覆盖 base.run.write_samples（默认：0）
#   --repeats <int>         *有效重复次数*（默认：10）
#   --warmup <int>          额外 warmup 次数（默认：1；绘图会自动排除）
#   --j_star <u64>          覆盖 base.run.j_star（默认：100000）
#
# 绘图相关：
#   --plot_mode <paper|full>    默认 paper（只生成 2 张主图）
#   --plot_t0 <int>             Δruntime 的 baseline t0（默认：10000）
#   --plot_error <auto|p95|stdev> 误差条策略（默认：auto）
#   --topk <int>                paper 图里仅绘制 top-k 最快方法 + ours（默认：6）
#   --exclude_methods <list>     逗号分隔；从图中排除某些方法（如 rejection）
#
# 输出目录：
#   - Temp（每次覆盖）：<repo_root>/run/temp/exp2/
#   - Final（成功后覆盖同步）：<repo_root>/results/raw/exp2/
#
# 注意：
#   - repeats 会被 patch 成 (repeats + warmup)，即 raw.csv 里会看到更多 rep。
#   - exp2_plot.py 会根据 --warmup_reps 自动忽略 rep < warmup。
#   - 这属于“实验口径修正”（排除冷启动），不改变算法语义。
#

set -euo pipefail
IFS=$'\n\t'

trap 'echo -e "[EXP2][FATAL] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2; exit 1' ERR

# ----------------------------
# Helpers
# ----------------------------
log() { echo -e "[EXP2] $*"; }
die() { echo -e "[EXP2][FATAL] $*" >&2; exit 1; }

need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

usage() {
  cat <<'EOF'
EXP-2 runner (Runtime vs t)

Usage:
  bash run/run_exp2.sh [options]

Options:
  --config <path>                 Sweep JSON (default: config/sweeps/sweep_t.json)
  --build_type <type>             Release|Debug|RelWithDebInfo|MinSizeRel (default: Release)
  --clean                         Remove build/<type> before building
  --no-build                      Skip build step
  --no-plot                       Skip plotting step
  --threads <int>                 Override base.sys.threads (default: 1)
  --write_samples <0|1>           Override base.run.write_samples (default: 0)
  --repeats <int>                 Effective repeats used for reporting/plots (default: 10)
  --warmup <int>                  Extra warmup repeats (excluded from plots) (default: 1)
  --j_star <u64>                  Override base.run.j_star (default: 100000)

Plot options:
  --plot_mode <paper|full>        Default: paper
  --plot_t0 <int>                 Baseline t0 for Δruntime (default: 10000)
  --plot_error <auto|p95|stdev>   Default: auto
  --topk <int>                    Paper plots: top-k fastest methods + ours (default: 6)
  --exclude_methods <csv>         Comma-separated methods to exclude in plots (default: empty)

  -h, --help                      Show help

Outputs (fixed):
  - Temp (always overwritten):      <repo_root>/run/temp/exp2/
  - Final (overwritten on success): <repo_root>/results/raw/exp2/
EOF
}

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

  if [[ -x "${build_dir}/${name}" ]]; then
    echo "${build_dir}/${name}"
    return
  fi
  if [[ -x "${build_dir}/apps/${name}" ]]; then
    echo "${build_dir}/apps/${name}"
    return
  fi

  local p
  p="$(find "${build_dir}" -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "${p}" ]] || die "Cannot find executable '${name}' under ${build_dir}. Did the build succeed?"
  echo "${p}"
}

to_bool_json() { [[ "$1" == "1" ]] && echo "true" || echo "false"; }

awk_col_idx() {
  local file="$1"
  local col="$2"
  awk -F',' -v c="${col}" 'NR==1{for(i=1;i<=NF;i++){if($i==c){print i; exit}}}' "${file}"
}

# ----------------------------
# Resolve paths (relative to this script)
# ----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ----------------------------
# Defaults (paper-friendly)
# ----------------------------
CONFIG="config/sweeps/sweep_t.json"
BUILD_TYPE="Release"
THREADS=1
WRITE_SAMPLES=0
REPEATS=10        # effective repeats (excluding warmup)
WARMUP=1          # extra repeats to warm cache; excluded from plots
J_STAR=100000

PLOT_MODE="paper"
PLOT_T0=10000
PLOT_ERROR="auto"
TOPK=6
EXCLUDE_METHODS=""

DO_CLEAN=0
DO_BUILD=1
DO_PLOT=1

# ----------------------------
# Parse args
# ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)          CONFIG="$2"; shift 2;;
    --build_type)      BUILD_TYPE="$2"; shift 2;;
    --threads)         THREADS="$2"; shift 2;;
    --write_samples)   WRITE_SAMPLES="$2"; shift 2;;
    --repeats)         REPEATS="$2"; shift 2;;
    --warmup)          WARMUP="$2"; shift 2;;
    --j_star)          J_STAR="$2"; shift 2;;
    --plot_mode)       PLOT_MODE="$2"; shift 2;;
    --plot_t0)         PLOT_T0="$2"; shift 2;;
    --plot_error)      PLOT_ERROR="$2"; shift 2;;
    --topk)            TOPK="$2"; shift 2;;
    --exclude_methods) EXCLUDE_METHODS="$2"; shift 2;;
    --clean)           DO_CLEAN=1; shift;;
    --no-build)        DO_BUILD=0; shift;;
    --no-plot)         DO_PLOT=0; shift;;
    -h|--help)         usage; exit 0;;
    *) die "Unknown argument: $1 (try --help)";;
  esac
done

if [[ "${CONFIG}" != /* ]]; then
  CONFIG="${REPO_ROOT}/${CONFIG}"
fi

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
need_cmd jq
need_cmd awk
need_cmd sort

[[ -f "${CONFIG}" ]] || die "Sweep config not found: ${CONFIG}"
[[ -f "${REPO_ROOT}/CMakeLists.txt" ]] || die "CMakeLists.txt not found at repo root: ${REPO_ROOT}"

if ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || [[ "${THREADS}" -le 0 ]]; then die "--threads must be a positive integer"; fi
if [[ "${WRITE_SAMPLES}" != "0" && "${WRITE_SAMPLES}" != "1" ]]; then die "--write_samples must be 0 or 1"; fi
if ! [[ "${REPEATS}" =~ ^[0-9]+$ ]] || [[ "${REPEATS}" -le 0 ]]; then die "--repeats must be a positive integer"; fi
if ! [[ "${WARMUP}" =~ ^[0-9]+$ ]] || [[ "${WARMUP}" -lt 0 ]]; then die "--warmup must be a non-negative integer"; fi
if ! [[ "${J_STAR}" =~ ^[0-9]+$ ]]; then die "--j_star must be a non-negative integer"; fi

if [[ "${PLOT_MODE}" != "paper" && "${PLOT_MODE}" != "full" ]]; then die "--plot_mode must be paper|full"; fi
if ! [[ "${PLOT_T0}" =~ ^[0-9]+$ ]] || [[ "${PLOT_T0}" -le 0 ]]; then die "--plot_t0 must be a positive integer"; fi
if [[ "${PLOT_ERROR}" != "auto" && "${PLOT_ERROR}" != "p95" && "${PLOT_ERROR}" != "stdev" ]]; then die "--plot_error must be auto|p95|stdev"; fi
if ! [[ "${TOPK}" =~ ^[0-9]+$ ]] || [[ "${TOPK}" -le 0 ]]; then die "--topk must be a positive integer"; fi

TOTAL_REPEATS=$((REPEATS + WARMUP))

# Clean temp (always overwrite temp).
rm -rf "${TEMP_DIR}"
mkdir -p "${TEMP_DIR}/logs"

TS="$(date +%Y%m%d_%H%M%S)"

log "Repo root     : ${REPO_ROOT}"
log "Config (in)   : ${CONFIG}"
log "Build type    : ${BUILD_TYPE}"
log "Build dir     : ${BUILD_DIR}"
log "Threads       : ${THREADS}"
log "write_samples : ${WRITE_SAMPLES}"
log "repeats(eff)  : ${REPEATS}"
log "warmup reps   : ${WARMUP}"
log "repeats(total): ${TOTAL_REPEATS}"
log "j_star        : ${J_STAR}"
log "plot_mode     : ${PLOT_MODE}"
log "plot_t0       : ${PLOT_T0}"
log "plot_error    : ${PLOT_ERROR}"
log "topk          : ${TOPK}"
log "exclude_methods: ${EXCLUDE_METHODS:-<none>}"
log "Temp dir      : ${TEMP_DIR}"
log "Final results : ${RESULT_DIR} (will be overwritten on success)"

# Force single-thread behavior for common libs.
export OMP_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export VECLIB_MAXIMUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"

# ----------------------------
# Patch sweep config (so overrides REALLY take effect)
# ----------------------------
ORIG_CFG_COPY="${TEMP_DIR}/sweep_t_original.json"
PATCHED_CFG="${TEMP_DIR}/sweep_t_used.json"

cp -f "${CONFIG}" "${ORIG_CFG_COPY}"

WRITE_SAMPLES_BOOL="$(to_bool_json "${WRITE_SAMPLES}")"

jq \
  --argjson repeats "${TOTAL_REPEATS}" \
  --argjson repeats_eff "${REPEATS}" \
  --argjson warmup_reps "${WARMUP}" \
  --argjson j_star "${J_STAR}" \
  --argjson threads "${THREADS}" \
  --argjson write_samples "${WRITE_SAMPLES_BOOL}" \
  --arg ts "${TS}" \
  '.base.run.repeats = $repeats
   | .base.run.j_star = $j_star
   | .base.sys.threads = $threads
   | .base.run.write_samples = $write_samples
   | .sweep.repeats = $repeats
   | .meta.patch = {
        "timestamp": $ts,
        "patched_by": "run/run_exp2.sh",
        "overrides": {
          "threads": $threads,
          "write_samples": $write_samples,
          "j_star": $j_star,
          "repeats_total": $repeats,
          "repeats_effective": $repeats_eff,
          "warmup_reps": $warmup_reps
        }
     }
  ' "${ORIG_CFG_COPY}" > "${PATCHED_CFG}"

# Save environment + manifest.
{
  echo "EXP-2 manifest"
  echo "timestamp=${TS}"
  echo "repo_root=${REPO_ROOT}"
  echo "config_in=${CONFIG}"
  echo "config_used=${PATCHED_CFG}"
  echo "build_type=${BUILD_TYPE}"
  echo "build_dir=${BUILD_DIR}"
  echo "threads=${THREADS}"
  echo "write_samples=${WRITE_SAMPLES}"
  echo "repeats_effective=${REPEATS}"
  echo "warmup_reps=${WARMUP}"
  echo "repeats_total=${TOTAL_REPEATS}"
  echo "j_star=${J_STAR}"
  echo "plot_mode=${PLOT_MODE}"
  echo "plot_t0=${PLOT_T0}"
  echo "plot_error=${PLOT_ERROR}"
  echo "topk=${TOPK}"
  echo "exclude_methods=${EXCLUDE_METHODS}"
  echo "temp_dir=${TEMP_DIR}"
  echo "result_dir=${RESULT_DIR}"
} > "${TEMP_DIR}/MANIFEST.txt"

{
  echo "timestamp=${TS}"
  echo "repo_root=${REPO_ROOT}"
  echo "config_in=${CONFIG}"
  echo "config_used=${PATCHED_CFG}"
  echo "build_dir=${BUILD_DIR}"
  echo "temp_dir=${TEMP_DIR}"
  echo "threads=${THREADS}"
  echo "write_samples=${WRITE_SAMPLES}"
  echo "repeats_effective=${REPEATS}"
  echo "warmup_reps=${WARMUP}"
  echo "repeats_total=${TOTAL_REPEATS}"
  echo "j_star=${J_STAR}"
  echo
  echo "uname:"; uname -a || true
  echo
  echo "compiler:"; (c++ --version || g++ --version || clang++ --version) 2>/dev/null || true
  echo
  echo "cmake:"; cmake --version 2>/dev/null || true
  echo
  echo "jq:"; jq --version 2>/dev/null || true
  echo
  if command -v git >/dev/null 2>&1 && [[ -d "${REPO_ROOT}/.git" ]]; then
    echo "git:"; (cd "${REPO_ROOT}" && git rev-parse HEAD && git status --porcelain) || true
  fi
} > "${TEMP_DIR}/logs/env.txt"

# ----------------------------
# Build
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
# Run sweep
# ----------------------------
OUT_DIR="${TEMP_DIR}"
RAW_FILE="sweep_raw.csv"
SUMMARY_FILE="sweep_summary.csv"

log "Running sweep (methods × variants × t × repeats_total=${TOTAL_REPEATS})..."
log "Config used: ${PATCHED_CFG}"

pushd "${REPO_ROOT}" >/dev/null
set +e
"${SJS_SWEEP}" \
  --config="${PATCHED_CFG}" \
  --out_dir="${OUT_DIR}" \
  --raw_file="${RAW_FILE}" \
  --summary_file="${SUMMARY_FILE}" \
  2>&1 | tee "${TEMP_DIR}/logs/sjs_sweep.log"
rc="${PIPESTATUS[0]}"
set -e
popd >/dev/null

if [[ "${rc}" -ne 0 ]]; then
  die "sjs_sweep failed with exit code ${rc}. See ${TEMP_DIR}/logs/sjs_sweep.log"
fi

[[ -f "${OUT_DIR}/${RAW_FILE}" ]] || die "Missing expected output: ${OUT_DIR}/${RAW_FILE}"
[[ -f "${OUT_DIR}/${SUMMARY_FILE}" ]] || die "Missing expected output: ${OUT_DIR}/${SUMMARY_FILE}"

log "Sweep finished OK."
log "Raw     : ${OUT_DIR}/${RAW_FILE}"
log "Summary : ${OUT_DIR}/${SUMMARY_FILE}"

# ----------------------------
# Post-checks (fast sanity)
# ----------------------------
SUM_PATH="${OUT_DIR}/${SUMMARY_FILE}"
RAW_PATH="${OUT_DIR}/${RAW_FILE}"

idx_ok="$(awk_col_idx "${SUM_PATH}" "ok_rate")"
idx_rep="$(awk_col_idx "${SUM_PATH}" "repeats")"
idx_cnt="$(awk_col_idx "${SUM_PATH}" "count_mean")"

if [[ -n "${idx_ok}" && -n "${idx_rep}" ]]; then
  uniq_repeats="$(awk -F',' -v okc="${idx_ok}" -v rc="${idx_rep}" 'NR>1 && $okc==1 {print $rc}' "${SUM_PATH}" | sort -u | tr '\n' ' ')"
  log "Summary ok points repeats values (includes warmup): ${uniq_repeats:-<none>}"
fi

if [[ -n "${idx_ok}" && -n "${idx_cnt}" ]]; then
  cnt_minmax="$(awk -F',' -v okc="${idx_ok}" -v cc="${idx_cnt}" 'NR>1 && $okc==1 {v=$cc; if(min==""||v<min)min=v; if(max==""||v>max)max=v} END{if(min!="")printf("min=%s max=%s",min,max)}' "${SUM_PATH}")"
  log "count_mean over ok points: ${cnt_minmax:-<none>}"
fi

# Adaptive branch stats (from raw)
idx_raw_ok="$(awk_col_idx "${RAW_PATH}" "ok")"
idx_raw_var="$(awk_col_idx "${RAW_PATH}" "variant")"
idx_raw_branch="$(awk_col_idx "${RAW_PATH}" "adaptive_branch")"

if [[ -n "${idx_raw_ok}" && -n "${idx_raw_var}" && -n "${idx_raw_branch}" ]]; then
  branch_stats="$(awk -F',' -v okc="${idx_raw_ok}" -v vc="${idx_raw_var}" -v bc="${idx_raw_branch}" 'NR>1 && $okc==1 && $vc=="adaptive" {cnt[$bc]++} END{for(k in cnt){printf("%s:%d ",k,cnt[k])}}' "${RAW_PATH}")"
  log "Adaptive branches (ok runs): ${branch_stats:-<none>}"
fi

# ----------------------------
# Plot
# ----------------------------
if [[ "${DO_PLOT}" -eq 1 ]]; then
  if command -v python3 >/dev/null 2>&1; then
    # Prefer new location; fall back for legacy repos.
    PLOT_SCRIPT="${REPO_ROOT}/run/include/exp2_plot.py"
    if [[ ! -f "${PLOT_SCRIPT}" ]]; then
      PLOT_SCRIPT="${REPO_ROOT}/run/plot_exp2.py"
    fi
    [[ -f "${PLOT_SCRIPT}" ]] || die "Missing plot script at run/include/exp2_plot.py or run/plot_exp2.py"

    log "Generating plots via: ${PLOT_SCRIPT}"
    python3 "${PLOT_SCRIPT}" \
      --out_dir "${OUT_DIR}" \
      --t0 "${PLOT_T0}" \
      --error "${PLOT_ERROR}" \
      --warmup_reps "${WARMUP}" \
      --mode "${PLOT_MODE}" \
      --topk "${TOPK}" \
      --exclude_methods "${EXCLUDE_METHODS}" \
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
