#!/usr/bin/env bash
# run.sh (repo root)
#
# Legacy 2D throughput 实验脚本
# --------------------------------
# 该脚本保留原 2D 实验计划（A*~D*）的一键复现流程。
# 对于当前仓库的通用 ND(2/3/4/5) 运行入口，请使用
#   ./build/sjs_run / ./build/sjs_sweep / ./build/sjs_verify / ./build/sjs_gen_dataset
# 并参考 README.md 与 EXTENSION_NOTES.md。
#
# 本版本用于“更大参数网络”的 TH（throughput）评估：
#   - 重点把 t 提升到 1e7~1e9 区间（SJS 主战场），并扫描 alpha_out / N / family。
#
# 主实验：
#   (A*) alpha_out 扫描（密度敏感性, TH）
#        N=1e6, t=1e7, d=2, alpha_out in {0.1,0.3,1,3,10,30,100,300,1000}
#        families: F0 + F1
#
#   (B*) N 扫描（输入规模可扩展性, TH）
#        alpha_out=100, t=1e7, d=2
#        N in {1e5,2e5,5e5,1e6,2e6,5e6}
#        family: F0（F1 可选：INCLUDE_B_F1=1）
#
#   (C*) t 扫描（输出规模可扩展性, TH）
#        N=1e6, alpha_out=100, d=2
#        t in {1e6,3e6,1e7,3e7,1e8,3e8,1e9}
#        family: F0（F1 可选：INCLUDE_C_F1=1）
#        额外汇报：samples/sec 与 ns/sample（基于 run_sample_ms）
#
#   (D*) family 敏感性（形状/尺度方差）
#        N=1e6, alpha_out=100, t=1e7, d=2
#        F0: fixed + shape_sigma=0.0
#        F1: lognormal(cv=1.0) + shape_sigma=1.0
#
# 模型：
#   1) KDTree          : --method=kd_tree --variant=sampling
#   2) SJS-Enum        : --method=ours    --variant=enum_sampling   (Framework I)
#   3) SJS-Sampling    : --method=ours    --variant=sampling        (两遍 sweep)
#
# 使用：
#   chmod +x run.sh
#   ./run.sh
#
# 常用覆盖（环境变量）：
#   BUILD_TYPE=Release|Debug|RelWithDebInfo|MinSizeRel
#   CLEAN_BUILD=0|1
#   CLEAN_TEMP=0|1
#   JOBS=8
#   THREADS=1
#   REPEATS=3
#   SEED=1
#
#   # 数据生成（Alacarte-rectgen）
#   PYTHON=python3
#   GEN_SEED=1
#   AUDIT_PAIRS=2000000
#   AUDIT_SEED=1
#
#   # EnumSampling 安全阈值（0 表示不设 cap；注意可能 OOM）
#   # 说明：B* 包含 N=5e6, alpha=100，期望 |J|≈1e9，若希望 enum_sampling 可运行，
#   #      需把 ENUM_CAP 提到 >= 1e9（并确保机器内存足够）。
#   ENUM_CAP=1200000000
#   ENUM_PRECHECK=1    # 1: 若预计 |J| > ENUM_CAP，直接标记 Not runnable 而不实际枚举
#
#   # 子集开关
#   RUN_A=1 RUN_B=1 RUN_C=1 RUN_D=1
#   INCLUDE_B_F1=0
#   INCLUDE_C_F1=0
#   INCLUDE_ALPHA_CTRL=0   # 1: A* 额外加入 ALPHAS_A_CTRL（默认空）
#   RUN_D_SHAPE_SWEEP=0    # 1: 跑 D* 的 shape_sigma 扩展
#
#   # 并行/健壮性
#   MAX_PARALLEL=2
#   MAX_RETRIES=2
#   TIMEOUT_SEC=0
#
#   # 输出/保留
#   RUN_LOG_TO_STDOUT=0
#   RESULT_COPY_DATASETS=0
#   PRUNE_TEMP_DATASETS=1

#   # TH 配置（throughput profile）
#   WRITE_SAMPLES=0
#   VERIFY=0
#
# 说明：
#   - 数据集生成采用 python 脚本 tools/alacarte_rectgen_generate.py（复制为 preset，不改源码），
#     preset 会把 Alacarte 返回的 info（含 coverage / alpha_expected_est / tune_history / seeds 等）
#     写入 report 的 "info" 字段，方便审计。
#   - 所有模型均复用同一份缓存数据集（dataset_source=binary）。
#   - EnumSampling 若被 cap / OOM 等限制，会被标记为 Not runnable（.nr），不会导致整次实验失败。

set -Eeuo pipefail
IFS=$' \t\n'

timestamp_now() { date -Is; }
log()  { printf '[%s][run.sh][INFO] %s\n' "$(timestamp_now)" "$*"; }
warn() { printf '[%s][run.sh][WARN] %s\n' "$(timestamp_now)" "$*" >&2; }
fatal_handler() {
  local line="$1"
  local cmd="$2"
  printf '[%s][run.sh][FATAL] Failed at line %s: %s\n' "$(timestamp_now)" "$line" "$cmd" >&2
}
trap 'fatal_handler "${LINENO}" "${BASH_COMMAND}"' ERR

# ----------------------------
# Small helpers
# ----------------------------

task_count() {
  local tf="$1"
  [[ -f "$tf" ]] || { echo 0; return 0; }
  awk 'BEGIN{c=0} !/^[[:space:]]*($|#)/ {c++} END{print c}' "$tf"
}

nproc_safe() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
  elif command -v getconf >/dev/null 2>&1; then
    getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4
  else
    echo 4
  fi
}

mem_gb_safe() {
  if [[ -r /proc/meminfo ]]; then
    awk '/MemTotal:/ {printf("%.0f", $2/1024/1024)}' /proc/meminfo
  else
    echo 16
  fi
}

build_subdir_for_type() {
  case "$1" in
    Release|release) echo "release" ;;
    Debug|debug) echo "debug" ;;
    RelWithDebInfo|relwithdebinfo) echo "relwithdebinfo" ;;
    MinSizeRel|minsizerel) echo "minsizerel" ;;
    *) echo "$(echo "$1" | tr '[:upper:]' '[:lower:]')" ;;
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
  local p
  p="$(find "${build_dir}" -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "${p}" && -x "${p}" ]] || return 1
  echo "${p}"
}

cmake_cache_internal_value() {
  local cache="$1"
  local key="$2"
  [[ -f "$cache" ]] || return 1
  local line
  line="$(grep -E "^${key}:INTERNAL=" "$cache" | head -n 1 || true)"
  [[ -n "$line" ]] || return 1
  echo "${line#*=}"
}

ensure_cmake_configured() {
  local cache="${BUILD_DIR}/CMakeCache.txt"

  if [[ ! -f "$cache" ]]; then
    log "Configuring CMake..."
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    return 0
  fi

  local cache_home=""
  local cache_dir=""
  cache_home="$(cmake_cache_internal_value "$cache" "CMAKE_HOME_DIRECTORY" || true)"
  cache_dir="$(cmake_cache_internal_value "$cache" "CMAKE_CACHEFILE_DIR" || true)"

  if [[ -n "$cache_home" && "$cache_home" != "$ROOT" ]]; then
    warn "Detected stale CMake cache source path: $cache_home (root is $ROOT)"
    warn "Recreating build dir: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
    log "Configuring CMake..."
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    return 0
  fi

  if [[ -n "$cache_dir" && "$cache_dir" != "$BUILD_DIR" ]]; then
    warn "Detected stale CMake cache build path: $cache_dir (build is $BUILD_DIR)"
    warn "Recreating build dir: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
    log "Configuring CMake..."
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    return 0
  fi
}

sanitize_token() {
  echo "$1" \
    | sed -e 's/-/m/g' -e 's/\./p/g' -e 's/[^A-Za-z0-9_]/_/g'
}

supports_wait_n() {
  local major="${BASH_VERSINFO[0]:-0}"
  local minor="${BASH_VERSINFO[1]:-0}"
  (( major > 4 || (major == 4 && minor >= 3) ))
}

wait_for_slot() {
  local max_parallel="$1"
  while true; do
    local running
    running="$(jobs -pr | wc -l | tr -d '[:space:]')"
    [[ -z "${running}" ]] && running=0
    if (( running < max_parallel )); then
      return 0
    fi
    if supports_wait_n; then
      wait -n || true
    else
      sleep 0.2
    fi
  done
}

csv_col_idx() {
  local file="$1"
  local col="$2"
  awk -F',' -v c="$col" 'NR==1{for(i=1;i<=NF;i++){if($i==c){print i; exit}}}' "$file"
}

csv_row_count() {
  local file="$1"
  awk 'END{print NR-1}' "$file" 2>/dev/null || echo 0
}

csv_ok_count() {
  local file="$1"
  local ok_idx
  ok_idx="$(csv_col_idx "$file" "ok")"
  [[ -n "${ok_idx}" ]] || { echo 0; return 0; }
  awk -F',' -v k="${ok_idx}" 'NR==1{next} {if($k==1) c++} END{print c+0}' "$file"
}

merge_csv_dir() {
  local dir="$1"
  local out_csv="$2"
  local first=""
  first="$(find "$dir" -type f -name run.csv | sort | head -n 1 || true)"
  if [[ -z "$first" ]]; then
    warn "No run.csv files found under: $dir"
    return 0
  fi
  mkdir -p "$(dirname "$out_csv")"
  head -n 1 "$first" > "$out_csv"
  while IFS= read -r f; do
    tail -n +2 "$f" >> "$out_csv" || true
  done < <(find "$dir" -type f -name run.csv | sort)
}

copy_dir_if_exists() {
  local src="$1"
  local dst="$2"
  [[ -d "$src" ]] || return 0
  cp -a "$src" "$dst"
}

dir_size_human() {
  local path="$1"
  [[ -e "$path" ]] || { echo "0"; return 0; }
  du -sh "$path" 2>/dev/null | awk '{print $1}'
}

# ----------------------------
# Resolve repo root
# ----------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ----------------------------
# Parameters (defaults align to 实验计划.md)
# ----------------------------
EXP_TAG="${EXP_TAG:-exp3_th_biggrid_2d}"

BUILD_TYPE="${BUILD_TYPE:-Release}"
CLEAN_BUILD="${CLEAN_BUILD:-0}"
CLEAN_TEMP="${CLEAN_TEMP:-0}"

THREADS="${THREADS:-1}"
REPEATS="${REPEATS:-3}"
SEED="${SEED:-1}"

# TH profile toggles
WRITE_SAMPLES="${WRITE_SAMPLES:-0}"
VERIFY="${VERIFY:-0}"

PYTHON="${PYTHON:-python3}"
GEN_SEED="${GEN_SEED:-1}"
AUDIT_PAIRS="${AUDIT_PAIRS:-2000000}"
AUDIT_SEED="${AUDIT_SEED:-$GEN_SEED}"

# EnumSampling safety
ENUM_CAP="${ENUM_CAP:-1200000000}"
ENUM_PRECHECK="${ENUM_PRECHECK:-1}"

# Which sweeps to run
RUN_A="${RUN_A:-1}"
RUN_B="${RUN_B:-1}"
RUN_C="${RUN_C:-1}"
RUN_D="${RUN_D:-1}"

INCLUDE_B_F1="${INCLUDE_B_F1:-0}"
INCLUDE_C_F1="${INCLUDE_C_F1:-0}"
INCLUDE_ALPHA_CTRL="${INCLUDE_ALPHA_CTRL:-0}"
RUN_D_SHAPE_SWEEP="${RUN_D_SHAPE_SWEEP:-0}"

# Grids
ALPHAS_A_BASE="${ALPHAS_A_BASE:-0.1 0.3 1 3 10 30 100 300 1000}"
ALPHAS_A_CTRL="${ALPHAS_A_CTRL:-}"
N_A="${N_A:-1000000}"
T_A="${T_A:-10000000}"

NS_B="${NS_B:-100000 200000 500000 1000000 2000000 5000000}"
ALPHA_B="${ALPHA_B:-100}"
T_B="${T_B:-10000000}"

TS_C="${TS_C:-1000000 3000000 10000000 30000000 100000000 300000000 1000000000}"
N_C="${N_C:-1000000}"
ALPHA_C="${ALPHA_C:-100}"

N_D="${N_D:-1000000}"
ALPHA_D="${ALPHA_D:-100}"
T_D="${T_D:-10000000}"

# Optional D extension: lognormal volume, shape_sigma sweep
D_SHAPE_SIGMAS="${D_SHAPE_SIGMAS:-0.0 0.5 1.0 1.5 2.0}"

# Robustness / parallelism
MAX_RETRIES="${MAX_RETRIES:-2}"
RETRY_BACKOFF_SEC="${RETRY_BACKOFF_SEC:-2}"
TIMEOUT_SEC="${TIMEOUT_SEC:-0}"

MAX_PARALLEL_CAP="${MAX_PARALLEL_CAP:-8}"
MEM_PER_TASK_GB="${MEM_PER_TASK_GB:-24}"
MEM_RESERVE_GB="${MEM_RESERVE_GB:-8}"

# Conservative split thresholds
N_LARGE_THRESHOLD="${N_LARGE_THRESHOLD:-2000000}"
T_LARGE_THRESHOLD="${T_LARGE_THRESHOLD:-30000000}"

# Logging + retention
RUN_LOG_TO_STDOUT="${RUN_LOG_TO_STDOUT:-0}"
RESULT_COPY_DATASETS="${RESULT_COPY_DATASETS:-0}"
PRUNE_TEMP_DATASETS="${PRUNE_TEMP_DATASETS:-1}"

# Generator scripts
RECTGEN_BASE_SCRIPT="${RECTGEN_BASE_SCRIPT:-$ROOT/tools/alacarte_rectgen_generate.py}"
ALACARTE_MODULE="${ALACARTE_MODULE:-$ROOT/third_party/Alacarte/alacarte_rectgen.py}"

# Families (Alacarte params)
F0_VOLUME_DIST="${F0_VOLUME_DIST:-fixed}"
F0_VOLUME_CV="${F0_VOLUME_CV:-0.0}"        # fixed 时该值通常无意义；设 0 更直观
F0_SHAPE_SIGMA="${F0_SHAPE_SIGMA:-0.0}"

F1_VOLUME_DIST="${F1_VOLUME_DIST:-lognormal}"
F1_VOLUME_CV="${F1_VOLUME_CV:-1.0}"
F1_SHAPE_SIGMA="${F1_SHAPE_SIGMA:-1.0}"

# Build parallelism
if [[ -n "${JOBS:-}" ]]; then
  JOBS="$JOBS"
else
  local_nproc="$(nproc_safe)"
  local_mem_gb="$(mem_gb_safe)"
  JOBS="$local_mem_gb"
  [[ -z "$JOBS" || "$JOBS" -lt 1 ]] && JOBS=1
  (( JOBS > local_nproc )) && JOBS="$local_nproc"
  (( JOBS > 16 )) && JOBS=16
fi

# ----------------------------
# Model set
# ----------------------------
MODELS_MAIN=(
  "ours enum_sampling"
  "ours sampling"
  "kd_tree sampling"
)

# ----------------------------
# Directories
# ----------------------------
BUILD_SUBDIR="$(build_subdir_for_type "$BUILD_TYPE")"
BUILD_DIR="${ROOT}/build/${BUILD_SUBDIR}"

TEMP_ROOT="${ROOT}/run/temp/${EXP_TAG}"
OUT_ROOT="${TEMP_ROOT}/out"
LOG_ROOT="${TEMP_ROOT}/logs"
STATUS_ROOT="${TEMP_ROOT}/status"
MANIFEST_DIR="${TEMP_ROOT}/manifest"

DATA_ROOT="${TEMP_ROOT}/datasets"
DATA_STATUS_ROOT="${TEMP_ROOT}/dataset_status"
DATA_LOG_ROOT="${TEMP_ROOT}/dataset_logs"
PRESET_DIR="${TEMP_ROOT}/rectgen_presets"

RESULT_ROOT="${ROOT}/results/raw/${EXP_TAG}"
RUN_LOG_DIR="${ROOT}/log"
RUN_LOG_STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_LOG_FILE="${RUN_LOG_DIR}/run_${EXP_TAG}_${RUN_LOG_STAMP}_pid$$.log"

setup_run_logger() {
  mkdir -p "$RUN_LOG_DIR"
  if [[ "$RUN_LOG_TO_STDOUT" == "1" ]]; then
    exec > >(tee -a "$RUN_LOG_FILE") 2>&1
  else
    exec >>"$RUN_LOG_FILE" 2>&1
  fi
}

setup_run_logger
log "Run log initialized: $RUN_LOG_FILE"

if [[ "$CLEAN_TEMP" == "1" ]]; then
  warn "CLEAN_TEMP=1: removing temp dir: $TEMP_ROOT"
  rm -rf "$TEMP_ROOT"
fi

mkdir -p "$TEMP_ROOT" "$OUT_ROOT" "$LOG_ROOT" "$STATUS_ROOT" "$MANIFEST_DIR" \
         "$DATA_ROOT" "$DATA_STATUS_ROOT" "$DATA_LOG_ROOT" "$PRESET_DIR"

# ----------------------------
# Build
# ----------------------------
log "Repo root  : $ROOT"
log "Build type : $BUILD_TYPE"
log "Build dir  : $BUILD_DIR"
log "Temp dir   : $TEMP_ROOT"
log "Results dir: $RESULT_ROOT"

if [[ "$CLEAN_BUILD" == "1" ]]; then
  log "Cleaning build dir: $BUILD_DIR"
  rm -rf "$BUILD_DIR"
fi

ensure_cmake_configured

log "Building... (JOBS=$JOBS)"
cmake --build "$BUILD_DIR" -j "$JOBS"

SJS_RUN="$(find_exe "$BUILD_DIR" sjs_run || true)"
if [[ -z "$SJS_RUN" ]]; then
  echo "[run.sh][FATAL] Could not locate sjs_run under $BUILD_DIR" >&2
  exit 2
fi
log "Using sjs_run: $SJS_RUN"

# ----------------------------
# Parallelism defaults
# ----------------------------
NPROC="$(nproc_safe)"
MEM_GB="$(mem_gb_safe)"

if [[ -n "${MAX_PARALLEL:-}" ]]; then
  MAX_PARALLEL_USER="$MAX_PARALLEL"
else
  cpu_per_task="$THREADS"
  (( cpu_per_task < 1 )) && cpu_per_task=1
  cpu_cap=$(( NPROC / cpu_per_task ))
  (( cpu_cap < 1 )) && cpu_cap=1
  (( cpu_cap > 1 )) && cpu_cap=$(( cpu_cap - 1 ))

  mem_budget=$(( MEM_GB - MEM_RESERVE_GB ))
  (( mem_budget < MEM_PER_TASK_GB )) && mem_budget="$MEM_PER_TASK_GB"
  mem_cap=$(( mem_budget / MEM_PER_TASK_GB ))
  (( mem_cap < 1 )) && mem_cap=1

  MAX_PARALLEL_USER="$cpu_cap"
  (( MAX_PARALLEL_USER > mem_cap )) && MAX_PARALLEL_USER="$mem_cap"
  (( MAX_PARALLEL_USER > MAX_PARALLEL_CAP )) && MAX_PARALLEL_USER="$MAX_PARALLEL_CAP"
fi
MAX_PARALLEL="$MAX_PARALLEL_USER"
(( MAX_PARALLEL < 1 )) && MAX_PARALLEL=1

# Per-sweep parallelism (override if needed)
PAR_A="${PAR_A:-$MAX_PARALLEL}"
PAR_B_SMALL="${PAR_B_SMALL:-$(( MAX_PARALLEL > 2 ? 2 : MAX_PARALLEL ))}"
PAR_B_LARGE="${PAR_B_LARGE:-1}"
PAR_C_SMALL="${PAR_C_SMALL:-$(( MAX_PARALLEL > 3 ? 3 : MAX_PARALLEL ))}"
PAR_C_LARGE="${PAR_C_LARGE:-1}"
PAR_D="${PAR_D:-$MAX_PARALLEL}"

log "Host: nproc=$NPROC, mem~${MEM_GB}GB"
log "Parallelism: max=$MAX_PARALLEL (A=$PAR_A, B_small=$PAR_B_SMALL, B_large=$PAR_B_LARGE, C_small=$PAR_C_SMALL, C_large=$PAR_C_LARGE, D=$PAR_D)"
log "Enum: cap=$ENUM_CAP, precheck=$ENUM_PRECHECK"

# ----------------------------
# Build run signature (for resumability)
# ----------------------------
build_run_signature() {
  {
    echo "exp_tag=$EXP_TAG"
    echo "build_type=$BUILD_TYPE"
    echo "threads=$THREADS"
    echo "repeats=$REPEATS"
    echo "seed=$SEED"

    echo "python=$PYTHON"
    echo "gen_seed=$GEN_SEED"
    echo "audit_pairs=$AUDIT_PAIRS"
    echo "audit_seed=$AUDIT_SEED"
    echo "rectgen_base_script=$RECTGEN_BASE_SCRIPT"
    echo "alacarte_module=$ALACARTE_MODULE"

    echo "F0=$F0_VOLUME_DIST,$F0_VOLUME_CV,$F0_SHAPE_SIGMA"
    echo "F1=$F1_VOLUME_DIST,$F1_VOLUME_CV,$F1_SHAPE_SIGMA"

    echo "enum_cap=$ENUM_CAP"
    echo "enum_precheck=$ENUM_PRECHECK"

    echo "RUN_A=$RUN_A RUN_B=$RUN_B RUN_C=$RUN_C RUN_D=$RUN_D"
    echo "INCLUDE_B_F1=$INCLUDE_B_F1 INCLUDE_C_F1=$INCLUDE_C_F1"
    echo "INCLUDE_ALPHA_CTRL=$INCLUDE_ALPHA_CTRL"
    echo "RUN_D_SHAPE_SWEEP=$RUN_D_SHAPE_SWEEP D_SHAPE_SIGMAS=$D_SHAPE_SIGMAS"

    echo "ALPHAS_A_BASE=$ALPHAS_A_BASE ALPHAS_A_CTRL=$ALPHAS_A_CTRL N_A=$N_A T_A=$T_A"
    echo "NS_B=$NS_B ALPHA_B=$ALPHA_B T_B=$T_B"
    echo "TS_C=$TS_C N_C=$N_C ALPHA_C=$ALPHA_C"
    echo "N_D=$N_D ALPHA_D=$ALPHA_D T_D=$T_D"

    echo "N_LARGE_THRESHOLD=$N_LARGE_THRESHOLD"
    echo "T_LARGE_THRESHOLD=$T_LARGE_THRESHOLD"
  } | cksum | awk '{print $1 "-" $2}'
}

RUN_SIGNATURE="$(build_run_signature)"
log "Run signature: $RUN_SIGNATURE"

# ----------------------------
# Prepare rectgen preset scripts (F0/F1 + optional D shape sweep)
# ----------------------------

if [[ ! -f "$RECTGEN_BASE_SCRIPT" ]]; then
  echo "[run.sh][FATAL] RECTGEN_BASE_SCRIPT not found: $RECTGEN_BASE_SCRIPT" >&2
  exit 2
fi
if [[ ! -f "$ALACARTE_MODULE" ]]; then
  echo "[run.sh][FATAL] ALACARTE_MODULE not found: $ALACARTE_MODULE" >&2
  exit 2
fi

# Patches:
#  - override volume_dist / volume_cv / shape_sigma
#  - inject `report["info"] = info` before report is written
make_rectgen_preset() {
  local src="$1"; local dst="$2"; local vol_dist="$3"; local vol_cv="$4"; local shape_sigma="$5"; local tag="$6"

  "$PYTHON" - "$src" "$dst" "$vol_dist" "$vol_cv" "$shape_sigma" "$tag" <<'PY'
import re
import sys
from pathlib import Path

src, dst, vol_dist, vol_cv, shape_sigma, tag = sys.argv[1:7]
text = Path(src).read_text(encoding='utf-8')

def sub_one(pattern, repl, name):
    global text
    new, n = re.subn(pattern, repl, text, count=1, flags=re.MULTILINE)
    if n != 1:
        raise SystemExit(f"[preset][FATAL] Expected to patch 1 occurrence of {name}, got {n}")
    text = new

# param overrides
sub_one(r'^\s*volume_dist\s*=\s*"[^"]+"\s*$', f'    volume_dist = "{vol_dist}"  # preset:{tag}', 'volume_dist')
sub_one(r'^\s*volume_cv\s*=\s*[-+0-9.eE]+\s*$', f'    volume_cv = {vol_cv}  # preset:{tag}', 'volume_cv')
sub_one(r'^\s*shape_sigma\s*=\s*[-+0-9.eE]+\s*$', f'    shape_sigma = {shape_sigma}  # preset:{tag}', 'shape_sigma')

# inject report["info"] = info once
needle = r'^\s*if\s+args\.report_path\s*:\s*$'
ins = '    # preset: inject full Alacarte info for auditing\n    report["info"] = info\n\n'
new, n = re.subn(needle, ins + r'\g<0>', text, count=1, flags=re.MULTILINE)
if n != 1:
    raise SystemExit(f"[preset][FATAL] Expected to inject info before report_path, got {n}")
text = new

Path(dst).write_text(text, encoding='utf-8')
PY

  chmod +x "$dst"
}

# Map family -> preset script
declare -A RECTGEN_PRESET

RECTGEN_F0="$PRESET_DIR/alacarte_rectgen_F0.py"
RECTGEN_F1="$PRESET_DIR/alacarte_rectgen_F1.py"
make_rectgen_preset "$RECTGEN_BASE_SCRIPT" "$RECTGEN_F0" "$F0_VOLUME_DIST" "$F0_VOLUME_CV" "$F0_SHAPE_SIGMA" "F0"
make_rectgen_preset "$RECTGEN_BASE_SCRIPT" "$RECTGEN_F1" "$F1_VOLUME_DIST" "$F1_VOLUME_CV" "$F1_SHAPE_SIGMA" "F1"
RECTGEN_PRESET["F0"]="$RECTGEN_F0"
RECTGEN_PRESET["F1"]="$RECTGEN_F1"

# Optional D shape sweep presets (lognormal + cv=1.0, varying shape_sigma)
if [[ "$RUN_D_SHAPE_SWEEP" == "1" ]]; then
  for sigma in $D_SHAPE_SIGMAS; do
    tag="F1sigma$(sanitize_token "$sigma")"
    dst="$PRESET_DIR/alacarte_rectgen_${tag}.py"
    make_rectgen_preset "$RECTGEN_BASE_SCRIPT" "$dst" "lognormal" "1.0" "$sigma" "$tag"
    RECTGEN_PRESET["$tag"]="$dst"
  done
fi

rectgen_for_family() {
  local fam="$1"
  if [[ -n "${RECTGEN_PRESET[$fam]+x}" ]]; then
    echo "${RECTGEN_PRESET[$fam]}"
  else
    echo "$RECTGEN_F0"
  fi
}

log "Rectgen presets under: $PRESET_DIR"
log "  F0 => ${RECTGEN_PRESET[F0]}"
log "  F1 => ${RECTGEN_PRESET[F1]}"
if [[ "$RUN_D_SHAPE_SWEEP" == "1" ]]; then
  log "  D-shape presets:"
  for sigma in $D_SHAPE_SIGMAS; do
    tag="F1sigma$(sanitize_token "$sigma")"
    log "    $tag => ${RECTGEN_PRESET[$tag]}"
  done
fi

# ----------------------------
# Dataset cache helpers
# ----------------------------

ds_name_for() {
  local fam="$1"; local N="$2"; local alpha="$3"
  echo "d2_${fam}_n${N}_a$(sanitize_token "$alpha")_gs${GEN_SEED}"
}

ds_dir_for() {
  local fam="$1"; echo "$DATA_ROOT/$fam"
}

ds_paths_for() {
  local fam="$1"; local N="$2"; local alpha="$3"
  local name; name="$(ds_name_for "$fam" "$N" "$alpha")"
  local dir; dir="$(ds_dir_for "$fam")"
  echo "$dir/${name}_R.bin"$'\t'"$dir/${name}_S.bin"$'\t'"$dir/${name}_gen_report.json"$'\t'"$name"
}

minify_json_inplace() {
  local json_path="$1"
  [[ -f "$json_path" ]] || return 0
  "$PYTHON" - "$json_path" <<'PY'
import json
import sys
from pathlib import Path

p = Path(sys.argv[1])
try:
    obj = json.loads(p.read_text(encoding='utf-8'))
except Exception:
    sys.exit(0)
# compact form
p.write_text(json.dumps(obj, separators=(',', ':'), ensure_ascii=False), encoding='utf-8')
PY
}

patch_run_csv_gen_report() {
  local csv_file="$1"
  local gen_json_path="$2"
  [[ -f "$csv_file" && -f "$gen_json_path" ]] || return 0

  "$PYTHON" - "$csv_file" "$gen_json_path" <<'PY'
import csv
import sys
from pathlib import Path

csv_path, gen_path = sys.argv[1], sys.argv[2]

gen_json = Path(gen_path).read_text(encoding='utf-8').strip() or "{}"

with open(csv_path, newline='') as f:
    reader = csv.reader(f)
    rows = list(reader)
if not rows:
    sys.exit(0)
header = rows[0]
try:
    idx = header.index('gen_report_json')
except ValueError:
    sys.exit(0)

out_rows = [header]
for row in rows[1:]:
    if len(row) < len(header):
        row = row + [''] * (len(header) - len(row))
    row[idx] = gen_json
    out_rows.append(row)

with open(csv_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(out_rows)
PY
}

# Fill note column for failed reps using per-task log (best-effort)
patch_run_csv_note_from_log() {
  local csv_file="$1"
  local log_file="$2"
  [[ -f "$csv_file" && -f "$log_file" ]] || return 0

  "$PYTHON" - "$csv_file" "$log_file" <<'PY'
import csv
import re
import sys
from pathlib import Path

csv_path, log_path = sys.argv[1], sys.argv[2]

# Extract errors: "Run failed (rep= X ): <msg>"
rep_err = {}
pat = re.compile(r"Run failed \(rep=\s*(\d+)\s*\)\s*:\s*(.*)")
for line in Path(log_path).read_text(encoding='utf-8', errors='ignore').splitlines():
    m = pat.search(line)
    if not m:
        continue
    rep = int(m.group(1))
    msg = m.group(2).strip()
    # keep first occurrence
    rep_err.setdefault(rep, msg)

with open(csv_path, newline='') as f:
    rows = list(csv.reader(f))
if not rows:
    sys.exit(0)

hdr = rows[0]
try:
    rep_i = hdr.index('rep')
    ok_i  = hdr.index('ok')
    note_i = hdr.index('note')
except ValueError:
    sys.exit(0)

out = [hdr]
changed = False
for r in rows[1:]:
    if len(r) < len(hdr):
        r = r + [''] * (len(hdr) - len(r))
    try:
        rep = int(r[rep_i])
    except Exception:
        rep = None
    ok = r[ok_i]
    if ok != '1':
        if (not r[note_i]) and (rep is not None) and (rep in rep_err):
            r[note_i] = rep_err[rep]
            changed = True
    out.append(r)

if changed:
    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerows(out)
PY
}

# ----------------------------
# Dataset generation (python rectgen presets)
# ----------------------------

ensure_dataset() {
  local fam="$1"; local N="$2"; local alpha="$3"

  local paths
  paths="$(ds_paths_for "$fam" "$N" "$alpha")"
  local path_r path_s rep name
  IFS=$'\t' read -r path_r path_s rep name <<< "$paths"

  mkdir -p "$(ds_dir_for "$fam")"

  local ok_flag="$DATA_STATUS_ROOT/${name}.ok"
  local meta_flag="$DATA_STATUS_ROOT/${name}.meta"
  local log_file="$DATA_LOG_ROOT/${name}.log"

  if [[ -f "$ok_flag" && -f "$path_r" && -f "$path_s" && -f "$rep" ]]; then
    if grep -Fxq "run_signature=$RUN_SIGNATURE" "$meta_flag" 2>/dev/null; then
      return 0
    fi
  fi

  local lock_dir="$DATA_STATUS_ROOT/${name}.lock"
  while ! mkdir "$lock_dir" 2>/dev/null; do
    if [[ -f "$ok_flag" && -f "$path_r" && -f "$path_s" && -f "$rep" ]]; then
      return 0
    fi
    sleep 0.5
  done

  rm -f "$ok_flag" "$meta_flag"

  {
    echo "run_signature=$RUN_SIGNATURE"
    echo "family=$fam"
    echo "N=$N"
    echo "alpha=$alpha"
    echo "gen_seed=$GEN_SEED"
    echo "rectgen_script=$(rectgen_for_family "$fam")"
    echo "alacarte_module=$ALACARTE_MODULE"
    echo "audit_pairs=$AUDIT_PAIRS"
    echo "audit_seed=$AUDIT_SEED"
  } > "$meta_flag"

  log "[dataset] Generating $name (fam=$fam N=$N alpha=$alpha)"

  local preset
  preset="$(rectgen_for_family "$fam")"

  local -a gen_cmd=(
    "$PYTHON" "$preset"
    "--nR=$N" "--nS=$N" "--d=2"
    "--alpha_out=$alpha"
    "--seed=$GEN_SEED"
    "--out_r=$path_r" "--out_s=$path_s"
    "--dataset_name=$name"
    "--report_path=$rep"
    "--audit_pairs=$AUDIT_PAIRS" "--audit_seed=$AUDIT_SEED"
    "--alacarte_module=$ALACARTE_MODULE"
  )

  {
    echo "[$(timestamp_now)][run.sh][dataset] start name=$name fam=$fam N=$N alpha=$alpha"
    echo "[$(timestamp_now)][run.sh][dataset] cmd=${gen_cmd[*]}"
  } > "$log_file"

  local rc=0
  if (( TIMEOUT_SEC > 0 )) && command -v timeout >/dev/null 2>&1; then
    timeout "$TIMEOUT_SEC" "${gen_cmd[@]}" >>"$log_file" 2>&1 || rc=$?
  else
    "${gen_cmd[@]}" >>"$log_file" 2>&1 || rc=$?
  fi
  echo "[$(timestamp_now)][run.sh][dataset] end rc=$rc name=$name" >> "$log_file"

  if [[ "$rc" -ne 0 ]]; then
    warn "[dataset] FAILED rc=$rc name=$name (see $log_file)"
    rmdir "$lock_dir" || true
    return 1
  fi

  if [[ ! -f "$path_r" || ! -f "$path_s" || ! -f "$rep" ]]; then
    warn "[dataset] FAILED: missing outputs for $name"
    warn "  R=$path_r"
    warn "  S=$path_s"
    warn "  rep=$rep"
    rmdir "$lock_dir" || true
    return 1
  fi

  # Minify JSON for easier CSV embedding
  minify_json_inplace "$rep" || true

  touch "$ok_flag"
  rmdir "$lock_dir" || true
  return 0
}

# ----------------------------
# Task registry (pre-generate datasets)
# ----------------------------

declare -A DATASET_SPEC
DATASET_ORDER=()

register_dataset() {
  local fam="$1"; local N="$2"; local alpha="$3"
  local name
  name="$(ds_name_for "$fam" "$N" "$alpha")"
  if [[ -z "${DATASET_SPEC[$name]+x}" ]]; then
    DATASET_SPEC[$name]="${fam}"$'\t'"${N}"$'\t'"${alpha}"
    DATASET_ORDER+=("$name")
  fi
}

# ----------------------------
# Task runner
# ----------------------------

write_not_runnable_csv() {
  local csv_file="$1"
  local ds_name="$2"
  local method="$3"
  local variant="$4"
  local N="$5"
  local t="$6"
  local reason="$7"

  "$PYTHON" - "$csv_file" "$ds_name" "$method" "$variant" "$N" "$t" "$REPEATS" "$SEED" "$ENUM_CAP" "$reason" <<'PY'
import csv
import json
import sys
from pathlib import Path

(csv_path, ds_name, method, variant, N, t, repeats, seed0, enum_cap, reason) = sys.argv[1:]
N = int(N)
t = int(t)
repeats = int(repeats)
seed0 = int(seed0)
enum_cap = int(enum_cap)

hdr = [
  "dataset","method","variant","rep","seed","n_r","n_s","t","ok","wall_ms",
  "count_value","count_exact","count_stderr","count_ci_low","count_ci_high",
  "used_enumeration","enum_truncated","enum_cap","adaptive_branch","adaptive_pilot_pairs",
  "phases_json","note","config_json","gen_report_json",
]

cfg = {
  "dataset": {"name": ds_name, "source": "binary", "dim": 2},
  "run": {"method": method, "variant": variant, "t": t, "seed": seed0, "repeats": repeats, "threads": 1, "enum_cap": enum_cap},
  "note": "not_runnable_stub"
}

Path(csv_path).parent.mkdir(parents=True, exist_ok=True)
with open(csv_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(hdr)
    for rep in range(repeats):
        seed = seed0 + rep
        row = [
          ds_name, method, variant, rep, seed, N, N, t,
          0, "", "", 0, "", "", "",
          1 if variant=="enum_sampling" else 0,
          0, enum_cap,
          "", 0,
          "{}",
          reason,
          json.dumps(cfg, separators=(',',':')),
          "{}",
        ]
        w.writerow(row)
PY
}

# Expected join size from alpha_out definition: |J| ≈ alpha_out*(|R|+|S|) = alpha*2N
expected_join_size() {
  local N="$1"; local alpha="$2"
  "$PYTHON" - "$N" "$alpha" <<'PY'
import sys
N = int(float(sys.argv[1]))
a = float(sys.argv[2])
# alpha_out is |J|/(|R|+|S|); assume |R|=|S|=N
print(int(round(a * (2.0 * N))))
PY
}

run_one_task() {
  local task_id="$1"
  local sweep="$2"
  local fam="$3"
  local N="$4"
  local alpha="$5"
  local t="$6"
  local method="$7"
  local variant="$8"
  local force="$9"
  local attempt="${10:-0}"

  local ok_flag="$STATUS_ROOT/${task_id}.ok"
  local nr_flag="$STATUS_ROOT/${task_id}.nr"
  local fail_flag="$STATUS_ROOT/${task_id}.fail"
  local exit_flag="$STATUS_ROOT/${task_id}.exit"
  local meta_flag="$STATUS_ROOT/${task_id}.meta"

  local out_dir="$OUT_ROOT/${sweep}/${task_id}"
  local log_file="$LOG_ROOT/${task_id}.log"
  local csv_file="$out_dir/run.csv"

  # fast path (already ok / nr with same run_signature)
  if [[ "$force" != "1" && -f "$meta_flag" ]]; then
    if grep -Fxq "run_signature=$RUN_SIGNATURE" "$meta_flag" 2>/dev/null; then
      if [[ -f "$ok_flag" && -f "$csv_file" ]]; then
        echo 0 > "$exit_flag"
        return 0
      fi
      if [[ -f "$nr_flag" && -f "$csv_file" ]]; then
        echo 0 > "$exit_flag"
        return 0
      fi
    fi
  fi

  rm -f "$ok_flag" "$nr_flag" "$fail_flag" "$exit_flag" "$meta_flag"
  rm -rf "$out_dir"
  mkdir -p "$out_dir"

  if ! ensure_dataset "$fam" "$N" "$alpha"; then
    echo "run_signature=$RUN_SIGNATURE" > "$meta_flag"
    echo "dataset_generation_failed=1" >> "$meta_flag"
    touch "$fail_flag"
    echo 1 > "$exit_flag"
    return 0
  fi

  local paths
  paths="$(ds_paths_for "$fam" "$N" "$alpha")"
  local path_r path_s rep_path ds_name
  IFS=$'\t' read -r path_r path_s rep_path ds_name <<< "$paths"

  # EnumSampling precheck: if expected |J| > ENUM_CAP, mark Not runnable (no actual run)
  if [[ "$method" == "ours" && "$variant" == "enum_sampling" && "$ENUM_PRECHECK" == "1" && "$ENUM_CAP" -gt 0 ]]; then
    local expJ
    expJ="$(expected_join_size "$N" "$alpha")"
    if [[ "$expJ" -gt "$ENUM_CAP" ]]; then
      local reason
      reason="Not runnable: expected |J|≈${expJ} > enum_cap=${ENUM_CAP} (precheck)"

      {
        echo "run_signature=$RUN_SIGNATURE"
        echo "attempt=$attempt"
        echo "task_id=$task_id"
        echo "sweep=$sweep"
        echo "family=$fam"
        echo "dataset=$ds_name"
        echo "N=$N"
        echo "alpha=$alpha"
        echo "t=$t"
        echo "method=$method"
        echo "variant=$variant"
        echo "precheck_expected_J=$expJ"
        echo "enum_cap=$ENUM_CAP"
        echo "reason=$reason"
        echo "path_r=$path_r"
        echo "path_s=$path_s"
        echo "gen_report=$rep_path"
      } > "$meta_flag"

      echo "[$(timestamp_now)][run.sh][task] NR(precheck) task_id=$task_id $reason" > "$log_file"

      write_not_runnable_csv "$csv_file" "$ds_name" "$method" "$variant" "$N" "$t" "$reason"
      patch_run_csv_gen_report "$csv_file" "$rep_path" || true

      touch "$nr_flag"
      echo 0 > "$exit_flag"
      return 0
    fi
  fi

  local -a cmd=(
    "$SJS_RUN"
    "--dataset_source=binary"
    "--dataset=$ds_name"
    "--dim=2"
    "--path_r=$path_r" "--path_s=$path_s"
    "--method=$method"
    "--variant=$variant"
    "--t=$t"
    "--seed=$SEED"
    "--repeats=$REPEATS"
    "--threads=$THREADS"
    "--enum_cap=$ENUM_CAP"
    "--write_samples=$WRITE_SAMPLES"
    "--verify=$VERIFY"
    "--out_dir=$out_dir"
    "--results_file=$csv_file"
    "--log_level=info"
    "--log_timestamp=1"
  )

  {
    echo "run_signature=$RUN_SIGNATURE"
    echo "attempt=$attempt"
    echo "task_id=$task_id"
    echo "sweep=$sweep"
    echo "family=$fam"
    echo "dataset=$ds_name"
    echo "N=$N"
    echo "alpha=$alpha"
    echo "t=$t"
    echo "method=$method"
    echo "variant=$variant"
    echo "enum_cap=$ENUM_CAP"
    echo "path_r=$path_r"
    echo "path_s=$path_s"
    echo "gen_report=$rep_path"
    echo "cmd=${cmd[*]}"
  } > "$meta_flag"

  {
    echo "[$(timestamp_now)][run.sh][task] start task_id=$task_id sweep=$sweep family=$fam N=$N alpha=$alpha t=$t method=$method variant=$variant attempt=$attempt"
    echo "[$(timestamp_now)][run.sh][task] cmd=${cmd[*]}"
  } > "$log_file"

  local rc=0
  if (( TIMEOUT_SEC > 0 )) && command -v timeout >/dev/null 2>&1; then
    timeout "$TIMEOUT_SEC" "${cmd[@]}" >>"$log_file" 2>&1 || rc=$?
  else
    "${cmd[@]}" >>"$log_file" 2>&1 || rc=$?
  fi
  echo "[$(timestamp_now)][run.sh][task] end task_id=$task_id rc=$rc" >> "$log_file"
  echo "$rc" > "$exit_flag"

  if [[ "$rc" -ne 0 ]]; then
    echo "nonzero_exit=$rc" >> "$meta_flag"
    # EnumSampling: treat crash as Not runnable; others: fail
    if [[ "$method" == "ours" && "$variant" == "enum_sampling" ]]; then
      local reason="Not runnable: process exit rc=$rc"
      write_not_runnable_csv "$csv_file" "$ds_name" "$method" "$variant" "$N" "$t" "$reason" || true
      patch_run_csv_gen_report "$csv_file" "$rep_path" || true
      touch "$nr_flag"
      return 0
    fi
    touch "$fail_flag"
    return 0
  fi
  if [[ ! -f "$csv_file" ]]; then
    echo "missing_csv=1" >> "$meta_flag"
    if [[ "$method" == "ours" && "$variant" == "enum_sampling" ]]; then
      local reason="Not runnable: missing run.csv (rc=0 but no output)"
      write_not_runnable_csv "$csv_file" "$ds_name" "$method" "$variant" "$N" "$t" "$reason" || true
      patch_run_csv_gen_report "$csv_file" "$rep_path" || true
      touch "$nr_flag"
      return 0
    fi
    touch "$fail_flag"
    return 0
  fi

  # Patch gen_report_json and (best-effort) note for failures
  patch_run_csv_gen_report "$csv_file" "$rep_path" || true
  patch_run_csv_note_from_log "$csv_file" "$log_file" || true

  local rows
  rows="$(csv_row_count "$csv_file")"
  if [[ "$rows" -lt "$REPEATS" ]]; then
    echo "short_csv_rows=$rows" >> "$meta_flag"
    if [[ "$method" == "ours" && "$variant" == "enum_sampling" ]]; then
      local reason="Not runnable: short CSV rows=$rows (<repeats=$REPEATS)"
      # keep produced CSV but mark NR
      touch "$nr_flag"
      echo "nr_reason=$reason" >> "$meta_flag"
      return 0
    fi
    touch "$fail_flag"
    return 0
  fi

  local okc
  okc="$(csv_ok_count "$csv_file")"

  if [[ "$okc" -eq "$REPEATS" ]]; then
    touch "$ok_flag"
    return 0
  fi

  # Not all ok
  echo "ok_count=$okc" >> "$meta_flag"

  if [[ "$method" == "ours" && "$variant" == "enum_sampling" && "$okc" -eq 0 ]]; then
    # EnumSampling not runnable (cap / OOM / other failure)
    touch "$nr_flag"
    return 0
  fi

  # otherwise treat as failure (will retry)
  touch "$fail_flag"
  return 0
}

run_task_file_with_parallel() {
  local task_file="$1"
  local parallel="$2"
  local attempt="$3"
  local force="$4"

  log "Running tasks: $(basename "$task_file") (attempt=$attempt, parallel=$parallel, force=$force)"

  while IFS=$'\t' read -r task_id sweep fam N alpha t method variant; do
    [[ -z "$task_id" ]] && continue
    [[ "$task_id" =~ ^# ]] && continue
    wait_for_slot "$parallel"
    (run_one_task "$task_id" "$sweep" "$fam" "$N" "$alpha" "$t" "$method" "$variant" "$force" "$attempt") &
  done < "$task_file"

  wait || true
}

collect_failures_from_task_file() {
  local task_file="$1"
  local out_fail_list="$2"
  : > "$out_fail_list"

  while IFS=$'\t' read -r task_id sweep fam N alpha t method variant; do
    [[ -z "$task_id" ]] && continue
    [[ "$task_id" =~ ^# ]] && continue
    if [[ ! -f "$STATUS_ROOT/${task_id}.ok" && ! -f "$STATUS_ROOT/${task_id}.nr" ]]; then
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$task_id" "$sweep" "$fam" "$N" "$alpha" "$t" "$method" "$variant" \
        >> "$out_fail_list"
    fi
  done < "$task_file"
}

write_manifest_csv() {
  local task_file="$1"
  local out_csv="$2"
  mkdir -p "$(dirname "$out_csv")"
  echo "task_id,sweep,family,N,alpha,t,method,variant,out_dir,log_file,status" > "$out_csv"
  while IFS=$'\t' read -r task_id sweep fam N alpha t method variant; do
    [[ -z "$task_id" ]] && continue
    [[ "$task_id" =~ ^# ]] && continue
    local out_dir="$OUT_ROOT/${sweep}/${task_id}"
    local log_file="$LOG_ROOT/${task_id}.log"
    local status=""
    if [[ -f "$STATUS_ROOT/${task_id}.ok" ]]; then status="OK"; fi
    if [[ -f "$STATUS_ROOT/${task_id}.nr" ]]; then status="NR"; fi
    if [[ -f "$STATUS_ROOT/${task_id}.fail" ]]; then status="FAIL"; fi
    echo "${task_id},${sweep},${fam},${N},${alpha},${t},${method},${variant},${out_dir},${log_file},${status}" >> "$out_csv"
  done < "$task_file"
}

run_with_retries() {
  local task_file="$1"
  local parallel="$2"
  local name="$3"

  run_task_file_with_parallel "$task_file" "$parallel" 0 0

  local tmp_fail="$MANIFEST_DIR/${name}_FAILURES.tsv"
  collect_failures_from_task_file "$task_file" "$tmp_fail"
  local failed_now
  failed_now="$(task_count "$tmp_fail")"

  local attempt=1
  while (( failed_now > 0 && attempt <= MAX_RETRIES )); do
    warn "$name: retrying failures (attempt $attempt / $MAX_RETRIES, failed=$failed_now)"
    if (( RETRY_BACKOFF_SEC > 0 )); then
      sleep $(( RETRY_BACKOFF_SEC * attempt ))
    fi
    run_task_file_with_parallel "$tmp_fail" "$parallel" "$attempt" 1
    collect_failures_from_task_file "$task_file" "$tmp_fail"
    failed_now="$(task_count "$tmp_fail")"
    attempt=$((attempt + 1))
  done

  if (( failed_now > 0 )); then
    warn "$name: still has failures after retries: $failed_now tasks"
  else
    log "$name: all tasks done (OK or NR)"
  fi
}

# ----------------------------
# Generate task manifests (and register datasets)
# ----------------------------

gen_task_id() {
  local sweep="$1"; local fam="$2"; local N="$3"; local alpha="$4"; local t="$5"; local method="$6"; local variant="$7"
  echo "${sweep}__${fam}__n${N}__a$(sanitize_token "$alpha")__t${t}__${method}__${variant}"
}

log "Generating task manifests..."

TASK_A="$MANIFEST_DIR/tasks_A_alpha.tsv"
TASK_B_SMALL="$MANIFEST_DIR/tasks_B_N_small.tsv"
TASK_B_LARGE="$MANIFEST_DIR/tasks_B_N_large.tsv"
TASK_C_SMALL="$MANIFEST_DIR/tasks_C_t_small.tsv"
TASK_C_LARGE="$MANIFEST_DIR/tasks_C_t_large.tsv"
TASK_D="$MANIFEST_DIR/tasks_D_family.tsv"
TASK_D_SHAPE="$MANIFEST_DIR/tasks_D_shape.tsv"

TASK_HEADER="# task_id\tsweep\tfamily\tN\talpha\tt\tmethod\tvariant"

: > "$TASK_A";      echo -e "$TASK_HEADER" > "$TASK_A"
: > "$TASK_B_SMALL";echo -e "$TASK_HEADER" > "$TASK_B_SMALL"
: > "$TASK_B_LARGE";echo -e "$TASK_HEADER" > "$TASK_B_LARGE"
: > "$TASK_C_SMALL";echo -e "$TASK_HEADER" > "$TASK_C_SMALL"
: > "$TASK_C_LARGE";echo -e "$TASK_HEADER" > "$TASK_C_LARGE"
: > "$TASK_D";      echo -e "$TASK_HEADER" > "$TASK_D"
: > "$TASK_D_SHAPE";echo -e "$TASK_HEADER" > "$TASK_D_SHAPE"

# Compose ALPHAS_A
ALPHAS_A="$ALPHAS_A_BASE"
if [[ "$INCLUDE_ALPHA_CTRL" == "1" ]]; then
  ALPHAS_A="$ALPHAS_A $ALPHAS_A_CTRL"
fi

# (A*) alpha sweep, F0+F1
if [[ "$RUN_A" == "1" ]]; then
  for alpha in $ALPHAS_A; do
    for fam in F0 F1; do
      register_dataset "$fam" "$N_A" "$alpha"
      for mv in "${MODELS_MAIN[@]}"; do
        read -r method variant <<< "$mv"
        sweep="A_alpha"
        task_id="$(gen_task_id "$sweep" "$fam" "$N_A" "$alpha" "$T_A" "$method" "$variant")"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "$task_id" "$sweep" "$fam" "$N_A" "$alpha" "$T_A" "$method" "$variant" \
          >> "$TASK_A"
      done
    done
  done
fi

# (B*) N sweep, default F0; F1 optional
if [[ "$RUN_B" == "1" ]]; then
  families_b=("F0")
  if [[ "$INCLUDE_B_F1" == "1" ]]; then
    families_b+=("F1")
  fi

  for fam in "${families_b[@]}"; do
    for N in $NS_B; do
      register_dataset "$fam" "$N" "$ALPHA_B"
      for mv in "${MODELS_MAIN[@]}"; do
        read -r method variant <<< "$mv"
        sweep="B_N"
        task_id="$(gen_task_id "$sweep" "$fam" "$N" "$ALPHA_B" "$T_B" "$method" "$variant")"
        if (( N >= N_LARGE_THRESHOLD )); then
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$task_id" "$sweep" "$fam" "$N" "$ALPHA_B" "$T_B" "$method" "$variant" \
            >> "$TASK_B_LARGE"
        else
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$task_id" "$sweep" "$fam" "$N" "$ALPHA_B" "$T_B" "$method" "$variant" \
            >> "$TASK_B_SMALL"
        fi
      done
    done
  done
fi

# (C*) t sweep, default F0; F1 optional
if [[ "$RUN_C" == "1" ]]; then
  families_c=("F0")
  if [[ "$INCLUDE_C_F1" == "1" ]]; then
    families_c+=("F1")
  fi

  for fam in "${families_c[@]}"; do
    register_dataset "$fam" "$N_C" "$ALPHA_C"
    for t in $TS_C; do
      for mv in "${MODELS_MAIN[@]}"; do
        read -r method variant <<< "$mv"
        sweep="C_t"
        task_id="$(gen_task_id "$sweep" "$fam" "$N_C" "$ALPHA_C" "$t" "$method" "$variant")"
        if (( t > T_LARGE_THRESHOLD )); then
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$task_id" "$sweep" "$fam" "$N_C" "$ALPHA_C" "$t" "$method" "$variant" \
            >> "$TASK_C_LARGE"
        else
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$task_id" "$sweep" "$fam" "$N_C" "$ALPHA_C" "$t" "$method" "$variant" \
            >> "$TASK_C_SMALL"
        fi
      done
    done
  done
fi

# (D*) family sensitivity (fixed point), F0+F1
if [[ "$RUN_D" == "1" ]]; then
  for fam in F0 F1; do
    register_dataset "$fam" "$N_D" "$ALPHA_D"
    for mv in "${MODELS_MAIN[@]}"; do
      read -r method variant <<< "$mv"
      sweep="D_family"
      task_id="$(gen_task_id "$sweep" "$fam" "$N_D" "$ALPHA_D" "$T_D" "$method" "$variant")"
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$task_id" "$sweep" "$fam" "$N_D" "$ALPHA_D" "$T_D" "$method" "$variant" \
        >> "$TASK_D"
    done
  done
fi

# (D ext) shape_sigma sweep, lognormal + cv=1.0
if [[ "$RUN_D_SHAPE_SWEEP" == "1" ]]; then
  for sigma in $D_SHAPE_SIGMAS; do
    fam="F1sigma$(sanitize_token "$sigma")"
    # preset already created in RECTGEN_PRESET
    register_dataset "$fam" "$N_D" "$ALPHA_D"
    for mv in "${MODELS_MAIN[@]}"; do
      read -r method variant <<< "$mv"
      sweep="D_shape"
      task_id="$(gen_task_id "$sweep" "$fam" "$N_D" "$ALPHA_D" "$T_D" "$method" "$variant")"
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$task_id" "$sweep" "$fam" "$N_D" "$ALPHA_D" "$T_D" "$method" "$variant" \
        >> "$TASK_D_SHAPE"
    done
  done
fi

write_manifest_csv "$TASK_A" "$MANIFEST_DIR/manifest_A_alpha.csv"
write_manifest_csv "$TASK_B_SMALL" "$MANIFEST_DIR/manifest_B_N_small.csv"
write_manifest_csv "$TASK_B_LARGE" "$MANIFEST_DIR/manifest_B_N_large.csv"
write_manifest_csv "$TASK_C_SMALL" "$MANIFEST_DIR/manifest_C_t_small.csv"
write_manifest_csv "$TASK_C_LARGE" "$MANIFEST_DIR/manifest_C_t_large.csv"
write_manifest_csv "$TASK_D" "$MANIFEST_DIR/manifest_D_family.csv"
write_manifest_csv "$TASK_D_SHAPE" "$MANIFEST_DIR/manifest_D_shape.csv"

CNT_A="$(task_count "$TASK_A")"
CNT_BS="$(task_count "$TASK_B_SMALL")"
CNT_BL="$(task_count "$TASK_B_LARGE")"
CNT_CS="$(task_count "$TASK_C_SMALL")"
CNT_CL="$(task_count "$TASK_C_LARGE")"
CNT_D="$(task_count "$TASK_D")"
CNT_DS="$(task_count "$TASK_D_SHAPE")"
TOTAL_TASKS=$(( CNT_A + CNT_BS + CNT_BL + CNT_CS + CNT_CL + CNT_D + CNT_DS ))

log "Task counts: A=$CNT_A, B_small=$CNT_BS, B_large=$CNT_BL, C_small=$CNT_CS, C_large=$CNT_CL, D=$CNT_D, D_shape=$CNT_DS (total=$TOTAL_TASKS)"
if (( TOTAL_TASKS == 0 )); then
  warn "No tasks generated. Check RUN_* flags and grids."
  exit 2
fi

# ----------------------------
# Pre-generate all required datasets (once)
# ----------------------------

log "Pre-generating datasets (unique count=${#DATASET_ORDER[@]})..."
for name in "${DATASET_ORDER[@]}"; do
  IFS=$'\t' read -r fam N alpha <<< "${DATASET_SPEC[$name]}"
  if ! ensure_dataset "$fam" "$N" "$alpha"; then
    warn "Dataset generation failed for $name. Inspect: $DATA_LOG_ROOT/${name}.log"
    exit 3
  fi
done
log "All datasets are ready."

# ----------------------------
# Execute sweeps
# ----------------------------

if (( CNT_A > 0 )); then
  run_with_retries "$TASK_A" "$PAR_A" "A_alpha"
fi
if (( CNT_BS > 0 )); then
  run_with_retries "$TASK_B_SMALL" "$PAR_B_SMALL" "B_N_small"
fi
if (( CNT_BL > 0 )); then
  run_with_retries "$TASK_B_LARGE" "$PAR_B_LARGE" "B_N_large"
fi
if (( CNT_CS > 0 )); then
  run_with_retries "$TASK_C_SMALL" "$PAR_C_SMALL" "C_t_small"
fi
if (( CNT_CL > 0 )); then
  run_with_retries "$TASK_C_LARGE" "$PAR_C_LARGE" "C_t_large"
fi
if (( CNT_D > 0 )); then
  run_with_retries "$TASK_D" "$PAR_D" "D_family"
fi
if (( CNT_DS > 0 )); then
  run_with_retries "$TASK_D_SHAPE" "$PAR_D" "D_shape"
fi

# ----------------------------
# Post: merge CSVs + summarize status
# ----------------------------

MERGED_DIR="$MANIFEST_DIR/merged"
mkdir -p "$MERGED_DIR"

while IFS= read -r d; do
  sweep="$(basename "$d")"
  merge_csv_dir "$d" "$MERGED_DIR/${sweep}_merged.csv"
done < <(find "$OUT_ROOT" -mindepth 1 -maxdepth 1 -type d | sort)

GLOBAL_FAIL="$MANIFEST_DIR/FAILURES.tsv"
GLOBAL_NR="$MANIFEST_DIR/NOT_RUNNABLE.tsv"

{
  echo -e "$TASK_HEADER"
  for tf in "$TASK_A" "$TASK_B_SMALL" "$TASK_B_LARGE" "$TASK_C_SMALL" "$TASK_C_LARGE" "$TASK_D" "$TASK_D_SHAPE"; do
    while IFS=$'\t' read -r task_id sweep fam N alpha t method variant; do
      [[ -z "$task_id" ]] && continue
      [[ "$task_id" =~ ^# ]] && continue
      if [[ ! -f "$STATUS_ROOT/${task_id}.ok" && ! -f "$STATUS_ROOT/${task_id}.nr" ]]; then
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "$task_id" "$sweep" "$fam" "$N" "$alpha" "$t" "$method" "$variant"
      fi
    done < "$tf"
  done
} > "$GLOBAL_FAIL"

{
  echo -e "$TASK_HEADER"
  for tf in "$TASK_A" "$TASK_B_SMALL" "$TASK_B_LARGE" "$TASK_C_SMALL" "$TASK_C_LARGE" "$TASK_D" "$TASK_D_SHAPE"; do
    while IFS=$'\t' read -r task_id sweep fam N alpha t method variant; do
      [[ -z "$task_id" ]] && continue
      [[ "$task_id" =~ ^# ]] && continue
      if [[ -f "$STATUS_ROOT/${task_id}.nr" ]]; then
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "$task_id" "$sweep" "$fam" "$N" "$alpha" "$t" "$method" "$variant"
      fi
    done < "$tf"
  done
} > "$GLOBAL_NR"

FAIL_N=0
if [[ -s "$GLOBAL_FAIL" ]]; then
  FAIL_N=$(( $(wc -l < "$GLOBAL_FAIL" | tr -d '[:space:]') - 1 ))
  (( FAIL_N < 0 )) && FAIL_N=0
fi
NR_N=0
if [[ -s "$GLOBAL_NR" ]]; then
  NR_N=$(( $(wc -l < "$GLOBAL_NR" | tr -d '[:space:]') - 1 ))
  (( NR_N < 0 )) && NR_N=0
fi

# ----------------------------
# Post: summary statistics (mean/median/CV + throughput for C*)
# ----------------------------

SUMMARY_DIR="$MANIFEST_DIR/summary"
mkdir -p "$SUMMARY_DIR"

"$PYTHON" - "$MERGED_DIR" "$SUMMARY_DIR" <<'PY'
import csv
import json
import os
import statistics
import sys
from collections import defaultdict
from pathlib import Path

merged_dir = Path(sys.argv[1])
out_dir = Path(sys.argv[2])
out_dir.mkdir(parents=True, exist_ok=True)

# metrics we care about (phases_json keys)
PHASE_KEYS = [
    "run_build_ms", "run_count_ms", "run_sample_ms",
    # SJS-Sampling fine-grained
    "build_events_ms", "build_y_domain_ms", "build_ranks_ms", "build_active_indices_ms",
    "phase1_count_ms", "phase2_plan_ms", "phase3_sample_ms",
    # KD-tree fine-grained
    "build_ms", "build_points_ms", "phase3_fill_ms",
    # EnumSampling fine-grained
    "phase1_enumerate_materialize_ms", "phase2_resample_ms",
]

# group key -> list of rows (ok=1 only)
# key: (sweep, family, N, alpha, t, method, variant)
rows_ok = defaultdict(list)
rows_all = defaultdict(list)

for csv_path in sorted(merged_dir.glob("*_merged.csv")):
    sweep = csv_path.name.replace("_merged.csv", "")
    with csv_path.open(newline='') as f:
        reader = csv.DictReader(f)
        for r in reader:
            # infer family/N/alpha from dataset name: d2_<fam>_n<N>_a<alpha>_gs<seed>
            ds = r.get('dataset','')
            fam = ''
            N = ''
            alpha = ''
            try:
                parts = ds.split('_')
                # ['d2', fam, 'n...', 'a...', 'gs...']
                if len(parts) >= 5 and parts[0] == 'd2':
                    fam = parts[1]
                    N = parts[2][1:]
                    alpha = parts[3][1:]
            except Exception:
                pass

            method = r.get('method','')
            variant = r.get('variant','')
            t = r.get('t','')

            key = (sweep, fam, N, alpha, t, method, variant)
            rows_all[key].append(r)
            if r.get('ok','0') == '1':
                rows_ok[key].append(r)

# helper: safe float

def f64(x):
    try:
        return float(x)
    except Exception:
        return None

# parse phases_json

def parse_phases(s):
    if not s:
        return {}
    try:
        return json.loads(s)
    except Exception:
        return {}

# summarize list -> mean/median/cv

def summarize(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None, None, None, 0
    mean = statistics.fmean(vals)
    med = statistics.median(vals)
    if len(vals) >= 2:
        stdev = statistics.pstdev(vals)  # population stdev
        cv = (stdev / mean) if mean != 0 else None
    else:
        cv = None
    return mean, med, cv, len(vals)

out_path = out_dir / 'summary_stats.csv'
with out_path.open('w', newline='') as f:
    fieldnames = [
        'sweep','family','N','alpha','t','method','variant',
        'ok_reps','total_reps','status',
        'wall_ms_mean','wall_ms_median','wall_ms_cv',
        # phases
    ]
    for k in PHASE_KEYS:
        fieldnames += [f'{k}_mean', f'{k}_median', f'{k}_cv']

    # derived throughput (based on run_sample_ms)
    fieldnames += ['samples_per_sec_mean','samples_per_sec_median','samples_per_sec_cv',
                   'ns_per_sample_mean','ns_per_sample_median','ns_per_sample_cv']

    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()

    for key in sorted(rows_all.keys()):
        sweep, fam, N, alpha, t, method, variant = key
        ok_list = rows_ok.get(key, [])
        all_list = rows_all.get(key, [])
        ok_reps = len(ok_list)
        total_reps = len(all_list)

        status = 'OK' if ok_reps == total_reps and total_reps > 0 else ('NR/FAIL' if ok_reps == 0 else 'PARTIAL')

        wall_vals = [f64(r.get('wall_ms','')) for r in ok_list]
        wall_mean, wall_med, wall_cv, _ = summarize(wall_vals)

        out_row = {
            'sweep': sweep,
            'family': fam,
            'N': N,
            'alpha': alpha,
            't': t,
            'method': method,
            'variant': variant,
            'ok_reps': ok_reps,
            'total_reps': total_reps,
            'status': status,
            'wall_ms_mean': wall_mean,
            'wall_ms_median': wall_med,
            'wall_ms_cv': wall_cv,
        }

        # collect phase metrics
        phase_vals = {k: [] for k in PHASE_KEYS}
        run_sample_ms_vals = []
        t_int = None
        try:
            t_int = int(float(t))
        except Exception:
            t_int = None

        for r in ok_list:
            ph = parse_phases(r.get('phases_json',''))
            for k in PHASE_KEYS:
                phase_vals[k].append(f64(ph.get(k,'')))
            rs = f64(ph.get('run_sample_ms',''))
            run_sample_ms_vals.append(rs)

        for k in PHASE_KEYS:
            m, md, cv, _ = summarize(phase_vals[k])
            out_row[f'{k}_mean'] = m
            out_row[f'{k}_median'] = md
            out_row[f'{k}_cv'] = cv

        # Throughput metrics based on run_sample_ms
        if t_int is not None and t_int > 0:
            sps = []
            nsp = []
            for rs in run_sample_ms_vals:
                if rs is None or rs <= 0:
                    continue
                sps.append(1000.0 * t_int / rs)
                nsp.append(1e6 * rs / t_int)
            m, md, cv, _ = summarize(sps)
            out_row['samples_per_sec_mean'] = m
            out_row['samples_per_sec_median'] = md
            out_row['samples_per_sec_cv'] = cv
            m, md, cv, _ = summarize(nsp)
            out_row['ns_per_sample_mean'] = m
            out_row['ns_per_sample_median'] = md
            out_row['ns_per_sample_cv'] = cv

        w.writerow(out_row)

# also write a compact NR/FAIL list (ok_reps==0)
fail_path = out_dir / 'summary_not_ok.csv'
with fail_path.open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['sweep','family','N','alpha','t','method','variant','total_reps'])
    for key in sorted(rows_all.keys()):
        if len(rows_ok.get(key, [])) == 0:
            sweep, fam, N, alpha, t, method, variant = key
            w.writerow([sweep, fam, N, alpha, t, method, variant, len(rows_all[key])])

print(f"[summary] wrote: {out_path}")
print(f"[summary] wrote: {fail_path}")
PY

# ----------------------------
# Write RUN_SUMMARY
# ----------------------------

{
  echo "# legacy EXP-2D runner output (align to 实验计划.md; TH/throughput profile)"
  echo "timestamp=$(date -Is)"
  echo "exp_tag=$EXP_TAG"
  echo "build_type=$BUILD_TYPE"
  echo "threads=$THREADS"
  echo "repeats=$REPEATS"
  echo "seed=$SEED"
  echo "write_samples=$WRITE_SAMPLES"
  echo "verify=$VERIFY"
  echo "python=$PYTHON"
  echo "gen_seed=$GEN_SEED"
  echo "audit_pairs=$AUDIT_PAIRS"
  echo "audit_seed=$AUDIT_SEED"
  echo "rectgen_base_script=$RECTGEN_BASE_SCRIPT"
  echo "alacarte_module=$ALACARTE_MODULE"
  echo "F0=$F0_VOLUME_DIST,$F0_VOLUME_CV,$F0_SHAPE_SIGMA"
  echo "F1=$F1_VOLUME_DIST,$F1_VOLUME_CV,$F1_SHAPE_SIGMA"
  echo "enum_cap=$ENUM_CAP"
  echo "enum_precheck=$ENUM_PRECHECK"
  echo "run_signature=$RUN_SIGNATURE"
  echo "run_log_file=$RUN_LOG_FILE"
  echo "result_copy_datasets=$RESULT_COPY_DATASETS"
  echo "prune_temp_datasets=$PRUNE_TEMP_DATASETS"
  echo
  echo "RUN_A=$RUN_A (ALPHAS_A=$ALPHAS_A N_A=$N_A T_A=$T_A INCLUDE_ALPHA_CTRL=$INCLUDE_ALPHA_CTRL)"
  echo "RUN_B=$RUN_B (NS_B=$NS_B ALPHA_B=$ALPHA_B T_B=$T_B INCLUDE_B_F1=$INCLUDE_B_F1)"
  echo "RUN_C=$RUN_C (TS_C=$TS_C N_C=$N_C ALPHA_C=$ALPHA_C INCLUDE_C_F1=$INCLUDE_C_F1)"
  echo "RUN_D=$RUN_D (N_D=$N_D ALPHA_D=$ALPHA_D T_D=$T_D)"
  echo "RUN_D_SHAPE_SWEEP=$RUN_D_SHAPE_SWEEP (D_SHAPE_SIGMAS=$D_SHAPE_SIGMAS)"
  echo
  echo "Parallelism: max=$MAX_PARALLEL (A=$PAR_A, B_small=$PAR_B_SMALL, B_large=$PAR_B_LARGE, C_small=$PAR_C_SMALL, C_large=$PAR_C_LARGE, D=$PAR_D)"
  echo "Not runnable tasks (NR)=$NR_N"
  echo "Failed tasks (FAIL)=$FAIL_N"
  echo
  echo "Outputs (temp workspace):"
  echo "  - Root run log    : $RUN_LOG_FILE"
  echo "  - Per-task logs   : $LOG_ROOT/<task_id>.log"
  echo "  - Per-task CSV    : $OUT_ROOT/<sweep>/<task_id>/run.csv"
  echo "  - Merged CSVs     : $MERGED_DIR/*_merged.csv"
  echo "  - Summary stats   : $SUMMARY_DIR/summary_stats.csv"
  echo "  - NR list         : $GLOBAL_NR"
  echo "  - Failures list   : $GLOBAL_FAIL"
  echo "  - Dataset logs    : $DATA_LOG_ROOT/*.log"
} > "$MANIFEST_DIR/RUN_SUMMARY.txt"

TMP_DATASET_SIZE="$(dir_size_human "$DATA_ROOT")"
TMP_OUT_SIZE="$(dir_size_human "$OUT_ROOT")"
TMP_TASK_LOG_SIZE="$(dir_size_human "$LOG_ROOT")"
TMP_DATASET_LOG_SIZE="$(dir_size_human "$DATA_LOG_ROOT")"

rm -rf "$RESULT_ROOT"
mkdir -p "$RESULT_ROOT"

copy_dir_if_exists "$MANIFEST_DIR" "$RESULT_ROOT/manifest"
copy_dir_if_exists "$OUT_ROOT" "$RESULT_ROOT/out"
copy_dir_if_exists "$LOG_ROOT" "$RESULT_ROOT/logs"
copy_dir_if_exists "$DATA_LOG_ROOT" "$RESULT_ROOT/dataset_logs"
copy_dir_if_exists "$STATUS_ROOT" "$RESULT_ROOT/status"
copy_dir_if_exists "$DATA_STATUS_ROOT" "$RESULT_ROOT/dataset_status"
copy_dir_if_exists "$PRESET_DIR" "$RESULT_ROOT/rectgen_presets"
cp -a "$RUN_LOG_FILE" "$RESULT_ROOT/run.log"

if [[ "$RESULT_COPY_DATASETS" == "1" ]]; then
  copy_dir_if_exists "$DATA_ROOT" "$RESULT_ROOT/datasets"
else
  log "Skipping dataset binaries in result snapshot (RESULT_COPY_DATASETS=0)."
fi

if [[ "$PRUNE_TEMP_DATASETS" == "1" ]]; then
  log "PRUNE_TEMP_DATASETS=1: removing temporary dataset binaries under $DATA_ROOT"
  rm -rf "$DATA_ROOT"
fi

RESULT_SIZE="$(dir_size_human "$RESULT_ROOT")"

log "Synced curated results to: $RESULT_ROOT"
log "Space usage: temp datasets=$TMP_DATASET_SIZE, temp out=$TMP_OUT_SIZE, temp task_logs=$TMP_TASK_LOG_SIZE, temp dataset_logs=$TMP_DATASET_LOG_SIZE, result_total=$RESULT_SIZE"
log "Merged CSVs under: $RESULT_ROOT/manifest/merged/"
log "Summary under    : $RESULT_ROOT/manifest/summary/"
log "Root run log     : $RUN_LOG_FILE"

if (( FAIL_N > 0 )); then
  warn "There are still $FAIL_N failed tasks after retries. See: $RESULT_ROOT/manifest/FAILURES.tsv"
  exit 1
fi

log "All tasks completed (OK + NR)."
