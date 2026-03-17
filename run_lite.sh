#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$' \t\n'

# run_lite.sh
#
# 轻量级 alpha≈100 全组评测脚本（默认对齐上一轮建议的主评测网格）
# ----------------------------------------------------------------------
# 默认会自动：
#   1) 配置并构建 build/sjs_sweep
#   2) 为 d=2/3/4/5 生成 4 份 sweep JSON
#   3) 逐维运行 sjs_sweep
#   4) 合并所有维度的 raw/summary CSV
#
# 默认评测网格：
#   - methods      : ours/sampling, kd_tree/sampling
#   - alpha        : 100
#   - dim          : 2 / 3 / 4 / 5
#   - n_r = n_s    : d2=10000, d3=2000, d4=500, d5=200
#   - t            : 10000, 50000, 100000
#   - gen_seed     : 7, 17
#   - run seed     : 1, 2, 3
#   - threads      : 1
#   - write_samples: false
#   - audit_pairs  : 20000
#
# 用法：
#   chmod +x run_lite.sh
#   ./run_lite.sh
#
# 常用覆盖（环境变量）：
#   BUILD_DIR=build_rel
#   BUILD_TYPE=Release
#   JOBS=8
#   CLEAN_BUILD=0|1
#   CLEAN_OUT=0|1
#   OUT_ROOT=out/alpha100_light
#   LOG_LEVEL=warn|info
#   THREADS=1
#   ALPHA=100
#   AUDIT_PAIRS=20000
#   T_LIST="10000 50000 100000"
#   GEN_SEEDS="7 17"
#   RUN_SEEDS="1 2 3"
#   METHODS="ours kd_tree"
#   VARIANTS="sampling"
#   N_D2=10000 N_D3=2000 N_D4=500 N_D5=200
#
# 可选：切到更轻的 smoke 版
#   PROFILE=smoke ./run_lite.sh
#   # 等价于：GEN_SEEDS="7"，T_LIST="50000 100000"
#
# 仅生成配置、不实际执行：
#   DRY_RUN=1 SKIP_BUILD=1 ./run_lite.sh

script_dir() {
  cd "$(dirname "${BASH_SOURCE[0]}")" && pwd
}
ROOT="$(script_dir)"

resolve_path() {
  local p="$1"
  if [[ -z "$p" ]]; then
    printf '%s' "$ROOT"
  elif [[ "$p" = /* ]]; then
    printf '%s' "$p"
  else
    printf '%s/%s' "$ROOT" "${p#./}"
  fi
}

timestamp_now() { date -Is; }
log()  { printf '[%s][run_lite][INFO] %s\n'  "$(timestamp_now)" "$*"; }
warn() { printf '[%s][run_lite][WARN] %s\n'  "$(timestamp_now)" "$*" >&2; }
die()  { printf '[%s][run_lite][FATAL] %s\n' "$(timestamp_now)" "$*" >&2; exit 1; }

print_usage() {
  sed -n '1,120p' "$0" | sed 's/^# \{0,1\}//'
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  print_usage
  exit 0
fi

nproc_safe() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
  elif command -v getconf >/dev/null 2>&1; then
    getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4
  else
    echo 4
  fi
}

BUILD_DIR="$(resolve_path "${BUILD_DIR:-build}")"
BUILD_TYPE="${BUILD_TYPE:-Release}"
JOBS="${JOBS:-$(nproc_safe)}"
OUT_ROOT="$(resolve_path "${OUT_ROOT:-out/alpha100_light}")"
CONFIG_DIR="$(resolve_path "${CONFIG_DIR:-$OUT_ROOT/_configs}")"
LOG_LEVEL="${LOG_LEVEL:-warn}"
THREADS="${THREADS:-1}"
ALPHA="${ALPHA:-100}"
AUDIT_PAIRS="${AUDIT_PAIRS:-20000}"
METHODS="${METHODS:-ours kd_tree}"
VARIANTS="${VARIANTS:-sampling}"
GEN_SEEDS="${GEN_SEEDS:-7 17}"
RUN_SEEDS="${RUN_SEEDS:-1 2 3}"
T_LIST="${T_LIST:-10000 50000 100000}"
PROFILE="${PROFILE:-main}"
CLEAN_BUILD="${CLEAN_BUILD:-0}"
CLEAN_OUT="${CLEAN_OUT:-0}"
SKIP_BUILD="${SKIP_BUILD:-0}"
DRY_RUN="${DRY_RUN:-0}"

N_D2="${N_D2:-10000}"
N_D3="${N_D3:-2000}"
N_D4="${N_D4:-500}"
N_D5="${N_D5:-200}"

if [[ "$PROFILE" == "smoke" ]]; then
  [[ "$GEN_SEEDS" == "7 17" ]] && GEN_SEEDS="7"
  [[ "$T_LIST" == "10000 50000 100000" ]] && T_LIST="50000 100000"
fi

json_num_array_from_words() {
  local s="$1"
  local -a arr=()
  read -r -a arr <<< "$s"
  local out="["
  local first=1
  local x
  for x in "${arr[@]}"; do
    [[ -z "$x" ]] && continue
    if (( first )); then
      out+="$x"
      first=0
    else
      out+=", $x"
    fi
  done
  out+="]"
  printf '%s' "$out"
}

json_str_array_from_words() {
  local s="$1"
  local -a arr=()
  read -r -a arr <<< "$s"
  local out="["
  local first=1
  local x escaped
  for x in "${arr[@]}"; do
    [[ -z "$x" ]] && continue
    escaped="${x//\\/\\\\}"
    escaped="${escaped//\"/\\\"}"
    if (( first )); then
      out+="\"$escaped\""
      first=0
    else
      out+=", \"$escaped\""
    fi
  done
  out+="]"
  printf '%s' "$out"
}

n_for_dim() {
  case "$1" in
    2) printf '%s' "$N_D2" ;;
    3) printf '%s' "$N_D3" ;;
    4) printf '%s' "$N_D4" ;;
    5) printf '%s' "$N_D5" ;;
    *) return 1 ;;
  esac
}

cmake_cache_internal_value() {
  local cache="$1"
  local key="$2"
  [[ -f "$cache" ]] || return 1
  local line
  line="$(grep -E "^${key}:INTERNAL=" "$cache" | head -n 1 || true)"
  [[ -n "$line" ]] || return 1
  printf '%s' "${line#*=}"
}

ensure_build_dir() {
  local cache="$BUILD_DIR/CMakeCache.txt"

  if [[ "$CLEAN_BUILD" == "1" ]]; then
    log "Removing build dir: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
  fi

  if [[ ! -f "$cache" ]]; then
    log "Configuring CMake in $BUILD_DIR"
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    return 0
  fi

  local cache_home cache_dir
  cache_home="$(cmake_cache_internal_value "$cache" "CMAKE_HOME_DIRECTORY" || true)"
  cache_dir="$(cmake_cache_internal_value "$cache" "CMAKE_CACHEFILE_DIR" || true)"

  if [[ -n "$cache_home" && "$cache_home" != "$ROOT" ]]; then
    warn "Stale source path in CMakeCache: $cache_home"
    rm -rf "$BUILD_DIR"
    log "Reconfiguring CMake in $BUILD_DIR"
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    return 0
  fi

  if [[ -n "$cache_dir" && "$cache_dir" != "$BUILD_DIR" ]]; then
    warn "Stale build path in CMakeCache: $cache_dir"
    rm -rf "$BUILD_DIR"
    log "Reconfiguring CMake in $BUILD_DIR"
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    return 0
  fi
}

find_sweep_app() {
  local -a candidates=(
    "${SWEEP_APP:-}"
    "$BUILD_DIR/sjs_sweep"
    "$BUILD_DIR/apps/sjs_sweep"
  )
  local c
  for c in "${candidates[@]}"; do
    [[ -z "$c" ]] && continue
    if [[ -x "$c" ]]; then
      printf '%s' "$c"
      return 0
    fi
  done
  local found
  found="$(find "$BUILD_DIR" -type f -name sjs_sweep -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "$found" ]] || return 1
  printf '%s' "$found"
}

build_sweep_app() {
  if [[ "$SKIP_BUILD" == "1" ]]; then
    return 0
  fi
  ensure_build_dir
  log "Building target sjs_sweep (BUILD_TYPE=$BUILD_TYPE, JOBS=$JOBS)"
  cmake --build "$BUILD_DIR" --target sjs_sweep -j "$JOBS"
}

write_sweep_config() {
  local dim="$1"
  local n="$2"
  local cfg_path="$3"
  local dim_out_dir="$4"

  mkdir -p "$(dirname "$cfg_path")" "$dim_out_dir"

  cat > "$cfg_path" <<JSON
{
  "base": {
    "dataset": {
      "source": "synthetic",
      "name": "alpha100_lite_d${dim}",
      "dim": ${dim},
      "synthetic": {
        "generator": "rectgen",
        "n_r": ${n},
        "n_s": ${n},
        "alpha": ${ALPHA},
        "seed": 7,
        "params": {
          "audit_pairs": ${AUDIT_PAIRS}
        }
      }
    },
    "run": {
      "method": "ours",
      "variant": "sampling",
      "t": 100000,
      "seed": 1,
      "repeats": 3,
      "write_samples": false
    },
    "output": {
      "out_dir": "${dim_out_dir}"
    },
    "logging": {
      "level": "${LOG_LEVEL}"
    },
    "sys": {
      "threads": ${THREADS}
    }
  },
  "sweep": {
    "alpha": [${ALPHA}],
    "gen_seed": $(json_num_array_from_words "$GEN_SEEDS"),
    "t": $(json_num_array_from_words "$T_LIST"),
    "method": $(json_str_array_from_words "$METHODS"),
    "variant": $(json_str_array_from_words "$VARIANTS"),
    "seed": $(json_num_array_from_words "$RUN_SEEDS")
  },
  "files": {
    "raw": "sweep_raw.csv",
    "summary": "sweep_summary.csv"
  }
}
JSON
}

merge_csvs() {
  local out_csv="$1"
  shift
  local first=1
  : > "$out_csv"
  local f
  for f in "$@"; do
    [[ -f "$f" ]] || continue
    if (( first )); then
      cat "$f" > "$out_csv"
      first=0
    else
      tail -n +2 "$f" >> "$out_csv"
    fi
  done
  if (( first )); then
    rm -f "$out_csv"
    return 1
  fi
  return 0
}

mkdir -p "$CONFIG_DIR"
if [[ "$CLEAN_OUT" == "1" ]]; then
  log "Removing output dir: $OUT_ROOT"
  rm -rf "$OUT_ROOT"
fi
mkdir -p "$OUT_ROOT"

build_sweep_app
SWEEP_APP="$(find_sweep_app || true)"
[[ -n "$SWEEP_APP" || "$DRY_RUN" == "1" ]] || die "Cannot find executable sjs_sweep under $BUILD_DIR"

manifest="$OUT_ROOT/manifest.tsv"
printf 'dim\tn\talpha\tt_list\tgen_seeds\trun_seeds\tconfig\tout_dir\n' > "$manifest"

summary_files=()
raw_files=()

for dim in 2 3 4 5; do
  n="$(n_for_dim "$dim")" || die "Unsupported dim: $dim"
  dim_out_dir="$OUT_ROOT/d${dim}"
  cfg_path="$CONFIG_DIR/alpha100_lite_d${dim}.json"

  write_sweep_config "$dim" "$n" "$cfg_path" "$dim_out_dir"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$dim" "$n" "$ALPHA" "$T_LIST" "$GEN_SEEDS" "$RUN_SEEDS" "$cfg_path" "$dim_out_dir" \
    >> "$manifest"

  summary_files+=("$dim_out_dir/sweep_summary.csv")
  raw_files+=("$dim_out_dir/sweep_raw.csv")

  if [[ "$DRY_RUN" == "1" ]]; then
    log "[DRY_RUN] $SWEEP_APP --config=$cfg_path"
  else
    log "Running d=$dim, n_r=n_s=$n"
    "$SWEEP_APP" --config="$cfg_path"
  fi
done

if [[ "$DRY_RUN" != "1" ]]; then
  if merge_csvs "$OUT_ROOT/all_dims_summary.csv" "${summary_files[@]}"; then
    log "Merged summary -> $OUT_ROOT/all_dims_summary.csv"
  else
    warn "No summary CSV merged"
  fi

  if merge_csvs "$OUT_ROOT/all_dims_raw.csv" "${raw_files[@]}"; then
    log "Merged raw -> $OUT_ROOT/all_dims_raw.csv"
  else
    warn "No raw CSV merged"
  fi
fi

cat <<EOF2

[DONE] run_lite.sh finished.

Root:      $ROOT
Build dir: $BUILD_DIR
Sweep app: ${SWEEP_APP:-<skipped in dry-run>}
Out root:  $OUT_ROOT
Manifest:  $manifest

Per-dimension outputs:
  - d2: $OUT_ROOT/d2/sweep_summary.csv
  - d3: $OUT_ROOT/d3/sweep_summary.csv
  - d4: $OUT_ROOT/d4/sweep_summary.csv
  - d5: $OUT_ROOT/d5/sweep_summary.csv

Merged outputs:
  - $OUT_ROOT/all_dims_summary.csv
  - $OUT_ROOT/all_dims_raw.csv

Key defaults:
  alpha      = $ALPHA
  gen_seed   = $GEN_SEEDS
  run_seed   = $RUN_SEEDS
  t          = $T_LIST
  methods    = $METHODS
  variants   = $VARIANTS
  threads    = $THREADS
  n(d=2..5)  = $N_D2 / $N_D3 / $N_D4 / $N_D5

EOF2
