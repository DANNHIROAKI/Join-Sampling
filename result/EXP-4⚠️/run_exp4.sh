#!/usr/bin/env bash
# ==============================================================================
# EXP-4: Scalability vs N (RQ4)
#
# Aligned with EXP-4_aligned.md.
#
# Integrity note:
#   This script aims to be FAIR and REPRODUCIBLE:
#   - same predicate (half-open), same resources (threads), same reporting (wall_ms).
#   - no baseline is hidden by default.
#
# Profiles:
#   RUN_PROFILE=balanced  -> alpha5 + alpha100, t=1e5, repeats=5 (p95-friendly)
#   RUN_PROFILE=ours      -> alpha5 + alpha200, t=1e6, repeats=5 (stress/high-sample)
#
# Optional:
#   KEEP_HISTORY=1 saves a copy to results/raw/exp4_runs/<RUN_ID>/ (in addition to overwriting results/raw/exp4).
#
# Quick usage:
#   chmod +x run/run_exp4.sh
#   RUN_PROFILE=balanced bash run/run_exp4.sh
#   RUN_PROFILE=ours     bash run/run_exp4.sh
#
# Common overrides:
#   REGIMES="alpha5:5 alpha200:200"      bash run/run_exp4.sh
#   T=1000000 T_LIST=""                  bash run/run_exp4.sh
#   VARIANTS="sampling"                  bash run/run_exp4.sh
#   METHODS="ours pbsm kd_tree"          bash run/run_exp4.sh
#   EXTRA_RUN_ARGS="--enum_cap=50000000" bash run/run_exp4.sh
#   TIMEOUT_SEC=3600                     bash run/run_exp4.sh
# ==============================================================================

set -euo pipefail
IFS=$' \n\t'

trap 'echo -e "[EXP4][FATAL] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

# --------------------------------------
# helpers
# --------------------------------------
msg() { echo "[$(date -Is)] [EXP4] $*"; }
warn() { echo "[$(date -Is)] [EXP4][WARN] $*" >&2; }
die() { echo "[$(date -Is)] [EXP4][FATAL] $*" >&2; exit 1; }

need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

sanitize_token() { echo "$1" | sed 's/[^A-Za-z0-9_.-]/_/g'; }

detect_jobs() {
  if [[ -n "${JOBS:-}" ]]; then echo "${JOBS}"; return; fi
  if command -v nproc >/dev/null 2>&1; then nproc; return; fi
  if [[ "$(uname -s)" == "Darwin" ]] && command -v sysctl >/dev/null 2>&1; then
    sysctl -n hw.ncpu 2>/dev/null || echo 4; return
  fi
  echo 4
}

build_subdir_for_type() {
  case "$1" in
    Release) echo "release" ;;
    Debug) echo "debug" ;;
    RelWithDebInfo) echo "relwithdebinfo" ;;
    MinSizeRel) echo "minsizerel" ;;
    *) echo "$(echo "$1" | tr '[:upper:]' '[:lower:]')" ;;
  esac
}

resolve_exe() {
  local name="$1"; local root="$2"
  local candidates=("$root/$name" "$root/apps/$name" "$root/bin/$name" "$root/src/apps/$name")
  for p in "${candidates[@]}"; do [[ -x "$p" ]] && { echo "$p"; return 0; }; done
  local found
  found="$(find "$root" -maxdepth 4 -type f -name "$name" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "$found" && -x "$found" ]] && { echo "$found"; return 0; }
  return 1
}

parse_rss_kb_from_timev() {
  local mem_log="$1"
  [[ -f "$mem_log" ]] || { echo ""; return 0; }
  awk -F: '/Maximum resident set size/ {gsub(/^[ \t]+/,"",$2); print $2}' "$mem_log" | tail -n 1 | tr -d '\r'
}

# FIXED: does not "exit" early (because AWK exit still runs END blocks).
parse_enum_truncated_any_from_run_csv() {
  local run_csv="$1"
  [[ -f "$run_csv" ]] || { echo ""; return 0; }
  awk -F, '
    NR==1 {
      for (i=1;i<=NF;i++) if ($i=="enum_truncated") col=i;
      next
    }
    {
      if (col>0 && $col+0==1) found=1
    }
    END {
      if (col>0) print (found?1:0)
    }
  ' "$run_csv" 2>/dev/null | tail -n 1
}

run_with_optional_timeout() {
  local timeout_sec="$1"; shift
  if [[ "${timeout_sec}" != "0" && "${timeout_sec}" != "" ]] && command -v timeout >/dev/null 2>&1; then
    timeout "${timeout_sec}" "$@"; return $?
  fi
  "$@"
}

# --------------------------------------
# locate repo root
# --------------------------------------
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --------------------------------------
# Standard directories
# --------------------------------------
BUILD_ROOT="${ROOT_DIR}/build"
TEMP_ROOT="${ROOT_DIR}/run/temp/exp4"
RESULT_ROOT="${ROOT_DIR}/results/raw/exp4"
DATA_ROOT="${ROOT_DIR}/data/synthetic/exp4"

# Optional history root
KEEP_HISTORY="${KEEP_HISTORY:-0}"
HISTORY_ROOT="${HISTORY_ROOT:-${ROOT_DIR}/results/raw/exp4_runs}"

# --------------------------------------
# RUN_PROFILE presets (defaults only)
# --------------------------------------
RUN_PROFILE="${RUN_PROFILE:-balanced}"

N_LIST_DEFAULT="50000 100000 200000 400000 800000"

# Make p95 meaningful by default (EXP-4.md recommends 5 when using p95).
REPEATS_DEFAULT="5"

T_DEFAULT="100000"
T_LIST_DEFAULT=""                 # if empty, use T_DEFAULT
DEFAULT_REGIMES_SPEC="alpha5:5 alpha100:100"

# Ours "stress/high-sample" preset:
# - higher alpha (200) is explicitly allowed as a clearer crossover check
# - t=1e6 makes sampling stage visible
if [[ "$RUN_PROFILE" == "ours" || "$RUN_PROFILE" == "stress" ]]; then
  DEFAULT_REGIMES_SPEC="alpha5:5 alpha200:200"
  T_LIST_DEFAULT="1000000"        # single-t by default; override with T_LIST="100000 1000000" if needed
  REPEATS_DEFAULT="5"
fi

# --------------------------------------
# Defaults (env overrides always win)
# --------------------------------------
N_LIST="${N_LIST:-$N_LIST_DEFAULT}"
T="${T:-$T_DEFAULT}"
T_LIST="${T_LIST:-$T_LIST_DEFAULT}"
REPEATS="${REPEATS:-$REPEATS_DEFAULT}"

SEED="${SEED:-1}"
GEN_SEED="${GEN_SEED:-1}"
THREADS="${THREADS:-1}"

METHODS="${METHODS:-"ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection tsunami"}"
VARIANTS="${VARIANTS:-"sampling enum_sampling adaptive"}"

BUILD_TYPE="${BUILD_TYPE:-Release}"
CLEAN_BUILD="${CLEAN_BUILD:-0}"
SKIP_BUILD="${SKIP_BUILD:-0}"
RUN_TESTS="${RUN_TESTS:-0}"
JOBS_DETECTED="$(detect_jobs)"

SKIP_GEN="${SKIP_GEN:-0}"
FORCE_REGEN="${FORCE_REGEN:-0}"
LINK_DATA_IN_TEMP="${LINK_DATA_IN_TEMP:-1}"

WRITE_SAMPLES="${WRITE_SAMPLES:-0}"

# Safety defaults (EXP-4.md strongly recommends enum_cap + timeout).
EXTRA_RUN_ARGS="${EXTRA_RUN_ARGS:-"--enum_cap=50000000"}"
TIMEOUT_SEC="${TIMEOUT_SEC:-3600}"

# IMPORTANT: sjs_run gates rejection×sampling behind a flag.
# Default: enabled here to avoid "all exit=2" columns; set to 0 to reproduce sjs_sweep default behavior.
ALLOW_REJECTION_SAMPLING="${ALLOW_REJECTION_SAMPLING:-1}"

# Regimes (alpha settings)
REGIMES_IN="${REGIMES:-""}"
ALPHA_LIST_IN="${ALPHA_LIST:-""}"
ALPHA_IN="${ALPHA:-""}"

if [[ -n "$REGIMES_IN" ]]; then
  REGIMES_SPEC="$REGIMES_IN"
elif [[ -n "$ALPHA_LIST_IN" ]]; then
  REGIMES_SPEC=""
  for a in $ALPHA_LIST_IN; do REGIMES_SPEC+=" alpha${a}:${a}"; done
  REGIMES_SPEC="${REGIMES_SPEC# }"
elif [[ -n "$ALPHA_IN" ]]; then
  REGIMES_SPEC="alpha${ALPHA_IN}:${ALPHA_IN}"
else
  REGIMES_SPEC="$DEFAULT_REGIMES_SPEC"
fi

# T list: if T_LIST is non-empty use it; otherwise single T.
if [[ -n "${T_LIST}" ]]; then
  T_LIST_SPEC="$T_LIST"
else
  T_LIST_SPEC="$T"
fi

# --------------------------------------
# Parse regimes
# --------------------------------------
declare -a REG_NAMES=(); declare -a REG_ALPHAS=(); declare -a REG_DIRS=()
for token in $REGIMES_SPEC; do
  if [[ "$token" == *:* ]]; then name="${token%%:*}"; alpha="${token#*:}"; else alpha="$token"; name="alpha${alpha}"; fi
  [[ -n "$name" && -n "$alpha" ]] || die "Bad regime token: '$token' (expected name:alpha)"
  reg_dir="$(sanitize_token "$name")"
  REG_NAMES+=("$name"); REG_ALPHAS+=("$alpha"); REG_DIRS+=("$reg_dir")
done
[[ "${#REG_NAMES[@]}" -ge 1 ]] || die "No regimes configured."

# --------------------------------------
# Thread fairness (FORCE caps to THREADS)
# EXP-4.md requires single-thread fairness by env + --threads.
# --------------------------------------
export OMP_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"

# --------------------------------------
# Preflight
# --------------------------------------
need_cmd cmake
if ! command -v g++ >/dev/null 2>&1 && ! command -v clang++ >/dev/null 2>&1; then
  die "Need a C++ compiler in PATH (g++ or clang++)"
fi
need_cmd awk; need_cmd sed; need_cmd find

# GNU time for RSS (EXP-4.md: use ONLY for RSS).
TIME_BIN=""
if command -v gtime >/dev/null 2>&1; then
  TIME_BIN="$(command -v gtime)"
elif [[ -x "/usr/bin/time" ]]; then
  TIME_BIN="/usr/bin/time"
else
  if command -v time >/dev/null 2>&1; then
    TIME_BIN="$(command -v time)"
  fi
fi
[[ -n "$TIME_BIN" ]] || die "Cannot find an external 'time' command (GNU time recommended)."

# Require -v support; otherwise RSS collection is not reliable.
set +e
"$TIME_BIN" -v true >/dev/null 2>&1
rc_v=$?
set -e
[[ "$rc_v" -eq 0 ]] || die "'$TIME_BIN' does not support '-v'. Please install GNU time (e.g., 'gtime' on macOS)."

TIME_SUPPORTS_O=0
tmp_time_out="$(mktemp)"
set +e
"$TIME_BIN" -v -o "$tmp_time_out" true >/dev/null 2>&1
rc_time=$?
set -e
rm -f "$tmp_time_out"
[[ "$rc_time" -eq 0 ]] && TIME_SUPPORTS_O=1

# --------------------------------------
# Clean temp (overwrite) + ensure dirs
# --------------------------------------
msg "Repo root:      $ROOT_DIR"
msg "Build root:     $BUILD_ROOT"
msg "Temp root:      $TEMP_ROOT"
msg "Results root:   $RESULT_ROOT (overwritten)"
msg "Data root:      $DATA_ROOT"
msg "RUN_PROFILE:    $RUN_PROFILE"
msg "Build type:     $BUILD_TYPE"
msg "Jobs:           $JOBS_DETECTED"
msg "Threads:        $THREADS"
msg "N_LIST:         $N_LIST"
msg "T_LIST:         $T_LIST_SPEC"
msg "REPEATS:        $REPEATS"
msg "SEED:           $SEED"
msg "GEN_SEED:       $GEN_SEED"
msg "REGIMES:        $REGIMES_SPEC"
msg "METHODS:        $METHODS"
msg "VARIANTS:       $VARIANTS"
msg "WRITE_SAMPLES:  $WRITE_SAMPLES"
msg "EXTRA_RUN_ARGS: $EXTRA_RUN_ARGS"
msg "TIMEOUT_SEC:    $TIMEOUT_SEC (timeout cmd: $(command -v timeout >/dev/null 2>&1 && echo yes || echo no))"
msg "GNU time:       $TIME_BIN (supports -o: $TIME_SUPPORTS_O)"
msg "ALLOW_REJECTION_SAMPLING: $ALLOW_REJECTION_SAMPLING"

rm -rf "$TEMP_ROOT"
mkdir -p "$TEMP_ROOT" "$DATA_ROOT" "$TEMP_ROOT/meta"
GEN_REPORT_DIR="$TEMP_ROOT/meta/gen_reports"; mkdir -p "$GEN_REPORT_DIR"
DATA_LINK_ROOT="$TEMP_ROOT/datasets"; mkdir -p "$DATA_LINK_ROOT"

# --------------------------------------
# Manifest + sysinfo
# --------------------------------------
RUN_ID="$(date +%Y%m%d_%H%M%S)"
{
  echo "EXP=exp4"
  echo "RUN_ID=$RUN_ID"
  echo "DATE=$(date -Is)"
  echo "RUN_PROFILE=$RUN_PROFILE"
  echo "ROOT_DIR=$ROOT_DIR"
  echo "BUILD_ROOT=$BUILD_ROOT"
  echo "DATA_ROOT=$DATA_ROOT"
  echo "TEMP_ROOT=$TEMP_ROOT"
  echo "RESULT_ROOT=$RESULT_ROOT"
  echo "BUILD_TYPE=$BUILD_TYPE"
  echo "JOBS=$JOBS_DETECTED"
  echo "THREADS=$THREADS"
  echo "N_LIST=$N_LIST"
  echo "T_LIST=$T_LIST_SPEC"
  echo "REPEATS=$REPEATS"
  echo "SEED=$SEED"
  echo "GEN_SEED=$GEN_SEED"
  echo "METHODS=$METHODS"
  echo "VARIANTS=$VARIANTS"
  echo "WRITE_SAMPLES=$WRITE_SAMPLES"
  echo "EXTRA_RUN_ARGS=$EXTRA_RUN_ARGS"
  echo "TIMEOUT_SEC=$TIMEOUT_SEC"
  echo "REGIMES_SPEC=$REGIMES_SPEC"
  echo "ALLOW_REJECTION_SAMPLING=$ALLOW_REJECTION_SAMPLING"
  echo "REGIMES_PARSED_COUNT=${#REG_NAMES[@]}"
  for i in "${!REG_NAMES[@]}"; do
    echo "REGIME_$i=${REG_NAMES[$i]} (dir=${REG_DIRS[$i]}) alpha=${REG_ALPHAS[$i]}"
  done
} > "$TEMP_ROOT/meta/manifest.txt"

{
  echo "date: $(date -Is)"
  echo "pwd:  $ROOT_DIR"
  uname -a || true
  if command -v lscpu >/dev/null 2>&1; then echo; lscpu; fi
  if command -v free  >/dev/null 2>&1; then echo; free -h; fi
  echo
  cmake --version || true
  if command -v g++ >/dev/null 2>&1; then echo; g++ --version; fi
  if command -v clang++ >/dev/null 2>&1; then echo; clang++ --version; fi
  if command -v git >/dev/null 2>&1 && [[ -d "$ROOT_DIR/.git" ]]; then echo; git rev-parse HEAD; fi
} > "$TEMP_ROOT/meta/sysinfo.txt"

# --------------------------------------
# Build once
# --------------------------------------
BUILD_DIR="${BUILD_ROOT}/$(build_subdir_for_type "$BUILD_TYPE")"
if [[ "$SKIP_BUILD" != "1" ]]; then
  [[ "$CLEAN_BUILD" == "1" ]] && { msg "[build] CLEAN_BUILD=1 -> rm -rf $BUILD_DIR"; rm -rf "$BUILD_DIR"; }
  mkdir -p "$BUILD_DIR"
  msg "[build] configure: $BUILD_DIR ($BUILD_TYPE)"
  cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE" 2>&1 | tee "$TEMP_ROOT/meta/cmake_configure.log"
  msg "[build] build: -j$JOBS_DETECTED"
  cmake --build "$BUILD_DIR" -j"$JOBS_DETECTED" 2>&1 | tee "$TEMP_ROOT/meta/cmake_build.log"
  if [[ "$RUN_TESTS" == "1" ]]; then
    msg "[build] RUN_TESTS=1 -> ctest"
    ctest --test-dir "$BUILD_DIR" --output-on-failure 2>&1 | tee "$TEMP_ROOT/meta/ctest.log"
  fi
else
  msg "[build] SKIP_BUILD=1 (skipping compilation)"
fi

# --------------------------------------
# Locate executables
# --------------------------------------
SJS_GEN="$(resolve_exe sjs_gen_dataset "$BUILD_DIR" || true)"
SJS_RUN="$(resolve_exe sjs_run "$BUILD_DIR" || true)"
[[ -x "$SJS_GEN" ]] || die "Cannot find executable: sjs_gen_dataset (looked under $BUILD_DIR)"
[[ -x "$SJS_RUN" ]] || die "Cannot find executable: sjs_run (looked under $BUILD_DIR)"
msg "[bin] sjs_gen_dataset = $SJS_GEN"
msg "[bin] sjs_run         = $SJS_RUN"

# global commands log
COMMANDS_LOG_GLOBAL="$TEMP_ROOT/commands.log"; : > "$COMMANDS_LOG_GLOBAL"

# ============================================================
# Run each regime (alpha)
# ============================================================
msg "[run] EXP-4 start. Artifacts -> $TEMP_ROOT"

for i in "${!REG_NAMES[@]}"; do
  REG_NAME="${REG_NAMES[$i]}"; REG_DIR="${REG_DIRS[$i]}"; ALPHA_CUR="${REG_ALPHAS[$i]}"
  msg "------------------------------------------------------------"
  msg "[regime] $REG_NAME (dir=$REG_DIR, alpha=$ALPHA_CUR)"
  msg "------------------------------------------------------------"

  # Offline datasets live in data/synthetic/exp4/<reg_dir>
  REG_DATA_DIR="$DATA_ROOT/$REG_DIR"
  mkdir -p "$REG_DATA_DIR"

  # 1) Generate datasets once per (regime, N)
  if [[ "$SKIP_GEN" != "1" ]]; then
    msg "[gen] generating datasets into: $REG_DATA_DIR"
    for N in $N_LIST; do
      DS="exp4_${REG_DIR}_N${N}"
      R_BIN="$REG_DATA_DIR/${DS}_R.bin"
      S_BIN="$REG_DATA_DIR/${DS}_S.bin"

      if [[ "$FORCE_REGEN" != "1" && -f "$R_BIN" && -f "$S_BIN" ]]; then
        msg "[gen] skip $DS (exists)"
        continue
      fi

      msg "[gen] $DS: n_r=n_s=$N, alpha=$ALPHA_CUR, gen_seed=$GEN_SEED"
      GEN_STDOUT="$TEMP_ROOT/meta/gen_${DS}.stdout.log"
      GEN_STDERR="$TEMP_ROOT/meta/gen_${DS}.stderr.log"

      set +e
      run_with_optional_timeout "$TIMEOUT_SEC" \
        "$SJS_GEN" \
        --dataset_source=synthetic --gen=stripe --dataset="$DS" \
        --n_r="$N" --n_s="$N" --alpha="$ALPHA_CUR" --gen_seed="$GEN_SEED" \
        --out_dir="$REG_DATA_DIR" --write_csv=0 \
        --threads="$THREADS" \
        1> "$GEN_STDOUT" 2> "$GEN_STDERR"
      rc_gen=$?
      set -e
      [[ "$rc_gen" -eq 0 ]] || die "Dataset generation failed for $DS (exit=$rc_gen). See: $GEN_STDERR"
      [[ -f "$R_BIN" && -f "$S_BIN" ]] || die "Missing generated binary files for $DS under $REG_DATA_DIR"

      if [[ "$LINK_DATA_IN_TEMP" == "1" ]]; then
        mkdir -p "$DATA_LINK_ROOT/$REG_DIR"
        ln -sf "$R_BIN" "$DATA_LINK_ROOT/$REG_DIR/$(basename "$R_BIN")" 2>/dev/null || true
        ln -sf "$S_BIN" "$DATA_LINK_ROOT/$REG_DIR/$(basename "$S_BIN")" 2>/dev/null || true
      fi

      for j in "$REG_DATA_DIR/${DS}"*.json; do
        [[ -f "$j" ]] && cp -f "$j" "$GEN_REPORT_DIR/" 2>/dev/null || true
      done
    done
  else
    msg "[gen] SKIP_GEN=1 (expect binaries under $REG_DATA_DIR)"
  fi

  # 2) For each t, run methods × variants
  for T_CUR in $T_LIST_SPEC; do
    msg "[run] regime=$REG_NAME alpha=$ALPHA_CUR t=$T_CUR outputs (temp) under $TEMP_ROOT/$REG_DIR/alpha${ALPHA_CUR}_t${T_CUR}"

    REG_OUT_BASE="$TEMP_ROOT/$REG_DIR/alpha${ALPHA_CUR}_t${T_CUR}"
    MEM_DIR="$REG_OUT_BASE/mem"
    mkdir -p "$REG_OUT_BASE" "$MEM_DIR"

    RSS_CSV="$REG_OUT_BASE/exp4_rss_peak_kb.csv"
    STATUS_CSV="$REG_OUT_BASE/exp4_status.csv"
    COMMANDS_LOG="$REG_OUT_BASE/commands.log"
    echo "run_id,regime,alpha,t,repeats,seed,gen_seed,N,method,variant,exit_code,rss_kb,enum_truncated_any,run_csv" > "$RSS_CSV"
    echo "run_id,regime,alpha,t,repeats,seed,gen_seed,N,method,variant,exit_code,out_dir,stderr_log,enum_truncated_any" > "$STATUS_CSV"
    : > "$COMMANDS_LOG"

    {
      echo "RUN_ID=$RUN_ID"
      echo "REG_NAME=$REG_NAME"
      echo "REG_DIR=$REG_DIR"
      echo "ALPHA=$ALPHA_CUR"
      echo "T=$T_CUR"
      echo "N_LIST=$N_LIST"
      echo "METHODS=$METHODS"
      echo "VARIANTS=$VARIANTS"
      echo "THREADS=$THREADS"
      echo "WRITE_SAMPLES=$WRITE_SAMPLES"
      echo "EXTRA_RUN_ARGS=$EXTRA_RUN_ARGS"
      echo "TIMEOUT_SEC=$TIMEOUT_SEC"
      echo "ALLOW_REJECTION_SAMPLING=$ALLOW_REJECTION_SAMPLING"
      echo "DATE=$(date -Is)"
    } > "$REG_OUT_BASE/regime_manifest.txt"

    for N in $N_LIST; do
      DS="exp4_${REG_DIR}_N${N}"
      R_BIN="$REG_DATA_DIR/${DS}_R.bin"
      S_BIN="$REG_DATA_DIR/${DS}_S.bin"
      [[ -f "$R_BIN" && -f "$S_BIN" ]] || die "Missing dataset files for $DS: $R_BIN / $S_BIN"

      for method in $METHODS; do
        for variant in $VARIANTS; do
          OUT_DIR="$REG_OUT_BASE/N${N}/${method}/${variant}"
          mkdir -p "$OUT_DIR"

          STDOUT_LOG="$OUT_DIR/stdout.log"
          STDERR_LOG="$OUT_DIR/stderr.log"
          MEM_LOG="$MEM_DIR/N${N}_${method}_${variant}.timev.log"

          CMD=(
            "$SJS_RUN"
            --dataset_source=binary --dataset="$DS"
            --path_r="$R_BIN" --path_s="$S_BIN"
            --method="$method" --variant="$variant"
            --t="$T_CUR" --seed="$SEED" --repeats="$REPEATS"
            --out_dir="$OUT_DIR"
            --write_samples="$WRITE_SAMPLES"
            --threads="$THREADS"
          )

          # Method-specific gate: rejection×sampling requires an explicit allow flag.
          if [[ "$method" == "rejection" && "$variant" == "sampling" ]]; then
            if [[ "$ALLOW_REJECTION_SAMPLING" == "1" ]]; then
              CMD+=(--allow_rejection_sampling=1)
            else
              echo "[SKIP] rejection×sampling is gated. Set ALLOW_REJECTION_SAMPLING=1 to enable." > "$STDERR_LOG"
              echo "$RUN_ID,$REG_DIR,$ALPHA_CUR,$T_CUR,$REPEATS,$SEED,$GEN_SEED,$N,$method,$variant,2,,,$OUT_DIR/run.csv" >> "$RSS_CSV"
              echo "$RUN_ID,$REG_DIR,$ALPHA_CUR,$T_CUR,$REPEATS,$SEED,$GEN_SEED,$N,$method,$variant,2,$OUT_DIR,$STDERR_LOG," >> "$STATUS_CSV"
              warn "skipped rejection×sampling (ALLOW_REJECTION_SAMPLING!=1)"
              continue
            fi
          fi

          if [[ -n "$EXTRA_RUN_ARGS" ]]; then
            # shellcheck disable=SC2206
            extra=($EXTRA_RUN_ARGS)
            CMD+=("${extra[@]}")
          fi

          msg "[run] alpha=$ALPHA_CUR t=$T_CUR N=$N method=$method variant=$variant"
          printf "%s\t%s\n" "$(date -Is)" "${CMD[*]}" >> "$COMMANDS_LOG"
          printf "%s\t%s\n" "$(date -Is)" "${CMD[*]}" >> "$COMMANDS_LOG_GLOBAL"

          set +e
          if [[ "$TIME_SUPPORTS_O" == "1" ]]; then
            run_with_optional_timeout "$TIMEOUT_SEC" \
              "$TIME_BIN" -v -o "$MEM_LOG" "${CMD[@]}" 1> "$STDOUT_LOG" 2> "$STDERR_LOG"
            rc=$?
          else
            run_with_optional_timeout "$TIMEOUT_SEC" \
              "$TIME_BIN" -v "${CMD[@]}" 1> "$STDOUT_LOG" 2> "$STDERR_LOG"
            rc=$?
            cp "$STDERR_LOG" "$MEM_LOG" 2>/dev/null || true
          fi
          set -e

          RUN_CSV="$OUT_DIR/run.csv"
          rss_kb="$(parse_rss_kb_from_timev "$MEM_LOG")"
          enum_trunc_any="$(parse_enum_truncated_any_from_run_csv "$RUN_CSV")"

          echo "$RUN_ID,$REG_DIR,$ALPHA_CUR,$T_CUR,$REPEATS,$SEED,$GEN_SEED,$N,$method,$variant,$rc,$rss_kb,$enum_trunc_any,$RUN_CSV" >> "$RSS_CSV"
          echo "$RUN_ID,$REG_DIR,$ALPHA_CUR,$T_CUR,$REPEATS,$SEED,$GEN_SEED,$N,$method,$variant,$rc,$OUT_DIR,$STDERR_LOG,$enum_trunc_any" >> "$STATUS_CSV"

          if [[ "$rc" -ne 0 ]]; then
            warn "failed (exit=$rc). See: $STDERR_LOG"
            continue
          fi
          [[ -f "$RUN_CSV" ]] || warn "run succeeded but run.csv missing: $RUN_CSV"
        done
      done
    done

    msg "[done] regime=$REG_NAME alpha=$ALPHA_CUR t=$T_CUR finished."
    msg "  Raw results:            $REG_OUT_BASE/N*/<method>/<variant>/run.csv"
    msg "  Peak RSS table (CSV):   $RSS_CSV"
    msg "  Run status table (CSV): $STATUS_CSV"
  done
done

# --------------------------------------
# Sync temp -> results/raw/exp4 (overwrite) + optional history
# --------------------------------------
if [[ "$KEEP_HISTORY" == "1" ]]; then
  mkdir -p "$HISTORY_ROOT"
  HIST_DIR="$HISTORY_ROOT/$RUN_ID"
  rm -rf "$HIST_DIR"
  mkdir -p "$HIST_DIR"
  cp -a "$TEMP_ROOT/." "$HIST_DIR/"
  msg "[history] Saved a copy to: $HIST_DIR"
fi

msg "[sync] Copy temp outputs -> results/raw/exp4 (overwrite)"
rm -rf "$RESULT_ROOT"
mkdir -p "$RESULT_ROOT"
cp -a "$TEMP_ROOT/." "$RESULT_ROOT/"

msg "[done] EXP-4 finished ✅"
msg "  Temp:    $TEMP_ROOT"
msg "  Results: $RESULT_ROOT"
msg "  Data:    $DATA_ROOT"
if [[ "$KEEP_HISTORY" == "1" ]]; then
  msg "  History: $HISTORY_ROOT/$RUN_ID"
fi
