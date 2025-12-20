#!/usr/bin/env bash
# ============================================================
# EXP-4: Scalability vs N (RQ4)
#
# Scan N (n_r = n_s = N): 50k -> 800k, while fixing t and alpha.
#
# This script runs EXP-4 end-to-end from scratch:
#   1) Configure + build (Release) via CMake
#   2) Generate offline binary datasets for each N (stripe/alpha-controlled)
#   3) Run each method × variant, with repeats, recording:
#        - Algorithm phase timings in <out_dir>/run.csv  (use for runtime/build-time curves)
#        - Peak RSS via GNU time -v                      (use for memory table)
#
# Outputs:
#   data/synthetic/exp4/*.bin
#   results/raw/exp4/<tag>/
#     ├─ alpha<alpha>_t<t>/
#     │   ├─ N<N>/<method>/<variant>/run.csv (+ stdout/stderr logs)
#     │   ├─ mem/*.timev.log
#     │   ├─ exp4_rss_peak_kb.csv
#     │   ├─ exp4_status.csv
#     │   └─ commands.log
#     └─ meta/{manifest.txt,sysinfo.txt,cmake_*.log,gen_*.log}
#
# Usage:
#   bash run/run_exp4.sh
#
# Override defaults (examples):
#   ALPHA=5 T=100000 REPEATS=3 SEED=1 GEN_SEED=1 bash run/run_exp4.sh
#   N_LIST="50000 100000 200000" METHODS="ours aabb" VARIANTS="sampling adaptive" bash run/run_exp4.sh
#   SKIP_BUILD=1 SKIP_GEN=1 bash run/run_exp4.sh
#
# Notes:
#   - EXP-4 timing in the paper should use the program's internal phase times
#     from run.csv (Build/Count/Enumerate/Sample). The external /usr/bin/time
#     numbers include IO and are only used here to extract peak RSS.
#   - Some method×variant combos may be unsupported or may fail at large |J|.
#     This script records failures in exp4_status.csv and continues.
# ============================================================

set -euo pipefail

# --- helpers ---
msg() { echo "[$(date -Is)] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

# Find an executable produced by CMake, regardless of output directory layout.
resolve_exe() {
  local name="$1"
  local root="$2"

  local candidates=(
    "$root/$name"
    "$root/apps/$name"
    "$root/bin/$name"
    "$root/src/apps/$name"
    "$root/../$name"
  )

  for p in "${candidates[@]}"; do
    if [[ -x "$p" ]]; then
      echo "$p"
      return 0
    fi
  done

  # fallback: search a bit
  local found
  found="$(find "$root" -maxdepth 4 -type f -name "$name" -perm -111 2>/dev/null | head -n 1 || true)"
  if [[ -n "$found" && -x "$found" ]]; then
    echo "$found"
    return 0
  fi

  return 1
}

# --- locate repo root ---
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- defaults (override via env vars) ---
N_LIST="${N_LIST:-"50000 100000 200000 400000 800000"}"
ALPHA="${ALPHA:-5}"
T="${T:-100000}"
REPEATS="${REPEATS:-3}"
SEED="${SEED:-1}"
GEN_SEED="${GEN_SEED:-1}"

# All methods listed in README (Dim=2). You can override METHODS="..."
METHODS="${METHODS:-"ours aabb interval_tree kd_tree r_tree range_tree pbsm tlsop sirs rejection tsunami"}"
VARIANTS="${VARIANTS:-"sampling enum_sampling adaptive"}"

# Build / output dirs
BUILD_DIR="${BUILD_DIR:-"$ROOT_DIR/build"}"
DATA_DIR="${DATA_DIR:-"$ROOT_DIR/data/synthetic/exp4"}"

EXP4_TAG="${EXP4_TAG:-"exp4_$(date +%Y%m%d_%H%M%S)"}"
OUT_BASE="${OUT_BASE:-"$ROOT_DIR/results/raw/exp4/$EXP4_TAG"}"

# Build parallelism
JOBS="${JOBS:-""}"
if [[ -z "$JOBS" ]]; then
  if command -v nproc >/dev/null 2>&1; then JOBS="$(nproc)"; else JOBS="4"; fi
fi

# Optional controls
SKIP_BUILD="${SKIP_BUILD:-0}"
SKIP_GEN="${SKIP_GEN:-0}"
FORCE_REGEN="${FORCE_REGEN:-0}"
WRITE_SAMPLES="${WRITE_SAMPLES:-0}"   # 1 => write per-repeat sample TSVs (can be huge)
EXTRA_RUN_ARGS="${EXTRA_RUN_ARGS:-""}" # extra flags forwarded to sjs_run, e.g. "--enum_cap=50000000"

# Threading guardrails (keep experiments single-threaded/fair)
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

# --- sanity checks ---
need_cmd cmake
# Need a compiler (g++ or clang++)
if ! command -v g++ >/dev/null 2>&1 && ! command -v clang++ >/dev/null 2>&1; then
  die "Need a C++ compiler in PATH (g++ or clang++)"
fi

mkdir -p "$BUILD_DIR" "$DATA_DIR" "$OUT_BASE"

META_DIR="$OUT_BASE/meta"
LOG_DIR="$OUT_BASE/meta"
MEM_DIR="$OUT_BASE/alpha${ALPHA}_t${T}/mem"
mkdir -p "$META_DIR" "$MEM_DIR"

# --- record manifest + sysinfo ---
{
  echo "EXP4_TAG=$EXP4_TAG"
  echo "ROOT_DIR=$ROOT_DIR"
  echo "BUILD_DIR=$BUILD_DIR"
  echo "DATA_DIR=$DATA_DIR"
  echo "OUT_BASE=$OUT_BASE"
  echo "N_LIST=$N_LIST"
  echo "ALPHA=$ALPHA"
  echo "T=$T"
  echo "REPEATS=$REPEATS"
  echo "SEED=$SEED"
  echo "GEN_SEED=$GEN_SEED"
  echo "METHODS=$METHODS"
  echo "VARIANTS=$VARIANTS"
  echo "WRITE_SAMPLES=$WRITE_SAMPLES"
  echo "EXTRA_RUN_ARGS=$EXTRA_RUN_ARGS"
  echo "JOBS=$JOBS"
  echo "DATE=$(date -Is)"
} > "$META_DIR/manifest.txt"

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
} > "$META_DIR/sysinfo.txt"

# --- build ---
if [[ "$SKIP_BUILD" != "1" ]]; then
  msg "[build] configure (Release) ..."
  pushd "$BUILD_DIR" >/dev/null
  cmake .. -DCMAKE_BUILD_TYPE=Release 2>&1 | tee "$META_DIR/cmake_configure.log"
  msg "[build] build (-j$JOBS) ..."
  cmake --build . -j"$JOBS" 2>&1 | tee "$META_DIR/cmake_build.log"
  popd >/dev/null
else
  msg "[build] SKIP_BUILD=1 (skipping compilation)"
fi

# --- locate executables ---
SJS_GEN="$(resolve_exe sjs_gen_dataset "$BUILD_DIR" || true)"
SJS_RUN="$(resolve_exe sjs_run        "$BUILD_DIR" || true)"

[[ -x "$SJS_GEN" ]] || die "Cannot find executable: sjs_gen_dataset (looked under $BUILD_DIR)"
[[ -x "$SJS_RUN" ]] || die "Cannot find executable: sjs_run (looked under $BUILD_DIR)"

msg "[bin] sjs_gen_dataset = $SJS_GEN"
msg "[bin] sjs_run         = $SJS_RUN"

# --- GNU time for RSS ---
TIME_BIN=""
if [[ -x "/usr/bin/time" ]]; then
  TIME_BIN="/usr/bin/time"
elif command -v gtime >/dev/null 2>&1; then
  TIME_BIN="$(command -v gtime)"
elif command -v time >/dev/null 2>&1; then
  # best effort (may be shell builtin on some systems)
  TIME_BIN="$(command -v time)"
fi
[[ -n "$TIME_BIN" ]] || die "Cannot find an external 'time' command to measure RSS (need GNU time recommended)"

TIME_SUPPORTS_O=0
tmp_time_out="$(mktemp)"
set +e
"$TIME_BIN" -v -o "$tmp_time_out" true >/dev/null 2>&1
rc_time=$?
set -e
rm -f "$tmp_time_out"
if [[ "$rc_time" -eq 0 ]]; then TIME_SUPPORTS_O=1; fi
msg "[time] using: $TIME_BIN (supports -o: $TIME_SUPPORTS_O)"

# --- generate datasets (offline, binary) ---
if [[ "$SKIP_GEN" != "1" ]]; then
  msg "[gen] generating datasets into: $DATA_DIR"
  for N in $N_LIST; do
    DS="exp4_N${N}"
    R_BIN="$DATA_DIR/${DS}_R.bin"
    S_BIN="$DATA_DIR/${DS}_S.bin"

    if [[ "$FORCE_REGEN" != "1" && -f "$R_BIN" && -f "$S_BIN" ]]; then
      msg "[gen] skip $DS (exists)"
      continue
    fi

    msg "[gen] $DS: n_r=n_s=$N, alpha=$ALPHA, gen_seed=$GEN_SEED"
    set +e
    "$SJS_GEN" \
      --dataset_source=synthetic --gen=stripe --dataset="$DS" \
      --n_r="$N" --n_s="$N" --alpha="$ALPHA" --gen_seed="$GEN_SEED" \
      --out_dir="$DATA_DIR" --write_csv=0 \
      1> "$META_DIR/gen_${DS}.stdout.log" 2> "$META_DIR/gen_${DS}.stderr.log"
    rc_gen=$?
    set -e

    if [[ "$rc_gen" -ne 0 ]]; then
      die "Dataset generation failed for $DS (exit=$rc_gen). See $META_DIR/gen_${DS}.stderr.log"
    fi
    [[ -f "$R_BIN" && -f "$S_BIN" ]] || die "Missing generated binary files for $DS under $DATA_DIR"
  done
else
  msg "[gen] SKIP_GEN=1 (skipping dataset generation)"
fi

# --- run experiments ---
EXP_DIR="$OUT_BASE/alpha${ALPHA}_t${T}"
mkdir -p "$EXP_DIR"

RSS_CSV="$EXP_DIR/exp4_rss_peak_kb.csv"
STATUS_CSV="$EXP_DIR/exp4_status.csv"
COMMANDS_LOG="$EXP_DIR/commands.log"

echo "tag,alpha,t,repeats,seed,gen_seed,N,method,variant,exit_code,rss_kb,run_csv" > "$RSS_CSV"
echo "tag,alpha,t,repeats,seed,gen_seed,N,method,variant,exit_code,out_dir,stderr_log" > "$STATUS_CSV"
: > "$COMMANDS_LOG"

msg "[run] EXP-4 start: outputs under $EXP_DIR"

for N in $N_LIST; do
  DS="exp4_N${N}"
  R_BIN="$DATA_DIR/${DS}_R.bin"
  S_BIN="$DATA_DIR/${DS}_S.bin"
  [[ -f "$R_BIN" && -f "$S_BIN" ]] || die "Missing dataset files for $DS: $R_BIN / $S_BIN"

  for method in $METHODS; do
    for variant in $VARIANTS; do
      OUT_DIR="$EXP_DIR/N${N}/${method}/${variant}"
      mkdir -p "$OUT_DIR"

      STDOUT_LOG="$OUT_DIR/stdout.log"
      STDERR_LOG="$OUT_DIR/stderr.log"
      MEM_LOG="$MEM_DIR/N${N}_${method}_${variant}.timev.log"

      CMD=(
        "$SJS_RUN"
        --dataset_source=binary --dataset="$DS"
        --path_r="$R_BIN" --path_s="$S_BIN"
        --method="$method" --variant="$variant"
        --t="$T" --seed="$SEED" --repeats="$REPEATS"
        --out_dir="$OUT_DIR"
        --write_samples="$WRITE_SAMPLES"
      )

      if [[ -n "$EXTRA_RUN_ARGS" ]]; then
        # shellcheck disable=SC2206
        extra=($EXTRA_RUN_ARGS)
        CMD+=("${extra[@]}")
      fi

      msg "[run] N=$N method=$method variant=$variant"
      printf "%s\t%s\n" "$(date -Is)" "${CMD[*]}" >> "$COMMANDS_LOG"

      # Execute, but do NOT abort EXP-4 on a single failure.
      set +e
      if [[ "$TIME_SUPPORTS_O" == "1" ]]; then
        "$TIME_BIN" -v -o "$MEM_LOG" "${CMD[@]}" 1> "$STDOUT_LOG" 2> "$STDERR_LOG"
        rc=$?
      else
        # fallback: time output will be mixed into stderr; still usable for RSS in many setups
        "$TIME_BIN" -v "${CMD[@]}" 1> "$STDOUT_LOG" 2> "$STDERR_LOG"
        rc=$?
        cp "$STDERR_LOG" "$MEM_LOG" 2>/dev/null || true
      fi
      set -e

      RUN_CSV="$OUT_DIR/run.csv"

      rss_kb=""
      if [[ -f "$MEM_LOG" ]]; then
        rss_kb="$(awk -F: '/Maximum resident set size/ {gsub(/^[ \t]+/,"",$2); print $2}' "$MEM_LOG" | tail -n 1 | tr -d '\r')"
      fi

      echo "$EXP4_TAG,$ALPHA,$T,$REPEATS,$SEED,$GEN_SEED,$N,$method,$variant,$rc,$rss_kb,$RUN_CSV" >> "$RSS_CSV"
      echo "$EXP4_TAG,$ALPHA,$T,$REPEATS,$SEED,$GEN_SEED,$N,$method,$variant,$rc,$OUT_DIR,$STDERR_LOG" >> "$STATUS_CSV"

      if [[ "$rc" -ne 0 ]]; then
        msg "  [WARN] failed (exit=$rc). See: $STDERR_LOG"
        continue
      fi
      if [[ ! -f "$RUN_CSV" ]]; then
        msg "  [WARN] run succeeded but run.csv missing: $RUN_CSV"
      fi
    done
  done
done

msg "[done] EXP-4 finished."
msg "  Raw results (per run):   $EXP_DIR/N*/<method>/<variant>/run.csv"
msg "  Peak RSS table (CSV):    $RSS_CSV"
msg "  Run status table (CSV):  $STATUS_CSV"
msg "  Commands log:            $COMMANDS_LOG"
msg "  Meta (sysinfo, build):   $META_DIR"
