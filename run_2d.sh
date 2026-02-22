#!/usr/bin/env bash
# run_2d.sh (repo root)
#
# Experiment II: Synthetic performance test (d=2, 7-model comparison)
#
# This script:
#   1) Builds C++ binaries (Release by default)
#   2) Generates/reuses all required alacarte-rectgen synthetic datasets
#   3) Runs the full Experiment II grid for d=2 (A_alpha, B_N, C_t)
#   4) Splits jobs into light/heavy queues and runs with bounded parallelism
#   5) Is resumable and will retry failed jobs across multiple attempts
#
# Key env knobs (safe defaults):
#   BUILD_TYPE=Release
#   PARALLEL_LIGHT=8
#   PARALLEL_HEAVY=1
#   HEAVY_T=30000000
#   TIME_LIMIT_SEC=7200
#   TIME_KILL_GRACE_SEC=30
#   REPEATS=5
#   MASTER_SEED=1
#   JOB_THREADS=1
#   MAX_ATTEMPTS=3
#   BUDGET=10000000
#   W_SMALL=1024
#   PREFETCH=1
#   CLEAN_RESULTS=0   # 1 => rm -rf results/exp2_d2 before running
#   CLEAN_DATA=0      # 1 => rm -rf data/synthetic/exp2_d2 before generating
#
# Notes:
#   - Framework name mapping (per SJS-HighDims.md):
#       Enum      -> --variant=enum_sampling
#       Adaptive  -> --variant=sampling
#       Sampling  -> --variant=adaptive   (budgeted cache/prefetch)
#   - Heavy jobs are: any enum_sampling OR any job with t >= HEAVY_T.

set -Eeuo pipefail
IFS=$'\n\t'

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

# ----------------------
# Logging helpers
# ----------------------
_ts() { date '+%F %T'; }
log() { echo "[$(_ts)] $*"; }
warn() { echo "[$(_ts)] [WARN] $*" >&2; }
err() { echo "[$(_ts)] [ERR ] $*" >&2; }

# ----------------------
# Config (env-overridable)
# ----------------------
BUILD_TYPE="${BUILD_TYPE:-Release}"
PARALLEL_LIGHT="${PARALLEL_LIGHT:-8}"
PARALLEL_HEAVY="${PARALLEL_HEAVY:-1}"
HEAVY_T="${HEAVY_T:-30000000}"
TIME_LIMIT_SEC="${TIME_LIMIT_SEC:-7200}"
TIME_KILL_GRACE_SEC="${TIME_KILL_GRACE_SEC:-30}"
REPEATS="${REPEATS:-5}"
MASTER_SEED="${MASTER_SEED:-1}"
JOB_THREADS="${JOB_THREADS:-1}"
MAX_ATTEMPTS="${MAX_ATTEMPTS:-3}"

BUDGET="${BUDGET:-10000000}"
W_SMALL="${W_SMALL:-1024}"
PREFETCH="${PREFETCH:-1}"

CLEAN_RESULTS="${CLEAN_RESULTS:-0}"
CLEAN_DATA="${CLEAN_DATA:-0}"

RESULTS_DIR="${RESULTS_DIR:-results/exp2_d2}"
DATA_DIR="${DATA_DIR:-data/synthetic/exp2_d2}"
DRIVER_DIR="${RESULTS_DIR}/_driver"

# ----------------------
# Internal helpers
# ----------------------
build_subdir_for_type() {
  case "$1" in
    Release|release) echo "release" ;;
    Debug|debug) echo "debug" ;;
    RelWithDebInfo|relwithdebinfo) echo "relwithdebinfo" ;;
    MinSizeRel|minsizerel) echo "minsizerel" ;;
    *) echo "$(echo "$1" | tr '[:upper:]' '[:lower:]')" ;;
  esac
}

alpha_tag() {
  # Turn e.g. "0.1" -> "0p1". Keep ints as-is.
  local a="$1"
  # Normalize possible scientific format (rare in our fixed grids).
  # If it contains 'e' or 'E', fall back to python formatting.
  if [[ "$a" == *e* || "$a" == *E* ]]; then
    python3 - <<PY
import math
x=float("$a")
# Prefer a short decimal representation when possible.
# We only need stable file names, not perfect reversibility.
if abs(x-round(x))<1e-12:
  s=str(int(round(x)))
else:
  s=("%.12g"%x)
print(s.replace('.','p'))
PY
    return
  fi
  echo "${a//./p}"
}

need_cmd() {
  local c="$1"
  command -v "$c" >/dev/null 2>&1 || {
    err "Missing required command: $c"
    exit 2
  }
}

ensure_dirs() {
  mkdir -p "$RESULTS_DIR" "$DATA_DIR" "$DRIVER_DIR"
}

maybe_clean() {
  if [[ "$CLEAN_RESULTS" == "1" ]]; then
    warn "CLEAN_RESULTS=1: removing $RESULTS_DIR"
    rm -rf "$RESULTS_DIR"
  fi
  if [[ "$CLEAN_DATA" == "1" ]]; then
    warn "CLEAN_DATA=1: removing $DATA_DIR"
    rm -rf "$DATA_DIR"
  fi
}

# ----------------------
# Dataset generation
# ----------------------

gen_dataset_one() {
  local N="$1"
  local alpha="$2"
  local seed="$3"

  local tag
  tag="$(alpha_tag "$alpha")"
  local name="rectgen_d2_N${N}_alpha${tag}_seed${seed}"

  local ddir="$DATA_DIR/$name"
  local rbin="$ddir/R.bin"
  local sbin="$ddir/S.bin"
  local rep="$ddir/gen_report.json"

  mkdir -p "$ddir"

  if [[ -s "$rbin" && -s "$sbin" && -s "$rep" ]]; then
    log "[DS] OK cached: $name"
    return 0
  fi

  log "[DS] Generating: $name (N=$N alpha_out=$alpha seed=$seed)"

  # If partially present, remove to avoid mixing old/new.
  rm -f "$rbin" "$sbin" "$rep"

  set +e
  python3 tools/alacarte_rectgen_generate.py \
    --nR "$N" --nS "$N" --d 2 --alpha_out "$alpha" --seed "$seed" \
    --out_r "$rbin" --out_s "$sbin" \
    --dataset_name "$name" \
    --report_path "$rep" \
    --audit_seed "$seed" \
    >/dev/null 2>"$ddir/gen.stderr.log"
  local ec=$?
  set -e

  if [[ $ec -ne 0 ]]; then
    err "[DS] FAILED: $name (exit=$ec). See: $ddir/gen.stderr.log"
    return $ec
  fi

  # Basic sanity: ensure files exist and non-empty.
  if [[ ! -s "$rbin" || ! -s "$sbin" || ! -s "$rep" ]]; then
    err "[DS] FAILED: $name produced missing/empty outputs"
    return 2
  fi

  log "[DS] DONE: $name"
  return 0
}

# ----------------------
# Job list (Experiment II, d=2)
# ----------------------

write_job_manifest() {
  local manifest="$1"

  # Parameter grids (d=2 locked).
  local -a ALPHAS=(0.1 0.3 1 3 10 30 100 300 1000)
  local -a NS=(100000 200000 500000 1000000 2000000 5000000)
  local -a TS=(100000 300000 1000000 3000000 10000000 30000000 100000000 300000000)

  # 7 models: label, method, variant (CLI).
  # IMPORTANT: naming follows paper/framework labels:
  #   - "*_adaptive" means paper's Adaptive framework (2-pass) => CLI variant=sampling.
  #   - "*_sampling" means paper's Sampling framework (budgeted) => CLI variant=adaptive.
  #
  # CRITICAL FIX:
  #   Use $'...' to embed REAL tab characters. DO NOT write "\t" inside "..."
  #   (that would be a literal backslash+t, and parsing would fail again).
  local -a MODELS=(
    $'sjs_enum\tours\tenum_sampling'
    $'sjs_adaptive\tours\tsampling'
    $'sjs_sampling\tours\tadaptive'
    $'rt_enum\trange_tree\tenum_sampling'
    $'rt_adaptive\trange_tree\tsampling'
    $'rt_sampling\trange_tree\tadaptive'
    $'kds_sampling\tkd_tree\tsampling'
  )

  rm -f "$manifest"
  mkdir -p "$(dirname "$manifest")"

  local lines=0

  # Helper to emit one job row (as real TSV columns).
  _emit() {
    local exp="$1" dataset="$2" t="$3" label="$4" method="$5" variant="$6"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$exp" "$dataset" "$t" "$label" "$method" "$variant" >>"$manifest"
    lines=$((lines+1))
  }

  # (A) vary alpha_out, fixed N=1e6, t=1e5
  local exp="A_alpha"
  local N_fixed=1000000
  local t_fixed=100000
  for alpha in "${ALPHAS[@]}"; do
    local tag name
    tag="$(alpha_tag "$alpha")"
    name="rectgen_d2_N${N_fixed}_alpha${tag}_seed${MASTER_SEED}"
    for mspec in "${MODELS[@]}"; do
      IFS=$'\t' read -r label method variant <<<"$mspec"
      _emit "$exp" "$name" "$t_fixed" "$label" "$method" "$variant"
    done
  done

  # (B) vary N, fixed alpha_out=200, t=1e5
  exp="B_N"
  local alpha_fixed=200
  for N in "${NS[@]}"; do
    local tag name
    tag="$(alpha_tag "$alpha_fixed")"
    name="rectgen_d2_N${N}_alpha${tag}_seed${MASTER_SEED}"
    for mspec in "${MODELS[@]}"; do
      IFS=$'\t' read -r label method variant <<<"$mspec"
      _emit "$exp" "$name" "$t_fixed" "$label" "$method" "$variant"
    done
  done

  # (C) vary t, fixed N=1e6, alpha_out=200
  exp="C_t"
  local N_ct=1000000
  local tag_ct name_ct
  tag_ct="$(alpha_tag "$alpha_fixed")"
  name_ct="rectgen_d2_N${N_ct}_alpha${tag_ct}_seed${MASTER_SEED}"
  for t in "${TS[@]}"; do
    for mspec in "${MODELS[@]}"; do
      IFS=$'\t' read -r label method variant <<<"$mspec"
      _emit "$exp" "$name_ct" "$t" "$label" "$method" "$variant"
    done
  done

  log "Wrote job manifest: $manifest (lines=$lines)"
}

# ----------------------
# Job execution
# ----------------------

out_dir_for_job() {
  local exp="$1" dataset="$2" t="$3" label="$4" method="$5" variant="$6"
  # Clean directory layout (no literal '\t' in names):
  #   results/exp2_d2/<exp>/<dataset>/t<t>/<label>/<method>/<variant>/
  echo "$RESULTS_DIR/$exp/$dataset/t${t}/$label/$method/$variant"
}

is_heavy_job() {
  local t="$1" variant="$2"
  if [[ "$variant" == "enum_sampling" ]]; then
    return 0
  fi
  if [[ "$t" -ge "$HEAVY_T" ]]; then
    return 0
  fi
  return 1
}

run_one_job() {
  local exp="$1" dataset="$2" t="$3" label="$4" method="$5" variant="$6" attempt="$7"

  local out_dir
  out_dir="$(out_dir_for_job "$exp" "$dataset" "$t" "$label" "$method" "$variant")"
  local done_marker="$out_dir/DONE"
  local status_file="$out_dir/status.txt"
  local log_file="$out_dir/run.log"

  # Resume: if DONE marker exists, skip.
  if [[ -f "$done_marker" ]]; then
    return 0
  fi

  mkdir -p "$out_dir"

  local rbin="$DATA_DIR/$dataset/R.bin"
  local sbin="$DATA_DIR/$dataset/S.bin"
  if [[ ! -s "$rbin" || ! -s "$sbin" ]]; then
    {
      echo "timestamp=$(_ts)"
      echo "attempt=$attempt"
      echo "exp=$exp"
      echo "dataset=$dataset"
      echo "t=$t"
      echo "label=$label"
      echo "method=$method"
      echo "variant=$variant"
      echo "exit_code=2"
      echo "reason=DATASET_MISSING"
      echo "out_dir=$out_dir"
      echo "path_r=$rbin"
      echo "path_s=$sbin"
    } >"$status_file"
    err "[JOB] Dataset missing for $exp/$dataset (t=$t) => $label/$method/$variant"
    return 0
  fi

  local -a cmd=(
    "$SJS_RUN"
    --dataset_source=binary
    "--dataset=$dataset"
    --dim=2
    "--path_r=$rbin"
    "--path_s=$sbin"
    "--method=$method"
    "--variant=$variant"
    "--t=$t"
    "--seed=$MASTER_SEED"
    "--repeats=$REPEATS"
    "--threads=$JOB_THREADS"
    "--out_dir=$out_dir"
    "--results_file=$out_dir/run.csv"
    --write_samples=0
    --verify=0
    --enum_cap=0
    "--budget=$BUDGET"
    "--w_small=$W_SMALL"
    "--prefetch=$PREFETCH"
  )

  set +e
  {
    /usr/bin/time -v \
      timeout -k "${TIME_KILL_GRACE_SEC}s" "${TIME_LIMIT_SEC}s" \
      "${cmd[@]}"
  } >"$log_file" 2>&1
  local ec=$?
  set -e

  local reason
  case "$ec" in
    0) reason="OK" ;;
    124) reason="TIMEOUT" ;;
    137) reason="KILLED" ;;
    139) reason="SEGFAULT" ;;
    *) reason="EXIT_${ec}" ;;
  esac

  {
    echo "timestamp=$(_ts)"
    echo "attempt=$attempt"
    echo "exp=$exp"
    echo "dataset=$dataset"
    echo "t=$t"
    echo "label=$label"
    echo "method=$method"
    echo "variant=$variant"
    echo "exit_code=$ec"
    echo "reason=$reason"
    echo "out_dir=$out_dir"
  } >"$status_file"

  if [[ $ec -eq 0 ]]; then
    : >"$done_marker"
  fi

  return 0
}

write_pending_lists() {
  local manifest="$1" pending_light="$2" pending_heavy="$3"

  rm -f "$pending_light" "$pending_heavy"
  : >"$pending_light"
  : >"$pending_heavy"

  local exp dataset t label method variant
  while IFS=$'\t' read -r exp dataset t label method variant; do
    [[ -z "${exp:-}" ]] && continue

    local out_dir
    out_dir="$(out_dir_for_job "$exp" "$dataset" "$t" "$label" "$method" "$variant")"
    if [[ -f "$out_dir/DONE" ]]; then
      continue
    fi

    if is_heavy_job "$t" "$variant"; then
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$exp" "$dataset" "$t" "$label" "$method" "$variant" >>"$pending_heavy"
    else
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$exp" "$dataset" "$t" "$label" "$method" "$variant" >>"$pending_light"
    fi
  done <"$manifest"
}

count_lines() {
  local f="$1"
  if [[ ! -f "$f" ]]; then
    echo 0
    return
  fi
  awk 'NF{c++} END{print c+0}' "$f"
}

run_jobs_from_tsv() {
  local tsv="$1" parallel="$2" attempt="$3"

  local total
  total="$(count_lines "$tsv")"
  if [[ "$total" -le 0 ]]; then
    return 0
  fi

  log "Running $(basename "$tsv") (jobs=$total, parallel=$parallel, attempt=$attempt)"

  # Throttled background execution (portable: no `wait -n` dependency).
  local -a pids=()
  local exp dataset t label method variant

  while IFS=$'\t' read -r exp dataset t label method variant; do
    [[ -z "${exp:-}" ]] && continue

    run_one_job "$exp" "$dataset" "$t" "$label" "$method" "$variant" "$attempt" &
    pids+=("$!")

    if [[ "${#pids[@]}" -ge "$parallel" ]]; then
      wait "${pids[0]}" || true
      pids=("${pids[@]:1}")
    fi
  done <"$tsv"

  local pid
  for pid in "${pids[@]}"; do
    wait "$pid" || true
  done

  return 0
}

write_final_summary() {
  local manifest="$1" out_path="$2"

  local total done failed
  total="$(count_lines "$manifest")"
  done=0
  failed=0

  {
    echo "Experiment II (d=2) - Final Summary"
    echo "timestamp=$(_ts)"
    echo "total_jobs=$total"
    echo "results_dir=$RESULTS_DIR"
    echo
  } >"$out_path"

  local exp dataset t label method variant
  while IFS=$'\t' read -r exp dataset t label method variant; do
    [[ -z "${exp:-}" ]] && continue
    local out_dir
    out_dir="$(out_dir_for_job "$exp" "$dataset" "$t" "$label" "$method" "$variant")"
    if [[ -f "$out_dir/DONE" ]]; then
      done=$((done+1))
    fi
  done <"$manifest"

  failed=$((total-done))

  {
    echo "done_jobs=$done"
    echo "failed_or_skipped_jobs=$failed"
    echo
  } >>"$out_path"

  if [[ "$failed" -le 0 ]]; then
    echo "All jobs DONE." >>"$out_path"
    return 0
  fi

  echo "FAILED / SKIPPED jobs (no DONE marker):" >>"$out_path"

  while IFS=$'\t' read -r exp dataset t label method variant; do
    [[ -z "${exp:-}" ]] && continue
    local out_dir status_file log_file
    out_dir="$(out_dir_for_job "$exp" "$dataset" "$t" "$label" "$method" "$variant")"
    if [[ -f "$out_dir/DONE" ]]; then
      continue
    fi
    status_file="$out_dir/status.txt"
    log_file="$out_dir/run.log"
    printf -- "- %s\t%s\tt=%s\t%s\t%s\t%s\tstatus=%s\tlog=%s\n" \
      "$exp" "$dataset" "$t" "$label" "$method" "$variant" "$status_file" "$log_file" >>"$out_path"
  done <"$manifest"

  return 0
}

# ----------------------
# Main
# ----------------------

main() {
  need_cmd python3
  need_cmd timeout
  need_cmd /usr/bin/time

  log "BUILD_TYPE=$BUILD_TYPE"
  log "PARALLEL_LIGHT=$PARALLEL_LIGHT PARALLEL_HEAVY=$PARALLEL_HEAVY HEAVY_T=$HEAVY_T"
  log "TIME_LIMIT_SEC=$TIME_LIMIT_SEC (kill grace ${TIME_KILL_GRACE_SEC}s)"
  log "REPEATS=$REPEATS (warmup 1) MASTER_SEED=$MASTER_SEED"

  maybe_clean
  ensure_dirs

  # 1) Build.
  log "Building C++ binaries (this may take a while the first time)..."
  BUILD_TYPE="$BUILD_TYPE" ./run.sh build

  local build_subdir
  build_subdir="$(build_subdir_for_type "$BUILD_TYPE")"
  local build_dir="$ROOT/build/$build_subdir"
  SJS_RUN="$build_dir/sjs_run"

  if [[ ! -x "$SJS_RUN" ]]; then
    err "Cannot find executable: $SJS_RUN"
    err "Did the build succeed?"
    exit 2
  fi

  log "Using sjs_run: ${SJS_RUN#$ROOT/}"

  # 2) Generate datasets.
  log "Generating (or reusing) all required synthetic datasets..."

  local -a ALPHAS=(0.1 0.3 1 3 10 30 100 300 1000)
  local -a NS=(100000 200000 500000 1000000 2000000 5000000)

  local N_fixed=1000000
  for alpha in "${ALPHAS[@]}"; do
    gen_dataset_one "$N_fixed" "$alpha" "$MASTER_SEED" || true
  done

  local alpha_fixed=200
  for N in "${NS[@]}"; do
    gen_dataset_one "$N" "$alpha_fixed" "$MASTER_SEED" || true
  done

  # 3) Write job manifest.
  local manifest="$DRIVER_DIR/jobs_all.tsv"
  write_job_manifest "$manifest"

  # 4) Attempts loop.
  log "Starting Experiment II (d=2) job execution..."

  for ((attempt=1; attempt<=MAX_ATTEMPTS; attempt++)); do
    local pending_light="$DRIVER_DIR/jobs_pending_light_attempt${attempt}.tsv"
    local pending_heavy="$DRIVER_DIR/jobs_pending_heavy_attempt${attempt}.tsv"

    write_pending_lists "$manifest" "$pending_light" "$pending_heavy"

    local n_light n_heavy
    n_light="$(count_lines "$pending_light")"
    n_heavy="$(count_lines "$pending_heavy")"

    log "Attempt ${attempt}/${MAX_ATTEMPTS}: pending light=${n_light} heavy=${n_heavy}"

    if [[ "$n_light" -eq 0 && "$n_heavy" -eq 0 ]]; then
      break
    fi

    # Run heavy first (more memory-hungry), then light.
    run_jobs_from_tsv "$pending_heavy" "$PARALLEL_HEAVY" "$attempt"
    run_jobs_from_tsv "$pending_light" "$PARALLEL_LIGHT" "$attempt"
  done

  # 5) Final summary.
  log "Collecting final status summary..."
  local summary="$DRIVER_DIR/final_summary.txt"
  write_final_summary "$manifest" "$summary"

  local total done
  total="$(count_lines "$manifest")"
  done=0
  local exp dataset t label method variant
  while IFS=$'\t' read -r exp dataset t label method variant; do
    [[ -z "${exp:-}" ]] && continue
    local out_dir
    out_dir="$(out_dir_for_job "$exp" "$dataset" "$t" "$label" "$method" "$variant")"
    if [[ -f "$out_dir/DONE" ]]; then
      done=$((done+1))
    fi
  done <"$manifest"

  log "DONE. Summary written to: $summary"
  log "Finished jobs: $done / $total"
  if [[ "$done" -lt "$total" ]]; then
    warn "There are $((total-done)) failed/skipped jobs. See: $summary"
    exit 0
  fi

  exit 0
}

SJS_RUN=""  # set in main() after build
main "$@"