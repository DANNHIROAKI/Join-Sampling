#!/usr/bin/env bash
# run_all.sh (place at repo root)
#
# Run EXP-0 .. EXP-7 once, in order.
#
# Conventions (aligned with your 5 global requirements):
#   - Each exp runner should write artifacts to:  run/temp/expX/  first
#     then overwrite-sync to:                     results/raw/expX/
#   - Builds live under:                          build/<type>/
#   - This script writes its own logs to:         run/temp/run_all/
#     then overwrites:                            results/raw/run_all/
#
# Usage:
#   bash ./run_all.sh
#
# Optional overrides (env vars):
#   RUN_ALL_BUILD_TYPE=Release            # default: Release
#   RUN_ALL_THREADS=1                     # empty => do not override
#   RUN_ALL_CONTINUE_ON_FAIL=1            # default: 0 (stop at first failure)
#   RUN_EXPS="exp0 exp1 exp2"            # subset (default: exp0..exp7)

set -euo pipefail
IFS=$'\n\t'

trap 'echo -e "[run_all][FATAL] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}"

# -------------------------
# Config
# -------------------------
BUILD_TYPE="${RUN_ALL_BUILD_TYPE:-Release}"
THREADS_OVERRIDE="${RUN_ALL_THREADS:-}"            # empty => do not override
CONTINUE_ON_FAIL="${RUN_ALL_CONTINUE_ON_FAIL:-0}"  # 1 => keep going
RUN_EXPS="${RUN_EXPS:-exp0 exp1 exp2 exp3 exp4 exp5 exp6 exp7}"

TEMP_DIR="${ROOT_DIR}/run/temp/run_all"
OUT_DIR="${ROOT_DIR}/results/raw/run_all"

# Thread caps for common libs (fairness / reproducibility)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# -------------------------
# Helpers
# -------------------------
log() { echo -e "[run_all] $*"; }
die() { echo -e "[run_all][FATAL] $*" >&2; exit 1; }

check_exp_token() {
  local e="$1"
  case "${e}" in
    exp0|exp1|exp2|exp3|exp4|exp5|exp6|exp7) return 0 ;;
    *) return 1 ;;
  esac
}

run_one() {
  local exp="$1"; shift
  local cmd_desc="$1"; shift

  local log_file="${TEMP_DIR}/${exp}.log"

  local start_iso end_iso t0 t1 elapsed rc
  start_iso="$(date -Is || date)"
  t0="$(date +%s)"

  {
    echo "[${start_iso}] [run_all] ===== START ${exp} ====="
    echo "[${start_iso}] [run_all] cmd: ${cmd_desc}"
  } | tee "${log_file}"

  # Run (capture exit code even with tee)
  set +e
  ( "$@" ) 2>&1 | tee -a "${log_file}"
  rc="${PIPESTATUS[0]}"
  set -e

  end_iso="$(date -Is || date)"
  t1="$(date +%s)"
  elapsed=$(( t1 - t0 ))

  echo "[${end_iso}] [run_all] ===== END ${exp} (exit=${rc}, elapsed=${elapsed}s) =====" | tee -a "${log_file}"
  echo "${exp},${start_iso},${end_iso},${elapsed},${rc},$(basename "${log_file}")" >> "${SUMMARY_CSV}"

  return "${rc}"
}

# -------------------------
# Preflight + init dirs
# -------------------------
log "Repo root        : ${ROOT_DIR}"
log "Build type       : ${BUILD_TYPE}"
log "Threads override : ${THREADS_OVERRIDE:-<none>}"
log "Continue on fail : ${CONTINUE_ON_FAIL}"
log "Run list         : ${RUN_EXPS}"

# Validate tokens + existence of per-exp scripts
for e in ${RUN_EXPS}; do
  check_exp_token "${e}" || die "Unknown exp token in RUN_EXPS: '${e}' (allowed: exp0..exp7)"
  script_path="${ROOT_DIR}/run/run_${e}.sh"
  [[ -f "${script_path}" ]] || die "Missing runner script: ${script_path}"
done

# Reset run_all temp/out (overwrite)
rm -rf "${TEMP_DIR}"
mkdir -p "${TEMP_DIR}"
rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

SUMMARY_CSV="${TEMP_DIR}/summary.csv"
echo "exp,start_iso,end_iso,elapsed_sec,exit_code,log_file" > "${SUMMARY_CSV}"

# Record manifest (reproducibility)
{
  echo "date=$(date -Is || date)"
  echo "repo_root=${ROOT_DIR}"
  echo "build_type=${BUILD_TYPE}"
  echo "threads_override=${THREADS_OVERRIDE:-}"
  echo "continue_on_fail=${CONTINUE_ON_FAIL}"
  echo "run_exps=${RUN_EXPS}"
  echo
  echo "uname:"; uname -a || true
  echo
  echo "git:";
  if command -v git >/dev/null 2>&1 && [[ -d "${ROOT_DIR}/.git" ]]; then
    git rev-parse HEAD || true
    git status --porcelain || true
  else
    echo "(no git repo detected)"
  fi
} > "${TEMP_DIR}/manifest.txt"

# -------------------------
# Run all experiments
# -------------------------
FAILED=0

for e in ${RUN_EXPS}; do
  log "------------------------------------------------------------"
  log "Running ${e} ..."

  set +e
  case "${e}" in
    exp0)
      if [[ -n "${THREADS_OVERRIDE}" ]]; then
        run_one "${e}" "env BUILD_TYPE=${BUILD_TYPE} THREADS=${THREADS_OVERRIDE} bash run/run_exp0.sh" \
          env BUILD_TYPE="${BUILD_TYPE}" THREADS="${THREADS_OVERRIDE}" \
          bash "${ROOT_DIR}/run/run_exp0.sh"
      else
        run_one "${e}" "env BUILD_TYPE=${BUILD_TYPE} bash run/run_exp0.sh" \
          env BUILD_TYPE="${BUILD_TYPE}" \
          bash "${ROOT_DIR}/run/run_exp0.sh"
      fi
      ;;

    exp1)
      if [[ -n "${THREADS_OVERRIDE}" ]]; then
        run_one "${e}" "env BUILD_TYPE=${BUILD_TYPE} THREADS=${THREADS_OVERRIDE} bash run/run_exp1.sh" \
          env BUILD_TYPE="${BUILD_TYPE}" THREADS="${THREADS_OVERRIDE}" \
          bash "${ROOT_DIR}/run/run_exp1.sh"
      else
        run_one "${e}" "env BUILD_TYPE=${BUILD_TYPE} bash run/run_exp1.sh" \
          env BUILD_TYPE="${BUILD_TYPE}" \
          bash "${ROOT_DIR}/run/run_exp1.sh"
      fi
      ;;

    exp2)
      # exp2 uses CLI flags for build_type / threads.
      if [[ -n "${THREADS_OVERRIDE}" ]]; then
        run_one "${e}" "bash run/run_exp2.sh --build_type ${BUILD_TYPE} --threads ${THREADS_OVERRIDE}" \
          bash "${ROOT_DIR}/run/run_exp2.sh" --build_type "${BUILD_TYPE}" --threads "${THREADS_OVERRIDE}"
      else
        run_one "${e}" "bash run/run_exp2.sh --build_type ${BUILD_TYPE}" \
          bash "${ROOT_DIR}/run/run_exp2.sh" --build_type "${BUILD_TYPE}"
      fi
      ;;

    exp3)
      # exp3 uses env vars EXP3_BUILD_TYPE / EXP3_THREADS.
      if [[ -n "${THREADS_OVERRIDE}" ]]; then
        run_one "${e}" "env EXP3_BUILD_TYPE=${BUILD_TYPE} EXP3_THREADS=${THREADS_OVERRIDE} bash run/run_exp3.sh" \
          env EXP3_BUILD_TYPE="${BUILD_TYPE}" EXP3_THREADS="${THREADS_OVERRIDE}" \
          bash "${ROOT_DIR}/run/run_exp3.sh"
      else
        run_one "${e}" "env EXP3_BUILD_TYPE=${BUILD_TYPE} bash run/run_exp3.sh" \
          env EXP3_BUILD_TYPE="${BUILD_TYPE}" \
          bash "${ROOT_DIR}/run/run_exp3.sh"
      fi
      ;;

    exp4)
      if [[ -n "${THREADS_OVERRIDE}" ]]; then
        run_one "${e}" "env BUILD_TYPE=${BUILD_TYPE} THREADS=${THREADS_OVERRIDE} bash run/run_exp4.sh" \
          env BUILD_TYPE="${BUILD_TYPE}" THREADS="${THREADS_OVERRIDE}" \
          bash "${ROOT_DIR}/run/run_exp4.sh"
      else
        run_one "${e}" "env BUILD_TYPE=${BUILD_TYPE} bash run/run_exp4.sh" \
          env BUILD_TYPE="${BUILD_TYPE}" \
          bash "${ROOT_DIR}/run/run_exp4.sh"
      fi
      ;;

    exp5)
      # exp5 uses env vars EXP5_BUILD_TYPE / EXP5_THREADS.
      if [[ -n "${THREADS_OVERRIDE}" ]]; then
        run_one "${e}" "env EXP5_BUILD_TYPE=${BUILD_TYPE} EXP5_THREADS=${THREADS_OVERRIDE} bash run/run_exp5.sh" \
          env EXP5_BUILD_TYPE="${BUILD_TYPE}" EXP5_THREADS="${THREADS_OVERRIDE}" \
          bash "${ROOT_DIR}/run/run_exp5.sh"
      else
        run_one "${e}" "env EXP5_BUILD_TYPE=${BUILD_TYPE} bash run/run_exp5.sh" \
          env EXP5_BUILD_TYPE="${BUILD_TYPE}" \
          bash "${ROOT_DIR}/run/run_exp5.sh"
      fi
      ;;

    exp6)
      # exp6 uses env vars EXP6_BUILD_TYPE / EXP6_THREADS.
      if [[ -n "${THREADS_OVERRIDE}" ]]; then
        run_one "${e}" "env EXP6_BUILD_TYPE=${BUILD_TYPE} EXP6_THREADS=${THREADS_OVERRIDE} bash run/run_exp6.sh" \
          env EXP6_BUILD_TYPE="${BUILD_TYPE}" EXP6_THREADS="${THREADS_OVERRIDE}" \
          bash "${ROOT_DIR}/run/run_exp6.sh"
      else
        run_one "${e}" "env EXP6_BUILD_TYPE=${BUILD_TYPE} bash run/run_exp6.sh" \
          env EXP6_BUILD_TYPE="${BUILD_TYPE}" \
          bash "${ROOT_DIR}/run/run_exp6.sh"
      fi
      ;;

    exp7)
      if [[ -n "${THREADS_OVERRIDE}" ]]; then
        run_one "${e}" "env BUILD_TYPE=${BUILD_TYPE} THREADS=${THREADS_OVERRIDE} bash run/run_exp7.sh" \
          env BUILD_TYPE="${BUILD_TYPE}" THREADS="${THREADS_OVERRIDE}" \
          bash "${ROOT_DIR}/run/run_exp7.sh"
      else
        run_one "${e}" "env BUILD_TYPE=${BUILD_TYPE} bash run/run_exp7.sh" \
          env BUILD_TYPE="${BUILD_TYPE}" \
          bash "${ROOT_DIR}/run/run_exp7.sh"
      fi
      ;;
  esac
  rc=$?
  set -e

  if [[ "${rc}" -ne 0 ]]; then
    FAILED=1
    log "[WARN] ${e} failed (exit=${rc}). Log: ${TEMP_DIR}/${e}.log"
    if [[ "${CONTINUE_ON_FAIL}" != "1" ]]; then
      log "Stop on first failure (RUN_ALL_CONTINUE_ON_FAIL=0)."
      break
    fi
  else
    log "${e} finished OK."
  fi

done

# -------------------------
# Sync run_all logs -> results/raw/run_all (overwrite)
# -------------------------
cp -a "${TEMP_DIR}/." "${OUT_DIR}/"

log "============================================================"
log "run_all finished."
log "Summary CSV : ${OUT_DIR}/summary.csv"
log "Manifest    : ${OUT_DIR}/manifest.txt"
log "Logs        : ${OUT_DIR}/exp*.log"
log "(Each experiment's own results are under results/raw/exp0 .. results/raw/exp7)"

if [[ "${FAILED}" -eq 0 ]]; then
  log "ALL EXPERIMENTS PASSED ✅"
  exit 0
else
  log "Some experiments FAILED ❗ (see summary.csv + individual logs)"
  exit 1
fi
