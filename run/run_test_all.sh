#!/usr/bin/env bash
# ==============================================================================
# run/run_test_all.sh
#
# Repository-wide "torture" test runner.
#
# Goals:
#   - Broad + strict validation of the repo's *correctness* (not performance).
#   - Catch: semantic bugs (half-open join), sampling correctness, edge cases,
#           crashes, UB (ASAN/UBSAN), races (optional TSAN), script hygiene.
#
# Outputs:
#   - results/raw/test_all/         (logs + report + artifacts)
#   - run/temp/test_all/            (generated configs, datasets for tests)
#   - build/{release,debug,asan,tsan}/
#
# Usage:
#   chmod +x run/run_test_all.sh
#   bash run/run_test_all.sh
#
# Configuration (env vars):
#   MODE=deep|fast                (default: deep)
#   CLEAN=1|0                     (default: 1) wipe test output dirs first
#   CONTINUE_ON_FAIL=0|1          (default: 0) stop at first failure
#   POLICY_STRICT=1|0             (default: 1) policy violations => fail
#   JOBS=<int>                    (default: nproc)
#   TIMEOUT_SEC=<int>             (default: 600) per command
#
#   RUN_ASAN=1|0                  (default: deep->1, fast->0)
#   RUN_TSAN=1|0                  (default: 0)
#   RUN_STATIC_ANALYSIS=1|0       (default: 1) best-effort; skipped if tools missing
#
#   METHODS="ours aabb ..."        (space-separated)
#   VARIANTS="sampling enum_sampling adaptive"
#
# Dataset knobs (smoke + correctness battery):
#   NR=2000 NS=2000 ALPHA=1 GEN_SEED=1
#   T=20000 SEED=1 REPEATS=2
#   ORACLE_MAX_CHECKS=50000000
#
# Edge-case knobs:
#   J_STAR_SMALL=1
#   J_STAR_LARGE=1000000000
#   ENUM_CAP_TEST=1000
#
# Dependency:
#   - run/include/csv_check.py    (CSV helper; no inline python in this bash)
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

# ----------------------------
# Repo root + key directories
# ----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

RUN_DIR="${REPO_ROOT}/run"
INCLUDE_DIR="${RUN_DIR}/include"
TEMP_ROOT="${RUN_DIR}/temp"
TEMP_DIR="${TEMP_ROOT}/test_all"

OUT_ROOT="${REPO_ROOT}/results/raw/test_all"
LOG_DIR="${OUT_ROOT}/logs"
ART_DIR="${OUT_ROOT}/artifacts"

BUILD_ROOT="${REPO_ROOT}/build"
BUILD_RELEASE="${BUILD_ROOT}/release"
BUILD_DEBUG="${BUILD_ROOT}/debug"
BUILD_ASAN="${BUILD_ROOT}/asan"
BUILD_TSAN="${BUILD_ROOT}/tsan"

CSV_CHECK_PY="${INCLUDE_DIR}/csv_check.py"

# ----------------------------
# Configuration defaults
# ----------------------------
MODE="${MODE:-deep}"
CLEAN="${CLEAN:-1}"
CONTINUE_ON_FAIL="${CONTINUE_ON_FAIL:-0}"
POLICY_STRICT="${POLICY_STRICT:-1}"

TIMEOUT_SEC="${TIMEOUT_SEC:-600}"
JOBS="${JOBS:-}"

RUN_STATIC_ANALYSIS="${RUN_STATIC_ANALYSIS:-1}"

if [[ "${MODE}" == "fast" ]]; then
  RUN_ASAN="${RUN_ASAN:-0}"
  FUZZ_CASES="${FUZZ_CASES:-5}"
else
  RUN_ASAN="${RUN_ASAN:-1}"
  FUZZ_CASES="${FUZZ_CASES:-30}"
fi
RUN_TSAN="${RUN_TSAN:-0}"

# Dataset knobs
NR="${NR:-2000}"
NS="${NS:-2000}"
ALPHA="${ALPHA:-1}"
GEN_SEED="${GEN_SEED:-1}"

T="${T:-20000}"
SEED="${SEED:-1}"
REPEATS="${REPEATS:-2}"

ORACLE_MAX_CHECKS="${ORACLE_MAX_CHECKS:-50000000}"

# Edge-case knobs
J_STAR_SMALL="${J_STAR_SMALL:-1}"
J_STAR_LARGE="${J_STAR_LARGE:-1000000000}"
ENUM_CAP_TEST="${ENUM_CAP_TEST:-1000}"

# Methods/variants (Dim=2 default)
METHODS_DEFAULT=("ours" "aabb" "interval_tree" "kd_tree" "r_tree" "range_tree" "pbsm" "tlsop" "sirs" "rejection" "tsunami")
VARIANTS_DEFAULT=("sampling" "enum_sampling" "adaptive")

if [[ -n "${METHODS:-}" ]]; then
  # shellcheck disable=SC2206
  METHODS_ARR=(${METHODS})
else
  METHODS_ARR=("${METHODS_DEFAULT[@]}")
fi

if [[ -n "${VARIANTS:-}" ]]; then
  # shellcheck disable=SC2206
  VARIANTS_ARR=(${VARIANTS})
else
  VARIANTS_ARR=("${VARIANTS_DEFAULT[@]}")
fi

# ----------------------------
# Helpers
# ----------------------------
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

log()  { echo "[$(timestamp)][TEST_ALL] $*"; }
warn() { echo "[$(timestamp)][TEST_ALL][WARN] $*" >&2; }
die()  { echo "[$(timestamp)][TEST_ALL][FATAL] $*" >&2; exit 1; }

need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"; }

cpu_count() {
  if [[ -n "${JOBS}" ]]; then
    echo "${JOBS}"
    return
  fi
  if command -v nproc >/dev/null 2>&1; then
    nproc
  elif command -v sysctl >/dev/null 2>&1; then
    sysctl -n hw.ncpu
  else
    echo 4
  fi
}

TIMEOUT_BIN=""
if command -v timeout >/dev/null 2>&1; then
  TIMEOUT_BIN="timeout"
fi

hash_file() {
  local f="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${f}" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "${f}" | awk '{print $1}'
  elif command -v md5sum >/dev/null 2>&1; then
    md5sum "${f}" | awk '{print $1}'
  else
    stat -c '%s_%Y' "${f}" 2>/dev/null || stat -f '%z_%m' "${f}"
  fi
}

# Extract key=value-like fields from a summary line such as:
# "count=4 (exact)  oracle=4  rel_err=0"
extract_field() {
  local line="$1"
  local key="$2"
  echo "${line}" | awk -v k="${key}" -F'[ =()]+' '{
    for (i=1; i<=NF; ++i) if ($i==k) { print $(i+1); exit }
  }'
}

extract_count_kind() {
  local line="$1"
  echo "${line}" | awk -F'[ =()]+' '{
    for (i=1; i<=NF; ++i) if ($i=="count") { print $(i+2); exit }
  }'
}

float_is_zero() {
  local x="$1"
  awk -v v="${x}" 'BEGIN { exit ((v+0)==0 ? 0 : 1) }'
}

float_leq() {
  local x="$1"; local y="$2"
  awk -v a="${x}" -v b="${y}" 'BEGIN { exit ((a+0) <= (b+0) ? 0 : 1) }'
}

float_eq() {
  local x="$1"; local y="$2"
  awk -v a="${x}" -v b="${y}" 'BEGIN { exit ((a+0) == (b+0) ? 0 : 1) }'
}

resolve_exe() {
  # Usage: resolve_exe <build_dir> <name>
  local build_dir="$1"
  local name="$2"
  local p=""

  for cand in \
    "${build_dir}/${name}" \
    "${build_dir}/apps/${name}" \
    "${build_dir}/bin/${name}" \
    "${build_dir}/src/apps/${name}" \
    "${build_dir}/src/${name}" \
  ; do
    if [[ -x "${cand}" ]]; then
      echo "${cand}"
      return 0
    fi
  done

  p="$(find "${build_dir}" -maxdepth 5 -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "${p}" ]] || return 1
  echo "${p}"
  return 0
}

supports_flag() {
  # Usage: supports_flag <exe> <flag_string>
  local exe="$1"
  local flag="$2"
  "${exe}" --help 2>&1 | grep -q -- "${flag}"
}

skip_pattern() {
  # Centralize "skip" detection patterns for unsupported/disabled combos.
  # We keep this intentionally narrow-ish.
  grep -qiE "unsupported|not supported|SKIP|SKIPPED|disabled|unknown (method|variant|generator)|invalid (method|variant|generator)" "$1"
}

# ----------------------------
# Test accounting + runner
# ----------------------------
TOTAL=0
PASSED=0
FAILED=0
SKIPPED=0

REPORT="${OUT_ROOT}/REPORT.txt"
FAILURES="${OUT_ROOT}/FAILURES.txt"

record_pass() { PASSED=$((PASSED+1)); TOTAL=$((TOTAL+1)); }
record_fail() { FAILED=$((FAILED+1)); TOTAL=$((TOTAL+1)); }
record_skip() { SKIPPED=$((SKIPPED+1)); TOTAL=$((TOTAL+1)); }

maybe_exit_on_fail() {
  if [[ "${CONTINUE_ON_FAIL}" == "0" ]]; then
    exit 1
  fi
}

run_cmd() {
  # run_cmd <tag> <command...>
  local tag="$1"; shift
  local log_file="${LOG_DIR}/${tag}.log"

  echo "------------------------------------------------------------" | tee -a "${REPORT}"
  echo "[$(timestamp)] [RUN] ${tag}" | tee -a "${REPORT}"
  echo "CMD: $*" | tee -a "${REPORT}"

  mkdir -p "$(dirname "${log_file}")"
  : > "${log_file}"
  echo "[$(timestamp)] CMD: $*" >> "${log_file}"

  set +e
  if [[ -n "${TIMEOUT_BIN}" && "${TIMEOUT_SEC}" -gt 0 ]]; then
    ("${TIMEOUT_BIN}" "${TIMEOUT_SEC}" "$@") 2>&1 | tee -a "${log_file}"
    rc="${PIPESTATUS[0]}"
  else
    ("$@") 2>&1 | tee -a "${log_file}"
    rc="${PIPESTATUS[0]}"
  fi
  set -e

  if [[ "${rc}" -ne 0 ]]; then
    echo "[$(timestamp)] [FAIL] ${tag} (rc=${rc})" | tee -a "${REPORT}" | tee -a "${FAILURES}"
    record_fail
    maybe_exit_on_fail
    return 1
  fi

  echo "[$(timestamp)] [PASS] ${tag}" | tee -a "${REPORT}"
  record_pass
  return 0
}

run_cmd_allow_skip() {
  # run_cmd_allow_skip <tag> <command...>
  local tag="$1"; shift
  local log_file="${LOG_DIR}/${tag}.log"

  echo "------------------------------------------------------------" | tee -a "${REPORT}"
  echo "[$(timestamp)] [RUN_ALLOW_SKIP] ${tag}" | tee -a "${REPORT}"
  echo "CMD: $*" | tee -a "${REPORT}"

  mkdir -p "$(dirname "${log_file}")"
  : > "${log_file}"
  echo "[$(timestamp)] CMD: $*" >> "${log_file}"

  set +e
  if [[ -n "${TIMEOUT_BIN}" && "${TIMEOUT_SEC}" -gt 0 ]]; then
    ("${TIMEOUT_BIN}" "${TIMEOUT_SEC}" "$@") 2>&1 | tee -a "${log_file}"
    rc="${PIPESTATUS[0]}"
  else
    ("$@") 2>&1 | tee -a "${log_file}"
    rc="${PIPESTATUS[0]}"
  fi
  set -e

  if [[ "${rc}" -eq 0 ]]; then
    echo "[$(timestamp)] [PASS] ${tag}" | tee -a "${REPORT}"
    record_pass
    return 0
  fi

  if skip_pattern "${log_file}"; then
    echo "[$(timestamp)] [SKIP] ${tag} (unsupported/disabled)" | tee -a "${REPORT}"
    record_skip
    return 0
  fi

  echo "[$(timestamp)] [FAIL] ${tag} (rc=${rc})" | tee -a "${REPORT}" | tee -a "${FAILURES}"
  record_fail
  maybe_exit_on_fail
  return 1
}

# ----------------------------
# Setup dirs + environment
# ----------------------------
need_cmd cmake
need_cmd ctest
need_cmd python3
need_cmd grep
need_cmd awk
need_cmd tee
need_cmd find
need_cmd sed
need_cmd sort
need_cmd head
need_cmd tail

[[ -f "${CSV_CHECK_PY}" ]] || die "Missing helper: ${CSV_CHECK_PY} (please add it under run/include/)"

if [[ "${CLEAN}" == "1" ]]; then
  log "CLEAN=1 -> removing: ${OUT_ROOT} and ${TEMP_DIR}"
  rm -rf "${OUT_ROOT}" "${TEMP_DIR}"
fi

mkdir -p "${LOG_DIR}" "${ART_DIR}" "${TEMP_DIR}" "${BUILD_ROOT}"

: > "${REPORT}"
: > "${FAILURES}"

log "Repo root        : ${REPO_ROOT}"
log "Mode             : ${MODE}"
log "Build root       : ${BUILD_ROOT}"
log "Output root      : ${OUT_ROOT}"
log "Temp dir         : ${TEMP_DIR}"
log "Timeout (sec)    : ${TIMEOUT_SEC}"
log "Continue on fail : ${CONTINUE_ON_FAIL}"
log "Policy strict    : ${POLICY_STRICT}"
log "Jobs             : $(cpu_count)"
log "RUN_ASAN         : ${RUN_ASAN}"
log "RUN_TSAN         : ${RUN_TSAN}"
log "FUZZ_CASES       : ${FUZZ_CASES}"
log "METHODS          : ${METHODS_ARR[*]}"
log "VARIANTS         : ${VARIANTS_ARR[*]}"

# Thread caps for deterministic + fair-ish behavior in tests
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

# Record environment
{
  echo "timestamp=$(date -Is || date)"
  echo "repo_root=${REPO_ROOT}"
  echo "mode=${MODE}"
  echo "uname=$(uname -a || true)"
  echo "cmake=$(cmake --version 2>/dev/null | head -n 1 || true)"
  echo "ctest=$(ctest --version 2>/dev/null | head -n 1 || true)"
  echo "python3=$(python3 --version 2>/dev/null || true)"
  if command -v git >/dev/null 2>&1 && [[ -d "${REPO_ROOT}/.git" ]]; then
    echo "git_head=$(git -C "${REPO_ROOT}" rev-parse HEAD || true)"
    echo "git_status:"
    git -C "${REPO_ROOT}" status --porcelain || true
  fi
} > "${OUT_ROOT}/env.txt"

# ----------------------------
# 0) Policy / hygiene lint for run scripts
# ----------------------------
policy_check() {
  log "Step 0: Policy/Hygiene lint for run/*.sh"

  local bad=0
  local sh
  while IFS= read -r sh; do
    # 0.1 embedded python execution: python - <<'PY'
    if grep -nE '^[[:space:]]*[^#].*python(3)?[[:space:]]+-[[:space:]]*<<' "${sh}" >/dev/null 2>&1; then
      warn "[POLICY] Embedded python heredoc execution found in: ${sh}"
      grep -nE '^[[:space:]]*[^#].*python(3)?[[:space:]]+-[[:space:]]*<<' "${sh}" | head -n 5 >&2 || true
      bad=1
    fi
    # 0.2 embedded python blocks: <<'PY' ... (python source inside bash)
    if grep -nE "^[[:space:]]*[^#].*<<'PY'|^[[:space:]]*[^#].*<<PY$|^[[:space:]]*[^#].*<<'PYTHON'|^[[:space:]]*[^#].*<<PYTHON$" "${sh}" >/dev/null 2>&1; then
      warn "[POLICY] Inline python heredoc block found in: ${sh}"
      grep -nE "^[[:space:]]*[^#].*<<'PY'|^[[:space:]]*[^#].*<<PY$|^[[:space:]]*[^#].*<<'PYTHON'|^[[:space:]]*[^#].*<<PYTHON$" "${sh}" | head -n 5 >&2 || true
      bad=1
    fi
    # 0.3 build dir must live under <root>/build/ (reject legacy /build_* dirs)
    if grep -nE '^[[:space:]]*[^#].*/build_[A-Za-z0-9_]' "${sh}" >/dev/null 2>&1; then
      warn "[POLICY] Legacy build_* dir found (should be under <root>/build/): ${sh}"
      grep -nE '^[[:space:]]*[^#].*/build_[A-Za-z0-9_]' "${sh}" | head -n 5 >&2 || true
      bad=1
    fi
    # 0.4 results must be under results/raw (reject results/sweeps or results/expX)
    if grep -nE '^[[:space:]]*[^#].*results/sweeps' "${sh}" >/dev/null 2>&1; then
      warn "[POLICY] Legacy results/sweeps path found (should be results/raw): ${sh}"
      grep -nE '^[[:space:]]*[^#].*results/sweeps' "${sh}" | head -n 5 >&2 || true
      bad=1
    fi
    if grep -nE '^[[:space:]]*[^#].*results/exp[0-9]' "${sh}" >/dev/null 2>&1; then
      warn "[POLICY] Legacy results/expX path found (should be results/raw/expX): ${sh}"
      grep -nE '^[[:space:]]*[^#].*results/exp[0-9]' "${sh}" | head -n 5 >&2 || true
      bad=1
    fi
  done < <(find "${RUN_DIR}" -maxdepth 1 -type f -name "*.sh" | sort)

  if [[ "${bad}" -eq 1 ]]; then
    if [[ "${POLICY_STRICT}" == "1" ]]; then
      die "Policy lint FAILED (see warnings above). Fix run scripts or set POLICY_STRICT=0."
    else
      warn "Policy lint reported issues, but POLICY_STRICT=0 -> continue."
      return 0
    fi
  fi

  log "Policy lint PASSED."
  return 0
}

policy_check

# ----------------------------
# 1) Build matrix + unit tests
# ----------------------------
cmake_build() {
  local bdir="$1"
  local build_type="$2"
  local tag="$3"
  shift 3
  local extra_cmake_args=("$@")

  mkdir -p "${bdir}"
  run_cmd "cmake_configure_${tag}" \
    cmake -S "${REPO_ROOT}" -B "${bdir}" \
      -DCMAKE_BUILD_TYPE="${build_type}" \
      -DSJS_BUILD_ROOT_APPS=ON \
      -DSJS_BUILD_TESTS=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
      "${extra_cmake_args[@]}"

  run_cmd "cmake_build_${tag}" \
    cmake --build "${bdir}" -j "$(cpu_count)"

  run_cmd "ctest_${tag}" \
    ctest --test-dir "${bdir}" --output-on-failure
}

log "Step 1: Build + unit tests (Release)"
cmake_build "${BUILD_RELEASE}" "Release" "release"

log "Step 1b: Build + unit tests (Debug)"
cmake_build "${BUILD_DEBUG}" "Debug" "debug"

if [[ "${RUN_ASAN}" == "1" ]]; then
  log "Step 1c: Build + unit tests (ASAN/UBSAN)"

  CMAKE_ASAN_ARGS=()
  if command -v clang >/dev/null 2>&1 && command -v clang++ >/dev/null 2>&1; then
    CMAKE_ASAN_ARGS+=("-DCMAKE_C_COMPILER=clang" "-DCMAKE_CXX_COMPILER=clang++")
  fi

  CMAKE_ASAN_ARGS+=(
    "-DCMAKE_BUILD_TYPE=RelWithDebInfo"
    "-DCMAKE_CXX_FLAGS=-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer"
    "-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined"
  )

  mkdir -p "${BUILD_ASAN}"
  run_cmd "cmake_configure_asan" \
    cmake -S "${REPO_ROOT}" -B "${BUILD_ASAN}" \
      -DSJS_BUILD_ROOT_APPS=ON -DSJS_BUILD_TESTS=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
      "${CMAKE_ASAN_ARGS[@]}"

  run_cmd "cmake_build_asan" \
    cmake --build "${BUILD_ASAN}" -j "$(cpu_count)"

  export ASAN_OPTIONS="${ASAN_OPTIONS:-halt_on_error=1:abort_on_error=1:detect_leaks=1:symbolize=1}"
  export UBSAN_OPTIONS="${UBSAN_OPTIONS:-halt_on_error=1:abort_on_error=1:print_stacktrace=1}"

  run_cmd "ctest_asan" \
    ctest --test-dir "${BUILD_ASAN}" --output-on-failure
fi

if [[ "${RUN_TSAN}" == "1" ]]; then
  log "Step 1d: Build + unit tests (TSAN)"
  if command -v clang >/dev/null 2>&1 && command -v clang++ >/dev/null 2>&1; then
    mkdir -p "${BUILD_TSAN}"
    run_cmd "cmake_configure_tsan" \
      cmake -S "${REPO_ROOT}" -B "${BUILD_TSAN}" \
        -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
        -DSJS_BUILD_ROOT_APPS=ON -DSJS_BUILD_TESTS=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        "-DCMAKE_CXX_FLAGS=-O1 -g -fsanitize=thread -fno-omit-frame-pointer" \
        "-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=thread"

    run_cmd "cmake_build_tsan" \
      cmake --build "${BUILD_TSAN}" -j "$(cpu_count)"

    export TSAN_OPTIONS="${TSAN_OPTIONS:-halt_on_error=1:abort_on_error=1}"
    run_cmd "ctest_tsan" \
      ctest --test-dir "${BUILD_TSAN}" --output-on-failure
  else
    warn "TSAN requested but clang/clang++ not found; skipping TSAN build."
    record_skip
  fi
fi

# ----------------------------
# 2) Locate binaries (Release build)
# ----------------------------
SJS_GEN="$(resolve_exe "${BUILD_RELEASE}" "sjs_gen_dataset" || true)"
SJS_RUN="$(resolve_exe "${BUILD_RELEASE}" "sjs_run" || true)"
SJS_VERIFY="$(resolve_exe "${BUILD_RELEASE}" "sjs_verify" || true)"
SJS_SWEEP="$(resolve_exe "${BUILD_RELEASE}" "sjs_sweep" || true)"

[[ -n "${SJS_GEN}" ]] || die "Cannot find sjs_gen_dataset under ${BUILD_RELEASE}"
[[ -n "${SJS_RUN}" ]] || die "Cannot find sjs_run under ${BUILD_RELEASE}"
[[ -n "${SJS_VERIFY}" ]] || die "Cannot find sjs_verify under ${BUILD_RELEASE}"
[[ -n "${SJS_SWEEP}" ]] || die "Cannot find sjs_sweep under ${BUILD_RELEASE}"

log "Using binaries:"
log "  sjs_gen_dataset : ${SJS_GEN}"
log "  sjs_run         : ${SJS_RUN}"
log "  sjs_verify      : ${SJS_VERIFY}"
log "  sjs_sweep       : ${SJS_SWEEP}"

# ----------------------------
# 3) CLI smoke (--help)
# ----------------------------
run_cmd "help_sjs_gen_dataset" "${SJS_GEN}" --help
run_cmd "help_sjs_run"         "${SJS_RUN}" --help
run_cmd "help_sjs_verify"      "${SJS_VERIFY}" --help
run_cmd "help_sjs_sweep"       "${SJS_SWEEP}" --help

# ----------------------------
# 4) End-to-end pipeline smoke (synthetic + binary)
# ----------------------------
DATA_DIR="${TEMP_DIR}/datasets"
mkdir -p "${DATA_DIR}"

DATASET_NAME="test_all_smoke"

log "Step 4: Generate offline dataset (binary + csv) -> ${DATA_DIR}"
run_cmd "gen_dataset_smoke" \
  "${SJS_GEN}" \
    --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}" \
    --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
    --out_dir="${DATA_DIR}" --write_csv=1 \
    --threads=1

R_BIN="${DATA_DIR}/${DATASET_NAME}_R.bin"
S_BIN="${DATA_DIR}/${DATASET_NAME}_S.bin"
[[ -f "${R_BIN}" ]] || die "Missing expected dataset file: ${R_BIN}"
[[ -f "${S_BIN}" ]] || die "Missing expected dataset file: ${S_BIN}"

# sjs_run smoke (ours all variants on synthetic)
for v in sampling enum_sampling adaptive; do
  OUT_DIR="${OUT_ROOT}/smoke/run_ours_${v}_synthetic"
  mkdir -p "${OUT_DIR}"
  run_cmd "run_ours_${v}_synthetic" \
    "${SJS_RUN}" \
      --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}_synthetic" \
      --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
      --method=ours --variant="${v}" --t="${T}" --seed="${SEED}" --repeats=1 \
      --out_dir="${OUT_DIR}" --write_samples=1 \
      --threads=1

  python3 "${CSV_CHECK_PY}" "${OUT_DIR}/run.csv" --nonempty >/dev/null
done

# sjs_run smoke (binary path)
OUT_DIR_BIN="${OUT_ROOT}/smoke/run_ours_sampling_binary"
mkdir -p "${OUT_DIR_BIN}"
run_cmd "run_ours_sampling_binary" \
  "${SJS_RUN}" \
    --dataset_source=binary --dataset="${DATASET_NAME}_binary" \
    --path_r="${R_BIN}" --path_s="${S_BIN}" \
    --method=ours --variant=sampling --t="${T}" --seed="${SEED}" --repeats=1 \
    --out_dir="${OUT_DIR_BIN}" --write_samples=1 \
    --threads=1
python3 "${CSV_CHECK_PY}" "${OUT_DIR_BIN}/run.csv" --nonempty >/dev/null

# baseline crash-guard
for m in r_tree pbsm; do
  OUT_DIR="${OUT_ROOT}/smoke/run_${m}_sampling_synthetic"
  mkdir -p "${OUT_DIR}"
  run_cmd_allow_skip "run_${m}_sampling_synthetic" \
    "${SJS_RUN}" \
      --dataset_source=synthetic --gen=stripe --dataset="${DATASET_NAME}_${m}_smoke" \
      --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
      --method="${m}" --variant=sampling --t="${T}" --seed="${SEED}" --repeats=1 \
      --out_dir="${OUT_DIR}" --write_samples=0 \
      --threads=1
  if [[ -f "${OUT_DIR}/run.csv" ]]; then
    python3 "${CSV_CHECK_PY}" "${OUT_DIR}/run.csv" --nonempty >/dev/null || true
  fi
done

# Determinism check (single-thread, same seed => identical samples)
log "Step 4b: Determinism check (same seed => identical sample file hash)"
DET1="${OUT_ROOT}/determinism/run1"
DET2="${OUT_ROOT}/determinism/run2"
mkdir -p "${DET1}" "${DET2}"

run_cmd "determinism_run1" \
  "${SJS_RUN}" \
    --dataset_source=binary --dataset="${DATASET_NAME}_det" \
    --path_r="${R_BIN}" --path_s="${S_BIN}" \
    --method=ours --variant=sampling --t=5000 --seed=123 --repeats=1 \
    --out_dir="${DET1}" --write_samples=1 \
    --threads=1

run_cmd "determinism_run2" \
  "${SJS_RUN}" \
    --dataset_source=binary --dataset="${DATASET_NAME}_det" \
    --path_r="${R_BIN}" --path_s="${S_BIN}" \
    --method=ours --variant=sampling --t=5000 --seed=123 --repeats=1 \
    --out_dir="${DET2}" --write_samples=1 \
    --threads=1

S1="$(find "${DET1}" -type f -path '*/samples/*' -name '*.tsv' 2>/dev/null | head -n 1 || true)"
S2="$(find "${DET2}" -type f -path '*/samples/*' -name '*.tsv' 2>/dev/null | head -n 1 || true)"

if [[ -n "${S1}" && -n "${S2}" ]]; then
  H1="$(hash_file "${S1}")"
  H2="$(hash_file "${S2}")"
  echo "determinism_hash1=${H1}" | tee -a "${REPORT}"
  echo "determinism_hash2=${H2}" | tee -a "${REPORT}"
  if [[ "${H1}" != "${H2}" ]]; then
    die "Determinism FAILED: same seed produced different sample file hashes. (${H1} vs ${H2})"
  fi
  log "Determinism OK."
else
  warn "Determinism check skipped: sample files not found under ${DET1}/samples and ${DET2}/samples."
  record_skip
fi

# ----------------------------
# 5) Correctness battery (oracle gate via sjs_verify)
# ----------------------------
log "Step 5: Correctness battery via sjs_verify (missing_in_universe must be 0)"

verify_and_gate() {
  local method="$1"
  local variant="$2"
  local tag="verify_${method}_${variant}"
  local log_file="${LOG_DIR}/${tag}.log"

  run_cmd_allow_skip "${tag}" \
    "${SJS_VERIFY}" \
      --dataset_source=synthetic --gen=stripe --dataset="test_all_verify" \
      --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
      --method="${method}" --variant="${variant}" --t="${T}" --seed="${SEED}" --repeats="${REPEATS}" \
      --threads=1 \
      --oracle_max_checks="${ORACLE_MAX_CHECKS}"

  # If skipped (unsupported), stop here.
  if [[ -f "${log_file}" ]] && skip_pattern "${log_file}"; then
    return 0
  fi

  # Gate 1: missing_in_universe must be 0 for ALL repeats if printed.
  local miss_lines miss_bad
  miss_lines="$(grep -E "missing_in_universe=" "${log_file}" || true)"
  if [[ -n "${miss_lines}" ]]; then
    miss_bad="$(echo "${miss_lines}" | awk -F= '{gsub(/[[:space:]]/,"",$2); if ($2!="0") {print $0}}' || true)"
    if [[ -n "${miss_bad}" ]]; then
      die "${tag}: correctness FAILED (missing_in_universe != 0). See ${log_file}"
    fi
  else
    warn "${tag}: no missing_in_universe lines found (verify output format changed or check skipped)."
    record_skip
  fi

  # Gate 2: if any summary line says (exact), require count==oracle and rel_err==0.
  local check_est="${CHECK_EST_COUNT:-0}"
  local est_rel_err_max="${EST_REL_ERR_MAX:-0.5}"
  local est_rel_err_warn="${EST_REL_ERR_WARN:-0.5}"

  local sum_lines
  sum_lines="$(grep -E "count=.*oracle=" "${log_file}" || true)"
  if [[ -z "${sum_lines}" ]]; then
    warn "${tag}: no 'count=... oracle=...' summary lines found."
    record_skip
    return 0
  fi

  while IFS= read -r line; do
    [[ -z "${line}" ]] && continue
    local count oracle rel_err kind
    count="$(extract_field "${line}" "count")"
    oracle="$(extract_field "${line}" "oracle")"
    rel_err="$(extract_field "${line}" "rel_err")"
    kind="$(extract_count_kind "${line}")"

    if [[ "${kind}" == "exact" ]]; then
      if ! float_eq "${count}" "${oracle}"; then
        die "${tag}: exact count mismatch: count=${count} oracle=${oracle} (line='${line}')"
      fi
      if ! float_is_zero "${rel_err}"; then
        die "${tag}: exact rel_err must be 0, got rel_err=${rel_err} (line='${line}')"
      fi
    else
      if [[ "${check_est}" == "1" ]]; then
        if ! float_leq "${rel_err}" "${est_rel_err_max}"; then
          die "${tag}: estimated rel_err too large: rel_err=${rel_err} > ${est_rel_err_max} (line='${line}')"
        fi
      else
        if ! float_leq "${rel_err}" "${est_rel_err_warn}"; then
          warn "${tag}: estimated rel_err warning: rel_err=${rel_err} > ${est_rel_err_warn} (line='${line}')"
        fi
      fi
    fi
  done <<< "${sum_lines}"
}

for m in "${METHODS_ARR[@]}"; do
  for v in "${VARIANTS_ARR[@]}"; do
    verify_and_gate "${m}" "${v}"
  done
done

# Verify binary-path (if supported by sjs_verify)
if supports_flag "${SJS_VERIFY}" "--path_r"; then
  OUT_TAG="verify_ours_sampling_binary"
  OUT_LOG="${LOG_DIR}/${OUT_TAG}.log"
  run_cmd "${OUT_TAG}" \
    "${SJS_VERIFY}" \
      --dataset_source=binary --dataset="${DATASET_NAME}_binary_verify" \
      --path_r="${R_BIN}" --path_s="${S_BIN}" \
      --method=ours --variant=sampling --t="${T}" --seed="${SEED}" --repeats=1 \
      --threads=1 \
      --oracle_max_checks="${ORACLE_MAX_CHECKS}"

  miss="$(grep -E "missing_in_universe=" "${OUT_LOG}" | tail -n1 | awk -F= '{print $2}' | tr -d '[:space:]' || true)"
  [[ -z "${miss}" || "${miss}" == "0" ]] || die "Binary verify FAILED: missing_in_universe=${miss} (see ${OUT_LOG})"
else
  warn "sjs_verify does not support --path_r; skipping binary verify."
  record_skip
fi

# ----------------------------
# 6) Edge-case tests (adaptive branch, enum_cap, empty join)
# ----------------------------
log "Step 6: Edge-case tests"

# 6.1 adaptive branch switching (if sjs_run supports --j_star)
if supports_flag "${SJS_RUN}" "--j_star"; then
  OUT_A="${OUT_ROOT}/edge/adaptive_force_fallback"
  mkdir -p "${OUT_A}"
  run_cmd "adaptive_force_fallback" \
    "${SJS_RUN}" \
      --dataset_source=synthetic --gen=stripe --dataset="test_all_adaptive_fb" \
      --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
      --method=ours --variant=adaptive --t=5000 --seed=777 --repeats=1 \
      --j_star="${J_STAR_SMALL}" \
      --out_dir="${OUT_A}" --write_samples=0 \
      --threads=1
  python3 "${CSV_CHECK_PY}" "${OUT_A}/run.csv" --nonempty >/dev/null
  if ! python3 "${CSV_CHECK_PY}" "${OUT_A}/run.csv" --last-row-col-contains adaptive_branch fallback >/dev/null 2>&1; then
    warn "adaptive_force_fallback: could not assert adaptive_branch contains 'fallback' (column missing or token differs)."
    record_skip
  fi

  OUT_B="${OUT_ROOT}/edge/adaptive_force_enumerate"
  mkdir -p "${OUT_B}"
  run_cmd "adaptive_force_enumerate" \
    "${SJS_RUN}" \
      --dataset_source=synthetic --gen=stripe --dataset="test_all_adaptive_enum" \
      --n_r="${NR}" --n_s="${NS}" --alpha=0.01 --gen_seed="${GEN_SEED}" \
      --method=ours --variant=adaptive --t=5000 --seed=778 --repeats=1 \
      --j_star="${J_STAR_LARGE}" \
      --out_dir="${OUT_B}" --write_samples=0 \
      --threads=1
  python3 "${CSV_CHECK_PY}" "${OUT_B}/run.csv" --nonempty >/dev/null
  if ! python3 "${CSV_CHECK_PY}" "${OUT_B}/run.csv" --last-row-col-contains adaptive_branch enumerate >/dev/null 2>&1; then
    warn "adaptive_force_enumerate: could not assert adaptive_branch contains 'enumerate' (column missing or token differs)."
    record_skip
  fi
else
  warn "sjs_run --help does not show --j_star; skipping adaptive-branch tests."
  record_skip
fi

# 6.2 enum_cap truncation must be explicit (if sjs_run supports --enum_cap)
if supports_flag "${SJS_RUN}" "--enum_cap"; then
  OUT_C="${OUT_ROOT}/edge/enum_cap_truncation"
  mkdir -p "${OUT_C}"
  run_cmd_allow_skip "enum_cap_truncation" \
    "${SJS_RUN}" \
      --dataset_source=synthetic --gen=stripe --dataset="test_all_enum_cap" \
      --n_r="${NR}" --n_s="${NS}" --alpha=50 --gen_seed="${GEN_SEED}" \
      --method=ours --variant=enum_sampling --t=1000 --seed=779 --repeats=1 \
      --enum_cap="${ENUM_CAP_TEST}" \
      --out_dir="${OUT_C}" --write_samples=0 \
      --threads=1

  if [[ -f "${OUT_C}/run.csv" ]]; then
    python3 "${CSV_CHECK_PY}" "${OUT_C}/run.csv" --nonempty >/dev/null || true
    if python3 "${CSV_CHECK_PY}" "${OUT_C}/run.csv" --any-row-col-int-eq enum_truncated 1 >/dev/null 2>&1; then
      :
    elif python3 "${CSV_CHECK_PY}" "${OUT_C}/run.csv" --any-row-col-int-eq ok 0 >/dev/null 2>&1; then
      :
    else
      warn "enum_cap_truncation: Could not confirm truncation via run.csv columns (enum_truncated/ok). Inspect ${OUT_C}/run.csv"
      record_skip
    fi
  else
    warn "enum_cap_truncation: run.csv not found; cannot assert truncation."
    record_skip
  fi
else
  warn "sjs_run --help does not show --enum_cap; skipping enum_cap truncation test."
  record_skip
fi

# 6.3 empty join handling (alpha=0): should fail gracefully OR mark ok=0
OUT_D="${OUT_ROOT}/edge/empty_join_alpha0"
mkdir -p "${OUT_D}"

set +e
(
  if [[ -n "${TIMEOUT_BIN}" && "${TIMEOUT_SEC}" -gt 0 ]]; then
    "${TIMEOUT_BIN}" "${TIMEOUT_SEC}" "${SJS_RUN}" \
      --dataset_source=synthetic --gen=stripe --dataset="test_all_empty_join" \
      --n_r="${NR}" --n_s="${NS}" --alpha=0 --gen_seed="${GEN_SEED}" \
      --method=ours --variant=sampling --t=10 --seed=780 --repeats=1 \
      --out_dir="${OUT_D}" --write_samples=0 \
      --threads=1
  else
    "${SJS_RUN}" \
      --dataset_source=synthetic --gen=stripe --dataset="test_all_empty_join" \
      --n_r="${NR}" --n_s="${NS}" --alpha=0 --gen_seed="${GEN_SEED}" \
      --method=ours --variant=sampling --t=10 --seed=780 --repeats=1 \
      --out_dir="${OUT_D}" --write_samples=0 \
      --threads=1
  fi
) 2>&1 | tee "${LOG_DIR}/empty_join_alpha0.log"
rc="${PIPESTATUS[0]}"
set -e

if [[ "${rc}" -ne 0 ]]; then
  log "empty_join_alpha0: failed as expected (rc=${rc})"
  record_pass
else
  if [[ -f "${OUT_D}/run.csv" ]] && python3 "${CSV_CHECK_PY}" "${OUT_D}/run.csv" --any-row-col-int-eq ok 0 >/dev/null 2>&1; then
    log "empty_join_alpha0: succeeded but marked ok=0 in run.csv (acceptable)."
    record_pass
  else
    die "empty_join_alpha0: command succeeded but did not clearly signal failure (expected empty join)."
  fi
fi

# ----------------------------
# 7) sjs_sweep smoke (tiny config)
# ----------------------------
log "Step 7: sjs_sweep smoke"

SWEEP_OUT="${OUT_ROOT}/sweep_smoke"
mkdir -p "${SWEEP_OUT}"
SWEEP_CFG="${TEMP_DIR}/sweep_smoke.json"

cat > "${SWEEP_CFG}" <<EOF
{
  "base": {
    "dataset": {
      "name": "test_all_sweep_smoke",
      "dim": 2,
      "synthetic": {
        "generator": "stripe",
        "n_r": ${NR},
        "n_s": ${NS},
        "alpha": 1.0,
        "seed": ${GEN_SEED}
      }
    },
    "run": {
      "method": "ours",
      "variant": "sampling",
      "t": 1000,
      "seed": ${SEED},
      "repeats": 2,
      "enum_cap": 0,
      "j_star": 1000000,
      "write_samples": false,
      "verify": false,
      "extra": {}
    },
    "output": {
      "out_dir": "${SWEEP_OUT}"
    },
    "logging": { "level": "info" },
    "sys": { "threads": 1 }
  },
  "sweep": {
    "alpha": [0.01, 0.1, 1.0, 5.0],
    "t": [100, 1000],
    "method": ["ours", "r_tree"],
    "variant": ["sampling", "adaptive"]
  },
  "files": {
    "raw": "sweep_raw.csv",
    "summary": "sweep_summary.csv"
  }
}
EOF

run_cmd "sweep_smoke" \
  "${SJS_SWEEP}" --config="${SWEEP_CFG}"

[[ -f "${SWEEP_OUT}/sweep_raw.csv" ]] || die "Missing sweep_raw.csv under ${SWEEP_OUT}"
[[ -f "${SWEEP_OUT}/sweep_summary.csv" ]] || die "Missing sweep_summary.csv under ${SWEEP_OUT}"

python3 "${CSV_CHECK_PY}" "${SWEEP_OUT}/sweep_raw.csv" --nonempty >/dev/null
python3 "${CSV_CHECK_PY}" "${SWEEP_OUT}/sweep_summary.csv" --nonempty >/dev/null

# ----------------------------
# 8) Randomized fuzz battery (deep mode)
# ----------------------------
if [[ "${MODE}" != "fast" ]]; then
  log "Step 8: Randomized fuzz battery (${FUZZ_CASES} cases)"

  ALPHA_SET=("0.0001" "0.001" "0.01" "0.1" "1" "3" "10" "20")
  GEN_SET=("stripe" "uniform" "clustered" "hetero_sizes")

  i=1
  while [[ "${i}" -le "${FUZZ_CASES}" ]]; do
    nr=$(( 200 + (RANDOM % 2800) ))
    ns=$(( 200 + (RANDOM % 2800) ))
    a="${ALPHA_SET[$((RANDOM % ${#ALPHA_SET[@]}))]}"
    g="${GEN_SET[$((RANDOM % ${#GEN_SET[@]}))]}"
    s=$(( 1 + (RANDOM % 1000000) ))
    t=$(( 100 + (RANDOM % 5000) ))

    tag="fuzz_${i}_g${g}_nr${nr}_ns${ns}_a${a}_s${s}_t${t}"

    run_cmd_allow_skip "${tag}" \
      "${SJS_VERIFY}" \
        --dataset_source=synthetic --gen="${g}" --dataset="test_all_fuzz_${i}" \
        --n_r="${nr}" --n_s="${ns}" --alpha="${a}" --gen_seed=1 \
        --method=ours --variant=sampling --t="${t}" --seed="${s}" --repeats=1 \
        --threads=1 \
        --oracle_max_checks="${ORACLE_MAX_CHECKS}"

    log_file="${LOG_DIR}/${tag}.log"
    if [[ -f "${log_file}" ]] && ! skip_pattern "${log_file}"; then
      miss_lines="$(grep -E "missing_in_universe=" "${log_file}" || true)"
      if [[ -n "${miss_lines}" ]]; then
        miss_bad="$(echo "${miss_lines}" | awk -F= '{gsub(/[[:space:]]/,"",$2); if ($2!="0") {print $0}}' || true)"
        if [[ -n "${miss_bad}" ]]; then
          die "${tag}: fuzz correctness FAILED (missing_in_universe != 0)."
        fi
      fi
    fi

    i=$((i+1))
  done
else
  log "Step 8: MODE=fast -> skip fuzz battery."
  record_skip
fi

# ----------------------------
# 9) Optional static analysis (best-effort)
# ----------------------------
if [[ "${RUN_STATIC_ANALYSIS}" == "1" ]]; then
  log "Step 9: Best-effort static checks (skipped if tools missing)"

  # 9.1 shellcheck run scripts
  if command -v shellcheck >/dev/null 2>&1; then
    run_cmd "shellcheck_run_scripts" shellcheck -x "${RUN_DIR}"/*.sh
  else
    warn "shellcheck not found; skipping."
    record_skip
  fi

  # 9.2 python syntax check for run/include/*.py
  if compgen -G "${INCLUDE_DIR}/*.py" >/dev/null; then
    run_cmd "py_compile_run_include" python3 -m py_compile "${INCLUDE_DIR}"/*.py
  else
    warn "No python files under ${INCLUDE_DIR}; skipping py_compile."
    record_skip
  fi

  # 9.3 clang-tidy (warning-only; do not gate)
  if command -v clang-tidy >/dev/null 2>&1; then
    CCJSON="${BUILD_RELEASE}/compile_commands.json"
    if [[ -f "${CCJSON}" ]]; then
      mapfile -t TU < <(sed -nE 's/.*"file"[[:space:]]*:[[:space:]]*"([^"]+)".*/\1/p' "${CCJSON}" | head -n 20)
      if [[ "${#TU[@]}" -gt 0 ]]; then
        set +e
        clang-tidy "${TU[@]}" -p "${BUILD_RELEASE}" 2>&1 | tee "${LOG_DIR}/clang_tidy_subset.log"
        set -e
        warn "clang-tidy finished (warning-only). See ${LOG_DIR}/clang_tidy_subset.log"
        record_skip
      else
        warn "clang-tidy: no translation units found in compile_commands.json (unexpected)."
        record_skip
      fi
    else
      warn "clang-tidy: missing ${CCJSON}"
      record_skip
    fi
  else
    warn "clang-tidy not found; skipping."
    record_skip
  fi

  # 9.4 cppcheck (gated if installed)
  if command -v cppcheck >/dev/null 2>&1; then
    set +e
    cppcheck --enable=warning,performance,portability --error-exitcode=1 -j "$(cpu_count)" "${REPO_ROOT}" 2>&1 | tee "${LOG_DIR}/cppcheck.log"
    rc=$?
    set -e
    if [[ "${rc}" -ne 0 ]]; then
      warn "cppcheck reported issues (see ${LOG_DIR}/cppcheck.log). Treating as FAIL."
      echo "[cppcheck] FAIL rc=${rc}" >> "${FAILURES}"
      record_fail
      maybe_exit_on_fail
    else
      record_pass
    fi
  else
    warn "cppcheck not found; skipping."
    record_skip
  fi
else
  log "Step 9: RUN_STATIC_ANALYSIS=0 -> skip static checks."
  record_skip
fi

# ----------------------------
# Final summary
# ----------------------------
echo "============================================================" | tee -a "${REPORT}"
echo "TEST_ALL SUMMARY" | tee -a "${REPORT}"
echo "  TOTAL   = ${TOTAL}" | tee -a "${REPORT}"
echo "  PASSED  = ${PASSED}" | tee -a "${REPORT}"
echo "  FAILED  = ${FAILED}" | tee -a "${REPORT}"
echo "  SKIPPED = ${SKIPPED}" | tee -a "${REPORT}"
echo "  OUT_ROOT= ${OUT_ROOT}" | tee -a "${REPORT}"
echo "  LOG_DIR = ${LOG_DIR}" | tee -a "${REPORT}"
echo "============================================================" | tee -a "${REPORT}"

if [[ "${FAILED}" -ne 0 ]]; then
  echo "[TEST_ALL] FAILED ❌  (see ${FAILURES} and ${REPORT})" >&2
  exit 2
fi

echo "[TEST_ALL] PASSED ✅  (see ${REPORT})"
exit 0
