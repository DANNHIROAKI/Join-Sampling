#!/usr/bin/env bash
# run.sh (repo root)
#
# HighDims build + experiment launcher.
#
# Typical usage:
#   chmod +x run.sh
#   ./run.sh build
#   ./run.sh test
#   ./run.sh sweep all
#   ./run.sh sweep dim
#   ./run.sh sweep config/sweeps/sweep_alpha.json
#
# You can also forward args:
#   ./run.sh run --help
#   ./run.sh gen --help
#   ./run.sh verify --help
#
# Environment variables:
#   BUILD_TYPE=Release|Debug|RelWithDebInfo|MinSizeRel   (default: Release)
#   CLEAN_BUILD=0|1                                     (default: 0)
#   JOBS=8                                              (default: auto)
#   THREADS=1                                           (default: 1)   # passed to apps if not overridden
#
# Note:
#  - Sweep output dirs are controlled by config/sweeps/*.json ("base.output.out_dir").
#  - This script runs from repo root so relative paths work as expected.

set -Eeuo pipefail
IFS=$' \t\n'

trap 'echo -e "[run.sh][FATAL] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

log()  { echo -e "[run.sh] $*"; }
warn() { echo -e "[run.sh][WARN] $*" >&2; }

usage() {
  cat <<'EOF'
Usage:
  ./run.sh build
  ./run.sh test
  ./run.sh sweep <all|alpha|scale|t|dim|path/to/sweep.json> [extra args...]
  ./run.sh run    [args passed to sjs_run...]
  ./run.sh gen    [args passed to sjs_gen_dataset...]
  ./run.sh verify [args passed to sjs_verify...]
  ./run.sh rectgen [args passed to tools/alacarte_rectgen_generate.py...]

Examples:
  ./run.sh
  ./run.sh sweep all
  ./run.sh sweep dim --threads=8
  ./run.sh sweep config/sweeps/sweep_alpha.json --out_dir=results/sweeps/alpha_override

  # Generate a dataset with alacarte-rectgen (Python) directly:
  ./run.sh rectgen --nR=100000 --nS=100000 --d=4 --alpha_out=10 --seed=1 \
    --out_r=data/synthetic/demo__d4_R.bin --out_s=data/synthetic/demo__d4_S.bin \
    --dataset_name=demo --report_path=data/synthetic/demo__d4_gen_report.json

Environment:
  BUILD_TYPE=Release|Debug|RelWithDebInfo|MinSizeRel
  CLEAN_BUILD=0|1
  JOBS=<n>
  THREADS=<n>
EOF
}

nproc_safe() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
  elif command -v sysctl >/dev/null 2>&1; then
    sysctl -n hw.ncpu 2>/dev/null || echo 4
  elif command -v getconf >/dev/null 2>&1; then
    getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4
  else
    echo 4
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

has_arg_prefix() {
  # has_arg_prefix "--threads" "${@:2}"
  local prefix="$1"; shift
  for a in "$@"; do
    if [[ "$a" == "${prefix}" || "$a" == "${prefix}="* ]]; then
      return 0
    fi
  done
  return 1
}

# ----------------------------
# Config / build paths
# ----------------------------
BUILD_TYPE="${BUILD_TYPE:-Release}"
CLEAN_BUILD="${CLEAN_BUILD:-0}"
THREADS="${THREADS:-1}"

if [[ -n "${JOBS:-}" ]]; then
  JOBS="${JOBS}"
else
  JOBS="$(nproc_safe)"
fi

BUILD_SUBDIR="$(build_subdir_for_type "$BUILD_TYPE")"
BUILD_DIR="${ROOT}/build/${BUILD_SUBDIR}"

cd "$ROOT"

do_build() {
  log "Repo root : $ROOT"
  log "Build type: $BUILD_TYPE"
  log "Build dir : $BUILD_DIR"
  log "JOBS      : $JOBS"

  if [[ "$CLEAN_BUILD" == "1" ]]; then
    warn "CLEAN_BUILD=1: removing build dir: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
  fi

  if [[ ! -f "${BUILD_DIR}/CMakeCache.txt" ]]; then
    log "Configuring CMake..."
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
  fi

  log "Building..."
  cmake --build "$BUILD_DIR" -j "$JOBS"
}

resolve_sweep_path() {
  local name="$1"
  case "$name" in
    all|alpha|scale|t|dim)
      echo "${ROOT}/config/sweeps/sweep_${name}.json"
      ;;
    *)
      # treat as a path
      if [[ "$name" = /* ]]; then
        echo "$name"
      else
        echo "${ROOT}/${name}"
      fi
      ;;
  esac
}

run_sweep_one() {
  local sweep_json="$1"; shift
  local sjs_sweep_exe="$1"; shift

  if [[ ! -f "$sweep_json" ]]; then
    warn "Sweep config not found: $sweep_json"
    return 2
  fi

  local extra_args=("$@")

  # If user didn't pass --threads, add a default.
  if ! has_arg_prefix "--threads" "${extra_args[@]}"; then
    extra_args+=("--threads=${THREADS}")
  fi

  log "Running sweep: $sweep_json"
  log "Cmd: ${sjs_sweep_exe} --config=${sweep_json} ${extra_args[*]}"
  "${sjs_sweep_exe}" --config="${sweep_json}" "${extra_args[@]}"
}

do_test() {
  do_build
  log "Running ctest..."
  (cd "$BUILD_DIR" && ctest --output-on-failure)
}

do_sweep() {
  do_build

  local SJS_SWEEP
  SJS_SWEEP="$(find_exe "$BUILD_DIR" sjs_sweep || true)"
  if [[ -z "${SJS_SWEEP}" ]]; then
    warn "Could not locate sjs_sweep under: $BUILD_DIR"
    return 2
  fi

  local which="$1"; shift || true

  if [[ -z "${which}" || "${which}" == "all" ]]; then
    run_sweep_one "$(resolve_sweep_path alpha)" "$SJS_SWEEP" "$@"
    run_sweep_one "$(resolve_sweep_path scale)" "$SJS_SWEEP" "$@"
    run_sweep_one "$(resolve_sweep_path t)"     "$SJS_SWEEP" "$@"
    run_sweep_one "$(resolve_sweep_path dim)"   "$SJS_SWEEP" "$@"
    return 0
  fi

  local sweep_json
  sweep_json="$(resolve_sweep_path "$which")"
  run_sweep_one "$sweep_json" "$SJS_SWEEP" "$@"
}

do_rectgen() {
  # Run the Python generator directly (no C++ build needed).
  local PY="${SJS_PYTHON:-python3}"
  local SCRIPT="${SJS_ALACARTE_RECTGEN_SCRIPT:-${ROOT}/tools/alacarte_rectgen_generate.py}"
  if [[ ! -f "${SCRIPT}" ]]; then
    warn "RectGen script not found: ${SCRIPT}"
    warn "Set SJS_ALACARTE_RECTGEN_SCRIPT or run from repo root."
    return 2
  fi
  log "Cmd: ${PY} ${SCRIPT} $*"
  "${PY}" "${SCRIPT}" "$@"
}

forward_app() {
  local exe_name="$1"; shift
  do_build
  local exe
  exe="$(find_exe "$BUILD_DIR" "$exe_name" || true)"
  if [[ -z "$exe" ]]; then
    warn "Could not locate ${exe_name} under: $BUILD_DIR"
    return 2
  fi
  log "Cmd: $exe $*"
  "$exe" "$@"
}

# ----------------------------
# Main
# ----------------------------
cmd="${1:-}"
if [[ -z "$cmd" ]]; then
  # default behavior
  do_sweep all
  exit $?
fi

shift || true
case "$cmd" in
  -h|--help|help)
    usage
    ;;
  build)
    do_build
    ;;
  test)
    do_test
    ;;
  sweep)
    do_sweep "${1:-all}" "${@:2}"
    ;;
  run)
    forward_app sjs_run "$@"
    ;;
  gen)
    forward_app sjs_gen_dataset "$@"
    ;;
  verify)
    forward_app sjs_verify "$@"
    ;;
  rectgen|alacarte)
    do_rectgen "$@"
    ;;
  *)
    warn "Unknown command: $cmd"
    usage
    exit 2
    ;;
esac
