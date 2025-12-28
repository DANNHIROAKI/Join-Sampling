#!/usr/bin/env bash
# run/run_exp2.sh  (paper+extended runner)
#
# EXP-2: Runtime vs t (RQ2)
#
# This is a *paper-ready* runner that:
#   - keeps the original EXP-2 behavior (sweep by JSON config)
#   - fixes the most common "paper figure" issues (cold start bias, inconsistent baseline t0)
#   - optionally runs an *extended t-range* to make asymptotic trends visible
#   - always preserves transparency (raw + summary kept; full plots optional)
#
# Key upgrades vs the previous runner:
#   1) Default Δruntime baseline t0=1000 (matches EXP-2 writeup).
#   2) Two profiles (paper / extended) via --t_profile={paper|extended|both}.
#   3) Optional patch of the sweep t-list from the runner (no manual JSON edits).
#   4) Plotter is called with --min_repeats so paper plots only use fully-successful repeats.
#   5) Optional paper method set control (topk OR explicit list) without hiding raw data.
#
# NOTE (integrity):
#   - This script does NOT delete "unfavorable" raw results.
#   - It only changes: (i) t-range coverage, (ii) plotting defaults, (iii) reproducibility metadata.
#   - If you choose to show a method subset in the main paper figure, keep the full outputs
#     (mode=full) as supplementary material.
#
# Outputs (fixed):
#   - Temp (always overwritten):      <repo_root>/run/temp/exp2/
#   - Final (overwritten on success): <repo_root>/results/raw/exp2/
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

Core options:
  --config <path>                 Sweep JSON (default: config/sweeps/sweep_t.json)
  --t_profile <paper|extended|both|full>  Which t-range(s) to run (default: full)

  --t_list_paper <csv>            Override t-list for paper profile
  --t_list_ext <csv>              Override t-list for extended profile
  --t_list_full <csv>             Override t-list for full profile (default: union(paper,ext))

Build options:
  --build_type <type>             Release|Debug|RelWithDebInfo|MinSizeRel (default: Release)
  --clean                         Remove build/<type> before building
  --no-build                      Skip build step

Run options:
  --threads <int>                 Override base.sys.threads (default: 1)
  --write_samples <0|1>           Override base.run.write_samples (default: 0)
  --j_star <u64>                  Override base.run.j_star (default: 10000)

Repeats:
  --repeats_paper <int>           Effective repeats for paper profile (default: 10)
  --warmup_paper <int>            Extra warmup repeats excluded from plots (default: 2)
  --repeats_ext <int>             Effective repeats for extended profile (default: 5)
  --warmup_ext <int>              Extra warmup repeats excluded from plots (default: 1)
  --repeats_full <int>            Effective repeats for full profile (default: 5)
  --warmup_full <int>             Extra warmup repeats excluded from plots (default: 1)

Plot options:
  --no-plot                       Skip plotting
  --plot_mode <paper|full>        Default: paper
  --plot_t0 <int>                 Baseline t0 for Δruntime (default: 1000)
  --plot_error <auto|p95|stdev>   Default: auto
  --paper_with_sample_phase <0|1> Also output sample-phase paper figure (default: 1)

Paper figure method selection:
  --topk <int>                    Paper plots: top-k fastest methods + always_include (default: 6)
  --always_include <csv>          Default: ours
  --paper_methods <csv>           If set, use exactly this comma-list (plus always_include)
  --exclude_methods <csv>         Methods to exclude from plots (e.g., rejection)

Misc:
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
  if [[ -x "${build_dir}/${name}" ]]; then echo "${build_dir}/${name}"; return; fi
  if [[ -x "${build_dir}/apps/${name}" ]]; then echo "${build_dir}/apps/${name}"; return; fi
  local p
  p="$(find "${build_dir}" -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  [[ -n "${p}" ]] || die "Cannot find executable '${name}' under ${build_dir}. Did the build succeed?"
  echo "${p}"
}

to_bool_json() { [[ "$1" == "1" ]] && echo "true" || echo "false"; }

trim() { local s="$*"; s="${s#"${s%%[![:space:]]*}"}"; s="${s%"${s##*[![:space:]]}"}"; printf "%s" "${s}"; }

csv_ints_to_json_array() {
  local csv="$(trim "$1")"
  [[ -n "${csv}" ]] || die "Empty CSV list"
  local out="["
  local first=1
  IFS=',' read -ra parts <<< "${csv}"
  for p in "${parts[@]}"; do
    p="$(trim "${p}")"
    [[ -z "${p}" ]] && continue
    [[ "${p}" =~ ^[0-9]+$ ]] || die "Non-integer value in list: '${p}'"
    if [[ "${first}" -eq 0 ]]; then out+=", "; fi
    out+="${p}"
    first=0
  done
  out+="]"
  echo "${out}"
}

# Merge two CSV integer lists into a sorted, unique CSV list.
merge_t_lists_csv() {
  local a="$(trim "$1")"
  local b="$(trim "$2")"
  local merged
  merged="$(
    printf "%s,%s\n" "${a}" "${b}" \
      | tr ',' '\n' \
      | awk 'NF {gsub(/^[ \t]+|[ \t]+$/, "", $0); if($0 ~ /^[0-9]+$/) print $0}' \
      | sort -n \
      | awk '!seen[$0]++ {out = out (out==""? $0 : "," $0)} END{print out}'
  )"
  [[ -n "${merged}" ]] || die "Failed to merge t-lists into a non-empty CSV."
  echo "${merged}"
}

awk_col_idx() {
  local file="$1"
  local col="$2"
  awk -F',' -v c="${col}" 'NR==1{for(i=1;i<=NF;i++){if($i==c){print i; exit}}}' "${file}"
}

# Extract the effective sweep t-list from a sweep config (best effort across common schemas).
extract_t_list_csv() {
  local cfg="$1"
  jq -r '
    def first_nonnull(a):
      reduce a[] as $x (null; . // $x);
    def axes_t_values:
      if (.sweep.axes? != null) then
        ([.sweep.axes[]? | select((.name?=="t") or (.key?=="t") or (.param?=="t")) | .values?] | .[0])
      else null end;

    def axis_values:
      if ((.sweep.axis? // "") == "t") then .sweep.values? else null end;

    def candidate:
      first_nonnull([
        .sweep.t_list?,
        .sweep.ts?,
        .sweep.t_values?,
        .sweep.t?,
        (.sweep.params.t?),
        axis_values,
        axes_t_values
      ]);

    (candidate // empty) | if type=="array" then map(tostring)|join(",") else empty end
  ' "${cfg}"
}

# Patch the sweep config:
#   - repeats, threads, write_samples, j_star
#   - t_list (best effort across common schemas)
patch_config() {
  local cfg_in="$1"
  local cfg_out="$2"
  local profile="$3"
  local repeats_total="$4"
  local repeats_eff="$5"
  local warmup_reps="$6"
  local threads="$7"
  local write_samples_bool="$8"
  local j_star="$9"
  local t_list_json="${10}"
  local ts="${11}"

  jq \
    --arg profile "${profile}" \
    --argjson repeats "${repeats_total}" \
    --argjson repeats_eff "${repeats_eff}" \
    --argjson warmup_reps "${warmup_reps}" \
    --argjson threads "${threads}" \
    --argjson write_samples "${write_samples_bool}" \
    --argjson j_star "${j_star}" \
    --argjson t_list "${t_list_json}" \
    --arg ts "${ts}" \
    '
    def set_if_exists(path; val):
      if (getpath(path) != null) then setpath(path; val) else . end;

    def patch_t(val):
      .
      | set_if_exists(["sweep","t_list"]; val)
      | set_if_exists(["sweep","ts"]; val)
      | set_if_exists(["sweep","t_values"]; val)
      | set_if_exists(["sweep","t"]; val)
      | (if (.sweep.params? != null and .sweep.params.t? != null) then .sweep.params.t = val else . end)
      | (if ((.sweep.axis? // "") == "t" and .sweep.values? != null) then .sweep.values = val else . end)
      | (if (.sweep.axes? != null) then
            .sweep.axes = (.sweep.axes | map(
              if ((.name?=="t") or (.key?=="t") or (.param?=="t")) then
                (if .values? != null then .values = val else . end)
              else .
              end))
         else . end);

    .
    | .base.run.repeats = $repeats
    | .base.run.j_star = $j_star
    | .base.sys.threads = $threads
    | .base.run.write_samples = $write_samples
    | .sweep.repeats = $repeats
    | patch_t($t_list)
    | .meta.patch = {
        "timestamp": $ts,
        "patched_by": "run/run_exp2.sh",
        "profile": $profile,
        "overrides": {
          "threads": $threads,
          "write_samples": $write_samples,
          "j_star": $j_star,
          "repeats_total": $repeats,
          "repeats_effective": $repeats_eff,
          "warmup_reps": $warmup_reps,
          "t_list": $t_list
        }
      }
    ' "${cfg_in}" > "${cfg_out}"
}

# ----------------------------
# Plotter capability probe (best effort)
# ----------------------------
PLOT_HELP=""

plot_supports_flag() {
  local flag="$1"
  [[ -n "${PLOT_HELP}" ]] || return 1
  [[ "${PLOT_HELP}" == *"${flag}"* ]]
}

should_pass_plot_flag() {
  local flag="$1"
  # If we couldn't probe help, keep legacy behavior: pass the flag (better compatibility).
  if [[ -z "${PLOT_HELP}" ]]; then
    return 0
  fi
  plot_supports_flag "${flag}"
}

run_sweep_and_plot() {
  local profile="$1"
  local t_list_csv="$2"
  local repeats_eff="$3"
  local warmup_reps="$4"

  local total_repeats=$((repeats_eff + warmup_reps))
  local out_dir="${TEMP_DIR}/${profile}"
  mkdir -p "${out_dir}/logs"

  log "------------------------------"
  log "PROFILE=${profile}"
  log "t_list=${t_list_csv}"
  log "repeats_effective=${repeats_eff} warmup=${warmup_reps} total=${total_repeats}"
  log "------------------------------"

  local orig_cfg_copy="${out_dir}/sweep_original.json"
  local patched_cfg="${out_dir}/sweep_used.json"
  cp -f "${CONFIG}" "${orig_cfg_copy}"

  local t_list_json
  t_list_json="$(csv_ints_to_json_array "${t_list_csv}")"

  local write_samples_bool
  write_samples_bool="$(to_bool_json "${WRITE_SAMPLES}")"

  patch_config "${orig_cfg_copy}" "${patched_cfg}" "${profile}" \
    "${total_repeats}" "${repeats_eff}" "${warmup_reps}" \
    "${THREADS}" "${write_samples_bool}" "${J_STAR}" \
    "${t_list_json}" "${TS}"

  # Validate t-list extraction worked and includes t0 (for consistent Δruntime).
  local t_list_used
  t_list_used="$(extract_t_list_csv "${patched_cfg}")"
  [[ -n "${t_list_used}" ]] || die "Could not extract t-list from patched config (${patched_cfg}). Please update patch_t() candidates."
  log "t_list_used(extracted)=${t_list_used}"
  if ! echo ",${t_list_used}," | grep -q ",${PLOT_T0},"; then
    die "plot_t0=${PLOT_T0} is not present in the sweep t-list for profile=${profile}. Fix --plot_t0 or --t_list_*."
  fi

  # Write per-profile manifest.
  {
    echo "EXP-2 profile manifest"
    echo "timestamp=${TS}"
    echo "profile=${profile}"
    echo "repo_root=${REPO_ROOT}"
    echo "config_in=${CONFIG}"
    echo "config_used=${patched_cfg}"
    echo "threads=${THREADS}"
    echo "write_samples=${WRITE_SAMPLES}"
    echo "j_star=${J_STAR}"
    echo "t_list_requested=${t_list_csv}"
    echo "t_list_used=${t_list_used}"
    echo "repeats_effective=${repeats_eff}"
    echo "warmup_reps=${warmup_reps}"
    echo "repeats_total=${total_repeats}"
    echo "plot_mode=${PLOT_MODE}"
    echo "plot_t0=${PLOT_T0}"
    echo "plot_error=${PLOT_ERROR}"
    echo "paper_with_sample_phase=${PAPER_WITH_SAMPLE_PHASE}"
    echo "topk=${TOPK}"
    echo "always_include=${ALWAYS_INCLUDE}"
    echo "paper_methods=${PAPER_METHODS}"
    echo "exclude_methods=${EXCLUDE_METHODS}"
    echo "out_dir=${out_dir}"
  } > "${out_dir}/MANIFEST.txt"

  # Run sweep
  log "Running sjs_sweep for profile=${profile} ..."
  pushd "${REPO_ROOT}" >/dev/null
  set +e
  "${SJS_SWEEP}" \
    --config="${patched_cfg}" \
    --out_dir="${out_dir}" \
    --raw_file="sweep_raw.csv" \
    --summary_file="sweep_summary.csv" \
    2>&1 | tee "${out_dir}/logs/sjs_sweep.log"
  local rc="${PIPESTATUS[0]}"
  set -e
  popd >/dev/null
  [[ "${rc}" -eq 0 ]] || die "sjs_sweep failed for profile=${profile} (exit=${rc}). See ${out_dir}/logs/sjs_sweep.log"

  [[ -f "${out_dir}/sweep_raw.csv" ]] || die "Missing ${out_dir}/sweep_raw.csv"
  [[ -f "${out_dir}/sweep_summary.csv" ]] || die "Missing ${out_dir}/sweep_summary.csv"

  # Fast sanity: count_mean should be constant across t (fixed dataset).
  local sum_path="${out_dir}/sweep_summary.csv"
  local idx_ok idx_cnt
  idx_ok="$(awk_col_idx "${sum_path}" "ok_rate")"
  idx_cnt="$(awk_col_idx "${sum_path}" "count_mean")"
  if [[ -n "${idx_ok}" && -n "${idx_cnt}" ]]; then
    local cnt_minmax
    cnt_minmax="$(awk -F',' -v okc="${idx_ok}" -v cc="${idx_cnt}" 'NR>1 && $okc==1 {v=$cc; if(min==""||v<min)min=v; if(max==""||v>max)max=v} END{if(min!="")printf("min=%s max=%s",min,max)}' "${sum_path}")"
    log "[${profile}] count_mean over ok points: ${cnt_minmax:-<none>}"
  fi

  # Plot (non-fatal)
  if [[ "${DO_PLOT}" -eq 1 ]]; then
    if command -v python3 >/dev/null 2>&1; then
      log "Plotting profile=${profile} ..."

      # If we have help text, and it doesn't even contain required flags, skip plotting.
      if [[ -n "${PLOT_HELP}" ]]; then
        if ! plot_supports_flag "--out_dir" || ! plot_supports_flag "--t0"; then
          log "[WARN] Plot script '${PLOT_SCRIPT}' does not appear to support --out_dir/--t0 (per --help). Skipping plot for profile=${profile}."
          return 0
        fi
      fi

      local plot_args=()
      plot_args+=(--out_dir "${out_dir}")
      plot_args+=(--t0 "${PLOT_T0}")

      # Map --error auto -> p95 if plotter doesn't mention "auto" in --help.
      local plot_error_eff="${PLOT_ERROR}"
      if [[ "${PLOT_ERROR}" == "auto" && -n "${PLOT_HELP}" && "${PLOT_HELP}" != *"auto"* ]]; then
        plot_error_eff="p95"
        log "[WARN] Plot script help does not mention 'auto' for --error; using --error ${plot_error_eff}"
      fi
      if should_pass_plot_flag "--error"; then
        plot_args+=(--error "${plot_error_eff}")
      fi

      if should_pass_plot_flag "--warmup_reps"; then
        plot_args+=(--warmup_reps "${warmup_reps}")
      fi
      if should_pass_plot_flag "--mode"; then
        plot_args+=(--mode "${PLOT_MODE}")
      fi

      if should_pass_plot_flag "--topk"; then
        plot_args+=(--topk "${TOPK}")
      fi
      if should_pass_plot_flag "--always_include"; then
        plot_args+=(--always_include "${ALWAYS_INCLUDE}")
      fi

      # Only pass paper_methods/exclude_methods if non-empty (avoid "empty string means empty set").
      if [[ -n "${PAPER_METHODS}" ]]; then
        if should_pass_plot_flag "--paper_methods"; then
          plot_args+=(--paper_methods "${PAPER_METHODS}")
        else
          log "[WARN] Plot script does not support --paper_methods; ignoring PAPER_METHODS for profile=${profile}."
        fi
      fi

      if [[ -n "${EXCLUDE_METHODS}" ]]; then
        if should_pass_plot_flag "--exclude_methods"; then
          plot_args+=(--exclude_methods "${EXCLUDE_METHODS}")
        else
          log "[WARN] Plot script does not support --exclude_methods; ignoring EXCLUDE_METHODS for profile=${profile}."
        fi
      fi

      if should_pass_plot_flag "--min_repeats"; then
        plot_args+=(--min_repeats "${repeats_eff}")
      fi
      if should_pass_plot_flag "--paper_with_sample_phase"; then
        plot_args+=(--paper_with_sample_phase "${PAPER_WITH_SAMPLE_PHASE}")
      fi

      set +e
      python3 "${PLOT_SCRIPT}" "${plot_args[@]}" 2>&1 | tee "${out_dir}/logs/plot_exp2.log"
      local rc_plot="${PIPESTATUS[0]}"
      set -e

      if [[ "${rc_plot}" -ne 0 ]]; then
        log "[WARN] Plot step failed for profile=${profile} (exit=${rc_plot}). Raw/summary are kept; see ${out_dir}/logs/plot_exp2.log"
      fi
    else
      log "python3 not found; skipping plot step for profile=${profile}."
    fi
  else
    log "Plot step disabled (--no-plot)."
  fi
}

# ----------------------------
# Resolve paths
# ----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ----------------------------
# Defaults
# ----------------------------
CONFIG="config/sweeps/sweep_t.json"
T_PROFILE="full"

T_LIST_PAPER="1000,3000,10000,30000,100000,300000,1000000"
T_LIST_EXT="1000,1000000,3000000,10000000,30000000"
T_LIST_FULL=""  # If empty, will be set to union(t_list_paper, t_list_ext)

BUILD_TYPE="Release"
THREADS=1
WRITE_SAMPLES=0
J_STAR=10000

REPEATS_PAPER=10
WARMUP_PAPER=2
REPEATS_EXT=5
WARMUP_EXT=1
REPEATS_FULL=5
WARMUP_FULL=1

PLOT_MODE="paper"
PLOT_T0=1000
PLOT_ERROR="auto"
PAPER_WITH_SAMPLE_PHASE=1

TOPK=6
ALWAYS_INCLUDE="ours"
PAPER_METHODS=""
EXCLUDE_METHODS=""

DO_CLEAN=0
DO_BUILD=1
DO_PLOT=1

# ----------------------------
# Parse args
# ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)            CONFIG="$2"; shift 2;;
    --t_profile)         T_PROFILE="$2"; shift 2;;
    --t_list_paper)      T_LIST_PAPER="$2"; shift 2;;
    --t_list_ext)        T_LIST_EXT="$2"; shift 2;;
    --t_list_full)       T_LIST_FULL="$2"; shift 2;;

    --build_type)        BUILD_TYPE="$2"; shift 2;;
    --threads)           THREADS="$2"; shift 2;;
    --write_samples)     WRITE_SAMPLES="$2"; shift 2;;
    --j_star)            J_STAR="$2"; shift 2;;

    --repeats_paper)     REPEATS_PAPER="$2"; shift 2;;
    --warmup_paper)      WARMUP_PAPER="$2"; shift 2;;
    --repeats_ext)       REPEATS_EXT="$2"; shift 2;;
    --warmup_ext)        WARMUP_EXT="$2"; shift 2;;
    --repeats_full)      REPEATS_FULL="$2"; shift 2;;
    --warmup_full)       WARMUP_FULL="$2"; shift 2;;

    --plot_mode)         PLOT_MODE="$2"; shift 2;;
    --plot_t0)           PLOT_T0="$2"; shift 2;;
    --plot_error)        PLOT_ERROR="$2"; shift 2;;
    --paper_with_sample_phase) PAPER_WITH_SAMPLE_PHASE="$2"; shift 2;;

    --topk)              TOPK="$2"; shift 2;;
    --always_include)    ALWAYS_INCLUDE="$2"; shift 2;;
    --paper_methods)     PAPER_METHODS="$2"; shift 2;;
    --exclude_methods)   EXCLUDE_METHODS="$2"; shift 2;;

    --clean)             DO_CLEAN=1; shift;;
    --no-build)          DO_BUILD=0; shift;;
    --no-plot)           DO_PLOT=0; shift;;

    -h|--help)           usage; exit 0;;
    *) die "Unknown argument: $1 (try --help)";;
  esac
done

if [[ "${CONFIG}" != /* ]]; then
  CONFIG="${REPO_ROOT}/${CONFIG}"
fi

if [[ -z "${T_LIST_FULL}" ]]; then
  T_LIST_FULL="$(merge_t_lists_csv "${T_LIST_PAPER}" "${T_LIST_EXT}")"
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
if ! [[ "${J_STAR}" =~ ^[0-9]+$ ]]; then die "--j_star must be a non-negative integer"; fi

if [[ "${T_PROFILE}" != "paper" && "${T_PROFILE}" != "extended" && "${T_PROFILE}" != "both" && "${T_PROFILE}" != "full" ]]; then
  die "--t_profile must be paper|extended|both|full"
fi

for v in "${REPEATS_PAPER}" "${WARMUP_PAPER}" "${REPEATS_EXT}" "${WARMUP_EXT}" "${REPEATS_FULL}" "${WARMUP_FULL}"; do
  if ! [[ "${v}" =~ ^[0-9]+$ ]]; then die "Repeats/warmup values must be integers"; fi
done

if [[ "${PLOT_MODE}" != "paper" && "${PLOT_MODE}" != "full" ]]; then die "--plot_mode must be paper|full"; fi
if ! [[ "${PLOT_T0}" =~ ^[0-9]+$ ]] || [[ "${PLOT_T0}" -le 0 ]]; then die "--plot_t0 must be a positive integer"; fi
if [[ "${PLOT_ERROR}" != "auto" && "${PLOT_ERROR}" != "p95" && "${PLOT_ERROR}" != "stdev" ]]; then die "--plot_error must be auto|p95|stdev"; fi
if ! [[ "${TOPK}" =~ ^[0-9]+$ ]] || [[ "${TOPK}" -le 0 ]]; then die "--topk must be a positive integer"; fi
if [[ "${PAPER_WITH_SAMPLE_PHASE}" != "0" && "${PAPER_WITH_SAMPLE_PHASE}" != "1" ]]; then die "--paper_with_sample_phase must be 0 or 1"; fi

rm -rf "${TEMP_DIR}"
mkdir -p "${TEMP_DIR}/logs"

TS="$(date +%Y%m%d_%H%M%S)"

log "Repo root      : ${REPO_ROOT}"
log "Config (in)    : ${CONFIG}"
log "t_profile      : ${T_PROFILE}"
log "t_list_paper   : ${T_LIST_PAPER}"
log "t_list_ext     : ${T_LIST_EXT}"
log "t_list_full    : ${T_LIST_FULL}"
log "Build type     : ${BUILD_TYPE}"
log "Build dir      : ${BUILD_DIR}"
log "Threads        : ${THREADS}"
log "write_samples  : ${WRITE_SAMPLES}"
log "j_star         : ${J_STAR}"
log "paper repeats  : eff=${REPEATS_PAPER} warmup=${WARMUP_PAPER}"
log "ext repeats    : eff=${REPEATS_EXT} warmup=${WARMUP_EXT}"
log "full repeats   : eff=${REPEATS_FULL} warmup=${WARMUP_FULL}"
log "plot_mode      : ${PLOT_MODE}"
log "plot_t0        : ${PLOT_T0}"
log "plot_error     : ${PLOT_ERROR}"
log "paper_with_sample_phase: ${PAPER_WITH_SAMPLE_PHASE}"
log "topk           : ${TOPK}"
log "always_include : ${ALWAYS_INCLUDE}"
log "paper_methods  : ${PAPER_METHODS:-<auto/topk>}"
log "exclude_methods: ${EXCLUDE_METHODS:-<none>}"
log "Temp dir       : ${TEMP_DIR}"
log "Final results  : ${RESULT_DIR} (will be overwritten on success)"

export OMP_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export VECLIB_MAXIMUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"

{
  echo "timestamp=${TS}"
  echo "repo_root=${REPO_ROOT}"
  echo "config_in=${CONFIG}"
  echo "build_type=${BUILD_TYPE}"
  echo "build_dir=${BUILD_DIR}"
  echo "t_profile=${T_PROFILE}"
  echo "t_list_paper=${T_LIST_PAPER}"
  echo "t_list_ext=${T_LIST_EXT}"
  echo "t_list_full=${T_LIST_FULL}"
  echo "threads=${THREADS}"
  echo "write_samples=${WRITE_SAMPLES}"
  echo "j_star=${J_STAR}"
  echo "repeats_paper_effective=${REPEATS_PAPER}"
  echo "warmup_paper=${WARMUP_PAPER}"
  echo "repeats_ext_effective=${REPEATS_EXT}"
  echo "warmup_ext=${WARMUP_EXT}"
  echo "repeats_full_effective=${REPEATS_FULL}"
  echo "warmup_full=${WARMUP_FULL}"
  echo "plot_mode=${PLOT_MODE}"
  echo "plot_t0=${PLOT_T0}"
  echo "plot_error=${PLOT_ERROR}"
  echo "paper_with_sample_phase=${PAPER_WITH_SAMPLE_PHASE}"
  echo "topk=${TOPK}"
  echo "always_include=${ALWAYS_INCLUDE}"
  echo "paper_methods=${PAPER_METHODS}"
  echo "exclude_methods=${EXCLUDE_METHODS}"
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

# Plot script location (prefer run/include/exp2_plot.py)
PLOT_SCRIPT="${REPO_ROOT}/run/include/exp2_plot.py"
if [[ ! -f "${PLOT_SCRIPT}" ]]; then
  PLOT_SCRIPT="${REPO_ROOT}/run/plot_exp2.py"
fi
[[ -f "${PLOT_SCRIPT}" ]] || die "Missing plot script at run/include/exp2_plot.py or run/plot_exp2.py"

# Probe plotter help once (best effort). If it fails, we keep legacy behavior.
if [[ "${DO_PLOT}" -eq 1 ]] && command -v python3 >/dev/null 2>&1; then
  PLOT_HELP="$(python3 "${PLOT_SCRIPT}" --help 2>&1 || true)"
fi

# ----------------------------
# Run profiles
# ----------------------------
if [[ "${T_PROFILE}" == "full" ]]; then
  run_sweep_and_plot "full" "${T_LIST_FULL}" "${REPEATS_FULL}" "${WARMUP_FULL}"
else
  if [[ "${T_PROFILE}" == "paper" || "${T_PROFILE}" == "both" ]]; then
    run_sweep_and_plot "paper" "${T_LIST_PAPER}" "${REPEATS_PAPER}" "${WARMUP_PAPER}"
  fi

  if [[ "${T_PROFILE}" == "extended" || "${T_PROFILE}" == "both" ]]; then
    run_sweep_and_plot "extended" "${T_LIST_EXT}" "${REPEATS_EXT}" "${WARMUP_EXT}"
  fi
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
