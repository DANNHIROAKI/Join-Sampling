#!/usr/bin/env bash
# run/run_exp1.sh
#
# EXP-1: Correctness & Sample Quality (RQ1)
# ---------------------------------------
# This script runs EXP-1 end-to-end, from zero:
#   1) Configure + build (Release by default)
#   2) (Optional) run ctest
#   3) Run sjs_verify for: method × variant × seeds (repeats)
#   4) Parse logs into CSV tables:
#        results/exp1/exp1_quality_raw.csv
#        results/exp1/exp1_quality_summary.csv
#
# Usage (from repo root):
#   bash run/run_exp1.sh
#
# Optional overrides (env vars):
#   CLEAN_BUILD=1            # delete build dir before building
#   RUN_TESTS=0              # skip ctest
#   BUILD_TYPE=Release
#   BUILD_DIR=build
#   JOBS=8
#   OUT_BASE=results/exp1
#
# EXP1 dataset / sampling knobs:
#   GEN=stripe              # synthetic generator
#   NR=2000                 # |R|
#   NS=2000                 # |S|
#   ALPHA=1                 # |J| ~= alpha*(NR+NS) for stripe_ctrl_alpha
#   GEN_SEED=1              # synthetic dataset seed
#   T=100000                # sample size per run
#   SEED0=1                 # first run seed
#   REPEATS=5               # number of runs: SEED0..SEED0+REPEATS-1
#   ORACLE_MAX_CHECKS=50000000
#   ORACLE_COLLECT_LIMIT=1000000
#   ORACLE_CAP=0
#
# Restrict threads for fairness / reproducibility:
#   OMP_NUM_THREADS=1 etc. are set by default below.
#
set -euo pipefail

########################################
# Helpers
########################################

die() { echo "[run_exp1] ERROR: $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"
}

detect_jobs() {
  if [[ -n "${JOBS:-}" ]]; then
    echo "${JOBS}"
    return
  fi
  if command -v nproc >/dev/null 2>&1; then
    nproc
    return
  fi
  if command -v sysctl >/dev/null 2>&1; then
    sysctl -n hw.ncpu 2>/dev/null || echo 4
    return
  fi
  echo 4
}

########################################
# Resolve paths
########################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BUILD_TYPE="${BUILD_TYPE:-Release}"
BUILD_DIR="${BUILD_DIR:-${REPO_ROOT}/build}"
OUT_BASE="${OUT_BASE:-${REPO_ROOT}/results/exp1}"
LOG_DIR="${OUT_BASE}/logs"

JOBS_DETECTED="$(detect_jobs)"

########################################
# EXP-1 parameters (defaults)
########################################

GEN="${GEN:-stripe}"
NR="${NR:-2000}"
NS="${NS:-2000}"
ALPHA="${ALPHA:-1}"
GEN_SEED="${GEN_SEED:-1}"
T="${T:-100000}"

SEED0="${SEED0:-1}"
REPEATS="${REPEATS:-5}"

ORACLE_MAX_CHECKS="${ORACLE_MAX_CHECKS:-50000000}"
ORACLE_COLLECT_LIMIT="${ORACLE_COLLECT_LIMIT:-1000000}"
ORACLE_CAP="${ORACLE_CAP:-0}"

RUN_TESTS="${RUN_TESTS:-1}"
CLEAN_BUILD="${CLEAN_BUILD:-0}"

# Methods / variants in this Dim=2 build (override if you want a subset)
METHODS_DEFAULT=("ours" "aabb" "interval_tree" "kd_tree" "r_tree" "range_tree" "pbsm" "tlsop" "sirs" "rejection" "tsunami")
VARIANTS_DEFAULT=("sampling" "enum_sampling" "adaptive")

# Allow overriding METHODS / VARIANTS as space-separated strings
if [[ -n "${METHODS:-}" ]]; then
  read -r -a METHODS_ARR <<< "${METHODS}"
else
  METHODS_ARR=("${METHODS_DEFAULT[@]}")
fi

if [[ -n "${VARIANTS:-}" ]]; then
  read -r -a VARIANTS_ARR <<< "${VARIANTS}"
else
  VARIANTS_ARR=("${VARIANTS_DEFAULT[@]}")
fi

########################################
# Preflight
########################################

need_cmd cmake
need_cmd python3

mkdir -p "${OUT_BASE}" "${LOG_DIR}"

# Thread caps (fairness + reproducibility)
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

########################################
# Build (from scratch)
########################################

echo "[run_exp1] Repo root: ${REPO_ROOT}"
echo "[run_exp1] Build dir: ${BUILD_DIR} (type=${BUILD_TYPE}, jobs=${JOBS_DETECTED})"
echo "[run_exp1] Output dir: ${OUT_BASE}"

if [[ "${CLEAN_BUILD}" == "1" ]]; then
  echo "[run_exp1] CLEAN_BUILD=1 -> rm -rf ${BUILD_DIR}"
  rm -rf "${BUILD_DIR}"
fi

cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
cmake --build "${BUILD_DIR}" -j "${JOBS_DETECTED}"

if [[ "${RUN_TESTS}" == "1" ]]; then
  echo "[run_exp1] Running tests (ctest)..."
  ctest --test-dir "${BUILD_DIR}" --output-on-failure
else
  echo "[run_exp1] RUN_TESTS=0 -> skip ctest"
fi

SJS_VERIFY="${SJS_VERIFY:-${BUILD_DIR}/sjs_verify}"
[[ -x "${SJS_VERIFY}" ]] || die "sjs_verify not found or not executable at: ${SJS_VERIFY}"

########################################
# Manifest (what we ran)
########################################

DATASET_NAME="exp1_${GEN}_nr${NR}_ns${NS}_a${ALPHA}_g${GEN_SEED}"
MANIFEST="${OUT_BASE}/MANIFEST.txt"

{
  echo "EXP-1 manifest"
  echo "-------------"
  echo "repo_root=${REPO_ROOT}"
  echo "build_dir=${BUILD_DIR}"
  echo "build_type=${BUILD_TYPE}"
  echo "jobs=${JOBS_DETECTED}"
  echo ""
  echo "dataset_source=synthetic"
  echo "gen=${GEN}"
  echo "dataset=${DATASET_NAME}"
  echo "n_r=${NR}"
  echo "n_s=${NS}"
  echo "alpha=${ALPHA}"
  echo "gen_seed=${GEN_SEED}"
  echo ""
  echo "t=${T}"
  echo "seed0=${SEED0}"
  echo "repeats=${REPEATS}"
  echo ""
  echo "oracle_max_checks=${ORACLE_MAX_CHECKS}"
  echo "oracle_collect_limit=${ORACLE_COLLECT_LIMIT}"
  echo "oracle_cap=${ORACLE_CAP}"
  echo ""
  echo "methods=${METHODS_ARR[*]}"
  echo "variants=${VARIANTS_ARR[*]}"
} > "${MANIFEST}"

########################################
# Run EXP-1 (method × variant)
########################################

FAIL_FILE="${OUT_BASE}/FAILURES.txt"
: > "${FAIL_FILE}"

echo "[run_exp1] Dataset label: ${DATASET_NAME}"
echo "[run_exp1] NR=${NR} NS=${NS} ALPHA=${ALPHA} T=${T} SEED0=${SEED0} REPEATS=${REPEATS}"
echo "[run_exp1] oracle_max_checks=${ORACLE_MAX_CHECKS} oracle_collect_limit=${ORACLE_COLLECT_LIMIT} oracle_cap=${ORACLE_CAP}"
echo "[run_exp1] Manifest written to: ${MANIFEST}"

total_jobs=0
failed_jobs=0

for m in "${METHODS_ARR[@]}"; do
  for v in "${VARIANTS_ARR[@]}"; do
    total_jobs=$((total_jobs+1))
    log="${LOG_DIR}/${m}__${v}.log"

    echo "[run_exp1] (${total_jobs}) method=${m} variant=${v}"
    echo "           log=${log}"

    # Run and capture log. Keep going even if one combo fails.
    if ! "${SJS_VERIFY}" \
      --dataset_source=synthetic --gen="${GEN}" --dataset="${DATASET_NAME}" \
      --n_r="${NR}" --n_s="${NS}" --alpha="${ALPHA}" --gen_seed="${GEN_SEED}" \
      --method="${m}" --variant="${v}" --t="${T}" --seed="${SEED0}" --repeats="${REPEATS}" \
      --oracle_max_checks="${ORACLE_MAX_CHECKS}" --oracle_collect_limit="${ORACLE_COLLECT_LIMIT}" --oracle_cap="${ORACLE_CAP}" \
      > "${log}" 2>&1
    then
      rc=$?
      failed_jobs=$((failed_jobs+1))
      echo "[run_exp1] FAIL rc=${rc} method=${m} variant=${v} log=${log}" | tee -a "${FAIL_FILE}"
      continue
    fi
  done
done

########################################
# Parse logs -> CSV tables
########################################

echo "[run_exp1] Parsing logs -> CSV ..."

python3 - <<'PY' "${LOG_DIR}" "${OUT_BASE}"
import csv, glob, os, re, statistics, sys

log_dir = sys.argv[1]
out_base = sys.argv[2]
fail_path = os.path.join(out_base, "FAILURES.txt")

logs = sorted(glob.glob(os.path.join(log_dir, "*.log")))

# Regex patterns that match sjs_verify output
re_run = re.compile(r"^---- run rep=(\d+)\s+seed=(\d+)\s+----$")
re_method = re.compile(r"^method=([^\s]+)\s+variant=([^\s]+)\s+t=(\d+)\s*$")
re_count = re.compile(r"^count=([0-9.eE+-]+)\s+\((exact|est)\)\s+oracle=([0-9.eE+-]+)\s+rel_err=([0-9.eE+-]+)\s*$")
re_samples = re.compile(r"^samples=(\d+)\s*$")
re_failed = re.compile(r"^FAILED:\s*(.*)$")
re_quality_skipped = re.compile(r"^quality:\s*skipped\b")
re_missing = re.compile(r"^\s+missing_in_universe=([0-9.eE+-]+)\s*$")
re_chi2 = re.compile(r"^\s+chi2_stat=.*\s+p_value=([0-9.eE+-]+)\s*$")
re_ac1 = re.compile(r"^\s+autocorr_hash_lag1=([0-9.eE+-]+)\s*$")
re_ks = re.compile(r"^\s+ks_hash_uniform01\s+D=.*\s+p=([0-9.eE+-]+)\s*$")

rows = []

def push_row(cur):
    if not cur:
        return
    cur.setdefault("status", "OK")
    rows.append(cur)

# Parse each log file
for path in logs:
    base = os.path.basename(path)
    # filename: <method>__<variant>.log
    m_guess, v_guess = "", ""
    if "__" in base:
        m_guess = base.split("__", 1)[0]
        v_guess = base.split("__", 1)[1].rsplit(".", 1)[0]

    cur = None
    saw_any_run = False

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")

            m = re_run.match(line)
            if m:
                saw_any_run = True
                # flush previous run record
                push_row(cur)
                cur = {
                    "method": m_guess,
                    "variant": v_guess,
                    "rep": int(m.group(1)),
                    "seed": int(m.group(2)),
                    "log_file": base,
                    "status": "OK",
                }
                continue

            if cur is None:
                continue  # ignore prologue before first run

            m = re_failed.match(line)
            if m:
                cur["status"] = "FAILED"
                cur["error"] = m.group(1).strip()
                continue

            m = re_method.match(line)
            if m:
                cur["method"] = m.group(1)
                cur["variant"] = m.group(2)
                cur["t"] = int(m.group(3))
                continue

            m = re_count.match(line)
            if m:
                cur["count"] = float(m.group(1))
                cur["count_kind"] = m.group(2)
                cur["oracle"] = float(m.group(3))
                cur["rel_err"] = float(m.group(4))
                continue

            m = re_samples.match(line)
            if m:
                cur["samples"] = int(m.group(1))
                continue

            if re_quality_skipped.match(line):
                cur["status"] = "QUALITY_SKIPPED"
                continue

            m = re_missing.match(line)
            if m:
                cur["missing_in_universe"] = float(m.group(1))
                continue

            m = re_chi2.match(line)
            if m:
                cur["chi2_p"] = float(m.group(1))
                continue

            m = re_ac1.match(line)
            if m:
                cur["autocorr_lag1"] = float(m.group(1))
                continue

            m = re_ks.match(line)
            if m:
                cur["ks_p"] = float(m.group(1))
                continue

    push_row(cur)

    # If a log has no "---- run ..." blocks at all, keep a placeholder row.
    if not saw_any_run:
        rows.append({
            "method": m_guess,
            "variant": v_guess,
            "rep": "",
            "seed": "",
            "log_file": base,
            "status": "NO_RUN_BLOCKS",
        })

# Incorporate top-level failures (non-zero exit codes) so they appear in the CSV too.
if os.path.exists(fail_path):
    with open(fail_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            m = re.search(r"rc=(\d+)\s+method=([^\s]+)\s+variant=([^\s]+)\s+log=(.+)$", line)
            if not m:
                continue
            rc = int(m.group(1))
            method = m.group(2)
            variant = m.group(3)
            log = os.path.basename(m.group(4))
            rows.append({
                "method": method,
                "variant": variant,
                "rep": "",
                "seed": "",
                "log_file": log,
                "status": "FAILED_CREATE_OR_ARGS",
                "exit_code": rc,
            })

# Write raw table
raw_csv = os.path.join(out_base, "exp1_quality_raw.csv")
os.makedirs(out_base, exist_ok=True)

all_fields = sorted({k for r in rows for k in r.keys()})
with open(raw_csv, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=all_fields)
    w.writeheader()
    for r in rows:
        w.writerow(r)

# Summary (median by method×variant over OK + QUALITY_SKIPPED)
def median(vals):
    vals = [v for v in vals if v is not None and v != "" and (not isinstance(v, float) or v == v)]
    return statistics.median(vals) if vals else ""

summary_csv = os.path.join(out_base, "exp1_quality_summary.csv")
keys = sorted(set((r.get("method",""), r.get("variant","")) for r in rows))

with open(summary_csv, "w", newline="", encoding="utf-8") as f:
    fields = [
        "method","variant","n_rows",
        "n_ok","n_failed","n_failed_create_or_args","n_quality_skipped","n_no_run_blocks",
        "t_median","oracle_median","rel_err_median",
        "missing_median","chi2_p_median","ks_p_median","autocorr_median"
    ]
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()

    for method, variant in keys:
        grp = [r for r in rows if r.get("method")==method and r.get("variant")==variant]
        n_rows = len(grp)
        n_ok = sum(1 for r in grp if r.get("status")=="OK")
        n_failed = sum(1 for r in grp if r.get("status")=="FAILED")
        n_fc = sum(1 for r in grp if r.get("status")=="FAILED_CREATE_OR_ARGS")
        n_qs = sum(1 for r in grp if r.get("status")=="QUALITY_SKIPPED")
        n_nrb = sum(1 for r in grp if r.get("status")=="NO_RUN_BLOCKS")

        w.writerow({
            "method": method,
            "variant": variant,
            "n_rows": n_rows,
            "n_ok": n_ok,
            "n_failed": n_failed,
            "n_failed_create_or_args": n_fc,
            "n_quality_skipped": n_qs,
            "n_no_run_blocks": n_nrb,
            "t_median": median([r.get("t") for r in grp]),
            "oracle_median": median([r.get("oracle") for r in grp]),
            "rel_err_median": median([r.get("rel_err") for r in grp]),
            "missing_median": median([r.get("missing_in_universe") for r in grp]),
            "chi2_p_median": median([r.get("chi2_p") for r in grp]),
            "ks_p_median": median([r.get("ks_p") for r in grp]),
            "autocorr_median": median([r.get("autocorr_lag1") for r in grp]),
        })

print("[run_exp1] Wrote:", raw_csv)
print("[run_exp1] Wrote:", summary_csv)

# Quick sanity: missing should be 0 when quality is computed.
bad_missing = [r for r in rows if r.get("status")=="OK" and float(r.get("missing_in_universe", 0.0) or 0.0) != 0.0]
if bad_missing:
    print("[run_exp1] WARNING: missing_in_universe != 0 found in", len(bad_missing), "OK runs (correctness failure).")
PY

########################################
# Final report
########################################

echo "[run_exp1] Done."
echo "[run_exp1] Logs:   ${LOG_DIR}/*.log"
echo "[run_exp1] Tables: ${OUT_BASE}/exp1_quality_raw.csv"
echo "               ${OUT_BASE}/exp1_quality_summary.csv"
echo "[run_exp1] Manifest: ${MANIFEST}"

if [[ "${failed_jobs}" -gt 0 ]]; then
  echo "[run_exp1] WARNING: ${failed_jobs}/${total_jobs} (method×variant) jobs failed to run. See ${FAIL_FILE}" >&2
  exit 2
fi

exit 0
