#!/usr/bin/env bash
# run/run_exp6.sh
#
# EXP-6 (RQ6): Adaptivity effectiveness on an (alpha, t) grid.
#
# This script is "from zero": it will
#   1) configure + build (Release by default)
#   2) generate an EXP-6 sweep JSON config (self-contained)
#   3) run sjs_sweep over (alpha, t) × {sampling, enum_sampling, adaptive}
#   4) produce two figures per method:
#        - phase diagram (winner per grid point)
#        - ratio heatmap = T(adaptive) / min(T(sampling), T(enum_sampling))
#
# Usage (from repo root):
#   bash run/run_exp6.sh
#
# Optional environment overrides (examples):
#   EXP6_N=200000 EXP6_REPEATS=5 bash run/run_exp6.sh
#   EXP6_METHODS="ours,kd_tree" EXP6_ALPHAS="0.01,0.1,1,10" bash run/run_exp6.sh
#
# Notes:
#   - The sweep uses the synthetic generator "stripe" (alias of stripe_ctrl_alpha),
#     so alpha controls |J| exactly (paper-friendly).
#   - By default we run single-thread (fairness) and disable correctness verification.
#   - Enum-based modes can be very slow at large alpha; adjust EXP6_ALPHAS if needed.

set -euo pipefail
IFS=$'\n\t'

# --------------------------
# Helpers
# --------------------------

die() { echo "[EXP6][ERROR] $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"
}

cpu_count() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
  elif command -v sysctl >/dev/null 2>&1; then
    sysctl -n hw.ncpu
  else
    echo 4
  fi
}

trim_spaces() { echo "$1" | tr -d '[:space:]'; }

json_array_numbers() {
  local s
  s="$(trim_spaces "$1")"
  IFS=',' read -r -a arr <<< "$s"
  local out="["
  local first=1
  for v in "${arr[@]}"; do
    [[ -z "$v" ]] && continue
    if [[ $first -eq 0 ]]; then out+=", "; fi
    out+="$v"
    first=0
  done
  out+="]"
  echo "$out"
}

json_array_strings() {
  local s
  s="$(trim_spaces "$1")"
  IFS=',' read -r -a arr <<< "$s"
  local out="["
  local first=1
  for v in "${arr[@]}"; do
    [[ -z "$v" ]] && continue
    if [[ $first -eq 0 ]]; then out+=", "; fi
    out+="\"$v\""
    first=0
  done
  out+="]"
  echo "$out"
}

json_bool() {
  case "$(trim_spaces "$1")" in
    1|true|TRUE|True|yes|YES|y|Y) echo "true" ;;
    0|false|FALSE|False|no|NO|n|N) echo "false" ;;
    *) echo "false" ;;
  esac
}

find_exe() {
  local build_dir="$1"
  local name="$2"

  if [[ -x "${build_dir}/${name}" ]]; then
    echo "${build_dir}/${name}"
    return 0
  fi

  local f=""
  f="$(find "${build_dir}" -maxdepth 4 -type f -name "${name}" -perm -111 2>/dev/null | head -n 1 || true)"
  if [[ -n "$f" ]]; then
    echo "$f"
    return 0
  fi

  return 1
}

# --------------------------
# Locate repo root
# --------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# --------------------------
# Dependencies
# --------------------------

need_cmd cmake
need_cmd python3
need_cmd awk
need_cmd sed
need_cmd find

# --------------------------
# Parameters (override via env)
# --------------------------

# Build
EXP6_BUILD_DIR="${EXP6_BUILD_DIR:-${ROOT_DIR}/build}"
EXP6_BUILD_TYPE="${EXP6_BUILD_TYPE:-Release}"
EXP6_BUILD_JOBS="${EXP6_BUILD_JOBS:-$(cpu_count)}"
EXP6_CLEAN_BUILD="${EXP6_CLEAN_BUILD:-0}"   # 1 -> rm -rf build dir before configure
EXP6_RUN_TESTS="${EXP6_RUN_TESTS:-0}"       # 1 -> ctest after build

# Dataset (synthetic stripe)
# If EXP6_N is provided, it overrides both n_r and n_s.
EXP6_N="${EXP6_N:-}"
EXP6_N_R="${EXP6_N_R:-100000}"
EXP6_N_S="${EXP6_N_S:-100000}"
if [[ -n "${EXP6_N}" ]]; then
  EXP6_N_R="${EXP6_N}"
  EXP6_N_S="${EXP6_N}"
fi

EXP6_GEN_SEED="${EXP6_GEN_SEED:-1}"

# Generator params (match stripe_ctrl_alpha defaults unless overridden)
EXP6_CONTROL_AXIS="${EXP6_CONTROL_AXIS:-1}"
EXP6_CORE_LO="${EXP6_CORE_LO:-0.45}"
EXP6_CORE_HI="${EXP6_CORE_HI:-0.55}"
EXP6_GAP_FACTOR="${EXP6_GAP_FACTOR:-0.10}"
EXP6_DELTA_FACTOR="${EXP6_DELTA_FACTOR:-0.25}"
EXP6_SHUFFLE_STRIPS="${EXP6_SHUFFLE_STRIPS:-true}"
EXP6_SHUFFLE_R="${EXP6_SHUFFLE_R:-false}"
EXP6_SWAP_SIDES="${EXP6_SWAP_SIDES:-false}"

# Sweep grid
EXP6_ALPHAS="${EXP6_ALPHAS:-0.01,0.03,0.1,0.3,1,3,10,30}"
EXP6_TS="${EXP6_TS:-1000,3000,10000,30000,100000,300000,1000000}"

# Which methods to run (comma-separated). Default: ours only (fastest to iterate).
EXP6_METHODS="${EXP6_METHODS:-ours}"

# Variants are fixed for EXP-6 (compare three modes)
EXP6_VARIANTS="${EXP6_VARIANTS:-sampling,enum_sampling,adaptive}"

# Run control
EXP6_REPEATS="${EXP6_REPEATS:-3}"
EXP6_RUN_SEED="${EXP6_RUN_SEED:-1}"
EXP6_THREADS="${EXP6_THREADS:-1}"

# Adaptive threshold & enumeration safety cap
EXP6_J_STAR="${EXP6_J_STAR:-1000000}"
EXP6_ENUM_CAP="${EXP6_ENUM_CAP:-0}"

# Output
EXP6_OUT_DIR="${EXP6_OUT_DIR:-${ROOT_DIR}/results/sweeps/exp6_alpha_t}"
EXP6_CLEAN_RESULTS="${EXP6_CLEAN_RESULTS:-0}"  # 1 -> remove EXP6_OUT_DIR before run

# --------------------------
# Prepare dirs
# --------------------------

GEN_DIR="${ROOT_DIR}/run/_generated"
mkdir -p "${GEN_DIR}"

if [[ "${EXP6_CLEAN_RESULTS}" == "1" ]]; then
  echo "[EXP6] Cleaning results dir: ${EXP6_OUT_DIR}"
  rm -rf "${EXP6_OUT_DIR}"
fi
mkdir -p "${EXP6_OUT_DIR}"
mkdir -p "${EXP6_OUT_DIR}/logs"
mkdir -p "${EXP6_OUT_DIR}/figs"

# Record manifest (for reproducibility)
{
  echo "timestamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  echo "root_dir=${ROOT_DIR}"
  echo "build_dir=${EXP6_BUILD_DIR}"
  echo "build_type=${EXP6_BUILD_TYPE}"
  echo "methods=${EXP6_METHODS}"
  echo "variants=${EXP6_VARIANTS}"
  echo "alphas=${EXP6_ALPHAS}"
  echo "t_list=${EXP6_TS}"
  echo "n_r=${EXP6_N_R}"
  echo "n_s=${EXP6_N_S}"
  echo "gen_seed=${EXP6_GEN_SEED}"
  echo "run_seed=${EXP6_RUN_SEED}"
  echo "repeats=${EXP6_REPEATS}"
  echo "threads=${EXP6_THREADS}"
  echo "j_star=${EXP6_J_STAR}"
  echo "enum_cap=${EXP6_ENUM_CAP}"
  if command -v git >/dev/null 2>&1 && git -C "${ROOT_DIR}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "git_sha=$(git -C "${ROOT_DIR}" rev-parse HEAD)"
    echo "git_dirty=$(git -C "${ROOT_DIR}" status --porcelain | wc -l | awk '{print $1}')"
  else
    echo "git_sha=unknown"
    echo "git_dirty=unknown"
  fi
} > "${EXP6_OUT_DIR}/manifest.txt"

# --------------------------
# Build (from zero)
# --------------------------

if [[ "${EXP6_CLEAN_BUILD}" == "1" ]]; then
  echo "[EXP6] Cleaning build dir: ${EXP6_BUILD_DIR}"
  rm -rf "${EXP6_BUILD_DIR}"
fi

echo "[EXP6] Configuring (CMake, ${EXP6_BUILD_TYPE}) ..."
cmake -S "${ROOT_DIR}" -B "${EXP6_BUILD_DIR}" -DCMAKE_BUILD_TYPE="${EXP6_BUILD_TYPE}"

echo "[EXP6] Building (jobs=${EXP6_BUILD_JOBS}) ..."
cmake --build "${EXP6_BUILD_DIR}" -j "${EXP6_BUILD_JOBS}"

if [[ "${EXP6_RUN_TESTS}" == "1" ]]; then
  echo "[EXP6] Running tests (ctest) ..."
  (cd "${EXP6_BUILD_DIR}" && ctest --output-on-failure)
fi

# Find sjs_sweep
SJS_SWEEP="$(find_exe "${EXP6_BUILD_DIR}" "sjs_sweep" || true)"
[[ -n "${SJS_SWEEP}" ]] || die "Could not find built executable: sjs_sweep under ${EXP6_BUILD_DIR}"
echo "[EXP6] Using sjs_sweep: ${SJS_SWEEP}"

# --------------------------
# Generate sweep config JSON (self-contained)
# --------------------------

ALPHAS_JSON="$(json_array_numbers "${EXP6_ALPHAS}")"
TS_JSON="$(json_array_numbers "${EXP6_TS}")"
METHODS_JSON="$(json_array_strings "${EXP6_METHODS}")"
VARIANTS_JSON="$(json_array_strings "${EXP6_VARIANTS}")"

SHUFFLE_STRIPS_JSON="$(json_bool "${EXP6_SHUFFLE_STRIPS}")"
SHUFFLE_R_JSON="$(json_bool "${EXP6_SHUFFLE_R}")"
SWAP_SIDES_JSON="$(json_bool "${EXP6_SWAP_SIDES}")"

CONFIG_PATH="${GEN_DIR}/exp6_alpha_t.json"

cat > "${CONFIG_PATH}" <<EOF
{
  "base": {
    "dataset": {
      "source": "synthetic",
      "name": "exp6_stripe",
      "dim": 2,
      "synthetic": {
        "generator": "stripe",
        "n_r": ${EXP6_N_R},
        "n_s": ${EXP6_N_S},
        "alpha": 0.1,
        "seed": ${EXP6_GEN_SEED},
        "params": {
          "control_axis": ${EXP6_CONTROL_AXIS},
          "core_lo": ${EXP6_CORE_LO},
          "core_hi": ${EXP6_CORE_HI},
          "gap_factor": ${EXP6_GAP_FACTOR},
          "delta_factor": ${EXP6_DELTA_FACTOR},
          "shuffle_strips": ${SHUFFLE_STRIPS_JSON},
          "shuffle_r": ${SHUFFLE_R_JSON},
          "swap_sides": ${SWAP_SIDES_JSON}
        }
      }
    },

    "run": {
      "method": "ours",
      "variant": "sampling",
      "t": 100000,
      "seed": ${EXP6_RUN_SEED},
      "repeats": ${EXP6_REPEATS},
      "enum_cap": ${EXP6_ENUM_CAP},
      "j_star": ${EXP6_J_STAR},
      "write_samples": false,
      "verify": false,
      "extra": {}
    },

    "output": {
      "out_dir": "${EXP6_OUT_DIR}"
    },

    "logging": {
      "level": "info"
    },

    "sys": {
      "threads": ${EXP6_THREADS}
    }
  },

  "sweep": {
    "alpha": ${ALPHAS_JSON},
    "t": ${TS_JSON},
    "method": ${METHODS_JSON},
    "variant": ${VARIANTS_JSON}
  },

  "files": {
    "raw": "sweep_raw.csv",
    "summary": "sweep_summary.csv"
  }
}
EOF

# Keep a copy next to results for provenance
cp -f "${CONFIG_PATH}" "${EXP6_OUT_DIR}/exp6_alpha_t.json"

echo "[EXP6] Generated sweep config: ${CONFIG_PATH}"
echo "[EXP6] Output dir: ${EXP6_OUT_DIR}"

# --------------------------
# Run sweep
# --------------------------

LOG_FILE="${EXP6_OUT_DIR}/logs/sjs_sweep_$(date +"%Y%m%d_%H%M%S").log"
echo "[EXP6] Running sweep ..."
echo "[EXP6] Log: ${LOG_FILE}"

# Run from ROOT_DIR so relative paths (if any) resolve as expected.
(
  cd "${ROOT_DIR}"
  "${SJS_SWEEP}" --config="${CONFIG_PATH}"
) 2>&1 | tee "${LOG_FILE}"

SUMMARY_CSV="${EXP6_OUT_DIR}/sweep_summary.csv"
RAW_CSV="${EXP6_OUT_DIR}/sweep_raw.csv"

[[ -f "${SUMMARY_CSV}" ]] || die "Missing expected output: ${SUMMARY_CSV}"
[[ -f "${RAW_CSV}" ]] || die "Missing expected output: ${RAW_CSV}"

echo "[EXP6] Sweep done."
echo "  raw:     ${RAW_CSV}"
echo "  summary: ${SUMMARY_CSV}"

# --------------------------
# Plot (phase diagram + ratio heatmap)
# --------------------------

PLOT_PY="${GEN_DIR}/plot_exp6.py"
cat > "${PLOT_PY}" <<'PY'
import csv
import math
import os
import sys
from collections import defaultdict

def ffloat(x, default=math.nan):
    try:
        return float(x)
    except Exception:
        return default

def fint(x, default=None):
    try:
        return int(x)
    except Exception:
        return default

def read_summary(path):
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        return list(r)

def unique_sorted(vals, key=None):
    s = sorted(set(vals), key=key)
    return s

def write_matrix_csv(path, row_labels, col_labels, matrix):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["alpha\\t"] + list(col_labels))
        for i, a in enumerate(row_labels):
            w.writerow([a] + list(matrix[i]))

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_exp6.py <sweep_summary.csv> <out_dir_figs>", file=sys.stderr)
        sys.exit(2)

    summary_csv = sys.argv[1]
    figs_dir = sys.argv[2]
    os.makedirs(figs_dir, exist_ok=True)

    rows = read_summary(summary_csv)
    if not rows:
        print(f"Empty summary: {summary_csv}", file=sys.stderr)
        sys.exit(1)

    # Expect these columns (from sjs_sweep):
    # alpha, t, method, variant, ok_rate, wall_median_ms
    # We'll be robust if extra columns exist.
    parsed = []
    for r in rows:
        a = ffloat(r.get("alpha"))
        t = fint(r.get("t"))
        method = r.get("method", "")
        variant = r.get("variant", "")
        ok_rate = ffloat(r.get("ok_rate"), 0.0)
        wall_median_ms = ffloat(r.get("wall_median_ms"))
        parsed.append((method, a, t, variant, ok_rate, wall_median_ms))

    methods = unique_sorted([m for (m, *_rest) in parsed if m])

    # Try importing matplotlib only if needed
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        print("[plot_exp6] matplotlib not available. Will write CSV matrices only.")
        plt = None

    for method in methods:
        sub = [(a, t, v, ok, tm) for (m, a, t, v, ok, tm) in parsed if m == method]
        if not sub:
            continue

        alphas = unique_sorted([a for (a, *_rest) in sub], key=float)
        ts = unique_sorted([t for (_a, t, *_rest) in sub], key=int)

        # time[(alpha,t,variant)] = median_ms if ok_rate==1 else NaN
        time = {}
        for a, t, v, ok, tm in sub:
            if ok >= 1.0:
                time[(a, t, v)] = tm
            else:
                time[(a, t, v)] = math.nan

        # Build matrices: winner (0/1/2) and ratio
        variants = ["sampling", "enum_sampling", "adaptive"]

        winner_mat = [[-1 for _ in ts] for _ in alphas]
        ratio_mat = [[math.nan for _ in ts] for _ in alphas]

        for i, a in enumerate(alphas):
            for j, t in enumerate(ts):
                vals = []
                for v in variants:
                    vals.append(time.get((a, t, v), math.nan))

                # winner among available finite vals
                best_idx = -1
                best_val = math.inf
                for k, x in enumerate(vals):
                    if x is None or math.isnan(x):
                        continue
                    if x < best_val:
                        best_val = x
                        best_idx = k
                winner_mat[i][j] = best_idx

                # ratio
                adaptive = vals[2]
                denom = math.inf
                for x in (vals[0], vals[1]):
                    if x is None or math.isnan(x):
                        continue
                    denom = min(denom, x)
                if adaptive is not None and not math.isnan(adaptive) and denom != math.inf:
                    ratio_mat[i][j] = adaptive / denom
                else:
                    ratio_mat[i][j] = math.nan

        # Write CSV matrices (always)
        winner_csv = os.path.join(figs_dir, f"exp6_phase_winner_{method}.csv")
        ratio_csv  = os.path.join(figs_dir, f"exp6_ratio_{method}.csv")
        write_matrix_csv(winner_csv, [str(a) for a in alphas], [str(t) for t in ts], winner_mat)
        write_matrix_csv(ratio_csv,  [str(a) for a in alphas], [str(t) for t in ts], ratio_mat)

        if plt is None:
            print(f"[plot_exp6] Wrote CSV only (no PNG): {winner_csv}, {ratio_csv}")
            continue

        # Phase diagram (imshow)
        plt.figure()
        plt.imshow(winner_mat, aspect="auto", interpolation="nearest")
        plt.yticks(range(len(alphas)), [str(a) for a in alphas])
        plt.xticks(range(len(ts)), [str(t) for t in ts], rotation=45, ha="right")
        plt.xlabel("t")
        plt.ylabel("alpha")
        plt.title(f"EXP-6 Phase diagram (winner: 0 sampling / 1 enum_sampling / 2 adaptive) — method={method}")
        plt.colorbar()
        plt.tight_layout()
        plt.savefig(os.path.join(figs_dir, f"exp6_phase_{method}.png"), dpi=200)

        # Ratio heatmap
        plt.figure()
        plt.imshow(ratio_mat, aspect="auto", interpolation="nearest")
        plt.yticks(range(len(alphas)), [str(a) for a in alphas])
        plt.xticks(range(len(ts)), [str(t) for t in ts], rotation=45, ha="right")
        plt.xlabel("t")
        plt.ylabel("alpha")
        plt.title(f"EXP-6 ratio = T(adaptive)/min(T(sampling),T(enum_sampling)) — method={method}")
        plt.colorbar()
        plt.tight_layout()
        plt.savefig(os.path.join(figs_dir, f"exp6_ratio_{method}.png"), dpi=200)

        print(f"[plot_exp6] Wrote PNGs and CSVs for method={method} into {figs_dir}")

    print("[plot_exp6] done.")

if __name__ == "__main__":
    main()
PY

echo "[EXP6] Plotting ..."
python3 "${PLOT_PY}" "${SUMMARY_CSV}" "${EXP6_OUT_DIR}/figs"

echo "[EXP6] DONE ✅"
echo "Results:"
echo "  ${EXP6_OUT_DIR}/sweep_raw.csv"
echo "  ${EXP6_OUT_DIR}/sweep_summary.csv"
echo "Figures:"
echo "  ${EXP6_OUT_DIR}/figs/exp6_phase_<method>.png"
echo "  ${EXP6_OUT_DIR}/figs/exp6_ratio_<method>.png"
