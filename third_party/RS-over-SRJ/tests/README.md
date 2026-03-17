# RS-over-SRJ third_party self-tests

These are lightweight, framework-free correctness tests for the vendored
RS-over-SRJ KDTree baseline.

They are not wired into the top-level CMake by default (to keep third_party
self-contained), but can be compiled manually.

From the **the main repository repository root**:

```bash
# KD-tree correctness on the legacy 4D self-test
g++ -std=c++17 -O2 \
  -Iinclude -Ithird_party/RS-over-SRJ \
  third_party/RS-over-SRJ/tests/test_kdtree_4d.cpp \
  -o /tmp/test_kdtree_4d && /tmp/test_kdtree_4d

# Baseline smoke test on a tiny rectangle dataset
g++ -std=c++17 -O2 \
  -Iinclude -Ithird_party/RS-over-SRJ \
  third_party/RS-over-SRJ/tests/test_rs_over_srj_baseline_small.cpp \
  -o /tmp/test_rs_over_srj_small && /tmp/test_rs_over_srj_small
```

Optional debug flags (useful when auditing correctness):

- `-DRS_OVER_SRJ_VALIDATE_SAMPLES=1` (fail fast if a non-intersecting sample is produced)
- `-DRS_OVER_SRJ_DISABLE_BATCHING=1` (mimic upstream per-sample control flow)
- `-DRS_OVER_SRJ_USE_UPSTREAM_ALIAS=1` (use self-contained Walker alias table)

