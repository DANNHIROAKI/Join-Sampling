# Modifications relative to the upstream RS-over-SRJ KDTree/KDS baseline

This directory vendors and adapts the **KDTree / KDS baseline** from the paper:

> Daichi Amagata, *Random Sampling over Spatial Range Joins*, ICDE 2025.

Upstream project link (reference):
- https://github.com/amgt-d1/RS-over-SRJ

The version we compared against is the **trimmed upstream snapshot** provided as `RS-over-SRJ-KDS.zip` (KDTree baseline only).

The intent is to keep the **algorithmic core** identical (KD-tree range counting + KD-tree range sampling + alias sampling), while:
1) fixing correctness bugs present in the trimmed upstream snapshot, and
2) upgrading the implementation from the upstream point range-join demo to the **(2*Dim)-dimensional lifting/embedding** required for box intersection joins in Dim=2/3/4/5.

This file documents every non-trivial deviation so that reviewers can audit fairness.

---

## What stays the same (algorithm)

The baseline still implements the same high-level sampler described in the RS-over-SRJ paper baseline (Section III-A):

1. Build a static KD-tree on the “S side” objects.
2. For each “R side” object, compute its exact overlap count (range counting).
3. Build a Walker alias table over R with weights equal to these counts.
4. For each output sample:
   - draw r ∼ w(r) / Σ w(r)
   - draw s uniformly from the overlap set S(r) via KD-tree range sampling

This yields **exact** uniform sampling over the join result with replacement.

---

## Changes (implementation / integration)

### 1) From a standalone demo to a library baseline

Upstream code is a standalone executable with global datasets/parameters and file I/O.

In this repository, the baseline is integrated as:

- `isrjs_kds/isrjs_kds.hpp` : `sjs::baselines::rs_over_srj::RSOverSRJKDTreeSamplingBaseline<2, Scalar>`
- Uses the repository’s `sjs::Dataset`, `sjs::IBaseline`, and `sjs::Rng`.

This removes hard-coded paths and avoids ODR issues from global variables.

### 2) From point range join to box intersection join via `2*Dim` lifting

The upstream KDTree baseline is demonstrated for a point/range setting.

Our experiments target a **Dim-dimensional box strict intersection join**, expressed as a **`2*Dim`-dimensional orthogonal range query** using the standard lifting:

- Embed each box `s ∈ S` as a `2*Dim` point:
  `p(s) = (s.lo[0], ..., s.lo[Dim-1], s.hi[0], ..., s.hi[Dim-1])`.
- For a query box `r`, strict intersection constraints become:
  `s.lo[a] < r.hi[a]` and `s.hi[a] > r.lo[a]` for every axis `a`.

Because the KD-tree code uses **closed** ranges (`qlo ≤ x ≤ qhi`), strict inequalities are converted using `std::nextafter`:

- `x < a  ⇔  x ≤ nextDown(a)`
- `x > a  ⇔  x ≥ nextUp(a)`

This matches the half-open box semantics used throughout the main repository:
`[lo, hi)` and strict overlap (touching edges is not an intersection).

### 3) KD-tree rewrite to fix subtree index invariants and sampling correctness

The trimmed upstream snapshot includes a pointer-based KD-tree that maintains
`left_idx/right_idx` intervals.

That implementation can violate subtree interval invariants after `nth_element`
partitioning (and when recursing on pointer offsets), which can lead to:

- incorrect range counts (range counting mismatch), and/or
- range sampling returning points outside the query range.

To make the baseline correct and extensible to general `D = 2*Dim`, we replaced the KD-tree with
a lightweight static KD-tree that:

- stores all points in a single array `points_`
- maps each subtree to a contiguous interval `[l, r)` inside that array
- stores the median point at `mid` and recurses on `[l, mid)` and `[mid+1, r)`
- computes bounding boxes bottom-up by merging children

The interface still supports exactly the two operations required by the paper baseline:

- `Count(qlo, qhi)` : exact orthogonal range counting
- `Search(qlo, qhi, out)` : returns a decomposition into
  - fully-contained subtrees (weight = subtree size)
  - single median points (weight = 1)

This decomposition is used for **exact uniform sampling** from the query result.

### 4) Alias table & RNG integration (with an opt-in upstream-like implementation)

Upstream uses a local alias-table implementation coupled to its own RNG stream.

For reproducibility across baselines, this repository defaults to using the
project’s alias table + RNG:

- `utils/weighted_sampling.hpp` wraps `sjs::sampling::AliasTable` and samples using `sjs::Rng*`.

To reduce “you changed their code too much” concerns, we also provide an
optional self-contained Walker alias implementation:

- define `RS_OVER_SRJ_USE_UPSTREAM_ALIAS=1` when compiling

Both variants implement the same sampling distribution.

### 5) Optional batching optimization (can be disabled)

`RSOverSRJKDTreeSamplingBaseline::Sample()` defaults to batching multiple output
slots that happened to draw the same `ridx`.

This amortizes a KD-tree query and per-query alias construction. It is a pure
performance optimization and does **not** change the sampling distribution.

If you want to mimic the upstream per-sample control flow more closely, compile with:

- `RS_OVER_SRJ_DISABLE_BATCHING=1`

### 6) Debug validation switch

For debugging and reviewer reassurance, you can compile with:

- `RS_OVER_SRJ_VALIDATE_SAMPLES=1`

This performs an explicit rectangle intersection check for every sampled pair
and fails fast if a non-intersecting pair is produced.

---

## Diff against upstream snapshot

A full diff against the provided upstream snapshot is included as:

- `UPSTREAM_DIFF.patch`

This is intended for auditing and is not required for building.

