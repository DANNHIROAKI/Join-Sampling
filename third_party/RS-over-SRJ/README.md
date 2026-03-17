https://github.com/amgt-d1/RS-over-SRJ or https://github.com/DANNHIROAKI/RS-over-SRJ/releases/tag/snap-20260305

> **Upstream reference (audit / fairness):**
> - Source: *RS-over-SRJ* (KDTree/KDS baseline) by Daichi Amagata.
> - This directory is a vendored and adapted copy of the official KDTree baseline (trimmed package).
> - See `MODIFICATIONS.md` for an itemized change list and `UPSTREAM_DIFF.patch` for a full diff
>   against the provided upstream snapshot.
>
> **Build-time switches:**
> - `RS_OVER_SRJ_DISABLE_BATCHING=1` to disable batching (mimics upstream per-sample control flow).
> - `RS_OVER_SRJ_USE_UPSTREAM_ALIAS=1` to use a self-contained Walker alias implementation.
> - `RS_OVER_SRJ_VALIDATE_SAMPLES=1` to runtime-check every sampled pair (debug only).

# KDTree-(2*Dim) baseline

A static lifted KD-tree baseline for **exact** uniform random sampling over Dim-dimensional box intersection joins.

This directory contains a KDS-style baseline (range counting + range sampling) and a lightweight KD-tree implementation that supports orthogonal range **COUNT** and range **SEARCH** in configurable dimensionality. In this repository the join predicate is Dim-dimensional box strict intersection, expressed as a `2*Dim` orthogonal range query via the standard lifting.

---

## Problem

Given two relations of axis-aligned boxes in dimension `Dim`:

- $R = \{r_0,\dots,r_{n-1}\}$
- $S = \{s_0,\dots,s_{m-1}\}$

the intersection join result set is:

$J = \{(r, s)\mid r\in R,\ s\in S,\ r \cap s \neq \emptyset\}$.

The sampler returns $t$ pairs from $J$ **with replacement**.  
Each returned pair is an **independent uniform** sample from $J$.

---

## Rectangle semantics

### Half-open rectangles

Each rectangle is half-open on every axis:

$[lo, hi)$.

The input must satisfy, for each box $b$ and each axis $d\in\{0,\dots,\mathrm{Dim}-1\}$:

$b.lo[d] < b.hi[d]$.

All rectangle coordinates are finite IEEE-754 floating-point values (no NaN).

### Strict overlap (touching is not an intersection)

Rectangles that only touch at edges or corners are **not** considered intersecting.

For rectangles $a$ and $b$ in 2D, $a \cap b \neq \emptyset$ iff for each axis $d\in\{0,1\}$:

$a.lo[d] < b.hi[d]\ \land\ b.lo[d] < a.hi[d]$.

---

## Algorithm

For each rectangle $r\in R$, define:

$w(r) = |\{ s \in S \mid r \cap s \neq \emptyset \}|$,

and:

$W = \sum_{r\in R} w(r) = |J|$.

The sampler is:

1. Build a static KD-tree over embedded points of $S$ (Section “`2*Dim` lifting”).
2. For every $r\in R$, compute $w(r)$ by exact range counting on the KD-tree.
3. Build a Walker alias table over $R$ with weights $\{w(r)\}$.
4. Repeat $t$ times:
   - sample $r$ from $R$ with probability $w(r)/W$
   - sample $s$ uniformly from $\{s\in S \mid r\cap s\neq\emptyset\}$ using KD-tree range sampling
   - output $(r,s)$

### Uniformity proof

For any $(r,s)\in J$:

$P[(r,s)] = P[r]\cdot P[s\mid r] = \frac{w(r)}{W}\cdot \frac{1}{w(r)} = \frac{1}{|J|}$.

Independence follows from sampling with replacement using a random number generator stream.

---

## `2*Dim` lifting / embedding

Each box $s$ is embedded as a `2*Dim` point:

$p(s) = (s.lo[0],\dots,s.lo[\mathrm{Dim}-1],\ s.hi[0],\dots,s.hi[\mathrm{Dim}-1])\in\mathbb{R}^{2\mathrm{Dim}}$.

For a fixed query box $r$, strict overlap constraints
$s.lo[a] < r.hi[a]$ and $s.hi[a] > r.lo[a]$ for every axis $a$ form a `2*Dim` orthogonal range query over points $p(s)$, with:

- upper bounds on the lifted `lo` coordinates
- lower bounds on the lifted `hi` coordinates

---

## Closed-range conversion via `std::nextafter`

`utils/kdtree.hpp` uses **closed** ranges:

$qlo[d] \le x[d] \le qhi[d]$.

Strict inequalities are converted to closed bounds using `std::nextafter` on floating-point values:

- $x < a \iff x \le \mathrm{nextDown}(a)$
- $x > a \iff x \ge \mathrm{nextUp}(a)$

For rectangle $r$, the embedded 4D query bounds are:

- $qlo = (-\infty,\ -\infty,\ \mathrm{nextUp}(r.lo.x),\ \mathrm{nextUp}(r.lo.y))$
- $qhi = (\mathrm{nextDown}(r.hi.x),\ \mathrm{nextDown}(r.hi.y),\ +\infty,\ +\infty)$

`nextUp/nextDown` are identity on $\pm\infty$.

---

## KD-tree range counting and uniform range sampling

### KD-tree representation

`utils/kdtree.hpp` implements a static KD-tree with:

- a single point array `points_`
- each subtree mapped to a contiguous interval `[l, r)` inside `points_`
- per-node bounding boxes `(lo, hi)` computed bottom-up by merging children
- closed range semantics: $qlo \le x \le qhi$

Each KD-tree node stores:

- `l, r, mid`: subtree interval and median index
- `left, right`: child node indices
- `axis`: split axis
- `lo, hi`: bounding box of the subtree (inclusive)

### Exact range counting

`Count(qlo, qhi)` returns the exact number of embedded points in the closed range.

A subtree is pruned when its bounding box is disjoint from the query.  
A subtree contributes its full size when its bounding box is contained in the query.

Worst-case time is $O(|S|)$.

### Range decomposition for uniform sampling

`Search(qlo, qhi, out)` returns a decomposition of the query result set into items:

- **subtree items** (`fully_contained=true`): a fully contained subtree block
- **point items** (`fully_contained=false`): the median point stored at a node that lies in range

This decomposition is a partition of the result set with no duplication.

Uniform sampling from the query result is implemented by:

1. assign an integer weight to each item:
   - subtree item: `weight = subtree_size`
   - point item: `weight = 1`
2. sample one item index via an alias table using these weights
3. if a subtree item is chosen, sample a uniform offset inside the subtree interval `[l, r)` and return that point
4. if a point item is chosen, return the node’s median point

Because each result point belongs to exactly one item, and each item is selected proportionally to its size contribution, the sampled point is uniform over the query result set.

---

## Baseline class and API

The join sampler is implemented in `isrjs_kds/isrjs_kds.hpp`:

` sjs::baselines::rs_over_srj::RSOverSRJKDTreeSamplingBaseline<2, sjs::Scalar> `

Properties:

- input dimension: `Dim in [2,5]` in this repository build
- embedding dimension: `D = 2 * Dim`
- `Build(ds, cfg)`: builds a KD-tree over embedded rectangles in `ds.S`
- `Count(cfg, rng, out)`: computes exact $|J|$ by summing $w(r)$ over `ds.R`
- `Sample(cfg, rng, out)`: produces `t = cfg.run.t` samples with replacement
- `Name()`: `"kd_tree_sampling"`
- `Enumerate(...)`: deterministic plane-sweep enumeration on axis 0 (verification / enumerator API)

Hard limits:

- `|S| < 2^32` (KD-tree indices and payloads are `uint32_t`)
- `t < 2^32` (sample count is stored as `uint32_t` inside `Sample()`)
- `|J| < 2^64` (otherwise `Count()` reports overflow)

Precondition:

- $|J| \ge 1$.

---

## Repository layout

```
RS-over-SRJ/
  isrjs_kds/
    isrjs_kds.hpp         # Baseline implementation (lifted embedding + sampling)
    main.cpp              # Sanity-check demo (guarded by RS_OVER_SRJ_DEMO_MAIN)
    parameter/
      dataset_id.txt      # Legacy file (not used by RSOverSRJKDTreeSamplingBaseline)
      range.txt           # Legacy file (not used by RSOverSRJKDTreeSamplingBaseline)
      sample_size.txt     # Legacy file (not used by RSOverSRJKDTreeSamplingBaseline)
      scalability.txt     # Legacy file (not used by RSOverSRJKDTreeSamplingBaseline)
  utils/
    kdtree.hpp            # Static KD-tree with closed-range count/search
    weighted_sampling.hpp # Wrapper around sjs::sampling::AliasTable
    utils.hpp             # Process RSS helper
    csv.h                 # Third-party CSV parser (BSD-3-Clause; license in file header)
    pcg_random.hpp        # Third-party RNG (Apache-2.0; license in file header)
    pcg_extras.hpp        # PCG support header (Apache-2.0; license in file header)
    pcg_uint128.hpp       # PCG support header (Apache-2.0; license in file header)
```

---

## Build and run: sanity-check demo

`isrjs_kds/main.cpp` is a minimal demo guarded by `RS_OVER_SRJ_DEMO_MAIN`.

This codebase is designed to be compiled inside the main repository source tree. The directory layout assumed by the command below is:

- repository root: `./`
- this directory: `./third_party/RS-over-SRJ/`

Build and run:

~~~bash
g++ -std=c++17 -O3 -DRS_OVER_SRJ_DEMO_MAIN \
  -I. -Ithird_party/RS-over-SRJ \
  third_party/RS-over-SRJ/isrjs_kds/main.cpp \
  -o rs_over_srj_demo

./rs_over_srj_demo
~~~

The executable constructs a small in-memory rectangle dataset, prints the exact join size, then prints `t = 10` sampled pairs as `(r_id, s_id)`.

---

## License

### KDTree-(2*Dim) baseline code

MIT License

~~~text
Copyright 2025 Daichi Amagata

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
~~~

### Third-party code

- `utils/csv.h`: BSD-3-Clause license text is included in the file header.
- `utils/pcg_random.hpp`, `utils/pcg_extras.hpp`, `utils/pcg_uint128.hpp`: Apache-2.0 license text is included in the file headers.

---

## Citation

~~~bibtex
@inproceedings{amagata2025random,
  title={Random Sampling over Spatial Range Joins},
  author={Amagata, Daichi},
  booktitle={ICDE},
  pages={2080--2093},
  year={2025}
}
~~~