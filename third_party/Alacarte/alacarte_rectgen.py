
"""
alacarte_rectgen.py

A practical, reproducible synthetic axis-aligned box/rectangle generator
in the spirit of "Benchmarking Spatial Joins À La Carte".

Goal:
- Generate two box sets R, S with specified |R|, |S|
- Control output density alpha_out = |J| / (|R|+|S|),
  where J = {(r,s): r in R, s in S, r intersects s}.

Key knob:
- coverage C := (sum of box volumes) / (universe volume)
  For a set of size n, mean box volume is v_bar = C * U / n.

We solve for coverage C that matches the requested alpha_out (in expectation),
then generate R and S with that coverage.

This implementation is dimension-agnostic (d>=1). In 2D it produces rectangles.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Dict, List, Optional, Tuple, Literal

import numpy as np

VolumeDist = Literal["fixed", "exponential", "normal", "lognormal"]


# ----------------------------
# Data container
# ----------------------------

@dataclass(frozen=True)
class BoxSet:
    """
    Stores n axis-aligned half-open boxes in d dimensions:
        box i = Π_k [lower[i,k], upper[i,k])
    """
    lower: np.ndarray  # shape (n, d)
    upper: np.ndarray  # shape (n, d)
    universe: np.ndarray  # shape (d, 2)

    @property
    def n(self) -> int:
        return int(self.lower.shape[0])

    @property
    def d(self) -> int:
        return int(self.lower.shape[1])

    def as_array(self) -> np.ndarray:
        """Return a (n, 2d) array: [L0..Ld-1, R0..Rd-1]."""
        return np.hstack([self.lower, self.upper])


# ----------------------------
# Helpers
# ----------------------------

def ensure_universe(universe: Optional[np.ndarray], d: int) -> np.ndarray:
    """
    universe: array-like of shape (d,2) giving [min, max] per dimension.
    Default is unit hypercube [0,1]^d.
    """
    if universe is None:
        U = np.array([[0.0, 1.0]] * d, dtype=np.float64)
        return U
    U = np.asarray(universe, dtype=np.float64)
    if U.shape != (d, 2):
        raise ValueError(f"universe must have shape (d,2); got {U.shape}")
    if np.any(U[:, 1] <= U[:, 0]):
        raise ValueError("universe upper must be > lower in every dim")
    return U


def sample_volumes(
    n: int,
    mean_vol: float,
    dist: VolumeDist,
    rng: np.random.Generator,
    cv: float = 0.25,
) -> np.ndarray:
    """
    Sample positive volumes with E[V]=mean_vol (approximately / exactly depends on dist).
    cv is coefficient of variation for normal/lognormal.
    """
    mean_vol = float(mean_vol)
    if mean_vol <= 0:
        return np.full(n, 1e-18, dtype=np.float64)

    if dist == "fixed":
        return np.full(n, mean_vol, dtype=np.float64)

    if dist == "exponential":
        # E = mean_vol
        return rng.exponential(scale=mean_vol, size=n).astype(np.float64)

    if dist == "normal":
        # Truncated normal: clip to positive
        sigma = cv * mean_vol
        vols = rng.normal(loc=mean_vol, scale=sigma, size=n).astype(np.float64)
        return np.clip(vols, a_min=mean_vol * 1e-12, a_max=None)

    if dist == "lognormal":
        # choose sigma so that cv matches: cv^2 = exp(sigma^2)-1
        sigma = math.sqrt(math.log1p(cv**2))
        mu = math.log(mean_vol) - 0.5 * sigma**2
        return rng.lognormal(mean=mu, sigma=sigma, size=n).astype(np.float64)

    raise ValueError(f"Unknown volume distribution: {dist}")


def sample_shape_factors(
    n: int,
    d: int,
    shape_sigma: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Shape factors g[i,k] with product over k equals 1.
    Then side lengths = (V^(1/d)) * g.
    - shape_sigma=0 => hypercubes (all g=1).
    - shape_sigma>0 => log-normal spread of aspect ratios.
    """
    if shape_sigma <= 0:
        return np.ones((n, d), dtype=np.float64)

    z = rng.normal(loc=0.0, scale=shape_sigma, size=(n, d)).astype(np.float64)
    g = np.exp(z)
    # Normalize per row so that Π_k g_k = 1
    g /= np.prod(g, axis=1, keepdims=True) ** (1.0 / d)
    return g


def sample_lengths_from_mean(
    n_samples: int,
    d: int,
    universe: np.ndarray,
    mean_vol: float,
    volume_dist: VolumeDist,
    volume_cv: float,
    shape_sigma: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Sample side lengths (n_samples, d) given mean volume.
    """
    spans = universe[:, 1] - universe[:, 0]
    vols = sample_volumes(n_samples, mean_vol, volume_dist, rng, cv=volume_cv)
    g = sample_shape_factors(n_samples, d, shape_sigma, rng)
    base = vols ** (1.0 / d)
    lengths = base[:, None] * g

    # Cap to slightly below spans to avoid degeneracy in formulas.
    eps = 1e-12
    lengths = np.minimum(lengths, spans * (1.0 - eps))
    return lengths


# ----------------------------
# Exact 1D overlap probability (given lengths)
# ----------------------------

def interval_overlap_prob(a: np.ndarray, b: np.ndarray, W: float) -> np.ndarray:
    """
    Exact probability that two half-open intervals overlap, when
      X ~ Uniform(0, W-a), Y ~ Uniform(0, W-b), independent,
      intervals are [X, X+a), [Y, Y+b).

    For a+b >= W: overlap probability is 1.
    For a+b <  W: p = 1 - (W-(a+b))^2 / ((W-a)(W-b))

    a, b can be arrays (broadcastable). W is scalar (span of universe).
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    W = float(W)

    a = np.clip(a, 0.0, W)
    b = np.clip(b, 0.0, W)

    c = W - (a + b)
    p = np.ones(np.broadcast(a, b).shape, dtype=np.float64)

    mask = c > 0.0
    if np.any(mask):
        denom = (W - a) * (W - b)  # positive when mask holds
        p[mask] = 1.0 - (c[mask] ** 2) / denom[mask]
    return p


def estimate_alpha_expected(
    nR: int,
    nS: int,
    alpha_out: float,
    d: int,
    universe: Optional[np.ndarray],
    coverage: float,
    volume_dist: VolumeDist,
    volume_cv: float,
    shape_sigma: float,
    num_samples: int,
    seed: int,
) -> Tuple[float, float]:
    """
    Estimate expected alpha_out for a given coverage using:
      p = E[ Π_k p1D(lenR_k, lenS_k) ]
    where p1D() is exact 1D overlap probability under uniform placement.

    Returns (alpha_est, p_est).
    """
    U = ensure_universe(universe, d)
    spans = U[:, 1] - U[:, 0]
    U_vol = float(np.prod(spans))

    mean_vol_R = coverage * U_vol / nR
    mean_vol_S = coverage * U_vol / nS

    rngR = np.random.default_rng(seed)
    rngS = np.random.default_rng(seed + 1)

    lenR = sample_lengths_from_mean(num_samples, d, U, mean_vol_R, volume_dist, volume_cv, shape_sigma, rngR)
    lenS = sample_lengths_from_mean(num_samples, d, U, mean_vol_S, volume_dist, volume_cv, shape_sigma, rngS)

    p = np.ones(num_samples, dtype=np.float64)
    for k in range(d):
        p *= interval_overlap_prob(lenR[:, k], lenS[:, k], spans[k])

    p_est = float(p.mean())
    alpha_est = p_est * (nR * nS) / (nR + nS)
    return alpha_est, p_est


def solve_coverage_for_alpha(
    nR: int,
    nS: int,
    alpha_target: float,
    *,
    d: int = 2,
    universe: Optional[np.ndarray] = None,
    volume_dist: VolumeDist = "fixed",
    volume_cv: float = 0.25,
    shape_sigma: float = 0.0,
    num_samples: int = 200_000,
    seed: int = 0,
    tol_rel: float = 0.02,
    max_iter: int = 30,
) -> Tuple[float, List[Dict[str, float]]]:
    """
    Find a coverage C such that expected alpha_out ~= alpha_target (relative tolerance).

    We use bracketing + log-space binary search. The objective is monotone in C.

    Returns (coverage, history), where history records alpha_est for tried coverages.
    """
    if alpha_target < 0:
        raise ValueError("alpha_target must be >= 0")
    if alpha_target == 0:
        return 0.0, [{"coverage": 0.0, "alpha_est": 0.0}]

    alpha_max = (nR * nS) / (nR + nS)  # since |J| <= nR*nS
    if alpha_target > alpha_max:
        raise ValueError(
            f"alpha_target={alpha_target} exceeds max possible {alpha_max:.6g} "
            f"for given |R|={nR}, |S|={nS}."
        )

    # Initial guess: for small rectangles and nR≈nS, alpha_out ≈ 2C
    C0 = max(1e-12, alpha_target / 2.0)

    def eval_alpha(C: float) -> float:
        a_est, _ = estimate_alpha_expected(
            nR=nR, nS=nS, alpha_out=alpha_target, d=d, universe=universe,
            coverage=C, volume_dist=volume_dist, volume_cv=volume_cv, shape_sigma=shape_sigma,
            num_samples=num_samples, seed=seed
        )
        return a_est

    alpha0 = eval_alpha(C0)
    history: List[Dict[str, float]] = [{"coverage": C0, "alpha_est": alpha0}]

    if abs(alpha0 - alpha_target) <= tol_rel * alpha_target:
        return C0, history

    # Bracket [lo, hi] s.t. alpha(lo) < target <= alpha(hi)
    lo = 0.0
    hi = C0
    if alpha0 < alpha_target:
        lo = C0
        hi = C0
        for _ in range(60):
            hi *= 2.0
            a_hi = eval_alpha(hi)
            history.append({"coverage": hi, "alpha_est": a_hi})
            if a_hi >= alpha_target:
                break
            lo = hi
        else:
            raise RuntimeError("Failed to bracket alpha_target; try increasing max bound or check parameters.")
    else:
        # alpha(0)=0 <= target < alpha(C0)
        lo = 0.0
        hi = C0

    # Log-space binary search
    lo_eps = 1e-15
    for _ in range(max_iter):
        lo_for_log = max(lo, lo_eps)
        mid = math.exp(0.5 * (math.log(lo_for_log) + math.log(hi)))
        a_mid = eval_alpha(mid)
        history.append({"coverage": mid, "alpha_est": a_mid})

        if abs(a_mid - alpha_target) <= tol_rel * alpha_target:
            return mid, history

        if a_mid < alpha_target:
            lo = mid
        else:
            hi = mid

    best = min(history, key=lambda h: abs(h["alpha_est"] - alpha_target))
    return float(best["coverage"]), history


# ----------------------------
# Final box generation
# ----------------------------

def generate_boxset(
    n: int,
    *,
    d: int = 2,
    universe: Optional[np.ndarray] = None,
    coverage: float = 1.0,
    volume_dist: VolumeDist = "fixed",
    volume_cv: float = 0.25,
    shape_sigma: float = 0.0,
    seed: int = 0,
    dtype: np.dtype = np.float32,
) -> BoxSet:
    """
    Generate a set of n boxes in the given universe with target coverage (in expectation).
    Placement: lower corner uniform so the box fits (same assumption as overlap formula).
    """
    U = ensure_universe(universe, d)
    spans = U[:, 1] - U[:, 0]
    U_vol = float(np.prod(spans))

    mean_vol = coverage * U_vol / n
    rng = np.random.default_rng(seed)

    lengths = sample_lengths_from_mean(n, d, U, mean_vol, volume_dist, volume_cv, shape_sigma, rng)
    # lower corner uniform so it fits (float64)
    u = rng.random((n, d), dtype=np.float64)
    lower64 = U[:, 0] + u * (spans - lengths)
    upper64 = lower64 + lengths
    
    # Cast to requested dtype and enforce strict lower < upper **in that dtype**.
    # This is important for float32: if boxes become extremely small (e.g., alpha_out=0),
    # rounding can make upper==lower, violating L<R invariants.
    U_d = U.astype(dtype, copy=False)
    lo = U_d[:, 0]
    hi = U_d[:, 1]
    
    lower = lower64.astype(dtype, copy=False)
    upper = upper64.astype(dtype, copy=False)

    # Keep coordinates inside the universe in the output dtype.
    lower = np.maximum(lower, lo)
    upper = np.minimum(upper, hi)

    # Ensure lower < hi so there exists a representable value above it.
    hi_prev = np.nextafter(hi, lo)  # largest representable value < hi (per dimension)
    lower = np.minimum(lower, hi_prev)

    # Fix any degeneracy caused by dtype rounding: enforce upper > lower by at least one ULP.
    hi_b = np.broadcast_to(hi, lower.shape)
    upper_fix = np.nextafter(lower, hi_b)
    mask = upper <= lower
    if np.any(mask):
        upper = np.where(mask, upper_fix, upper)

    return BoxSet(lower=lower, upper=upper, universe=U_d)



def make_rectangles_R_S(
    nR: int,
    nS: int,
    alpha_out: float,
    *,
    d: int = 2,
    universe: Optional[np.ndarray] = None,
    volume_dist: VolumeDist = "fixed",
    volume_cv: float = 0.25,
    shape_sigma: float = 0.0,
    tune_samples: int = 200_000,
    tune_tol_rel: float = 0.02,
    tune_max_iter: int = 30,
    seed: int = 0,
    dtype: np.dtype = np.float32,
) -> Tuple[BoxSet, BoxSet, Dict]:
    """
    Main entry:
      - solve coverage C for target alpha_out (expected)
      - generate R and S with that coverage

    Returns (R, S, info) where info includes coverage and tuning history.
    """
    C, history = solve_coverage_for_alpha(
        nR=nR, nS=nS, alpha_target=alpha_out,
        d=d, universe=universe,
        volume_dist=volume_dist, volume_cv=volume_cv, shape_sigma=shape_sigma,
        num_samples=tune_samples, seed=seed + 10_000,
        tol_rel=tune_tol_rel, max_iter=tune_max_iter,
    )

    R = generate_boxset(
        nR, d=d, universe=universe, coverage=C,
        volume_dist=volume_dist, volume_cv=volume_cv, shape_sigma=shape_sigma,
        seed=seed + 1, dtype=dtype,
    )
    S = generate_boxset(
        nS, d=d, universe=universe, coverage=C,
        volume_dist=volume_dist, volume_cv=volume_cv, shape_sigma=shape_sigma,
        seed=seed + 2, dtype=dtype,
    )

    # report expected alpha using the tuned C (for logging)
    alpha_est, p_est = estimate_alpha_expected(
        nR=nR, nS=nS, alpha_out=alpha_out, d=d, universe=universe,
        coverage=C, volume_dist=volume_dist, volume_cv=volume_cv, shape_sigma=shape_sigma,
        num_samples=min(200_000, tune_samples), seed=seed + 20_000,
    )

    info = {
        "coverage": float(C),
        "alpha_target": float(alpha_out),
        "alpha_expected_est": float(alpha_est),
        "pair_intersection_prob_est": float(p_est),
        "tune_history": history,
        "params": {
            "nR": int(nR),
            "nS": int(nS),
            "d": int(d),
            "universe": None if universe is None else np.asarray(universe).tolist(),
            "volume_dist": volume_dist,
            "volume_cv": float(volume_cv),
            "shape_sigma": float(shape_sigma),
            "seed": int(seed),
        },
    }
    return R, S, info


# ----------------------------
# Optional: quick sanity check (pair sampling on generated sets)
# ----------------------------

def estimate_alpha_by_pair_sampling(
    R: BoxSet,
    S: BoxSet,
    *,
    num_pairs: int = 1_000_000,
    seed: int = 0,
) -> Tuple[float, float]:
    """
    Estimate realized alpha_out by random pair sampling on the concrete generated sets.
    Warning: for very large N and very small alpha_out, true p is tiny (≈O(1/N)),
    so you may need a lot of pairs to observe enough intersections.
    """
    rng = np.random.default_rng(seed)
    idxR = rng.integers(0, R.n, size=num_pairs)
    idxS = rng.integers(0, S.n, size=num_pairs)

    lower_max = np.maximum(R.lower[idxR], S.lower[idxS])
    upper_min = np.minimum(R.upper[idxR], S.upper[idxS])
    hits = np.all(lower_max < upper_min, axis=1).sum()

    p_hat = hits / num_pairs
    alpha_hat = p_hat * (R.n * S.n) / (R.n + S.n)
    return float(alpha_hat), float(p_hat)


if __name__ == "__main__":
    # Example
    nR, nS = 5000, 5000
    alpha = 10.0

    R, S, info = make_rectangles_R_S(
        nR=nR, nS=nS, alpha_out=alpha,
        d=2,
        volume_dist="fixed",
        shape_sigma=0.0,
        seed=0,
        tune_samples=200_000,
        tune_tol_rel=0.01,
    )
    print("Solved coverage:", info["coverage"])
    print("Expected alpha (MC over lengths):", info["alpha_expected_est"])

    alpha_ps, p_ps = estimate_alpha_by_pair_sampling(R, S, num_pairs=2_000_000, seed=123)
    print("Alpha by pair-sampling on generated sets:", alpha_ps, "p:", p_ps)
