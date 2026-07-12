"""Bayesian optimisation of RNAmotifs Discovery-mode parameters.

Replaces the exhaustive grid search over (hw, ew) with a Gaussian-process-
based Bayesian optimisation loop using Expected Improvement (EI).

Implementation
    The BO loop is self-contained (see ``optimise_rbp``): a Latin-hypercube
    initial design (custom NumPy routine), a Gaussian-process surrogate re-fit
    each iteration, and an
    Expected-Improvement acquisition maximised by multi-start scipy L-BFGS-B.
    Only the GP surrogate comes from a library; the design, acquisition and
    optimisation loop are all hand-rolled here.

    GP surrogate: ``sklearn.gaussian_process.GaussianProcessRegressor`` with a
    Matern-5/2 kernel that assigns n and e their own length-scales (per-dimension
    scaling — "ARD" in the 2-D sense, i.e. one length-scale per axis, not
    high-dimensional feature selection).  scikit-learn is the only hard
    requirement (>= 1.0).

    scikit-optimize (``skopt``) is supported as an optional drop-in for the
    surrogate class only (``skopt.learning.GaussianProcessRegressor``); it is
    NOT used as an optimiser and is not required.  When skopt is absent the code
    uses scikit-learn — which was the case for the reported benchmark.  A GP
    backend was preferred over BoTorch/GPyOpt to avoid a PyTorch dependency and
    because a lightweight surrogate is sufficient for this 2-D, deterministic,
    expensive objective.
"""

from __future__ import annotations

import logging
import warnings
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize as scipy_minimize
from scipy.stats import norm

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# GP surrogate -- try skopt first, fall back to sklearn
# ---------------------------------------------------------------------------

_BACKEND: Optional[str] = None

try:
    from skopt.learning import GaussianProcessRegressor as _SkoptGPR
    from skopt.learning.gaussian_process.kernels import Matern as _SkoptMatern
    _BACKEND = "skopt"
except ImportError:
    pass

if _BACKEND is None:
    try:
        from sklearn.gaussian_process import GaussianProcessRegressor as _SklearnGPR
        from sklearn.gaussian_process.kernels import (
            Matern as _SklearnMatern,
            ConstantKernel,
            WhiteKernel,
        )
        _BACKEND = "sklearn"
    except ImportError:
        pass

if _BACKEND is None:
    raise ImportError(
        "Bayesian optimisation requires scikit-optimize or scikit-learn.  "
        "Install one:  pip install scikit-optimize   OR   pip install scikit-learn"
    )


def _build_gp(n_dims: int = 2, random_state: int = 0):
    """Return a fresh GP regressor with Matern-5/2 kernel and ARD."""
    if _BACKEND == "skopt":
        kernel = _SkoptMatern(
            length_scale=[1.0] * n_dims,
            length_scale_bounds=[(1e-3, 1e3)] * n_dims,
            nu=2.5,
        )
        return _SkoptGPR(
            kernel=kernel,
            alpha=1e-6,
            normalize_y=True,
            n_restarts_optimizer=5,
            random_state=random_state,
        )
    # sklearn fallback
    kernel = (
        ConstantKernel(1.0, (1e-3, 1e3))
        * _SklearnMatern(
            length_scale=[1.0] * n_dims,
            length_scale_bounds=[(1e-3, 1e3)] * n_dims,
            nu=2.5,
        )
        + WhiteKernel(noise_level=1e-5, noise_level_bounds=(1e-10, 1e-1))
    )
    return _SklearnGPR(
        kernel=kernel,
        alpha=0.0,
        normalize_y=True,
        n_restarts_optimizer=5,
        random_state=random_state,
    )


# ---------------------------------------------------------------------------
# Latin Hypercube Sampling
# ---------------------------------------------------------------------------

def _lhs(n_samples: int, bounds: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Simple Latin Hypercube Sampling in ``[lo, hi]`` per dimension.

    Parameters
    ----------
    n_samples : int
        Number of design points.
    bounds : ndarray, shape (d, 2)
        Lower and upper bounds per dimension.
    rng : Generator
        Numpy random generator.

    Returns
    -------
    ndarray, shape (n_samples, d)
    """
    d = bounds.shape[0]
    result = np.empty((n_samples, d))
    for j in range(d):
        lo, hi = bounds[j]
        perm = rng.permutation(n_samples)
        cuts = np.linspace(0, 1, n_samples + 1)
        u = cuts[:-1] + rng.uniform(size=n_samples) * (cuts[1:] - cuts[:-1])
        result[:, j] = lo + u[perm] * (hi - lo)
    return result


# ---------------------------------------------------------------------------
# Expected Improvement
# ---------------------------------------------------------------------------

def _expected_improvement(
    X: np.ndarray,
    gp,
    y_best: float,
    xi: float = 0.01,
) -> np.ndarray:
    """Compute Expected Improvement at points *X*.

    Parameters
    ----------
    X : ndarray, shape (n, d)
    gp : fitted GP regressor
    y_best : float
        Best observed objective value so far.
    xi : float
        Exploration–exploitation trade-off.

    Returns
    -------
    ndarray, shape (n,)
        EI values (non-negative).
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mu, sigma = gp.predict(X, return_std=True)
    sigma = np.maximum(sigma, 1e-9)
    z = (mu - y_best - xi) / sigma
    ei = (mu - y_best - xi) * norm.cdf(z) + sigma * norm.pdf(z)
    ei[sigma < 1e-9] = 0.0
    return ei


# ---------------------------------------------------------------------------
# EI maximisation via multi-start L-BFGS-B
# ---------------------------------------------------------------------------

def _maximise_ei(
    gp,
    y_best: float,
    bounds: np.ndarray,
    evaluated: np.ndarray,
    xi: float = 0.01,
    n_restarts: int = 1000,
    rng: np.random.Generator = None,
) -> np.ndarray:
    """Find the point that maximises EI, rounded to integers, avoiding duplicates.

    Parameters
    ----------
    gp : fitted GP
    y_best : float
    bounds : ndarray, shape (d, 2)
    evaluated : ndarray, shape (m, d)
        Already-evaluated integer points (to skip).
    xi : float
    n_restarts : int
        Random starting points for L-BFGS-B.
    rng : Generator

    Returns
    -------
    ndarray, shape (d,)  — integer candidate.
    """
    if rng is None:
        rng = np.random.default_rng()

    d = bounds.shape[0]
    scipy_bounds = list(zip(bounds[:, 0], bounds[:, 1]))

    # Pre-screen: evaluate EI on random points, pick top starts
    X_random = rng.uniform(bounds[:, 0], bounds[:, 1], size=(n_restarts, d))
    ei_random = _expected_improvement(X_random, gp, y_best, xi)
    top_idx = np.argsort(ei_random)[-min(20, n_restarts):]

    best_ei = -1.0
    best_x = None

    for idx in top_idx:
        x0 = X_random[idx]
        try:
            res = scipy_minimize(
                lambda x: -_expected_improvement(x.reshape(1, -1), gp, y_best, xi)[0],
                x0,
                bounds=scipy_bounds,
                method="L-BFGS-B",
            )
            if -res.fun > best_ei:
                best_ei = -res.fun
                best_x = res.x
        except Exception:
            continue

    if best_x is None:
        # Fallback: best random point
        best_x = X_random[np.argmax(ei_random)]

    # Round to integer and ensure within bounds
    candidate = np.rint(best_x).astype(int)
    candidate = np.clip(candidate, bounds[:, 0].astype(int), bounds[:, 1].astype(int))

    # If already evaluated, try nearby integer grid points
    if evaluated.shape[0] > 0:
        max_attempts = 50
        for _ in range(max_attempts):
            if not any(np.array_equal(candidate, e) for e in evaluated):
                break
            # Perturb by +-1 in each dimension
            perturbation = rng.integers(-2, 3, size=d)
            candidate = np.rint(best_x + perturbation).astype(int)
            candidate = np.clip(
                candidate, bounds[:, 0].astype(int), bounds[:, 1].astype(int)
            )
        else:
            # Last resort: random unevaluated point
            for _ in range(500):
                candidate = rng.integers(
                    bounds[:, 0].astype(int),
                    bounds[:, 1].astype(int) + 1,
                    size=d,
                )
                if not any(np.array_equal(candidate, e) for e in evaluated):
                    break

    return candidate


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def optimise_rbp(
    target: str,
    eval_fn: Callable[[int, int], float],
    search_space: Optional[Dict[str, Tuple[int, int]]] = None,
    budget: Optional[Tuple[int, int]] = None,
    seed: int = 42,
    xi: float = 0.01,
    n_restarts: int = 1000,
) -> Dict:
    """Bayesian optimisation of (hw, ew) for a single RBP.

    Parameters
    ----------
    target : str
        RBP name (used only for logging).
    eval_fn : callable(hw: int, ew: int) -> float
        Black-box objective that returns AUROC for a given (hw, ew).
        Must use the **same** AUROC computation as the grid search.
    search_space : dict, optional
        ``{'hw': (lo, hi), 'ew': (lo, hi)}``.
        Default: ``{'hw': (5, 100), 'ew': (20, 400)}``.
    budget : (n_init, n_iter), optional
        Number of initial LHS points and BO iterations.
        Default: ``(8, 30)``.
    seed : int
        Random seed for reproducibility.
    xi : float
        EI exploration parameter (default 0.01).
    n_restarts : int
        Multi-start count for EI maximisation (default 1000).

    Returns
    -------
    dict
        ``n_star``    : int   — optimal hw
        ``e_star``    : int   — optimal ew
        ``auroc_star``: float — best AUROC observed
        ``trace``     : list of (hw, ew, auroc) tuples in evaluation order
        ``gp``        : fitted GP surrogate after final iteration
    """
    if search_space is None:
        search_space = {"hw": (5, 35), "ew": (30, 300)}
    if budget is None:
        budget = (8, 30)

    n_init, n_iter = budget
    rng = np.random.default_rng(seed)

    bounds = np.array([
        [search_space["hw"][0], search_space["hw"][1]],
        [search_space["ew"][0], search_space["ew"][1]],
    ], dtype=float)

    # ── Step 1: Latin Hypercube initial design ──────────────────────────
    X_init = _lhs(n_init, bounds, rng)
    X_init = np.rint(X_init).astype(int)
    X_init = np.clip(X_init, bounds[:, 0].astype(int), bounds[:, 1].astype(int))

    # De-duplicate initial points
    seen = set()
    unique_rows = []
    for row in X_init:
        key = (int(row[0]), int(row[1]))
        if key not in seen:
            seen.add(key)
            unique_rows.append(row)
    X_init = np.array(unique_rows)

    trace: List[Tuple[int, int, float]] = []

    log.info("[%s] BO initial design: %d points", target, len(X_init))

    Y_init = np.empty(len(X_init))
    for i, (hw, ew) in enumerate(X_init):
        hw_int, ew_int = int(hw), int(ew)
        auroc = eval_fn(hw_int, ew_int)
        Y_init[i] = auroc
        trace.append((hw_int, ew_int, auroc))
        log.info("[%s] init %d/%d: hw=%d ew=%d -> AUROC=%.4f",
                 target, i + 1, len(X_init), hw_int, ew_int, auroc)

    X_all = X_init.astype(float)
    Y_all = Y_init.copy()

    # ── Step 2-3: BO iterations ─────────────────────────────────────────
    for it in range(1, n_iter + 1):
        gp = _build_gp(n_dims=2, random_state=seed + it)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gp.fit(X_all, Y_all)

        y_best = float(Y_all.max())
        candidate = _maximise_ei(
            gp, y_best, bounds,
            evaluated=X_all.astype(int),
            xi=xi, n_restarts=n_restarts, rng=rng,
        )

        hw_int, ew_int = int(candidate[0]), int(candidate[1])
        auroc = eval_fn(hw_int, ew_int)
        trace.append((hw_int, ew_int, auroc))

        X_all = np.vstack([X_all, candidate.astype(float).reshape(1, -1)])
        Y_all = np.append(Y_all, auroc)

        improved = " *" if auroc > y_best else ""
        log.info("[%s] iter %d/%d: hw=%d ew=%d -> AUROC=%.4f (best=%.4f)%s",
                 target, it, n_iter, hw_int, ew_int, auroc,
                 max(y_best, auroc), improved)

    # ── Step 4: Return best ─────────────────────────────────────────────
    best_idx = int(np.argmax(Y_all))
    n_star = int(X_all[best_idx, 0])
    e_star = int(X_all[best_idx, 1])
    auroc_star = float(Y_all[best_idx])

    # Final GP fit for the returned surrogate
    gp_final = _build_gp(n_dims=2, random_state=seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gp_final.fit(X_all, Y_all)

    log.info("[%s] BO complete: hw*=%d ew*=%d AUROC*=%.4f (%d evaluations)",
             target, n_star, e_star, auroc_star, len(trace))

    return {
        "n_star": n_star,
        "e_star": e_star,
        "auroc_star": auroc_star,
        "trace": trace,
        "gp": gp_final,
    }
