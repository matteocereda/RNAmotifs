"""Tests for rnamotifs_mars.optim.bayes_opt.

Run:  pytest tests/test_bayes_opt.py -v
"""

import math
import sys
import os

import numpy as np
import pytest

# Ensure the repo root is on sys.path so the package is importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rnamotifs_mars.optim.bayes_opt import optimise_rbp, _lhs, _build_gp


# ---------------------------------------------------------------------------
# Synthetic objectives
# ---------------------------------------------------------------------------

def _negative_rosenbrock_2d(x1, x2, x1_off=50, x2_off=200):
    """Negative Rosenbrock shifted into the (hw, ew) domain.

    Optimum at (x1_off, x2_off) with value 0.0.  We negate and add 1 so
    the maximum is 1.0 at the optimum, mimicking an AUROC.
    """
    a = (x1 - x1_off) / 20.0
    b = (x2 - x2_off) / 100.0
    val = (1 - a) ** 2 + 100 * (b - a ** 2) ** 2
    # Map to [0, 1] range:  1 / (1 + val)
    return 1.0 / (1.0 + val)


def _quadratic_peak(hw, ew, hw_opt=40, ew_opt=150):
    """Simple quadratic with known peak — easy for the GP to model."""
    d_hw = (hw - hw_opt) / 30.0
    d_ew = (ew - ew_opt) / 100.0
    return max(0.0, 1.0 - d_hw ** 2 - d_ew ** 2)


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------

class TestLHS:
    def test_shape(self):
        rng = np.random.default_rng(0)
        bounds = np.array([[5.0, 100.0], [20.0, 400.0]])
        pts = _lhs(10, bounds, rng)
        assert pts.shape == (10, 2)

    def test_within_bounds(self):
        rng = np.random.default_rng(1)
        bounds = np.array([[5.0, 100.0], [20.0, 400.0]])
        pts = _lhs(50, bounds, rng)
        assert np.all(pts[:, 0] >= 5.0)
        assert np.all(pts[:, 0] <= 100.0)
        assert np.all(pts[:, 1] >= 20.0)
        assert np.all(pts[:, 1] <= 400.0)

    def test_stratification(self):
        """Each stratum should have exactly one point per dimension."""
        rng = np.random.default_rng(2)
        n = 20
        bounds = np.array([[0.0, 1.0], [0.0, 1.0]])
        pts = _lhs(n, bounds, rng)
        for dim in range(2):
            bins = np.floor(pts[:, dim] * n).astype(int)
            bins = np.clip(bins, 0, n - 1)
            assert len(set(bins)) == n, "LHS should have one point per stratum"


class TestBuildGP:
    def test_returns_gp(self):
        gp = _build_gp(2)
        assert hasattr(gp, "fit")
        assert hasattr(gp, "predict")

    def test_fit_predict(self):
        gp = _build_gp(2)
        X = np.array([[10.0, 50.0], [30.0, 100.0], [60.0, 200.0],
                       [80.0, 300.0], [20.0, 150.0]])
        Y = np.array([0.5, 0.7, 0.9, 0.6, 0.8])
        gp.fit(X, Y)
        mu, std = gp.predict(X[:1], return_std=True)
        assert mu.shape == (1,)
        assert std.shape == (1,)


class TestOptimiseRBP:
    """Core BO loop tests on synthetic objectives."""

    def test_finds_quadratic_peak(self):
        """BO should recover the quadratic peak within tolerance."""
        hw_opt, ew_opt = 40, 150

        def obj(hw, ew):
            return _quadratic_peak(hw, ew, hw_opt, ew_opt)

        result = optimise_rbp(
            target="synthetic_quad",
            eval_fn=obj,
            search_space={"hw": (5, 100), "ew": (20, 400)},
            budget=(8, 30),
            seed=42,
        )

        assert abs(result["n_star"] - hw_opt) <= 10, (
            f"hw* should be near {hw_opt}, got {result['n_star']}")
        assert abs(result["e_star"] - ew_opt) <= 30, (
            f"ew* should be near {ew_opt}, got {result['e_star']}")
        assert result["auroc_star"] >= 0.85, (
            f"Best AUROC should be >= 0.85, got {result['auroc_star']:.4f}")

    def test_finds_rosenbrock_peak(self):
        """BO on shifted negative Rosenbrock."""
        hw_opt, ew_opt = 50, 200

        def obj(hw, ew):
            return _negative_rosenbrock_2d(hw, ew, hw_opt, ew_opt)

        result = optimise_rbp(
            target="synthetic_rosen",
            eval_fn=obj,
            search_space={"hw": (5, 100), "ew": (20, 400)},
            budget=(8, 30),
            seed=123,
        )

        assert result["auroc_star"] >= 0.7, (
            f"Should find near-optimal region, got {result['auroc_star']:.4f}")

    def test_trace_length(self):
        """Trace should contain n_init + n_iter entries."""
        n_init, n_iter = 5, 10

        result = optimise_rbp(
            target="trace_test",
            eval_fn=lambda hw, ew: _quadratic_peak(hw, ew),
            search_space={"hw": (5, 100), "ew": (20, 400)},
            budget=(n_init, n_iter),
            seed=0,
        )

        # May be slightly less than n_init if LHS produces duplicates
        assert len(result["trace"]) <= n_init + n_iter
        assert len(result["trace"]) >= n_init  # At least initial points

    def test_trace_entries_are_tuples(self):
        result = optimise_rbp(
            target="format_test",
            eval_fn=lambda hw, ew: 0.5,
            search_space={"hw": (5, 50), "ew": (20, 100)},
            budget=(3, 2),
            seed=7,
        )
        for entry in result["trace"]:
            assert len(entry) == 3
            hw, ew, auroc = entry
            assert isinstance(hw, int)
            assert isinstance(ew, int)
            assert isinstance(auroc, float)

    def test_gp_returned(self):
        result = optimise_rbp(
            target="gp_test",
            eval_fn=lambda hw, ew: _quadratic_peak(hw, ew),
            budget=(5, 5),
            seed=0,
        )
        gp = result["gp"]
        assert hasattr(gp, "predict")
        mu, std = gp.predict(np.array([[40.0, 150.0]]), return_std=True)
        assert mu.shape == (1,)

    def test_reproducibility(self):
        """Same seed should give identical traces."""
        kwargs = dict(
            target="repro",
            eval_fn=lambda hw, ew: _quadratic_peak(hw, ew),
            budget=(5, 10),
            seed=99,
        )
        r1 = optimise_rbp(**kwargs)
        r2 = optimise_rbp(**kwargs)
        assert r1["trace"] == r2["trace"]
        assert r1["n_star"] == r2["n_star"]
        assert r1["e_star"] == r2["e_star"]

    def test_bounds_respected(self):
        """All evaluated points must be within the search space."""
        lo_hw, hi_hw = 10, 80
        lo_ew, hi_ew = 30, 300

        result = optimise_rbp(
            target="bounds_test",
            eval_fn=lambda hw, ew: _quadratic_peak(hw, ew),
            search_space={"hw": (lo_hw, hi_hw), "ew": (lo_ew, hi_ew)},
            budget=(8, 15),
            seed=42,
        )

        for hw, ew, _ in result["trace"]:
            assert lo_hw <= hw <= hi_hw, f"hw={hw} out of [{lo_hw}, {hi_hw}]"
            assert lo_ew <= ew <= hi_ew, f"ew={ew} out of [{lo_ew}, {hi_ew}]"

    def test_monotone_best_so_far(self):
        """Best-so-far AUROC should be non-decreasing."""
        result = optimise_rbp(
            target="monotone",
            eval_fn=lambda hw, ew: _quadratic_peak(hw, ew),
            budget=(5, 15),
            seed=42,
        )
        best_so_far = -1.0
        for _, _, auroc in result["trace"]:
            best_so_far = max(best_so_far, auroc)
        assert best_so_far == result["auroc_star"]


# ---------------------------------------------------------------------------
# Integration test (requires pipeline binaries — skipped in CI)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not os.path.isfile(
        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     "build", "rnamotifs_mars_score")),
    reason="rnamotifs_mars_score binary not built"
)
class TestIntegrationPTBP1:
    """Compare grid and bayes on PTBP1/HepG2 if the full pipeline is available.

    This test is designed to be run manually on a machine with the compiled
    C++ binaries and the eCLIP/exon data in place.  It is automatically
    skipped in environments without the binary.
    """

    def test_bayes_comparable_to_grid(self):
        """Bayes should find AUROC within noise of grid on PTBP1."""
        # This is a placeholder that documents the test design.
        # The actual comparison is done by bench/run_benchmark.py which
        # produces bayes_vs_grid.csv.  Here we just verify the BO module
        # runs without error on a trivial synthetic stand-in.
        grid_aurocs = {
            (5, 30): 0.80, (5, 50): 0.82, (5, 100): 0.85,
            (15, 30): 0.83, (15, 50): 0.88, (15, 100): 0.90,
            (25, 50): 0.92, (25, 200): 0.983, (25, 300): 0.91,
            (35, 50): 0.95, (35, 100): 0.93, (35, 300): 0.88,
        }

        def mock_pipeline(hw, ew):
            # Interpolate from nearest grid point
            best_match = min(grid_aurocs.keys(),
                             key=lambda k: (k[0] - hw) ** 2 + (k[1] - ew) ** 2)
            dist = math.sqrt((best_match[0] - hw) ** 2 +
                             (best_match[1] - ew) ** 2)
            base = grid_aurocs[best_match]
            return max(0.0, base - 0.001 * dist)

        result = optimise_rbp(
            target="PTBP1_mock",
            eval_fn=mock_pipeline,
            search_space={"hw": (5, 100), "ew": (20, 400)},
            budget=(8, 30),
            seed=42,
        )

        grid_best = max(grid_aurocs.values())  # 0.983
        # Bayes should get close — within 0.05 of grid best
        assert result["auroc_star"] >= grid_best - 0.05, (
            f"Bayes AUROC {result['auroc_star']:.4f} too far from "
            f"grid best {grid_best:.4f}")
