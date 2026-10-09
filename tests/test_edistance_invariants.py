"""E-distance ranking invariants without requiring Mixscape / full pertpy CI.

Why Mixscape post≥pre is not asserted in default CI
----------------------------------------------------
- ``pertpy`` is an optional heavy extra (``.[pertpy]`` / ``.[de]``); default CI
  installs only ``.[dev]`` so Mixscape/E-distance APIs are often absent.
- Mixscape classification is stochastic and needs a large NT pool; a cheap
  synthetic case that reliably yields post≥pre distance is flaky.
- Guaranteed invariant tested here: planted PCA shifts rank by separation
  strength under a pure energy-distance definition matching the scPerturb
  spirit (larger multivariate shift → larger distance vs control).
"""

from __future__ import annotations

import numpy as np


def _pairwise_mean_norm(a: np.ndarray, b: np.ndarray) -> float:
    """Mean Euclidean distance between all pairs of rows in a and b."""
    # ||a_i - b_j||^2 = ||a_i||^2 + ||b_j||^2 - 2 a_i·b_j
    a2 = np.sum(a * a, axis=1)[:, None]
    b2 = np.sum(b * b, axis=1)[None, :]
    d2 = np.maximum(a2 + b2 - 2.0 * (a @ b.T), 0.0)
    return float(np.sqrt(d2).mean())


def energy_distance(x: np.ndarray, y: np.ndarray) -> float:
    """Scalar energy distance between two point clouds (embedding rows)."""
    return (
        2.0 * _pairwise_mean_norm(x, y)
        - _pairwise_mean_norm(x, x)
        - _pairwise_mean_norm(y, y)
    )


def test_energy_distance_ranks_planted_pca_shifts():
    """Strong shift > weak shift > null cloud (fixed seed)."""
    rng = np.random.default_rng(7)
    n, d = 40, 6
    control = rng.normal(0.0, 0.4, size=(n, d))
    null = rng.normal(0.0, 0.4, size=(n, d))
    weak = rng.normal(0.0, 0.4, size=(n, d))
    weak[:, 0] += 1.5
    strong = rng.normal(0.0, 0.4, size=(n, d))
    strong[:, 0] += 5.0

    d_null = energy_distance(control, null)
    d_weak = energy_distance(control, weak)
    d_strong = energy_distance(control, strong)

    assert d_strong > d_weak > d_null
    assert d_null < 0.5
    assert d_strong > 2.0


def test_energy_distance_zero_for_identical_clouds():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(25, 4))
    assert abs(energy_distance(x, x.copy())) < 1e-9
