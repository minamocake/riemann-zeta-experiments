"""Kernel and window-functional checks against independent routes."""
from fractions import Fraction as Fr

import numpy as np
import pytest

from riemann_research import kernels as K
from riemann_research.montgomery import window_functional


def test_fejer_kernel_gives_montgomery_four_thirds_exactly():
    W = [Fr(1)]
    C = K.kernel_exact(W[0], W[0], Fr(0)) + 2 * K.weighted_integral_exact(W[0], W[0])
    assert C == Fr(4, 3)


@pytest.mark.parametrize("L,M", [(1.0, 1.0), (0.5, 1.5), (0.1, 1.9), (2.5, 0.3), (1e-3, 1.0), (1.0, 3.0)])
def test_float_kernel_matches_exact_rational_kernel(L, M):
    for a in [0.0, abs(L - M) / 2, (L + M) / 2, 0.37, 1.0, 1.73, 3.2]:
        exact = float(K.kernel_exact(Fr(L), Fr(M), Fr(a)))
        assert K.kernel_matrix([L, M], a)[0, 1] == pytest.approx(exact, abs=1e-13)
        assert K.overlap(L, M, a) / (L * M) == pytest.approx(exact, abs=1e-13)


def test_weighted_integral_matches_dense_quadrature():
    rng = np.random.default_rng(0)
    for _ in range(20):
        L, M = rng.uniform(0.05, 2.5, 2)
        a = np.linspace(0, 1, 200001)
        vals = a * np.array([K.overlap(L, M, x) for x in a]) / (L * M)
        assert float(K.weighted_integral_exact(Fr(L), Fr(M))) == pytest.approx(np.trapezoid(vals, a), abs=1e-8)


def test_correlation_matches_brute_force_convolution_of_h():
    W = [0.3, 0.7, 1.0, 1.4, 2.0]
    q = np.array([0.8, -0.5, 1.2, -0.9, 0.4])
    q /= q.sum()
    dt = 2e-5
    t = np.arange(-2.0, 2.0, dt)
    def h(s):
        return sum(qj * (np.abs(s) <= L / 2) / L for qj, L in zip(q, W, strict=True))
    for a in [0.0, 0.65, 1.17, 1.85]:
        brute = np.sum(h(t) * h(t + a)) * dt
        assert K.correlation(W, q, a) == pytest.approx(brute, abs=1e-6)


def test_window_functional_matches_physical_side_integral():
    """C_window = int f(x) mu_hat(x) dx with mu = delta + |a| 1_[-1,1]: a completely different route."""
    W = list(np.linspace(1 / 12, 1, 12))
    q = np.random.default_rng(1).dirichlet(np.ones(12))
    X, dx = 2000.0, 4e-4
    x = np.arange(dx / 2, X, dx)
    H = sum(qj * np.sinc(L * x) for qj, L in zip(q, W, strict=True))
    w = 2 * np.pi * x
    g = 2 * (np.sin(w) / w + (np.cos(w) - 1) / w ** 2)
    tail = sum(qj ** 2 / (np.pi ** 2 * L ** 2 * X) for qj, L in zip(q, W, strict=True))
    C_phys = 2 * np.sum(H * H * (1 + g)) * dx + tail
    assert window_functional(W, q) == pytest.approx(C_phys, abs=2e-6)


def test_nonnegative_mixtures_have_nonincreasing_correlation():
    W = np.linspace(2 / 12, 2, 12)
    a = np.linspace(0, 2.2, 2001)
    rng = np.random.default_rng(2)
    for _ in range(30):
        q = rng.dirichlet(np.ones(12))
        R = np.array([K.correlation(W, q, x) for x in a])
        assert np.max(np.diff(R)) <= 1e-15


def test_knot_evaluation_is_never_beaten_by_a_dense_grid():
    W = np.linspace(2 / 12, 2, 12)
    rng = np.random.default_rng(3)
    for _ in range(20):
        q = rng.normal(size=12)
        q /= q.sum()
        for lo, hi in [(1.0, 1.1), (1.23, 1.71), (1.0, 2.0)]:
            kn = max(K.correlation(W, q, k) for k in K.knots(W, lo, hi))
            grid = max(K.correlation(W, q, x) for x in np.linspace(lo, hi, 1501))
            assert grid <= kn + 1e-14
