"""Numerical checks of the sign-only degeneracy proposition (docs/derivations/sign-only-degeneracy.md)."""
from fractions import Fraction as Fr

import numpy as np

from riemann_research import kernels as K
from conftest import requires_cvxpy

MONTGOMERY_TAYLOR_P = 0.6725  # literature value, rounded


def exact_max_R_on_1_S(W, X):
    n = len(W)
    S = max(W)
    def R(a):
        return sum(X[i][j] * K.kernel_exact(W[i], W[j], a) for i in range(n) for j in range(n) if X[i][j])
    return max(R(k) for k in K.knots_exact(W, Fr(1), S))


def test_any_weight_on_a_box_wider_than_one_makes_R_positive_somewhere_on_1_S():
    rng = np.random.default_rng(0)
    W = [Fr(k, 6) for k in range(1, 13)]
    for _ in range(25):
        n_comp = rng.integers(1, 4)
        X = [[Fr(0)] * 12 for _ in range(12)]
        for _ in range(n_comp):
            v = [Fr(int(x)) for x in rng.integers(-5, 6, size=12)]
            k = int(rng.integers(6, 12))   # force a nonzero entry on some width > 1
            v[k] = Fr(int(rng.integers(1, 6)))
            for i in range(12):
                for j in range(12):
                    X[i][j] += v[i] * v[j]
        assert exact_max_R_on_1_S(W, X) > 0


def test_legacy_sign_only_candidates_put_only_rounding_noise_on_widths_above_one():
    import signed_sensitivity as legacy_signed
    for A in [1.5, 2.0]:
        prog = legacy_signed.make_program(A, 12)
        opt = legacy_signed.optimize_signed(prog, mode="sign_only", restarts=4)
        big = prog.model.widths > 1 + 1e-12
        assert np.max(np.abs(opt["q"][big])) < 1e-9


@requires_cvxpy
def test_best_psd_mixture_of_boxes_of_width_at_most_one_stays_below_montgomery_taylor():
    from riemann_research.optimize import sdp
    for n in [6, 8, 11, 12]:
        W = np.linspace(1 / n, 1, n)
        res = sdp(W, "sign-only")
        assert 2 - res["objective"] < MONTGOMERY_TAYLOR_P + 1e-4
