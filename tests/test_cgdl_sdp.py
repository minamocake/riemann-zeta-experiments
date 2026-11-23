"""Research B2-mini: small CGdL-style SDP, formulations A (two hard normalisations, fixed width)
and B (one normalisation, dilation-invariant Z), basis degree m = 0 (deg_q_max = 1) and m = 1
(deg_q_max = 3). Solver outputs are numerical candidates; nothing here is a certificate.
"""
import math

import numpy as np
import pytest

from conftest import requires_cvxpy
from riemann_research import cgdl_sdp as sdp

Z_GAUSS = 1.4525472652249223          # B1a: (1 - y^2) e^{-pi y^2} normalised by dilation
TOL = 1e-8


@pytest.fixture(scope="module")
def runs():
    return {(F, m): sdp.solve(F, m) for F in "AB" for m in (0, 1)}


def test_degree_bookkeeping():
    assert [sdp.deg_q_max(m) for m in range(4)] == [1, 3, 5, 7]
    with pytest.raises(ValueError):
        sdp.deg_q_max(-1)


@requires_cvxpy
def test_dimensions_of_solutions(runs):
    for (_, m), r in runs.items():
        assert len(r.q) == sdp.deg_q_max(m) + 1 == len(r.tq)
        for name in ("X2", "X3", "X4"):
            assert np.asarray(r.grams[name]).shape == (m + 1, m + 1)


def test_known_point_m0_exactly_feasible_and_benchmark():
    """m = 0: the only feasible point of A is q = 1 - t, X3 = 1 - 1/(2 pi), X4 = 1 (numpy re-check)."""
    ev = sdp.evaluate_point("A", 0, [[1]], [[1 - 1 / (2 * math.pi)]], [[1]])
    assert ev["q"] == [1, -1]
    assert abs(ev["fourier_match_max_abs"]) < 1e-15 and abs(ev["a0_tq0_minus_1"]) < 1e-15
    assert abs(ev["objective"] - Z_GAUSS) < 1e-14


def test_lower_degree_point_embeds_in_m1():
    """deg_q_max is a maximum: the m = 0 point (padded with zeros) is feasible for m = 1."""
    ev = sdp.evaluate_point("A", 1, np.diag([1, 0]), np.diag([1 - 1 / (2 * math.pi), 0]), np.diag([1, 0]))
    assert ev["q"] == [1, -1, 0, 0]
    assert abs(ev["fourier_match_max_abs"]) < 1e-15 and abs(ev["q0_minus_1"]) < 1e-15
    assert abs(ev["objective"] - Z_GAUSS) < 1e-14


@requires_cvxpy
def test_m0_solutions_reproduce_the_gaussian(runs):
    a, b = runs["A", 0], runs["B", 0]
    assert a.status == b.status == "optimal"
    assert np.allclose(a.q, [1, -1], atol=1e-9)
    assert np.allclose(a.grams["X3"], [[1 - 1 / (2 * math.pi)]], atol=1e-9)
    assert abs(a.objective - Z_GAUSS) < TOL and abs(b.objective - Z_GAUSS) < TOL


@requires_cvxpy
def test_residuals_normalisation_and_psd(runs):
    for r in runs.values():
        assert r.numerically_acceptable, r.notes
        assert abs(r.primal_residuals["fourier_match_max_abs"]) < 1e-9
        assert r.primal_residuals["min_gram_eigenvalue"] > -1e-9
        for k, v in r.primal_residuals.items():
            if k.endswith("minus_1"):
                assert abs(v) < 1e-9, k
        for name in ("X2", "X3", "X4"):
            assert r.dual_psd_min_eigenvalues[name] > -1e-8
            assert abs(r.complementarity[name]) < 1e-8


@requires_cvxpy
def test_numpy_reimplementation_agrees_with_solver(runs):
    for (F, m), r in runs.items():
        ev = sdp.evaluate_point(F, m, r.grams["X2"], r.grams["X3"], r.grams["X4"])
        assert np.allclose(ev["q"], r.q, atol=1e-12)
        assert abs(ev["objective"] - r.objective) < 1e-9


@requires_cvxpy
def test_independent_objective_and_sign_checks(runs):
    for r in runs.values():
        ind = sdp.independent_checks(r)
        assert ind["sign_radius_admissible"]
        assert abs(ind["r_minus_rho"]) < 1e-8          # last sign change is at rho
        assert ind["z_routes_abs_diff"] < 1e-12          # closed form versus quadrature
        assert ind["z_vs_solver_abs_diff"] < TOL
        assert ind["fourier_sign"]["passed"] in (True, None)    # never a found violation
        if ind["fourier_sign"]["passed"] is None:
            assert ind["fourier_sign"]["rigor"] == "exploratory-float"


@requires_cvxpy
def test_map_a_to_b_preserves_objective_and_constraints(runs):
    """q_B = a0 q_A: (T q_B)(0) = 1, q_B(0) = a0, Z_inv(q_B) = Z_A."""
    for m in (0, 1):
        a = runs["A", m]
        qb = sdp.map_a_to_b(a.q)
        tq = sdp.transform_matrix(len(qb)).T @ np.array(qb)
        assert abs(tq[0] - 1) < 1e-9
        assert abs(qb[0] - sdp.A0) < 1e-9
        assert abs(sdp.objective_b(qb) - a.objective) < 1e-9


@requires_cvxpy
def test_value_ordering(runs):
    for m in (0, 1):
        assert runs["B", m].objective <= runs["A", m].objective + TOL
    assert runs["A", 1].objective <= Z_GAUSS + TOL
    assert runs["B", 1].objective <= runs["B", 0].objective + TOL


@requires_cvxpy
def test_characterization_values_m1(runs):
    """Characterization (recorded from this run; paired with the independent checks above)."""
    assert abs(runs["A", 1].objective - 1.4525472652591453) < 1e-7
    assert abs(runs["B", 1].objective - 1.328646166101496) < 1e-7


@requires_cvxpy
def test_infeasible_instance_is_reported_without_a_value():
    r = sdp.solve("A", 0, a0=2.0)   # m = 0 forces (Tq)(0) = 1 - 1/(2 pi), so 2 (Tq)(0) != 1
    assert r.status in ("infeasible", "infeasible_inaccurate")
    assert r.objective is None and not r.numerically_acceptable and r.notes
    with pytest.raises(ValueError):
        sdp.independent_checks(r)


@requires_cvxpy
def test_numerical_uncertainty_is_flagged():
    r = sdp.solve("B", 1, tol=1e-30)
    assert r.objective is not None and not r.numerically_acceptable
    assert any("numerically uncertain" in n for n in r.notes)


def test_invalid_arguments():
    with pytest.raises(ValueError):
        sdp.solve("C", 0)
    with pytest.raises(ValueError):
        sdp.solve("A", 0, rho=1.5)
    with pytest.raises(ValueError):
        sdp.evaluate_point("B", 1, [[1]], [[1]], [[1]])


@requires_cvxpy
def test_duals_satisfy_the_homogeneity_identity(runs):
    """Scaling every normalisation right-hand side by c scales V by c, so dV/dc = V at c = 1.

    With CVXPY's sign convention this reads: -(sum of normalisation duals) = V. A check of the
    dual sign convention needed by C0 (not an adjoint experiment).
    """
    for (F, m), r in runs.items():
        names = ("q0", "a0_tq0") if F == "A" else ("tq0",)
        assert abs(-sum(r.duals[n] for n in names) - r.objective) < 1e-7, (F, m)
