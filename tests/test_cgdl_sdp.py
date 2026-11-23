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


# --------------------------------------------------------------------------- review of PR #2

def _saved_b_m1_q():
    import json
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "experiments/cgdl_reproduction/b2_mini/results/b2_mini_results.json"
    runs = json.loads(path.read_text(encoding="utf-8"))["runs"]
    return next(e for e in runs if e["result"]["formulation"] == "B" and e["result"]["m"] == 1)["result"]["q"]


def test_fourier_near_contact_from_saved_coefficients():
    """numerical-high-precision: continuous minimum of (Tq)(u) for the saved B, m = 1 candidate.

    About 4.0928e-9 at u = 1.0668272166, far below the sampled grid minimum 4.8554e-7. This is
    neither a positivity certificate nor a proof of an exact double zero.
    """
    from fractions import Fraction

    from riemann_research import cgdl_fourier as cf
    q = _saved_b_m1_q()
    mn = cf.even_poly_minimum(cf.transform_laguerre(cf.even_poly(Fraction(c) for c in q)))
    assert not mn["unbounded_below"] and mn["rigor"] == "numerical-high-precision"
    assert abs(float(mn["min_value"]) - 4.092813391e-9) < 1e-17
    assert abs(float(mn["argmin_u"]) - 1.0668272166) < 1e-9
    assert float(mn["min_value"]) < 4.8554e-7 / 100


@requires_cvxpy
def test_stationary_minimum_reported_for_fresh_runs(runs):
    ind = sdp.independent_checks(runs["B", 1])["fourier_stationary_minimum"]
    assert 0 < ind["min_value"] < 1e-7 and abs(ind["argmin_u"] - 1.0668272) < 1e-5
    ind_a = sdp.independent_checks(runs["A", 1])["fourier_stationary_minimum"]
    assert ind_a["min_value"] > 0.8 and ind_a["argmin_u"] == 0   # A is far from Fourier contact


def test_q0_positivity_is_checked_explicitly():
    """PSD only gives q(0) = rho^2 X2[0,0] >= 0; a zero q(0) must be reported as not normalisable."""
    ev = sdp.evaluate_point("B", 1, np.diag([0.0, 1.0]), np.diag([1.0, 1.0]), np.diag([1.0, 1.0]))
    assert ev["q"][0] == 0 and ev["q0_positive"] is False
    ev = sdp.evaluate_point("B", 0, [[1.0]], [[1.0]], [[1.0]])
    assert ev["q0_positive"] is True


@requires_cvxpy
def test_q0_positive_for_all_candidates(runs):
    for r in runs.values():
        assert r.q0_positive is True
        assert sdp.independent_checks(r)["q0_exact_binary_value_positive"] is True


@requires_cvxpy
def test_kkt_stationarity_and_duality_gap(runs):
    """Explicit KKT stationarity residual and primal-dual gap from a numpy model of the SDP."""
    for (F, m), r in runs.items():
        k = r.kkt
        assert k["stationarity_max_abs"] < 1e-8, (F, m, k)
        assert abs(k["duality_gap"]) < 1e-8, (F, m, k)
        assert abs(k["dual_objective"] - r.objective) < 1e-8


@requires_cvxpy
def test_kkt_diagnostic_detects_wrong_duals(runs):
    """Corrupting a dual multiplier must show up in the stationarity residual."""
    r = runs["B", 1]
    bad = dict(r.duals)
    bad["tq0"] = r.duals["tq0"] * 1.01
    k = sdp.kkt_diagnostics("B", 1, 1.0, sdp.A0, r.grams, bad, r.objective)
    assert k["stationarity_max_abs"] > 1e-4 and abs(k["duality_gap"]) > 1e-3
