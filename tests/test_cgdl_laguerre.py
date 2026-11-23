"""Research B2-full S1: formulation B in the Laguerre basis (cgdl_laguerre).

Solver outputs are numerical candidates; nothing here is a certificate. Oracles: tests/ORACLES.md.
"""
import json
import os
from fractions import Fraction
from pathlib import Path

import mpmath
import numpy as np
import pytest

from conftest import requires_cvxpy
from riemann_research import cgdl_fourier as cf
from riemann_research import cgdl_laguerre as L

Z_B2MINI_M1 = 1.328646166101496       # B2-mini, formulation B, m = 1, rho = 1 (monomial basis)
Z_GAUSS = 1.4525472652249223
SWEEP = Path(__file__).resolve().parents[1] / "experiments/cgdl_reproduction/b2_full/results/s1_sweep.json"


def w_as_even_poly(n):
    """w_n(2 pi x^2) = L_n^(-1/2)(2 pi t) as an exact B1a EvenPoly in powers of t = x^2."""
    return cf.even_poly(cf.PiLaurent.pi_power(j, c * 2 ** j) for j, c in enumerate(L.laguerre_y(n)))


def test_transform_is_diagonal_on_w_basis_exactly():
    """exact-symbolic: B1a's T (Laguerre monomial formula) maps w_n to (-1)^n w_n."""
    for n in range(16):
        w = w_as_even_poly(n)
        assert cf.transform_laguerre(w) == cf.poly_scale(w, (-1) ** n), n


def test_multiply_by_y_matrix_is_exact():
    """exact-rational: y * w_n reproduced by the three-term recurrence used in the model."""
    N = 12
    Y = L.multiply_by_y_matrix(N)
    for n in range(N):
        lhs = L.to_w_basis([Fraction(0)] + list(L.laguerre_y(n)), N + 1)
        rhs = [Fraction(Y[k, n]) for k in range(N + 1)]
        assert lhs == rhs, n


def test_w_basis_round_trip_and_paper_basis():
    p = [Fraction(3, 7), Fraction(-2), Fraction(5, 11), Fraction(1, 13)]
    coords = L.to_w_basis(p)
    back = [Fraction(0)] * len(p)
    for n, c in enumerate(coords):
        for j, wj in enumerate(L.laguerre_y(n)):
            back[j] += c * wj
    assert back == p
    # paper basis: v_k(y) = L_k^(-1/2)(y/2), i.e. L_k^(-1/2)(pi t)
    assert L.gram_basis_y("paper", 2) == tuple(c * Fraction(1, 2) ** j for j, c in enumerate(L.laguerre_y(2)))


def test_product_tensor_matches_exact_products():
    for basis in ("paper", "hermite"):
        m = 3
        P = L.product_tensor(basis, m)
        for i in range(m + 1):
            for j in range(m + 1):
                exact = L.to_w_basis(L._pmul(list(L.gram_basis_y(basis, i)), list(L.gram_basis_y(basis, j))), 2 * m + 1)
                assert np.allclose(P[i, j], [float(c) for c in exact], rtol=1e-15, atol=0)


def test_objective_integrals_against_quadrature():
    """numerical-high-precision: J_n(rho) = int_0^rho x w_n(2 pi x^2) e^(-pi x^2) dx."""
    for rho in (0.9, 1.2):
        J = L.objective_integrals(21, rho)
        with mpmath.workdps(40):
            for n in (0, 1, 5, 12, 20):
                w = [mpmath.mpf(c.numerator) / c.denominator for c in L.laguerre_y(n)]

                def f(x, w=w):
                    return x * mpmath.polyval(w, 2 * mpmath.pi * x * x, asc=True) * mpmath.exp(-mpmath.pi * x * x)

                q = mpmath.quad(f, mpmath.linspace(0, rho, 9))
                assert abs(J[n] - float(q)) <= 1e-13 * max(1.0, abs(float(q))), (rho, n)


@requires_cvxpy
def test_m0_and_m1_reproduce_earlier_values():
    for basis in ("paper", "hermite"):
        r0 = L.solve(0, 1.0, basis)
        assert L.acceptable(r0) and abs(r0.objective - Z_GAUSS) < 1e-10
        r1 = L.solve(1, 1.0, basis)
        assert L.acceptable(r1) and abs(r1.objective - Z_B2MINI_M1) < 1e-8, (basis, r1.objective)


@requires_cvxpy
def test_bases_agree_and_values_decrease_with_degree():
    vals = []
    for m in (1, 2, 3, 4, 5):
        rp = L.solve(m, 1.0, "paper")
        assert L.acceptable(rp), (m, rp.notes)
        vals.append(rp.objective)
        if m <= 3:
            rh = L.solve(m, 1.0, "hermite")
            assert abs(rh.objective - rp.objective) < 1e-8, m
    assert all(b <= a + 1e-9 for a, b in zip(vals[:-1], vals[1:], strict=True))       # nested families


@requires_cvxpy
def test_independent_checks_and_physical_radius():
    r = L.solve(3, 1.0, "paper")
    ind = L.independent_check(r)
    assert ind["z_routes_abs_diff"] < 1e-12 and ind["z_vs_solver_abs_diff"] < 1e-8
    assert abs(ind["tq0_b1a"] - 1) < 1e-8 and abs(ind["q_at_rho_sq"]) < 1e-8
    hp = L.highprec_check(r)
    assert hp["matching_defect_fn_max"] < 1e-10 and hp["z_vs_solver_abs_diff"] < 1e-8
    assert abs(hp["R_physical"] - r.rho * r.q0 / r.tq0) < 1e-8
    assert r.R_physical >= 1                                          # Cohn-Elkies diagnostic


@requires_cvxpy
def test_kkt_diagnostics_and_settings_are_recorded():
    r = L.solve(2, 1.0, "paper")
    assert r.kkt_stationarity_max_abs < 1e-7 and abs(r.duality_gap) < 1e-7
    assert any(n.startswith("settings ") for n in r.notes)


@pytest.mark.skipif(not SWEEP.exists(), reason="raw sweep output not generated yet (run run_s1.py)")
def test_saved_sweep_is_monotone_and_documented():
    d = json.loads(SWEEP.read_text(encoding="utf-8"))
    acc = [(e["m"], e["best"]["objective"]) for e in d["runs"] if e["best_acceptable_1e-6"]]
    assert [m for m, _ in acc][:5] == [1, 2, 3, 5, 8]
    vals = [v for _, v in acc]
    assert all(b <= a + 1e-8 for a, b in zip(vals[:-1], vals[1:], strict=True))
    assert all(v > 1.3207 for v in vals)
    for e in d["runs"]:
        if e["best"]["objective"] is not None and e["highprec_check"]:
            assert abs(e["highprec_check"]["R_physical"] - e["best"]["R_physical"]) < 1e-6


@pytest.mark.skipif(not os.environ.get("CGDL_Z40"), reason="set CGDL_Z40 to the path of the published Z-40.txt")
def test_model_accepts_published_witness():
    import subprocess
    import sys
    script = Path(__file__).resolve().parents[1] / "experiments/cgdl_reproduction/b2_full/check_published_witness.py"
    out = json.loads(subprocess.run([sys.executable, str(script), os.environ["CGDL_Z40"]], capture_output=True,
                                    text=True, check=True).stdout)
    assert out["max_abs_matching_residual"] < 1e-12 and abs(out["tq0_minus_1"]) < 1e-12
    assert abs(out["objective_minus_published_1.3207644973215569"]) < 1e-12
