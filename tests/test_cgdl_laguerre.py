"""Research B2-full S1: formulation B in the Laguerre basis (cgdl_laguerre).

Solver outputs are numerical candidates; nothing here is a certificate. Oracles: tests/ORACLES.md.
"""
import copy
import json
import os
import subprocess
import sys
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
B2FULL = Path(__file__).resolve().parents[1] / "experiments/cgdl_reproduction/b2_full"
SWEEP = B2FULL / "results/s1_sweep.json"
WITNESSES = B2FULL / "results/s1_witnesses.json"


def _reject_constant(name):
    raise ValueError(f"non-standard JSON constant {name}")


def load_strict(path):
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=_reject_constant)


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
    hp = L.highprec_check(r, fourier_grid_points=4001)
    assert hp["matching_defect_fn_max"] < 1e-10 and hp["z_vs_solver_abs_diff"] < 1e-8
    assert abs(hp["R_physical"] - r.rho * r.q0 / r.tq0) < 1e-8
    assert r.R_physical >= 1                                          # Cohn-Elkies diagnostic
    f = hp["fourier"]
    assert f["reliable"] and abs(f["unweighted_P"]["value_at_t0"] - 1) < 1e-8
    # the weighted minimum is P(t*) exp(-pi t*) at its own stationary point, not the unweighted minimum reweighted
    w = f["weighted_P_exp_minus_pi_t"]
    assert w["infimum"] <= 0 or not w["infimum_attained"]


@requires_cvxpy
def test_kkt_diagnostics_and_settings_are_recorded():
    r = L.solve(2, 1.0, "paper")
    assert r.kkt_stationarity_max_abs < 1e-7 and abs(r.duality_gap) < 1e-7
    assert min(r.dual_min_eig.values()) >= -1e-7 and abs(r.complementarity) < 1e-7
    assert any(n.startswith("settings ") for n in r.notes)
    assert r.solver.endswith("(settings strict)")


# ----------------------------------------------------------------------------- saved sweep (strict JSON, M2/M3)

@pytest.mark.skipif(not SWEEP.exists(), reason="raw sweep output not generated yet (run run_s1.py)")
def test_saved_sweep_is_strict_and_truthful():
    """No benchmark floor: a legitimate lower value must never fail this test (only structure and consistency)."""
    d = load_strict(SWEEP)                           # rejects NaN / Infinity
    assert d["nonfinite_values_replaced_by_null"] == []
    tol = d["acceptance"]["residual_max_abs_max"]
    acc = []
    for e in d["runs"]:
        b = e["best_accepted"]
        if b is None:                                # every solve rejected: nothing reported as a candidate
            assert e["highprec_check"] is None and e["n_accepted"] == 0
            assert all(not ok for _, _, ok, _ in e["rho_evaluations"])
            continue
        r = L.S1Result(**{k: b[k] for k in b})
        assert L.acceptable(r, tol), e["m"]
        assert b["objective"] == min(v for _, v, ok, _ in e["rho_evaluations"] if ok)
        assert abs(e["highprec_check"]["R_physical"] - b["R_physical"]) < 1e-6
        assert abs(e["highprec_check"]["z_vs_solver_abs_diff"]) < 1e-8
        acc.append((e["m"], b["objective"]))
    assert [m for m, _ in acc][:5] == [1, 2, 3, 5, 8]
    # nested families: the true optimum is non-increasing in m; the searched values are checked for it as a
    # consistency diagnostic of the rho search (a failure flags the search, not the mathematics)
    vals = [v for _, v in acc]
    assert all(b <= a + 1e-8 for a, b in zip(vals[:-1], vals[1:], strict=True))


# ----------------------------------------------------------------------------- witnesses and replay (M4)

def _saved_witness(m):
    return copy.deepcopy(next(w for w in load_strict(WITNESSES)["witnesses"] if w["m"] == m))


needs_witnesses = pytest.mark.skipif(not (SWEEP.exists() and WITNESSES.exists()), reason="witness file not generated")


@needs_witnesses
def test_replay_script_is_solver_free_and_passes():
    """Runs replay_s1.py in a fresh interpreter in which importing cvxpy or clarabel fails."""
    code = ("import sys, runpy; sys.modules['cvxpy'] = None; sys.modules['clarabel'] = None; "
            f"sys.argv = ['replay_s1.py', {str(SWEEP)!r}]; runpy.run_path({str(B2FULL / 'replay_s1.py')!r}, "
            "run_name='__main__')")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert p.returncode == 0, p.stdout[-2000:] + p.stderr[-2000:]
    rep = json.loads(p.stdout)
    assert rep["ok"] and rep["witness_sha256_matches"] and len(rep["witnesses"]) >= 2
    assert all(w["replay_ok"] and w["matches_sweep_best_accepted"] and w["q_w_matches_sweep"]
               and w["status_and_tolerance_ok"] for w in rep["witnesses"])


@needs_witnesses
def test_replay_detects_corruption():
    """Each corruption of a saved witness or of its recorded diagnostics makes the replay fail (no exception)."""
    base = _saved_witness(1)
    assert L.replay_witness(base)["ok"]

    def bad(mutate):
        w = copy.deepcopy(base)
        mutate(w)
        out = L.replay_witness(w)
        assert not out["ok"] and out["errors"]
        return out

    def bump(M, i, j, d):
        M[i][j] += d

    bad(lambda w: [bump(w["X"]["X2"], 0, 0, 1e-6)])                                 # primal entry
    out = bad(lambda w: [bump(w["X"]["X3"], 0, 1, 1e-9)])                           # breaks symmetry
    assert not out["checks"]["symmetric"]
    out = bad(lambda w: w["eq_duals"].__setitem__(0, w["eq_duals"][0] + 1e-5))      # equality dual
    assert not out["checks"]["kkt_stationarity"]
    out = bad(lambda w: [bump(w["psd_duals"]["X4"], k, k, -1e-3) for k in range(len(w["psd_duals"]["X4"]))])
    assert not out["checks"]["dual_psd"] and not out["checks"]["kkt_stationarity"]
    bad(lambda w: w["recorded"].__setitem__("objective", w["recorded"]["objective"] + 1e-9))
    bad(lambda w: w["recorded"]["dual_min_eig"].__setitem__("X2", -1.0))
    bad(lambda w: w.__setitem__("rho", w["rho"] * (1 + 1e-9)))
    bad(lambda w: w["recorded"].__setitem__("accepted", False))
    bad(lambda w: [bump(w["X"]["X2"], 0, 0, float("nan"))])
    bad(lambda w: w["X"].pop("X4"))
    bad(lambda w: w["X"].__setitem__("X2", w["X"]["X2"][:-1]))
    bad(lambda w: w.__setitem__("rho", float("inf")))


@needs_witnesses
def test_replay_detects_witness_file_tampering(tmp_path):
    d = load_strict(SWEEP)
    (tmp_path / "s1_sweep.json").write_text(json.dumps(d), encoding="utf-8")
    raw = WITNESSES.read_text(encoding="utf-8")
    (tmp_path / d["witness_file"]["path"]).write_text(raw.replace('"m": 1,', '"m": 1 ,', 1), encoding="utf-8")
    rep = _replay_module().replay(tmp_path / "s1_sweep.json")
    assert not rep["witness_sha256_matches"] and not rep["ok"]


def _replay_module():
    import importlib.util
    spec = importlib.util.spec_from_file_location("replay_s1", B2FULL / "replay_s1.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _replay_pair(tmp_path, degrees=None, mutate_sweep=None, mutate_witnesses=None, fourier=False):
    """Copy the saved sweep and witness file to tmp_path (optionally only some witness degrees), apply the
    mutations, re-record the witness SHA-256 and degrees so that only the targeted field differs, and replay."""
    import hashlib
    sweep, wdoc = load_strict(SWEEP), load_strict(WITNESSES)
    if degrees is not None:
        wdoc["witnesses"] = [w for w in wdoc["witnesses"] if w["m"] in degrees]
    if mutate_witnesses:
        mutate_witnesses(wdoc["witnesses"])
    raw = (json.dumps(wdoc, indent=1, allow_nan=False) + "\n").encode("utf-8")
    (tmp_path / sweep["witness_file"]["path"]).write_bytes(raw)
    sweep["witness_file"]["sha256"] = hashlib.sha256(raw).hexdigest()
    sweep["witness_file"]["degrees"] = [w["m"] for w in wdoc["witnesses"]]
    if mutate_sweep:
        mutate_sweep(sweep, {e["m"]: e for e in sweep["runs"]})
    (tmp_path / "s1_sweep.json").write_text(json.dumps(sweep, allow_nan=False), encoding="utf-8")
    return _replay_module().replay(tmp_path / "s1_sweep.json", fourier)


@needs_witnesses
def test_replay_compares_saved_q_w_status_and_tolerance(tmp_path):
    """Corruptions of the sweep's saved q coefficients, status or acceptance tolerance (or of the witness status
    or tolerance) must fail the replay; the unmodified pair passes."""
    assert _replay_pair(tmp_path)["ok"]

    def scale_first(s, runs):
        runs[1]["best_accepted"]["q_w"][0] *= 1.01

    rep = _replay_pair(tmp_path, mutate_sweep=scale_first)
    e = next(x for x in rep["witnesses"] if x["m"] == 1)
    assert not rep["ok"] and not e["q_w_matches_sweep"]
    assert e["replay_ok"] and e["matches_sweep_best_accepted"] and e["status_and_tolerance_ok"]   # isolated
    assert not _replay_pair(tmp_path, mutate_sweep=lambda s, r: r[1]["best_accepted"]["q_w"].pop())["ok"]
    rep = _replay_pair(tmp_path, mutate_sweep=lambda s, r: r[10]["best_accepted"]["q_w"].__setitem__(3, None))
    assert not rep["ok"]

    def sweep_status(s, runs):
        runs[30]["best_accepted"]["status"] = "optimal_inaccurate"

    rep = _replay_pair(tmp_path, mutate_sweep=sweep_status)
    assert not rep["ok"] and not next(x for x in rep["witnesses"] if x["m"] == 30)["matches_sweep_best_accepted"]

    def witness_status(ws):
        ws[0]["status"] = "optimal_inaccurate"

    rep = _replay_pair(tmp_path, mutate_witnesses=witness_status)
    assert not rep["ok"] and not rep["witnesses"][0]["status_and_tolerance_ok"]

    def witness_tol(ws):
        ws[0]["tol"] = 1e-7

    rep = _replay_pair(tmp_path, mutate_witnesses=witness_tol)
    assert not rep["ok"] and not rep["witnesses"][0]["status_and_tolerance_ok"]

    def sweep_kkt_threshold(s, runs):
        s["acceptance"]["kkt_stationarity_max_abs_max"] = 1e-3

    rep = _replay_pair(tmp_path, mutate_sweep=sweep_kkt_threshold)
    assert not rep["ok"] and not any(x["status_and_tolerance_ok"] for x in rep["witnesses"])


@needs_witnesses
def test_replay_compares_fourier_infimum(tmp_path):
    """With Fourier replay, a changed weighted infimum or attainment flag in the sweep must fail (m = 1 only)."""
    assert _replay_pair(tmp_path, degrees=(1,), fourier=True)["ok"]
    key = "weighted_P_exp_minus_pi_t"

    def infimum(s, runs):
        runs[1]["highprec_check"]["fourier"][key]["infimum"] = -1e-12

    def attained(s, runs):
        f = runs[1]["highprec_check"]["fourier"][key]
        f["infimum_attained"] = not f["infimum_attained"]

    for mut in (infimum, attained):
        rep = _replay_pair(tmp_path, degrees=(1,), mutate_sweep=mut, fourier=True)
        assert not rep["ok"] and not rep["witnesses"][0]["fourier_reproduced"]
        assert rep["witnesses"][0]["q_w_matches_sweep"] and rep["witnesses"][0]["replay_ok"]       # isolated


@requires_cvxpy
def test_witness_round_trip_through_strict_json():
    r = L.solve(2, 1.0, "paper")
    w = json.loads(json.dumps(L.witness_record(r, 1e-7), allow_nan=False))
    out = L.replay_witness(w)
    assert out["ok"], out["errors"]


# ----------------------------------------------------------------------------- acceptance and error paths (M1)

def test_build_model_rejects_non_finite_or_invalid_inputs():
    for m, rho in ((1, float("nan")), (1, float("inf")), (1, -1.0), (1, 0.0), (-1, 1.0), (1.5, 1.0), (True, 1.0)):
        with pytest.raises(ValueError):
            L.build_model(m, rho, "paper")
    with pytest.raises(ValueError):
        L.build_model(1, 1.0, "nope")


def _good_result():
    eig = {b: 1e-12 for b in L.BLOCKS}
    return L.S1Result(2, 1.0, "paper", "optimal", "x", objective=1.3, q0=0.9, tq0=1.0, R_physical=0.9,
                      gram_min_eig=dict(eig), dual_min_eig=dict(eig), residual_max_abs=1e-12,
                      kkt_stationarity_max_abs=1e-12, duality_gap=1e-12)


def test_acceptable_requires_every_diagnostic_finite_and_within_threshold():
    assert L.acceptable(_good_result())
    for field_name, value in (("objective", float("nan")), ("objective", None), ("residual_max_abs", None),
                              ("duality_gap", float("inf")), ("kkt_stationarity_max_abs", None), ("q0", float("nan")),
                              ("tq0", None), ("status", "optimal_inaccurate"), ("status", "solver-error"),
                              ("kkt_stationarity_max_abs", 2e-7), ("residual_max_abs", 2e-7), ("duality_gap", -2e-7),
                              ("q0", 1e-9)):
        r = _good_result()
        setattr(r, field_name, value)
        assert not L.acceptable(r), (field_name, value)
    for attr in ("gram_min_eig", "dual_min_eig"):
        for broken in ({}, {"X2": 0.0}, {"X2": 0.0, "X3": 0.0, "X4": float("nan")},
                       {"X2": 0.0, "X3": 0.0, "X4": -2e-7}):
            r = _good_result()
            setattr(r, attr, broken)
            assert not L.acceptable(r), (attr, broken)
    assert not L.acceptable(_good_result(), float("nan"))


@requires_cvxpy
def test_failed_rung_never_relabels_an_earlier_solution():
    """An OPTIMAL_INACCURATE rung followed by a rung that raises SolverError: the reported point, status and
    settings label are those of the rung that produced the point."""
    ladder = (("tiny-iters", dict(max_iter=8)), ("nan-reg", dict(static_regularization_constant=float("nan"))))
    r = L.solve(3, 1.0, "paper", ladder=ladder)
    assert r.notes[1].startswith("settings nan-reg: solver error")
    assert r.status == "optimal_inaccurate" and r.solver.endswith("(settings tiny-iters)")
    assert not L.acceptable(r)


@requires_cvxpy
def test_all_rungs_failing_gives_solver_error():
    ladder = (("nan-reg", dict(static_regularization_constant=float("nan"))),)
    r = L.solve(2, 1.0, "paper", ladder=ladder)
    assert r.status == "solver-error" and r.objective is None and not r.grams and not L.acceptable(r)


# ----------------------------------------------------------------------------- published witness SHA-256 (M5)

def test_published_witness_checker_rejects_other_files(tmp_path):
    fake = tmp_path / "Z-40.txt"
    fake.write_text("1.0\n[1 0; 0 1]\n[1 0; 0 1]\n[1 0; 0 1]\n", encoding="utf-8")
    p = subprocess.run([sys.executable, str(B2FULL / "check_published_witness.py"), str(fake)],
                       capture_output=True, text=True)
    assert p.returncode == 2 and "does not match the expected" in p.stderr and p.stdout == ""


# ----------------------------------------------------------------------------- Fourier-side minima (M6)

def test_fourier_minima_unweighted_and_weighted_are_distinct():
    """P(t) = (t - 1)^2 + eps: unweighted minimum eps at t = 1; weighted minimum of P(t) e^(-pi t) at the
    root of P' - pi P near 1 (u = t - 1 solves pi u^2 - 2u + pi eps = 0), with a different value."""
    with mpmath.workdps(60):
        eps = mpmath.mpf(10) ** -6
        f = L.fourier_side_minima([1 + eps, -2, 1], dps=60, grid_points=4001)
        pi = mpmath.pi
        u = (1 - mpmath.sqrt(1 - pi ** 2 * eps)) / pi
        g = (u ** 2 + eps) * mpmath.exp(-pi * (1 + u))
    assert f["reliable"] and f["polyroots_converged"]
    unw, wtd = f["unweighted_P"], f["weighted_P_exp_minus_pi_t"]
    assert abs(unw["min_value"] - 1e-6) < 1e-15 and abs(unw["argmin_t"] - 1) < 1e-15
    assert abs(wtd["argmin_t"] - float(1 + u)) < 1e-15 and abs(wtd["min_value_attained"] - float(g)) < 1e-20
    assert wtd["infimum"] == 0.0 and not wtd["infimum_attained"]           # tail -> 0+, not attained


def test_fourier_minima_boundary_and_tail_cases():
    f = L.fourier_side_minima([-1, 1], dps=40, grid_points=2001)           # P = t - 1
    assert f["unweighted_P"]["min_value"] == -1.0 and f["unweighted_P"]["argmin_t"] == 0.0     # boundary t = 0
    w = f["weighted_P_exp_minus_pi_t"]
    assert w["min_value_attained"] == -1.0 and w["argmin_t"] == 0.0 and w["infimum_attained"]
    assert w["n_positive_stationary_points"] == 1                          # t = 1 + 1/pi is a maximum of g
    f = L.fourier_side_minima([1, 0, -1], dps=40, grid_points=2001)        # P = 1 - t^2
    assert f["unweighted_P"]["unbounded_below"] and f["unweighted_P"]["min_value"] is None
    w = f["weighted_P_exp_minus_pi_t"]
    assert w["min_value_attained"] < 0 and "- side" in w["tail"]           # g -> 0- but attains a minimum
    with mpmath.workdps(40):
        t = mpmath.findroot(lambda t: -2 * t - mpmath.pi * (1 - t * t), 1.5)
    assert abs(w["argmin_t"] - float(t)) < 1e-14
    f = L.fourier_side_minima([2], dps=40, grid_points=101)                # constant
    assert f["unweighted_P"]["min_value"] == 2.0 and f["weighted_P_exp_minus_pi_t"]["infimum"] == 0.0


def test_fourier_minima_report_unreliable_instead_of_a_global_claim(monkeypatch):
    """If the root finder misses stationary points, the grid cross-check must flag it."""
    monkeypatch.setattr(L, "_positive_real_roots", lambda a, dps: ([], True, 0.0))
    f = L.fourier_side_minima([1 + 1e-6, -2, 1], dps=40, grid_points=4001)
    assert not f["reliable"] and "global minimum NOT established" in f["limitations"]


@pytest.mark.skipif(not os.environ.get("CGDL_Z40"), reason="set CGDL_Z40 to the path of the published Z-40.txt")
def test_model_accepts_published_witness():
    import subprocess
    import sys
    script = Path(__file__).resolve().parents[1] / "experiments/cgdl_reproduction/b2_full/check_published_witness.py"
    out = json.loads(subprocess.run([sys.executable, str(script), os.environ["CGDL_Z40"]], capture_output=True,
                                    text=True, check=True).stdout)
    assert out["sha256_verified"]
    assert out["file_sha256"] == "a268555389798aa132a46d90057f2dfa28bc1dcf81d212282b057349d2cbd9e9"
    assert out["max_abs_matching_residual"] < 1e-12 and abs(out["tq0_minus_1"]) < 1e-12
    assert abs(out["objective_minus_published_1.3207644973215569"]) < 1e-12
