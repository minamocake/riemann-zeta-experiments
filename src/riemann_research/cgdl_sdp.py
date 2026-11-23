"""Research B2-mini: a small CGdL-style semidefinite program on the Gaussian-polynomial family.

Exploratory numerical optimisation. Solver outputs are CANDIDATES: nothing here certifies
feasibility, Fourier non-negativity or a bound on zeta zeros. Every use for zeta zeros assumes RH.

Variables and degree bookkeeping
--------------------------------
Dimensionless variable y, t = y^2, and p(y) = q(t). With a monomial basis v(t) = [1, t, ..., t^m]
and Gram blocks X_2, X_3, X_4 of size (m + 1):

    s_i(t) = v(t)^T X_i v(t)                    (degree <= 2m in t)
    q(t)   = (rho^2 - t) s_2(t)                 (degree <= deg_q_max = 2m + 1; X_1 = 0, CGdL 4.1.1)
    (Tq)(t) = s_3(t) + t s_4(t)                 (degree <= 2m + 1; T preserves the degree in t)

deg_q_max is a maximum: lower-degree solutions (e.g. X_2 = diag(1, 0, ...)) stay feasible.
q <= 0 for t >= rho^2 holds by construction once X_2 is PSD; Fourier non-negativity holds once
X_3, X_4 are PSD AND the coefficient-matching identity T q = s_3 + t s_4 holds exactly. The
solver satisfies both only up to tolerances, so neither is certified here (that is B4).

CGdL write p of degree 2d with Gram blocks of size d + 1 (arXiv:1810.08843v2, Section 4.1); read
literally, (x^2 - R^2) s_2(x^2) can have degree 4d + 2 in x. This module does not use CGdL's d:
it uses m and deg_q_max as above. The discrepancy in the paper's notation is not resolved here.

Two formulations
----------------
(A) Fixed width and both hard normalisations (rho = 1, so R = a0 in x):
        q(0) = 1,  a0 (Tq)(0) = 1,   minimise  Z_A = a0 [q(0) + 2 sum_k q_k M_k(1)],
    with a0 = 1 / (1 - 1/(2 pi)) by default and M_k(r) = int_0^r y^(2k+1) exp(-pi y^2) dy.
    On the feasible set q(0) = 1, so Z_A equals a0 [1 + 2 int_0^1 y p(y) exp(-pi y^2) dy].
(B) One normalisation and the dilation-invariant functional:
        (Tq)(0) = 1,  q(0) free,   minimise  Z_inv = rho q(0) + (2 / rho) sum_k q_k M_k(rho).
Map A -> B: q_B = a0 q_A satisfies (T q_B)(0) = 1 and q_B(0) = a0, and Z_inv(q_B) = Z_A. So
(A) is (B) plus the extra constraint q(0) = a0 (at rho = 1), and V_B <= V_A. The dual multiplier
of that extra constraint is a LOCAL sensitivity (derivative of V_A with respect to its right-hand
side at the current value); it does not by itself account for the finite gap V_A - V_B, which
can involve changes of the active set (zero eigenvalues of X_i, contact points) along the way.
In formulation B, q(0) = rho^2 X_2[0, 0] >= 0 follows from PSD, but normalisation needs
q(0) > 0 strictly; this is checked explicitly (`q0_positive`).

The meaningful sensitivity parameter is rho = R / (Gaussian width): Z_inv is invariant under a
single common dilation, which moves R and the width together. CGdL search R >= 1 (rho >= 1), so
rho = 1 is a boundary point and only one-sided sensitivities are meaningful there; V(rho) may
be non-differentiable where the active set (zero eigenvalues of X_i) changes.

Linearity remark (recorded for C0): f_lam = (1 - x^2)(1 + lam x^2) e^(-pi x^2) satisfies
f_lam + f_(-lam) = 2 f_0, and the same identity after Fourier transformation. This is linearity
only; it does not make both perturbations feasible (lam = +0.5 violates the Fourier condition,
lam = -0.5 the real-space condition; see test_cgdl_fourier.py).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

import mpmath
import numpy as np

from . import cgdl_fourier as cf

A0 = 1.0 / (1.0 - 1.0 / (2.0 * np.pi))
POSITIVITY_MARGIN = 1e-8      # q(0) must exceed this for a formulation-B candidate to count as normalisable
A0_INV_EXACT = 1 - cf.PiLaurent.pi_power(-1, Fraction(1, 2))     # 1/a0 = 1 - 1/(2 pi), exactly


def deg_q_max(m: int) -> int:
    """Maximum degree of q in t for basis degree m (X_1 = 0 construction)."""
    if m < 0:
        raise ValueError("m must be >= 0")
    return 2 * m + 1


def transform_matrix(n: int) -> np.ndarray:
    """T as a float matrix: (T q)_j = sum_k M[k, j] q_k for degrees 0..n-1 in t (from B1a)."""
    M = np.zeros((n, n))
    for k in range(n):
        for j, c in enumerate(cf.laguerre_monomial_transform(k)):
            M[k, j] = float(c)
    return M


def moments(n: int, r: float, dps: int = 40) -> np.ndarray:
    """M_k(r) = int_0^r y^(2k+1) exp(-pi y^2) dy, k < n, via B1a's stable gaussian_moment."""
    with mpmath.workdps(dps):
        return np.array([float(cf.gaussian_moment(k, mpmath.mpf(r))) for k in range(n)])


def _antidiag(X, k: int, m: int):
    lo, hi = max(0, k - m), min(k, m)
    terms = [X[i, k - i] for i in range(lo, hi + 1)]
    return sum(terms[1:], terms[0]) if terms else 0


@dataclass
class B2Result:
    """Solver output and diagnostics for one instance (numerical candidate only)."""

    formulation: str
    m: int
    rho: float
    a0: float
    status: str
    solver: str
    solver_version: str
    objective: float | None = None
    q: list[float] = field(default_factory=list)
    tq: list[float] = field(default_factory=list)
    grams: dict[str, list[list[float]]] = field(default_factory=dict)
    gram_eigenvalues: dict[str, list[float]] = field(default_factory=dict)
    primal_residuals: dict[str, float] = field(default_factory=dict)
    duals: dict[str, object] = field(default_factory=dict)
    dual_psd_min_eigenvalues: dict[str, float] = field(default_factory=dict)
    complementarity: dict[str, float] = field(default_factory=dict)
    kkt: dict[str, float] = field(default_factory=dict)
    q0_positive: bool | None = None
    numerically_acceptable: bool = False
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return dict(self.__dict__)


def solve(formulation: str, m: int, rho: float = 1.0, a0: float = A0, tol: float = 1e-7,
          solver_tol: float = 1e-10) -> B2Result:
    """Build and solve formulation 'A' or 'B' with basis degree m. Requires cvxpy + Clarabel."""
    import clarabel
    import cvxpy as cp

    if formulation not in ("A", "B"):
        raise ValueError("formulation must be 'A' or 'B'")
    if formulation == "A" and rho != 1.0:
        raise ValueError("formulation A is defined at rho = 1 (R = a0 in x)")
    if rho <= 0:
        raise ValueError("rho must be positive")
    n, N = m + 1, deg_q_max(m) + 1
    X = {name: cp.Variable((n, n), symmetric=True, name=name) for name in ("X2", "X3", "X4")}
    s2 = [_antidiag(X["X2"], k, m) for k in range(2 * m + 1)]
    s3 = [_antidiag(X["X3"], k, m) for k in range(2 * m + 1)]
    s4 = [_antidiag(X["X4"], k, m) for k in range(2 * m + 1)]
    q = [rho ** 2 * (s2[k] if k <= 2 * m else 0) - (s2[k - 1] if k >= 1 else 0) for k in range(N)]
    T = transform_matrix(N)
    tq = [sum(T[k, j] * q[k] for k in range(N)) for j in range(N)]
    rhs = [(s3[j] if j <= 2 * m else 0) + (s4[j - 1] if j >= 1 else 0) for j in range(N)]
    Mk = moments(N, rho)

    cons = {name: X[name] >> 0 for name in X}
    match = cp.hstack([tq[j] - rhs[j] for j in range(N)])
    cons["fourier_match"] = match == 0
    if formulation == "A":
        cons["q0"] = q[0] == 1
        cons["a0_tq0"] = a0 * tq[0] == 1
        objective = a0 * (q[0] + 2 * sum(Mk[k] * q[k] for k in range(N)))
    else:
        cons["tq0"] = tq[0] == 1
        objective = rho * q[0] + (2 / rho) * sum(Mk[k] * q[k] for k in range(N))

    prob = cp.Problem(cp.Minimize(objective), list(cons.values()))
    res = B2Result(formulation, m, rho, a0, "not-run", "CLARABEL", clarabel.__version__)
    try:
        prob.solve(solver=cp.CLARABEL, tol_gap_abs=solver_tol, tol_gap_rel=solver_tol, tol_feas=solver_tol)
    except cp.error.SolverError as exc:  # pragma: no cover - depends on solver internals
        res.status = "solver-error"
        res.notes.append(f"solver error: {exc}")
        return res
    res.status = str(prob.status)
    if prob.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        res.notes.append("no solution: objective and variables are not reported")
        return res
    if prob.status == cp.OPTIMAL_INACCURATE:
        res.notes.append("solver reported OPTIMAL_INACCURATE")

    res.objective = float(prob.value)
    qv = np.array([float(np.asarray(getattr(e, "value", e))) for e in q])
    res.q = qv.tolist()
    res.tq = (T.T @ qv).tolist()
    res.grams = {k: X[k].value.tolist() for k in X}
    res.gram_eigenvalues = {k: np.linalg.eigvalsh(X[k].value).tolist() for k in X}

    # Primal residuals recomputed from the returned variables (not taken from the solver).
    rhs_v = np.array([float(np.asarray(getattr(e, "value", e))) for e in rhs])
    res.primal_residuals["fourier_match_max_abs"] = float(np.max(np.abs(T.T @ qv - rhs_v)))
    if formulation == "A":
        res.primal_residuals["q0_minus_1"] = float(qv[0] - 1)
        res.primal_residuals["a0_tq0_minus_1"] = float(a0 * (T.T @ qv)[0] - 1)
    else:
        res.primal_residuals["tq0_minus_1"] = float((T.T @ qv)[0] - 1)
    res.primal_residuals["min_gram_eigenvalue"] = float(min(min(v) for v in res.gram_eigenvalues.values()))

    for name, c in cons.items():
        dv = c.dual_value
        res.duals[name] = np.asarray(dv).tolist() if dv is not None else None
        if name in X and dv is not None:
            S = np.asarray(dv)
            res.dual_psd_min_eigenvalues[name] = float(np.linalg.eigvalsh((S + S.T) / 2).min())
            res.complementarity[name] = float(np.trace(S @ X[name].value))

    res.kkt = kkt_diagnostics(formulation, m, rho, a0, res.grams, res.duals, res.objective)
    # q(0) = rho^2 X2[0, 0] >= 0 follows from PSD; normalisation needs q(0) > 0 strictly.
    res.q0_positive = bool(qv[0] > POSITIVITY_MARGIN)
    if formulation == "B" and not res.q0_positive:
        res.notes.append(f"q(0) = {qv[0]:.3e} is not > {POSITIVITY_MARGIN}: candidate is not normalisable")

    worst = max(abs(res.primal_residuals["fourier_match_max_abs"]),
                *(abs(v) for k, v in res.primal_residuals.items() if k.endswith("minus_1")))
    res.numerically_acceptable = (prob.status == cp.OPTIMAL and worst <= tol
                                  and res.primal_residuals["min_gram_eigenvalue"] >= -tol
                                  and res.q0_positive
                                  and abs(res.kkt["stationarity_max_abs"]) <= 100 * tol
                                  and abs(res.kkt["duality_gap"]) <= 100 * tol)
    if not res.numerically_acceptable:
        res.notes.append(f"residuals exceed tol = {tol} or PSD violation; treat as numerically uncertain")
    return res


def _affine_parts(formulation: str, m: int, rho: float, a0: float, X2, X3, X4):
    """Objective value and constraint expressions (lhs - rhs) at given Gram matrices (numpy)."""
    ev = evaluate_point(formulation, m, X2, X3, X4, rho, a0)
    N = deg_q_max(m) + 1
    q = np.array(ev["q"])
    tq = transform_matrix(N).T @ q
    s3 = np.zeros(N)
    n = m + 1
    for i in range(n):
        for j in range(n):
            s3[i + j] += np.asarray(X3)[i, j]
            s3[i + j + 1] += np.asarray(X4)[i, j]
    cons = {"fourier_match": tq - s3}
    if formulation == "A":
        cons["q0"] = np.array([q[0] - 1])
        cons["a0_tq0"] = np.array([a0 * tq[0] - 1])
    else:
        cons["tq0"] = np.array([tq[0] - 1])
    return ev["objective"], cons


def kkt_diagnostics(formulation: str, m: int, rho: float, a0: float, grams: dict, duals: dict,
                    primal_value: float) -> dict:
    """KKT stationarity residual and primal-dual gap, from an explicit numpy model of the SDP.

    The SDP is: minimise <C, X> subject to A_i(X) - b_i = 0 and X_b PSD (blocks b = X2, X3, X4).
    Coefficient matrices are obtained by probing the affine maps with symmetric basis matrices
    (independent of CVXPY's internal canonicalisation). With CVXPY's convention (Lagrangian
    <C, X> + sum_i y_i (A_i(X) - b_i) - sum_b <S_b, X_b>), stationarity reads
    C_b + sum_i y_i A_ib - S_b = 0 and the dual objective is -sum_i y_i b_i. Reported:
    max |stationarity residual| over all blocks, primal and dual objectives, and their gap.
    """
    n = m + 1
    zero = [np.zeros((n, n))] * 3
    names = ("X2", "X3", "X4")
    c0, cons0 = _affine_parts(formulation, m, rho, a0, *zero)       # c0 = 0, cons0 = -b
    C = {b: np.zeros((n, n)) for b in names}
    A = {k: {b: np.zeros((len(v), n, n)) for b in names} for k, v in cons0.items()}
    for bi, b in enumerate(names):
        for j in range(n):
            for k in range(j, n):
                E = np.zeros((n, n))
                E[j, k] = E[k, j] = 1.0
                mats = list(zero)
                mats[bi] = E
                cv, cons = _affine_parts(formulation, m, rho, a0, *mats)
                w = 1.0 if j == k else 0.5          # <M, E> = 2 M_jk off the diagonal
                C[b][j, k] = C[b][k, j] = (cv - c0) * w
                for name in cons:
                    d = (cons[name] - cons0[name]) * w
                    A[name][b][:, j, k] = d
                    A[name][b][:, k, j] = d
    worst = 0.0
    for b in names:
        S = np.asarray(duals[b])
        R = C[b] - S
        for name in cons0:
            y = np.atleast_1d(np.asarray(duals[name], dtype=float))
            R = R + np.tensordot(y, A[name][b], axes=1)
        worst = max(worst, float(np.max(np.abs(R))))
    dual_value = 0.0
    for name, v in cons0.items():             # v = -b
        y = np.atleast_1d(np.asarray(duals[name], dtype=float))
        dual_value += float(y @ v)            # -sum y_i b_i
    return {"stationarity_max_abs": worst, "primal_objective": float(primal_value),
            "dual_objective": dual_value, "duality_gap": float(primal_value) - dual_value}


def map_a_to_b(q_a: list[float], a0: float = A0) -> list[float]:
    """q_B = a0 q_A: maps a formulation-A point to a formulation-B point with the same objective."""
    return [a0 * c for c in q_a]


def objective_b(q: list[float], rho: float = 1.0) -> float:
    """Z_inv = rho q(0) + (2/rho) sum_k q_k M_k(rho) for given coefficients (formulation B)."""
    Mk = moments(len(q), rho)
    return float(rho * q[0] + (2 / rho) * float(np.dot(Mk, q)))


def independent_checks(res: B2Result, dps: int = 40) -> dict:
    """Re-evaluate a solution with B1a routines on the exact binary values of its coefficients.

    Uses the actual sign radius of the returned q (not the prescribed rho), the closed-form and
    quadrature routes for Z_inv, and B1a's Fourier sign check. Numerical-high-precision /
    exploratory only: no certificate.
    """
    if res.objective is None:
        raise ValueError("no solution to check")
    p = cf.even_poly(Fraction(c) for c in res.q)
    inv = A0_INV_EXACT if (res.formulation == "A" and res.a0 == A0) else cf.PiLaurent.const(Fraction(1 / res.a0))
    f = cf.GaussianPolynomial(p, inv if res.formulation == "A" else cf.PiLaurent.const(1))
    sr = cf.sign_radius(p)
    out: dict = {"sign_radius_admissible": sr.admissible, "sign_radius_reason": sr.reason}
    if sr.admissible:
        with mpmath.workdps(dps):
            r = sr.radius_mp()
            out["r_of_q"] = float(r)
            out["r_minus_rho"] = float(r - res.rho)
        out["z_closed_form"] = float(cf.z_closed_form(f, dps))
        out["z_quadrature"] = float(cf.z_quadrature(f, dps))
        out["z_routes_abs_diff"] = abs(out["z_closed_form"] - out["z_quadrature"])
        out["z_vs_solver_abs_diff"] = abs(out["z_closed_form"] - res.objective)
    fc = cf.fourier_sign_check(f)
    out["fourier_sign"] = {"passed": fc.passed, "rigor": fc.rigor, "detail": fc.detail}
    mn = cf.even_poly_minimum(f.fourier_poly(), dps)
    out["fourier_stationary_minimum"] = {
        "unbounded_below": mn["unbounded_below"],
        "min_value": None if mn["min_value"] is None else float(mn["min_value"]),
        "argmin_u": None if mn["argmin_u"] is None else float(mn["argmin_u"]),
        "rigor": mn["rigor"],
        "note": "continuous minimum of (Tq)(u) over u >= 0 from stationary points; not a certificate"}
    out["q0_exact_binary_value_positive"] = Fraction(res.q[0]) > 0
    with mpmath.workdps(dps):
        out["f0"] = float(f.f(0))
        out["fhat0"] = float(f.fhat(0))
    return out


def evaluate_point(formulation: str, m: int, X2, X3, X4, rho: float = 1.0, a0: float = A0) -> dict:
    """Residuals and objective of given Gram matrices, with numpy only (no solver involved).

    Used to check that lower-degree points embed into higher-degree problems and to re-check
    solver output with a second implementation of the same linear maps.
    """
    mats = [np.asarray(Y, dtype=float) for Y in (X2, X3, X4)]
    n = m + 1
    if any(Y.shape != (n, n) for Y in mats):
        raise ValueError(f"Gram matrices must be {n}x{n}")
    N = deg_q_max(m) + 1

    def coeffs(Y):  # s(t) = v^T Y v with v = [1, t, ..., t^m]
        c = np.zeros(2 * m + 1)
        for i in range(n):
            for j in range(n):
                c[i + j] += Y[i, j]
        return c

    s2, s3, s4 = (coeffs(Y) for Y in mats)
    q = np.zeros(N)
    q[: 2 * m + 1] += rho ** 2 * s2
    q[1:] -= s2
    tq = transform_matrix(N).T @ q
    rhs = np.zeros(N)
    rhs[: 2 * m + 1] += s3
    rhs[1:] += s4
    Mk = moments(N, rho)
    out = {"q": q.tolist(), "fourier_match_max_abs": float(np.max(np.abs(tq - rhs))),
           "q0_positive": bool(q[0] > POSITIVITY_MARGIN),
           "min_gram_eigenvalue": float(min(np.linalg.eigvalsh(Y).min() for Y in mats))}
    if formulation == "A":
        out["q0_minus_1"] = float(q[0] - 1)
        out["a0_tq0_minus_1"] = float(a0 * tq[0] - 1)
        out["objective"] = float(a0 * (q[0] + 2 * Mk @ q))
    elif formulation == "B":
        out["tq0_minus_1"] = float(tq[0] - 1)
        out["objective"] = float(rho * q[0] + (2 / rho) * Mk @ q)
    else:
        raise ValueError("formulation must be 'A' or 'B'")
    return out
