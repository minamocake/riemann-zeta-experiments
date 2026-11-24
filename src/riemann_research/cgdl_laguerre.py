"""Research B2-full S1: formulation B in a Laguerre (Hermite-function) basis.

Numerical candidates only. Nothing here certifies feasibility, Fourier non-negativity or a bound on
zeta zeros; every use for zeta zeros assumes RH (AGENTS.md).

Variables and bases
-------------------
x is the dimensionless variable of the undilated function f(x) = q(x^2) exp(-pi x^2), t = x^2, and
y = 2 pi t. Output basis: w_n(y) = L_n^(-1/2)(y). Because L_n^(-1/2)(2 pi x^2) exp(-pi x^2) are even
Hermite functions, the Fourier operator T of B1a is diagonal there:

    T w_n = (-1)^n w_n                                   (exact; checked against B1a in the tests).

Gram basis v_k(t) = L_k^(-1/2)(sigma t), k = 0..m, with sigma = 2 pi ("hermite", v_k = w_k) or
sigma = pi ("paper", the basis of CGdL's ZetaBounds.jl, v_k(y) = L_k^(-1/2)(y / 2)). All connection,
product and multiplication coefficients are rational and are computed exactly with Fractions; the
only irrational inputs are 1/(2 pi) (from t = y / (2 pi)), rho and the objective integrals.

Formulation B (as in B2-mini, deg_q_max = 2m + 1):
    q(t) = (rho^2 - t) s_2(t),  (Tq)(t) = s_3(t) + t s_4(t),  s_i = v^T X_i v,  X_i PSD,
    (Tq)(0) = 1,  minimise  Z_inv = rho q(0) + (2 / rho) int_0^rho x q(x^2) exp(-pi x^2) dx.
rho is dimensionless. The normalised physical sign radius is R = rho q(0) / (Tq)(0) (= rho q(0) here),
and an exactly admissible function has R >= 1 (CGdL Section 4: one-dimensional Cohn-Elkies bound).
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache

import mpmath
import numpy as np

from . import cgdl_fourier as cf

ALPHA = Fraction(-1, 2)
BASES = {"hermite": Fraction(1), "paper": Fraction(1, 2)}   # sigma / (2 pi)


# =========================================================================== exact rational algebra in y

def _half_binom(a: Fraction, k: int) -> Fraction:
    out = Fraction(1)
    for i in range(k):
        out = out * (a - i) / (i + 1)
    return out


@lru_cache(maxsize=None)
def laguerre_y(n: int) -> tuple[Fraction, ...]:
    """Coefficients (low degree first) of L_n^(-1/2)(y): (-1)^j binom(n - 1/2, n - j) / j!."""
    return tuple((-1) ** j * _half_binom(n + ALPHA, n - j) / math.factorial(j) for j in range(n + 1))


def w_at_zero(n: int) -> Fraction:
    """L_n^(-1/2)(0) = binom(n - 1/2, n)."""
    return laguerre_y(n)[0]


def _pmul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def to_w_basis(poly: list[Fraction], size: int | None = None) -> list[Fraction]:
    """Exact coordinates of a polynomial in y in the basis {w_n}, by elimination from the top degree."""
    rem = list(poly)
    deg = len(rem) - 1
    coords = [Fraction(0)] * (size if size is not None else deg + 1)
    for n in range(deg, -1, -1):
        if rem[n] == 0:
            continue
        w = laguerre_y(n)
        c = rem[n] / w[n]
        coords[n] = c
        for j in range(n + 1):
            rem[j] -= c * w[j]
    return coords


@lru_cache(maxsize=None)
def gram_basis_y(basis: str, k: int) -> tuple[Fraction, ...]:
    """v_k as a polynomial in y: L_k^(-1/2)(lambda y), lambda = sigma / (2 pi)."""
    lam = BASES[basis]
    return tuple(c * lam ** j for j, c in enumerate(laguerre_y(k)))


@lru_cache(maxsize=None)
def product_tensor(basis: str, m: int) -> np.ndarray:
    """P[i, j, n] = n-th w-coordinate of v_i v_j (exact rationals rounded once to float64)."""
    N = 2 * m + 1
    P = np.zeros((m + 1, m + 1, N))
    for i in range(m + 1):
        for j in range(i, m + 1):
            coords = to_w_basis(_pmul(list(gram_basis_y(basis, i)), list(gram_basis_y(basis, j))), N)
            P[i, j, :] = P[j, i, :] = [float(c) for c in coords]
    return P


def multiply_by_y_matrix(N: int) -> np.ndarray:
    """Y with (y * sum_n a_n w_n) = sum_n (Y a)_n w_n on degrees < N + 1 (output length N + 1).

    Three-term recurrence for alpha = -1/2: y w_n = -(n+1) w_{n+1} + (2n + 1/2) w_n - (n - 1/2) w_{n-1}.
    """
    Y = np.zeros((N + 1, N))
    for n in range(N):
        Y[n + 1, n] = -(n + 1)
        Y[n, n] = 2 * n + 0.5
        if n >= 1:
            Y[n - 1, n] = -(n - 0.5)
    return Y


def objective_integrals(N: int, rho: float, dps: int | None = None) -> np.ndarray:
    """J_n(rho) = int_0^rho x w_n(2 pi x^2) exp(-pi x^2) dx, n < N, via B1a's stable Gaussian moments.

    Expansion in monomials cancels; the working precision grows with n to absorb it.
    """
    dps = dps if dps is not None else 40 + 2 * N
    out = np.zeros(N)
    with mpmath.workdps(dps):
        r = mpmath.mpf(rho)
        mom = [cf.gaussian_moment(j, r) for j in range(N)]
        tp = 2 * mpmath.pi
        for n in range(N):
            out[n] = float(mpmath.fsum(mpmath.mpf(c.numerator) / c.denominator * tp ** j * mom[j]
                                       for j, c in enumerate(laguerre_y(n))))
    return out


# =========================================================================== the SDP (explicit affine data)

@dataclass
class LinearModel:
    """minimise sum_b <C_b, X_b> s.t. sum_b <A_kb, X_b> = b_k, X_b PSD (blocks X2, X3, X4)."""

    m: int
    rho: float
    basis: str
    C: dict
    A: list            # list of (name, {block: matrix}, rhs)
    qmap: np.ndarray   # a = qmap @ vec(X2) gives the w-coordinates of q (length 2m + 2)


def build_model(m: int, rho: float, basis: str = "hermite") -> LinearModel:
    if isinstance(m, bool) or not isinstance(m, (int, np.integer)) or m < 0:
        raise ValueError("m must be a non-negative integer")
    if isinstance(rho, bool) or not isinstance(rho, (int, float, np.integer, np.floating)) \
            or not math.isfinite(rho) or rho <= 0:
        raise ValueError("rho must be a finite number > 0")
    if basis not in BASES:
        raise ValueError(f"unknown basis {basis!r}")
    m, rho = int(m), float(rho)
    N = 2 * m + 2                      # coordinates of q and Tq: degrees 0..2m+1
    P = product_tensor(basis, m)       # (n, n, 2m + 1)
    Y = multiply_by_y_matrix(2 * m + 1)
    inv2pi = 1.0 / (2.0 * math.pi)
    # Coordinates as linear maps of the Gram matrices: S(X)_k = sum_ij X_ij P[i, j, k].
    # q = rho^2 S2 - inv2pi * Y S2 (length N);  rhs of matching = S3 + inv2pi * Y S4.
    E = np.zeros((N, 2 * m + 1))
    E[: 2 * m + 1, :] = np.eye(2 * m + 1)
    Mq = rho ** 2 * E - inv2pi * Y                     # N x (2m+1), acting on S2
    Mf = inv2pi * Y                                    # N x (2m+1), acting on S4
    sign = np.array([(-1.0) ** k for k in range(N)])
    w0 = np.array([float(w_at_zero(k)) for k in range(N)])
    J = objective_integrals(N, rho)

    def block(vec_over_S):            # matrix G with <G, X> = sum_k vec_k S(X)_k
        return np.tensordot(P, vec_over_S, axes=([2], [0]))

    A = []
    for k in range(N):                # (Tq)_k - S3_k - (inv2pi Y S4)_k = 0
        e3 = np.zeros(2 * m + 1)
        if k < 2 * m + 1:
            e3[k] = 1.0
        A.append((f"match_{k}", {"X2": block(sign[k] * Mq[k]), "X3": block(-e3), "X4": block(-Mf[k])}, 0.0))
    A.append(("tq0", {"X2": block((sign * w0) @ Mq)}, 1.0))
    C = {"X2": block(rho * (w0 @ Mq) + (2.0 / rho) * (J @ Mq))}
    return LinearModel(m, rho, basis, C, A, Mq)


@dataclass
class S1Result:
    m: int
    rho: float
    basis: str
    status: str
    solver: str
    objective: float | None = None
    q_w: list = field(default_factory=list)        # q in the w basis
    q0: float | None = None
    tq0: float | None = None
    R_physical: float | None = None
    gram_min_eig: dict = field(default_factory=dict)
    gram_max_eig: dict = field(default_factory=dict)
    dual_min_eig: dict = field(default_factory=dict)   # minimum eigenvalue of each PSD dual block
    residual_max_abs: float | None = None
    kkt_stationarity_max_abs: float | None = None
    duality_gap: float | None = None
    complementarity: float | None = None             # sum_b <X_b, S_b>
    seconds: float = 0.0
    notes: list = field(default_factory=list)
    grams: dict = field(default_factory=dict)        # primal Gram matrices X2, X3, X4
    duals: dict = field(default_factory=dict)        # {"eq": [y_k], "psd": {block: S_b}}

    def as_dict(self, with_witness: bool = False) -> dict:
        d = dict(self.__dict__)
        if not with_witness:
            d.pop("grams")
            d.pop("duals")
        return d


BLOCKS = ("X2", "X3", "X4")


def witness_diagnostics(model: LinearModel, X: dict, y=None, S: dict | None = None) -> dict:
    """Diagnostics of a primal (and optionally dual) point from the explicit affine data alone (numpy, no solver).

    Dual convention (CVXPY): L = <C, X> + sum_k y_k (A_k(X) - b_k) - sum_b <S_b, X_b>, so stationarity is
    C_b + sum_k y_k A_kb - S_b = 0, the dual objective is -sum_k y_k b_k, and dual feasibility is S_b PSD.
    """
    out = {}
    out["objective"] = float(sum(np.sum(M * X[b]) for b, M in model.C.items()))
    out["residual_max_abs"] = float(max(abs(sum(np.sum(M * X[b]) for b, M in mats.items()) - rhs)
                                        for _, mats, rhs in model.A))
    S2 = np.einsum("ij,ijk->k", X["X2"], product_tensor(model.basis, model.m))
    a = model.qmap @ S2
    w0 = np.array([float(w_at_zero(k)) for k in range(len(a))])
    out["q_w"] = a.tolist()
    out["q0"] = float(w0 @ a)
    out["tq0"] = float((np.array([(-1.0) ** k for k in range(len(a))]) * w0) @ a)
    out["R_physical"] = (model.rho * out["q0"] / out["tq0"]
                         if math.isfinite(out["tq0"]) and out["tq0"] != 0 else None)
    out["gram_min_eig"], out["gram_max_eig"] = {}, {}
    for b in BLOCKS:
        ev = np.linalg.eigvalsh(X[b])
        out["gram_min_eig"][b], out["gram_max_eig"][b] = float(ev[0]), float(ev[-1])
    if y is not None and S is not None:
        y = np.asarray(y, dtype=float)
        worst = 0.0
        for b in BLOCKS:
            Rm = model.C.get(b, 0) - S[b]
            for k, (_, mats, _) in enumerate(model.A):
                if b in mats:
                    Rm = Rm + y[k] * mats[b]
            worst = max(worst, float(np.max(np.abs(Rm))))
        out["kkt_stationarity_max_abs"] = worst
        out["duality_gap"] = out["objective"] - float(-sum(y[k] * rhs for k, (_, _, rhs) in enumerate(model.A)))
        out["dual_min_eig"] = {b: float(np.linalg.eigvalsh(S[b])[0]) for b in BLOCKS}
        out["complementarity"] = float(sum(np.sum(X[b] * S[b]) for b in BLOCKS))
    return out


# Clarabel settings tried in order. Interior-point iterations on these near-degenerate SDPs (optimal
# Gram matrices are rank deficient, Fourier side near contact) often stop with NumericalError at the
# default static regularisation; a larger static regularisation lets them converge. Every attempt is
# recorded, and acceptance is decided from recomputed residuals, not from the solver status alone.
SETTINGS_LADDER = (
    ("strict", dict(tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10)),
    ("reg1e-7-noequil", dict(tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10,
                             equilibrate_enable=False, static_regularization_constant=1e-7)),
    ("reg1e-6", dict(static_regularization_constant=1e-6)),
)


def solve(m: int, rho: float = 1.0, basis: str = "paper", ladder=SETTINGS_LADDER) -> S1Result:
    """Solve formulation B with Clarabel (float64). Returns a numerical candidate with diagnostics.

    The rungs of `ladder` are tried in order until one reports OPTIMAL. The reported point is the one
    returned by the last rung that produced a solution (OPTIMAL or OPTIMAL_INACCURATE), and `solver`
    names that rung; a rung that raises a solver error never overwrites it. If no rung produced a
    solution, the status is "solver-error" and no point is reported.
    """
    import clarabel
    import cvxpy as cp

    t0 = time.time()
    model = build_model(m, rho, basis)
    n = m + 1
    X = {b: cp.Variable((n, n), symmetric=True, name=b) for b in BLOCKS}

    def lin(mats):
        return sum(cp.sum(cp.multiply(M, X[b])) for b, M in mats.items())

    cons = [X[b] >> 0 for b in BLOCKS]
    eq = [lin(mats) == rhs for _, mats, rhs in model.A]
    prob = cp.Problem(cp.Minimize(lin(model.C)), cons + eq)
    res = S1Result(m, rho, basis, "not-run", f"CLARABEL {clarabel.__version__}")
    snap = None        # (rung, status, primal, eq duals, psd duals) of the last rung that produced a solution
    for name, opts in ladder:
        try:
            prob.solve(solver=cp.CLARABEL, **opts)
        except cp.error.SolverError:
            res.notes.append(f"settings {name}: solver error (NumericalError)")
            continue
        status = prob.status
        res.notes.append(f"settings {name}: {status}")
        if status in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE) and all(X[b].value is not None for b in BLOCKS):
            ys = [c.dual_value for c in eq]
            Ss = [c.dual_value for c in cons]
            duals_ok = all(v is not None for v in ys + Ss)
            snap = (name, status, {b: np.array(X[b].value, dtype=float) for b in BLOCKS},
                    [float(np.asarray(v)) for v in ys] if duals_ok else None,
                    {b: np.array(Ss[i], dtype=float) for i, b in enumerate(BLOCKS)} if duals_ok else None)
        if status == cp.OPTIMAL:
            break
    res.seconds = time.time() - t0
    if snap is None:
        res.status = "solver-error"
        res.notes.append("no settings rung produced a solution")
        return res
    name, status, Xraw, y, Sraw = snap
    res.solver = f"CLARABEL {clarabel.__version__} (settings {name})"
    res.status = str(status)
    Xv = {b: (Xraw[b] + Xraw[b].T) / 2 for b in BLOCKS}
    Sv = None if Sraw is None else {b: (Sraw[b] + Sraw[b].T) / 2 for b in BLOCKS}
    res.grams = {b: Xv[b].tolist() for b in BLOCKS}
    if Sv is not None:
        res.duals = {"eq": list(y), "psd": {b: Sv[b].tolist() for b in BLOCKS}}
    else:
        res.notes.append("solver returned no dual values: KKT diagnostics not computed")
    dg = witness_diagnostics(model, Xv, y, Sv)
    for key, val in dg.items():
        setattr(res, key, val)
    if status != cp.OPTIMAL:
        res.notes.append("solver reported OPTIMAL_INACCURATE")
    if not (math.isfinite(res.q0) and res.q0 > Q0_MIN):
        res.notes.append(f"q(0) is not > {Q0_MIN}: not normalisable")
    if res.R_physical is not None and res.R_physical < 1 - 1e-6:
        res.notes.append("R < 1: impossible for an exactly admissible function (numerical artefact)")
    return res


Q0_MIN = 1e-8


def acceptance_thresholds(tol: float = 1e-7) -> dict:
    """The numerical acceptance rule of `acceptable`, written out (recorded with every sweep)."""
    return {"solver_status": "optimal", "residual_max_abs_max": tol, "abs_duality_gap_max": tol,
            "kkt_stationarity_max_abs_max": tol, "gram_min_eig_min": -tol, "dual_min_eig_min": -tol,
            "q0_min_exclusive": Q0_MIN, "all_diagnostics_present_and_finite": True}


def _finite(v) -> bool:
    return isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool) and math.isfinite(v)


def acceptable(r: S1Result, tol: float = 1e-7) -> bool:
    """Numerical acceptance used by the sweep (not a certificate); see `acceptance_thresholds`.

    Never raises on missing or non-finite diagnostics: such a candidate is simply not accepted.
    """
    if r.status != "optimal" or not _finite(tol) or tol <= 0:
        return False
    scalars = (r.objective, r.residual_max_abs, r.duality_gap, r.kkt_stationarity_max_abs, r.q0, r.tq0)
    if not all(_finite(v) for v in scalars):
        return False
    for eigs in (r.gram_min_eig, r.dual_min_eig):
        if not isinstance(eigs, dict) or set(eigs) != set(BLOCKS) or not all(_finite(v) for v in eigs.values()):
            return False
    return (r.residual_max_abs <= tol and abs(r.duality_gap) <= tol and r.kkt_stationarity_max_abs <= tol
            and min(r.gram_min_eig.values()) >= -tol and min(r.dual_min_eig.values()) >= -tol
            and r.q0 > Q0_MIN)


# =========================================================================== saved witnesses and solver-free replay

WITNESS_RECORDED_KEYS = ("objective", "residual_max_abs", "kkt_stationarity_max_abs", "duality_gap", "q0", "tq0",
                         "R_physical", "complementarity", "gram_min_eig", "dual_min_eig")
REPLAY_AGREEMENT_ABS = 1e-12      # recomputed versus recorded diagnostics (same affine data, same floats)


def witness_record(r: S1Result, tol: float) -> dict:
    """Primal and dual witness of a solve, with the diagnostics recorded at solve time (JSON-ready)."""
    if not r.grams or not r.duals:
        raise ValueError("result has no primal/dual witness")
    return {"m": r.m, "rho": r.rho, "basis": r.basis, "solver": r.solver, "status": r.status, "tol": tol,
            "X": r.grams, "eq_duals": r.duals["eq"], "psd_duals": r.duals["psd"],
            "recorded": {**{k: getattr(r, k) for k in WITNESS_RECORDED_KEYS}, "accepted": acceptable(r, tol)}}


def replay_witness(w: dict) -> dict:
    """Re-check a saved witness without any solver: rebuild the affine data from (m, rho, basis), recompute
    every diagnostic with numpy, re-apply the acceptance thresholds and compare with the recorded values.

    Returns {"ok": bool, "checks": {name: bool}, "errors": [...], "recomputed": {...}}. Malformed or
    non-finite input gives ok = False (never an exception). Numerical check only, not a certificate.
    """
    checks, errors = {}, []

    def fail(msg):
        errors.append(msg)
        return {"ok": False, "checks": checks, "errors": errors, "recomputed": None}

    try:
        m, rho, basis, tol = w["m"], w["rho"], w["basis"], w["tol"]
        model = build_model(m, rho, basis)
        if not _finite(tol) or tol <= 0:
            return fail("tol must be finite and > 0")
        X, S, y = {}, {}, np.asarray(w["eq_duals"], dtype=float)
        for b in BLOCKS:
            X[b] = np.asarray(w["X"][b], dtype=float)
            S[b] = np.asarray(w["psd_duals"][b], dtype=float)
        rec = w["recorded"]
    except (KeyError, TypeError, ValueError) as exc:
        return fail(f"malformed witness: {type(exc).__name__}: {exc}")
    n = m + 1
    checks["shapes"] = (y.shape == (len(model.A),)
                        and all(X[b].shape == (n, n) and S[b].shape == (n, n) for b in BLOCKS))
    if not checks["shapes"]:
        return fail("wrong array shapes")
    checks["finite"] = bool(np.all(np.isfinite(y)) and all(np.all(np.isfinite(X[b])) and np.all(np.isfinite(S[b]))
                                                            for b in BLOCKS))
    if not checks["finite"]:
        return fail("non-finite entries")
    checks["symmetric"] = all(np.array_equal(X[b], X[b].T) and np.array_equal(S[b], S[b].T) for b in BLOCKS)
    dg = witness_diagnostics(model, X, y, S)
    dg.pop("q_w")
    checks["residual"] = dg["residual_max_abs"] <= tol
    checks["duality_gap"] = abs(dg["duality_gap"]) <= tol
    checks["kkt_stationarity"] = dg["kkt_stationarity_max_abs"] <= tol
    checks["primal_psd"] = min(dg["gram_min_eig"].values()) >= -tol
    checks["dual_psd"] = min(dg["dual_min_eig"].values()) >= -tol
    checks["q0_positive"] = dg["q0"] > Q0_MIN
    agree = True
    for key in WITNESS_RECORDED_KEYS:
        a, b = dg[key], rec.get(key) if isinstance(rec, dict) else None
        if isinstance(a, dict):
            ok = isinstance(b, dict) and set(b) == set(a) and all(
                _finite(b[k]) and abs(a[k] - b[k]) <= REPLAY_AGREEMENT_ABS for k in a)
        else:
            ok = _finite(a) and _finite(b) and abs(a - b) <= REPLAY_AGREEMENT_ABS
        if not ok:
            agree = False
            errors.append(f"recomputed {key} disagrees with the recorded value")
    checks["agrees_with_recorded"] = agree
    accepted_now = all(checks[k] for k in ("residual", "duality_gap", "kkt_stationarity", "primal_psd", "dual_psd",
                                           "q0_positive"))
    checks["recorded_acceptance_flag"] = isinstance(rec, dict) and rec.get("accepted") is accepted_now
    for k, v in checks.items():
        if not v and not any(k in e for e in errors):
            errors.append(f"check failed: {k}")
    return {"ok": all(checks.values()), "checks": checks, "errors": errors, "recomputed": dg}


# =========================================================================== independent re-evaluation (B1a)

def q_as_even_poly(q_w: list[float]) -> cf.EvenPoly:
    """Exact q(x^2) in powers of x^2 from the binary floats of the w-coordinates: coefficients in Q[pi]."""
    coeffs: dict[int, Fraction] = {}
    for n, a in enumerate(q_w):
        fa = Fraction(a)
        for j, c in enumerate(laguerre_y(n)):
            coeffs[j] = coeffs.get(j, Fraction(0)) + fa * c * 2 ** j
    return cf.even_poly(cf.PiLaurent.pi_power(j, coeffs.get(j, Fraction(0))) for j in range(len(q_w)))


def independent_check(r: S1Result, dps: int = 60, fourier_minimum: bool = True) -> dict:
    """Re-evaluate a candidate with B1a: exact T (Laguerre formula) and two objective routes.

    Uses the same w-coordinates but none of the S1 linear maps. numerical-high-precision only.
    """
    p = q_as_even_poly(r.q_w)
    tp = cf.transform_laguerre(p)
    out = {}
    with mpmath.workdps(dps + 4 * len(r.q_w)):
        rho = mpmath.mpf(r.rho)
        q0 = cf.poly_value_at_zero(p).to_mp()
        tq0 = cf.poly_value_at_zero(tp).to_mp()
        integral = mpmath.fsum(c.to_mp() * cf.gaussian_moment(k, rho) for k, c in enumerate(p))
        z_moments = (q0 * rho + 2 * integral / rho) / tq0
        fx = lambda x: x * cf.poly_eval_mp(p, x) * mpmath.exp(-mpmath.pi * x * x)  # noqa: E731
        Iq = mpmath.quad(fx, mpmath.linspace(0, rho, 17))
        z_quad = (q0 * rho + 2 * Iq / rho) / tq0
        out["z_moments"] = float(z_moments)
        out["z_quadrature"] = float(z_quad)
        out["z_routes_abs_diff"] = float(abs(z_moments - z_quad))
        out["z_vs_solver_abs_diff"] = float(abs(z_moments - mpmath.mpf(r.objective)))
        out["tq0_b1a"] = float(tq0)
        out["q_at_rho_sq"] = float(cf.poly_eval_mp(p, rho))       # should be ~0: built-in sign change
    if fourier_minimum:
        # Unweighted polynomial (Tq)(xi^2) over xi >= 0 (not fhat = (Tq)(xi^2) exp(-pi xi^2)).
        mn = cf.even_poly_minimum(tp, dps)
        out["tq_unweighted_min_value"] = None if mn["min_value"] is None else float(mn["min_value"])
        out["tq_unweighted_argmin_xi"] = None if mn["argmin_u"] is None else float(mn["argmin_u"])
        out["tq_unbounded_below"] = mn["unbounded_below"]
    return out


# =========================================================================== Fourier side: two different minima

def _fujiwara_bound(a: list) -> mpmath.mpf:
    """Fujiwara's bound on the moduli of all complex roots of sum_k a_k t^k (a ascending, a[-1] != 0)."""
    n = len(a) - 1
    if n <= 0:
        return mpmath.mpf(0)
    an = abs(a[-1])
    terms = [(abs(a[n - k]) / an) ** (mpmath.mpf(1) / k) for k in range(1, n)]
    terms.append((abs(a[0]) / (2 * an)) ** (mpmath.mpf(1) / n))
    return 2 * max(terms)


def _positive_real_roots(a: list, dps: int) -> tuple[list, bool, float]:
    """Positive real roots of sum_k a_k t^k (ascending) via mpmath.polyroots; (roots, converged, max |residual|).

    The residual of each accepted root is relative to sum_k |a_k| t^k.
    """
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    if len(a) < 2:
        return [], True, 0.0
    try:
        roots = mpmath.polyroots(a, maxsteps=2000, extraprec=4 * dps, asc=True)
    except mpmath.libmp.NoConvergence:
        return [], False, float("inf")
    thr = mpmath.mpf(10) ** (-(dps // 4))
    out, worst = [], mpmath.mpf(0)
    for z in roots:
        z = mpmath.mpc(z)
        if abs(z.imag) <= thr * (1 + abs(z)) and z.real > 0:
            t = z.real
            scale = mpmath.fsum(abs(c) * t ** k for k, c in enumerate(a))
            worst = max(worst, abs(mpmath.polyval(a, t, asc=True)) / scale)
            out.append(t)
    return sorted(out), True, float(worst)


def fourier_side_minima(P: list, dps: int = 200, grid_points: int = 20001) -> dict:
    """Two distinct minima of the Fourier side of a candidate, over the half-line t = xi^2 >= 0.

    P is (Tq)(t) as ascending coefficients in t, so fhat(xi) = P(xi^2) exp(-pi xi^2). Reported separately:

    * unweighted: min of P(t) on t >= 0. Candidates t = 0 and the positive roots of P'(t) = 0; tail: P -> +inf
      if the leading coefficient is positive, otherwise P is unbounded below (no minimum).
    * weighted: min of g(t) = P(t) exp(-pi t) on t >= 0. Candidates t = 0 and the positive roots of
      P'(t) - pi P(t) = 0 (g' = (P' - pi P) exp(-pi t)); tail: g -> 0, from the side of the sign of the leading
      coefficient. The infimum is min(candidate values, 0); if every candidate value is positive the
      infimum 0 is approached in the tail and not attained.

    Reliability (numerical, not rigorous): polyroots must converge, accepted real roots must have small
    relative residuals, and a dense grid on [0, B] (B = Fujiwara bound of both stationary polynomials, so it
    contains every stationary point) must not find a value below the stationary-point minimum. If any of
    these fails, `reliable` is False and no global minimum is claimed. This is a numerical diagnostic of the
    float candidate, never a certificate of Fourier non-negativity.
    """
    with mpmath.workdps(dps):
        c = [mpmath.mpf(x) for x in P]
        while len(c) > 1 and c[-1] == 0:
            c.pop()
        n = len(c) - 1
        pi = mpmath.pi
        lead = c[-1]
        dP = [k * c[k] for k in range(1, n + 1)]
        G = [(k + 1) * c[k + 1] - pi * c[k] for k in range(n)] + [-pi * c[n]]   # P' - pi P
        r_dP, conv_dP, res_dP = _positive_real_roots(dP, dps)
        r_G, conv_G, res_G = _positive_real_roots(G, dps)

        def Pv(t):
            return mpmath.polyval(c, t, asc=True)

        def gv(t):
            return Pv(t) * mpmath.exp(-pi * t)

        limitations = []
        # unweighted
        if n >= 1 and lead < 0:
            unw = {"min_value": None, "argmin_t": None, "argmin_xi": None, "unbounded_below": True,
                   "tail": "-inf (negative leading coefficient)"}
        else:
            cand = [mpmath.mpf(0)] + r_dP
            v, t = min((Pv(t), t) for t in cand)
            unw = {"min_value": float(v), "argmin_t": float(t), "argmin_xi": float(mpmath.sqrt(t)),
                   "unbounded_below": False, "tail": "+inf" if n >= 1 else "constant"}
        unw.update({"value_at_t0": float(c[0]), "leading_coefficient": float(lead),
                    "n_positive_stationary_points": len(r_dP)})
        # weighted
        cand = [mpmath.mpf(0)] + r_G
        v, t = min((gv(t), t) for t in cand)
        tail_sign = "+" if lead > 0 else ("-" if lead < 0 else "0")
        wtd = {"min_value_attained": float(v), "argmin_t": float(t), "argmin_xi": float(mpmath.sqrt(t)),
               "infimum": float(min(v, 0)), "infimum_attained": bool(v <= 0),
               "tail": f"-> 0 from the {tail_sign} side (leading coefficient sign)",
               "value_at_t0": float(c[0]), "n_positive_stationary_points": len(r_G)}
        # cross-check on a dense grid covering every stationary point
        B = max(_fujiwara_bound(dP) if len(dP) > 1 else mpmath.mpf(0),
                _fujiwara_bound(G) if len(G) > 1 else mpmath.mpf(0), mpmath.mpf(1))
        xi_max = mpmath.sqrt(B * mpmath.mpf("1.05"))
        g_min_grid, p_min_grid = mpmath.inf, mpmath.inf
        for k in range(grid_points):
            t = (xi_max * k / (grid_points - 1)) ** 2
            pv = Pv(t)
            p_min_grid = min(p_min_grid, pv)
            g_min_grid = min(g_min_grid, pv * mpmath.exp(-pi * t))
        slack_w = mpmath.mpf(10) ** (-(dps // 2)) + mpmath.mpf("1e-9") * abs(wtd["min_value_attained"])
        grid_ok_w = g_min_grid >= wtd["min_value_attained"] - slack_w
        if unw["unbounded_below"]:
            grid_ok_u = True
        else:
            slack_u = mpmath.mpf(10) ** (-(dps // 2)) + mpmath.mpf("1e-9") * abs(unw["min_value"])
            grid_ok_u = p_min_grid >= unw["min_value"] - slack_u
        res_ok = res_dP <= 1e-20 and res_G <= 1e-20
        if not (conv_dP and conv_G):
            limitations.append("polyroots did not converge: stationary points incomplete")
        if not res_ok:
            limitations.append("a real root has a large relative residual")
        if not (grid_ok_u and grid_ok_w):
            limitations.append("grid found a value below the stationary-point minimum: a stationary point was missed")
        reliable = conv_dP and conv_G and res_ok and grid_ok_u and grid_ok_w
        if not reliable:
            limitations.append("global minimum NOT established")
        return {"variable": "t = xi^2 >= 0; fhat(xi) = P(xi^2) exp(-pi xi^2), P = Tq",
                "unweighted_P": unw, "weighted_P_exp_minus_pi_t": wtd,
                "stationary_point_bound_t": float(B),
                "grid_cross_check": {"xi_max": float(xi_max), "points": grid_points,
                                     "min_P_on_grid": float(p_min_grid), "min_weighted_on_grid": float(g_min_grid),
                                     "consistent_unweighted": bool(grid_ok_u), "consistent_weighted": bool(grid_ok_w)},
                "root_residual_max_rel": {"P_prime": res_dP, "P_prime_minus_pi_P": res_G},
                "polyroots_converged": bool(conv_dP and conv_G), "reliable": bool(reliable),
                "limitations": limitations, "dps": dps,
                "rigor": "numerical-high-precision diagnostic of the float candidate; not a certificate"}


def golden_minimise(fun, lo: float, hi: float, iters: int = 14):
    """Golden-section search on [lo, hi] for a function that may be non-smooth (bracketing only)."""
    g = (math.sqrt(5) - 1) / 2
    c, d = hi - g * (hi - lo), lo + g * (hi - lo)
    fc, fd = fun(c), fun(d)
    evals = [(c, fc), (d, fd)]
    for _ in range(iters):
        if fc <= fd:
            hi, d, fd = d, c, fc
            c = hi - g * (hi - lo)
            fc = fun(c)
            evals.append((c, fc))
        else:
            lo, c, fc = c, d, fd
            d = lo + g * (hi - lo)
            fd = fun(d)
            evals.append((d, fd))
    return min(evals, key=lambda t: t[1]), evals


def highprec_check(r: S1Result, dps: int = 200, grid: int = 4001, xi_max: float = 4.0,
                   fourier_grid_points: int = 20001) -> dict:
    """Re-evaluate a candidate from its Gram matrices alone, independently of the w-basis linear maps.

    The binary floats of X2, X3, X4 are taken as exact data. In mpmath at `dps` digits:
    s_b(t) = v(t)^T X_b v(t) with v_k(t) = L_k^(-1/2)(sigma t) expanded in monomials of t,
    q = (rho^2 - t) s_2, Tq by B1a's monomial formula T(t^k) = (k!/pi^k) L_k^(-1/2)(pi xi^2), and
    E = s_3 + t s_4 - Tq (the coefficient-matching defect). Reports max |E|, (Tq)(0), q(0), R, Z by
    Gaussian moments, and `fourier_side_minima` of Tq: the unweighted minimum of Tq(t) and the weighted
    minimum of Tq(t) exp(-pi t) over the half-line t = xi^2 >= 0, kept separate (the Fourier side of this
    float candidate; a numerical diagnostic, not a certificate). The matching defect is sampled on a grid
    xi in [0, xi_max].
    """
    m = r.m
    with mpmath.workdps(dps):
        pi = mpmath.pi
        lam = mpmath.mpf(BASES[r.basis].numerator) / BASES[r.basis].denominator
        sig = 2 * pi * lam

        def lag_t(n, scale):           # L_n^(-1/2)(scale * t) in monomials of t
            return [mpmath.mpf(c.numerator) / c.denominator * scale ** j for j, c in enumerate(laguerre_y(n))]

        V = [lag_t(k, sig) for k in range(m + 1)]

        def pmul(a, b):
            out = [mpmath.mpf(0)] * (len(a) + len(b) - 1)
            for i, x in enumerate(a):
                for j, y in enumerate(b):
                    out[i + j] += x * y
            return out

        def padd(a, b):
            n = max(len(a), len(b))
            return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]

        def quad(Xb):
            Xm = [[mpmath.mpf(Fraction(v).numerator) / Fraction(v).denominator for v in row] for row in Xb]
            out = [mpmath.mpf(0)]
            for i in range(m + 1):
                wi = [mpmath.mpf(0)]
                for j in range(m + 1):
                    wi = padd(wi, [Xm[i][j] * c for c in V[j]])
                out = padd(out, pmul(V[i], wi))
            return out

        s2, s3, s4 = (quad(r.grams[b]) for b in ("X2", "X3", "X4"))
        rho = mpmath.mpf(r.rho)
        q = pmul([rho * rho, mpmath.mpf(-1)], s2)
        tq = [mpmath.mpf(0)]
        for k, c in enumerate(q):
            tq = padd(tq, [c * mpmath.factorial(k) / pi ** k * x for x in lag_t(k, pi)])
        rhs = padd(s3, [mpmath.mpf(0)] + s4)
        E = padd(rhs, [-x for x in tq])
        integral = mpmath.fsum(c * cf.gaussian_moment(k, rho) for k, c in enumerate(q))
        z = (q[0] * rho + 2 * integral / rho) / tq[0]
        Ef = lambda u: mpmath.polyval(E, u, asc=True)  # noqa: E731
        xs = [mpmath.mpf(xi_max) * k / (grid - 1) for k in range(grid)]
        defect_fn = max(abs(Ef(x * x)) * mpmath.exp(-pi * x * x) for x in xs)
        fourier = fourier_side_minima(tq, dps, fourier_grid_points)
        return {"dps": dps, "matching_defect_fn_max": float(defect_fn), "matching_defect_xi_max": xi_max,
                "tq0": float(tq[0]), "q0": float(q[0]), "R_physical": float(rho * q[0] / tq[0]),
                "z_highprec": float(z), "z_vs_solver_abs_diff": float(abs(z - mpmath.mpf(r.objective))),
                "fourier": fourier,
                "note": ("Tq of the float candidate itself, so negative values measure its Fourier-side "
                         "infeasibility; matching_defect_fn_max = max over xi in [0, matching_defect_xi_max] of "
                         "|E(xi^2)| exp(-pi xi^2)")}
