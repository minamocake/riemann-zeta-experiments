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
    if m < 0 or rho <= 0 or basis not in BASES:
        raise ValueError("need m >= 0, rho > 0 and a known basis")
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
    residual_max_abs: float | None = None
    kkt_stationarity_max_abs: float | None = None
    duality_gap: float | None = None
    seconds: float = 0.0
    notes: list = field(default_factory=list)
    grams: dict = field(default_factory=dict)

    def as_dict(self, with_grams: bool = False) -> dict:
        d = dict(self.__dict__)
        if not with_grams:
            d.pop("grams")
        return d


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
    """Solve formulation B with Clarabel (float64). Returns a numerical candidate with diagnostics."""
    import clarabel
    import cvxpy as cp

    t0 = time.time()
    model = build_model(m, rho, basis)
    n = m + 1
    X = {b: cp.Variable((n, n), symmetric=True, name=b) for b in ("X2", "X3", "X4")}

    def lin(mats):
        return sum(cp.sum(cp.multiply(M, X[b])) for b, M in mats.items())

    cons = [X[b] >> 0 for b in X]
    eq = [lin(mats) == rhs for _, mats, rhs in model.A]
    prob = cp.Problem(cp.Minimize(lin(model.C)), cons + eq)
    res = S1Result(m, rho, basis, "not-run", f"CLARABEL {clarabel.__version__}")
    for name, opts in ladder:
        try:
            prob.solve(solver=cp.CLARABEL, **opts)
        except cp.error.SolverError:
            res.notes.append(f"settings {name}: solver error (NumericalError)")
            continue
        res.notes.append(f"settings {name}: {prob.status}")
        if prob.status == cp.OPTIMAL:
            break
    else:
        if prob.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
            res.status = "solver-error"
            res.seconds = time.time() - t0
            return res
    res.solver = f"CLARABEL {clarabel.__version__} ({res.notes[-1].split(':')[0]})"
    res.status = str(prob.status)
    res.seconds = time.time() - t0
    if prob.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        res.notes.append("no solution reported")
        return res
    Xv = {b: (X[b].value + X[b].value.T) / 2 for b in X}
    res.grams = {b: Xv[b].tolist() for b in Xv}
    res.objective = float(sum(np.sum(M * Xv[b]) for b, M in model.C.items()))
    viol = [abs(sum(np.sum(M * Xv[b]) for b, M in mats.items()) - rhs) for _, mats, rhs in model.A]
    res.residual_max_abs = float(max(viol))
    S2 = np.einsum("ij,ijk->k", Xv["X2"], product_tensor(basis, m))
    a = model.qmap @ S2
    res.q_w = a.tolist()
    w0 = np.array([float(w_at_zero(k)) for k in range(len(a))])
    res.q0 = float(w0 @ a)
    res.tq0 = float((np.array([(-1.0) ** k for k in range(len(a))]) * w0) @ a)
    res.R_physical = rho * res.q0 / res.tq0
    for b in Xv:
        ev = np.linalg.eigvalsh(Xv[b])
        res.gram_min_eig[b], res.gram_max_eig[b] = float(ev[0]), float(ev[-1])
    # KKT from the explicit affine data (CVXPY convention: L = <C,X> + sum y_k (A_k(X) - b_k) - sum <S,X>).
    y = np.array([float(np.asarray(c.dual_value)) for c in eq])
    worst = 0.0
    for bi, b in enumerate(("X2", "X3", "X4")):
        Rm = model.C.get(b, 0) - np.asarray(cons[bi].dual_value)
        for k, (_, mats, _) in enumerate(model.A):
            if b in mats:
                Rm = Rm + y[k] * mats[b]
        worst = max(worst, float(np.max(np.abs(Rm))))
    res.kkt_stationarity_max_abs = worst
    res.duality_gap = res.objective - float(-sum(y[k] * rhs for k, (_, _, rhs) in enumerate(model.A)))
    if prob.status != cp.OPTIMAL:
        res.notes.append("solver reported OPTIMAL_INACCURATE")
    if res.q0 <= 1e-8:
        res.notes.append("q(0) is not > 1e-8: not normalisable")
    if res.R_physical < 1 - 1e-6:
        res.notes.append("R < 1: impossible for an exactly admissible function (numerical artefact)")
    return res


def acceptable(r: S1Result, tol: float = 1e-7) -> bool:
    """Numerical acceptance used by the sweep (not a certificate)."""
    return (r.status == "optimal" and r.residual_max_abs is not None and r.residual_max_abs <= tol
            and abs(r.duality_gap) <= tol and r.kkt_stationarity_max_abs <= 1e3 * tol
            and min(r.gram_min_eig.values()) >= -tol and r.q0 > 1e-8)


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
        mn = cf.even_poly_minimum(tp, dps)
        out["fourier_min_value"] = None if mn["min_value"] is None else float(mn["min_value"])
        out["fourier_min_u"] = None if mn["argmin_u"] is None else float(mn["argmin_u"])
        out["fourier_unbounded_below"] = mn["unbounded_below"]
    return out


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


def highprec_check(r: S1Result, dps: int = 200, grid: int = 4001, xi_max: float = 4.0) -> dict:
    """Re-evaluate a candidate from its Gram matrices alone, independently of the w-basis linear maps.

    The binary floats of X2, X3, X4 are taken as exact data. In mpmath at `dps` digits:
    s_b(t) = v(t)^T X_b v(t) with v_k(t) = L_k^(-1/2)(sigma t) expanded in monomials of t,
    q = (rho^2 - t) s_2, Tq by B1a's monomial formula T(t^k) = (k!/pi^k) L_k^(-1/2)(pi xi^2), and
    E = s_3 + t s_4 - Tq (the coefficient-matching defect). Reports max |E|, (Tq)(0), q(0), R, Z by
    Gaussian moments, and the minimum of Tq over a grid in xi refined by golden-section search
    (the Fourier side of this float candidate; exploratory, not a certificate).
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
        P = lambda u: mpmath.polyval(tq, u, asc=True)  # noqa: E731
        Ef = lambda u: mpmath.polyval(E, u, asc=True)  # noqa: E731
        xs = [mpmath.mpf(xi_max) * k / (grid - 1) for k in range(grid)]
        damp = [mpmath.exp(-pi * x * x) for x in xs]
        fh = [P(x * x) * dmp for x, dmp in zip(xs, damp, strict=True)]          # fhat(xi) of the candidate
        defect_fn = max(abs(Ef(x * x)) * dmp for x, dmp in zip(xs, damp, strict=True))
        vals = fh
        mins = []
        g = (mpmath.sqrt(5) - 1) / 2
        for i in range(1, grid - 1):
            if vals[i] <= vals[i - 1] and vals[i] <= vals[i + 1]:
                a, b = xs[i - 1], xs[i + 1]
                F = lambda x: P(x * x) * mpmath.exp(-pi * x * x)  # noqa: E731
                for _ in range(120):
                    c1, c2 = b - g * (b - a), a + g * (b - a)
                    if F(c1) < F(c2):
                        b = c2
                    else:
                        a = c1
                x0 = (a + b) / 2
                mins.append((float(x0), float(P(x0 * x0) * mpmath.exp(-pi * x0 * x0)), float(P(x0 * x0))))
        mins.sort(key=lambda t: t[1])
        return {"dps": dps, "matching_defect_fn_max": float(defect_fn),
                "tq0": float(tq[0]), "q0": float(q[0]), "R_physical": float(rho * q[0] / tq[0]),
                "z_highprec": float(z), "z_vs_solver_abs_diff": float(abs(z - mpmath.mpf(r.objective))),
                "fhat_local_minima_lowest3": mins[:3],
                "fhat_min_on_grid": float(min(fh)),
                "note": ("xi in units of the undilated function; entries (xi, fhat(xi), Tq(xi^2)); fhat of the "
                         "float candidate itself, so negative values measure its Fourier-side infeasibility; "
                         "matching_defect_fn_max = max over the grid of |E(xi^2)| exp(-pi xi^2)")}
