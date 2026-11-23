"""Independent re-implementation for the audit (written from the math, not from the audited code).

Test function family:  h(t) = sum_j q_j * (1/L_j) 1_{[-L_j/2, L_j/2]}(t),  sum_j q_j = 1.
  H(x) = FT h(x) = sum_j q_j sinc(L_j x)          (sinc(y) = sin(pi y)/(pi y))
  f(x) = |H(x)|^2 >= 0,  f(0) = (sum q)^2 = 1
  R(a) = FT f (a) = (h * h~)(a) = q^T K(a) q,     K_jk(a) = |I_j ∩ (I_k + a)| / (L_j L_k)
Montgomery (RH):  N*(T) <= (1+o(1)) N(T) * [ R(0) + int_{-1}^{1} |a| R(a) da + 2 int_{1}^{inf} F(a) R(a) da ]
"""
from fractions import Fraction as Fr
import numpy as np


# ---------- exact rational arithmetic ----------
def overlap_exact(L, M, x):
    x = abs(x)
    s = (L + M) / 2
    d = abs(L - M) / 2
    if x >= s:
        return Fr(0)
    if x <= d:
        return min(L, M)
    return s - x


def K_exact(L, M, x):
    return overlap_exact(L, M, x) / (L * M)


def xK_integral_exact(L, M, upper=Fr(1)):
    """int_0^upper x K_LM(x) dx, exact, via antiderivatives on each linear piece."""
    s = (L + M) / 2
    d = abs(L - M) / 2
    m = min(L, M)
    tot = Fr(0)
    a = min(d, upper)                       # flat piece [0, a]: value m
    tot += m * a * a / 2
    lo, hi = min(d, upper), min(s, upper)    # linear piece: s - x
    if hi > lo:
        F = lambda x: s * x * x / 2 - x ** 3 / 3
        tot += F(hi) - F(lo)
    return tot / (L * M)


def R_exact(q, W, x):
    n = len(W)
    return sum(q[i] * q[j] * K_exact(W[i], W[j], x) for i in range(n) for j in range(n))


def C_window_exact(q, W):
    n = len(W)
    r0 = sum(q[i] * q[j] * K_exact(W[i], W[j], Fr(0)) for i in range(n) for j in range(n))
    I = sum(q[i] * q[j] * xK_integral_exact(W[i], W[j]) for i in range(n) for j in range(n))
    return r0 + 2 * I


def knots_exact(W, lo, hi):
    pts = {lo, hi}
    for L in W:
        for M in W:
            for p in (abs(L - M) / 2, (L + M) / 2):
                if lo < p < hi:
                    pts.add(p)
    return sorted(pts)


# ---------- float, vectorised ----------
def Kmat(W, x):
    W = np.asarray(W, float)
    L, M = W[:, None], W[None, :]
    s, d = (L + M) / 2, np.abs(L - M) / 2
    x = abs(x)
    ov = np.where(x <= d, np.minimum(L, M), np.where(x >= s, 0.0, s - x))
    return ov / (L * M)


def xKmat(W):
    W = [Fr(w) for w in W]  # exact on the binary floats
    return np.array([[float(xK_integral_exact(a, b)) for b in W] for a in W])


def window_matrix(W):
    return Kmat(W, 0.0) + 2 * xKmat(W)


def knots_float(W, lo, hi):
    W = np.asarray(W, float)
    p = np.concatenate([np.abs(W[:, None] - W[None, :]).ravel() / 2, (W[:, None] + W[None, :]).ravel() / 2])
    p = p[(p > lo) & (p < hi)]
    return np.unique(np.concatenate([[lo, hi], p]))


# ---------- x-space (physical side) cross-check, completely different route ----------
def C_window_xspace(q, W, X=3000.0, dx=2e-4):
    """C_window = int f(x) * mu_hat(x) dx with mu = delta_0 + |a| 1_{[-1,1]}.
    mu_hat(x) = 1 + 2[ sin w / w + (cos w - 1)/w^2 ],  w = 2 pi x.
    Returns (R0_numeric, C_numeric)."""
    x = np.arange(dx / 2, X, dx)  # even integrand: integrate over (0,X) and double
    H = sum(qj * np.sinc(Lj * x) for qj, Lj in zip(q, W))
    f = H * H
    w = 2 * np.pi * x
    g = 2 * (np.sin(w) / w + (np.cos(w) - 1) / w ** 2)
    tail = sum(qj ** 2 / (np.pi ** 2 * Lj ** 2 * X) for qj, Lj in zip(q, W))  # mean of f tail, both sides
    R0 = 2 * np.sum(f) * dx + tail
    C = R0 + 2 * np.sum(f * g) * dx
    return R0, C


def R_direct_convolution(q, W, alphas, dt=1e-4):
    """R(a) = int h(t) h(t+a) dt by brute-force quadrature on a grid (independent of K formula)."""
    A = max(W)
    t = np.arange(-A / 2 - 1, A / 2 + 1, dt)
    h = lambda tt: sum(qj * (np.abs(tt) <= Lj / 2) / Lj for qj, Lj in zip(q, W))
    ht = h(t)
    return np.array([np.sum(ht * h(t + a)) * dt for a in alphas])
