"""Centered-box test-function family: correlation kernels in float and exact rational arithmetic.

Family
------
    h(t) = sum_j q_j * (1/L_j) * 1_{[-L_j/2, L_j/2]}(t),        sum_j q_j = 1
    H(x) = (FT h)(x) = sum_j q_j * sinc(L_j x),                   sinc(y) = sin(pi y)/(pi y)
    f(x) = |H(x)|^2 >= 0,   f(0) = 1
    R(a) = (FT f)(a) = (h * h~)(a) = q^T K(a) q

More generally a positive semidefinite Gram matrix X (with 1^T X 1 = 1) defines the admissible
test function f(x) = v(x)^T X v(x), v_j(x) = sinc(L_j x), whose Fourier transform is
R(a) = tr(K(a) X).  The rank-one case X = q q^T is the family used by the legacy code.

Kernel
------
    K_LM(a) = |[-L/2, L/2] ∩ [a - M/2, a + M/2]| / (L M)

is piecewise linear in a with breakpoints |L-M|/2 and (L+M)/2, and vanishes for |a| >= (L+M)/2.
Hence every R is piecewise linear with those breakpoints ("knots"), and the extrema of R on a
closed interval are attained at knots or endpoints.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Sequence

import numpy as np


# --------------------------------------------------------------------------- float versions

def overlap(L: float, M: float, a: float) -> float:
    """Length of the overlap of two centered intervals of lengths L and M, one shifted by a."""
    if L <= 0 or M <= 0:
        raise ValueError("box widths must be positive")
    x = abs(a)
    s, d = (L + M) / 2, abs(L - M) / 2
    if x >= s:
        return 0.0
    if x <= d:
        return min(L, M)
    return s - x


def kernel_matrix(widths: Sequence[float], a: float) -> np.ndarray:
    """K(a) for all pairs of widths (vectorised)."""
    w = np.asarray(widths, dtype=float)
    if np.any(w <= 0):
        raise ValueError("box widths must be positive")
    L, M = w[:, None], w[None, :]
    s, d = (L + M) / 2, np.abs(L - M) / 2
    x = abs(float(a))
    ov = np.where(x <= d, np.minimum(L, M), np.where(x >= s, 0.0, s - x))
    return ov / (L * M)


def weighted_integral_matrix(widths: Sequence[float], upper: float = 1.0) -> np.ndarray:
    """Matrix of int_0^upper a K_LM(a) da, computed exactly on the binary floats then rounded."""
    W = [Fraction(float(x)) for x in widths]
    U = Fraction(float(upper))
    return np.array([[float(weighted_integral_exact(p, q, U)) for q in W] for p in W])


def window_matrix(widths: Sequence[float]) -> np.ndarray:
    """B with q^T B q = R(0) + 2 int_0^1 a R(a) da  (the Montgomery window functional)."""
    return kernel_matrix(widths, 0.0) + 2.0 * weighted_integral_matrix(widths, 1.0)


def knots(widths: Sequence[float], lo: float, hi: float) -> np.ndarray:
    """Sorted breakpoints of all kernels strictly inside (lo, hi), plus lo and hi."""
    if not hi > lo:
        raise ValueError("need lo < hi")
    w = np.asarray(widths, dtype=float)
    p = np.concatenate([(np.abs(w[:, None] - w[None, :]) / 2).ravel(),
                        ((w[:, None] + w[None, :]) / 2).ravel()])
    p = p[(p > lo) & (p < hi)]
    return np.unique(np.concatenate([[lo, hi], p]))


def as_gram(q_or_X) -> np.ndarray:
    """Accept a coefficient vector q (rank one) or a Gram matrix X; return X."""
    arr = np.asarray(q_or_X, dtype=float)
    if arr.ndim == 1:
        return np.outer(arr, arr)
    if arr.ndim == 2 and arr.shape[0] == arr.shape[1]:
        return (arr + arr.T) / 2
    raise ValueError("expected a vector q or a square Gram matrix X")


def correlation(widths: Sequence[float], q_or_X, a: float) -> float:
    """R(a) = tr(K(a) X)."""
    return float(np.sum(kernel_matrix(widths, a) * as_gram(q_or_X)))


# --------------------------------------------------------------------------- exact versions

def to_fraction(x) -> Fraction:
    """Exact conversion. Strings such as '1/6' are parsed exactly; floats are taken at their binary value."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, (int, str)):
        return Fraction(x)
    return Fraction(float(x))


def overlap_exact(L: Fraction, M: Fraction, a: Fraction) -> Fraction:
    x = abs(a)
    s, d = (L + M) / 2, abs(L - M) / 2
    if x >= s:
        return Fraction(0)
    if x <= d:
        return min(L, M)
    return s - x


def kernel_exact(L: Fraction, M: Fraction, a: Fraction) -> Fraction:
    return overlap_exact(L, M, a) / (L * M)


def weighted_integral_exact(L: Fraction, M: Fraction, upper: Fraction = Fraction(1)) -> Fraction:
    """int_0^upper a K_LM(a) da, by exact antiderivatives on the two linear pieces."""
    if upper <= 0:
        return Fraction(0)
    s, d, m = (L + M) / 2, abs(L - M) / 2, min(L, M)
    total = Fraction(0)
    a1 = min(d, upper)                       # flat piece on [0, a1], value m
    total += m * a1 * a1 / 2
    lo, hi = min(d, upper), min(s, upper)    # linear piece s - a on [lo, hi]
    if hi > lo:
        def prim(t):
            return s * t * t / 2 - t ** 3 / 3
        total += prim(hi) - prim(lo)
    return total / (L * M)


def knots_exact(widths: Iterable[Fraction], lo: Fraction, hi: Fraction) -> list[Fraction]:
    W = list(widths)
    pts = {lo, hi}
    for L in W:
        for M in W:
            for p in (abs(L - M) / 2, (L + M) / 2):
                if lo < p < hi:
                    pts.add(p)
    return sorted(pts)
