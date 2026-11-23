"""Optimisers. Their output is a CANDIDATE test function; validity is decided by `certify`.

* `positive_qp`: q >= 0, sum q = 1. For q >= 0 the correlation R is nonincreasing on [0, inf),
  so sup_{B_i} R = R(left_i) and the capped functional is the quadratic form q^T Q q.
  If Q is positive definite on {sum q = 0} (checked in an ORTHONORMAL basis) the problem is
  convex and a KKT point is the global optimum within the family.
* `sdp`: lift X = q q^T to X PSD (requires cvxpy). Every PSD X is itself an admissible test
  function, so the SDP optimum is a valid member of the framework and a lower bound for the
  rank-one problem. Solver output is floating point and must be certified separately.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
from scipy.optimize import minimize

from .kernels import kernel_matrix, knots, window_matrix
from .montgomery import Partition
from .status import Optimality


def positive_qp_matrix(widths, partition: Partition) -> np.ndarray:
    Q = window_matrix(widths).copy()
    for left, u in zip(partition.edges[:-1], partition.caps, strict=True):
        Q += 2.0 * u * kernel_matrix(widths, left)
    return (Q + Q.T) / 2


def tangent_min_eigenvalue(Q: np.ndarray) -> float:
    n = Q.shape[0]
    Z = np.linalg.qr(np.c_[np.ones(n), np.eye(n)[:, : n - 1]])[0][:, 1:]
    return float(np.linalg.eigvalsh(Z.T @ Q @ Z).min())


def positive_qp(widths: Sequence[float], partition: Partition, restarts: int = 6, seed: int = 0) -> dict:
    Q = positive_qp_matrix(widths, partition)
    n = Q.shape[0]
    rng = np.random.default_rng(seed)
    starts = [np.full(n, 1 / n), np.eye(n)[-1]] + [rng.dirichlet(np.ones(n)) for _ in range(restarts)]
    best = None
    for q0 in starts:
        r = minimize(lambda q: q @ Q @ q, q0, jac=lambda q: 2 * Q @ q, method="SLSQP",
                     bounds=[(0, 1)] * n,
                     constraints=[{"type": "eq", "fun": lambda q: q.sum() - 1, "jac": lambda q: np.ones(n)}],
                     options={"ftol": 1e-14, "maxiter": 2000})
        q = np.maximum(r.x, 0)
        q /= q.sum()
        val = float(q @ Q @ q)
        if best is None or val < best[0]:
            best = (val, q)
    C, q = best
    g = 2 * Q @ q
    supp = q > 1e-9
    mu = float(np.mean(g[supp]))
    kkt = float(max(np.ptp(g[supp]), max(0.0, mu - g[~supp].min()) if (~supp).any() else 0.0))
    lam = tangent_min_eigenvalue(Q)
    opt = Optimality.GLOBAL_CONVEX if (lam > 0 and kkt < 1e-5) else Optimality.LOCAL
    return {"C": C, "P": 2 - C, "q": q, "tangent_min_eigenvalue": lam, "kkt_residual": kkt,
            "optimality": opt}


def sdp(widths: Sequence[float], mode: str, partition: Partition | None = None) -> dict:
    """mode 'sign-only': R <= 0 on [1, S];  mode 'caps': capped functional with the given partition."""
    import cvxpy as cp  # optional dependency

    W = np.asarray(widths, float)
    n = W.size
    S = float(W.max())
    B = window_matrix(W)
    X = cp.Variable((n, n), PSD=True)
    cons = [np.ones(n) @ X @ np.ones(n) == 1]
    if mode == "sign-only":
        if S > 1:
            cons += [cp.trace(kernel_matrix(W, k) @ X) <= 0 for k in knots(W, 1.0, S)]
        obj = cp.trace(B @ X)
    elif mode == "caps":
        if partition is None:
            raise ValueError("caps mode needs a validated Partition")
        t = cp.Variable(partition.caps.size, nonneg=True)
        for j, (lo, hi) in enumerate(zip(partition.edges[:-1], partition.edges[1:], strict=True)):
            if lo < S:
                cons += [t[j] >= cp.trace(kernel_matrix(W, k) @ X) for k in knots(W, lo, min(hi, S))]
        obj = cp.trace(B @ X) + 2 * partition.caps @ t
    else:
        raise ValueError(mode)
    prob = cp.Problem(cp.Minimize(obj), cons)
    prob.solve(solver="CLARABEL")
    Xv = (X.value + X.value.T) / 2
    ev = np.sort(np.linalg.eigvalsh(Xv))[::-1]
    return {"status": prob.status, "objective": float(prob.value), "X": Xv,
            "eigenvalue_ratio": float(ev[1] / ev[0]) if ev[0] > 0 else float("nan")}
