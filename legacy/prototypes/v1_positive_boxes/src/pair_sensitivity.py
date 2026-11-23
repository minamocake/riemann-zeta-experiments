"""Pair-correlation band sensitivity, proof-controlled compact-support test family.

Research prototype; the out-of-window interval-mass hypotheses used in the
experiments are NOT established zeta-theoretic estimates.

For RH-conditional Montgomery pair correlation and a nonnegative f(x)=|H(x)|^2
with Fourier transform R_h(alpha)=int h(t) h(t+alpha) dt, normalized f(0)=1,
one obtains a conditional upper bound on the multiplicity sum with coefficient

    C = R_h(0) + 2 int_0^1 alpha R_h(alpha) d alpha
        + 2 sum_i U_i max(0,sup_{alpha in B_i} R_h(alpha)).

Here U_i is a *hypothetical or independently proven* asymptotic upper bound
on int_{B_i} F(alpha,T) d alpha; the function h is compactly supported in an
interval of length A and may be signed. In this baseline version we restrict
h to a nonnegative mixture of centered uniform boxes. This restriction makes
R_h(alpha) nonnegative and nonincreasing for alpha>=0, so
max_{B_i} R_h(alpha) = R_h(left_i) EXACTLY, without checking a grid.

The objective is a finite-dimensional quadratic program over the probability
simplex, solved numerically via SLSQP with restarts (not a global certificate).
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import numpy as np
from scipy.optimize import minimize


def overlap_kernel(L: float, M: float, shift: float) -> float:
    """Exact correlation of two centered uniform densities of widths L,M."""
    if L <= 0 or M <= 0:
        raise ValueError('Box widths must be positive')
    x = abs(float(shift))
    half_difference = abs(L - M) / 2
    half_sum = (L + M) / 2
    if x <= half_difference:
        overlap = min(L, M)
    elif x >= half_sum:
        overlap = 0.0
    else:
        overlap = half_sum - x
    return overlap / (L * M)


def weighted_kernel_integral(L: float, M: float, upper: float = 1.0) -> float:
    """Exact integral_0^upper x K_LM(x) dx, piecewise-linear integration."""
    if upper <= 0:
        return 0.0
    endpoints = [0.0, upper, abs(L-M)/2, (L+M)/2]
    endpoints = sorted(set(np.clip(endpoints, 0.0, upper)))
    value = 0.0
    for a, b in zip(endpoints[:-1], endpoints[1:]):
        if b <= a:
            continue
        # Two-point Gauss-Legendre is exact: x times a linear kernel is quadratic.
        midpoint, half_width = (a+b)/2, (b-a)/2
        x1, x2 = midpoint - half_width/np.sqrt(3), midpoint + half_width/np.sqrt(3)
        value += half_width * (x1 * overlap_kernel(L,M,x1)
                             + x2 * overlap_kernel(L,M,x2))
    return float(value)


@dataclass(frozen=True)
class Model:
    A: float
    widths: np.ndarray
    K0: np.ndarray
    known: np.ndarray

    def K(self, alpha: float) -> np.ndarray:
        w = self.widths
        return np.asarray([[overlap_kernel(a,b,alpha) for b in w] for a in w])

    def R(self, q: np.ndarray, alpha: float) -> float:
        return float(q @ self.K(alpha) @ q)


def make_model(A: float, basis_count: int = 12) -> Model:
    if A < 1 or basis_count < 2:
        raise ValueError('Require A>=1 and basis_count>=2')
    widths = np.linspace(A/basis_count, A, basis_count)
    K0 = np.asarray([[overlap_kernel(a,b,0) for b in widths] for a in widths])
    known = np.asarray([[weighted_kernel_integral(a,b) for b in widths] for a in widths])
    return Model(A, widths, K0, known)


def uniform_edges(start: float = 1.0, end: float = 2.0, count: int = 10) -> np.ndarray:
    return np.linspace(start, end, count + 1)


def coefficient_matrix(model: Model, edges: Sequence[float], upper_masses: Sequence[float]) -> np.ndarray:
    edges, masses = np.asarray(edges, float), np.asarray(upper_masses, float)
    if len(edges) != len(masses)+1 or np.any(np.diff(edges)<=0):
        raise ValueError('Strictly increasing edges and exactly one mass per bin required')
    if np.any(edges < 1) or np.any(masses < 0):
        raise ValueError('Only bins >= 1 and nonnegative mass upper bounds allowed')
    Q = model.K0 + 2 * model.known
    # For nonnegative symmetric unimodal h, R is monotone decreasing on alpha>=0.
    for left, U in zip(edges[:-1], masses):
        Q = Q + 2 * U * model.K(left)
    return (Q + Q.T) / 2


def optimize(model: Model, edges: Sequence[float], masses: Sequence[float],
             seed: int = 42, restarts: int = 5):
    Q = coefficient_matrix(model, edges, masses)
    n = len(model.widths)
    starts = [np.eye(n)[-1], np.full(n,1/n)]
    rng = np.random.default_rng(seed)
    for _ in range(max(0,restarts-2)):
        starts.append(rng.dirichlet(np.ones(n)))
    solutions = []
    for q0 in starts:
        opt = minimize(lambda q: float(q @ Q @ q), q0, jac=lambda q: 2*Q@q,
                       bounds=[(0,1)]*n,
                       constraints=[{'type':'eq', 'fun':lambda q: float(q.sum()-1),
                                     'jac':lambda q: np.ones(n)}],
                       method='SLSQP',
                       options={'ftol':1e-12,'maxiter':800})
        if np.max(-opt.x) > 1e-7 or abs(sum(opt.x)-1)>1e-7:
            continue
        solutions.append((float(opt.x@Q@opt.x), opt))
    if not solutions:
        raise RuntimeError('No feasible numerical solution')
    obj, result = min(solutions, key=lambda pair: pair[0])
    q = np.maximum(result.x, 0.0)
    q /= q.sum()
    return {'C':float(q@Q@q), 'P':float(2-q@Q@q),
            'q':q, 'optimizer_success': bool(result.success),
            'optimizer_message':str(result.message),
            'restart_spread':float(max(x[0] for x in solutions)-min(x[0] for x in solutions)),
            'num_feasible_starts':len(solutions),
            'min_tangent_eigenvalue':float(np.linalg.eigvalsh((np.eye(n)[:,:-1]-np.eye(n)[:,-1:]) .T @ Q @ (np.eye(n)[:,:-1]-np.eye(n)[:,-1:]))[0])}


def local_sensitivities(model: Model, q: np.ndarray, edges: Sequence[float]) -> np.ndarray:
    """Analytic derivatives dC/dU_i for a fixed optimizer; envelope at smooth optima."""
    return np.array([2*model.R(q, left) for left in np.asarray(edges)[:-1]])


def finite_difference_sensitivity(model, edges, masses, bin_idx: int, step=1e-5, seed=42):
    masses = np.asarray(masses, float)
    a,b = masses.copy(), masses.copy()
    a[bin_idx] += step
    b[bin_idx] -= step
    if b[bin_idx] < -1e-14:
        raise ValueError('Mass upper bound too small for central finite difference')
    forward = optimize(model, edges, a,seed=seed)
    backward = optimize(model, edges, b,seed=seed)
    return (forward['C']-backward['C'])/(2*step)


def refine_edges(edges: Sequence[float], selected: Sequence[int]) -> np.ndarray:
    edges = np.asarray(edges, float)
    idx = set(selected)
    if not idx.issubset(set(range(len(edges)-1))):
        raise ValueError('Selected bin outside edges')
    out=[edges[0]]
    for i,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
        if i in idx:
            out.append((a+b)/2)
        out.append(b)
    return np.array(out)


def adaptive_scan(A: float, refinement_rounds: int = 2, top_k: int = 3,
                  basis_count: int = 12, seed: int = 42):
    """Counterfactual experiment: U_i = width(B_i), NOT a proved bound for zeta."""
    model = make_model(A,basis_count=basis_count)
    edges = uniform_edges()
    stages = []
    for stage in range(refinement_rounds+1):
        upper_masses = np.diff(edges)  # Counterfactual Poisson-density upper mass.
        fit = optimize(model, edges, upper_masses, seed=seed, restarts=6)
        sens = local_sensitivities(model, fit['q'], edges)
        ranked = [int(i) for i in np.argsort(-sens,kind='stable') if sens[i] > 1e-10]
        stages.append({'stage':stage,'A':A,'edges':edges.copy(),
                       'upper_masses':upper_masses.copy(),
                       'C':fit['C'],'P':fit['P'],'q':fit['q'].copy(),
                       'sensitivity':sens.copy(),
                       'top_bins':ranked[:top_k],
                       'optimizer_success':fit['optimizer_success'],
                       'optimizer_message':fit['optimizer_message'],
                       'restart_spread':fit['restart_spread'],
                       'min_tangent_eigenvalue':fit['min_tangent_eigenvalue']})
        if stage < refinement_rounds:
            edges=refine_edges(edges,ranked[:top_k])
    return stages


def exact_band_extrema(model: Model, q: np.ndarray, left: float, right: float):
    """Exact extrema of R_h on a closed band, even if coefficients q are signed.

    The cross-correlations of centered rectangles are piecewise linear with
    kink positions |L-M|/2 and (L+M)/2. Therefore, evaluating all those kinks
    plus endpoints evaluates the entire continuum, not a discretized sample.
    """
    if right <= left or left < 0:
        raise ValueError('Require 0<=left<right')
    w = model.widths
    knots = {float(left),float(right)}
    for L in w:
        for M in w:
            for x in (abs(L-M)/2, (L+M)/2):
                if left < x < right:
                    knots.add(float(x))
    points = np.array(sorted(knots))
    values = np.array([model.R(np.asarray(q), x) for x in points])
    return {'min':float(values.min()),'max':float(values.max()),
            'argmin':float(points[np.argmin(values)]),
            'argmax':float(points[np.argmax(values)]),
            'knots':points,'values':values}
