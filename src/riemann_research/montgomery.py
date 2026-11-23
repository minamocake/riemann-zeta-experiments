"""Montgomery's pair-correlation inequality for multiplicities, in the centered-box family.

Mathematical setting (all statements conditional on RH)
------------------------------------------------------
Let F(a) = F(a, T) be Montgomery's function. Facts used (cited, not re-proved here):
  * F is real, even and F >= 0. Nonnegativity holds because the weight w(u) = 4/(4+u^2) is the
    Fourier transform of exp(-2|t|) >= 0, so F is an average of squared absolute values. (Being the
    Fourier transform of a nonnegative measure alone would NOT imply pointwise nonnegativity.)
  * F(a) = (1 + o(1)) T^{-2|a|} log T + |a| + o(1) uniformly for |a| <= 1
    (Montgomery 1973; uniformity up to |a| = 1 by Goldston-Montgomery).
For f >= 0 with f(0) = 1 and Fourier transform R supported in [-S, S]:

    N*(T) <= (1 + o(1)) N(T) * [ R(0) + 2 int_0^1 a R(a) da + 2 int_1^S F(a) R(a) da ]

and the proportion of simple zeros satisfies N_s >= 2 N - N*.

The outside term is controlled by a partition of [1, S] into bins B_i with caps U_i such that
int_{B_i} F <= U_i + o(1):

    2 int_1^S F R <= 2 sum_i U_i * max(0, sup_{B_i} R).

Whether a cap U_i is a theorem or an unproven hypothesis is NOT decided here; it is recorded
separately in the result registry (see riemann_research.status).

Audit fix (2026-11-23)
----------------------
The legacy code accepted partitions that did not cover [1, S]: bins starting above 1 or ending
below the support S were silently ignored, which yields invalid (too optimistic) bounds.
`Partition` now rejects such input.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .kernels import as_gram, kernel_matrix, knots, window_matrix

# Coverage of [1, S] is checked exactly (no floating tolerance): a tolerance would let a bin
# partition omit a sliver of [1, S] (2026-11-23 independent audit, finding F2).


class PartitionError(ValueError):
    """Raised when a bin partition cannot support a valid outside-window bound."""


@dataclass(frozen=True)
class Partition:
    """Bins [e_0, e_1], ..., [e_{m-1}, e_m] with caps U_i, validated against the support S."""
    edges: np.ndarray
    caps: np.ndarray
    support: float

    @staticmethod
    def build(edges: Sequence[float], caps: Sequence[float], support: float) -> "Partition":
        e = np.asarray(edges, dtype=float)
        u = np.asarray(caps, dtype=float)
        support = float(support)
        if not np.isfinite(support) or support <= 0:
            raise PartitionError(f"support must be finite and positive (got {support!r})")
        if e.ndim != 1 or e.size < 2:
            raise PartitionError("need at least two edges")
        if u.shape != (e.size - 1,):
            raise PartitionError("exactly one cap per bin is required")
        if not np.all(np.isfinite(e)) or not np.all(np.isfinite(u)):
            raise PartitionError("edges and caps must be finite")
        if np.any(np.diff(e) <= 0):
            raise PartitionError("edges must be strictly increasing")
        if np.any(u < 0):
            raise PartitionError("caps must be nonnegative")
        if e[0] != 1.0:
            raise PartitionError(f"partition must start at 1 (got {e[0]!r}); "
                                 "the interval [1, e_0) would otherwise be ignored")
        if e[-1] < support:
            raise PartitionError(f"partition ends at {e[-1]!r} but R is supported up to {support!r}; "
                                 "the interval (e_m, S] would otherwise be ignored")
        return Partition(e, u, support)


def support_of(widths: Sequence[float]) -> float:
    """R vanishes for |a| >= max width."""
    return float(np.max(np.asarray(widths, dtype=float)))


def window_functional(widths: Sequence[float], q_or_X) -> float:
    """C_window = R(0) + 2 int_0^1 a R(a) da."""
    return float(np.sum(window_matrix(widths) * as_gram(q_or_X)))


def band_positive_maxima(widths: Sequence[float], q_or_X, partition: Partition) -> np.ndarray:
    """max(0, sup_{B_i} R) for every bin, evaluated at all knots (exact up to float rounding)."""
    X = as_gram(q_or_X)
    S = partition.support
    out = np.zeros(partition.caps.size)
    for i, (lo, hi) in enumerate(zip(partition.edges[:-1], partition.edges[1:], strict=True)):
        if lo >= S:
            continue  # R == 0 there
        vals = [float(np.sum(kernel_matrix(widths, k) * X)) for k in knots(widths, lo, min(hi, S))]
        out[i] = max(0.0, max(vals))
    return out


def max_positive_part_outside(widths: Sequence[float], q_or_X) -> float:
    """max(0, sup_{[1, S]} R): zero means R <= 0 on [1, infinity) (up to float rounding)."""
    S = support_of(widths)
    if S <= 1.0:
        return 0.0
    X = as_gram(q_or_X)
    return max(0.0, max(float(np.sum(kernel_matrix(widths, k) * X)) for k in knots(widths, 1.0, S)))


def capped_functional(widths: Sequence[float], q_or_X, partition: Partition) -> float:
    """C = C_window + 2 sum_i U_i max(0, sup_{B_i} R)   (float evaluation, not a certificate)."""
    if partition.support < support_of(widths):
        raise PartitionError("partition was validated for a smaller support than these widths")
    return window_functional(widths, q_or_X) + 2.0 * float(
        partition.caps @ band_positive_maxima(widths, q_or_X, partition))


def proportion_from_C(C: float) -> float:
    """Lower bound for the proportion of simple zeros implied by N* <= C N."""
    return 2.0 - C


def sensitivities(widths: Sequence[float], q_or_X, partition: Partition) -> np.ndarray:
    """dC*/dU_i = 2 max(0, sup_{B_i} R) at an optimiser (Danskin; valid where the optimiser is unique)."""
    return 2.0 * band_positive_maxima(widths, q_or_X, partition)
