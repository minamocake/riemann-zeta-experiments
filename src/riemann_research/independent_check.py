"""Second, mathematically independent verifier for box-family certificates (exact rationals).

Deliberately shares NO code with kernels.py or certify.py (enforced by tests/test_independent_check.py):

* h is rebuilt as a step function on the partition of the real line by all endpoints +-L_j/2;
* R(a) = int h(t) h(t+a) dt is computed by summing overlaps of the step intervals, not by the
  centered-box kernel formula |L-M|/2, (L+M)/2;
* the knots of R are taken as differences of step endpoints, not from the kernel breakpoints;
* int_0^1 a R(a) da uses Simpson's rule on each linear piece of R (exact for the quadratic a R(a)),
  not the antiderivative formula;
* f(0) = sum_k lam_k (int h_k)^2 is computed from the integral of each step function.

Both verifiers check the claimed bound independently. Agreement of both is required for the
computational status "numerically-certified".
Like certify.verify, a pass establishes only a finite arithmetic statement.
"""
from __future__ import annotations

from fractions import Fraction


def _frac(x) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(str(x))


def _steps(widths, vector):
    """Step-function representation of h = sum v_j (1/L_j) 1_[-L_j/2, L_j/2]: list of (left, right, value)."""
    pts = sorted({p for L in widths for p in (-L / 2, L / 2)})
    out = []
    for a, b in zip(pts[:-1], pts[1:], strict=True):
        mid = (a + b) / 2
        val = sum((v / L for L, v in zip(widths, vector, strict=True) if abs(mid) < L / 2), Fraction(0))
        if val:
            out.append((a, b, val))
    return out


def _overlap(a1, b1, a2, b2) -> Fraction:
    lo, hi = max(a1, a2), min(b1, b2)
    return hi - lo if hi > lo else Fraction(0)


class _Correlation:
    """R(a) = sum_k lam_k int h_k(t) h_k(t + a) dt, combined into one interval Gram matrix."""

    def __init__(self, widths, components):
        pts = sorted({p for L in widths for p in (-L / 2, L / 2)})
        self.cells = list(zip(pts[:-1], pts[1:], strict=True))
        index = {c: i for i, c in enumerate(self.cells)}
        n = len(self.cells)
        self.G = [[Fraction(0)] * n for _ in range(n)]
        self.f0 = Fraction(0)
        for lam, vec in components:
            st = _steps(widths, vec)
            self.f0 += lam * sum(((b - a) * v for a, b, v in st), Fraction(0)) ** 2
            for a1, b1, v1 in st:
                for a2, b2, v2 in st:
                    self.G[index[(a1, b1)]][index[(a2, b2)]] += lam * v1 * v2
        self.endpoints = pts

    def __call__(self, a: Fraction) -> Fraction:
        total = Fraction(0)
        for i, (a1, b1) in enumerate(self.cells):
            row = self.G[i]
            for j, (a2, b2) in enumerate(self.cells):
                if row[j]:
                    total += row[j] * _overlap(a1, b1, a2 - a, b2 - a)
        return total

    def knots(self, lo: Fraction, hi: Fraction):
        diffs = {p - q for p in self.endpoints for q in self.endpoints}
        return sorted({lo, hi} | {d for d in diffs if lo < d < hi})


def check(cert: dict) -> dict:
    """Recompute C, P, the outside term and max R on [1, S] for a certificate dict."""
    W = [_frac(w) for w in cert["widths"]]
    comps = [(_frac(c["weight"]), [_frac(x) for x in c["vector"]]) for c in cert["components"]]
    if any(w <= 0 for w in W) or any(lam < 0 for lam, _ in comps):
        raise ValueError("invalid widths or weights")
    R = _Correlation(W, comps)
    if R.f0 <= 0:
        raise ValueError("f(0) = 0")
    S = max(W)
    # window: R(0) + 2 int_0^1 a R(a) da, Simpson on each linear piece of R
    window_integral = Fraction(0)
    ks = R.knots(Fraction(0), Fraction(1))
    for a, b in zip(ks[:-1], ks[1:], strict=True):
        m = (a + b) / 2
        window_integral += (b - a) / 6 * (a * R(a) + 4 * m * R(m) + b * R(b))
    C_window = (R(Fraction(0)) + 2 * window_integral) / R.f0
    max_R = max(R(k) for k in R.knots(Fraction(1), S)) / R.f0 if S > 1 else Fraction(0)
    out = cert["outside"]
    outside = Fraction(0)
    if out["mode"] == "support-in-window":
        if S > 1:
            raise ValueError("support exceeds [-1, 1]")
    elif out["mode"] == "sign-only":
        if max_R > 0:
            raise ValueError("R > 0 on [1, S]")
    elif out["mode"] == "caps":
        edges = [_frac(e) for e in out["edges"]]
        caps = [_frac(u) for u in out["caps"]]
        if edges[0] != 1 or edges[-1] < S or len(caps) != len(edges) - 1:
            raise ValueError("partition must cover [1, S]")
        for (lo, hi), u in zip(zip(edges[:-1], edges[1:], strict=True), caps, strict=True):
            if hi <= lo or u < 0:
                raise ValueError("invalid bin")
            if lo < S:
                m = max(R(k) for k in R.knots(lo, min(hi, S))) / R.f0
                outside += 2 * u * max(Fraction(0), m)
    else:
        raise ValueError("unknown mode")
    C = C_window + outside
    P = 2 - C
    # The claim is checked here too, so this verifier does not rely on certify.verify for it
    # (2026-11-23 independent audit, finding F3).
    if "claimed" in cert and P < _frac(cert["claimed"]["P_lower_bound"]):
        raise ValueError(f"claimed P >= {cert['claimed']['P_lower_bound']} but exact P = {P}")
    return {"C_window": C_window, "outside_term": outside, "C": C, "P": P, "max_R_on_1_S": max_R}
