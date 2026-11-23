"""Exact rational certificates for box-family test functions.

A certificate fixes a test function exactly (rational widths, and a Gram matrix given as a
nonnegative rational combination of rank-one terms, hence PSD by construction) together with the
treatment of the outside window. `verify` recomputes everything in exact rational arithmetic:

  * f(0) = 1^T X 1 is normalised to 1 exactly;
  * R is piecewise linear with known knots, so knot-wise checks are checks on the continuum;
  * C = R(0) + 2 int_0^1 a R + (outside term) and P = 2 - C are exact rationals.

What a passing certificate does and does not establish
------------------------------------------------------
It establishes a finite arithmetic inequality. Turning it into "N_s(T) >= (P + o(1)) N(T)"
additionally requires RH, Montgomery's theorem (cited), and, for mode "caps", the cap hypotheses,
whose status is recorded separately in results/assumptions.yaml. It is not a formal proof.
The full chain is spelled out in docs/derivations/what-a-certificate-establishes.md.
An independent second verifier lives in independent_check.py.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np

from .kernels import kernel_exact, knots_exact, to_fraction, weighted_integral_exact

SCHEMA = "riemann-research/box-certificate/v1"
MODES = ("support-in-window", "sign-only", "caps")


class CertificateError(ValueError):
    pass


def _gram_exact(components, n):
    X = [[Fraction(0)] * n for _ in range(n)]
    for wt, v in components:
        for i in range(n):
            if v[i] == 0:
                continue
            for j in range(n):
                X[i][j] += wt * v[i] * v[j]
    return X


def verify(cert: dict) -> dict:
    """Re-check a certificate in exact arithmetic. Raises CertificateError if it does not hold."""
    if cert.get("schema") != SCHEMA:
        raise CertificateError("unknown certificate schema")
    W = [to_fraction(w) for w in cert["widths"]]
    n = len(W)
    if n == 0 or any(w <= 0 for w in W):
        raise CertificateError("widths must be positive")
    comps = [(to_fraction(c["weight"]), [to_fraction(x) for x in c["vector"]]) for c in cert["components"]]
    if any(wt < 0 for wt, _ in comps) or any(len(v) != n for _, v in comps):
        raise CertificateError("weights must be nonnegative and vectors must match the widths")
    s = sum(wt * sum(v) ** 2 for wt, v in comps)
    if s <= 0:
        raise CertificateError("f(0) = 0; cannot normalise")
    comps = [(wt / s, v) for wt, v in comps]
    X = _gram_exact(comps, n)
    S = max(W)

    def R(a):
        return sum(X[i][j] * kernel_exact(W[i], W[j], a) for i in range(n) for j in range(n) if X[i][j])

    r0 = sum(X[i][j] * kernel_exact(W[i], W[j], Fraction(0)) for i in range(n) for j in range(n))
    integral = sum(X[i][j] * weighted_integral_exact(W[i], W[j]) for i in range(n) for j in range(n))
    C_window = r0 + 2 * integral

    out = cert["outside"]
    mode = out["mode"]
    if mode not in MODES:
        raise CertificateError(f"unknown mode {mode!r}")
    outside = Fraction(0)
    max_R_outside = Fraction(0)
    if S > 1:
        kn = knots_exact(W, Fraction(1), S)
        max_R_outside = max(R(k) for k in kn)
    if mode == "support-in-window":
        if S > 1:
            raise CertificateError("support exceeds the Montgomery window [-1, 1]")
    elif mode == "sign-only":
        if max_R_outside > 0:
            raise CertificateError(f"R > 0 somewhere on [1, S]: max {float(max_R_outside):.3e}")
    else:
        edges = [to_fraction(e) for e in out["edges"]]
        caps = [to_fraction(u) for u in out["caps"]]
        if len(caps) != len(edges) - 1 or any(b <= a for a, b in zip(edges[:-1], edges[1:], strict=True)):
            raise CertificateError("invalid partition")
        if edges[0] != 1 or edges[-1] < S or any(u < 0 for u in caps):
            raise CertificateError("partition must cover [1, S] with nonnegative caps")
        for (lo, hi), u in zip(zip(edges[:-1], edges[1:], strict=True), caps, strict=True):
            if lo >= S:
                continue
            m = max(R(k) for k in knots_exact(W, lo, min(hi, S)))
            outside += 2 * u * max(Fraction(0), m)
    C = C_window + outside
    P = 2 - C
    claimed = to_fraction(cert["claimed"]["P_lower_bound"])
    if P < claimed:
        raise CertificateError(f"claimed P >= {float(claimed)} but exact P = {float(P)}")
    return {"C_window": C_window, "outside_term": outside, "C": C, "P": P,
            "max_R_on_1_S": max_R_outside, "support": S, "claimed": claimed}


def verify_file(path) -> dict:
    return verify(json.loads(Path(path).read_text()))


def _rat(x: float, den: int) -> Fraction:
    return Fraction(float(x)).limit_denominator(den)


def build(cert_id: str, widths: Sequence, X: np.ndarray, mode: str, edges=None, caps=None,
          denominator: int = 10**12, digits: int = 10, notes: str = "", provenance: dict | None = None) -> dict:
    """Round a float Gram matrix to an exact certificate and set the claim to the exact P rounded down."""
    X = (np.asarray(X, float) + np.asarray(X, float).T) / 2
    lam, V = np.linalg.eigh(X)
    comps = []
    for k, eig in enumerate(lam):
        if eig > 1e-13 * max(lam.max(), 1e-300):
            comps.append({"weight": str(_rat(eig, denominator)),
                          "vector": [str(_rat(x, denominator)) for x in V[:, k]]})
    cert = {"schema": SCHEMA, "id": cert_id, "family": "centered-boxes",
            "widths": [str(to_fraction(w)) for w in widths], "components": comps,
            "outside": {"mode": mode}, "claimed": {"P_lower_bound": "0"}, "notes": notes}
    if provenance is not None:
        cert["provenance"] = provenance
    if mode == "caps":
        cert["outside"]["edges"] = [str(to_fraction(e)) for e in edges]
        cert["outside"]["caps"] = [str(to_fraction(u)) for u in caps]
    exact = verify(cert)
    scale = 10 ** digits
    cert["claimed"]["P_lower_bound"] = str(Fraction(math.floor(exact["P"] * scale), scale))
    cert["claimed"]["P_lower_bound_decimal"] = float(Fraction(math.floor(exact["P"] * scale), scale))
    verify(cert)
    return cert
