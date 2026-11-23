"""Gaussian-polynomial Fourier foundations for reproducing Chirre-Gonçalves-de Laat (Research B1a).

Reference: A. Chirre, F. Gonçalves, D. de Laat, "Pair correlation estimates for the zeros of the
zeta function via semidefinite programming", Adv. Math. 361 (2020) 106926, arXiv:1810.08843v2,
Sections 3.1 and 4.1 ("CGdL" below). This module concerns test functions only; every statement
about zeta zeros that uses it assumes RH (see AGENTS.md).

Conventions
-----------
Fourier transform: fhat(xi) = int f(x) exp(-2 pi i x xi) dx, so exp(-pi x^2) is self-dual.

Class A_LP (CGdL): even, continuous, integrable f with f(0) = fhat(0) = 1, fhat >= 0 on R and
f eventually non-positive. Last sign change: r(f) = inf{r > 0 : f(x) <= 0 for |x| >= r}.
Functional (CGdL Lemma 8; assuming RH, N*(T) <= (Z(f) + o(1)) N(T)):

    Z(f) = r(f) + (2 / r(f)) int_0^{r(f)} x f(x) dx.

The function that enters Montgomery's pair-correlation sum is NOT f but

    g(x) = fhat(x / r(f)) / r(f),      with transform  ghat(xi) = f(r(f) xi)

(see `pair_correlation_g`). Confusing f and g swaps the real-space and Fourier-space sign
conditions.

Family
------
f(x) = p(b x) exp(-pi b^2 x^2), with p an even polynomial and b = 1/a > 0 ("inv_scale"), so f is
p(x) exp(-pi x^2) dilated by a. Coefficients of p (powers x^(2k)) and b are stored exactly in
Q[pi, 1/pi] (`PiLaurent`). Then

    fhat(xi) = a (Tp)(a xi) exp(-pi a^2 xi^2),    T(x^(2k)) = (k! / pi^k) L_k^(-1/2)(pi xi^2)

(Laguerre form: CGdL Section 4.1). `transform_by_differentiation` is an independent second
implementation of T, from FT[x^2 h] = -(1 / (4 pi^2)) (hhat)''.

Normalisation and invariance
----------------------------
If f(0) > 0 and fhat(0) > 0, then F(x) = c f(x / a) with c = 1 / f(0) and a = f(0) / fhat(0)
satisfies F(0) = Fhat(0) = 1, r(F) = a r(f) and

    Z(F) = (f(0) r(f) + (2 / r(f)) int_0^{r(f)} x f(x) dx) / fhat(0)  =:  Z_inv(f).

Z_inv is invariant under f -> c f(x / a) for all c, a > 0, and equals Z on normalised functions.
Consequences recorded for B2 (not used to change any registry entry):
  * the Gaussian width is not an independent design parameter for Z: a dilation changes the
    family member but not Z_inv;
  * imposing BOTH f(0) = 1 and fhat(0) = 1 as hard constraints on the fixed-width family
    p(x) exp(-pi x^2) is a genuine restriction at fixed degree (for example
    c (R^2 - x^2) exp(-pi x^2) can satisfy at most one of them), whereas one normalisation
    together with Z_inv loses nothing.

Rigor labels
------------
exact-symbolic            exact arithmetic in Q[pi, 1/pi]; a sign is decided only when the exact
                          rational enclosure from PI_LOWER < pi < PI_UPPER excludes zero.
exact-rational-sturm      exact rational polynomial algebra (square-free decomposition, Sturm
                          sequences, rational bisection).
numerical-high-precision  mpmath at the stated working precision; agreement is evidence only.
exploratory-float         sampling on a finite grid; never a proof of a sign on an infinite domain.
Passing checks are software checks, not mathematical proofs (docs/status-policy.md). Rigorous
Fourier-side certification of general members of the family is deferred to B4 (interval
arithmetic).

Degree convention (open; settle before B2)
------------------------------------------
CGdL write p as an even polynomial of degree 2d (d = 40) and use four PSD blocks X_1..X_4 of size
d + 1 with s_i(y) = v(y)^T X_i v(y), v a basis of the polynomials of degree <= d, and
f(x) = -(s_1(x^2) + (x^2 - R^2) s_2(x^2)) exp(-pi x^2),
fhat(x) = (s_3(x^2) + x^2 s_4(x^2)) exp(-pi x^2).
Read literally, (x^2 - R^2) s_2(x^2) can have degree 4d + 2 in x rather than 2d. The degree
bookkeeping (and the basis used for v) must be checked against the paper before B2; nothing in
this module depends on it.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Callable, Iterable, Sequence

import mpmath

RIGOR_LABELS = ("exact-symbolic", "exact-rational-sturm", "numerical-high-precision", "exploratory-float")

# Rational enclosure of pi used for exact sign decisions. Source: the decimal expansion
# pi = 3.14159265358979323846... (OEIS A000796). Reproducible rigorous verification:
# `pi_enclosure_machin` computes exact rational bounds from Machin's formula, and the test suite
# checks PI_LOWER < lower bound and upper bound < PI_UPPER in exact rational arithmetic.
PI_LOWER = Fraction("3.14159265358979")
PI_UPPER = Fraction("3.14159265358980")

Rational = int | Fraction


class InadmissibleError(ValueError):
    """Raised when a quantity is requested that is undefined for the given function."""


def _arctan_inverse_bounds(m: int, terms: int) -> tuple[Fraction, Fraction]:
    """Exact rational bounds for arctan(1/m), m >= 2, from the alternating Taylor series.

    The terms (1/m)^(2n+1)/(2n+1) decrease strictly, so consecutive partial sums bracket the value.
    """
    if m < 2 or terms < 1:
        raise ValueError("need m >= 2 and terms >= 1")
    x = Fraction(1, m)
    partial = [Fraction(0)]
    for n in range(terms + 1):
        partial.append(partial[-1] + (-1) ** n * x ** (2 * n + 1) / (2 * n + 1))
    a, b = partial[-2], partial[-1]
    return min(a, b), max(a, b)


def pi_enclosure_machin(terms: int = 30) -> tuple[Fraction, Fraction]:
    """Rigorous rational (lo, hi) with lo < pi < hi from Machin's formula (1706):

        pi = 16 arctan(1/5) - 4 arctan(1/239).

    Exact rational arithmetic only; no floating point is involved.
    """
    a_lo, a_hi = _arctan_inverse_bounds(5, terms)
    b_lo, b_hi = _arctan_inverse_bounds(239, terms)
    return 16 * a_lo - 4 * b_hi, 16 * a_hi - 4 * b_lo


# =========================================================================== exact Q[pi, 1/pi]

@dataclass(frozen=True)
class PiLaurent:
    """Exact element sum_e c_e pi^e of Q[pi, 1/pi]; `terms` is sorted, without zero coefficients."""

    terms: tuple[tuple[int, Fraction], ...] = ()

    @staticmethod
    def from_dict(d: dict[int, Fraction]) -> PiLaurent:
        return PiLaurent(tuple(sorted((e, Fraction(c)) for e, c in d.items() if c != 0)))

    @staticmethod
    def const(c: Rational) -> PiLaurent:
        return PiLaurent.from_dict({0: Fraction(c)})

    @staticmethod
    def pi_power(e: int, c: Rational = 1) -> PiLaurent:
        return PiLaurent.from_dict({e: Fraction(c)})

    @staticmethod
    def coerce(x: PiLaurent | Rational) -> PiLaurent:
        return x if isinstance(x, PiLaurent) else PiLaurent.const(x)

    def as_dict(self) -> dict[int, Fraction]:
        return dict(self.terms)

    def __add__(self, other: PiLaurent | Rational) -> PiLaurent:
        d = self.as_dict()
        for e, c in PiLaurent.coerce(other).terms:
            d[e] = d.get(e, Fraction(0)) + c
        return PiLaurent.from_dict(d)

    __radd__ = __add__

    def __neg__(self) -> PiLaurent:
        return PiLaurent(tuple((e, -c) for e, c in self.terms))

    def __sub__(self, other: PiLaurent | Rational) -> PiLaurent:
        return self + (-PiLaurent.coerce(other))

    def __rsub__(self, other: PiLaurent | Rational) -> PiLaurent:
        return PiLaurent.coerce(other) - self

    def __mul__(self, other: PiLaurent | Rational) -> PiLaurent:
        o = PiLaurent.coerce(other)
        d: dict[int, Fraction] = {}
        for e1, c1 in self.terms:
            for e2, c2 in o.terms:
                d[e1 + e2] = d.get(e1 + e2, Fraction(0)) + c1 * c2
        return PiLaurent.from_dict(d)

    __rmul__ = __mul__

    def is_zero(self) -> bool:
        return not self.terms

    def is_rational(self) -> bool:
        return all(e == 0 for e, _ in self.terms)

    def rational_value(self) -> Fraction:
        if not self.is_rational():
            raise ValueError("element involves pi; it is not rational")
        return self.terms[0][1] if self.terms else Fraction(0)

    def enclosure(self) -> tuple[Fraction, Fraction]:
        """Exact rational interval containing the value, from PI_LOWER < pi < PI_UPPER."""
        lo = hi = Fraction(0)
        for e, c in self.terms:
            a, b = PI_LOWER ** e, PI_UPPER ** e
            plo, phi = (a, b) if e >= 0 else (b, a)
            lo += c * plo if c > 0 else c * phi
            hi += c * phi if c > 0 else c * plo
        return lo, hi

    def sign(self) -> int | None:
        """+1, -1 or 0 when decided exactly; None when the enclosure contains zero (undecided)."""
        if self.is_zero():
            return 0
        lo, hi = self.enclosure()
        if lo > 0:
            return 1
        if hi < 0:
            return -1
        return None

    def to_mp(self):
        return mpmath.fsum(mpmath.mpf(c.numerator) / c.denominator * mpmath.pi ** e for e, c in self.terms)

    def __float__(self) -> float:
        return float(self.to_mp())


EvenPoly = tuple[PiLaurent, ...]   # coefficient j multiplies x^(2j)


def even_poly(coeffs: Iterable[PiLaurent | Rational]) -> EvenPoly:
    """Even polynomial sum_j coeffs[j] x^(2j), trailing zero coefficients removed."""
    out = [PiLaurent.coerce(c) for c in coeffs]
    while out and out[-1].is_zero():
        out.pop()
    return tuple(out)


def poly_add(p: EvenPoly, q: EvenPoly) -> EvenPoly:
    n = max(len(p), len(q))
    z = PiLaurent()
    return even_poly((p[j] if j < len(p) else z) + (q[j] if j < len(q) else z) for j in range(n))


def poly_scale(p: EvenPoly, c: PiLaurent | Rational) -> EvenPoly:
    return even_poly(c * a for a in p)


def poly_mul(p: EvenPoly, q: EvenPoly) -> EvenPoly:
    if not p or not q:
        return ()
    out = [PiLaurent()] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] = out[i + j] + a * b
    return even_poly(out)


def poly_value_at_zero(p: EvenPoly) -> PiLaurent:
    return p[0] if p else PiLaurent()


def poly_eval_mp(p: EvenPoly, x):
    """p(x) in mpmath at the current working precision."""
    x2 = mpmath.mpf(x) ** 2
    acc = mpmath.mpf(0)
    for c in reversed(p):
        acc = acc * x2 + c.to_mp()
    return acc


def poly_enclosure_at(p: EvenPoly, x: Rational) -> tuple[Fraction, Fraction]:
    """Exact rational enclosure of p(x) at a rational point x (via the pi enclosure)."""
    x2 = Fraction(x) ** 2
    lo = hi = Fraction(0)
    for j, c in enumerate(p):
        clo, chi = c.enclosure()
        w = x2 ** j
        lo += clo * w
        hi += chi * w
    return lo, hi


# =========================================================================== the operator T

def _half_binomial(alpha: Fraction, m: int) -> Fraction:
    out = Fraction(1)
    for i in range(m):
        out = out * (alpha - i) / (i + 1)
    return out


@lru_cache(maxsize=None)
def laguerre_monomial_transform(k: int) -> EvenPoly:
    """T(x^(2k)) = (k!/pi^k) L_k^(-1/2)(pi xi^2) as an even polynomial in xi (CGdL Section 4.1).

    Coefficient of xi^(2j): (-1)^j (k!/j!) binom(k - 1/2, k - j) pi^(j - k).
    """
    if k < 0:
        raise ValueError("k must be >= 0")
    alpha = Fraction(k) - Fraction(1, 2)
    return even_poly(
        PiLaurent.pi_power(j - k, (-1) ** j * Fraction(math.factorial(k), math.factorial(j))
                           * _half_binomial(alpha, k - j))
        for j in range(k + 1))


def transform_laguerre(p: EvenPoly) -> EvenPoly:
    """T p via the Laguerre closed form; (T p)(xi) exp(-pi xi^2) is the transform of p(x) exp(-pi x^2)."""
    out: EvenPoly = ()
    for k, c in enumerate(p):
        out = poly_add(out, poly_scale(laguerre_monomial_transform(k), c))
    return out


@lru_cache(maxsize=None)
def _differentiation_monomial(k: int) -> EvenPoly:
    # q_{k+1}(xi) = -(1/(4 pi^2)) (q_k'' - 4 pi xi q_k' + (4 pi^2 xi^2 - 2 pi) q_k), from
    # FT[x^2 h] = -(1/(4 pi^2)) (hhat)'' and hhat = q_k(xi) exp(-pi xi^2). Full (not only even)
    # coefficients are tracked so that the vanishing of odd powers is checked, not assumed.
    full = {0: PiLaurent.const(1)}
    for _ in range(k):
        new: dict[int, PiLaurent] = {}
        top = max(full) + 2
        for n in range(top + 1):
            a_n = full.get(n, PiLaurent())
            term = (PiLaurent.const((n + 2) * (n + 1)) * full.get(n + 2, PiLaurent())
                    - PiLaurent.pi_power(1, 4 * n) * a_n
                    + PiLaurent.pi_power(2, 4) * full.get(n - 2, PiLaurent())
                    - PiLaurent.pi_power(1, 2) * a_n)
            v = PiLaurent.pi_power(-2, Fraction(-1, 4)) * term
            if not v.is_zero():
                new[n] = v
        full = new
    if any(n % 2 for n in full):
        raise AssertionError("odd powers appeared in the transform of an even function")
    top = max(full) if full else 0
    return even_poly(full.get(2 * j, PiLaurent()) for j in range(top // 2 + 1))


def transform_by_differentiation(p: EvenPoly) -> EvenPoly:
    """T p computed independently of the Laguerre formula (repeated differentiation identity)."""
    out: EvenPoly = ()
    for k, c in enumerate(p):
        out = poly_add(out, poly_scale(_differentiation_monomial(k), c))
    return out


# =========================================================================== exact rational polynomials in t = x^2

RPoly = list[Fraction]   # low degree first


def _rp(coeffs: Iterable[Rational]) -> RPoly:
    out = [Fraction(c) for c in coeffs]
    while out and out[-1] == 0:
        out.pop()
    return out


def _rp_eval(p: RPoly, t: Fraction) -> Fraction:
    acc = Fraction(0)
    for c in reversed(p):
        acc = acc * t + c
    return acc


def _rp_deriv(p: RPoly) -> RPoly:
    return _rp(i * c for i, c in enumerate(p) if i)


def _rp_divmod(a: RPoly, b: RPoly) -> tuple[RPoly, RPoly]:
    if not b:
        raise ZeroDivisionError("polynomial division by zero")
    a = list(a)
    q = [Fraction(0)] * max(len(a) - len(b) + 1, 0)
    while len(a) >= len(b) and a:
        s = len(a) - len(b)
        f = a[-1] / b[-1]
        q[s] = f
        for i, c in enumerate(b):
            a[s + i] -= f * c
        a = _rp(a)
    return _rp(q), a


def _rp_monic(p: RPoly) -> RPoly:
    return [c / p[-1] for c in p] if p else p


def _rp_gcd(a: RPoly, b: RPoly) -> RPoly:
    while b:
        a, b = b, _rp_divmod(a, b)[1]
    return _rp_monic(a)


def _rp_sub(a: RPoly, b: RPoly) -> RPoly:
    n = max(len(a), len(b))
    return _rp((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(n))


def _rp_mul(a: RPoly, b: RPoly) -> RPoly:
    if not a or not b:
        return []
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return _rp(out)


def squarefree_decomposition(p: RPoly) -> list[tuple[RPoly, int]]:
    """Yun's algorithm over Q: monic square-free a_i with p = const * prod a_i^i (a_i == 1 dropped)."""
    p = _rp(p)
    if len(p) <= 1:
        return []
    dp = _rp_deriv(p)
    a0 = _rp_gcd(p, dp)
    b = _rp_divmod(p, a0)[0]
    c = _rp_divmod(dp, a0)[0]
    d = _rp_sub(c, _rp_deriv(b))
    out, i = [], 1
    while len(b) > 1:
        a = _rp_gcd(b, d)
        b = _rp_divmod(b, a)[0]
        c = _rp_divmod(d, a)[0]
        d = _rp_sub(c, _rp_deriv(b))
        if len(a) > 1:
            out.append((_rp_monic(a), i))
        i += 1
    return out


def _sturm_sequence(p: RPoly) -> list[RPoly]:
    seq = [p, _rp_deriv(p)]
    while seq[-1]:
        r = _rp_divmod(seq[-2], seq[-1])[1]
        if not r:
            break
        seq.append([-c for c in r])
    return seq


def _sign_changes(seq: list[RPoly], t: Fraction) -> int:
    vals = [v for v in (_rp_eval(s, t) for s in seq) if v != 0]
    return sum(1 for u, v in zip(vals[:-1], vals[1:], strict=True) if (u > 0) != (v > 0))


def _cauchy_bound(p: RPoly) -> Fraction:
    lead = abs(p[-1])
    return 1 + max((abs(c) / lead for c in p[:-1]), default=Fraction(0))


def _largest_positive_root(p: RPoly, width: Fraction) -> tuple[Fraction, Fraction, bool]:
    """Enclosure (lo, hi] of the largest root in (0, inf) of a square-free p with p(0) != 0.

    Returns (lo, hi, exact); exact means hi is itself the root. Raises if there is no such root.
    """
    seq = _sturm_sequence(p)
    lo, hi = Fraction(0), _cauchy_bound(p)
    if _sign_changes(seq, lo) - _sign_changes(seq, hi) < 1:
        raise ValueError("no positive root")
    if _rp_eval(p, hi) == 0:
        return hi, hi, True
    while hi - lo > width:
        mid = (lo + hi) / 2
        if _sign_changes(seq, mid) - _sign_changes(seq, hi) >= 1:
            lo = mid
        else:
            if _rp_eval(p, mid) == 0:
                return mid, mid, True
            hi = mid
    if _rp_eval(p, hi) == 0:
        return lo, hi, True
    # Detect a rational root with a small denominator that bisection did not land on exactly.
    for bound in (10 ** 3, 10 ** 6, 10 ** 12):
        cand = ((lo + hi) / 2).limit_denominator(bound)
        if lo < cand <= hi and _rp_eval(p, cand) == 0:
            return cand, cand, True
    return lo, hi, False


# =========================================================================== sign radius (real space)

@dataclass(frozen=True)
class SignRadius:
    """Last sign change of p(x) exp(-pi x^2) in the variable of p (before any dilation).

    admissible: f(0) > 0 and f eventually non-positive. If admissible, r^2 lies in
    [t_lo, t_hi] (t_lo == t_hi when exact). If not, `witness_x` is a rational x with p(x) > 0
    beyond every root (when the reason is eventual positivity).
    """

    admissible: bool
    reason: str
    t_lo: Fraction | None = None
    t_hi: Fraction | None = None
    exact: bool = False
    witness_x: Fraction | None = None
    witness_value: Fraction | None = None
    rigor: str = "exact-rational-sturm"

    def radius_mp(self):
        if not self.admissible:
            raise InadmissibleError(f"no finite sign radius: {self.reason}")
        return mpmath.sqrt((mpmath.mpf(self.t_lo.numerator) / self.t_lo.denominator
                            + mpmath.mpf(self.t_hi.numerator) / self.t_hi.denominator) / 2)


def sign_radius(p: EvenPoly, width: Fraction = Fraction(1, 2 ** 240)) -> SignRadius:
    """Exact last sign change of p(x) exp(-pi x^2) for p with rational coefficients.

    With q(t) = p(sqrt(t)): r^2 = sup{t >= 0 : q(t) > 0}, which is the largest positive root of
    odd multiplicity of q (q < 0 beyond its largest root when the leading coefficient is negative).
    """
    if not all(c.is_rational() for c in p):
        raise ValueError("exact sign radius needs rational coefficients; use B4 interval methods")
    q = _rp(c.rational_value() for c in p)
    if not q:
        return SignRadius(False, "f is identically zero")
    if q[0] <= 0:
        return SignRadius(False, "f(0) <= 0, so f cannot be normalised to f(0) = 1")
    if len(q) == 1:
        return SignRadius(False, "f is a positive multiple of exp(-pi x^2): never non-positive")
    if q[-1] > 0:
        x = Fraction(math.isqrt(math.ceil(_cauchy_bound(q))) + 1)
        return SignRadius(False, "f(x) > 0 for all sufficiently large |x| (positive leading coefficient)",
                          witness_x=x, witness_value=_rp_eval(q, x * x))
    odd: RPoly = [Fraction(1)]
    for a, mult in squarefree_decomposition(q):
        if mult % 2:
            odd = _rp_mul(odd, a)
    lo, hi, exact = _largest_positive_root(odd, width)
    return SignRadius(True, "eventually non-positive", lo, hi, exact)


# =========================================================================== the family

@dataclass(frozen=True)
class GaussianPolynomial:
    """f(x) = p(b x) exp(-pi b^2 x^2), b = inv_scale > 0 (dilation by a = 1/b)."""

    p: EvenPoly
    inv_scale: PiLaurent = PiLaurent.const(1)

    def __post_init__(self):
        if self.inv_scale.sign() != 1:
            raise ValueError("inv_scale must be provably positive")

    @staticmethod
    def from_rational(coeffs: Sequence[Rational], scale: Rational = 1) -> GaussianPolynomial:
        """p with rational coefficients of x^(2j), dilated by a rational scale a > 0."""
        return GaussianPolynomial(even_poly(coeffs), PiLaurent.const(1 / Fraction(scale)))

    def fourier_poly(self) -> EvenPoly:
        """T p (transform polynomial in the undilated variable)."""
        return transform_laguerre(self.p)

    def value_at_zero(self) -> PiLaurent:
        """f(0) = p(0), exactly."""
        return poly_value_at_zero(self.p)

    def fourier_at_zero_ratio(self) -> tuple[PiLaurent, PiLaurent]:
        """fhat(0) = N / D exactly, with N = (T p)(0) and D = b."""
        return poly_value_at_zero(self.fourier_poly()), self.inv_scale

    def f(self, x):
        b = self.inv_scale.to_mp()
        u = b * mpmath.mpf(x)
        return poly_eval_mp(self.p, u) * mpmath.exp(-mpmath.pi * u * u)

    def fhat(self, xi):
        b = self.inv_scale.to_mp()
        u = mpmath.mpf(xi) / b
        return poly_eval_mp(self.fourier_poly(), u) * mpmath.exp(-mpmath.pi * u * u) / b

    def sign_radius_undilated(self) -> SignRadius:
        return sign_radius(self.p)


def normalize(p: EvenPoly) -> GaussianPolynomial:
    """Rescale p(x) exp(-pi x^2) to c f(x / a) with f(0) = fhat(0) = 1, exactly.

    c = 1 / p(0) must be rational; a = 1 / (T p_c)(0) with (T p_c)(0) provably positive.
    """
    p0 = poly_value_at_zero(p)
    if not p0.is_rational() or p0.rational_value() <= 0:
        raise InadmissibleError("normalisation needs a rational p(0) > 0")
    pc = poly_scale(p, 1 / p0.rational_value())
    t0 = poly_value_at_zero(transform_laguerre(pc))
    if t0.sign() != 1:
        raise InadmissibleError("fhat(0) is not provably positive")
    return GaussianPolynomial(pc, t0)


def lowest_degree_normalisation_gap(R: Rational, c: Rational = 1) -> PiLaurent:
    """f(0) - fhat(0) for f(x) = c (R^2 - x^2) exp(-pi x^2); equals c / (2 pi) for every R."""
    f = GaussianPolynomial(even_poly([c * Fraction(R) ** 2, -Fraction(c)]))
    num, den = f.fourier_at_zero_ratio()
    assert den == PiLaurent.const(1)
    return f.value_at_zero() - num


# =========================================================================== objective

def _lower_gamma_int(k: int, z):
    """Lower incomplete gamma gamma(k + 1, z) for integer k >= 0 and z >= 0, via mpmath.gammainc.

    The finite-sum identity gamma(k + 1, z) = k! (1 - exp(-z) sum_{m <= k} z^m / m!) is exact but
    numerically useless for large k / small z: the bracket is about z^(k+1) / (k + 1)! and is lost
    to cancellation (k = 80, z = pi gives 0 at 60 digits instead of 1.0309459049511504e37).
    """
    if k < 0:
        raise ValueError("k must be >= 0")
    return mpmath.gammainc(k + 1, 0, z)


def gaussian_moment(k: int, r):
    """int_0^r x^(2k+1) exp(-pi x^2) dx = gamma(k + 1, pi r^2) / (2 pi^(k + 1)).

    Rigor: numerical-high-precision (mpmath at the current working precision).
    """
    r = mpmath.mpf(r)
    return _lower_gamma_int(k, mpmath.pi * r * r) / (2 * mpmath.pi ** (k + 1))


def z_closed_form(f: GaussianPolynomial, dps: int = 50):
    """Z_inv(f) from incomplete-gamma moments; equals Z(f) when f is normalised.

    Uses only the undilated data: Z_inv = (p(0) r0 + (2 / r0) I0) / (T p)(0) with
    I0 = int_0^{r0} x p(x) exp(-pi x^2) dx = sum_k c_k gamma(k + 1, pi r0^2) / (2 pi^(k + 1)).
    Each moment is computed stably (`gaussian_moment`). The sum over k can still cancel when the
    coefficients of p have mixed signs and large size (high degree in the monomial basis); B2
    must choose a basis or precision accordingly. Rigor: numerical-high-precision.
    """
    sr = f.sign_radius_undilated()
    with mpmath.workdps(dps + 10):
        r0 = sr.radius_mp()
        i0 = mpmath.fsum(c.to_mp() * gaussian_moment(k, r0) for k, c in enumerate(f.p))
        num = poly_value_at_zero(f.p).to_mp() * r0 + 2 * i0 / r0
        val = num / poly_value_at_zero(f.fourier_poly()).to_mp()
    with mpmath.workdps(dps):
        return +val


def z_quadrature(f: GaussianPolynomial, dps: int = 50):
    """Z_inv(f) from the definition in the dilated variable, by mpmath quadrature.

    Independent of `z_closed_form` except for the sign radius. Rigor: numerical-high-precision.
    """
    sr = f.sign_radius_undilated()
    with mpmath.workdps(dps + 10):
        r = sr.radius_mp() / f.inv_scale.to_mp()
        integral = mpmath.quad(lambda x: x * f.f(x), mpmath.linspace(0, r, 9))
        val = (f.f(0) * r + 2 * integral / r) / f.fhat(0)
    with mpmath.workdps(dps):
        return +val


def fhat_quadrature(f: Callable, xi, dps: int = 40, support: float | None = None):
    """fhat(xi) = 2 int_0^inf f(x) cos(2 pi x xi) dx for even f, by mpmath quadrature.

    `support` bounds the integration range (the integrand must be negligible beyond it).
    Rigor: numerical-high-precision.
    """
    with mpmath.workdps(dps + 10):
        L = mpmath.mpf(support if support is not None else 12)
        n = max(16, int(4 * abs(float(xi)) * float(L)) + 16)
        val = 2 * mpmath.quad(lambda x: f(x) * mpmath.cos(2 * mpmath.pi * x * xi), mpmath.linspace(0, L, n))
    with mpmath.workdps(dps):
        return +val


# =========================================================================== baseline: triangle

def triangle(x):
    """H(x) = max(1 - |x|, 0). Admissible baseline; not a member of the Gaussian family."""
    x = abs(mpmath.mpf(x))
    return 1 - x if x < 1 else mpmath.mpf(0)


def triangle_hat(xi):
    """Hhat(xi) = (sin(pi xi) / (pi xi))^2 (Fejér kernel)."""
    xi = mpmath.mpf(xi)
    return mpmath.mpf(1) if xi == 0 else (mpmath.sin(mpmath.pi * xi) / (mpmath.pi * xi)) ** 2


def triangle_z_exact() -> Fraction:
    """Z(H) = 1 + 2 int_0^1 x (1 - x) dx = 4/3 exactly (r(H) = 1, H(0) = Hhat(0) = 1)."""
    r = Fraction(1)
    integral = Fraction(1, 2) - Fraction(1, 3)
    return r + 2 * integral / r


# =========================================================================== f versus g

def pair_correlation_g(fhat: Callable, r) -> Callable:
    """g(x) = fhat(x / r) / r: the function used in Montgomery's sum (CGdL Section 3.1).

    For even f, ghat(xi) = f(r xi): g >= 0 comes from fhat >= 0, and ghat <= 0 for |xi| >= 1 comes
    from f <= 0 beyond r. Do not substitute f itself into the pair-correlation sum.
    """
    r = mpmath.mpf(r)
    return lambda x: fhat(mpmath.mpf(x) / r) / r


# =========================================================================== diagnostics

@dataclass(frozen=True)
class Check:
    """One diagnostic. passed is True, False, or None (undetermined)."""

    name: str
    category: str     # normalization | real-space | fourier-space | objective | tail
    passed: bool | None
    rigor: str
    detail: str

    def __post_init__(self):
        if self.rigor not in RIGOR_LABELS:
            raise ValueError(f"unknown rigor label {self.rigor!r}")


def fourier_sign_check(f: GaussianPolynomial, grid_max: float = 6.0, grid_n: int = 601) -> Check:
    """Decide fhat >= 0 when possible; otherwise report undetermined (never infer a sign from samples)."""
    tp = f.fourier_poly()
    signs = [c.sign() for c in tp]
    if tp and all(s in (0, 1) for s in signs):
        return Check("fourier-nonnegative", "fourier-space", True, "exact-symbolic",
                     "all coefficients of T p in xi^2 are provably >= 0, so fhat >= 0 on R")
    if tp and signs[-1] == -1:
        return Check("fourier-nonnegative", "fourier-space", False, "exact-symbolic",
                     "leading coefficient of T p is provably negative: fhat(xi) < 0 for large |xi|")
    with mpmath.workdps(30):
        candidates = [k * grid_max / (grid_n - 1) for k in range(grid_n)]
        negative = [u for u in candidates if poly_eval_mp(tp, u) < 0]
    for u in negative:
        ur = Fraction(u).limit_denominator(10 ** 6)
        lo, hi = poly_enclosure_at(tp, ur)
        if hi < 0:
            return Check("fourier-nonnegative", "fourier-space", False, "exact-symbolic",
                         f"(T p)({ur}) lies in [{float(lo):.10g}, {float(hi):.10g}] < 0 "
                         f"(undilated variable; fhat(a * {ur}) < 0)")
    return Check("fourier-nonnegative", "fourier-space", None, "exploratory-float",
                 f"no negative value of T p found on {grid_n} points in [0, {grid_max}]; "
                 "this is not a proof (deferred to B4 interval arithmetic)")


def diagnose(f: GaussianPolynomial, dps: int = 40) -> list[Check]:
    """Classified checks: normalisation, real-space sign, Fourier sign, tail, objective."""
    out: list[Check] = []
    f0 = f.value_at_zero()
    out.append(Check("f(0) = 1", "normalization", f0 == PiLaurent.const(1), "exact-symbolic",
                     f"f(0) = {float(f0):.15g}"))
    num, den = f.fourier_at_zero_ratio()
    out.append(Check("fhat(0) = 1", "normalization", num == den, "exact-symbolic",
                     f"fhat(0) = {float(num.to_mp() / den.to_mp()):.15g}"))
    if all(c.is_rational() for c in f.p):
        sr = sign_radius(f.p)
        if sr.admissible:
            out.append(Check("eventually non-positive", "real-space", True, sr.rigor,
                             f"r(f)^2 * b^2 in [{float(sr.t_lo):.15g}, {float(sr.t_hi):.15g}]"
                             + (" (exact)" if sr.exact else "")))
        else:
            wit = (f" witness: p({sr.witness_x}) = {sr.witness_value} > 0"
                   if sr.witness_x is not None else "")
            out.append(Check("eventually non-positive", "real-space", False, sr.rigor, sr.reason + wit))
    else:
        out.append(Check("eventually non-positive", "tail", None, "exploratory-float",
                         "p has coefficients involving pi: the tail sign is not decided exactly here"))
    out.append(fourier_sign_check(f))
    real = out[2]
    if real.passed:
        zc, zq = z_closed_form(f, dps), z_quadrature(f, dps)
        agree = abs(zc - zq) <= mpmath.mpf(10) ** (-(dps - 10))
        out.append(Check("Z: closed form = quadrature", "objective", bool(agree), "numerical-high-precision",
                         f"closed form {mpmath.nstr(zc, 20)}, quadrature {mpmath.nstr(zq, 20)}"))
    return out
