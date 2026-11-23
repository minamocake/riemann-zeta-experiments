"""Research B1a: Gaussian-polynomial Fourier foundations for the CGdL reproduction.

Each test names its rigor label (see `riemann_research.cgdl_fourier`). Oracles are listed in
tests/ORACLES.md. Passing tests are software checks, not mathematical proofs.
"""
from fractions import Fraction

import mpmath
import pytest

from riemann_research import cgdl_fourier as cf
from riemann_research.cgdl_fourier import PiLaurent as PL

PI_INV = PL.pi_power(-1)


def poly_mp(p, x):
    with mpmath.workdps(40):
        return cf.poly_eval_mp(p, x)


def lam_family(lam):
    """f_lam(x) = (1 - x^2)(1 + lam x^2) exp(-pi x^2)."""
    lam = Fraction(lam)
    return cf.GaussianPolynomial(cf.even_poly([1, lam - 1, -lam]))


# --------------------------------------------------------------------------- B1-A: the operator T

def test_laguerre_formula_equals_differentiation_route_exactly():
    """exact-symbolic: Laguerre closed form (CGdL 4.1) == repeated FT[x^2 h] = -(hhat)''/(4 pi^2)."""
    for k in range(13):
        mono = cf.even_poly([0] * k + [1])
        assert cf.transform_laguerre(mono) == cf.transform_by_differentiation(mono), k


def test_transform_closed_forms_exactly():
    """exact-symbolic: hand-derived transforms of the required examples."""
    T = cf.transform_laguerre
    assert T(cf.even_poly([1])) == cf.even_poly([1])
    assert T(cf.even_poly([0, 1])) == cf.even_poly([PL.pi_power(-1, Fraction(1, 2)), -1])
    assert T(cf.even_poly([0, 0, 1])) == cf.even_poly(
        [PL.pi_power(-2, Fraction(3, 4)), PL.pi_power(-1, -3), 1])
    # (1 - x^2) exp(-pi x^2)  ->  (1 - 1/(2 pi) + xi^2) exp(-pi xi^2)
    assert T(cf.even_poly([1, -1])) == cf.even_poly([1 - PL.pi_power(-1, Fraction(1, 2)), 1])


@pytest.mark.parametrize("coeffs", [[1], [0, 1], [0, 0, 1], [1, -1]])
def test_transform_matches_quadrature(coeffs):
    """numerical-high-precision: closed form versus direct Fourier integration (40 digits)."""
    f = cf.GaussianPolynomial.from_rational(coeffs)
    with mpmath.workdps(40):
        for xi in ["0", "0.3", "1.1", "2.5"]:
            xi = mpmath.mpf(xi)
            num = cf.fhat_quadrature(f.f, xi, dps=40)
            assert abs(num - f.fhat(xi)) < mpmath.mpf("1e-30"), (coeffs, xi)


def test_transform_is_an_involution_on_even_gaussian_polynomials():
    """exact-symbolic: T(T p) = p (Fourier inversion for even functions)."""
    T = cf.transform_laguerre
    for k in range(11):
        mono = cf.even_poly([0] * k + [1])
        assert T(T(mono)) == mono, k
    p = cf.even_poly([Fraction(3, 7), -2, PL.pi_power(1, 5), 0, Fraction(-1, 9), 4])
    assert T(T(p)) == p


# --------------------------------------------------------------------------- B1-B: scaling and normalisation

def test_scaling_rule_against_quadrature():
    """numerical-high-precision: f_a(x) = f(x/a) has fhat_a(xi) = a fhat(a xi)."""
    a = Fraction(3, 2)
    f = cf.GaussianPolynomial.from_rational([1, Fraction(-1, 3), Fraction(1, 5)], scale=a)
    base = cf.GaussianPolynomial.from_rational([1, Fraction(-1, 3), Fraction(1, 5)])
    with mpmath.workdps(40):
        for xi in ["0", "0.4", "1.3"]:
            xi = mpmath.mpf(xi)
            assert abs(f.fhat(xi) - mpmath.mpf(3) / 2 * base.fhat(xi * 3 / 2)) < mpmath.mpf("1e-35")
            assert abs(cf.fhat_quadrature(f.f, xi, dps=40, support=18) - f.fhat(xi)) < mpmath.mpf("1e-30")


def test_gaussian_example_normalisation_radius_and_z():
    """f_raw = (1 - x^2) exp(-pi x^2), normalised by x -> x/a, a = 1/(1 - 1/(2 pi))."""
    f = cf.normalize(cf.even_poly([1, -1]))
    # exact-symbolic: both normalisations hold exactly, and b = 1/a = 1 - 1/(2 pi).
    assert f.value_at_zero() == PL.const(1)
    num, den = f.fourier_at_zero_ratio()
    assert num == den == 1 - PL.pi_power(-1, Fraction(1, 2))
    with mpmath.workdps(30):
        a = 1 / f.inv_scale.to_mp()
        assert abs(a - mpmath.mpf("1.1892797511079254")) < 1e-15
    # exact-rational-sturm: r(p) = 1 exactly, so r(f_a) = a.
    sr = f.sign_radius_undilated()
    assert sr.admissible and sr.exact and sr.t_lo == sr.t_hi == 1
    # exact-symbolic: fhat_a >= 0 (coefficients 1 - 1/(2 pi) > 0 and 1 in xi^2).
    assert cf.fourier_sign_check(f).passed is True
    assert cf.fourier_sign_check(f).rigor == "exact-symbolic"
    # numerical-high-precision: Z by two routes, and the expected value.
    zc, zq = cf.z_closed_form(f, 50), cf.z_quadrature(f, 50)
    assert abs(zc - zq) < mpmath.mpf("1e-40")
    assert abs(zc - mpmath.mpf("1.4525472652249223")) < 1e-15
    assert all(c.passed for c in cf.diagnose(f))


def test_unnormalised_example_fails_normalisation_but_has_same_invariant_z():
    """The undilated f_raw has fhat(0) = 1 - 1/(2 pi) != 1 (exact), yet Z_inv is the same."""
    raw = cf.GaussianPolynomial(cf.even_poly([1, -1]))
    checks = {c.name: c for c in cf.diagnose(raw)}
    assert checks["f(0) = 1"].passed is True
    assert checks["fhat(0) = 1"].passed is False
    assert abs(cf.z_closed_form(raw) - cf.z_closed_form(cf.normalize(raw.p))) < mpmath.mpf("1e-45")


def test_invariant_z_under_amplitude_and_dilation():
    """numerical-high-precision: Z_inv(c f(x/a)) = Z_inv(f), checked by quadrature in x."""
    base = cf.GaussianPolynomial.from_rational([1, -1, Fraction(-1, 4)])
    z0 = cf.z_closed_form(base, 40)
    for c, a in [(Fraction(3), Fraction(1, 2)), (Fraction(2, 7), Fraction(5, 3))]:
        g = cf.GaussianPolynomial(cf.poly_scale(base.p, c), PL.const(1 / a))
        assert abs(cf.z_quadrature(g, 40) - z0) < mpmath.mpf("1e-30")


def test_lowest_degree_ansatz_cannot_meet_both_hard_normalisations():
    """exact-symbolic: for c (R^2 - x^2) exp(-pi x^2), f(0) - fhat(0) = c/(2 pi) for every R."""
    for R in [Fraction(1), Fraction(3, 2), Fraction(7, 5), Fraction(10)]:
        for c in [Fraction(1), Fraction(-2, 3)]:
            gap = cf.lowest_degree_normalisation_gap(R, c)
            assert gap == PL.pi_power(-1, c / 2)
            assert gap.sign() != 0          # so f(0) = fhat(0) = 1 is impossible
    # With one normalisation and Z_inv (dilation allowed) the same shape is usable (R = 1: f_raw).
    assert abs(cf.z_closed_form(cf.GaussianPolynomial(cf.even_poly([1, -1])))
               - mpmath.mpf("1.4525472652249223")) < 1e-15


# --------------------------------------------------------------------------- Gaussian moments (review fix)

def test_lower_incomplete_gamma_has_no_cancellation_at_order_80():
    """Regression (review of PR #1): the finite-sum formula returned 0 for k = 80, z = pi at 60 digits."""
    with mpmath.workdps(60):
        z = mpmath.pi
        value = cf._lower_gamma_int(80, z)
        assert value != 0
        assert abs(value / mpmath.mpf("1.03094590495115e37") - 1) < mpmath.mpf("1e-14")
        # Independent series: gamma(s, z) = z^s e^(-z) sum_n z^n / (s (s+1) ... (s+n)), positive terms.
        s, term, total, n = 81, mpmath.mpf(1) / 81, mpmath.mpf(1) / 81, 0
        while term > total * mpmath.mpf(10) ** -70:
            n += 1
            term = term * z / (s + n)
            total += term
        series = z ** s * mpmath.exp(-z) * total
        assert abs(value / series - 1) < mpmath.mpf("1e-50")


@pytest.mark.parametrize("r", ["0.1", "0.5", "1", "2"])
def test_gaussian_moments_through_order_80_match_quadrature(r):
    """numerical-high-precision: relative error of every moment k = 0..80 against quadrature.

    Includes very small values (r = 0.1, k = 80 is about 1e-238). Quadrature in u = x/r with
    breakpoints accumulating at u = 1, where x^(2k+1) concentrates.
    """
    with mpmath.workdps(60):
        r = mpmath.mpf(r)
        pts = [0] + [1 - mpmath.mpf(2) ** -j for j in range(1, 12)] + [1]
        worst = mpmath.mpf(0)
        for k in range(81):
            exact = cf.gaussian_moment(k, r)
            quad = r ** (2 * k + 2) * mpmath.quad(
                lambda u, k=k: u ** (2 * k + 1) * mpmath.exp(-mpmath.pi * r * r * u * u), pts)
            assert exact > 0
            worst = max(worst, abs(quad / exact - 1))
        assert worst < mpmath.mpf("1e-45"), worst


# --------------------------------------------------------------------------- B1-C: triangle baseline

def test_triangle_baseline_exact_and_numerical():
    """exact: Z(H) = 4/3. numerical-high-precision: Hhat = Fejer kernel, Z(H) by quadrature."""
    assert cf.triangle_z_exact() == Fraction(4, 3)
    with mpmath.workdps(40):
        assert cf.triangle(0) == 1 and cf.triangle_hat(0) == 1
        for xi in ["0", "0.25", "1", "1.7", "3.5"]:
            xi = mpmath.mpf(xi)
            num = 2 * mpmath.quad(lambda x, xi=xi: (1 - x) * mpmath.cos(2 * mpmath.pi * x * xi), [0, 1])
            assert abs(num - cf.triangle_hat(xi)) < mpmath.mpf("1e-35")
        z = 1 + 2 * mpmath.quad(lambda x: x * cf.triangle(x), [0, 1])
        assert abs(z - mpmath.mpf(4) / 3) < mpmath.mpf("1e-35")


# --------------------------------------------------------------------------- B1-D: sign diagnostics

def test_lambda_plus_half_is_a_fourier_space_violation():
    """lam = +0.5: (T p)(2) ~ -3.2077 (polynomial) and fhat(2) ~ -1.1186e-5 (full transform)."""
    f = lam_family(Fraction(1, 2))
    tp = f.fourier_poly()
    assert abs(poly_mp(tp, 2) - mpmath.mpf("-3.2077135983")) < 1e-9
    with mpmath.workdps(30):
        assert abs(f.fhat(2) - mpmath.mpf("-1.11863955e-5")) < 1e-13
        assert abs(f.fhat(2) - poly_mp(tp, 2) * mpmath.exp(-4 * mpmath.pi)) < 1e-25
    lo, hi = cf.poly_enclosure_at(tp, 2)          # exact-symbolic witness
    assert hi < 0
    check = cf.fourier_sign_check(f)
    assert check.passed is False and check.rigor == "exact-symbolic"
    sr = f.sign_radius_undilated()                # real space is fine: r = 1
    assert sr.admissible and sr.exact and sr.t_hi == 1


def test_lambda_minus_half_is_inadmissible_in_real_space():
    """lam = -0.5: p(3) = 28 > 0 and f > 0 for large |x|; no finite radius is reported."""
    f = lam_family(Fraction(-1, 2))
    sr = f.sign_radius_undilated()
    assert not sr.admissible and sr.t_lo is None and sr.t_hi is None
    assert sr.witness_x == 3 and sr.witness_value == 28
    with pytest.raises(cf.InadmissibleError):
        sr.radius_mp()
    with pytest.raises(cf.InadmissibleError):
        cf.z_closed_form(f)
    checks = {c.name: c for c in cf.diagnose(f)}
    assert checks["eventually non-positive"].passed is False
    assert "objective" not in {c.category for c in cf.diagnose(f)}


def test_lambda_zero_is_admissible_after_scaling():
    f = cf.normalize(lam_family(0).p)
    assert all(c.passed for c in cf.diagnose(f))


def test_sign_radius_ignores_even_multiplicity_roots():
    """exact-rational-sturm: q(t) = (2 - t)(t - 1)^2 -> r^2 = 2; q(t) = (1 - t)(t - 2)^2 -> r^2 = 1."""
    # (2 - t)(t^2 - 2t + 1) = 2 - 5t + 4t^2 - t^3
    sr = cf.sign_radius(cf.even_poly([2, -5, 4, -1]))
    assert sr.admissible and sr.exact and sr.t_hi == 2
    # (1 - t)(t^2 - 4t + 4) = 4 - 8t + 5t^2 - t^3
    sr = cf.sign_radius(cf.even_poly([4, -8, 5, -1]))
    assert sr.admissible and sr.exact and sr.t_hi == 1


def test_irrational_squared_radius_is_enclosed_not_rounded():
    """exact-rational-sturm: q(t) = 2 - t^2 has r^2 = sqrt(2), enclosed to width < 2^-200."""
    sr = cf.sign_radius(cf.even_poly([2, 0, -1]))
    assert sr.admissible and not sr.exact
    assert sr.t_lo ** 2 < 2 <= sr.t_hi ** 2 and sr.t_hi - sr.t_lo < Fraction(1, 2 ** 200)


def test_irrational_radius_objective_routes_agree():
    """numerical-high-precision: Z_inv with irrational r = 2^(1/4) (q = 2 - t^2) by both routes."""
    f = cf.GaussianPolynomial(cf.even_poly([2, 0, -1]))
    assert abs(cf.z_closed_form(f, 40) - cf.z_quadrature(f, 40)) < mpmath.mpf("1e-30")


def test_undecided_fourier_sign_and_uncertain_tail_are_reported_as_such():
    """p = T(1 - xi^2 + xi^4) has pi-dependent coefficients; fhat > 0 but not termwise."""
    p = cf.transform_laguerre(cf.even_poly([1, -1, 1]))
    f = cf.GaussianPolynomial(p)
    fourier = cf.fourier_sign_check(f)
    assert fourier.passed is None and fourier.rigor == "exploratory-float"
    tail = [c for c in cf.diagnose(f) if c.category == "tail"]
    assert len(tail) == 1 and tail[0].passed is None and tail[0].rigor == "exploratory-float"
    with pytest.raises(ValueError):
        cf.sign_radius(p)


def test_no_finite_radius_for_pure_gaussian_or_nonpositive_center():
    assert not cf.sign_radius(cf.even_poly([1])).admissible
    assert not cf.sign_radius(cf.even_poly([-1, 1])).admissible
    assert not cf.sign_radius(cf.even_poly([])).admissible


def test_pi_enclosure_is_rigorous_and_decides_signs():
    with mpmath.workdps(40):
        assert cf.PI_LOWER < Fraction(mpmath.nstr(mpmath.pi, 35)) < cf.PI_UPPER
    assert (PL.pi_power(1) - Fraction(314159, 100000)).sign() == 1
    assert (1 - PL.pi_power(-1, Fraction(1, 2))).sign() == 1
    assert (PL.pi_power(1) - Fraction("3.141592653589795")).sign() is None   # too close: undecided


def test_unknown_rigor_label_is_rejected():
    with pytest.raises(ValueError):
        cf.Check("x", "objective", True, "proven", "labels must come from RIGOR_LABELS")


# --------------------------------------------------------------------------- f versus g

def test_pair_correlation_g_is_not_f():
    """numerical-high-precision: ghat(xi) = f(r xi); g >= 0 while ghat <= 0 for |xi| >= 1."""
    f = cf.normalize(cf.even_poly([1, -1]))
    with mpmath.workdps(40):
        r = f.sign_radius_undilated().radius_mp() / f.inv_scale.to_mp()
        g = cf.pair_correlation_g(f.fhat, r)
        for xi in ["0", "0.6", "1", "1.4"]:
            xi = mpmath.mpf(xi)
            ghat = cf.fhat_quadrature(g, xi, dps=40, support=20)
            assert abs(ghat - f.f(r * xi)) < mpmath.mpf("1e-30")
        assert g(0) > 0 and f.f(r * mpmath.mpf("1.4")) < 0
        assert g(0) != f.f(0)          # g(0) = fhat(0)/r = 1/a, not f(0) = 1
