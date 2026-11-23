# Pair Correlation Sensitivity: initial numerical experiment

**Date:** 2026-11-23  
**Research status:** executable prototype, tested; not a new number-theoretic theorem or a verified improvement on published bounds.

## Question

Which parts of the previously unproved Fourier-frequency window `1 <= alpha <= 2` would deliver the largest **local change in a simple-zero lower-bound functional**, *if* a corresponding interval-integral upper bound could be proven?

## Mathematically explicit conditional experiment

Take a nonnegative mixture `h` of centered normalized box functions of widths `<= A`, with `∫h=1`. Let `f(y)=|FT(h)(y)|² >= 0`, `R(alpha)=FT(f)(alpha)=∫h(t)h(t+alpha)dt`. Then `R>=0`, is even and nonincreasing for positive alpha, and vanishes for `|alpha|>=A`.

The known-window contribution is exactly

`R(0) + 2 ∫_0^1 alpha R(alpha) d alpha`.

On bands `B_i=[b_i,b_(i+1)]` with the **hypothetical** additional asymptotic upper bounds `limsup_T∫_(B_i) F(alpha,T)dalpha <= U_i`, the unknown contribution is at most `2 Σ U_i R(b_i)`. This last inequality is valid pointwise across the **entire continuum**, because every positive centered-box correlation decreases with alpha.

The candidate upper coefficient for the multiplicity ratio is

`C(q,U) = R(0) + 2∫_0^1 alpha R(alpha) d alpha + 2Σ_i U_i R(b_i)`.

Conditional on RH and on ALL displayed `U_i`, the corresponding simple-zero proportion has a lower bound `P>= 2-C`. This is a **mathematical implication**, not a new proven zeta estimate: the numerical runs below set `U_i=width(B_i)` (the `F=1` conjectural scenario), which is currently unavailable as a proved estimate for each bin.

For fixed `q`, `dC/dU_i=2R(b_i)`. At a smooth optimum it is also the envelope derivative. Numerical derivatives are checked against central differences of **reoptimized** problems.

## Actual executed results

### Established-window sanity check, no assumptions outside |alpha|<=1

With `A=1` and 12 positive centered-box basis functions, minimization gives

- `C_candidate = 1.3275509027565`
- `P_candidate = 0.6724490972435`, or 67.24491%

This closely reproduces the ~67.25% Montgomery–Taylor pair-correlation value. It is **not** a reproduction of Chirre–Gonçalves–de Laat 67.92%, which uses a broader function class permitting negative Fourier weights outside `|alpha|=1`.

### Hypothetical F=1 outside the known window

| Window length A | 10 coarse bins, P | 2 adaptive refinements, P |
|---|---:|---:|
| 1.1 | 72.6363% | 73.0734% |
| 1.5 | 84.7045% | 85.5412% |
| 2.0 | 90.7741% | 91.4430% |

**ALL NUMBERS IN THIS TABLE ARE COUNTERFACTUAL**, predicated on `U_i <= width(B_i)` for each band. They are not known lower bounds for the proportion of simple zeros, and they do not surpass any established record.

For initial 0.1-wide bins, the derivative `dC*/dU_i` is:

| Band | A=1.1 | A=1.5 | A=2.0 |
|---|---:|---:|---:|
| 1.00–1.10 | 0.10980 | 0.35318 | 0.47672 |
| 1.10–1.20 | 0 | 0.25807 | 0.40358 |
| 1.20–1.30 | 0 | 0.17352 | 0.33419 |
| 1.30–1.40 | 0 | 0.10267 | 0.26959 |
| 1.40–1.50 | 0 | 0.04576 | 0.21062 |
| 1.50–1.60 | 0 | 0 | 0.15786 |
| 1.60–1.70 | 0 | 0 | 0.11330 |
| 1.70–1.80 | 0 | 0 | 0.07543 |
| 1.80–1.90 | 0 | 0 | 0.04413 |
| 1.90–2.00 | 0 | 0 | 0.01912 |

The most sensitive bands are nearest 1.0. **This ranking is predetermined by the deliberately constrained family of positive centered boxes** (its autocorrelation is monotone). It is *not* evidence that number theory's true optimal choice must be near 1.0. The adaptive pass confirms the mechanics of refinement but cannot discover non-monotonic sensitivities in this family.

The magnitude of the numerical `P` improvement from refinement also partly reflects stronger hypothetical information: splitting one bin and imposing the mass upper bound on *each half* is stricter than merely imposing the bound on the original whole bin.

### Numerical verification

- `pytest -q`: **11 passed**.
- `A=1.5`, first-bin derivative: envelope formula `0.353184835`, central difference after reoptimization `0.353184837`, absolute error `2.30e-09`.
- Multiple SLSQP random restarts agreed to approximately `1e-12` in the tested cases. Computed tangent Hessian eigenvalues were strictly positive at tested window sizes, consistent with a convex quadratic optimization on the simplex *for these instances*; no machine-checkable global optimality certificate was produced.
- The `overlap_kernel()` evaluates uniform-box overlap exactly. `weighted_kernel_integral()` uses 2-point Gaussian quadrature, exactly integrating the piecewise quadratic integrand between every known breakpoint.
- A specific **signed-box** counterexample is included in the tests: at `alpha=1.0` and `1.1` its correlation is negative, but at an interior point near `1.090909` it becomes positive. The `exact_band_extrema()` method identifies this missed positive spike by enumerating **all** analytic breakpoints, not merely a fixed numerical mesh. This directly tests against the off-grid loophole reported in prior experiments.

## Interpretation and next experiment

What was validated: mathematical construction of a Fourier kernel with continuous support control, exact conditional band contributions, reference-value reproduction, runnable sensitivity/finite-difference validation, and an off-grid sign-change regression test.

What remains unproved: any new actual upper bound for `F(alpha,T)` outside `|alpha|<=1`, superiority to the 2020 SDP result, globally optimal performance across all allowable Fourier functions, and relevance of a monotone sensitivity ranking to the unrestricted problem.

**Next meaningful experiment:** permit signed `h` while maintaining the inherently nonnegative test function `f=|FT(h)|²`, and calculate each bin's `max(0,sup R(alpha))` using `exact_band_extrema()`. A nonmonotone sensitivity distribution would then be possible. Optimize this enlarged family with explicit continuous inequalities, compare against the published 2020 SDP objective, and separate numerical candidates from rigorous certificates.

## References

Chirre, Gonçalves, de Laat, *Pair Correlation Estimates for the Zeros of the Zeta Function via Semidefinite Programming*, Adv. Math. **361** (2020), 106926. https://arxiv.org/html/1810.08843

Carneiro, Milinovich, Ramos, *Fourier optimization and Montgomery's pair correlation conjecture*, arXiv:2310.01913, https://arxiv.org/abs/2310.01913

Reused research workflow: `minamocake/adjoint_waveform` (HH forward model, Jacobian, adjoint integration, alignment/sensitivity validation). The original repository was not modified.
