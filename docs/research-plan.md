# Research plan (persistent)

Last updated: 2026-11-23 (v0.1). Update this file when a question is resolved or added.

## Current position

* Best RH-only value here: P = 0.677634962 (`pc-provencap-A2`), `proof-pending`, relies on the
  unreviewed lemma `lemma-intF-8-3`. Pair-correlation benchmark under RH: 0.6792 (Chirre-Gonçalves-de Laat 2020; not the best known value by all methods).
* The centered-box family cannot use the sign-only relaxation (`neg-signonly-degeneracy`).
* All caps-mode values above 0.68 depend on unproven hypotheses.

## Research B: reproducing the CGdL pair-correlation benchmark (in progress)

* B1a (`src/riemann_research/cgdl_fourier.py`): exact Gaussian-polynomial Fourier operator, scaling and
  normalisation, sign-radius and sign-violation diagnostics with rigor labels, and the functional Z by two
  routes. Infrastructure only; no new result and no registry change.
* Open before B2: (a) the degree bookkeeping in CGdL Section 4.1 (p of degree 2d versus SOS terms in v(x^2));
  (b) the formulation: Z_inv(f) = (f(0) r + (2/r) int_0^r x f) / fhat(0) is invariant under f -> c f(x/a), so a
  single common scaling transformation is redundant and one normalisation plus Z_inv covers every common Gaussian
  width, while imposing both f(0) = 1 and fhat(0) = 1 on the fixed-width family is a restriction at fixed degree.
  Decide which one B2 uses. Changing the relative widths of several Gaussians in one family is NOT covered by
  this invariance and remains a genuine design direction;
  (c) numerical stability of Z for high-degree p: individual Gaussian moments are now computed stably (checked to
  order 80), but the sum over a monomial basis with mixed-sign coefficients can still cancel.
* Fourier-side signs of general members are only explored numerically in B1a; certification is B4.
* B2-mini (`src/riemann_research/cgdl_sdp.py`, `experiments/cgdl_reproduction/b2_mini/`): question (b) answered
  numerically at deg_q_max = 3: with both hard normalisations at fixed width (A) the optimum stays at
  1.45255, with one normalisation and Z_inv (B) it is 1.32865 (candidate, not certified). B2-full should use
  formulation B. Question (a) is sidestepped by m / deg_q_max = 2m + 1, not resolved for the paper's notation.
* New open questions: certifying near-active Fourier constraints (B4); one-sided sensitivity in
  rho = R / width at rho = 1 (C0).
* Findings to explain in Phase C (dual and adjoint sensitivity), recorded from B2-mini:
  - Why A's extra constraint q(0) / (Tq)(0) = a0 blocks all improvement at deg_q_max = 3 (V_A stays at
    1.45255 while V_B = 1.32865). Goal: an algebraic explanation via the active constraints and the dual
    price of that constraint, not only the numerical gap.
  - Fourier-side near-contact in B at m = 1: (Tq)(u) has a grid minimum of about 4.9e-7 at u = 1.067
    (exploratory-float), with eigenvalues of order 1e-9 in X_3 and X_4. Identify whether this is a true
    double zero of the optimal Fourier polynomial and how the contact point moves with rho and degree.
  - Feasible improvement directions: which perturbations of the optimum keep all constraints satisfied.

## Open mathematical claims awaiting human review

| Claim | Where | Needed |
|---|---|---|
| Deduction from certificate to `N_s ≥ (P − o(1)) N` | `derivations/what-a-certificate-establishes.md` | review of steps C1-C7 |
| `∫_{b−1/2}^{b+1/2} F ≤ 8/3 + o(1)` under RH | `derivations/int-F-cap.md` | review |
| Sign-only degeneracy for box mixtures | `derivations/sign-only-degeneracy.md` | review |
| Refinement strengthens cap hypotheses | `derivations/refinement.md` | review |

## Next research questions (in priority order)

1. Get the four arguments above reviewed by a human mathematician; record reviews in the registry.
2. Test-function families whose Fourier side is not compactly supported, to approach 0.6792 with
   certificates (needs validated interval arithmetic instead of rationals).
3. Replace the 8/3 cap with sharper published bounds (Carneiro-Chandee-Chirre-Milinovich), with
   the exact statement recorded as a cited assumption.
4. Dual certificates: lower bounds on the optimum within a family, so "optimal in family" can be certified.
5. Adjoint and duality methods (planned directories under `experiments/`).

## Not planned for now

Autonomous research agents, multi-agent orchestration, large-scale compute infrastructure.
