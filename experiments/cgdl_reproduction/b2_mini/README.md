# B2-mini: first small SDP

Code: `src/riemann_research/cgdl_sdp.py`. Tests: `tests/test_cgdl_sdp.py`. Raw output:
`results/b2_mini_results.json` (regenerate with `python experiments/cgdl_reproduction/b2_mini/run_b2_mini.py`).

Statuses: computation — numerical (Clarabel, tolerances 1e-10) with independent re-evaluation by the
B1a routines; mathematical status — none (no certificate, Fourier non-negativity not proved);
assumptions — RH for any statement about zeta zeros. Benchmark (registry `literature`): 0.6792,
i.e. Z = 1.3208. Nothing here reaches it.

## Model

y dimensionless, t = y^2, p(y) = q(t), v(t) = [1, ..., t^m], Gram blocks X_2, X_3, X_4 of size m + 1:
q(t) = (rho^2 - t) s_2(t), (Tq)(t) = s_3(t) + t s_4(t), deg q <= deg_q_max = 2m + 1, rho = 1.

* (A) fixed width a0 = 1/(1 - 1/(2 pi)), q(0) = 1, a0 (Tq)(0) = 1, minimise a0 [q(0) + 2 sum q_k M_k(1)].
* (B) (Tq)(0) = 1, q(0) free, minimise Z_inv = rho q(0) + (2/rho) sum q_k M_k(rho).
* q_B = a0 q_A maps A into B with the same objective, so V_B <= V_A.

## Results (m = 0: deg_q_max = 1; m = 1: deg_q_max = 3)

| | V_A | V_B | V_A - V_B |
|---|---|---|---|
| m = 0 | 1.4525472652249225 | 1.4525472652249220 | 4e-16 |
| m = 1 | 1.4525472652591453 | 1.3286461661014960 | 0.1239 |

* m = 0 has a single feasible point and reproduces the B1a Gaussian (plumbing check).
* m = 1, A: no improvement over m = 0 within solver tolerance. The extra hard normalisation pins
  q(0)/(Tq)(0) = a0, and that constraint is binding.
* m = 1, B: Z_inv = 1.32865 (would correspond to 67.13% if certified; below the 0.6792 benchmark and
  below the registry's best RH-only value). Numerical candidate only.
* m = 1, B: Fourier-side near-contact. The continuous minimum of (Tq)(u) over u >= 0, from its
  stationary points, is about 4.0928e-9 at u = 1.0668272166 (numerical-high-precision). The first
  report used an 8001-point grid and gave 4.8554e-7 at u = 1.067 (exploratory-float); the grid missed
  the bottom of the valley by a factor of about 120. X_3 and X_4 each have an eigenvalue of order 1e-9.
  This is not a positivity certificate and not a proof of an exact double zero. Certifying such a
  point (B4) needs a strictly feasible perturbation, as in CGdL Section 4.1.1.
* In A the Fourier polynomial stays far from contact (minimum 0.8408 at u = 0).
* Residuals: coefficient matching <= 5e-15, normalisations <= 3e-14, complementarity <= 3e-10.
  Independent Z (closed form = quadrature) agrees with the solver to <= 4e-14. The actual last sign
  change of the returned q equals rho = 1 in every run.
* Duals satisfy the homogeneity identity -(sum of normalisation duals) = V (tested), which fixes
  the sign convention for C0.
* KKT diagnostics from an explicit numpy model of the SDP (`kkt_diagnostics`): stationarity residual
  <= 1.4e-10 and primal-dual gap <= 5.3e-11 in absolute value in every run.
* q(0) > 0 is checked explicitly (PSD only gives q(0) >= 0); q(0) = 1.006836 for B, m = 1.

## Notes for C0 (not implemented here)

* Sensitivity parameter: the dimensionless rho (sign radius of the undilated q(y^2) e^(-pi y^2) in
  units of its Gaussian width), using formulation B. It is not the normalised physical sign radius
  R = rho q(0) / (Tq)(0), which is rho q(0) in B: R = 1.006836 for B at m = 1 and R = a0 = 1.189280 at
  m = 0. So rho = 1 is not the boundary R = 1 of CGdL's search, and B is defined for every rho > 0.
  V(rho) may be non-differentiable where the active set changes; compare left and right differences.
* The derivative of a moving Fourier contact point u*(rho) must include the derivative of the
  optimising function itself: from d/du (Tq_rho)(u*) = 0, du*/drho involves dq_rho/drho, not only the
  explicit rho-dependence of the model.
* The dual multiplier of A's extra constraint is a local sensitivity of V_A. It does not by itself
  explain the finite gap V_A - V_B = 0.1239: relaxing the constraint that far can change the active
  set (eigenvalues of X_i reaching zero, the Fourier contact point appearing near u = 1.067).
* Linearity remark: f_lam + f_(-lam) = 2 f_0 for f_lam = (1 - x^2)(1 + lam x^2) e^(-pi x^2), also after
  Fourier transformation. Linearity only; it does not make both perturbations feasible.
