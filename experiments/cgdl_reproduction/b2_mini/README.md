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
* m = 1, B: (Tq)(u) has a grid minimum of about 4.9e-7 at u = 1.067 (exploratory-float), and X_3, X_4
  each have an eigenvalue of order 1e-9: the Fourier condition is (nearly) active there. Certifying
  such a point (B4) needs a strictly feasible perturbation, as in CGdL Section 4.1.1.
* Residuals: coefficient matching <= 5e-15, normalisations <= 3e-14, complementarity <= 3e-10.
  Independent Z (closed form = quadrature) agrees with the solver to <= 4e-14. The actual last sign
  change of the returned q equals rho = 1 in every run.
* Duals satisfy the homogeneity identity -(sum of normalisation duals) = V (tested), which fixes
  the sign convention for C0.

## Notes for C0 (not implemented here)

* Sensitivity parameter: rho = R / (Gaussian width), using formulation B. CGdL search rho >= 1, so
  at rho = 1 only one-sided sensitivities are meaningful. V(rho) may be non-differentiable where the
  active set changes (near-zero eigenvalues of X_3, X_4 above; the contact point u ~ 1.067 may move).
* The A-B gap is itself a sensitivity question: the dual of A's extra constraint prices it.
* Linearity remark: f_lam + f_(-lam) = 2 f_0 for f_lam = (1 - x^2)(1 + lam x^2) e^(-pi x^2), also after
  Fourier transformation. Linearity only; it does not make both perturbations feasible.
