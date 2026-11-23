# B2-full S1: formulation B in the Laguerre basis

Code: `src/riemann_research/cgdl_laguerre.py`. Tests: `tests/test_cgdl_laguerre.py`.
Scripts: `run_s1.py` (bounded degree sweep with a rho search; writes `results/s1_sweep.json`),
`check_published_witness.py` (evaluates this model at CGdL's published d = 40 witness; the witness file is
not stored here; output of a local run against the owner-supplied file is `results/published_witness_check.json`).

Raw outputs regenerated at commit 56e8acc (numerics identical to the pre-commit run).

Statuses: computation — numerical (Clarabel 0.11.1 via CVXPY 1.9.3, float64) with independent re-evaluation in
mpmath; mathematical status — none (no certificate; Fourier non-negativity not proved; the float candidates
violate it by up to ~5e-11 in fhat); assumptions — RH for any statement about zeta zeros. Benchmark in
`results/registry.yaml` (literature): 0.6792, Z = 1.3208. Registry unchanged.

## Model

Formulation B as in B2-mini: q(t) = (rho^2 - t) s_2(t), (Tq)(t) = s_3(t) + t s_4(t), s_i = v^T X_i v, X_i PSD,
(Tq)(0) = 1, minimise Z_inv = rho q(0) + (2/rho) int_0^rho x q(x^2) e^(-pi x^2) dx. deg q <= 2m + 1.
* Output basis w_n(y) = L_n^(-1/2)(y), y = 2 pi t: even Hermite functions, so T w_n = (-1)^n w_n (exact,
  tested against B1a). Multiplication by y is the exact three-term recurrence; t = y / (2 pi).
* Gram basis v_k(t) = L_k^(-1/2)(pi t) ("paper", as in CGdL's ZetaBounds.jl) or L_k^(-1/2)(2 pi t) ("hermite").
  All connection and product coefficients are exact rationals, rounded once to float64.
* rho is dimensionless; the normalised physical sign radius is R = rho q(0) / (Tq)(0). R >= 1 is a diagnostic
  (CGdL Section 4: an exactly admissible function has r >= 1).

## Validation (what was run)

* Exact: T diagonal on w_n (n < 16) via B1a; recurrence matrix; basis round trips; product tensor.
* Objective integrals J_n(rho) against quadrature (n <= 20, relative 1e-13).
* m = 0 and m = 1 reproduce earlier values in both bases (m = 1, rho = 1: 1.3286461659 versus B2-mini
  1.3286461661, difference 3e-10, within solver tolerance).
* Paper and hermite bases agree to 1e-8 for m <= 3. For larger m the hermite basis is badly scaled (product
  coefficients up to 2e17 at m = 20; equality matrix rank-deficient in float64 at m = 15) and is not used.
* Published-witness oracle: our affine maps evaluated at CGdL's d = 40 witness (file sha256 a2685553...cbd9e9)
  give matching residual 7.5e-15, (Tq)(0) - 1 = 4e-16, q(0) = 1 - 9e-16 and Z_inv = 1.3207644973215573,
  i.e. the published value 1.3207644973215569 to 4e-16. This checks the linear maps, the paper basis and the
  objective integrals at full size (m = 40) against independent data.
* Per run: residuals, KKT stationarity and duality gap from the explicit affine data, Gram eigenvalues, q(0),
  R; at the best rho, `highprec_check` re-evaluates the candidate from its Gram matrices alone with B1a's
  monomial transform at 200 digits (matching defect as a function, Z by moments, fhat minima on a grid).

## Results (formulation B, paper basis; rho searched in a bounded bracket)

| m | deg q | rho* | V_B | R | Clarabel settings | residual | KKT | gap | min fhat (grid) |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 3 | 1.0664 | 1.3267511514 | 1.0080416 | strict | 6e-15 | 2e-11 | -6e-12 | +1.6e-18 |
| 2 | 5 | 0.9021 | 1.3234265282 | 1.0058574 | strict | 6e-16 | 5e-11 | 5e-12 | +1.3e-16 |
| 3 | 7 | 1.1582 | 1.3216712411 | 1.0046500 | strict | 3e-16 | 5e-12 | 9e-13 | +4.9e-15 |
| 5 | 11 | 0.8823 | 1.3212327196 | 1.0043676 | strict | 2e-13 | 2e-11 | 5e-13 | +6.7e-12 |
| 8 | 17 | 1.1188 | 1.3211448822 | 1.0042734 | reg1e-7-noequil | 3e-16 | 2e-11 | 2e-12 | +2.8e-11 |
| 10 | 21 | 1.0937 | 1.3209248519 | 1.0040663 | strict | 1e-14 | 2e-11 | -2e-13 | +4.1e-11 |
| 12 | 25 | 1.1117 | 1.3209175301 | 1.0040490 | strict | 3e-15 | 4e-11 | -1e-12 | -5.2e-11 |
| 15 | 31 | 1.0394 | 1.3209140863 | 1.0040445 | reg1e-7-noequil | 1e-15 | 4e-11 | -1e-12 | -2.9e-11 |
| 20 | 41 | 1.2310 | 1.3208280402 | 1.0039433 | reg1e-7-noequil | 1e-15 | 5e-11 | -8e-12 | -2.6e-11 |
| 25 | 51 | 1.2497 | 1.3208128207 | 1.0039144 | reg1e-7-noequil | 1e-13 | 5e-11 | -6e-12 | +1.3e-11 |
| 30 | 61 | 1.1853 | 1.3207724006 | 1.0038961 | reg1e-6 | 2e-12 | 5e-12 | -8e-13 | -2.8e-11 |
| 40 | 81 | — | not accepted | — | reg1e-6 | | | | |

* V_B is non-increasing in m (nested families), as it must be. At m = 30, V_B = 1.3207724 (2 - Z = 0.6792276):
  the published 1.3208 / 67.92% is reproduced at the stated 4-digit precision. Numerical candidate only.
* Reference: the published witness has Z = 1.3207644973 at d = 40 and was stored after adding a 1e-6 slack to
  the optimum (ZetaBounds.jl `rstore`), so CGdL's optimum is about 1.3207635. At m = 30 we are 8.9e-6 above it.
* m = 40: every rho tried returned OPTIMAL_INACCURATE (values 1.320758-1.320775, e.g. 1.3207634549 at rho = 1
  with R = 1.0038740, Fourier violation -2.7e-8 in fhat). Not accepted; whether B goes below CGdL's value at
  d = 40 (allowed, since B contains CGdL's reduced model at rho = R) is unresolved in float64.
* R converges towards CGdL's R = 1.0038743 as m grows (1.00804 -> 1.00390), while rho* does not converge:
  V_B(rho) becomes very flat in rho at high degree (spread ~1e-5 over rho in [0.9, 1.2] at m = 20-30), and the
  bracket edge is hit for m = 20, 25. The physical radius, not rho, is the stable quantity.
* The Clarabel settings ladder matters: from m = 8 the default settings often stop with NumericalError at
  the first iteration; a larger static regularisation converges. Acceptance is decided from recomputed
  residuals, KKT and gap, not from solver status alone.

## Limitations

* float64 interior point: the d = 40 optimum (Gram eigenvalues ~1e-9, Fourier contacts) is beyond the
  accepted range here; the last digits need multiprecision (CGdL used SDPA-GMP at 3000 bits).
* The candidates are not feasible exactly: fhat of the float candidate dips to about -5e-11 (m >= 12), and
  some Gram matrices have eigenvalues down to -3e-11. Certification would need a strictly feasible
  perturbation and interval arithmetic (B4).
* The rho search is a bounded golden-section search on a possibly non-smooth, very flat function; rho* is not
  resolved for m >= 20.
