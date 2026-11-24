# B2-full S1: formulation B in the Laguerre basis

Code: `src/riemann_research/cgdl_laguerre.py`. Tests: `tests/test_cgdl_laguerre.py`.
Scripts:
* `run_s1.py`: bounded degree sweep with a rho search; writes `results/s1_sweep.json` (strict JSON) and
  `results/s1_witnesses.json` (primal and dual witnesses for m = 1, 10, 30; its SHA-256 is recorded in the sweep).
* `replay_s1.py [--fourier]`: solver-free replay of the saved witnesses (no CVXPY, no Clarabel).
* `check_published_witness.py <Z-40.txt>`: evaluates this model at CGdL's published d = 40 witness. The file is not
  stored here and must have SHA-256 `a268555389798aa132a46d90057f2dfa28bc1dcf81d212282b057349d2cbd9e9`; any other
  file is rejected. `results/published_witness_check.json` is the output of a local run against the owner-supplied
  file, made before the digest was enforced; its recorded digest is the pinned one.

Raw outputs regenerated at commit af6336d on a clean worktree (`environment.commit`, `environment.worktree_clean`).

Statuses: computation — numerical (Clarabel 0.11.1 via CVXPY 1.9.3, float64) with independent re-evaluation in
mpmath; mathematical status — none (no certificate; Fourier non-negativity not proved, and the float candidates
violate it, see "Fourier side"); assumptions — RH for any statement about zeta zeros. Benchmark in
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

## Numerical acceptance (not a certificate)

`acceptable(r, tol)` with tol = 1e-6 in the sweep (written out in `s1_sweep.json` under `acceptance`): solver status
OPTIMAL; every diagnostic present and finite; recomputed equality residual <= tol; |duality gap| <= tol; KKT
stationarity <= tol; minimum eigenvalue of every primal Gram block and of every PSD dual block >= -tol;
q(0) > 1e-8. Diagnostics are recomputed from the explicit affine data, not taken from the solver.
`best_accepted` is the lowest accepted objective among all solves for that degree, or null if none is accepted;
rejected solves appear only in `rho_evaluations` and `lowest_rejected`.

## Validation (what was run)

* Exact: T diagonal on w_n (n < 16) via B1a; recurrence matrix; basis round trips; product tensor.
* Objective integrals J_n(rho) against quadrature (n <= 20, relative 1e-13).
* m = 0 and m = 1 reproduce earlier values in both bases (m = 1, rho = 1: 1.3286461659 versus B2-mini
  1.3286461661, difference 3e-10, within solver tolerance).
* Paper and hermite bases agree to 1e-8 for m <= 3. For larger m the hermite basis is badly scaled (product
  coefficients up to 2e17 at m = 20; equality matrix rank-deficient in float64 at m = 15) and is not used.
* Published-witness oracle: our affine maps evaluated at CGdL's d = 40 witness give matching residual 7.5e-15,
  (Tq)(0) - 1 = 4e-16, q(0) = 1 - 9e-16 and Z_inv = 1.3207644973215573, i.e. the published value
  1.3207644973215569 to 4e-16. This checks the linear maps, the paper basis and the objective integrals at full
  size (m = 40) against independent data. The digest of the original arXiv file has not been re-derived here.
* Per run: residuals, KKT stationarity, duality gap, complementarity, primal and dual Gram eigenvalues, q(0), R.
  At the best accepted rho, `highprec_check` re-evaluates the candidate from its Gram matrices alone with B1a's
  monomial transform at 200 digits (matching defect, Z by moments, Fourier-side minima).
* Solver-free replay (`replay_s1.py --fourier`, run at 6a02741 with cvxpy and clarabel imports blocked): for
  m = 1, 10, 30 every check passes (witness SHA-256, shapes, finiteness, symmetry, acceptance re-applied, agreement
  with the recorded diagnostics to 1e-12, same point as `best_accepted`, Fourier minima reproduced). The tests
  check that corrupting a primal entry, the symmetry, an equality dual, a dual PSD block, a recorded diagnostic,
  rho, the acceptance flag, the array shapes or the witness file makes the replay fail.

## Results (formulation B, paper basis; rho searched in a bounded bracket)

| m | deg q | rho* | V_B | R | Clarabel settings | residual | KKT | gap | min eig X | min eig S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 3 | 1.0664 | 1.3267511514 | 1.0080416 | strict | 6e-15 | 2e-11 | -6e-12 | +4e-12 | +1e-12 |
| 2 | 5 | 0.9021 | 1.3234265282 | 1.0058574 | strict | 6e-16 | 5e-11 | 5e-12 | -2e-11 | +4e-12 |
| 3 | 7 | 1.1582 | 1.3216712411 | 1.0046500 | strict | 3e-16 | 5e-12 | 9e-13 | -2e-12 | +2e-13 |
| 5 | 11 | 0.8823 | 1.3212327196 | 1.0043676 | strict | 2e-13 | 2e-11 | 5e-13 | -1e-11 | +1e-12 |
| 8 | 17 | 1.1188 | 1.3211448822 | 1.0042734 | reg1e-7-noequil | 3e-16 | 2e-11 | 2e-12 | -1e-11 | +7e-12 |
| 10 | 21 | 1.0937 | 1.3209248519 | 1.0040663 | strict | 1e-14 | 2e-11 | -2e-13 | -2e-11 | +4e-12 |
| 12 | 25 | 1.1117 | 1.3209175301 | 1.0040490 | strict | 3e-15 | 4e-11 | -1e-12 | -3e-11 | +3e-11 |
| 15 | 31 | 1.0394 | 1.3209140863 | 1.0040445 | reg1e-7-noequil | 1e-15 | 4e-11 | -1e-12 | -3e-11 | +2e-11 |
| 20 | 41 | 1.2310 | 1.3208280402 | 1.0039433 | reg1e-7-noequil | 1e-15 | 5e-11 | -8e-12 | -3e-11 | +3e-11 |
| 25 | 51 | 1.2497 | 1.3208128207 | 1.0039144 | reg1e-7-noequil | 1e-13 | 5e-11 | -6e-12 | -3e-11 | +3e-11 |
| 30 | 61 | 1.1853 | 1.3207724006 | 1.0038961 | reg1e-6 | 2e-12 | 5e-12 | -8e-13 | -3e-11 | +3e-13 |
| 40 | 81 | — | no accepted solve | — | reg1e-6 | | | | | |

* The true optimum of B is non-increasing in m (nested families, for each rho and hence after minimising over rho).
  The searched values above are non-increasing too; because rho is searched separately and only in a bounded bracket
  for each m, this is a consistency check of the search, not a consequence of the theory. At m = 30,
  V_B = 1.3207724 (2 - Z = 0.6792276): the published 1.3208 / 67.92% is reproduced at the stated 4-digit precision.
  Numerical candidate only.
* Reference: the published witness has Z = 1.3207644973 at d = 40 and was stored after adding a 1e-6 slack to
  the optimum (ZetaBounds.jl `rstore`), so CGdL's optimum is about 1.3207635. At m = 30 we are 8.9e-6 above it.
* m = 40: all six grid points returned OPTIMAL_INACCURATE (values 1.3207582-1.3207650; lowest 1.3207581964 at
  rho = 1.2). None is accepted, so no candidate is reported, no high-precision check is run and the bracket search
  is skipped. Whether B goes below CGdL's value at d = 40 (allowed, since B contains CGdL's reduced model at
  rho = R) is unresolved in float64.
* R converges towards CGdL's R = 1.0038743 as m grows (1.00804 -> 1.00390), while rho* does not converge:
  V_B(rho) becomes very flat in rho at high degree (spread ~1e-5 over rho in [0.9, 1.2] at m = 20-30), and the
  bracket edge is hit for m = 20, 25. The physical radius, not rho, is the stable quantity.
* The Clarabel settings ladder matters: from m = 8 the default settings often stop with NumericalError at
  the first iteration; a larger static regularisation converges. The reported settings are those of the rung
  that produced the point.

## Fourier side of the float candidates (numerical diagnostics, not certificates)

With t = xi^2 and P = Tq, fhat(xi) = P(xi^2) exp(-pi xi^2). Two different minima over the half-line t >= 0 are
reported by `fourier_side_minima` (in `highprec_check`, 200 digits, from the Gram matrices):
* unweighted: min P(t); candidates t = 0 and the positive roots of P'(t) = 0; P -> +inf in the tail here (positive
  leading coefficient for every m below);
* weighted: min P(t) exp(-pi t), i.e. the minimum of fhat; candidates t = 0 and the positive roots of
  P'(t) - pi P(t) = 0; the weighted function tends to 0 in the tail, so its infimum is min(candidate values, 0).
The two have different stationary points and are not interchangeable.

| m | min P (unweighted) | at xi | inf of P e^(-pi t) (= min fhat) | at xi | reliable |
|---|---|---|---|---|---|
| 1 | +2.5e-10 | 0.995 | 0 (not attained; all values > 0) | tail | yes |
| 2 | -6.8e-10 | 1.168 | -9.4e-12 | 1.168 | yes |
| 3 | -2.3e-9 | 1.733 | -1.9e-13 | 1.733 | yes |
| 5 | -5.2e-5 | 2.262 | -5.5e-12 | 2.262 | yes |
| 8 | -0.15 | 2.712 | -1.4e-11 | 2.712 | yes |
| 10 | -0.91 | 2.770 | -3.1e-11 | 2.770 | yes |
| 12 | -6.3e7 | 3.635 | -6.9e-11 | 2.726 | yes |
| 15 | -3.8e10 | 3.920 | -6.6e-11 | 2.916 | yes |
| 20 | -1.4e22 | 4.883 | -1.0e-10 | 4.072 | yes |
| 25 | -6.1e30 | 5.500 | -1.1e-10 | 4.810 | yes |
| 30 | -5.5e36 | 5.907 | -7.7e-11 | 4.258 | yes |

* Only m = 1 is numerically non-negative on the whole half-line. For m >= 2 fhat of the float candidate is negative
  somewhere, down to about -1e-10 (m = 20, 25), and for m >= 20 the lowest point lies beyond xi = 4.
* In the weighted sense the violation is small (|fhat| <= 1.1e-10), but the polynomial itself is not: relative to
  (Tq)(0) = 1, P reaches -0.91 at m = 10 and -6e7 to -6e36 for m >= 12, at xi where exp(-pi xi^2) is about 1e-18 to
  1e-48. The float candidates are therefore far from non-negative as polynomials; repairing them is a B4 question.
* "reliable" means: polyroots converged, every accepted real root has relative residual <= 1e-20, and a dense grid
  (20001 points in xi) on [0, B], where B is a Fujiwara bound containing every stationary point of both problems,
  found no value below the stationary-point minimum. The grid check can only reveal a missed minimum that a grid
  point falls below, so this is strong numerical evidence, not a proof that the reported minimum is global.
  Where these criteria fail, `reliable` is false and no global minimum is claimed.

## Limitations

* float64 interior point: the d = 40 optimum (Gram eigenvalues ~1e-9, Fourier contacts) is beyond the
  accepted range here; the last digits need multiprecision (CGdL used SDPA-GMP at 3000 bits).
* The candidates are not feasible exactly: fhat is negative as above, and primal Gram matrices have eigenvalues down
  to -3e-11. Certification would need a strictly feasible perturbation and interval arithmetic (B4).
* The rho search is a bounded golden-section search on a possibly non-smooth, very flat function; rho* is not
  resolved for m >= 20.
* Witnesses are stored for m = 1, 10, 30 only; the other degrees keep their diagnostics and q coefficients.
