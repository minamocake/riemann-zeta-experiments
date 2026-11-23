# Independent audit: pair-correlation band-sensitivity code (2026-11-23)

**Auditor:** Claude (Anthropic), acting as an independent third party at the request of the
repository owner. **Object:** the pre-audit code, preserved unchanged in
`legacy/pair_correlation_sensitivity_v0/` (SHA-256 manifest in `legacy/MANIFEST.sha256`).

The original audit was delivered in Japanese. This file is the English record, translated for
the repository's English-only policy. A short addendum at the end records what was learned when
the findings were turned into this repository.

## 0. What was actually run

* The package contained two source files and two run scripts. It had **no README and no tests**:
  `pytest -q` reported "no tests ran". No claim of passing tests could be checked.
* Both run scripts were executed on a copy (outputs archived in `legacy_run_outputs/`). The
  original files were not modified.
* Six independent check scripts were written and run (`scripts/`, outputs in `script_outputs/`).
  They use exact rational arithmetic, a physical-side integral as a second route, a convex QP
  solver, and an SDP relaxation. Only the import paths of these scripts were changed when they
  were copied into this repository.

## 1. What the code claims to compute

Under RH, for `h = Σ q_j (box of width L_j)`, `f = |ĥ|²`, `R = f̂ = h * h̃` and `Σ q = 1`:

    N*(T) ≤ (C + o(1)) N(T),   C = R(0) + 2 ∫₀¹ α R(α) dα + 2 ∫₁^A F(α) R(α) dα,

and the proportion of simple zeros is at least `P = 2 − C`. The code computes (a) the `A = 1`
baseline, (b) `P` under per-bin caps `U_i` on `∫_{B_i} F` (caps mode), (c) `P` with `R ≤ 0` for
`α ≥ 1` (sign-only mode), and (d) the sensitivities `∂C/∂U_i = 2 max(0, sup_{B_i} R)`.

## 2. Mathematically guaranteed (given RH and cited theorems)

* **Framework and implementation are correct.** `F ≥ 0` and Montgomery's asymptotic holds
  uniformly for `|α| ≤ 1` (Goldston-Montgomery), so the inequality is valid. Kernels and the
  integral `∫₀¹ x K dx` agree with exact rationals to 5e-14 and with the physical-side integral
  `∫ f μ̂` to 7 digits. The Fejér kernel reproduces 4/3.
* For `q ≥ 0`, `R` is nonincreasing, so a bin maximum is the value at the left edge (correct).
* `R` is piecewise linear, so knot evaluation finds the extrema (correct; a dense grid never beats it).
* The sensitivity formula follows from Danskin's theorem (for fixed `q` the functional is affine in
  `U`, so `C*` is concave in `U`).
* **RH-only values, re-checked in exact rationals:** A=1 baseline 0.6724490972; sign-only
  A=1.1: 0.6724406356, A=1.5: 0.6723649810, A=2.0: 0.6707507066. For A=1.5 and 2.0 the original
  float solutions had `R ≈ +1e-25` at some knots after rounding (not exactly admissible); setting
  the ~1e-12 weights on widths > 1 to zero makes them exactly admissible, with |ΔP| < 1e-13.
  All values are below the Montgomery-Taylor optimum (`N* ≤ 1.3275 N`, P ≈ 0.6725). The
  reference cited in the code, arXiv:1810.08843 (Chirre, Gonçalves, de Laat), already uses the same
  relaxation (`f̂ ≤ 0` for `|x| ≥ 1`) and obtains `N* < 1.3208 N`, i.e. more than 67.92% simple zeros.
  **Nothing here is new.**
* **Proposition found during the audit: sign-only degenerates** (`docs/derivations/sign-only-degeneracy.md`).
  For the largest width `w_k > 1`, on `((w_k + w_{k−1})/2, w_k)` only the self-term contributes,
  `R(α) = X_kk (w_k − α)/w_k² > 0`. So `R ≤ 0` on `[1, A]` forces zero weight on every box wider
  than 1; the mode reduces to the classical support-in-`[−1, 1]` problem and can never exceed 0.6725.
  The decrease of `P` with `A` comes from fewer boxes of width ≤ 1 (11, 8, 6).
* **A proven-cap variant the code does not consider.** Because `F` is the Fourier transform of a
  nonnegative measure, `∫F(α) R(α − b) dα ≤ ∫ F R` for any autocorrelation `R`; with
  `R = 2(1 − |α|)₊` this gives `∫ F ≤ 8/3 + o(1)` over any interval of length 1 under RH
  (`docs/derivations/int-F-cap.md`, unreviewed). Paying the positive part of `R` with this cap gives
  an exactly certified `P ≥ 0.6776349629` in the same 12-box family (a rank-3 PSD mixture): above
  Montgomery-Taylor, below 67.92%, not new.

## 3. Numerically reproduced

* **Nonnegative family:** all printed stages reproduced independently. The minimum eigenvalue on
  the tangent space in an orthonormal basis is 1.5e-2 to 3.1e-2, so the QP is convex and the
  reported values are global optima within the family (stronger than the code's own "local" claim).
* **Signed caps, 12 cases:** in 11 cases the SDP relaxation is numerically rank-one tight, so the
  legacy local optima are global even among PSD mixtures. At A=2, scale=8 a rank-2 mixture is better:
  0.68085 vs 0.68052 (40 seeds of the legacy optimiser did not improve on 0.68052).
* **Sensitivities:** agree with globally re-optimised finite differences to 2.8e-7.

## 4. Unproven assumptions

* RH, throughout.
* **Caps `U_i = scale × width` (all caps-mode values 0.70 to 0.93 depend on them).** The best known
  RH bounds for `∫_b^{b+ℓ} F` tend to 1, not 0, as `ℓ → 0` (Carneiro-Chandee-Chirre-Milinovich).
  `U ≤ 0.8` on a bin of width 0.1 is not known even for scale 8. Scale 1 is essentially the strong
  pair correlation conjecture on `[1, 2]`.
* **Refinement raises `P` because it strengthens the hypothesis** (the refined hypothesis implies the
  coarse one). For nonnegative boxes, A=2, scale 1, the limit of infinitely fine bins is
  "`F ≤ 1` a.e. on `[1, 2]`" with P = 0.932; the stages 0.908, 0.913, 0.914 approach it.

## 5. Errors and unverified items

1. **[Major, design]** The sign-only mode is degenerate (Section 2); the "signed kernels with larger
   A" extension cannot help in that mode.
2. **[Bug]** Partition coverage of `[1, A]` is not checked. With A=2.5 and default bins `[1, 2]`,
   `(2, 2.5]` is silently dropped: nonnegative case P=0.9707 vs 0.9356 with covering bins. Sign-only
   at A=2.5 returns `valid=True`, P=0.674637, although `R` reaches +1.5e-4 on `(2, 2.5]`.
3. **[Bug]** Partitions starting above 1 are accepted: A=1.5 with bins from 1.5 gives P=0.963
   (correct: 0.847).
4. **[Method]** Validity is a float tolerance (violation < 2e-7). A candidate with large weights on
   widths > 1 and P=0.673185 passes it. With the 8/3 cap the error is at most 1e-6, so it is
   quantitatively harmless, but the code's own logic does not prove it.
5. **[Minor]** `min_tangent_eigenvalue` uses a non-orthonormal basis (sign meaningful, magnitude not);
   width grids are not nested, so basis count 12 scores below 10; the "target ~0.6725" exceeds
   what the family reaches (0.672449); the restart filter (|q| ≤ 7) does not match the bounds (±8).
6. **Not checked:** the GitHub link in the code was not opened; Montgomery's theorem itself is taken
   from the literature; SDP values are floating point (only the explicit certificates are exact).

## Addendum (repository construction, 2026-11-23)

* An SDP solved with float tolerance for the sign-only mode at A=2 reported P=0.67291. It was first
  suspected to be spurious. The exact analysis showed that a violation of 2.6e-9 had cascaded into a
  weight of about 0.19 on a width > 1 (`check3b`): spurious as a sign-only claim, but such a candidate
  is legitimate once its positive part is paid for with the 8/3 cap. This led to the proven-cap
  variant above. `scripts/check3_exact.py` part (3) still ends with a solver error; the error is kept
  on purpose because it was the first sign of the degeneracy (margin constraints are infeasible).
* The audit-time certificate for the A=2 derived-cap result gave 0.6776349629. The repository
  regenerates it from a fresh SDP solve with different rounding, giving 0.6776349620
  (`certificates/pc-provencap-A2.json`); the registry records the committed value.
* Defects 2 to 4 are fixed in `src/riemann_research/` and covered by regression tests; the legacy
  tree is kept unchanged, and `tests/test_legacy_known_defects.py` documents its behaviour.

## Sources

* A. Chirre, F. Gonçalves, D. de Laat, *Pair correlation estimates for the zeros of the zeta
  function via semidefinite programming*, Adv. Math. 361 (2020) 106926, arXiv:1810.08843.
* E. Carneiro, V. Chandee, A. Chirre, M. B. Milinovich, *On Montgomery's pair correlation
  conjecture: a tale of three integrals*, arXiv:2108.09258.
* H. L. Montgomery (1973); D. A. Goldston and H. L. Montgomery (1987), as cited in the above.
