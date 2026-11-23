# What a certificate establishes, and what it does not

**Mathematical status of this argument:** `proof-pending` (written 2026-11-23, not reviewed by a
human mathematician). Results that cite it can be no more than `proof-pending` until it is reviewed.

## A. What the two verifiers check (finite statements, exact rational arithmetic)

A certificate fixes rational widths `L_1, ..., L_n`, nonnegative rational weights `λ_k` and
rational vectors `v_k`. Both `certify.verify` and the independent `independent_check.check` confirm:

* **A1.** `X = Σ_k λ_k v_k v_kᵀ` is positive semidefinite (by construction: `λ_k ≥ 0`).
* **A2.** With `h_k = Σ_j v_{k,j} L_j⁻¹ 1_{[-L_j/2, L_j/2]}` and `s = Σ_k λ_k (∫h_k)²`, the raw
  components need not satisfy `f(0) = 1`; both verifiers check `s > 0` and divide by `s` exactly.
  After this exact normalisation, `f = s⁻¹ Σ_k λ_k |ĥ_k|²` has `f(0) = 1`, and its Fourier transform
  `R = f̂ = s⁻¹ Σ_k λ_k h_k * h̃_k` satisfies the mode condition:
  * `support-in-window`: `max L_j ≤ 1`, so `R = 0` for `|α| ≥ 1`;
  * `sign-only`: `R(α) ≤ 0` for all `α ∈ [1, S]`, `S = max L_j` (checked at every breakpoint of the
    piecewise-linear `R`, which is equivalent to checking the continuum);
  * `caps`: the bins `[e_{i-1}, e_i]` cover `[1, S]` and the outside term is
    `O = 2 Σ_i U_i max(0, max_{[e_{i-1}, e_i]} R)`.
* **A3.** `C = R(0) + 2∫₀¹ α R(α) dα + O` and `P = 2 − C` are computed exactly, and each verifier
  separately checks that `P` is at least the claimed value.

The two verifiers share no code: one uses the centered-box kernel formula and antiderivatives, the
other rebuilds `h` as a step function, sums interval overlaps, takes knots from endpoint
differences and integrates with Simpson's rule. They must agree exactly.

## B. What they do not check

* That RH holds, or that the cited theorems are correctly stated and applicable.
* That a cap `U_i` is true. The registry checker only verifies that each cap lies within the
  *scope* of a cited assumption (for example bin length ≤ 1 and cap ≥ 8/3 for `lemma-intF-8-3`).
* Optimality of the test function (irrelevant for validity).
* The deduction in Section C. It is written out but has not been formally verified.

## C. Deduction from a certificate to a statement about zeros (to be reviewed)

Conventions: `ĝ(ξ) = ∫ g(x) e^{-2πixξ} dx`; `F(α) = (T log T / 2π)⁻¹ Σ_{0<γ,γ'≤T} T^{iα(γ-γ')} w(γ-γ')`,
`w(u) = 4/(4+u²)`; zeros counted with multiplicity; `N(T) ~ (T/2π) log T`.

1. **Test function.** `f ≥ 0` and, after the exact normalisation of A2, `f(0) = 1`; `f ∈ L¹`, and `f̂ = R` is continuous with support
   in `[-S, S]`. (Each `|ĥ_k|²` is integrable because `ĥ_k(x) = O(1/|x|)`.)
2. **Explicit formula for pairs** (Montgomery): for such `f`,
   `Σ_{γ,γ'} f((γ-γ') log T / 2π) w(γ-γ') = (T log T / 2π) ∫ F(α) R(α) dα`.
3. **Diagonal.** Under RH, pairs with `γ = γ'` contribute `Σ_ρ m_ρ = N*(T)` (each zero of
   multiplicity `m` gives `m²` ordered pairs, summed over distinct zeros), each with weight
   `f(0) w(0) = 1`. Off-diagonal terms are `≥ 0` because `f ≥ 0` and `w > 0`. Hence
   `N*(T) ≤ (T log T / 2π) ∫ F R`.
4. **Window** (assumption `montgomery-theorem`): on `|α| ≤ 1`,
   `F(α) = (1+o(1)) T^{-2|α|} log T + |α| + o(1)` uniformly, so
   `∫_{-1}^{1} F R = R(0) + 2∫₀¹ α R(α) dα + o(1)` (R is Lipschitz; `∫ T^{-2|α|} log T dα = 1 − T^{-2}`).
5. **Outside the window** (assumption `F-nonnegative`, plus the mode):
   * `support-in-window`: the term vanishes;
   * `sign-only`: `F ≥ 0` and `R ≤ 0` on `[1, S]` give `2∫₁^S F R ≤ 0`;
   * `caps`: `2∫₁^S F R ≤ 2 Σ_i max(0, sup_{B_i} R) ∫_{B_i} F ≤ O + o(1)`, using the cap assumptions.
6. **Normalisation.** Steps 3 to 5 give `N*(T) ≤ (C + o(1)) (T/2π) log T = (C + o(1)) N(T)`.
7. **Simple zeros** (assumption `simple-zero-inequality`): `N_s(T) ≥ 2N(T) − N*(T) ≥ (P − o(1)) N(T)`.

## D. Checklist for a human reviewer (to move a result to `independently-reviewed`)

* Steps C2 and C4: exact statements and hypotheses in Montgomery (1973) and Goldston-Montgomery
  (1987), including uniformity up to `|α| = 1` and the admissible class of test functions.
* Step C3: multiplicity bookkeeping and the use of RH.
* Fourier and `F` normalisations match the code (`kernels.py`: `ĥ(x) = Σ v_j sinc(L_j x)`).
* For `caps` results: the cited cap lemma or hypothesis and its scope (e.g. `int-F-cap.md`).
* Record the review in `results/registry.yaml` under `reviews:` with `kind: human`, `reviewer`,
  `date`, `scope` (and `accepted: true` only for `proven`). Tools never set these fields.
