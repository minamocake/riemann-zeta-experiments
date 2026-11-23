# Lemma: a cap of 8/3 for the integral of F over any interval of length 1

**Status:** `derived-unreviewed`. Written during the 2026-11-23 audit. Elementary, but not peer
reviewed. Sharper published bounds exist (Carneiro, Chandee, Chirre, Milinovich, *On Montgomery's
pair correlation conjecture: a tale of three integrals*); this repository relies only on the
argument below.

## Statement

Assume RH. For every real `b`,

    ∫_{b-1/2}^{b+1/2} F(α, T) dα  ≤  8/3 + o(1)     (T → ∞).

## Proof

1. **F is the Fourier transform of a nonnegative measure.** By definition

       F(α) = (2π / (T log T)) Σ_{0<γ,γ'≤T} T^{iα(γ-γ')} w(γ-γ'),   w(u) = 4/(4+u²) > 0,

   so `F(α) = ∫ e^{2πiαx} dμ(x)` with `μ = (2π/(T log T)) Σ w(γ-γ') δ_{(γ-γ') log T / 2π} ≥ 0`.
   This is used only for the translation inequality in step 2. Pointwise `F ≥ 0`, used in step 4,
   is a separate fact (assumption `F-nonnegative`); it does not follow from `μ ≥ 0` alone.

2. **Translation does not increase the pairing with a positive-definite kernel.** Let
   `R₀(α) = 2 (1 - |α|)₊`. Its inverse Fourier transform is `f₀(x) = 2 (sin πx / πx)² ≥ 0`. Then

       ∫ F(α) R₀(α - b) dα = ∫ e^{2πibx} f₀(x) dμ(x)  ≤  ∫ f₀(x) dμ(x) = ∫ F(α) R₀(α) dα.

   The left side is real because F and R₀ are real; the inequality uses `f₀ ≥ 0` and `μ ≥ 0`.

3. **Evaluate the right side with Montgomery's theorem.** `R₀` is supported in `[-1, 1]`, where
   `F(α) = (1+o(1)) T^{-2|α|} log T + |α| + o(1)` uniformly (Montgomery; Goldston-Montgomery).
   The first term integrates against a continuous function `g` to `g(0) + o(1)`, so

       ∫ F R₀ = R₀(0) + ∫_{-1}^{1} |α| R₀(α) dα + o(1) = 2 + 2/3 + o(1) = 8/3 + o(1).

4. **Conclude.** `R₀(α - b) ≥ 1` for `|α - b| ≤ 1/2`, and `F ≥ 0`, `R₀ ≥ 0`. Hence
   `∫_{b-1/2}^{b+1/2} F ≤ ∫ F(α) R₀(α - b) dα ≤ 8/3 + o(1)`. ∎

## How it is used

For a test function whose Fourier transform `R` is supported in `[-S, S]` with `1 < S ≤ 2`,

    2 ∫_1^S F R  ≤  2 · max(0, sup_{[1,S]} R) · ∫_1^S F  ≤  (16/3) · max(0, sup_{[1,S]} R) + o(1).

This is mode `caps` with the single bin `[1, S]` and cap `8/3`, used by
`experiments/pair_correlation/proven_cap/`.
