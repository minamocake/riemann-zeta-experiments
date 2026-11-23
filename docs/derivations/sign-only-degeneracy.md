# Proposition: the sign-only mode degenerates in the centered-box family

**Status:** `derived-unreviewed` (audit 2026-11-23). Checked numerically by
`tests/test_signonly_degeneracy.py`.

## Setting

Widths `0 < w_1 < ... < w_n = S`, a positive semidefinite Gram matrix `X` (rank one: `X = q qᵀ`),
and `R(α) = Σ_{i,j} X_ij K_ij(α)` with the box kernel `K_ij` of `src/riemann_research/kernels.py`.
The legacy "proven-from-positivity-only" mode requires `R(α) ≤ 0` for all `α ≥ 1`.

## Proposition

If `R(α) ≤ 0` for all `α ∈ [1, S]`, then `X_kk = 0` for every `k` with `w_k > 1`. Consequently the
corresponding rows and columns of `X` vanish, `R` is supported in `[-1, 1]`, and the problem is the
classical Montgomery problem restricted to boxes of width at most 1. Its value can therefore never
exceed the Montgomery-Taylor optimum (`N* ≤ 1.3275 N`, i.e. `P ≤ 0.6725`), whatever `S` is.

## Proof

Suppose not, and let `k` be the largest index with `w_k > 1` and `X_kk > 0`. For `j > k` we have
`X_jj = 0`, so row `j` of `X` vanishes (`|X_ij|² ≤ X_ii X_jj`). Hence `R = Σ_{i,j ≤ k} X_ij K_ij`.

`K_ij(α) = 0` once `α ≥ (w_i + w_j)/2`. Every pair `(i, j) ≠ (k, k)` with `i, j ≤ k` has
`(w_i + w_j)/2 ≤ (w_k + w_{k-1})/2 < w_k` (for `k = 1` there is no such pair). Choose
`α ∈ (max(1, (w_k + w_{k-1})/2), w_k)`, a nonempty interval contained in `[1, S]`. There

    R(α) = X_kk K_kk(α) = X_kk (w_k - α) / w_k²  >  0,

a contradiction. ∎

## Remarks

* **Why boxes fail here.** For even `h`, `R(S - u)` near the end of the support is the
  autoconvolution of the edge profile of `h`. For boxes the edge profile is constant, so `R > 0`
  just below `S`. Escaping this needs a different family (the 2020 pair-correlation benchmark of Chirre, Gonçalves and
  de Laat uses functions whose Fourier side is not compactly supported).
* **Fragility under tolerances.** If the constraint is only enforced up to `ε > 0`, the diagonal
  entries can grow in a cascade (roughly `ε^{1/2^m}` down the widths), so an `ε = 2·10⁻⁹`
  SDP solution carried weight about 0.19 on a width > 1 and reported a spurious gain of 0.0006.
  Float tolerances must not be read as feasibility. Such a candidate is still usable if its small
  positive part is paid for with a proven cap (see `int-F-cap.md`); this is what the `caps` mode
  with cap `8/3` does.
* **Legacy numbers.** The legacy sign-only values decrease with `A` (0.67244, 0.67236, 0.67075)
  because the number of boxes of width ≤ 1 drops (11, 8, 6), not because larger support hurts.
