# Bin refinement with caps proportional to width strengthens the hypothesis

**Status:** `derived-unreviewed` (audit 2026-11-23). The argument is one line.

Let `𝒫` be a partition of `[1, 2]` and `H_s(𝒫)` the hypothesis `∫_B F ≤ s|B| + o(1)` for every
bin `B ∈ 𝒫`. If `𝒫'` refines `𝒫`, every `B ∈ 𝒫` is a disjoint union of bins `B'_i ∈ 𝒫'`, so

    ∫_B F = Σ_i ∫_{B'_i} F ≤ Σ_i s|B'_i| + o(1) = s|B| + o(1).

Hence `H_s(𝒫') ⇒ H_s(𝒫)`: the refined hypothesis is **stronger**. Under it the bound can only
improve, because the coarse bound remains available. The legacy adaptive scan therefore compares
results under different hypotheses; the rising `P` across stages (A=2, s=1: 0.90774, 0.91305,
0.91443) is not a methodological gain.

In the limit of arbitrarily fine bins the hypothesis becomes "`F ≤ s` in the weak sense on
`[1, 2]`", for which the nonnegative box family gives `P = 0.93212` at `A = 2`, `s = 1`
(registry id `pc-positive-A2-refinement-limit`).

The sensitivities `dC*/dU_i` are still correct local derivatives (Danskin's theorem); they measure
the value of tightening one cap, not the value of refining a bin.
