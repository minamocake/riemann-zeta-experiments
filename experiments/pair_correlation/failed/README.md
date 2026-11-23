# Failed and negative results (kept on purpose)

Each item is also a registry entry (`kind: negative-result`).

## 1. Signed kernels with larger support in the sign-only mode (`neg-signonly-degeneracy`)

**Hope:** allowing signed coefficients and support `A > 1` while requiring `R ≤ 0` on `[1, A]`
would beat the classical problem using only RH.
**Outcome:** impossible in the box family. Any weight on a box wider than 1 makes `R > 0` just
below its width, so all such weights vanish (`docs/derivations/sign-only-degeneracy.md`). The legacy
values 0.67244 / 0.67236 / 0.67075 (A = 1.1 / 1.5 / 2.0) are classical-problem values with 11 / 8 / 6
usable boxes.

## 2. Adaptive refinement of cap bins (`neg-refinement-strengthens-hypothesis`)

**Hope:** refining high-sensitivity bins improves the bound.
**Outcome:** with caps proportional to width, refinement strengthens the assumption; the improvement
(A=2: 0.90774, 0.91305, 0.91443) is not a methodological gain. Limit of infinitely fine bins: 0.93212.

## 3. Rank-one signed optimisation (`neg-rank-one-not-global`)

At A=2, scale 8, the best rank-one value over 40 seeds is 0.68052; a rank-2 PSD mixture gives
0.68085. Restricting to a single `H` loses value.

## 4. Float-tolerance SDP for sign-only (audit addendum)

An SDP with solver tolerances reported 0.67291 at A=2. The exact check found a 2.6e-9 violation that
had cascaded into weight 0.19 on a box wider than 1. Not admissible as a sign-only result; usable
only after paying its positive part with a proven cap (which led to `proven_cap/`).

## 5. No improvement over the literature (`neg-no-improvement-over-literature`)

Best RH-only value here: 0.67763 (with an unreviewed lemma). Pair-correlation benchmark: 0.6792.
