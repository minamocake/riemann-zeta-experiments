# Signed Fourier Window Sensitivity — Execution and Validation Report

**Date:** 2026-11-23  
**Scope:** Frequencies outside Montgomery's known pair-correlation window `|α| <= 1`, especially `α in [1, 2]`.  
**Research classification:** Numerical optimization prototype. **Not a new number-theoretic theorem or an improvement on a published bound.**

> **Historical translation:** This document is an English translation of `signed_research_notes.md` from the original v2 archive. It records statements made *before* independent auditing; some original procedures had later-identified validation defects. It must not be cited as current certification. Consult the main repository's audit reports and corrected core. Original source-file SHA-256 values are recorded in `legacy/prototypes/SOURCE_MANIFEST.json`.

## Research question

If Fourier weights outside the known window are allowed to change sign, which frequency bands' unknown upper bounds most strongly influence the optimal lower bound on the proportion of simple zeros? In particular, can a remote second peak emerge in the sensitivity map, which is impossible for the monotone sensitivity profile of a mixture with only nonnegative coefficients?

## Precise mathematical distinctions

Let `h` be a real-valued, compactly supported function normalized by `∫h = 1`, expressed as a **signed** combination of centered box densities of widths `L_j <= A`: `h = Σ q_j b_j`, `Σ q_j = 1`. Then `f(y) = |FT h(y)|² >= 0` and `f(0) = 1` hold exactly. Its autocorrelation `R_h(α) = FT f(α) = ∫h(t)h(t+α)dt` is even, compactly supported, and piecewise linear, and **it may be negative**.

Assuming RH, the known-window contribution is

    C_known(q) = R_h(0) + 2∫_0^1 α R_h(α) dα.

If an independent upper bound `limsup_T ∫_{B_i} F(α,T)dα <= U_i` can be **proved** for each interval `B_i` outside the known window, using `F >= 0` gives

    C(q,U) = C_known(q) + 2 Σ U_i max(0, sup_{α∈B_i} R_h(α)).

The proposed argument uses this coefficient as an upper bound on `N*(T)/N(T)`, hence a simple-zero proportion `P >= 2-C`. Every such conclusion is conditional on RH and on the stated bounds `U_i` actually holding.

**The experimental values `U_i = λ |B_i|`, for `λ ∈ {1,2,4,8}`, were hypothetical bandwise bounds, not bounds proved at the time.** Therefore the tabulated `P` values are conditional optimization outputs, not established proportions of simple zeros.

A **separate experiment** imposed `R_h(α) <= 0` throughout `1 <= α <= A` without assuming additional bounds `U_i`. This can use only Montgomery's known window and `F >= 0`, but the floating-point optimizer's candidates did not yet have exact analytical, rational, or interval certificates. They were therefore numerical candidates, not established results.

## Coarse scan: effect of allowing signs

For each `A`, the range `[1,2]` was initially divided into ten bins of width 0.1. The experiment used 12 basis functions, local nonlinear optimization, and multiple initial points.

| A | λ | Nonnegative weights P | Signed weights P | Gain (percentage points) | Sign changes? |
|---|---:|---:|---:|---:|---|
| 1.1 | 4 | 69.956% | 69.956% | 0.000 | No |
| 1.5 | 4 | 70.797% | **70.935%** | **+0.139** | Yes |
| 2.0 | 2 | 77.371% | **77.436%** | **+0.065** | Yes |
| 2.0 | 4 | 70.682% | **71.239%** | **+0.557** | Yes |
| 2.0 | 8 | 67.515% | **68.052%** | **+0.537** | Yes |

In particular, for `A=2`, `λ=4`, the smallest coefficient was approximately `-0.0818`; `R_h(α)` became negative near `1.2–1.6` and formed a **very small positive lobe** again near `1.7–1.9`. The positive centered-box mixture in the previous experiment could not show this behavior. Thus signed degrees of freedom affected the numerical optimization result.

**Audit qualification:** These numbers and the claimed feasible solutions are historical experimental outputs under hypothetical assumptions. The subsequent audit found input-validation and sign-checking issues in the original programs; consult the corrected core and certified registry rather than treating this paragraph as a valid bound.

### Local sensitivity

If the local derivative exists and the optimum changes smoothly, the envelope theorem gives

    S_i = ∂C*/∂U_i = 2 max(0, sup_{B_i} R_{h*}(α)).

The original `A=2`, `λ=4` example reported:

| Band | S_i |
|---|---:|
| 1.00–1.10 | 0.122120 |
| 1.10–1.20 | 0.043360 |
| 1.20–1.30 | 0.000000 |
| 1.70–1.80 | 0.000634 |
| 1.80–1.90 | 0.000634 |
| 1.90–2.00 | 0.000350 |

Sensitivity remained strongest near the frequency-window boundary. The second peak was only **approximately 0.52%** of the largest sensitivity. Identifying a region worth exploring is not the same as proving that it is the best target for future number-theoretic estimates.

## Adaptive refinement

For `A=2`, `λ=4`, starting from bins of width 0.1, the three highest-ranking intervals were refined. The initial selection was `1.00–1.10`, `1.10–1.20`, and `1.70–1.80`. Subsequently, selecting only by highest sensitivity caused the search to concentrate near the boundary.

A second strategy therefore maintained **two top-ranked bins plus one remote bin (`α >= 1.5`) in every round**. This kept `1.70–1.75` among the next targets and allowed comparison of the refined sensitivity structure.

Caution: imposing the original hypothesis `U_i = λ |B_i|` independently on the sub-bins after subdivision **strengthens the assumptions**. An increase in `P` after refinement cannot be attributed solely to improved numerical precision.

## Numerical stability and tests

- The original `pytest -q` run reported **23 passed**.
- For the signed candidate `A=2`, `λ=4`, the original computation reported satisfaction of interval constraints at floating-point accuracy; the maximum continuum violation was on the order of `1e-13` or smaller.
- The original sensitivity computation was compared with central finite differences **after reoptimization** for perturbed upper bounds.
- Since `R_h` is analytically piecewise linear, enumerating all breakpoints of box cross-correlations yields its interval maxima without depending only on a sampling mesh. A regression test also reproduced the earlier failure mode where the endpoints were negative but an interior point was positive.
- Sign changes and a second sensitivity peak appeared with 8, 10, 12, 16, and 20 basis elements. **The values were basis-dependent**, and small variations among local minima occurred in the multistart computations with 20 basis functions.
- SLSQP is a local optimizer; **there was no guarantee of global optimality across the entire admissible function space**. Small numerical gains and comparisons across bases have no automatic force as mathematical proofs.

**Audit qualification:** The statements above describe the original prototype's testing, not the corrected implementation's certification or a current audit conclusion.

## Results using only proved input information

In the original numerical run with no extra integral caps and `R_h <= 0` over the whole unknown window, the 12-element basis gave the following candidates:

| A | Numerical P candidate |
|---|---:|
| 1.1 | 67.244% |
| 1.5 | 67.229% |
| 2.0 | 67.075% |

**This function family did not reproduce the 2020 RH-conditional 67.92% benchmark.** Neither this mode nor its numerical candidates established a new number-theoretic bound. The family of simple symmetric box functions might lack adequate expressiveness; subsequent independent audit established stronger structural restrictions within the centered-box family.

## Questions for future work

1. Can the **independent bandwise upper-bound hypotheses** be replaced with proved cumulative bounds or dual SDP constraints?
2. If the basis is extended to shifted box functions, B-splines, or Gaussian-polynomial families, will the second sensitivity peak survive?
3. Can exact certificates be built, using rational coefficients or interval arithmetic to check `R <= 0` over an entire interval?
4. Can band sensitivity and dual SDP solutions be used to reverse-engineer the **minimal additional inequality** that would have to be proved to reach 70.37%?

## Reproduction of the original prototype

```bash
cd pair_correlation_sensitivity
pip install -r requirements.txt
python run_signed_experiment.py
pytest -q
```

CSV, JSON, and PNG outputs are in `results/`. `run_signed_experiment.py` does not authenticate to or upload to GitHub. The original `adjoint_waveform` repository was not modified. The **research workflow** was informed by its Hodgkin–Huxley sensitivity and numerical-derivative validation approach; the HH differential equations were not imported directly.

**Prior literature:** Chirre, Gonçalves, and de Laat, *Pair Correlation Estimates for the Zeros of the Zeta Function via Semidefinite Programming*, *Advances in Mathematics* 361 (2020), 106926. https://arxiv.org/html/1810.08843
