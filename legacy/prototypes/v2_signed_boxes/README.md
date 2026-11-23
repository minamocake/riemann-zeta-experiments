# Pair-correlation sensitivity experiment (research prototype)

## Status

This repo is a **conditional, counterfactual research prototype**, not a new
proved estimate for the zeta function. It **does not reproduce** the 67.92%
Chirre–Gonçalves–de Laat SDP bound. It **does** reproduce Montgomery's 4/3
bound with an exact continuous Fourier-support construction and reports the
sensitivity to additional *hypothetical* bandwise pair-correlation mass bounds.

**User's prior project:** [minamocake/adjoint_waveform](https://github.com/minamocake/adjoint_waveform), especially adjoint / analytical Jacobian / validation ideas. No existing GitHub files are modified; the HH neuronal ODE is not physically used here.

## Key definitions and proof-safe properties

Choose `h` as a nonnegative mixture of centered, normalized boxes with widths
up to `A`, and `∫h=1`. Let `f(y)=|∫h(t) exp(2π i yt) dt|²` and
`R_h(alpha)=∫h(t)h(t+alpha)dt`. Then:

* `f>=0` on the real line, `f(0)=1`.
* `R_h>=0` and `R_h(alpha)=0` outside `[-A,A]`.
* Each pairwise box correlation is piecewise linear, symmetric, and
  nonincreasing on positive `alpha`; so is the full mixture correlation.
* The upper bound for `R_h` on a whole band `[left,right]` equals **exactly**
  its value at `left`, independently of the sampling grid.

Assuming RH and, for every band `B_i` beyond frequency 1, a **proved** upper
bound `U_i >= limsup_T ∫_{B_i} F(alpha,T) d alpha`, the standard Montgomery
approach gives the multiplicity coefficient

    C <= R_h(0) + 2 ∫_0^1 alpha R_h(alpha) d alpha
         + 2 Σ_i U_i R_h(left_i).

The corresponding lower bound for simple zeros is `P >= 2 - C` (assuming RH).
The code treats the RHS as the quantity **C to minimize**.

**Crucial:** In experiments, `U_i = right_i-left_i`, corresponding to a
hypothetical Poisson `F=1` above frequency 1. These bounds are **not known**
for individual bins. Therefore the displayed conditional P numbers are not
established bounds and must not be reported as advances over the literature.
Refining a band into two strengthens the hypothetical assumptions.

## Numerics

The model uses a small positive nested-box basis (12 widths). Its quadratic
objective is minimized over the unit simplex with SLSQP and multiple starts.
Unless a positive-semidefinite Hessian on the tangent space is verified, the
reported minimum is only a candidate within the finite basis, not a certified
optimal bound. In this positive symmetric class the sensitivities always
decrease with frequency by construction; seeing the first bands selected is
**not a novel numerical discovery**. The next meaningful extension is to
allow signed `h` (while retaining `f=|FT(h)|²>=0`) and certify sign and band
suprema using exact piecewise breakpoints or interval arithmetic.

## Run

```sh
pip install -r requirements.txt
pytest -q
python run_experiment.py
```

Outputs are written to `results/`. See `research_notes.md`.

## Literature

Chirre, Gonçalves and de Laat, *Pair Correlation Estimates for the Zeros of the
Zeta Function via Semidefinite Programming*, Advances in Mathematics 361
(2020), 106926. https://arxiv.org/html/1810.08843 . Its RH-conditional
constant is `C<=1.3208` (thus `P>=0.6792`). The 2020 paper's class of
admissible functions is more expressive than the nested-box family here.

### Signed-function continuum validation (additional test)

`exact_band_extrema(model, q, left, right)` supports **signed** box weights and
computes exact continuous band extrema by enumerating every analytic breakpoint
of the piecewise linear cross-correlations. Its regression test constructs a
correlation which is negative at **both sampled band endpoints** yet positive
at an interior breakpoint. This is precisely the class of sign-checking error
that a simple sampling mesh can miss. The current optimizer remains restricted
to nonnegative box mixtures; signed optimization is a **next step**, not an
implemented scientific result.

## v2: Signed Fourier weights and adaptive window search — IMPLEMENTED

**Run:** `python run_signed_experiment.py && pytest -q`

The signed extension is implemented in `src/signed_sensitivity.py`, with
`tests/test_signed_sensitivity.py` and `signed_research_notes.md`.
Unlike the v1 positive-box optimizer, signed coefficients can produce
negative autocorrelation in parts of `1<=alpha<=A` while `f=|FT(h)|²>=0`
still holds. The maxima of signed autocorrelation on each frequency band are
evaluated by **enumerating every analytic breakpoint**, not relying on
mesh-based sign testing. SLSQP solves an epigraph formulation of the
conditional upper-bound objective. Local candidates only; no global or
machine-checkable number-theoretic certificate is asserted.

Two experiment tracks are separate:

* **Counterfactual upper mass caps:** `U_i = lambda * width(B_i)` for
  `lambda=1,2,4,8`. These caps are **not proved zeta-theoretic bounds**. Their
  signed coefficients produce tiny secondary outer-frequency sensitivity
  peaks for some settings; see `signed_research_notes.md`.
* **Known-information sign-only:** requires `R_h(alpha)<=0` for **all**
  `alpha>=1` (until compact support ends), so the unknown positive pair
  correlation can only reduce the multiplicity bound. Numerically feasible
  candidates in our restricted box basis did **not** improve the published
  2020 SDP result.

The adaptive scan supports `outer_quota=1` to force continuing to examine
one remote band beyond `alpha=1.5`, instead of greedily concentrating every
refinement near `alpha=1`.

The v1 introductory `Numerics` and old `Signed-function continuum validation`
sections above describe **v1**, not the expanded v2 model.
