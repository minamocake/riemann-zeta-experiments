# Where test expectations come from

A test is only as good as its oracle. Each expectation below is classified by its source.
"Characterization" values were recorded from a run of code and protect against regressions; they
are NOT independent evidence of correctness and are always paired with an independent check where
one is practical.

## Independent oracles (formula, specification, or a different computational route)

| Test file | Oracle |
|---|---|
| `test_kernels.py` | Closed form: Fejér kernel gives C = 4/3 exactly (Montgomery 1973). Brute-force convolution of `h` on a grid. Physical-side integral `∫ f μ̂` (Parseval; a different mathematical route). Dense trapezoid quadrature. Theorem: nonnegative mixtures give nonincreasing R. Piecewise linearity: a dense grid never beats knot evaluation. |
| `test_kernels.py::test_float_kernel_matches_exact_rational_kernel` | Partial only: float and exact kernels implement the same overlap formula; this catches floating-point bugs, not formula errors. The formula itself is checked by the brute-force and physical-side tests. |
| `test_independent_check.py` | Second verifier sharing no code (step functions, interval overlaps, Simpson). Fejér 2/3 by that route. A deliberately corrupted kernel is caught. |
| `test_certificates.py` | Specification of a valid certificate (modes, coverage, nonnegative weights, exact claims); exact scaling invariance; Fejér 2/3. |
| `test_signonly_degeneracy.py` | The degeneracy proposition (exact rational property test) and the literature value 0.6725 (Montgomery-Taylor). |
| `test_partition_validation.py` | Specification: bins must cover [1, S]. |
| `test_status.py`, `test_registry.py` | Specification in `docs/status-policy.md`; pair-correlation benchmark 0.6792 (literature). |
| `test_reproduction.py::test_sensitivity_formula...` | Danskin's theorem versus re-optimised finite differences. |
| `test_historical_failures.py` | Specification (each stored input must be rejected) and exact evaluation by both verifiers. |
| `test_english_only.py` | Unicode script ranges. |

## Characterization values (recorded from runs) and their independent pairing

| Value | Source | Independent pairing |
|---|---|---|
| `LEGACY_POSITIVE` stage values (`test_reproduction.py`) | legacy `run_experiment.py` output | different optimiser (convex QP, curvature checked) and exact evaluation of the same test function by both verifiers |
| 0.935599, 0.970667, 0.962963 | audit runs | exact evaluation by both verifiers of the corrected and legacy test functions (`test_historical_failures.py`) |
| Signed caps values in the registry | legacy rank-one optimiser | independent SDP solve (different algorithm and solver) |
| `legacy_run_outputs/` comparisons (`test_legacy_preservation.py`) | archived legacy outputs | none needed: these tests assert unchanged legacy behaviour by design |
| `tests/data/historical_failures.json` | legacy code (deterministic generator) | the tests check rejection by the core and both exact verifiers |

Optimality claims ("global in family") rest on convexity checks or SDP rank-one tightness and are
numerical; they do not affect the validity of any bound.
