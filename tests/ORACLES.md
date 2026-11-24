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
| `test_cgdl_fourier.py` | Two implementations of the Gaussian-polynomial transform (Laguerre closed form of CGdL Section 4.1 versus the differentiation identity FT[x^2 h] = -(hhat)''/(4 pi^2)), compared exactly in Q[pi, 1/pi]. Fourier inversion T(T p) = p (exact). Hand-derived closed forms. Direct Fourier quadrature at 40 digits. Fejér kernel and Z(H) = 4/3 (Montgomery 1973). Z by incomplete-gamma moments versus quadrature in the dilated variable. Exact square-free decomposition and Sturm sequences for sign radii. Gaussian moments k = 0..80 (`mpmath.gammainc`) versus quadrature, relative error, including values near 1e-238; a positive-term series for gamma(81, pi). Rigorous pi enclosure from Machin's formula in exact rationals. Mutation checks during development: corrupting one Laguerre coefficient makes two tests fail; restoring the finite-sum incomplete gamma makes five tests fail. |
| `test_cgdl_sdp.py` | Analytic m = 0 point (q = 1 - t, X3 = 1 - 1/(2 pi), X4 = 1) and the B1a Gaussian Z. A numpy re-implementation of the linear maps (`evaluate_point`) versus the CVXPY model. Embedding of the m = 0 point into m = 1. The map q_B = a0 q_A (objective and constraints preserved). V_B <= V_A (theorem: A is B plus one constraint). Z re-evaluated by both B1a routes using the actual sign radius of the returned q. Dual homogeneity identity -(sum of normalisation duals) = V. KKT stationarity residual and primal-dual gap from a numpy model of the SDP built by probing the affine maps (independent of CVXPY's canonicalisation); a corrupted dual is detected. Explicit q(0) > 0 check. Continuous Fourier minimum from stationary points (analytic cases: (u^2 - 1)^2 + eps, the Gaussian example). A constructed infeasible instance. |
| `test_cgdl_laguerre.py` | Exact: B1a's T is diagonal (+-1) on L_n^(-1/2)(2 pi t); the three-term recurrence; basis round trips and exact products. Objective integrals versus quadrature. B2-mini (monomial basis) values at m = 0, 1. Agreement of two Gram bases. Nested families (V non-increasing in m). Independent re-evaluation from the Gram matrices with B1a's monomial transform at 200 digits. Fourier-side minima: closed-form cases for the unweighted (P' = 0) and weighted (P' - pi P = 0) problems, boundary t = 0 and both tail behaviours, and a sabotaged root finder that must be reported as unreliable. Acceptance with missing or non-finite diagnostics; a failing settings rung must not relabel an earlier solution. Solver-free replay of the saved witnesses with cvxpy/clarabel imports blocked; mutation checks: corrupting a primal entry, symmetry, an equality dual, a dual PSD block, a recorded diagnostic, rho, the acceptance flag, shapes, the witness file, the sweep's saved q coefficients, a status, the acceptance tolerance or the reported Fourier infimum makes the replay fail. The published-witness checker rejects a file with another SHA-256. Optional: CGdL's published d = 40 witness (set CGDL_Z40). |

## Characterization values (recorded from runs) and their independent pairing

| Value | Source | Independent pairing |
|---|---|---|
| `LEGACY_POSITIVE` stage values (`test_reproduction.py`) | legacy `run_experiment.py` output | different optimiser (convex QP, curvature checked) and exact evaluation of the same test function by both verifiers |
| 0.935599, 0.970667, 0.962963 | audit runs | exact evaluation by both verifiers of the corrected and legacy test functions (`test_historical_failures.py`) |
| Signed caps values in the registry | legacy rank-one optimiser | independent SDP solve (different algorithm and solver) |
| `legacy_run_outputs/` comparisons (`test_legacy_preservation.py`) | archived legacy outputs | none needed: these tests assert unchanged legacy behaviour by design |
| `tests/data/historical_failures.json` | legacy code (deterministic generator) | the tests check rejection by the core and both exact verifiers |
| Z = 1.4525472652249223, (T p)(2) = -3.2077135983, fhat(2) = -1.11863955e-5 (`test_cgdl_fourier.py`) | separate scipy and mpmath computations made while writing the B1a specification | closed form versus quadrature inside the same tests; exact symbolic transforms |
| V_A(m=1) = 1.4525472652591453, V_B(m=1) = 1.328646166101496 (`test_cgdl_sdp.py`) | Clarabel run of B2-mini | independent Z by closed form and quadrature; numpy re-evaluation; V_B <= V_A and V_A(m=1) <= Z of the embedded m = 0 point |
| Fourier near-contact 4.092813391e-9 at u = 1.0668272166 (`test_cgdl_sdp.py`) | stationary points of the saved B, m = 1 coefficients; matches the independent reviewer's value | analytic test cases of `even_poly_minimum`; consistent with (larger) grid minimum |
| S1 sweep values in `experiments/cgdl_reproduction/b2_full/results/s1_sweep.json` (`test_cgdl_laguerre.py`) | Clarabel runs | acceptance re-applied to every saved candidate; solver-free replay of the m = 1, 10, 30 witnesses; monotonicity in m (consistency of the rho search); high-precision re-evaluation; published-witness evaluation of the same model |

Optimality claims ("global in family") rest on convexity checks or SDP rank-one tightness and are
numerical; they do not affect the validity of any bound.
