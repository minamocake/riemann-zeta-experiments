# riemann-zeta-experiments

Computational research on Riemann zeta zeros, pair correlation, Fourier optimization, and adjoint methods.

Repository: `minamocake/riemann-zeta-experiments`. The Python package keeps its name, `riemann_research`.

> **Scope and status.** This repository does not claim a proof of the Riemann Hypothesis or any
> new theorem. Every result assumes RH. Each result records its computational status, its
> mathematical status (set only by people; nothing is reviewed or proven yet) and its assumptions.
> Passing tests and exact arithmetic re-checks are not mathematical proofs.
> See [docs/status-policy.md](docs/status-policy.md).

## What is here (v0.1)

The first topic is Montgomery's pair-correlation method for bounding the multiplicity of zeta
zeros: for a test function `f = |H|² ≥ 0` with `f(0) = 1` and Fourier transform `R`,

    N*(T) ≤ (C + o(1)) N(T),   C = R(0) + 2∫₀¹ α R(α) dα + 2∫₁^∞ F(α) R(α) dα,

so at least `P = 2 − C` of the zeros are simple. The repository contains

* a shared core (`src/riemann_research/`) with box-function kernels in float and exact rational
  arithmetic, partition validation, an exact certificate checker, a second verifier that shares no
  code with it, provenance recording, and the status logic;
* two historical experiments (v1 nonnegative boxes and v2 signed boxes), preserved separately in
  [`legacy/prototypes/`](legacy/README.md), a frozen pre-audit executable subset, and independent
  audits (`audits/`);
* certificates, a result registry with three status dimensions, and negative results kept on purpose.

Contributors and AI agents: start with [AGENTS.md](AGENTS.md). The persistent research plan and
the list of claims awaiting review are in [docs/research-plan.md](docs/research-plan.md).

## Results at a glance

P is a lower bound for the proportion of simple zeros. The pair-correlation benchmark is **0.6792**
(Chirre, Gonçalves, de Laat, Adv. Math. 2020): the best published value we know of obtained with
Montgomery's pair-correlation method under RH. It is not the best known proportion of simple zeros
under RH by all methods. Larger numerical values below depend on unproven hypotheses.

| Mathematical status | Computation | Conditional on | P | Result |
|---|---|---|---|---|
| proof-pending | numerically-certified | RH | 0.6724490972 | Nonnegative boxes, Fourier support in [−1, 1] (classical problem) |
| proof-pending | numerically-certified | RH + unreviewed lemma | 0.677634962 | PSD mixture of 12 boxes; outside mass paid with a derived cap ∫F ≤ 8/3 |
| proof-pending | numerically-certified | RH + **UNPROVEN** hypotheses | 0.7123931138 | Signed boxes with caps ∫_B F ≤ 4·\|B\| on bins of [1, 2] |

"proof-pending" means a written argument and an exactly re-checked certificate exist, but no human
has reviewed the argument. The full, generated table is [docs/results-table.md](docs/results-table.md).
Values that depend on unproven hypotheses (such as the 71.24% figure above) must not be quoted as
proven proportions.

## Status dimensions

* **Computational status:** `exploratory`, `reproduced`, or `numerically-certified` (exact rational
  re-check by two independent verifiers).
* **Mathematical status:** `conjectural`, `proof-pending`, `independently-reviewed`, or `proven`.
  Written by people only; the last two require a recorded human review. AI audits do not count.
* **Assumptions:** ids in `results/assumptions.yaml` (RH, cited theorems, unreviewed lemmas,
  unproven hypotheses), each with its own status and, for caps, a machine-checked scope.

What a certificate does and does not establish is spelled out in
[docs/derivations/what-a-certificate-establishes.md](docs/derivations/what-a-certificate-establishes.md).

## Quick start

```bash
python3.13 -m venv .venv && source .venv/bin/activate
pip install -r requirements-lock.txt && pip install --no-deps -e .   # exact versions used for v0.1
ruff check .
python -m pytest                           # software checks (not proofs)
python tools/verify_certificates.py        # two independent exact verifiers
python tools/check_registry.py --check     # registry consistency, cap scopes, generated table
python tools/check_english_only.py
```

Expensive regeneration runs are in the manual `research` workflow, not in pull-request checks.

## Layout

```
src/riemann_research/   shared core: kernels, montgomery, certify, independent_check, optimize, provenance, status
experiments/            pair_correlation/{baselines, signed_kernels, proven_cap, sensitivity, failed}
certificates/           exact certificates (JSON, rational numbers as strings)
results/                registry.yaml (results), assumptions.yaml (status of assumptions)
legacy/                 v1/v2 historical snapshots, frozen audited subset, SHA-256 manifests
audits/                 audit reports, scripts and outputs
docs/                   status policy, derivations, design decisions, known limitations
tests/                  regression and consistency tests; tests/data holds historical failing inputs
tools/                  certificate, registry and English-only checkers
```

Planned areas (adjoint methods, duality certificates, GRH variants) will be added as new
directories under `experiments/` that reuse the core.

## History and negative results

The historical v1/v2 prototypes archived under `legacy/` had input-validation defects that produced invalid bounds; they are
fixed in the core, covered by regression tests, and documented as `known-defect` registry entries.
The audit also showed that the prototype's sign-only mode cannot exceed the classical problem in
the box family. Failed approaches are explained in
[experiments/pair_correlation/failed/](experiments/pair_correlation/failed/README.md).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). New results must enter through the registry; certificates
are welcome for anything claimed beyond `Numerical`.

## Citation and license

Cite via [CITATION.cff](CITATION.cff). Code is released under the MIT License ([LICENSE](LICENSE)).

## Acknowledgements

The 2026-11-23 audit and the v0.1 foundation were prepared with Claude (Anthropic) acting as an
independent reviewer; see [audits/](audits/2026-11-23-pair-correlation-audit/REPORT.md). A further
independent review of the hardening changes and the response to it are in
[audits/2026-11-23-independent-review-of-changed-files/](audits/2026-11-23-independent-review-of-changed-files/RESPONSE.md).
