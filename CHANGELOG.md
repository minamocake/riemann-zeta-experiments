# Changelog

## 0.1.0 (unreleased)

### Fixes from the independent review of the changed files (2026-11-23)

* F2: exact coverage check in `Partition.build` (no 1e-12 tolerance); finite positive support required.
* F3: `independent_check.check` verifies the claimed bound on its own.
* F4: the manual research workflow fails on any audit-script error except the documented one.
* F5: values above the published record are allowed only with an explicit human audit flag.
* F6, F7: documentation of exact normalisation and of why F >= 0 holds.
* Reviewer's report and script archived under `audits/2026-11-23-independent-review-of-changed-files/`.

### Hardening (post-generation review)

* Status model: three independent dimensions (computation, mathematical status set by hand,
  assumptions). Tools never assign `proven`; reviewed/proven statuses need a recorded human review.
  Registry values unchanged; field values migrated (see docs/decisions/0005).
* Second exact verifier `independent_check` sharing no code with `certify`; both must agree.
* Machine-checked scopes for cap assumptions (8/3 lemma: bins of length <= 1, cap >= 8/3).
* Provenance in every certificate; certificates regenerated with identical mathematical content.
* Falsification-first regression tests on stored historical failing inputs (tests/data).
* Reproducibility and security: lock file, Python 3.13 pin, Ruff, read-only CI pinned to commit
  SHAs, manual research workflow, CODEOWNERS, documented repository settings.
* Behaviour-preserving lint fixes (`zip(..., strict=True)`, no lambda assignment).

### Initial foundation

* Shared core `riemann_research`: box kernels (float and exact), Montgomery functional,
  partition validation, exact certificate checker, two-axis status logic, optimisers.
* Fixed (relative to the legacy prototype): partitions that do not start at 1 or do not reach the
  support are rejected; validity is decided by exact rational certificates, not float tolerances.
* Result registry with separate computation and assumption status; generated results table.
* Certificates: classical baseline, legacy sign-only candidates (repaired), derived-cap results,
  and one conditional result to show that exact arithmetic does not remove hypotheses.
* Pre-audit prototype preserved byte for byte under `legacy/`; audit report, scripts and outputs.
* Negative results and known defects recorded in the registry.
