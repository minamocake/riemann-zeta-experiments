# Independent audit — `changed-files-only.zip`

Date: 2026-11-23. Scope: the 55 uploaded files only, **not** the complete repository. This is a software and finite-arithmetic audit, not a verification of a theorem about zeta zeros.

## Results reproduced

- ZIP entries: 55. The archive was extracted to a separate scratch directory; the source archive was not modified.
- Seven JSON certificates were recomputed using **both** `riemann_research.certify.verify` and `riemann_research.independent_check.check`. Their exact rational values of `C`, `P`, outside term, and outside maximum matched in all seven cases.
- Exact `P`: baseline A1 `0.6724490972435007`; hypothetical caps A2 scale4 `0.7123931138919648` (requires unproved cap hypotheses); derived-cap A2 `0.6776349620452484` (depends on **unreviewed** analytic 8/3 lemma).
- In the isolated scratch tree, after supplying minimal stand-in files for two unchanged missing pieces (`src/riemann_research/__init__.py` with `__version__`, and `tests/conftest.py` with a `root` fixture), the 26 tests in `test_independent_check.py`, `test_status.py`, and `test_kernels.py` all passed. The stand-ins are **not** evidence that the complete delivered repository works.
- Independent checker’s use of Simpson's rule is mathematically exact for this particular use: the implementation partitions at all step-endpoint-difference knots; on each interval `R` is linear, so `a*R(a)` is quadratic, and Simpson is exact over rationals.
- Mutation checks: a caps partition beginning at 3/2, a truncated partition, and a negative PSD weight were rejected by both verifiers. A false claimed bound was rejected by the main verifier but not by the independent verifier (see finding F3).
- The YAML files parse using PyYAML's BaseLoader and specify a PR trigger and a push-to-main trigger, read-only `contents: read`, plus immutable-SHA `uses:` entries.

## Findings

### F1 — Incomplete archive prevents a full go/no-go decision (scope blocker)

The submitted archive omits files used by the actual test suite and CI, including `src/riemann_research/__init__.py`, `tests/conftest.py`, `tests/test_certificates.py`, `tests/test_partition_validation.py`, `tests/test_signonly_degeneracy.py`, `tools/check_english_only.py`, `docs/derivations/int-F-cap.md`, `legacy/`, `audits/`, `LICENSE` and `CITATION.cff`. These omissions may be entirely intentional because the ZIP is named `changed-files-only.zip`; **do not infer that they are absent from the complete repository**. An untouched run against the extracted partial ZIP stops on import/collection errors, and the registry checker reports missing referenced files. Request the complete repository ZIP (excluding `.git`, virtual environments, caches and large generated files), then rerun under the locked Python environment.

### F2 — Candidate partition validator still accepts tiny uncovered gaps (medium)

`src/riemann_research/montgomery.py::Partition.build` accepts an edge starting at `1 + 5e-13` and a final edge ending at `S - 5e-13` because `EDGE_TOL=1e-12`. It also accepts `support=NaN`. These violate the literal mathematical specification of complete coverage. Importantly, **the exact certificate checker rejects such gaps**, so the seven checked certificates remain safe. For numerically produced candidates, reject nonfinite support and use a coverage rule that cannot omit any part of `[1,S]`, rather than a symmetric floating tolerance.

### F3 — Independent checker does not independently check the claimed bound (medium)

`independent_check.check(cert)` recomputes `P`, but never checks `cert['claimed']['P_lower_bound']`. When the claimed bound in `pc-provencap-A2.json` is changed to 1, `certify.verify` rejects it while `independent_check.check` accepts. The overall driver **does reject it** because the main verifier is called first; no false claim was accepted by the two-verifier driver. Still, the phrasing “both independent verifiers accepted the whole certificate” is too strong. Make the second checker validate the claim separately if that level of independence is desired.

### F4 — Manual research workflow hides most audit-script failures (medium)

`.github/workflows/research.yml` runs every `check*.py` with `|| echo "$s exited non-zero (check3_exact is expected to)"`. This turns **every** audit-script failure into a successful workflow step, not only the documented expected solver error. Use an explicit allowlist for the known expected failure and fail on any other nonzero exit.

### F5 — A test forbids future improvements over the published record (design issue)

`tests/test_registry.py::test_no_rh_only_result_exceeds_the_published_record` asserts every RH-only candidate remains below 0.6792. A genuine future improvement would make CI fail. Replace this as a hard invariant with a comparison that requires explicit human review or marks potential breakthroughs for audit. This historical record is a benchmark, not a mathematical upper bound.

### F6 — Documentation overstates raw input normalization (low)

`docs/derivations/what-a-certificate-establishes.md` says the certificate has `f(0)=1`. In fact, raw rational certificate components generally have `f(0) != 1` exactly; the verifiers **normalize by the exact rational total**. The computed `P` is therefore not affected, but the documentation should explicitly state “after exact normalization”.

### F7 — Small documentation inconsistencies (low)

`README.md` still refers to “two status axes” although the implementation has three dimensions. Update the sentence. Also avoid treating “Fourier transform of a nonnegative measure” **by itself** as proof of pointwise nonnegativity of `F`; pointwise nonnegativity needs the special positive-definite weight argument.

## Unresolved mathematical boundary

The analytic `8/3` cap lemma's proof (`docs/derivations/int-F-cap.md`) was not included in the ZIP. The seven verified arithmetic results do **not** establish the lemma or the implication to zeta-zero statistics. Treat `P=0.6776349620` as conditional on RH **and the unreviewed lemma** until the argument has been independently checked. `P=0.7123931138` additionally assumes currently unproved cap hypotheses.

## Checks not independently reproduced

- The full author-reported 84-test pass, the `legacy/` and `audits/` preservation checks, Ruff, actionlint, `check_english_only.py`, the manual research workflow, and GitHub-hosted CI.
- Reproduction on the exact claimed Python 3.13.18 and dependency versions in `requirements-lock.txt`. This run used the available Python 3.13.5 environment.
- The derivation of the `8/3` bound and the zero-theoretic steps C1–C7.

## Recommended next action

Obtain a **full** source snapshot, rerun the original CI commands on a clean checkout, fix F2/F4 and the misleading record-bound test, and independently review the analytic `8/3` lemma before describing any new unconditional-on-caps zeta bound. These are small targeted changes, not a reason to rebuild the repository.
