# Response to the independent review of `changed-files-only.zip` (2026-11-23)

`REPORT.md` and `audit_riemann_v01.py` in this folder are the reviewer's files, kept verbatim
(SHA-256 `1e6837e3...` and `e10a7cd8...`). Every finding was first reproduced on the repository,
then fixed with a minimal change and a regression test.

| Finding | Reproduced? | Change | Regression test |
|---|---|---|---|
| F1 Incomplete archive | Yes (by design of that ZIP) | The complete repository is delivered as one archive. The reviewer's script was run on it: 40 PASS, 5 WARN, 0 FAIL (before the `float` fix below). | n/a |
| F2 Partition tolerance accepts slivers and NaN support | Yes: `1+5e-13`, `S-5e-13` and `support=NaN` were accepted | `Partition.build` checks coverage exactly (`e_0 == 1`, `e_m >= S`) and requires finite positive support; the 1e-12 tolerance was removed. | `test_partition_validation.py::test_no_floating_tolerance_and_no_nonfinite_support` |
| F3 Second verifier does not check the claim | Yes: claim 1 accepted by `independent_check.check` | `independent_check.check` now rejects `P < claimed` on its own. | `test_independent_check.py::test_independent_verifier_checks_the_claim_itself` |
| F4 Research workflow masks audit-script failures | Yes (`|| echo`) | Only `check3_exact` may fail, and only with the documented `SolverError`; any other nonzero exit fails the job. | `test_registry.py::test_research_workflow_does_not_mask_audit_script_failures`, actionlint |
| F5 Test forbids exceeding the published record | Yes (design) | The hard invariant was replaced: a value above the record needs a human-written `record_comparison: {flagged_for_audit: true, note}` and is marked in the table. | `test_registry.py::test_values_above_the_published_record_need_an_explicit_audit_flag` |
| F6 Normalisation wording | Yes | `what-a-certificate-establishes.md` now states that raw components need not have `f(0) = 1` and both verifiers divide by `s` exactly. | documentation |
| F7 Wording ("two status axes"; nonnegativity) | Yes | README updated. `F ≥ 0` is now justified by the positive Fourier transform of the weight `w`, and the documents state that "Fourier transform of a nonnegative measure" alone does not imply it (montgomery.py, assumptions.yaml, int-F-cap.md). | documentation |

Reviewer-script warnings that remain (see `docs/known-limitations.md`): push runs on `main` only
(intended); no hash pinning (deferred to v0.2); a generic Simpson's-rule reminder (the reviewer's
own report confirms exactness for this use: R is linear between knots, so `a R(a)` is quadratic);
the 8/3 lemma is proof-pending. The `float` token warning was addressed by printing the exact
rational in `independent_check.py`.

Not changed: the reviewer's mathematical boundary statement is correct and already reflected in the
registry (`pc-provencap-A2` is conditional on RH and the unreviewed lemma; `pc-legacy-caps-A2-scale4`
also on unproven cap hypotheses).
