# 0004. English-only repository

**Status:** accepted (v0.1)

**Decision.** All files, names, comments, documentation, test names, registry entries and CI
messages are in English. `tools/check_english_only.py` (run in CI) rejects Japanese, Chinese and
Korean script. Mathematical symbols are allowed.

**Translation record for v0.1.**
* The legacy code was scanned for non-English text before import: none was found (all comments,
  docstrings and identifiers were already English). No translation edits were made, so the
  legacy tree is byte-identical to the audited package (`legacy/MANIFEST.sha256`).
* Behavioural regression check: `tests/test_legacy_preservation.py` re-runs legacy functions and
  compares them with the outputs archived during the audit.
* The audit report, originally delivered in Japanese, was translated into English
  (`audits/2026-11-23-pair-correlation-audit/REPORT.md`).
* Bug fixes were made only in new code (`src/riemann_research/`), never mixed with translation.

**For future imports.** Translate comments and human-readable strings in a separate commit from
any fix, keep formulas, identifiers and APIs unchanged, review the diff, and run the regression
tests before and after.
