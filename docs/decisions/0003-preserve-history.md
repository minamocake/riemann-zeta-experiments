# 0003. Keep pre-audit code, failed experiments and defects

**Status:** accepted (v0.1)

**Decision.**
* The pre-audit code is kept byte for byte in `legacy/pair_correlation_sensitivity_v0/`, with a
  SHA-256 manifest checked by `tests/test_legacy_preservation.py`. Defects are fixed in the core,
  never in the legacy tree.
* Negative results (`kind: negative-result`) and known defects (`kind: known-defect`) are registry
  entries, not deleted notes. `experiments/pair_correlation/failed/` explains them.
* Recommended first history: commit 1 contains only `legacy/` (the pre-audit state); commit 2 adds
  the audit, core, tests and documentation. This makes the pre-audit state a real Git revision
  rather than only a copied directory.

**Consequences.** Readers can reproduce both the original claims and their corrections.
