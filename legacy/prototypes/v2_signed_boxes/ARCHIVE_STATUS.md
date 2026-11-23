# v2 — Signed centered-box sensitivity prototype

**Historical extension of v1; not an independently established zeta-zero improvement.** v2 retains the positive-box experiment and adds signed coefficients, sign-only feasibility experiments, band sensitivities, adaptive subdivision, and numerical plots.

- Main added entry point: `run_signed_experiment.py`.
- Historical tests: `tests/test_pair_sensitivity.py` and `tests/test_signed_sensitivity.py`.
- All original code, historical tests, logs, CSV/JSON data, and PNG plots are byte-identical to the original v2 source ZIP.
- `signed_research_notes.md` is **translated into English** to comply with the repository language policy; its source and translated SHA-256 hashes are recorded in `../SOURCE_MANIFEST.json`. The original Japanese note is retained separately outside this English-only repository.
- Values assuming `U_i = λ |B_i|` are **conditional on unproven hypotheses**; they are not proven simple-zero proportions. The original sign-only feasibility tests also had defects identified by later independent review.
- Audit limitations and fixes: consult `../../../audits/` and corrected shared implementation in `../../../src/riemann_research/`.
- Do not run the original script directly inside the archive: it writes into its own `results/` folder; use a temporary copy.
