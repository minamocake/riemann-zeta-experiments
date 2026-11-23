# v1 — Nonnegative centered-box sensitivity prototype

**Historical experiment; not a verified new zeta-zero theorem.** The first prototype used nonnegative mixtures of centered box densities and hypothetical interval-mass caps beyond the classical Fourier window. It contains the original executable program, tests, research notes, data tables, and PNG plots.

- Original entry point: `run_experiment.py`.
- Historical test suite: `tests/test_pair_sensitivity.py`.
- Source and all 13 original files are preserved byte for byte.
- Outputs that assume `F = 1` beyond the known window are **counterfactual**, not unconditional results or advances over published bounds.
- Audit limitations and fixes: consult `../../../audits/` and the corrected shared implementation in `../../../src/riemann_research/`.
- Do not run the original script directly inside the archive: it writes into its own `results/` folder; use a temporary copy.
