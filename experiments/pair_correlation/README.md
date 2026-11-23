# Pair-correlation experiments

All scripts run from a plain checkout (`python experiments/pair_correlation/<dir>/<script>.py`) and
use the shared core. The statuses of their results are in `docs/results-table.md`.

| Directory | Script | Produces | Extra assumptions |
|---|---|---|---|
| `baselines/` | `run_baseline.py` | `certificates/pc-baseline-A1.json` | none beyond RH |
| `signed_kernels/` | `certify_sign_only.py` | legacy sign-only certificates; conditional A=2 scale-4 certificate | none / `caps-scale-4` |
| `proven_cap/` | `run_proven_cap.py` (needs cvxpy) | `certificates/pc-provencap-*.json` | `lemma-intF-8-3` (unreviewed) |
| `sensitivity/` | `reproduce_caps.py` | `outputs/caps_reproduction.json` | `caps-scale-*` (unproven) |
| `failed/` | notes only | negative results kept for the record | |

Each certificate records its provenance (commit, generator, inputs, seed, precision, solver and
versions, verifier versions). Regenerating with the locked environment reproduces the committed
mathematical content exactly (`tools/compare_certificates.py`, used by the manual `research` workflow).

## Historical experiment snapshots

The original **v1 positive-box** and **v2 signed-box** runs are separate, complete snapshots under
[`legacy/prototypes/`](../../legacy/README.md). These are archival experiments, not new theorem
certificates; current experiments above use the corrected shared core.
