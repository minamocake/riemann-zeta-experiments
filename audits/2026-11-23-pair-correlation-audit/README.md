# Audit 2026-11-23: pair-correlation band-sensitivity code

| Path | Content |
|---|---|
| `REPORT.md` | The audit report (English record of the original Japanese delivery, plus addendum) |
| `scripts/` | The six independent check scripts used in the audit. Only import paths were changed when copying them here. |
| `script_outputs/` | Their outputs when re-run from this repository on 2026-11-23 (`check3_exact` ends with a solver error by design; see the report addendum) |
| `legacy_run_outputs/` | Outputs of the legacy run scripts, produced on a copy during the audit |

Run a script from the repository root, for example
`python audits/2026-11-23-pair-correlation-audit/scripts/check4_edge_cases.py`.
`check2` to `check5` need the optional `sdp` extra (cvxpy).
