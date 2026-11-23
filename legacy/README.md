# Legacy experiments and historical prototypes

Two original research prototypes are preserved as **separate experiment snapshots**:

| Experiment | Snapshot | Description |
|---|---|---|
| v1 | [`prototypes/v1_positive_boxes/`](prototypes/v1_positive_boxes/ARCHIVE_STATUS.md) | Nonnegative centered-box mixtures; sensitivities under hypothetical band caps |
| v2 | [`prototypes/v2_signed_boxes/`](prototypes/v2_signed_boxes/ARCHIVE_STATUS.md) | Adds signed mixtures, sign-only tests, adaptive refinement and a remote sensitivity lobe |

**Version relationship:** v2 extends v1. They are two successive experiments, not two independent proofs. Every historical code file, test, log, plot and table from each source ZIP is archived in its corresponding snapshot. Only the v2 Japanese `signed_research_notes.md` has been translated to English; original and archived hashes are recorded in [`prototypes/SOURCE_MANIFEST.json`](prototypes/SOURCE_MANIFEST.json). Original Japanese source notes are retained in a separate, unmodified source backup outside the English-only repository.

The earlier frozen **pre-audit executable subset** remains at `pair_correlation_sensitivity_v0/`, exactly as it was, to keep the existing SHA-256 manifest and regression tests valid. It is not a third experiment: its five files are an audited executable subset of the v2 prototype. The corresponding manifest is `MANIFEST.sha256`.

**Known defects:** Both prototypes contain historical validation problems documented in `tests/test_legacy_known_defects.py`, the result registry, and audit reports. The corrected implementation is in `src/riemann_research/`. Do **not** present historical conditional numeric outputs as new proved zero-proportion bounds. Run original scripts only in scratch copies; they write into `results/`.
