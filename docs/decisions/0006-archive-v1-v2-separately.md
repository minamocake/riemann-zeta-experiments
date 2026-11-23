# 0006. Preserve v1 and v2 as separate historical experiments

**Status:** accepted for v0.1 history completion (2026-11-23)

The original v1 positive-box prototype and the v2 signed-box extension were both supplied as
self-contained ZIPs. v2 extends v1 but includes additional analysis, tests, and outputs; preserving
only their shared executable code loses historical research context.

We therefore archive each ZIP's full contents in `legacy/prototypes/`, without rewriting its
programs, tests, logs, figures, or numeric data. The sole textual change to an original file is an
English translation of v2's Japanese `signed_research_notes.md`, with original and archived hashes
recorded in `legacy/prototypes/SOURCE_MANIFEST.json`. Unmodified original source archives remain
outside this English-only repository.

The original five-file `legacy/pair_correlation_sensitivity_v0/` subtree remains untouched because
existing regression tests and `legacy/MANIFEST.sha256` verify that exact copy. It is a frozen
subset of v2, not a separate line of research. Historic outputs are explicitly **not** certificates
of proved estimates; the audit-verified replacement is the shared core.

The CI pytest collection remains confined to root `tests/`; historical tests are preserved for
inspection and should run only in isolated scratch environments.
