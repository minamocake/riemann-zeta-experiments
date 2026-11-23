# 0002. Separate computation status from assumption status

**Status:** superseded by 0005 (kept for the record)

**Context.** A single "Verified / Numerical / Conjectural" field mixes two questions: was the
number checked exactly, and are the hypotheses true? The 71.239% figure shows why: its arithmetic
can be certified exactly, yet it rests on an unproven cap hypothesis.

**Decision.** `results/registry.yaml` records `computation` (and `optimality`) per result;
`results/assumptions.yaml` records the status of each assumption. The display label is derived
(`src/riemann_research/status.py`) and hand-written labels are rejected by the checker.

**Consequences.** Upgrading an assumption (for example, a lemma passing review) relabels every
dependent result automatically. See `docs/status-policy.md`.
