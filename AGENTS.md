# AGENTS.md

Shared instructions for anyone, human or AI agent, changing this repository. Keep it short;
details live in the linked documents.

## Purpose

Computational research on Montgomery's pair-correlation method for zeta zeros. Every result
assumes RH. Reliability of claims matters more than speed or infrastructure.

## Commands (same as CI)

```bash
pip install -r requirements-lock.txt && pip install --no-deps -e .   # Python 3.13
ruff check .
python -m pytest
python tools/verify_certificates.py
python tools/check_registry.py --check     # use --write after editing results/registry.yaml
python tools/check_english_only.py
```

## Rules

1. **English only** in every file and name.
2. **Never upgrade a mathematical status.** Do not set `independently-reviewed` or `proven` in
   `results/registry.yaml`; only a person may, with a human review record. Passing tests or
   certificate checks are not proofs. See `docs/status-policy.md`.
3. **Every claim needs its assumptions listed** (`results/assumptions.yaml`). Unproven caps stay
   `unproven-hypothesis`.
4. **Certify before claiming.** Anything above `conjectural` needs a certificate accepted by both
   verifiers (`certify` and `independent_check`). Optimisers produce candidates only.
5. **Keep evidence.** Do not delete failed experiments, `legacy/`, `audits/` or
   `tests/data/historical_failures.json`. Record negative results in the registry.
6. **Do not edit `legacy/`** (SHA-256 manifest). Fix defects in `src/riemann_research/` with a regression test.
7. **Minimal, separate changes.** Translation or formatting, bug fixes and new results go in separate commits.
8. **No commits, pushes, merges, releases or GitHub setting changes without the owner's approval.**

## Where things are

* Status policy: `docs/status-policy.md`; what certificates establish: `docs/derivations/what-a-certificate-establishes.md`
* Research plan and open claims: `docs/research-plan.md`
* Results: `results/registry.yaml` (generated view: `docs/results-table.md`)
* Known limitations: `docs/known-limitations.md`; design decisions: `docs/decisions/`
