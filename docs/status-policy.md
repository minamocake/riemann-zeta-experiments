# Status policy

Every recorded result carries **three independent dimensions**. None of them is derived from the
others, and no tool ever upgrades one.

## 1. Computational status (`computation`)

| Value | Meaning |
|---|---|
| `exploratory` | One implementation, not independently reproduced. |
| `reproduced` | Floating point, reproduced by an independent implementation. |
| `numerically-certified` | A certificate in `certificates/` fixes the test function exactly, and two independent exact verifiers (`certify.verify`, `independent_check.check`) re-check the finite inequality in rational arithmetic and agree. |

A separate field, `optimality`, says whether the value is optimal within its finite family. It
does not affect validity: any admissible test function gives a valid bound.

## 2. Mathematical status (`mathematical_status`), set by hand

| Value | Meaning | Checker requirements |
|---|---|---|
| `conjectural` | Supported by computation only, or the argument is incomplete. | none |
| `proof-pending` | A complete written argument (`argument:`) and a numerically certified computation exist; no human review yet. | `numerically-certified` (for numeric claims), argument files exist, `open_issues` empty |
| `independently-reviewed` | A qualified human other than the author checked the argument. | a `reviews:` entry with `kind: human`, `reviewer`, `date`, `scope` |
| `proven` | Accepted after human review (for example published). | as above, plus `accepted: true` |

* Tools never assign or upgrade this field. `tools/check_registry.py` only rejects statuses that
  overstate the evidence; `tests/test_status.py` checks that no code path produces `proven`.
* AI audits (including the 2026-11-23 audit) and passing tests are evidence, not reviews.
* Changes to `results/` are merged manually by the maintainer after review of the diff; no automation
  merges them (`docs/repository-settings.md`).

## 3. Assumptions (`assumptions`)

Ids from `results/assumptions.yaml`, each with its own status: `theorem`, `theorem-under-RH`,
`riemann-hypothesis`, `derived-unreviewed`, `unproven-hypothesis`, `open-conjecture`. Cap
assumptions carry a `scope`; the checker rejects any certificate that uses a cap outside the scope
of the assumptions it cites (for example a bin longer than 1 with the 8/3 lemma).

A result can be `proof-pending` as a *conditional* statement and still depend on an unproven
hypothesis: `pc-legacy-caps-A2-scale4` (P = 0.71239) is "RH and caps-scale-4 imply P ≥ 0.71239".
The generated table therefore always shows what each row is conditional on, and that figure must
never be quoted as a proven proportion of simple zeros.

## What certificates and tests establish

See [derivations/what-a-certificate-establishes.md](derivations/what-a-certificate-establishes.md).
In short: a certificate establishes a finite arithmetic inequality; the deduction to a statement
about zeta zeros uses RH, cited theorems and the listed assumptions, and is not formally verified.
Passing tests show that software checks passed. Neither is a proof.

## Comparison with the literature

RH-only values are compared with the pair-correlation benchmark, currently 0.6792 (Chirre-Gonçalves-de
Laat 2020): the best published value we know of from Montgomery's pair-correlation method. It is not
the best known proportion of simple zeros under RH by all methods, and it is not an upper bound. A value above it (and not conditional on
unproven hypotheses) is accepted only with a human-written
`record_comparison: {flagged_for_audit: true, note: ...}`, and the generated table marks it
"EXCEEDS PAIR-CORRELATION BENCHMARK - flagged for audit". Such a value must be audited before it is reported.
