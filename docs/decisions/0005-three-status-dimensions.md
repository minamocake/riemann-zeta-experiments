# 0005. Three status dimensions; mathematical status set only by people

**Status:** accepted (v0.1 hardening). Supersedes 0002.

**Context.** The first draft derived a display label ("Verified (conditional on RH)", ...) from
computation and assumptions. A derived "Verified" label can be read as a proof, and it conflated
checking a computation with checking an argument.

**Decision.**
* Record `computation` (exploratory / reproduced / numerically-certified), `mathematical_status`
  (conjectural / proof-pending / independently-reviewed / proven) and `assumptions` separately.
* `mathematical_status` is written by a person. Tools only reject overstatements; `proven` and
  `independently-reviewed` need a recorded human review, and `proven` also `accepted: true`.
* `numerically-certified` requires agreement of two verifiers that share no code.
* Cap assumptions carry machine-checked scopes.
* Values were not changed by the migration. Mapping: exact-certificate -> numerically-certified,
  numerical-reproduced -> reproduced, numerical-single -> exploratory.

**Consequences.** All certified results start as `proof-pending`; the rest are `conjectural`.
Nothing is reviewed or proven until a human review is recorded.
