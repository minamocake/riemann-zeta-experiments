"""Three independent status dimensions for every recorded result.

1. Computational status (`Computation`): how the number was obtained and checked.
2. Mathematical status (`Mathematical`): how far the mathematical argument has been checked.
3. Assumptions (`Assumption`, per hypothesis in results/assumptions.yaml).

Rules
-----
* The mathematical status is always set explicitly by a person in results/registry.yaml.
  No function in this package computes, derives or upgrades it, and nothing ever assigns PROVEN.
* `consistency_errors` only REJECTS combinations that overstate the evidence; it never fills a
  status in.
* Passing tests, or an exact certificate re-check, establishes a finite computation. It is not a
  proof of a statement about zeta zeros. See docs/derivations/what-a-certificate-establishes.md.
"""
from __future__ import annotations

from enum import Enum
from typing import Iterable, Mapping, Sequence


class Computation(str, Enum):
    EXPLORATORY = "exploratory"                      # one implementation, not independently reproduced
    REPRODUCED = "reproduced"                        # floating point, reproduced by an independent implementation
    NUMERICALLY_CERTIFIED = "numerically-certified"  # finite inequality re-checked in exact rational arithmetic
                                                     # by two independent verifiers (certify + independent_check)


class Mathematical(str, Enum):
    CONJECTURAL = "conjectural"                      # supported by computation only; no complete written argument
    PROOF_PENDING = "proof-pending"                  # complete written argument + certificate; not reviewed by a human
    INDEPENDENTLY_REVIEWED = "independently-reviewed"  # checked by a qualified human other than the author
    PROVEN = "proven"                                # accepted after human review (e.g. published); set by hand only


class Optimality(str, Enum):
    """Whether the value is optimal within its finite test-function family.
    Irrelevant for the validity of a bound (any admissible test function gives a valid bound)."""
    GLOBAL_CONVEX = "global-in-family-convex"
    GLOBAL_SDP_NUMERICAL = "global-in-family-sdp"
    LOCAL = "local"
    NOT_APPLICABLE = "not-applicable"


class Assumption(str, Enum):
    THEOREM = "theorem"
    THEOREM_UNDER_RH = "theorem-under-RH"
    RH = "riemann-hypothesis"
    DERIVED_UNREVIEWED = "derived-unreviewed"
    UNPROVEN_HYPOTHESIS = "unproven-hypothesis"
    OPEN_CONJECTURE = "open-conjecture"


HUMAN_REVIEW_REQUIRED = {Mathematical.INDEPENDENTLY_REVIEWED, Mathematical.PROVEN}


def consistency_errors(computation: Computation, mathematical: Mathematical,
                       reviews: Sequence[Mapping] = (), open_issues: Sequence[str] = (),
                       has_argument: bool = False, numeric_claim: bool = True) -> list[str]:
    """Return reasons why the recorded statuses overstate the evidence (empty list = acceptable).

    numeric_claim=False is for propositions (e.g. negative results) whose argument is a proof, not a
    certified number; for them the certified-computation requirement does not apply."""
    errors = []
    if mathematical in HUMAN_REVIEW_REQUIRED:
        human = [r for r in reviews if r.get("kind") == "human" and r.get("reviewer") and r.get("date")
                 and r.get("scope")]
        if not human:
            errors.append(f"'{mathematical.value}' requires a recorded human review (kind: human, reviewer, "
                          "date, scope); AI audits and passing tests do not count")
        if mathematical is Mathematical.PROVEN and not any(r.get("accepted") for r in human):
            errors.append("'proven' requires a human review that records 'accepted: true'")
    if mathematical is not Mathematical.CONJECTURAL:
        if numeric_claim and computation is not Computation.NUMERICALLY_CERTIFIED:
            errors.append(f"'{mathematical.value}' requires a numerically certified computation")
        if not has_argument:
            errors.append(f"'{mathematical.value}' requires a written argument (field 'argument')")
        if open_issues:
            errors.append(f"'{mathematical.value}' is not allowed while open issues remain: {list(open_issues)}")
    return errors


def conditional_on(assumption_statuses: Iterable[Assumption]) -> str:
    """Plain-language summary of what a claim is conditional on (display only)."""
    a = set(assumption_statuses)
    if a & {Assumption.UNPROVEN_HYPOTHESIS, Assumption.OPEN_CONJECTURE}:
        return "RH + UNPROVEN hypotheses"
    if Assumption.DERIVED_UNREVIEWED in a:
        return "RH + unreviewed lemma"
    return "RH"


def summary(computation: Computation, mathematical: Mathematical, assumption_statuses: Iterable[Assumption]) -> str:
    """Display string that repeats the recorded statuses verbatim; it never upgrades them."""
    return f"{mathematical.value} | {computation.value} | conditional on {conditional_on(assumption_statuses)}"


DISCLAIMER = ("Statuses describe computational checking, the state of the mathematical argument, and the "
              "assumptions used. Passing tests or certificate re-checks are not mathematical proofs.")
