"""The status logic may reject overstated statuses but must never assign or upgrade one."""
import ast
import inspect
import itertools

from riemann_research import status
from riemann_research.status import Assumption, Computation, Mathematical, consistency_errors, summary

HUMAN = {"kind": "human", "reviewer": "A. Reviewer", "date": "2027-01-01", "scope": "steps C1-C7"}


def test_no_code_path_in_the_package_produces_proven():
    """'proven' may appear only in the enum definition; nothing computes or returns it."""
    tree = ast.parse(inspect.getsource(status))
    uses = [n for n in ast.walk(tree) if isinstance(n, ast.Attribute) and n.attr == "PROVEN"]
    # allowed uses: the HUMAN_REVIEW_REQUIRED set and the 'accepted' check in consistency_errors
    assert len(uses) == 2
    for comp, assumptions in itertools.product(Computation, [[Assumption.RH], [Assumption.UNPROVEN_HYPOTHESIS]]):
        for m in Mathematical:
            assert summary(comp, m, assumptions).startswith(m.value)  # display repeats, never upgrades


def test_reviewed_and_proven_require_a_recorded_human_review():
    for m in (Mathematical.INDEPENDENTLY_REVIEWED, Mathematical.PROVEN):
        assert consistency_errors(Computation.NUMERICALLY_CERTIFIED, m, [], [], True)
        ai = [{"kind": "ai-audit", "reviewer": "Claude", "date": "2026-11-23", "scope": "all"}]
        assert consistency_errors(Computation.NUMERICALLY_CERTIFIED, m, ai, [], True)
    assert not consistency_errors(Computation.NUMERICALLY_CERTIFIED, Mathematical.INDEPENDENTLY_REVIEWED,
                                  [HUMAN], [], True)


def test_proven_additionally_requires_explicit_acceptance():
    assert consistency_errors(Computation.NUMERICALLY_CERTIFIED, Mathematical.PROVEN, [HUMAN], [], True)
    assert not consistency_errors(Computation.NUMERICALLY_CERTIFIED, Mathematical.PROVEN,
                                  [dict(HUMAN, accepted=True)], [], True)


def test_passing_tests_or_numerics_cannot_raise_a_numeric_claim_above_conjectural():
    for comp in (Computation.EXPLORATORY, Computation.REPRODUCED):
        for m in Mathematical:
            if m is not Mathematical.CONJECTURAL:
                assert consistency_errors(comp, m, [dict(HUMAN, accepted=True)], [], True)


def test_open_issues_or_missing_argument_block_proof_pending():
    certified, pending = Computation.NUMERICALLY_CERTIFIED, Mathematical.PROOF_PENDING
    assert consistency_errors(certified, pending, [], ["normalisation"], True)
    assert consistency_errors(Computation.NUMERICALLY_CERTIFIED, Mathematical.PROOF_PENDING, [], [], False)
    assert not consistency_errors(Computation.NUMERICALLY_CERTIFIED, Mathematical.PROOF_PENDING, [], [], True)


def test_conjectural_is_always_allowed():
    for comp in Computation:
        assert not consistency_errors(comp, Mathematical.CONJECTURAL)


def test_unproven_hypotheses_are_always_visible_in_the_summary():
    s = summary(Computation.NUMERICALLY_CERTIFIED, Mathematical.PROOF_PENDING,
                [Assumption.RH, Assumption.UNPROVEN_HYPOTHESIS])
    assert "UNPROVEN" in s
