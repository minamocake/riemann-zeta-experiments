"""Falsification-first regression tests built from the failures found in the 2026-11-23 audit.

Every input here was ACCEPTED by the pre-audit code although it should not have been. The inputs
are stored in tests/data/historical_failures.json (regenerate with make_historical_failures.py).
Each test checks two things:
  1. the evidence still reproduces with the frozen legacy code, and
  2. the shared core and BOTH exact verifiers reject it.
"""
import json
from fractions import Fraction

import numpy as np
import pytest

import pair_sensitivity as legacy_pos
import signed_sensitivity as legacy_signed
from riemann_research import certify, independent_check
from riemann_research.montgomery import Partition, PartitionError, max_positive_part_outside
from riemann_research.optimize import positive_qp


@pytest.fixture(scope="module")
def evidence(root):
    return json.loads((root / "tests" / "data" / "historical_failures.json").read_text())


def cert_from(widths, components, mode, edges=None, caps=None):
    """Exact certificate from float data (binary floats taken exactly, no rounding)."""
    c = {"schema": certify.SCHEMA, "widths": [str(Fraction(float(w))) for w in widths],
         "components": [{"weight": str(Fraction(float(lam))), "vector": [str(Fraction(float(x))) for x in v]}
                        for lam, v in components],
         "outside": {"mode": mode}, "claimed": {"P_lower_bound": "0"}}
    if mode == "caps":
        c["outside"]["edges"] = [str(Fraction(float(e))) for e in edges]
        c["outside"]["caps"] = [str(Fraction(u)) if isinstance(u, str) else str(Fraction(float(u))) for u in caps]
    return c


def assert_both_verifiers_reject(cert, match_a=None):
    with pytest.raises(certify.CertificateError, match=match_a):
        certify.verify(cert)
    with pytest.raises(ValueError):
        independent_check.check(cert)


# 1. incomplete coverage of [1, A] -----------------------------------------------------------------
def test_incomplete_band_coverage_A_2_5(evidence):
    e = evidence["coverage_A2.5_positive"]
    m = legacy_pos.make_model(e["A"], 12)
    fit = legacy_pos.optimize(m, e["edges"], e["caps"])
    assert fit["P"] == pytest.approx(e["legacy_reported_P"], abs=1e-12)          # evidence reproduces
    with pytest.raises(PartitionError, match="supported up to"):
        Partition.build(e["edges"], e["caps"], support=e["A"])
    assert_both_verifiers_reject(cert_from(e["widths"], [(1, e["q"])], "caps", e["edges"], e["caps"]), "cover")
    edges = np.linspace(1, 2.5, 16)
    fit = positive_qp(e["widths"], Partition.build(edges, np.diff(edges), support=e["A"]))
    correct = fit["P"]
    assert correct == pytest.approx(e["correct_P_with_covering_bins"], abs=2e-6)
    assert correct < e["legacy_reported_P"] - 0.03
    # Independent oracle: the corrected value evaluated exactly by both verifiers.
    cert = cert_from(e["widths"], [(1, fit["q"])], "caps", edges, np.diff(edges))
    assert float(certify.verify(cert)["P"]) == pytest.approx(correct, abs=1e-9)
    assert certify.verify(cert)["P"] == independent_check.check(cert)["P"]
    # And the legacy optimum, evaluated honestly with covering bins, is worse than the legacy claim.
    honest = certify.verify(cert_from(e["widths"], [(1, e["q"])], "caps", edges, np.diff(edges)))["P"]
    assert float(honest) < e["legacy_reported_P"] - 0.03


# 2. band partition starting strictly above 1 -------------------------------------------------------
def test_partition_starting_above_one(evidence):
    e = evidence["partition_starting_at_1.5"]
    m = legacy_pos.make_model(e["A"], 12)
    assert legacy_pos.optimize(m, e["edges"], e["caps"])["P"] == pytest.approx(e["legacy_reported_P"], abs=1e-12)
    with pytest.raises(PartitionError, match="must start at 1"):
        Partition.build(e["edges"], e["caps"], support=e["A"])
    assert_both_verifiers_reject(cert_from(e["widths"], [(1, e["q"])], "caps", e["edges"], e["caps"]), "cover")


# 3. A > 2 in the sign-only ("theorem-safe") mode ---------------------------------------------------
def test_sign_only_with_support_beyond_two(evidence):
    e = evidence["sign_only_A2.5"]
    prog = legacy_signed.make_program(e["A"], 12)
    assert prog.validate_full_continuum(np.array(e["q"]), require_negative=True) < 2e-7  # legacy says valid
    assert max_positive_part_outside(e["widths"], e["q"]) > 1e-4                      # core checks all of [1, S]
    assert_both_verifiers_reject(cert_from(e["widths"], [(1, e["q"])], "sign-only"), "R > 0")


# 4. float tolerance accepting positive autocorrelation as nonpositive -----------------------------
def test_tolerance_accepted_positive_correlation(evidence):
    e = evidence["tolerance_accepted_positive_R"]
    prog = legacy_signed.make_program(2.0, 12)
    v = prog.validate_full_continuum(np.array(e["q"]), require_negative=True)
    assert 0 < v < e["legacy_tolerance"]                                                # legacy accepts
    sign_only = cert_from(e["widths"], [(1, e["q"])], "sign-only")
    assert_both_verifiers_reject(sign_only, "R > 0")
    # The same function IS usable once its positive part is paid with the derived cap (8/3 on [1, 2]);
    # the honest value is then below what the legacy code reported.
    capped = cert_from(e["widths"], [(1, e["q"])], "caps", [1, 2], ["8/3"])
    a, b = certify.verify(capped), independent_check.check(capped)
    assert a["P"] == b["P"]
    assert float(a["P"]) < e["legacy_reported_P"]


# 5. sign-only degeneracy: solver-tolerance SDP "improvement" ---------------------------------------
def test_sdp_cascade_is_not_a_sign_only_result(evidence, root):
    e = evidence["sdp_cascade_sign_only_A2"]
    d = json.loads((root / e["evidence_file"]).read_text())
    node = d
    for k in e["path"]:
        node = node[k]
    X, W = np.array(node["X"]), np.array(node["W"])
    assert node["P_sdp"] == pytest.approx(e["reported_P"], abs=1e-7)
    assert np.max(np.diag(X)[W > 1 + 1e-12]) > 0.1          # large weight on a box wider than 1
    lam, V = np.linalg.eigh(X)
    comps = [(lam[i], V[:, i]) for i in range(len(lam)) if lam[i] > 0]
    assert_both_verifiers_reject(cert_from(W, comps, "sign-only"), "R > 0")


def test_no_admissible_sign_only_function_puts_weight_on_wide_boxes():
    """Falsification attempt: search for a sign-only counterexample with weight on a box wider than 1."""
    rng = np.random.default_rng(7)
    W = [Fraction(k, 6) for k in range(1, 13)]
    for _ in range(40):
        v = [Fraction(int(x)) for x in rng.integers(-9, 10, size=12)]
        v[int(rng.integers(6, 12))] = Fraction(int(rng.integers(1, 9)))
        cert = {"schema": certify.SCHEMA, "widths": [str(w) for w in W],
                "components": [{"weight": "1", "vector": [str(x) for x in v]}],
                "outside": {"mode": "sign-only"}, "claimed": {"P_lower_bound": "-100"}}
        if sum(v) == 0:
            continue
        assert_both_verifiers_reject(cert, "R > 0")
