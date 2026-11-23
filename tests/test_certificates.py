"""Exact certificate verification, including adversarial (tampered) certificates."""
import copy
import json
from fractions import Fraction as Fr

import numpy as np
import pytest

from riemann_research import certify


def load(root, name):
    return json.loads((root / "certificates" / name).read_text())


def all_certificates(root):
    return sorted((root / "certificates").glob("*.json"))


def test_every_committed_certificate_passes_the_exact_recheck(root):
    paths = all_certificates(root)
    assert len(paths) >= 7
    for p in paths:
        certify.verify_file(p)


def test_overclaimed_value_is_rejected(root):
    c = load(root, "pc-baseline-A1.json")
    c["claimed"]["P_lower_bound"] = "0.6725"
    with pytest.raises(certify.CertificateError, match="claimed"):
        certify.verify(c)


def test_support_beyond_window_is_rejected_in_window_mode(root):
    c = load(root, "pc-legacy-signonly-A2.0.json")
    c["outside"] = {"mode": "support-in-window"}
    with pytest.raises(certify.CertificateError, match="support exceeds"):
        certify.verify(c)


def test_tiny_positive_correlation_fails_sign_only_exactly():
    """A candidate the legacy 2e-7 tolerance would accept is rejected by the exact check."""
    W = [Fr(k, 6) for k in range(1, 13)]
    v = [Fr(0)] * 12
    v[5] = Fr(1)            # width 1
    v[11] = Fr(2, 10**7)    # width 2, tiny weight: R ~ +1e-7, inside the legacy 2e-7 tolerance
    cert = {"schema": certify.SCHEMA, "widths": [str(w) for w in W],
            "components": [{"weight": "1", "vector": [str(x) for x in v]}],
            "outside": {"mode": "sign-only"}, "claimed": {"P_lower_bound": "0"}}
    with pytest.raises(certify.CertificateError, match="R > 0"):
        certify.verify(cert)


def test_caps_certificate_with_partition_not_starting_at_one_is_rejected(root):
    c = load(root, "pc-provencap-A2.json")
    c["outside"]["edges"] = ["3/2", "2"]
    with pytest.raises(certify.CertificateError, match="cover"):
        certify.verify(c)


def test_caps_certificate_not_reaching_the_support_is_rejected(root):
    c = load(root, "pc-provencap-A2.json")
    c["outside"]["edges"] = ["1", "3/2"]
    with pytest.raises(certify.CertificateError, match="cover"):
        certify.verify(c)


def test_negative_weights_are_rejected(root):
    c = copy.deepcopy(load(root, "pc-provencap-A2.json"))
    c["components"][0]["weight"] = "-" + c["components"][0]["weight"]
    with pytest.raises(certify.CertificateError):
        certify.verify(c)


def test_scaling_the_gram_matrix_does_not_change_the_value(root):
    """Normalisation f(0) = 1 is done exactly by the verifier."""
    c = load(root, "pc-baseline-A1.json")
    base = certify.verify(c)["P"]
    c2 = copy.deepcopy(c)
    for comp in c2["components"]:
        comp["weight"] = str(Fr(comp["weight"]) * 7)
    assert certify.verify(c2)["P"] == base


def test_fejer_certificate_is_exactly_two_thirds():
    cert = certify.build("fejer", [Fr(1)], np.array([[1.0]]), "support-in-window")
    assert certify.verify(cert)["P"] == Fr(2, 3)
