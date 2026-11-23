"""The second verifier must be independent of the first and agree with it exactly."""
import ast
import json
from fractions import Fraction

import pytest

from riemann_research import certify, independent_check, kernels


def test_independent_verifier_imports_neither_kernels_nor_certify(root):
    tree = ast.parse((root / "src" / "riemann_research" / "independent_check.py").read_text())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            imported |= {node.module or ""} | {a.name for a in node.names}
    assert not {"kernels", "certify", "riemann_research.kernels", "riemann_research.certify",
                "montgomery", "optimize"} & imported
    assert imported <= {"__future__", "annotations", "fractions", "Fraction"}


def test_both_verifiers_agree_exactly_on_every_certificate(root):
    for p in sorted((root / "certificates").glob("*.json")):
        c = json.loads(p.read_text())
        a, b = certify.verify(c), independent_check.check(c)
        for k in ("C_window", "outside_term", "C", "P", "max_R_on_1_S"):
            assert a[k] == b[k], (p.name, k)


def test_fejer_closed_form_by_the_independent_route():
    c = {"widths": ["1"], "components": [{"weight": "3", "vector": ["2"]}], "outside": {"mode": "support-in-window"}}
    assert independent_check.check(c)["P"] == Fraction(2, 3)


def test_a_bug_in_the_shared_kernel_is_caught_by_the_independent_route(root, monkeypatch):
    c = json.loads((root / "certificates" / "pc-provencap-A2.json").read_text())
    honest = independent_check.check(c)["P"]
    original = kernels.kernel_exact
    monkeypatch.setattr(certify, "kernel_exact", lambda L, M, a: original(L, M, a) * Fraction(1001, 1000))
    try:
        broken = certify.verify(c)["P"]
    except certify.CertificateError:
        return  # also acceptable: the broken verifier rejects
    assert broken != honest
    assert independent_check.check(c)["P"] == honest


@pytest.mark.parametrize("mode", ["support-in-window", "sign-only"])
def test_independent_route_rejects_what_the_mode_forbids(mode):
    c = {"widths": ["1", "2"], "components": [{"weight": "1", "vector": ["0", "1"]}], "outside": {"mode": mode}}
    with pytest.raises(ValueError):
        independent_check.check(c)


def test_verification_report_records_results_and_provenance(root, tmp_path):
    import subprocess
    import sys
    out = tmp_path / "report.json"
    subprocess.run([sys.executable, str(root / "tools" / "verify_certificates.py"), "--report", str(out)],
                   check=True, capture_output=True, text=True)
    rep = json.loads(out.read_text())
    assert {"commit", "python", "libraries", "verifiers"} <= set(rep["provenance"])
    assert len(rep["results"]) >= 7 and all(r["result"] == "passed" for r in rep["results"])


def test_independent_verifier_checks_the_claim_itself(root):
    """2026-11-23 independent audit, finding F3: a false claim is rejected by BOTH verifiers on their own."""
    c = json.loads((root / "certificates" / "pc-provencap-A2.json").read_text())
    c["claimed"]["P_lower_bound"] = "1"
    with pytest.raises(certify.CertificateError, match="claimed"):
        certify.verify(c)
    with pytest.raises(ValueError, match="claimed"):
        independent_check.check(c)
