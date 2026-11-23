"""The shared core reproduces the legacy numbers independently (numerical agreement, not proof)."""
import numpy as np
import pytest

from fractions import Fraction

from riemann_research import certify, independent_check
from riemann_research.montgomery import Partition, sensitivities
from riemann_research.optimize import positive_qp, sdp
from conftest import requires_cvxpy

LEGACY_POSITIVE = {  # from legacy run_experiment.py output (audit run 2026-11-23)
    1.1: [0.726363, 0.729280, 0.730734],
    1.5: [0.847045, 0.853599, 0.855412],
    2.0: [0.907741, 0.913049, 0.914430],
}


def exact_caps_certificate(widths, q, edges, caps):
    """Certificate taking the float inputs at their exact binary values (no rounding)."""
    fr = lambda x: str(Fraction(float(x)))  # noqa: E731
    return {"schema": certify.SCHEMA, "widths": [fr(w) for w in widths],
            "components": [{"weight": "1", "vector": [fr(x) for x in q]}],
            "outside": {"mode": "caps", "edges": [fr(e) for e in edges], "caps": [fr(u) for u in caps]},
            "claimed": {"P_lower_bound": "0"}}


def refine(edges, sel):
    out = [edges[0]]
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        if i in sel:
            out.append((a + b) / 2)
        out.append(b)
    return np.array(out)


@pytest.mark.parametrize("A", [1.1, 1.5, 2.0])
def test_positive_family_adaptive_scan_matches_legacy_and_is_globally_optimal(A):
    W = np.linspace(A / 12, A, 12)
    edges = np.linspace(1, 2, 11)
    for stage in range(3):
        part = Partition.build(edges, np.diff(edges), support=A)
        fit = positive_qp(W, part)
        assert fit["P"] == pytest.approx(LEGACY_POSITIVE[A][stage], abs=1e-6)   # characterization (legacy output)
        assert fit["optimality"].value == "global-in-family-convex"
        # Independent oracle for the value itself: exact evaluation of the same test function by both verifiers.
        cert = exact_caps_certificate(W, fit["q"], edges, np.diff(edges))
        a, b = certify.verify(cert), independent_check.check(cert)
        assert a["P"] == b["P"]
        assert float(a["P"]) == pytest.approx(fit["P"], abs=1e-9)
        s = sensitivities(W, fit["q"], part)
        edges = refine(edges, [int(i) for i in np.argsort(-s, kind="stable") if s[i] > 1e-10][:3])


def test_sensitivity_formula_matches_reoptimised_finite_differences():
    W = np.linspace(1.5 / 12, 1.5, 12)
    edges = np.linspace(1, 2, 11)
    U = np.diff(edges)
    fit = positive_qp(W, Partition.build(edges, U, support=1.5))
    s = sensitivities(W, fit["q"], Partition.build(edges, U, support=1.5))
    h = 1e-5
    for i in [0, 1, 4]:
        up, dn = U.copy(), U.copy()
        up[i] += h
        dn[i] -= h
        fd = (positive_qp(W, Partition.build(edges, up, support=1.5))["C"]
              - positive_qp(W, Partition.build(edges, dn, support=1.5))["C"]) / (2 * h)
        assert fd == pytest.approx(s[i], abs=1e-6)


@requires_cvxpy
def test_signed_caps_registry_values_are_reproduced_by_an_independent_sdp(root):
    import yaml
    reg = yaml.safe_load((root / "results" / "registry.yaml").read_text())
    vals = {e["id"]: e["value"] for e in reg["results"] if e["id"].startswith("pc-legacy-caps-A")}
    edges = np.linspace(1, 2, 11)
    for A in [1.1, 1.5, 2.0]:
        w = np.linspace(A / 12, A, 12)
        if not np.any(np.isclose(w, 1.0, atol=1e-12)):
            w[np.argmin(abs(w - 1))] = 1.0
        w = np.sort(w)
        for s in [1, 2, 4, 8]:
            key = f"pc-legacy-caps-A{A}-scale{s}" if not (A == 2.0 and s == 4) else "pc-legacy-caps-A2-scale4"
            P = 2 - sdp(w, "caps", Partition.build(edges, s * np.diff(edges), support=A))["objective"]
            if A == 2.0 and s == 8:
                assert P > vals[key] + 1e-4   # PSD mixture beats the legacy rank-one optimum here
            else:
                assert P == pytest.approx(vals[key], abs=3e-8)
