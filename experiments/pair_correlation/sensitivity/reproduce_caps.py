"""Independent reproduction of the legacy counterfactual-cap results with the shared core.

The caps U_i = scale * (bin width) on [1, 2] are an UNPROVEN hypothesis (see results/assumptions.yaml).
Refining bins with U_i proportional to width makes the hypothesis STRONGER, so P increasing with
refinement is not a methodological gain. Writes outputs/caps_reproduction.json.
Signed (PSD) results need cvxpy; they are skipped if it is missing.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_common"))
import paths  # noqa: E402

from riemann_research import provenance  # noqa: E402
from riemann_research.montgomery import Partition, sensitivities  # noqa: E402
from riemann_research.optimize import positive_qp, sdp  # noqa: E402


def refine(edges, selected):
    out = [edges[0]]
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        if i in selected:
            out.append((a + b) / 2)
        out.append(b)
    return np.array(out)


def legacy_signed_widths(A, n=12):
    w = np.linspace(A / n, A, n)
    if not np.any(np.isclose(w, 1.0, rtol=0, atol=1e-12)):
        w[np.argmin(abs(w - 1.0))] = 1.0
    return np.sort(w)


def positive_adaptive(A, rounds=2, top_k=3):
    widths = np.linspace(A / 12, A, 12)
    edges = np.linspace(1.0, 2.0, 11)
    stages = []
    for stage in range(rounds + 1):
        part = Partition.build(edges, np.diff(edges), support=A)
        fit = positive_qp(widths, part)
        sens = sensitivities(widths, fit["q"], part)
        stages.append({"stage": stage, "bins": int(edges.size - 1), "P": fit["P"],
                       "optimality": fit["optimality"].value,
                       "tangent_min_eigenvalue": fit["tangent_min_eigenvalue"]})
        ranked = [int(i) for i in np.argsort(-sens, kind="stable") if sens[i] > 1e-10]
        edges = refine(edges, ranked[:top_k])
    return stages


def main():
    out = {"hypothesis": "int_{B_i} F <= scale * |B_i| + o(1) on bins of [1,2] (UNPROVEN)",
           "positive_adaptive_scale1": {}, "signed_psd": {}}
    for A in [1.1, 1.5, 2.0]:
        out["positive_adaptive_scale1"][str(A)] = positive_adaptive(A)
    try:
        import cvxpy  # noqa: F401
        edges = np.linspace(1.0, 2.0, 11)
        for A in [1.1, 1.5, 2.0]:
            for s in [1, 2, 4, 8]:
                w = legacy_signed_widths(A)
                res = sdp(w, "caps", Partition.build(edges, s * np.diff(edges), support=A))
                out["signed_psd"][f"A{A}_scale{s}"] = {"P": 2 - res["objective"],
                                                       "eigenvalue_ratio": res["eigenvalue_ratio"]}
    except ImportError:
        out["signed_psd"] = "skipped: cvxpy not installed"
    out["provenance"] = provenance.collect(
        "experiments/pair_correlation/sensitivity/reproduce_caps.py",
        inputs={"A": [1.1, 1.5, 2.0], "scales": [1, 2, 4, 8], "edges": "linspace(1,2,11)", "rounds": 2, "top_k": 3},
        seed=0, precision="float64", solver="SLSQP (positive) and CLARABEL via cvxpy (signed PSD)",
        solver_package="scipy")
    paths.OUTPUTS.mkdir(exist_ok=True)
    (paths.OUTPUTS / "caps_reproduction.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
