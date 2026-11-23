"""Regenerate tests/data/historical_failures.json from the frozen legacy code (deterministic).

Each entry is a concrete input that the pre-audit code accepted although it should not have.
The JSON is committed as permanent evidence; tests/test_historical_failures.py checks that the
shared core and both exact verifiers reject every entry, and that the legacy behaviour still
reproduces. Run from the repository root: python tests/data/make_historical_failures.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "legacy" / "pair_correlation_sensitivity_v0" / "src"))
sys.path.insert(0, str(ROOT / "src"))

import pair_sensitivity as legacy_pos  # noqa: E402
import signed_sensitivity as legacy_signed  # noqa: E402

from riemann_research.kernels import kernel_matrix, knots, window_matrix  # noqa: E402


def tolerance_candidate():
    """Search (seed 1, as in audit check4 part C) for a rank-one q that the legacy 2e-7 tolerance accepts."""
    prog = legacy_signed.make_program(2.0, 12)
    W = prog.model.widths
    n = W.size
    B = window_matrix(W)
    Ks = np.stack([kernel_matrix(W, k) for k in knots(W, 1.0, 2.0)])
    rng = np.random.default_rng(1)
    eps, best = 1.9e-7, None
    for _ in range(300):
        q0 = rng.normal(0, 1, n)
        if abs(q0.sum()) <= 0.2:
            continue
        q0 /= q0.sum()
        res = minimize(lambda q: q @ B @ q, q0, jac=lambda q: 2 * B @ q, method="SLSQP",
                       constraints=[{"type": "eq", "fun": lambda q: q.sum() - 1, "jac": lambda q: np.ones(n)},
                                    {"type": "ineq", "fun": lambda q: eps - np.einsum("i,kij,j->k", q, Ks, q),
                                     "jac": lambda q: -2 * np.einsum("kij,j->ki", Ks, q)}],
                       bounds=[(-8, 8)] * n, options={"maxiter": 1000, "ftol": 1e-14})
        q = res.x
        v = prog.validate_full_continuum(q, require_negative=True)
        if v < 2e-7 and abs(q.sum() - 1) < 2e-7:
            P = 2 - q @ B @ q
            if best is None or P > best[0]:
                best = (float(P), float(v), q)
    return {"widths": W.tolist(), "q": best[2].tolist(), "legacy_reported_P": best[0],
            "legacy_violation": best[1], "legacy_tolerance": 2e-7}


def main():
    out = {"generated_by": "tests/data/make_historical_failures.py", "source": "legacy/pair_correlation_sensitivity_v0"}
    m = legacy_pos.make_model(2.5, 12)
    e = legacy_pos.uniform_edges()
    fit = legacy_pos.optimize(m, e, np.diff(e))
    out["coverage_A2.5_positive"] = {"A": 2.5, "widths": m.widths.tolist(), "edges": e.tolist(),
                                     "caps": np.diff(e).tolist(), "q": fit["q"].tolist(),
                                     "legacy_reported_P": fit["P"], "correct_P_with_covering_bins": 0.935599}
    m = legacy_pos.make_model(1.5, 12)
    e = np.linspace(1.5, 2.0, 6)
    fit = legacy_pos.optimize(m, e, np.diff(e))
    out["partition_starting_at_1.5"] = {"A": 1.5, "widths": m.widths.tolist(), "edges": e.tolist(),
                                        "caps": np.diff(e).tolist(), "q": fit["q"].tolist(),
                                        "legacy_reported_P": fit["P"], "correct_P_with_bins_from_1": 0.847045}
    prog = legacy_signed.make_program(2.5, 12)
    so = legacy_signed.optimize_signed(prog, mode="sign_only", restarts=4)
    out["sign_only_A2.5"] = {"A": 2.5, "widths": prog.model.widths.tolist(), "q": so["q"].tolist(),
                             "legacy_valid": bool(so["valid"]), "legacy_reported_P": so["P"],
                             "legacy_reported_violation": so["continuum_violation"]}
    out["tolerance_accepted_positive_R"] = tolerance_candidate()
    out["sdp_cascade_sign_only_A2"] = {
        "evidence_file": "audits/2026-11-23-pair-correlation-audit/scripts/check2_results.json",
        "path": ["sign_only", "2.0"], "reported_P": 0.67290647,
        "note": "Solver-tolerance SDP solution: 2.6e-9 violation cascaded into weight ~0.19 on width 7/6."}
    path = ROOT / "tests" / "data" / "historical_failures.json"
    path.write_text(json.dumps(out, indent=1))
    print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
