"""Research B2-mini: solve formulations A and B at m = 0 and m = 1 and save raw outputs.

Run from a checkout:  python experiments/cgdl_reproduction/b2_mini/run_b2_mini.py
Outputs are numerical candidates, not certificates. See README.md in this directory.
"""
import json
import platform
import subprocess
import sys
from pathlib import Path

import clarabel
import cvxpy
import mpmath
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from riemann_research import cgdl_fourier as cf  # noqa: E402
from riemann_research import cgdl_sdp as sdp  # noqa: E402

OUT = Path(__file__).resolve().parent / "results" / "b2_mini_results.json"


def fourier_grid_minimum(q, upper=8.0, n=8001):
    """Exploratory-float: minimum of (T q)(u) on a grid. Not a proof of anything."""
    from fractions import Fraction
    tq = cf.transform_laguerre(cf.even_poly(Fraction(c) for c in q))
    with mpmath.workdps(30):
        grid = [upper * k / (n - 1) for k in range(n)]
        vals = [float(cf.poly_eval_mp(tq, u)) for u in grid]
    i = int(np.argmin(vals))
    return {"min_value": vals[i], "argmin_u": grid[i], "grid": f"{n} points on [0, {upper}]",
            "rigor": "exploratory-float"}


def main():
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                                text=True, check=True).stdout.strip()
    except Exception:  # noqa: BLE001
        commit = "unavailable"
    record = {
        "description": "B2-mini: formulations A and B, basis degree m in {0, 1} (deg_q_max in {1, 3}), rho = 1",
        "status": "numerical candidates only; no certificate; no registry change; assumes RH for any zeta use",
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "cvxpy": cvxpy.__version__, "clarabel": clarabel.__version__,
                        "mpmath": mpmath.__version__, "commit": commit},
        "a0": sdp.A0,
        "runs": [],
    }
    values = {}
    for F in "AB":
        for m in (0, 1):
            r = sdp.solve(F, m)
            entry = {"result": r.as_dict(), "deg_q_max": sdp.deg_q_max(m)}
            if r.objective is not None:
                entry["independent_checks"] = sdp.independent_checks(r)
                entry["fourier_grid_minimum"] = fourier_grid_minimum(r.q)
                values[F, m] = r.objective
            record["runs"].append(entry)
    record["comparison"] = {
        f"m={m}": {"V_A": values["A", m], "V_B": values["B", m], "V_A_minus_V_B": values["A", m] - values["B", m]}
        for m in (0, 1)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=2, default=str) + "\n", encoding="utf-8")
    for k, v in record["comparison"].items():
        print(k, v)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
