"""Research B2-full S1: formulation B in the Laguerre basis, bounded degree sweep with a rho search.

Run from a checkout:  python experiments/cgdl_reproduction/b2_full/run_s1.py [--quick]
Outputs numerical candidates only (no certificate, no registry change). See README.md.
"""
import json
import platform
import subprocess
import sys
import time
import warnings
from pathlib import Path

import clarabel
import cvxpy
import mpmath
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from riemann_research import cgdl_laguerre as L  # noqa: E402

warnings.filterwarnings("ignore", message="Solution may be inaccurate")
OUT = Path(__file__).resolve().parent / "results" / "s1_sweep.json"
DEGREES = [1, 2, 3, 5, 8, 10, 12, 15, 20, 25, 30, 40]
RHO_GRID = [0.90, 0.95, 1.00, 1.05, 1.10, 1.20]
PUBLISHED_Z = 1.3207644973215569        # CGdL Z-40 witness (stored with +1e-6 objective slack)


def main(quick=False):
    degrees = [1, 2, 3, 5] if quick else DEGREES
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    rec = {"description": "B2-full S1: formulation B, paper Laguerre basis L_k^(-1/2)(pi t), rho search per degree",
           "status": "numerical candidates only; no certificate; no registry change; assumes RH for zeta use",
           "environment": {"python": platform.python_version(), "numpy": np.__version__, "cvxpy": cvxpy.__version__,
                           "clarabel": clarabel.__version__, "mpmath": mpmath.__version__, "commit": commit},
           "settings_ladder": [[n, o] for n, o in L.SETTINGS_LADDER],
           "published_reference_Z": PUBLISHED_Z, "runs": []}
    for m in degrees:
        t0 = time.time()
        cache = {}

        def V(rho, m=m, cache=cache):
            key = round(rho, 12)
            if key not in cache:
                cache[key] = L.solve(m, rho, "paper")
            r = cache[key]
            return r.objective if (r.objective is not None and L.acceptable(r, 1e-6)) else float("inf")

        grid = [(rho, V(rho)) for rho in RHO_GRID]
        best_rho = min(grid, key=lambda t: t[1])[0]
        lo, hi = max(0.5, best_rho - 0.05), best_rho + 0.05
        (rho_star, v_star), _ = L.golden_minimise(V, lo, hi, iters=10 if m >= 25 else 14)
        r = cache[round(rho_star, 12)]
        hp = L.highprec_check(r) if r.objective is not None else None
        evals = sorted((k, c.objective, L.acceptable(c, 1e-6), c.status) for k, c in cache.items())
        entry = {"m": m, "deg_q_max": 2 * m + 1, "rho_grid": grid, "rho_evaluations": evals,
                 "best": r.as_dict(), "best_acceptable_1e-7": L.acceptable(r),
                 "best_acceptable_1e-6": L.acceptable(r, 1e-6),
                 "highprec_check": hp, "rho1": cache.get(1.0).as_dict() if 1.0 in cache else None,
                 "seconds": round(time.time() - t0, 1)}
        rec["runs"].append(entry)
        print(f"m={m:2d} rho*={rho_star:.6f} V={v_star!r} R={r.R_physical} acc={L.acceptable(r)} "
              f"fhat_min={hp['fhat_min_on_grid'] if hp else None:.2e} ({entry['seconds']} s)", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out = OUT if not quick else OUT.with_name("s1_sweep_quick.json")
    out.write_text(json.dumps(rec, indent=1, default=str) + "\n", encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
