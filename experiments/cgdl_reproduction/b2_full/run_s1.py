"""Research B2-full S1: formulation B in the Laguerre basis, bounded degree sweep with a rho search.

Run from a clean checkout:  python experiments/cgdl_reproduction/b2_full/run_s1.py [--quick]
Outputs numerical candidates only (no certificate, no registry change). See README.md.

Writes strict JSON (no NaN/Infinity). Per degree, `best_accepted` is the lowest objective among the solves
that pass `acceptable(r, ACCEPT_TOL)`, or null when none does; rejected solves are listed in
`rho_evaluations` and summarised in `lowest_rejected`, never reported as a candidate. For WITNESS_DEGREES the
primal and dual witness of `best_accepted` is written to s1_witnesses.json (replay: replay_s1.py).
"""
import hashlib
import json
import math
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
RESULTS = Path(__file__).resolve().parent / "results"
DEGREES = [1, 2, 3, 5, 8, 10, 12, 15, 20, 25, 30, 40]
RHO_GRID = [0.90, 0.95, 1.00, 1.05, 1.10, 1.20]
WITNESS_DEGREES = (1, 10, 30)           # representative low and high degree witnesses
ACCEPT_TOL = 1e-6
PUBLISHED_Z = 1.3207644973215569        # CGdL Z-40 witness (stored with +1e-6 objective slack)


def strict_json(obj, replaced, path="$"):
    """Copy of obj with every non-finite float replaced by None; the replaced paths are appended to `replaced`."""
    if isinstance(obj, dict):
        return {k: strict_json(v, replaced, f"{path}.{k}") for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [strict_json(v, replaced, f"{path}[{i}]") for i, v in enumerate(obj)]
    if isinstance(obj, (float, np.floating)):
        if not math.isfinite(obj):
            replaced.append(path)
            return None
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    return obj


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()


def main(quick=False):
    degrees = [1, 2, 3, 5] if quick else DEGREES
    rec = {"description": "B2-full S1: formulation B, paper Laguerre basis L_k^(-1/2)(pi t), rho search per degree",
           "status": "numerical candidates only; no certificate; no registry change; assumes RH for zeta use",
           "environment": {"python": platform.python_version(), "numpy": np.__version__, "cvxpy": cvxpy.__version__,
                           "clarabel": clarabel.__version__, "mpmath": mpmath.__version__,
                           "commit": git("rev-parse", "HEAD"), "worktree_clean": git("status", "--porcelain") == ""},
           "settings_ladder": [[n, o] for n, o in L.SETTINGS_LADDER],
           "acceptance": L.acceptance_thresholds(ACCEPT_TOL),
           "published_reference_Z": PUBLISHED_Z, "runs": []}
    witnesses = []
    for m in degrees:
        t0 = time.time()
        cache = {}

        def V(rho, m=m, cache=cache):
            key = round(rho, 12)
            if key not in cache:
                cache[key] = L.solve(m, rho, "paper")
            r = cache[key]
            return r.objective if L.acceptable(r, ACCEPT_TOL) else math.inf      # inf only inside the search

        grid = [(rho, V(rho)) for rho in RHO_GRID]
        accepted_grid = [g for g in grid if math.isfinite(g[1])]
        notes = []
        if accepted_grid:
            best_rho = min(accepted_grid, key=lambda t: t[1])[0]
            lo, hi = max(0.5, best_rho - 0.05), best_rho + 0.05
            L.golden_minimise(V, lo, hi, iters=10 if m >= 25 else 14)
        else:
            notes.append("no accepted solve on the rho grid: bracket search skipped")
        evals = sorted((k, c.objective, L.acceptable(c, ACCEPT_TOL), c.status) for k, c in cache.items())
        accepted = [(k, c) for k, c in cache.items() if L.acceptable(c, ACCEPT_TOL)]
        rejected = [(k, c) for k, c in cache.items() if not L.acceptable(c, ACCEPT_TOL)]
        best = hp = None
        if accepted:
            _, best = min(accepted, key=lambda kc: kc[1].objective)
            hp = L.highprec_check(best)
            if m in WITNESS_DEGREES:
                witnesses.append(L.witness_record(best, ACCEPT_TOL))
        with_obj = [(k, c) for k, c in rejected if c.objective is not None and math.isfinite(c.objective)]
        lowest_rejected = None
        if with_obj:
            k, c = min(with_obj, key=lambda kc: kc[1].objective)
            lowest_rejected = {"rho": k, "objective": c.objective, "status": c.status, "solver": c.solver,
                               "R_physical": c.R_physical, "note": "rejected by the acceptance rule; not a candidate"}
        entry = {"m": m, "deg_q_max": 2 * m + 1,
                 "rho_grid": [[rho, v if math.isfinite(v) else None] for rho, v in grid],
                 "rho_evaluations": evals, "n_accepted": len(accepted), "n_rejected": len(rejected),
                 "best_accepted": best.as_dict() if best else None, "highprec_check": hp,
                 "lowest_rejected": lowest_rejected,
                 "rho1": cache[1.0].as_dict() if 1.0 in cache else None,
                 "notes": notes, "seconds": round(time.time() - t0, 1)}
        rec["runs"].append(entry)
        if best:
            f = hp["fourier"]
            print(f"m={m:2d} rho*={best.rho:.6f} V={best.objective!r} R={best.R_physical} "
                  f"minP={f['unweighted_P']['min_value']} "
                  f"min(P e^-pi t)={f['weighted_P_exp_minus_pi_t']['infimum']:.2e} "
                  f"reliable={f['reliable']} ({entry['seconds']} s)", flush=True)
        else:
            print(f"m={m:2d} no accepted solve; lowest rejected {lowest_rejected} ({entry['seconds']} s)", flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    suffix = "_quick" if quick else ""
    wit_path = RESULTS / f"s1_witnesses{suffix}.json"
    replaced_w = []
    wit_doc = strict_json({"description": "Primal (X2, X3, X4) and dual (eq_duals, psd_duals) witnesses of selected "
                                          "accepted S1 candidates, with the diagnostics recorded at solve time. "
                                          "Numerical data, not certificates. Replay: replay_s1.py.",
                           "commit": rec["environment"]["commit"], "witnesses": witnesses}, replaced_w)
    if replaced_w:
        raise RuntimeError(f"non-finite values in witnesses: {replaced_w}")
    wit_bytes = (json.dumps(wit_doc, indent=1, allow_nan=False) + "\n").encode("utf-8")
    wit_path.write_bytes(wit_bytes)
    rec["witness_file"] = {"path": wit_path.name, "sha256": hashlib.sha256(wit_bytes).hexdigest(),
                           "degrees": [w["m"] for w in witnesses]}
    replaced = []
    doc = strict_json(rec, replaced)
    doc["nonfinite_values_replaced_by_null"] = replaced
    out = RESULTS / f"s1_sweep{suffix}.json"
    out.write_text(json.dumps(doc, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT)} and {wit_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
