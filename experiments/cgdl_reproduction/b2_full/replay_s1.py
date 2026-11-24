"""Solver-free replay of the saved S1 witnesses (no CVXPY, no Clarabel).

    python experiments/cgdl_reproduction/b2_full/replay_s1.py [--fourier] [sweep.json]

For every witness in the witness file named by the sweep (default results/s1_sweep.json):
* the witness file's SHA-256 must equal the one recorded in the sweep;
* `replay_witness` rebuilds the affine data from (m, rho, basis), recomputes residuals, objective, q(0), (Tq)(0),
  R, primal and dual Gram eigenvalues, KKT stationarity, duality gap and complementarity with numpy, re-applies
  the acceptance thresholds and compares with the values recorded at solve time;
* the sweep's `best_accepted` for that degree must be the same point (rho, objective, diagnostics).
With --fourier, `highprec_check` is re-run from the saved Gram matrices and its Fourier-side minima are
compared with the sweep. Prints JSON; exit status 0 only if every check passes. Numerical checks only.
"""
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from riemann_research import cgdl_laguerre as L  # noqa: E402

DEFAULT_SWEEP = Path(__file__).resolve().parent / "results" / "s1_sweep.json"


def _reject_constant(name):
    raise ValueError(f"non-standard JSON constant {name}")


def load_strict(path: Path):
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=_reject_constant)


def replay(sweep_path: Path, fourier: bool = False) -> dict:
    sweep = load_strict(sweep_path)
    wf = sweep["witness_file"]
    wit_path = sweep_path.parent / wf["path"]
    raw = wit_path.read_bytes()
    report = {"sweep": str(sweep_path.name), "witness_file": wf["path"],
              "witness_sha256_matches": hashlib.sha256(raw).hexdigest() == wf["sha256"], "witnesses": []}
    witnesses = json.loads(raw.decode("utf-8"), parse_constant=_reject_constant)["witnesses"]
    report["witness_degrees_match"] = [w["m"] for w in witnesses] == wf["degrees"]
    runs = {e["m"]: e for e in sweep["runs"]}
    for w in witnesses:
        rp = L.replay_witness(w)
        entry = {"m": w["m"], "rho": w["rho"], "replay_ok": rp["ok"], "checks": rp["checks"], "errors": rp["errors"]}
        best = runs.get(w["m"], {}).get("best_accepted")
        entry["matches_sweep_best_accepted"] = bool(
            best is not None and best["rho"] == w["rho"]
            and all(best[k] == w["recorded"][k] for k in L.WITNESS_RECORDED_KEYS))
        if rp["ok"]:
            rc = rp["recomputed"]
            entry["recomputed"] = {k: rc[k] for k in ("objective", "residual_max_abs", "kkt_stationarity_max_abs",
                                                      "duality_gap", "complementarity", "R_physical")}
            entry["recomputed"]["gram_min_eig_min"] = min(rc["gram_min_eig"].values())
            entry["recomputed"]["dual_min_eig_min"] = min(rc["dual_min_eig"].values())
        if fourier:
            r = L.S1Result(w["m"], w["rho"], w["basis"], w["status"], w["solver"],
                           objective=w["recorded"]["objective"], grams=w["X"])
            f_new = L.highprec_check(r)["fourier"]
            f_old = runs[w["m"]]["highprec_check"]["fourier"]
            keys = [("unweighted_P", "min_value"), ("unweighted_P", "argmin_t"),
                    ("weighted_P_exp_minus_pi_t", "min_value_attained"), ("weighted_P_exp_minus_pi_t", "argmin_t")]
            entry["fourier_reproduced"] = f_new["reliable"] == f_old["reliable"] and all(
                f_new[a][b] == f_old[a][b] for a, b in keys)
        report["witnesses"].append(entry)
    report["ok"] = (report["witness_sha256_matches"] and report["witness_degrees_match"] and bool(witnesses)
                    and all(e["replay_ok"] and e["matches_sweep_best_accepted"] and e.get("fourier_reproduced", True)
                            for e in report["witnesses"]))
    return report


def main(argv):
    fourier = "--fourier" in argv
    paths = [a for a in argv if not a.startswith("--")]
    report = replay(Path(paths[0]) if paths else DEFAULT_SWEEP, fourier)
    print(json.dumps(report, indent=1, allow_nan=False))
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
