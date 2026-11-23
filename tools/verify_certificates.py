"""Re-check every certificate in certificates/ with two independent exact verifiers.

Usage: python tools/verify_certificates.py [--report PATH]
With --report, a JSON record of the verification results is written together with provenance
(commit, environment, verifier versions), so research runs keep their verification evidence.

A pass means the stated finite inequality holds exactly. It is not a proof of a statement about
zeta zeros: that deduction also needs RH, cited theorems and (for mode 'caps') the cap hypotheses.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from riemann_research import certify, independent_check, provenance  # noqa: E402


def main() -> int:
    report_path = sys.argv[sys.argv.index("--report") + 1] if "--report" in sys.argv else None
    records = []
    failures = 0
    for path in sorted((ROOT / "certificates").glob("*.json")):
        try:
            cert = json.loads(path.read_text())
            res = certify.verify(cert)
            other = independent_check.check(cert)
            if any(res[k] != other[k] for k in ("C", "P", "outside_term", "max_R_on_1_S")):
                raise certify.CertificateError("the two independent verifiers disagree")
            records.append({"certificate": path.name, "result": "passed", "P": str(res["P"]),
                            "P_decimal": float(res["P"]), "max_R_on_1_S": str(res["max_R_on_1_S"]),
                            "claimed": str(res["claimed"]), "mode": cert["outside"]["mode"]})
            print(f"arithmetic re-check passed (2 verifiers)  {path.name:36s} exact P = {float(res['P']):.12f}  "
                  f"(claimed >= {float(res['claimed']):.10f}, max R on [1,S] = {float(res['max_R_on_1_S']):.3e})")
        except (certify.CertificateError, ValueError) as exc:
            failures += 1
            records.append({"certificate": path.name, "result": "failed", "reason": str(exc)})
            print(f"arithmetic re-check FAILED  {path.name}: {exc}")
    print("Reminder: passing re-checks are not mathematical proofs; see docs/status-policy.md.")
    if report_path:
        report = {"provenance": provenance.collect("tools/verify_certificates.py", precision="exact rationals",
                                                   solver=None),
                  "verifiers": ["certify.verify", "independent_check.check"], "results": records,
                  "note": "Arithmetic re-checks only; not mathematical proofs."}
        Path(report_path).parent.mkdir(parents=True, exist_ok=True)
        Path(report_path).write_text(json.dumps(report, indent=1))
        print(f"wrote verification report {report_path}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
