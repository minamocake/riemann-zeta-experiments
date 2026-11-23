"""Outside-window mass paid with a cap derived (not assumed) under RH.

Lemma (docs/derivations/int-F-cap.md, elementary, NOT peer reviewed):
    for every b, int_{b-1/2}^{b+1/2} F(a, T) da <= 8/3 + o(1)   (assuming RH).
Hence a single bin [1, S] (S <= 2) with cap 8/3 is admissible without any unproven hypothesis.
The SDP optimum over PSD mixtures of 12 boxes is rounded to an exact certificate.
Writes certificates/pc-provencap-A1.5.json and pc-provencap-A2.json. Requires cvxpy.
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_common"))
import paths  # noqa: E402

from riemann_research import certify, provenance  # noqa: E402
from riemann_research.montgomery import Partition  # noqa: E402
from riemann_research.optimize import sdp  # noqa: E402

CAP = Fraction(8, 3)


def widths_for(A: Fraction):
    w = [A * k / 12 for k in range(1, 13)]
    if Fraction(1) not in w:
        w[min(range(12), key=lambda i: abs(w[i] - 1))] = Fraction(1)
    return sorted(w)


def main():
    paths.CERTIFICATES.mkdir(exist_ok=True)
    for tag, A in [("A1.5", Fraction(3, 2)), ("A2", Fraction(2))]:
        W = widths_for(A)
        wf = [float(x) for x in W]
        part = Partition.build([1.0, float(A)], [float(CAP)], support=float(A))
        res = sdp(wf, "caps", part)
        cert = certify.build(f"pc-provencap-{tag}", W, res["X"], "caps", edges=[1, A], caps=[CAP],
                             notes="PSD mixture of centered boxes; single bin [1,S] with derived cap 8/3.",
                             provenance=provenance.collect(
                                 "experiments/pair_correlation/proven_cap/run_proven_cap.py",
                                 inputs={"widths": f"{A}*k/12, k=1..12 (closest to 1 replaced by 1)",
                                         "edges": ["1", str(A)], "caps": ["8/3"], "mode": "caps"},
                                 seed=None, precision={"denominator": 10**12, "claim_digits": 10},
                                 solver="CLARABEL via cvxpy (SDP)", solver_package="clarabel"))
        (paths.CERTIFICATES / f"pc-provencap-{tag}.json").write_text(json.dumps(cert, indent=1))
        print(f"{tag}: SDP status={res['status']} objective P={2 - res['objective']:.8f} "
              f"certified P >= {cert['claimed']['P_lower_bound_decimal']}  rank={len(cert['components'])}")


if __name__ == "__main__":
    main()
