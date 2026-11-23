"""Baseline: classical Montgomery problem (Fourier support in [-1, 1]) in the nonnegative box family.

Widths k/12, k = 1..12 (exact rationals). With support <= 1 no outside-window term arises, so the
only assumptions are RH and Montgomery's theorem. Writes certificates/pc-baseline-A1.json.
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_common"))
import paths  # noqa: E402

from riemann_research import certify, provenance  # noqa: E402
from riemann_research.montgomery import Partition  # noqa: E402
from riemann_research.optimize import positive_qp  # noqa: E402


def main():
    widths = [Fraction(k, 12) for k in range(1, 13)]
    wf = [float(w) for w in widths]
    part = Partition.build([1.0, 2.0], [0.0], support=1.0)  # no outside mass: support is [-1, 1]
    fit = positive_qp(wf, part)
    cert = certify.build("pc-baseline-A1", widths, fit["q"][:, None] @ fit["q"][None, :],
                         "support-in-window",
                         notes="Nonnegative centered boxes of widths k/12; classical support condition.",
                         provenance=provenance.collect(
                             "experiments/pair_correlation/baselines/run_baseline.py",
                             inputs={"widths": "k/12, k=1..12", "family": "nonnegative centered boxes",
                                     "mode": "support-in-window"},
                             seed=0, precision={"denominator": 10**12, "claim_digits": 10},
                             solver="SLSQP via riemann_research.optimize.positive_qp", solver_package="scipy"))
    paths.CERTIFICATES.mkdir(exist_ok=True)
    (paths.CERTIFICATES / "pc-baseline-A1.json").write_text(json.dumps(cert, indent=1))
    print(f"float P = {fit['P']:.10f}  optimality = {fit['optimality'].value}  "
          f"certified P >= {cert['claimed']['P_lower_bound_decimal']}")


if __name__ == "__main__":
    main()
