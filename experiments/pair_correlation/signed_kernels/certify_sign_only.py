"""Exact certificates for the legacy sign-only candidates (R <= 0 on [1, A]).

The legacy optimiser returns float coefficients whose weights on widths > 1 are ~1e-12 and whose
correlation is ~+1e-16 at some knots, so they are not exactly admissible. By the degeneracy
proposition (docs/derivations/sign-only-degeneracy.md) those weights must be exactly zero; the
repair below sets them to zero and certifies the result. The change in P is below 1e-13.
Also writes a conditional (counterfactual-cap) certificate for the legacy A=2, scale=4 result,
to illustrate that exact arithmetic does not remove an unproven hypothesis.
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_common"))
import paths  # noqa: E402

sys.path.insert(1, str(paths.LEGACY_SRC))
import signed_sensitivity as legacy_signed  # noqa: E402

from riemann_research import certify, provenance  # noqa: E402


def main():
    paths.CERTIFICATES.mkdir(exist_ok=True)
    for A in [1.1, 1.5, 2.0]:
        prog = legacy_signed.make_program(A, 12)
        opt = legacy_signed.optimize_signed(prog, mode="sign_only", restarts=4)
        w = prog.model.widths
        q = np.where(w > 1 + 1e-12, 0.0, opt["q"])
        cert = certify.build(f"pc-legacy-signonly-A{A}", [Fraction(float(x)) for x in w], np.outer(q, q),
                             "sign-only",
                             notes="Legacy rank-one candidate; weights on widths > 1 set to exactly 0 "
                                   "(they were ~1e-12). Widths are the legacy binary floats.",
                             provenance=provenance.collect(
                                 "experiments/pair_correlation/signed_kernels/certify_sign_only.py",
                                 inputs={"A": A, "basis_count": 12, "mode": "sign-only",
                                         "source": "legacy signed_sensitivity.optimize_signed(restarts=4)"},
                                 seed=42, precision={"denominator": 10**12, "claim_digits": 10},
                                 solver="SLSQP via legacy signed_sensitivity", solver_package="scipy"))
        (paths.CERTIFICATES / f"pc-legacy-signonly-A{A}.json").write_text(json.dumps(cert, indent=1))
        print(f"A={A}: legacy float P={opt['P']:.12f}  certified P >= {cert['claimed']['P_lower_bound_decimal']}")

    edges = np.linspace(1.0, 2.0, 11)
    prog = legacy_signed.make_program(2.0, 12, edges, masses=4 * np.diff(edges))
    fit = legacy_signed.optimize_signed(prog, "caps", restarts=3, seed=42)
    cert = certify.build("pc-legacy-caps-A2-scale4", [Fraction(float(x)) for x in prog.model.widths],
                         np.outer(fit["q"], fit["q"]), "caps",
                         edges=[Fraction(float(e)) for e in edges],
                         caps=[Fraction(float(u)) for u in 4 * np.diff(edges)],
                         notes="Arithmetic is exact, but the caps U_i = 4 * width are an UNPROVEN hypothesis.",
                         provenance=provenance.collect(
                             "experiments/pair_correlation/signed_kernels/certify_sign_only.py",
                             inputs={"A": 2.0, "basis_count": 12, "mode": "caps", "edges": "linspace(1,2,11)",
                                     "caps": "4 * bin width", "source": "legacy optimize_signed(restarts=3)"},
                             seed=42, precision={"denominator": 10**12, "claim_digits": 10},
                             solver="SLSQP via legacy signed_sensitivity", solver_package="scipy"))
    (paths.CERTIFICATES / "pc-legacy-caps-A2-scale4.json").write_text(json.dumps(cert, indent=1))
    print(f"A=2 scale=4 (conditional): legacy float P={fit['P']:.10f}  "
          f"exact P >= {cert['claimed']['P_lower_bound_decimal']}")


if __name__ == "__main__":
    main()
