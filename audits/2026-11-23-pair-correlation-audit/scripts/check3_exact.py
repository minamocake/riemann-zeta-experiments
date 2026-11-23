"""Check 3: exact-rational certificates for the RH-only (hypothesis-free) claims.

For a test function given by PSD X (rank one: X = q q^T) on widths W:
  admissible  <=>  X PSD, 1^T X 1 = 1, and R(a) = tr(K(a) X) <= 0 at every knot in [1, A]
  (R is piecewise linear with knots |L-M|/2, (L+M)/2, so knot-wise <= 0 is equivalent to <= 0 on [1,A];
   R == 0 for a >= A.)   Then, under RH, N_s >= (2 - C + o(1)) N with C = tr(B X) computed exactly.
"""
import sys, json
from fractions import Fraction as Fr
import numpy as np
import cvxpy as cp
import os
HERE = os.path.dirname(os.path.abspath(__file__))
LEGACY_SRC = os.path.join(HERE, '..', '..', '..', 'legacy', 'pair_correlation_sensitivity_v0', 'src')
sys.path.insert(0, HERE)
sys.path.insert(1, LEGACY_SRC)
import indep as I


def rat(x, den=10**12):
    return Fr(x).limit_denominator(den)


def certify(comps, W, A):
    """comps: list of (lam, v) with lam>=0 rational, v rational vectors. X = sum lam v v^T (PSD by construction)."""
    n = len(W)
    s = sum(l * sum(v) ** 2 for l, v in comps)          # 1^T X 1
    comps = [(l / s, v) for l, v in comps]               # exact normalisation f(0)=1
    X = [[sum(l * v[i] * v[j] for l, v in comps) for j in range(n)] for i in range(n)]
    def R(a):
        return sum(X[i][j] * I.K_exact(W[i], W[j], a) for i in range(n) for j in range(n))
    kn = I.knots_exact(W, Fr(1), A)
    vals = [(k, R(k)) for k in kn]
    worst = max(vals, key=lambda kv: kv[1])
    r0 = sum(X[i][j] * I.K_exact(W[i], W[j], Fr(0)) for i in range(n) for j in range(n))
    integ = sum(X[i][j] * I.xK_integral_exact(W[i], W[j]) for i in range(n) for j in range(n))
    C = r0 + 2 * integ
    return dict(C=C, P=2 - C, max_R_on_knots=worst[1], argmax=worst[0], n_knots=len(kn),
                admissible=worst[1] <= 0)


def report(tag, c):
    print(f"{tag}: admissible(exact)={c['admissible']}  max R at {c['n_knots']} knots = {float(c['max_R_on_knots']):.3e} "
          f"(at a={float(c['argmax']):.6f})   P = 2 - C = {float(c['P']):.10f}")


if __name__ == '__main__':
    import pair_sensitivity as PS
    import signed_sensitivity as SS

    # --- (1) audited A=1 nonnegative baseline: R supported in [-1,1] so admissible automatically
    m = PS.make_model(1.0, 12)
    fit = PS.optimize(m, PS.uniform_edges(), np.diff(PS.uniform_edges()))
    W = [Fr(w) for w in m.widths]          # the exact binary floats used by the code
    q = [rat(x) for x in fit['q']]
    report('audited A=1.0 baseline (rank 1, q>=0)', certify([(Fr(1), q)], W, Fr(1)))

    # --- (2) audited sign-only candidates (rank one, signed q)
    for A in [1.1, 1.5, 2.0]:
        pr = SS.make_program(A, 12)
        opt = SS.optimize_signed(pr, mode='sign_only', restarts=4)
        W = [Fr(w) for w in pr.model.widths]
        q = [rat(x) for x in opt['q']]
        c = certify([(Fr(1), q)], W, Fr(A) if A != 1.1 else Fr(pr.model.widths.max()))
        report(f'audited sign-only A={A} (rank 1)', c)

    # --- (3) audit's own improved candidate: SDP with a safety margin, then exact rounding
    print('\nSDP-based candidates on exact rational widths k*A/12 (my construction, not from the audited code):')
    out = {}
    for A, Aq in [(1.5, Fr(3, 2)), (2.0, Fr(2))]:
        Wq = [Aq * k / 12 for k in range(1, 13)]
        if Fr(1) not in Wq:
            Wq[min(range(12), key=lambda i: abs(Wq[i] - 1))] = Fr(1); Wq.sort()
        Wf = np.array([float(w) for w in Wq])
        n = 12
        B = I.window_matrix(Wf)
        kn = [k for k in I.knots_float(Wf, 1.0, float(Aq)) if k < float(Aq) - 1e-12]
        for margin in [1e-7, 1e-6, 1e-5]:
            X = cp.Variable((n, n), PSD=True)
            cons = [np.ones(n) @ X @ np.ones(n) == 1]
            cons += [cp.trace(I.Kmat(Wf, k) @ X) <= -margin * (float(Aq) - k) for k in kn]
            cp.Problem(cp.Minimize(cp.trace(B @ X)), cons).solve(solver='CLARABEL')
            Xv = (X.value + X.value.T) / 2
            lam, V = np.linalg.eigh(Xv)
            comps = [(rat(l, 10**9), [rat(x, 10**9) for x in V[:, i]]) for i, l in enumerate(lam) if l > 1e-12]
            c = certify(comps, Wq, Aq)
            report(f'  A={A} margin={margin:g} rank={len(comps)}', c)
            if c['admissible']:
                out[str(A)] = dict(P=str(c['P']), P_float=float(c['P']), rank=len(comps), margin=margin)
                break
    json.dump(out, open(os.path.join(HERE, 'check3_certified.json'), 'w'), indent=1)
