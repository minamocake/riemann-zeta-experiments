"""Check 2: independent reproduction of the main numbers + global-optimality certificates.

* nonnegative-box QP: convexity on the simplex's affine hull checked with an ORTHONORMAL tangent basis,
  then solved by an independent convex solver (cvxpy/CLARABEL) -> global optimum within the family.
* signed boxes: Shor/SDP lift X = q q^T  ->  X >= 0 (PSD), 1^T X 1 = 1.  Any PSD X is itself an
  admissible test function f = v^T X v = sum_k lam_k |H_k|^2 >= 0, f(0)=1, supp f^ in [-A,A]; so the
  SDP optimum is (i) a lower bound for the rank-one problem the audited code solves and (ii) a valid
  (often better) member of the same Montgomery framework.
"""
import sys, json
import numpy as np
import cvxpy as cp
import os
HERE = os.path.dirname(os.path.abspath(__file__))
LEGACY_SRC = os.path.join(HERE, '..', '..', '..', 'legacy', 'pair_correlation_sensitivity_v0', 'src')
sys.path.insert(0, HERE)
import indep as I

EDGES0 = np.linspace(1.0, 2.0, 11)


def widths_positive(A, n):
    return np.linspace(A / n, A, n)


def widths_signed(A, n):
    w = np.linspace(A / n, A, n)
    if not np.any(np.isclose(w, 1.0, rtol=0, atol=1e-12)):
        w[np.argmin(abs(w - 1.0))] = 1.0
    return np.sort(w)


def refine(edges, sel):
    out = [edges[0]]
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        if i in sel:
            out.append((a + b) / 2)
        out.append(b)
    return np.array(out)


# ---------------- nonnegative family ----------------
def positive_global(W, edges, U):
    B = I.window_matrix(W)
    Q = B + sum(2 * u * I.Kmat(W, l) for l, u in zip(edges[:-1], U))
    Q = (Q + Q.T) / 2
    n = len(W)
    # orthonormal basis of {sum q = 0}
    Z = np.linalg.qr(np.c_[np.ones(n), np.eye(n)[:, :n - 1]])[0][:, 1:]
    lam = np.linalg.eigvalsh(Z.T @ Q @ Z).min()
    q = cp.Variable(n)
    prob = cp.Problem(cp.Minimize(cp.quad_form(q, cp.psd_wrap(Q))), [q >= 0, cp.sum(q) == 1])
    prob.solve(solver='CLARABEL')
    qv = np.maximum(q.value, 0); qv /= qv.sum()
    # KKT check: gradient g = 2Qq ; on support g_i = mu, off support g_i >= mu
    g = 2 * Q @ qv
    supp = qv > 1e-8
    mu = g[supp].mean()
    kkt = max(np.ptp(g[supp]), max(0, mu - g[~supp].min()) if (~supp).any() else 0)
    return float(qv @ Q @ qv), qv, lam, kkt


def positive_adaptive(A, n=12, rounds=2, top_k=3):
    W = widths_positive(A, n); edges = EDGES0.copy(); out = []
    for s in range(rounds + 1):
        U = np.diff(edges)
        C, q, lam, kkt = positive_global(W, edges, U)
        sens = np.array([2 * q @ I.Kmat(W, l) @ q for l in edges[:-1]])
        out.append(dict(stage=s, bins=len(U), C=C, P=2 - C, lam_min_orth=lam, kkt=kkt))
        ranked = [int(i) for i in np.argsort(-sens, kind='stable') if sens[i] > 1e-10]
        edges = refine(edges, ranked[:top_k])
    return out


# ---------------- signed family via SDP ----------------
def band_knots(W, edges, A):
    rows = []
    for j, (l, r) in enumerate(zip(edges[:-1], edges[1:])):
        if l >= A - 1e-12:
            continue
        for k in I.knots_float(W, l, min(r, A)):
            rows.append((j, k))
    return rows


def sdp(W, edges, U, A, mode):
    n = len(W); B = I.window_matrix(W)
    X = cp.Variable((n, n), PSD=True)
    one = np.ones(n)
    cons = [one @ X @ one == 1]
    rows = band_knots(W, edges, A)
    if mode == 'sign_only':
        cons += [cp.trace(I.Kmat(W, k) @ X) <= 0 for _, k in rows]
        obj = cp.trace(B @ X)
    else:
        t = cp.Variable(len(U), nonneg=True)
        cons += [t[j] >= cp.trace(I.Kmat(W, k) @ X) for j, k in rows]
        obj = cp.trace(B @ X) + 2 * np.asarray(U) @ t
    prob = cp.Problem(cp.Minimize(obj), cons)
    prob.solve(solver='CLARABEL', tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10, max_iter=500)
    Xv = (X.value + X.value.T) / 2
    ev = np.sort(np.linalg.eigvalsh(Xv))[::-1]
    # exact value of the (possibly higher-rank) test function, recomputed independently of the solver objective
    Rk = {k: float(np.trace(I.Kmat(W, k) @ Xv)) for _, k in rows}
    if mode == 'sign_only':
        val = float(np.trace(B @ Xv)); viol = max([0.0] + list(Rk.values()))
    else:
        mx = np.zeros(len(U))
        for j, k in rows:
            mx[j] = max(mx[j], Rk[k])
        val = float(np.trace(B @ Xv) + 2 * np.asarray(U) @ mx); viol = 0.0
    return dict(status=prob.status, C=val, P=2 - val, eig_ratio=float(ev[1] / ev[0]) if ev[0] > 0 else None,
                top_eigs=[float(e) for e in ev[:3]], viol=viol, X=Xv)


if __name__ == '__main__':
    res = {}
    print('=== nonnegative boxes, independent global solve (U_i = bin width, counterfactual) ===')
    C1, q1, lam1, kkt1 = positive_global(widths_positive(1.0, 12), EDGES0, np.diff(EDGES0))
    print(f'A=1.0 baseline: C={C1:.9f} P={2-C1:.9f}  lam_min(orthonormal tangent)={lam1:.3e} KKT resid={kkt1:.1e}')
    res['pos_A1'] = dict(C=C1, q=q1.tolist())
    for A in [1.1, 1.5, 2.0]:
        for r in positive_adaptive(A):
            print(f"A={A} stage={r['stage']} bins={r['bins']} C={r['C']:.6f} P={r['P']:.6f} "
                  f"lam_min_orth={r['lam_min_orth']:.3e} KKT={r['kkt']:.1e}")

    print('\n=== signed boxes: SDP (global over PSD-mixtures) vs audited rank-one local results ===')
    audited = {(1.1, 1): .72636305, (1.1, 2): .71615702, (1.1, 4): .69956055, (1.1, 8): .67697978,
               (1.5, 1): .84704480, (1.5, 2): .77241378, (1.5, 4): .70935274, (1.5, 8): .67941026,
               (2.0, 1): .90774070, (2.0, 2): .77436457, (2.0, 4): .71239311, (2.0, 8): .68052404}
    res['caps'] = {}
    for A in [1.1, 1.5, 2.0]:
        W = widths_signed(A, 12)
        for s in [1, 2, 4, 8]:
            r = sdp(W, EDGES0, s * np.diff(EDGES0), A, 'caps')
            print(f"A={A} scale={s}: SDP P={r['P']:.8f}  audited P={audited[(A,s)]:.8f}  "
                  f"diff={100*(r['P']-audited[(A,s)]):+.5f}pp  eig2/eig1={r['eig_ratio']:.1e} [{r['status']}]")
            res['caps'][f'{A}_{s}'] = dict(P_sdp=r['P'], P_audited=audited[(A, s)], eig_ratio=r['eig_ratio'])

    print('\n=== sign-only (RH-only, no hypothesis): SDP vs audited ===')
    aud_so = {1.1: 0.6724406356261092, 1.5: 0.6723649809728278, 2.0: 0.6707507066269873}
    res['sign_only'] = {}
    for A in [1.1, 1.5, 2.0]:
        W = widths_signed(A, 12)
        r = sdp(W, EDGES0, np.diff(EDGES0), A, 'sign_only')
        print(f"A={A}: SDP P={r['P']:.8f}  audited P={aud_so[A]:.8f}  eig2/eig1={r['eig_ratio']:.1e} "
              f"max R+ at knots={r['viol']:.1e} [{r['status']}]")
        res['sign_only'][str(A)] = dict(P_sdp=r['P'], eig_ratio=r['eig_ratio'], X=r['X'].tolist(), W=W.tolist())
    json.dump(res, open(os.path.join(HERE, 'check2_results.json'), 'w'), indent=1)
