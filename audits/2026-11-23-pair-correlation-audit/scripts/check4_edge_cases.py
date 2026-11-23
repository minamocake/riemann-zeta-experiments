"""Check 4: boundary conditions / counterexamples against the audited code, plus a proven-cap variant."""
import sys, warnings
from fractions import Fraction as Fr
import numpy as np
import cvxpy as cp
from scipy.optimize import minimize
warnings.filterwarnings('ignore')
import os
HERE = os.path.dirname(os.path.abspath(__file__))
LEGACY_SRC = os.path.join(HERE, '..', '..', '..', 'legacy', 'pair_correlation_sensitivity_v0', 'src')
sys.path.insert(0, HERE)
sys.path.insert(1, LEGACY_SRC)
import indep as I
from check3_exact import certify, rat
from check2_reproduce import widths_signed, sdp, EDGES0
import pair_sensitivity as PS
import signed_sensitivity as SS

print('=== (A) A > 2 with default edges [1,2]: out-of-window mass on (2,A] silently dropped ===')
A = 2.5
m = PS.make_model(A, 12)
fit = PS.optimize(m, PS.uniform_edges(), np.diff(PS.uniform_edges()))
e_ok = np.linspace(1, 2.5, 16)
fit_ok = PS.optimize(m, e_ok, np.diff(e_ok))
q = fit['q']
print(f'positive, A=2.5: code C={fit["C"]:.6f} (P={fit["P"]:.6f}) ; same hypothesis with bins covering [1,2.5]: C={fit_ok["C"]:.6f} (P={fit_ok["P"]:.6f})')
print(f'   R(2.0) for code optimum = {m.R(q, 2.0):.4e} > 0  -> contribution of (2,2.5] ignored, no error raised')
pr = SS.make_program(2.5, 12)
so = SS.optimize_signed(pr, mode='sign_only', restarts=4)
xs = np.linspace(2.0, 2.5, 2001)
Rs = np.array([pr.model.R(so['q'], x) for x in xs])
print(f'sign-only ("theorem-safe") A=2.5: code returns valid=True, P={so["P"]:.6f}, reported violation={so["continuum_violation"]:.1e}')
print(f'   but max R on (2,2.5] = {Rs.max():.4e} at a={xs[Rs.argmax()]:.3f}; the R<=0 requirement is unchecked there')
print(f'   with F==1 there (Montgomery conjecture) the omitted term 2*int_2^2.5 R = {2*np.trapezoid(Rs, xs):+.4e}')

print('\n=== (B) edges not starting at 1 are accepted: [1,1.5) silently dropped ===')
m = PS.make_model(1.5, 12)
e_bad = np.linspace(1.5, 2.0, 6)
fb = PS.optimize(m, e_bad, np.diff(e_bad))
fg = PS.optimize(m, PS.uniform_edges(), np.diff(PS.uniform_edges()))
print(f'A=1.5 edges from 1.5: C={fb["C"]:.6f} (P={fb["P"]:.6f})  vs edges from 1.0: C={fg["C"]:.6f} (P={fg["P"]:.6f})')
print(f'   (with edges from 1.5 the out-of-window term is identically 0: it is just the window functional)')

print('\n=== (C) tolerance-based "validity" in sign_only: violations up to 2e-7 accepted ===')
W = widths_signed(2.0, 12); n = 12
prg = SS.make_program(2.0, 12)
B = I.window_matrix(W)
kn = [k for k in I.knots_float(W, 1.0, 2.0)]
Ks = np.stack([I.Kmat(W, k) for k in kn])
best = None
rng = np.random.default_rng(1)
eps = 1.9e-7
for trial in range(300):
    q0 = rng.normal(0, 1, n); q0 /= q0.sum() if abs(q0.sum()) > 0.2 else 1
    if abs(q0.sum() - 1) > 1e-9:
        continue
    res = minimize(lambda q: q @ B @ q, q0, jac=lambda q: 2 * B @ q, method='SLSQP',
                   constraints=[{'type': 'eq', 'fun': lambda q: q.sum() - 1, 'jac': lambda q: np.ones(n)},
                                {'type': 'ineq', 'fun': lambda q: eps - np.einsum('i,kij,j->k', q, Ks, q),
                                 'jac': lambda q: -2 * np.einsum('kij,j->ki', Ks, q)}],
                   bounds=[(-8, 8)] * n, options={'maxiter': 1000, 'ftol': 1e-14})
    q = res.x
    v = prg.validate_full_continuum(q, require_negative=True)
    if v < 2e-7 and abs(q.sum() - 1) < 2e-7:
        P = 2 - q @ B @ q
        if best is None or P > best[0]:
            best = (P, v, q)
P, v, q = best
print(f'rank-1 q accepted by the code\'s own validity test (violation {v:.2e} < 2e-7): P = {P:.6f}')
print(f'   vs exactly-admissible optimum in this family (widths<=1 only): 0.670751 (rank 1) / 0.672295 (PSD mix)')
print(f'   weights on widths>1: {np.array2string(q[W > 1 + 1e-12], precision=3)}  -> large, not rounding noise')
print(f'   rigorous repair: pay R+ with the proven cap int_1^2 F <= 8/3 : P_rig >= {P - (16/3)*v:.6f}')

print('\n=== (D) proven-cap variant the code omits: single bin [1,2], U = 8/3 (RH only, no hypothesis) ===')
for Aw in [1.5, 2.0]:
    Wq = [Fr(Aw).limit_denominator(10) * k / 12 for k in range(1, 13)]
    if Fr(1) not in Wq:
        Wq[min(range(12), key=lambda i: abs(Wq[i] - 1))] = Fr(1); Wq.sort()
    Wf = np.array([float(w) for w in Wq])
    r = sdp(Wf, np.array([1.0, 2.0]), np.array([8 / 3]), Aw, 'caps')
    lam, V = np.linalg.eigh(r['X'])
    comps = [(rat(l, 10**12), [rat(x, 10**12) for x in V[:, i]]) for i, l in enumerate(lam) if l > 0]
    c = certify(comps, Wq, Fr(Aw).limit_denominator(10))
    Rplus = max(Fr(0), c['max_R_on_knots'])
    Pcert = c['P'] - Fr(16, 3) * Rplus
    print(f'A={Aw}: SDP P={r["P"]:.8f} rank~{np.sum(lam > 1e-6*lam.max())};  exact certificate: '
          f'P_window={float(c["P"]):.10f}, max R+ on [1,{Aw}]={float(Rplus):.3e}, P_cert={float(Pcert):.10f}')
print('   reference: Montgomery-Taylor/Cheer-Goldston 0.6725 ; Chirre-Goncalves-de Laat (2020) 0.6792')

print('\n=== (E) rank-one vs PSD at A=2, scale=8 (caps): is the audited local optimum global? ===')
prog = SS.make_program(2.0, 12, PS.uniform_edges(), masses=8 * np.diff(PS.uniform_edges()))
bestP = -1
for seed in range(40):
    try:
        o = SS.optimize_signed(prog, 'caps', restarts=6, seed=seed)
        bestP = max(bestP, o['P'])
    except RuntimeError:
        pass
r = sdp(widths_signed(2.0, 12), EDGES0, 8 * np.diff(EDGES0), 2.0, 'caps')
print(f'best rank-1 P over 40 seeds x 14 starts (audited optimizer) = {bestP:.8f} ; SDP (rank>1) P = {r["P"]:.8f}; '
      f'top eigenvalues {np.round(r["top_eigs"], 4)}')

print('\n=== (F) sensitivity = 2 max(0, sup_B R) vs re-optimised finite differences (global SDP), A=2 scale=4 ===')
W = widths_signed(2.0, 12); U = 4 * np.diff(EDGES0)
base = sdp(W, EDGES0, U, 2.0, 'caps')
lam, V = np.linalg.eigh(base['X']); qv = V[:, -1] * np.sqrt(lam[-1]); qv *= np.sign(qv.sum())
p2 = SS.make_program(2.0, 12, EDGES0, U)
S_formula = 2 * p2.exact_positive_max(qv)
h = 1e-4; fd = []
for i in range(10):
    up = U.copy(); up[i] += h; dn = U.copy(); dn[i] -= h
    fd.append((sdp(W, EDGES0, up, 2.0, 'caps')['C'] - sdp(W, EDGES0, dn, 2.0, 'caps')['C']) / (2 * h))
fd = np.array(fd)
print('formula :', np.array2string(S_formula, precision=5))
print('FD (SDP):', np.array2string(fd, precision=5), f'  max abs diff {np.max(np.abs(fd - S_formula)):.1e}')

print('\n=== (G) refinement strengthens the hypothesis: limit = "F <= scale a.e. on [1,2]" ===')
Wp = np.linspace(2 / 12, 2, 12)
al = np.linspace(1, 2, 20001)
Kint = np.trapezoid(np.stack([I.Kmat(Wp, a) for a in al]), al, axis=0)
Q = I.window_matrix(Wp) + 2 * Kint
qv = cp.Variable(12)
cp.Problem(cp.Minimize(cp.quad_form(qv, cp.psd_wrap((Q + Q.T) / 2))), [qv >= 0, cp.sum(qv) == 1]).solve()
print(f'positive family A=2, scale 1: stages 0/1/2 P = 0.907741 / 0.913049 / 0.914430 ; infinitely fine bins limit P = {2 - qv.value @ Q @ qv.value:.6f}')
