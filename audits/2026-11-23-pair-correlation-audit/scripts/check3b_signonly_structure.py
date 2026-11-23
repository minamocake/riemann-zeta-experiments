"""Check 3b: sign-only mode degenerates to the classical support-in-[-1,1] problem (box family).

Claim (audit): if R = tr(K X) <= 0 on [1, A] with X PSD, then X_kk = 0 for every width w_k > 1.
Proof sketch: let w_k be the largest width >1 with X_kk>0; on ((w_k+w_{k-1})/2, w_k) only the pair (k,k)
has (w_i+w_j)/2 > alpha, so R(alpha) = X_kk (w_k - alpha)/w_k^2 > 0 there; contradiction. PSD => row k = 0.
Numerical consequences tested here:
  (i)  audited sign-only q: weights on widths > 1 are ~0 (only rounding);
  (ii) zeroing them exactly gives an EXACTLY admissible test function with P equal to the audited value;
  (iii) SDP restricted to widths <= 1 reproduces the audited values (so they are global in the family);
  (iv) the SDP 'gain' at A=2 seen in check2 disappears when constraints are enforced exactly.
"""
import sys
from fractions import Fraction as Fr
import numpy as np
import cvxpy as cp
import os
HERE = os.path.dirname(os.path.abspath(__file__))
LEGACY_SRC = os.path.join(HERE, '..', '..', '..', 'legacy', 'pair_correlation_sensitivity_v0', 'src')
sys.path.insert(0, HERE)
sys.path.insert(1, LEGACY_SRC)
import indep as I
from check3_exact import certify, report, rat
import signed_sensitivity as SS

for A in [1.1, 1.5, 2.0]:
    pr = SS.make_program(A, 12)
    opt = SS.optimize_signed(pr, mode='sign_only', restarts=4)
    w = pr.model.widths; q = opt['q']
    big = w > 1 + 1e-12
    print(f'A={A}: widths>1 -> {np.round(w[big],4).tolist()}  audited q there = {np.array2string(q[big], precision=2)}')
    # (ii) exact repair
    qq = [Fr(0) if b else rat(x) for x, b in zip(q, big)]
    c = certify([(Fr(1), qq)], [Fr(x) for x in w], Fr(w.max()))
    report(f'   repaired (q=0 on widths>1)', c)
    print(f'   audited P = {opt["P"]:.12f};  repaired exact P - audited = {float(c["P"]) - opt["P"]:+.2e}')
    # (iii) SDP restricted to widths <= 1 (support of R in [-1,1]; no outside constraint needed)
    ws = w[~big]; n = len(ws); B = I.window_matrix(ws)
    X = cp.Variable((n, n), PSD=True)
    cp.Problem(cp.Minimize(cp.trace(B @ X)), [np.ones(n) @ X @ np.ones(n) == 1]).solve(solver='CLARABEL')
    ev = np.sort(np.linalg.eigvalsh(X.value))[::-1]
    print(f'   SDP over {n} boxes of width<=1: P = {2 - np.trace(B @ X.value):.10f}  (eig2/eig1={ev[1]/ev[0]:.1e})')

# (iv) the SDP 'improvement' at A=2 from check2, examined exactly
print('\nSDP solution from check2 (A=2, sign-only), exact evaluation after PSD rounding:')
import json
d = json.load(open(os.path.join(HERE, 'check2_results.json')))['sign_only']['2.0']
Xv = np.array(d['X']); W = np.array(d['W'])
lam, V = np.linalg.eigh(Xv)
comps = [(rat(l, 10**12), [rat(x, 10**12) for x in V[:, i]]) for i, l in enumerate(lam) if l > 0]
c = certify(comps, [Fr(x) for x in W], Fr(2))
report('   SDP X (check2) exact', c)
print('   diag(X) on widths>1:', np.array2string(np.diag(Xv)[W > 1 + 1e-12], precision=2))
# How much would the outside positive part cost even under the most optimistic proven-type cap (U=1)?
kn = I.knots_float(W, 1.0, 2.0)
Rpos = max(0, max(np.trace(I.Kmat(W, k) @ Xv) for k in kn))
print(f'   max R+ on [1,2] = {Rpos:.2e};  required |offdiag| from tiny diag is what produces the spurious gain')
