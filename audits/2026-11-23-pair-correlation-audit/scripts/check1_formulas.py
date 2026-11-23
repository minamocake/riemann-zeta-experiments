"""Check 1: formula-level verification of the audited kernels against independent routes."""
import sys
from fractions import Fraction as Fr
import numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__))
LEGACY_SRC = os.path.join(HERE, '..', '..', '..', 'legacy', 'pair_correlation_sensitivity_v0', 'src')
sys.path.insert(0, HERE)
sys.path.insert(1, LEGACY_SRC)
import indep as I
import pair_sensitivity as PS

# (a) Fejer / Montgomery: single box of width 1 => R = triangle, C = 1 + 2*int_0^1 a(1-a) = 4/3
C = I.C_window_exact([Fr(1)], [Fr(1)])
print('(a) Fejer exact C =', C, '=> P =', 2 - C)

# (b) audited overlap_kernel / weighted_kernel_integral vs exact rationals, random + adversarial widths
rng = np.random.default_rng(0)
maxerr_K = maxerr_I = 0.0
cases = [(1.0, 1.0), (0.5, 1.5), (2.0, 2.0), (0.1, 1.9), (1.0, 3.0), (2.5, 0.3), (1e-3, 1.0)]
cases += [tuple(rng.uniform(0.01, 3, 2)) for _ in range(300)]
for L, M in cases:
    for x in list(rng.uniform(0, 3.5, 20)) + [abs(L - M) / 2, (L + M) / 2, 0.0, 1.0]:
        e = float(I.K_exact(Fr(L), Fr(M), Fr(x)))
        maxerr_K = max(maxerr_K, abs(PS.overlap_kernel(L, M, x) - e))
    for up in [1.0, 0.37, 2.2]:
        e = float(I.xK_integral_exact(Fr(L), Fr(M), Fr(up)))
        maxerr_I = max(maxerr_I, abs(PS.weighted_kernel_integral(L, M, up) - e) / max(1, abs(e)))
print(f'(b) max |overlap_kernel - exact| = {maxerr_K:.2e};  max rel err weighted_kernel_integral = {maxerr_I:.2e}')

# (c) R(a) from K formula vs brute-force convolution of h (signed q)
W = [0.3, 0.7, 1.0, 1.4, 2.0]
q = np.array([0.8, -0.5, 1.2, -0.9, 0.4]); q = q / q.sum()
al = np.array([0.0, 0.2, 0.65, 1.0, 1.17, 1.5, 1.85, 2.0, 2.3])
Rk = np.array([q @ I.Kmat(W, a) @ q for a in al])
Rc = I.R_direct_convolution(q, W, al, dt=2e-5)
print('(c) max |R_K - R_conv| =', f'{np.max(np.abs(Rk - Rc)):.2e}', ' (grid dt=2e-5)')

# (d) C_window from frequency side (exact) vs physical side int f * mu_hat  (completely different route)
for W, q in [([1.0], [1.0]),
             ([0.3, 0.7, 1.0, 1.4, 2.0], list(q)),
             (list(np.linspace(1/12, 1, 12)), list(rng.dirichlet(np.ones(12))))]:
    Ce = float(I.C_window_exact([Fr(x) for x in q], [Fr(w) for w in W]))
    R0e = float(sum(Fr(a) * Fr(b) * I.K_exact(Fr(u), Fr(v), Fr(0)) for a, u in zip(q, W) for b, v in zip(q, W)))
    R0n, Cn = I.C_window_xspace(q, W)
    Qc = PS.make_model(max(max(W), 1.0), 2)  # just to access functions; build audited matrix directly:
    Qa = np.array([[PS.overlap_kernel(a, b, 0) + 2 * PS.weighted_kernel_integral(a, b) for b in W] for a in W])
    Ca = float(np.array(q) @ Qa @ np.array(q))
    print(f'(d) W={len(W)} boxes: exact C={Ce:.10f}  audited C={Ca:.10f}  x-space C={Cn:.7f}  (R0 exact {R0e:.8f} vs x-space {R0n:.7f})')

# (e) claim: for q>=0 R is nonincreasing on a>=0 (so band sup = R(left)); and for signed q it can increase
a = np.linspace(0, 2.2, 4001)
Wd = np.linspace(2/12, 2, 12)
worst_inc = 0
for _ in range(200):
    qq = rng.dirichlet(np.ones(12))
    R = np.array([qq @ I.Kmat(Wd, x) @ qq for x in a])
    worst_inc = max(worst_inc, np.max(np.diff(R)))
print(f'(e) nonneg q: largest increase of R on grid over 200 random q = {worst_inc:.2e} (should be <=0)')
