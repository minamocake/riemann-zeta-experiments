"""Check 5: cross-check the proven-cap certificate (A=2) by an independent route + knot-sufficiency + extrema."""
import sys, warnings
from fractions import Fraction as Fr
import numpy as np
warnings.filterwarnings('ignore')
import os
HERE = os.path.dirname(os.path.abspath(__file__))
LEGACY_SRC = os.path.join(HERE, '..', '..', '..', 'legacy', 'pair_correlation_sensitivity_v0', 'src')
sys.path.insert(0, HERE)
sys.path.insert(1, LEGACY_SRC)
import indep as I
from check2_reproduce import sdp
import pair_sensitivity as PS
import signed_sensitivity as SS

Wq = [Fr(k, 6) for k in range(1, 13)]
W = np.array([float(w) for w in Wq])
r = sdp(W, np.array([1.0, 2.0]), np.array([8 / 3]), 2.0, 'caps')
X = r['X']; X /= np.ones(12) @ X @ np.ones(12)
lam, V = np.linalg.eigh(X)
# physical-side window functional for f = sum_k lam_k H_k^2
Cx = 0.0; R0x = 0.0
for l, v in zip(lam, V.T):
    if l > 1e-14:
        r0, c = I.C_window_xspace(list(v), list(W))
        Cx += l * c; R0x += l * r0
Cf = float(np.trace(I.window_matrix(W) @ X))
print(f'window functional: frequency side {Cf:.9f}  physical side int f*mu_hat {Cx:.7f}')
# f >= 0 sanity and f(0)=1
xs = np.linspace(-40, 40, 200001)
fx = sum(l * (sum(vj * np.sinc(wj * xs) for vj, wj in zip(v, W))) ** 2 for l, v in zip(lam, V.T) if l > 0)
print(f'f(0) = {fx[len(xs)//2]:.12f}, min f on grid = {fx.min():.3e}')
# knot sufficiency: dense grid max of R on [1,2] vs knot max
al = np.linspace(1, 2, 200001)
Rg = np.array([np.trace(I.Kmat(W, a) @ X) for a in al[::50]])
kn = I.knots_float(W, 1.0, 2.0)
Rk = np.array([np.trace(I.Kmat(W, k) @ X) for k in kn])
print(f'max R on [1,2]: dense grid {Rg.max():.6e} | knots {Rk.max():.6e}  (grid must not exceed knots)')

# audited exact_band_extrema vs dense grid, random signed q, several bands
m = SS.signed_model(2.0, 12); rng = np.random.default_rng(3); worst = 0
for _ in range(50):
    q = rng.normal(size=12); q /= q.sum()
    for a, b in [(1.0, 1.1), (1.23, 1.71), (0.0, 2.0), (1.9, 2.0)]:
        ex = PS.exact_band_extrema(m, q, a, b)
        g = np.linspace(a, b, 3001); vals = np.array([m.R(q, x) for x in g])
        worst = max(worst, vals.max() - ex['max'], ex['min'] - vals.min())
print(f'audited exact_band_extrema: largest amount a dense grid beats it = {worst:.2e} (<=0 means correct)')
