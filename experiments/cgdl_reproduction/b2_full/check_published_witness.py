"""Oracle: evaluate the S1 formulation-B model (paper basis) at CGdL's published d = 40 witness.

The witness file is NOT stored in this repository. Pass the path of Z-40.txt from arXiv:1810.08843v2/anc:
    python experiments/cgdl_reproduction/b2_full/check_published_witness.py /path/to/Z-40.txt
Mapping: CGdL's X1 (f side) -> our X2; CGdL's X2, X3 (fhat side) -> our X3, X4; rho = R (Gaussian width 1);
CGdL's d = our m = 40. Float64 evaluation of our affine maps; the witness itself was verified with Arb
in the separate S0 audit. Prints JSON.

The file's SHA-256 must equal EXPECTED_SHA256 (the digest of the owner-supplied file used for the S0 audit
and for results/published_witness_check.json); any other file is rejected before parsing (exit status 2).
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from riemann_research import cgdl_laguerre as L  # noqa: E402

EXPECTED_SHA256 = "a268555389798aa132a46d90057f2dfa28bc1dcf81d212282b057349d2cbd9e9"


def main(path):
    raw = Path(path).read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != EXPECTED_SHA256:
        print(f"rejected: SHA-256 {digest} of {path} does not match the expected {EXPECTED_SHA256}", file=sys.stderr)
        return 2
    lines = [x for x in raw.decode().split("\n") if x.strip()]
    R = float(lines[0])
    mats = [np.array([[float(v) for v in row.split()] for row in lines[k].strip()[1:-1].split(";")]) for k in (1, 2, 3)]
    m = mats[0].shape[0] - 1
    model = L.build_model(m, R, "paper")
    X = {"X2": mats[0], "X3": mats[1], "X4": mats[2]}
    resid = {name: float(sum(np.sum(M * X[b]) for b, M in d.items()) - rhs) for name, d, rhs in model.A}
    obj = float(sum(np.sum(M * X[b]) for b, M in model.C.items()))
    S2 = np.einsum("ij,ijk->k", X["X2"], L.product_tensor("paper", m))
    a = model.qmap @ S2
    w0 = np.array([float(L.w_at_zero(k)) for k in range(len(a))])
    out = {"file_sha256": digest, "sha256_verified": True, "m": m, "R": R,
           "max_abs_matching_residual": max(abs(v) for k, v in resid.items() if k.startswith("match")),
           "tq0_minus_1": resid["tq0"], "q0": float(w0 @ a), "objective_Z_inv": obj,
           "objective_minus_published_1.3207644973215569": obj - 1.3207644973215569,
           "min_eig": {b: float(np.linalg.eigvalsh(M).min()) for b, M in X.items()}}
    print(json.dumps(out, indent=1, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
