"""Confirm the in-class (teacher-student) transition with several random targets
and a larger restart budget, to distinguish expressibility vs trainability by depth."""
import time, numpy as np
from scipy.optimize import minimize
from pqc_common import ansatz_ops, unitary_from_params, _cost_grad, DIM
I32 = np.eye(DIM, dtype=complex)

def best_delta(target, ops, npar, restarts, seed):
    rng = np.random.default_rng(seed); best = np.inf
    for _ in range(restarts):
        x0 = rng.uniform(-np.pi, np.pi, npar)
        r = minimize(lambda x: _cost_grad(ops, x, target, I32), x0, jac=True,
                     method="L-BFGS-B", options=dict(maxiter=5000, ftol=1e-15, gtol=1e-13))
        best = min(best, float(r.fun))
    return best

for B, ntgt, R in [(12, 3, 30), (16, 4, 30), (20, 3, 30)]:
    ops, npar = ansatz_ops(B)
    ds = []
    t = time.time()
    for k in range(ntgt):
        thstar = np.random.default_rng(5000 + 13*B + k).uniform(-np.pi, np.pi, npar)
        Ustar = unitary_from_params(ops, thstar)
        ds.append(best_delta(Ustar, ops, npar, R, seed=90 + k))
    ds = np.array(ds)
    print(f"B={B:>2} (npar={npar}, {R} restarts): in-class delta over {ntgt} targets = "
          f"[{', '.join(f'{d:.1e}' for d in ds)}]  min={ds.min():.1e} "
          f"({time.time()-t:.0f}s)", flush=True)
