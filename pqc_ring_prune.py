"""
CNOT pruning of the learned 18-CNOT gadget compilation.
=======================================================

Start from the exact 18-CNOT solution (pqc_ring_L3_params.npy) and greedily remove
CNOTs one at a time — keeping ALL rotation layers so the parameter vector and its
warm-start stay aligned — retraining from the current solution after each removal.
A CNOT is "removable" if the circuit still compiles U to delta < 1e-6.  This finds a
(locally) minimal CNOT count for this ansatz family and tests whether we can go below
the structured Step-4 count of 14.

Run:  python pqc_ring_prune.py
"""
import json
import numpy as np
from scipy.optimize import minimize
from pqc_ring_ansatz import GADGET_PAIRS, cost_grad, unitary, _I32, _obs_err, PARAMS_FILE
from pqc_common import U_TARGET, DIM

SEQ = GADGET_PAIRS * 3          # 18 ordered CNOT pairs
NPAR = 15 * (1 + len(SEQ))      # 285: initial rot + one rot layer after each of 18 slots


def ansatz_masked(mask):
    """Rotation layers always present (initial + one after each slot); CNOT at slot i
    only if mask[i].  Param indexing is independent of mask (285 params)."""
    ops, p = [], 0

    def rotlayer():
        nonlocal p
        for w in range(5):
            for k in ("rx", "ry", "rz"):
                ops.append(("g", k, w, p)); p += 1
    rotlayer()
    for i, pair in enumerate(SEQ):
        if mask[i]:
            ops.append(("cnot", pair[0], pair[1]))
        rotlayer()
    return ops


def _delta(ops, x):
    return cost_grad(ops, x, U_TARGET, _I32)[0]


def train_masked(mask, warm, rng, n_random=1, maxiter=6000):
    ops = ansatz_masked(mask)
    starts = [warm] + [rng.uniform(-np.pi, np.pi, NPAR) for _ in range(n_random)]
    best = (np.inf, None)
    for x0 in starts:
        r = minimize(lambda x: cost_grad(ops, x, U_TARGET, _I32), x0, jac=True,
                     method="L-BFGS-B", options=dict(maxiter=maxiter, ftol=1e-15, gtol=1e-13))
        if r.fun < best[0]:
            best = (float(r.fun), r.x)
    return best


def main(tol=1e-6):
    theta = np.load(PARAMS_FILE)
    mask = [True] * len(SEQ)
    d0 = _delta(ansatz_masked(mask), theta)
    print(f"start: {sum(mask)} CNOTs, delta = {d0:.2e}\n", flush=True)
    cur = theta.copy()
    rng = np.random.default_rng(0)
    while True:
        active = [i for i in range(len(SEQ)) if mask[i]]
        cand = []
        for i in active:
            m = mask.copy(); m[i] = False
            d, x = train_masked(m, cur, rng, n_random=1)
            cand.append((d, i, x))
        cand.sort(key=lambda z: z[0])
        d, i, x = cand[0]
        if d < tol:
            mask[i] = False
            cur = x
            print(f"  removed slot {i:>2} {str(SEQ[i]):>7}  ->  {sum(mask)} CNOTs, "
                  f"delta = {d:.2e}", flush=True)
        else:
            print(f"  no further removal keeps delta<{tol:g} (best try = {d:.2e})", flush=True)
            break
    ncx = sum(mask)
    remaining = [SEQ[i] for i in range(len(SEQ)) if mask[i]]
    oe = _obs_err(ansatz_masked(mask), cur)
    print(f"\n  MINIMAL (greedy) CNOT count = {ncx}   (vs Step-4 14, textbook 16, full 18)")
    print(f"  remaining CNOTs: {remaining}")
    print(f"  final delta = {_delta(ansatz_masked(mask), cur):.2e},  observable err = {oe:.2e}")
    np.save("pqc_ring_pruned_params.npy", cur)
    json.dump({"ncx": ncx, "mask": mask, "remaining": remaining},
              open("pqc_ring_pruned.json", "w"), indent=2)
    print("  saved pqc_ring_pruned_params.npy, pqc_ring_pruned.json")


if __name__ == "__main__":
    main()
