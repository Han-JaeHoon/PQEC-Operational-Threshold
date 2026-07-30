"""
Step 5b with the successful 5a ansatz: compile the ancilla-|0> ISOMETRY, then prune.
===================================================================================

Reuse the gadget-matched, per-CNOT-rotation ansatz that solved 5a, but target the
ancilla-|0> isometry U0 = U E0 (32x16) instead of the full unitary. Since 5b <= 5a in
difficulty, we ask: does it compile at fewer CNOTs, and does greedy pruning reach a
lower CNOT floor than 5a's 14?

Run:  python pqc_ring_5b.py
"""
import numpy as np
from scipy.optimize import minimize
from pqc_ring_ansatz import ansatz_percnot, GADGET_PAIRS, cost_grad, unitary
from pqc_ring_prune import ansatz_masked, SEQ, NPAR
from pqc_common import U0_TARGET, _ANC0_ISO

TGT, E = U0_TARGET, _ANC0_ISO          # 5b isometry target (32x16), dc=16


def train(ops, npar, restarts, seed, warm=None):
    rng = np.random.default_rng(seed)
    starts = ([warm] if warm is not None else []) + \
             [rng.uniform(-np.pi, np.pi, npar) for _ in range(restarts)]
    best = (np.inf, None)
    for x0 in starts:
        r = minimize(lambda x: cost_grad(ops, x, TGT, E), x0, jac=True, method="L-BFGS-B",
                     options=dict(maxiter=6000, ftol=1e-15, gtol=1e-13))
        if r.fun < best[0]:
            best = (float(r.fun), r.x)
    return best


def _delta(ops, x):
    return cost_grad(ops, x, TGT, E)[0]


def main(tol=1e-6):
    print("=" * 78)
    print(" Step 5b (isometry) with the 5a ansatz: compile + prune")
    print("=" * 78)
    print(f"\n compile scan (ansatz_percnot, GADGET_PAIRS):")
    print(f"   {'L':>2} {'CNOTs':>6} {'params':>7} {'delta_iso':>11}")
    for L in [1, 2, 3]:
        ops, npar = ansatz_percnot(GADGET_PAIRS, L)
        d, _ = train(ops, npar, 20, 10 + L)
        print(f"   {L:>2} {6*L:>6} {npar:>7} {d:>11.2e}", flush=True)

    # prune from an 18-CNOT isometry solution (L=3)
    ops3, npar3 = ansatz_percnot(GADGET_PAIRS, 3)
    d3, theta = train(ops3, npar3, 24, 33)
    print(f"\n prune start: 18 CNOTs, delta_iso = {d3:.2e}")
    mask = [True] * len(SEQ)
    cur = theta.copy()
    rng = np.random.default_rng(0)
    while True:
        active = [i for i in range(len(SEQ)) if mask[i]]
        cand = []
        for i in active:
            m = mask.copy(); m[i] = False
            ops = ansatz_masked(m)
            d, x = train(ops, NPAR, 1, 800 + i, warm=cur)
            cand.append((d, i, x))
        cand.sort(key=lambda z: z[0])
        d, i, x = cand[0]
        if d < tol:
            mask[i] = False; cur = x
            print(f"   removed slot {i:>2} {str(SEQ[i]):>7} -> {sum(mask)} CNOTs, "
                  f"delta_iso = {d:.2e}", flush=True)
        else:
            print(f"   no further removal keeps delta<{tol:g} (best try = {d:.2e})", flush=True)
            break
    ncx = sum(mask)
    print(f"\n  5b (isometry) minimal (greedy) CNOT count = {ncx}   (vs 5a full-unitary 14)")
    print(f"  remaining CNOTs: {[SEQ[i] for i in range(len(SEQ)) if mask[i]]}")
    print(f"  final delta_iso = {_delta(ansatz_masked(mask), cur):.2e}")
    np.save("pqc_ring_5b_params.npy", cur)


if __name__ == "__main__":
    main()
