"""Auxiliary (isometry): does the ISOMETRY need fewer CNOTs than the full unitary?

Start from the 14-CNOT solution that compiles the full U (which therefore compiles the
isometry exactly), and greedily prune further using the ISOMETRY cost. If it drops
below 14, the isometry genuinely needs fewer CNOTs (isometry < unitary); if it stops at
14, both share the floor. This avoids the greedy path-dependence of starting from a
separate 18-CNOT isometry solution.
"""
import json, numpy as np
from scipy.optimize import minimize
from pqc_ring_prune import ansatz_masked, SEQ, NPAR
from pqc_ring_ansatz import cost_grad
from pqc_common import U0_TARGET, _ANC0_ISO, U_TARGET, DIM

TGT, E = U0_TARGET, _ANC0_ISO


def iso_delta(ops, x):
    return cost_grad(ops, x, TGT, E)[0]


def train(mask, warm, seed, n_random=2):
    ops = ansatz_masked(mask)
    rng = np.random.default_rng(seed)
    starts = [warm] + [rng.uniform(-np.pi, np.pi, NPAR) for _ in range(n_random)]
    best = (np.inf, None)
    for x0 in starts:
        r = minimize(lambda x: cost_grad(ops, x, TGT, E), x0, jac=True, method="L-BFGS-B",
                     options=dict(maxiter=7000, ftol=1e-15, gtol=1e-13))
        if r.fun < best[0]:
            best = (float(r.fun), r.x)
    return best


def main(tol=1e-6):
    mask = json.load(open("pqc_ring_pruned.json"))["mask"]     # 14 active (full-U floor)
    theta = np.load("pqc_ring_pruned_params.npy")
    ops = ansatz_masked(mask)
    dfull = cost_grad(ops, theta, U_TARGET, np.eye(DIM, dtype=complex))[0]
    diso = iso_delta(ops, theta)
    print(f"14-CNOT (Step 5) solution: delta_full = {dfull:.2e}, delta_iso = {diso:.2e}")
    print(f"  (isometry compiled by the 14-CNOT full-U circuit as expected)\n", flush=True)

    cur = theta.copy()
    rng = np.random.default_rng(7)
    while True:
        active = [i for i in range(len(SEQ)) if mask[i]]
        cand = []
        for i in active:
            m = mask.copy(); m[i] = False
            d, x = train(m, cur, 800 + i, n_random=2)
            cand.append((d, i, x))
        cand.sort(key=lambda z: z[0])
        d, i, x = cand[0]
        if d < tol:
            mask[i] = False; cur = x
            print(f"  removed slot {i:>2} {str(SEQ[i]):>7} -> {sum(mask)} CNOTs, "
                  f"delta_iso = {d:.2e}", flush=True)
        else:
            print(f"  no further removal keeps delta_iso<{tol:g} (best try = {d:.2e})", flush=True)
            break
    print(f"\n  isometry floor (pruned from the 14-CNOT full-U solution) = {sum(mask)} CNOTs")
    print(f"  (Step 5 full-unitary floor = 14)")


if __name__ == "__main__":
    main()
