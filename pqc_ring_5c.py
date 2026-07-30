"""
Step 5c on the 14-CNOT circuit: prune with the OBSERVABLE cost (ancilla-parity).
===============================================================================

Take the exact 14-CNOT circuit and prune it under the weaker 5c requirement -- only the
purified observable via the ancilla-parity read-out F = <Z_a (x) O>/<Z_a>. We match the
correlators <Z_a> -> Tr(rho^2) and <Z_a (x) O> -> Tr(O rho^2) (anchored, so no
denominator degeneracy) for O in {Phi+, ZZ} over an eps grid. Greedy pruning then finds
how few CNOTs reproduce the observable -- the 5c rung of the relaxation ladder on the
same ansatz+prune footing as 5a/5b (both 14).

Exact analytic gradient of the expectation cost (verified vs finite differences).

Run:  python pqc_ring_5c.py
"""
import json, numpy as np
from scipy.optimize import minimize
from pqc_ring_ansatz import _rmat, _drmat, _perm
from pqc_ring_prune import ansatz_masked, SEQ, NPAR
from pqc_common import (DIM, _apply_1q_left, _apply_1q_right, _embed_1q,
                        _obs_ZO, _OBS_Z, _rho_in, O_PHI_PLUS)
from noisy_bell_state import rho_eps_analytic

_Z = np.array([[1, 0], [0, -1]], dtype=complex)
O_ZZ = np.kron(_Z, _Z)
EPS_T = [0.2, 0.4, 0.6]


def _tr(O, r2):
    return float(np.real(np.trace(O @ r2)))


def build_terms():
    """terms: list of (M 32x32, rho_in 32x32, target float)."""
    terms = []
    for e in EPS_T:
        r = rho_eps_analytic(e); r2 = r @ r
        rin = _rho_in(e)
        terms.append((_OBS_Z, rin, _tr(np.eye(4), r2)))                  # <Z_a>=Tr(rho^2)
        terms.append((_obs_ZO(O_PHI_PLUS), rin, _tr(O_PHI_PLUS, r2)))    # <Z_a O_Phi>
        terms.append((_obs_ZO(O_ZZ), rin, _tr(O_ZZ, r2)))               # <Z_a O_ZZ>
    return terms


# left/right application for ring ops
def apply_left(op, M, params):
    if op[0] == "g":
        return _apply_1q_left(M, _rmat(op[1], params[op[3]]), op[2])
    return M[_perm(op[1], op[2]), :]


def apply_right(op, M, params):
    if op[0] == "g":
        return _apply_1q_right(M, _rmat(op[1], params[op[3]]), op[2])
    return M[:, _perm(op[1], op[2])]


def obs_cost_grad(ops, params, terms):
    params = np.asarray(params, float)
    m = len(ops)
    P = [np.eye(DIM, dtype=complex)]
    for op in ops:
        P.append(apply_left(op, P[-1], params))
    V = P[m]; Vd = V.conj().T
    exps = [float(np.real(np.trace(M @ V @ rin @ Vd))) for (M, rin, _) in terms]
    cost = sum((e - tau) ** 2 for e, (_, _, tau) in zip(exps, terms))
    Suf = [None] * (m + 1); Suf[m] = np.eye(DIM, dtype=complex)
    for i in range(m - 1, -1, -1):
        Suf[i] = apply_right(ops[i], Suf[i + 1], params)
    Bracket = np.zeros((DIM, DIM), dtype=complex)
    for e, (M, rin, tau) in zip(exps, terms):
        Bracket += 2 * (e - tau) * (rin @ Vd @ M)
    grad = np.zeros_like(params)
    for i, op in enumerate(ops):
        if op[0] != "g":
            continue
        A = P[i] @ Bracket @ Suf[i + 1]
        grad[op[3]] = 2 * np.real(np.trace(A @ _embed_1q(_drmat(op[1], params[op[3]]), op[2])))
    return cost, grad


def train(mask, terms, warm, seed, n_random=1):
    ops = ansatz_masked(mask)
    rng = np.random.default_rng(seed)
    starts = [warm] + [rng.uniform(-np.pi, np.pi, NPAR) for _ in range(n_random)]
    best = (np.inf, None)
    for x0 in starts:
        r = minimize(lambda x: obs_cost_grad(ops, x, terms), x0, jac=True, method="L-BFGS-B",
                     options=dict(maxiter=1200, ftol=1e-14, gtol=1e-12))
        if r.fun < best[0]:
            best = (float(r.fun), r.x)
    return best


def main(tol=1e-6):
    print("=" * 78)
    print(" Step 5c on the 14-CNOT circuit: prune with the OBSERVABLE cost")
    print("=" * 78)
    terms = build_terms()
    mask = json.load(open("pqc_ring_pruned.json"))["mask"]
    theta = np.load("pqc_ring_pruned_params.npy")
    c0 = obs_cost_grad(ansatz_masked(mask), theta, terms)[0]
    print(f"\n start: {sum(mask)} CNOTs, obs cost = {c0:.2e}\n", flush=True)
    cur = theta.copy()
    while True:
        active = [i for i in range(len(SEQ)) if mask[i]]
        cand = []
        for i in active:
            m = mask.copy(); m[i] = False
            d, x = train(m, terms, cur, 800 + i, n_random=0)
            cand.append((d, i, x))
        cand.sort(key=lambda z: z[0])
        d, i, x = cand[0]
        if d < tol:
            mask[i] = False; cur = x
            print(f"  removed slot {i:>2} {str(SEQ[i]):>7} -> {sum(mask)} CNOTs, "
                  f"obs cost = {d:.2e}", flush=True)
        else:
            print(f"  no further removal keeps obs cost<{tol:g} (best try = {d:.2e})", flush=True)
            break
    ncx = sum(mask)
    print(f"\n  5c observable floor (ancilla-parity read-out) = {ncx} CNOTs")
    print(f"  (5a full-unitary = 14, 5b isometry = 14; Step 4b destructive = 2)")
    print(f"  remaining CNOTs: {[SEQ[i] for i in range(len(SEQ)) if mask[i]]}")
    np.save("pqc_ring_5c_params.npy", cur)
    json.dump({"ncx": ncx, "mask": mask}, open("pqc_ring_5c.json", "w"))


if __name__ == "__main__":
    main()
