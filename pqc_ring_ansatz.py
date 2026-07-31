"""
Step 5 (success) — RX-RY-RZ PQC that DOES compile the gadget.
=============================================================

After the generic-ansatz attempt (pqc_compile.py) failed to compile U = H_a . CSWAP(0;1,3) .
CSWAP(0;2,4) . H_a with generic ansaetze, we searched what an RX-RY-RZ + CNOT ansatz
needs to succeed. Three knobs, tested with exact analytic gradients + a teacher-student
(in-class recompilation) test that separates EXPRESSIBILITY (U not reachable) from
TRAINABILITY (optimizer fails even on a provably-reachable target):

  1. connectivity  — a linear chain 0-1-2-3-4 (ansatz_ring) cannot express U even at
     20 CNOTs (in-class solves to 1e-15, U stays at delta~0.6): the ancilla must reach
     BOTH registers and the swap pairs (1,3),(2,4).  GADGET_PAIRS supplies exactly
     those links.
  2. rotation placement — a full RX,RY,RZ layer after EVERY CNOT (ansatz_percnot),
     not just between CNOT blocks, is what makes U representable at a practical depth.
  3. depth — enough CNOTs (>= ~14).

Result (ansatz_percnot(GADGET_PAIRS, L)):

  | L | CNOTs | params | delta(U) |
  |---|-------|--------|----------|
  | 1 |   6   |  105   | 0.86  (expressibility: too few CNOTs) |
  | 2 |  12   |  195   | 0.44  (expressibility) |
  | 3 |  18   |  285   | ~5e-15  --  EXACT COMPILATION |

So a generic RX-RY-RZ PQC DOES reproduce the gadget unitary (to machine precision) at
18 CNOTs with the right connectivity and per-CNOT rotations. (At L=3 the landscape is
already in the barren-plateau onset — random in-class targets start to fail — but the
structured target U is still found within a modest number of restarts.)

main() retrains the L=3 winner, reports the restart success rate, verifies exact
compilation and observable reproduction, and saves the parameters to
pqc_ring_L3_params.npy (reload with load_and_verify()).

Run:  python pqc_ring_ansatz.py
"""

import numpy as np
from scipy.optimize import minimize

from pqc_common import (
    DIM, U_TARGET, _X, _Y, _Z, _apply_1q_left, _apply_1q_right, _embed_1q,
    _cnot_perm, F_from_unitary, F_exact,
)

CHAIN = [(0, 1), (1, 2), (2, 3), (3, 4)]     # linear chain (ring minus 4->0): 4 CNOTs
_PAULI = {"rx": _X, "ry": _Y, "rz": _Z}
_CACHE = {}


# --- ansatz schedule -------------------------------------------------------
def ansatz_pairs(pairs, L=1):
    """(ops, n_params): L blocks of [RX,RY,RZ per qubit; CNOTs on `pairs`], then a
    final RX,RY,RZ layer.  n_params = 15*(L+1), CNOTs = len(pairs)*L."""
    ops, p = [], 0
    for _ in range(L):
        for w in range(5):
            for kind in ("rx", "ry", "rz"):
                ops.append(("g", kind, w, p)); p += 1
        for c, t in pairs:
            ops.append(("cnot", c, t))
    for w in range(5):                        # final rotation layer
        for kind in ("rx", "ry", "rz"):
            ops.append(("g", kind, w, p)); p += 1
    return ops, p


def ansatz_ring(L):
    """(ops, n_params) for L layers of the linear CNOT chain (ring minus 4->0)."""
    return ansatz_pairs(CHAIN, L)


# gadget-native connectivity: ancilla to both registers + the swap pairs (1,3),(2,4)
GADGET_PAIRS = [(0, 1), (0, 3), (1, 3), (0, 2), (0, 4), (2, 4)]


def ansatz_percnot(pairs, L=1):
    """Maximum-expressibility placement: an initial full RX,RY,RZ layer, then a full
    RX,RY,RZ layer after EVERY CNOT.  n_params = 15*(1 + len(pairs)*L),
    CNOTs = len(pairs)*L."""
    ops, p = [], 0
    for w in range(5):
        for kind in ("rx", "ry", "rz"):
            ops.append(("g", kind, w, p)); p += 1
    for _ in range(L):
        for c, t in pairs:
            ops.append(("cnot", c, t))
            for w in range(5):
                for kind in ("rx", "ry", "rz"):
                    ops.append(("g", kind, w, p)); p += 1
    return ops, p


# --- single-qubit gate matrices and derivatives ----------------------------
def _rmat(kind, a):
    c, s = np.cos(a / 2), np.sin(a / 2)
    if kind == "rx":
        return np.array([[c, -1j * s], [-1j * s, c]], dtype=complex)
    if kind == "ry":
        return np.array([[c, -s], [s, c]], dtype=complex)
    return np.array([[np.exp(-0.5j * a), 0], [0, np.exp(0.5j * a)]], dtype=complex)


def _drmat(kind, a):
    return (-0.5j) * _PAULI[kind] @ _rmat(kind, a)   # d/da exp(-i a P/2) = -i/2 P R


def _perm(c, t):
    if (c, t) not in _CACHE:
        _CACHE[(c, t)] = _cnot_perm(c, t)
    return _CACHE[(c, t)]


def unitary(ops, params, cols=None):
    V = np.eye(DIM, dtype=complex) if cols is None else np.asarray(cols, dtype=complex)
    for op in ops:
        if op[0] == "g":
            _, kind, w, p = op
            V = _apply_1q_left(V, _rmat(kind, params[p]), w)
        else:
            V = V[_perm(op[1], op[2]), :]
    return V


# --- exact (cost, grad) for the HS cost 1 - |Tr(target^dag V E)|^2/dc^2 -----
def cost_grad(ops, params, target, E):
    params = np.asarray(params, float)
    m = len(ops)
    dc = target.shape[1]
    P = [np.eye(DIM, dtype=complex)]
    for op in ops:
        if op[0] == "g":
            P.append(_apply_1q_left(P[-1], _rmat(op[1], params[op[3]]), op[2]))
        else:
            P.append(P[-1][_perm(op[1], op[2]), :])
    V = P[m]
    g = np.vdot(target, V @ E)
    cost = 1.0 - abs(g) ** 2 / dc ** 2
    Suf = [None] * (m + 1)
    Suf[m] = np.eye(DIM, dtype=complex)
    for i in range(m - 1, -1, -1):
        op = ops[i]
        if op[0] == "g":
            Suf[i] = _apply_1q_right(Suf[i + 1], _rmat(op[1], params[op[3]]), op[2])
        else:
            Suf[i] = Suf[i + 1][:, _perm(op[1], op[2])]
    Tdag = target.conj().T
    pref = -2.0 / dc ** 2
    grad = np.zeros_like(params)
    for i, op in enumerate(ops):
        if op[0] != "g":
            continue
        _, kind, w, p = op
        A = (P[i] @ E) @ Tdag @ Suf[i + 1]
        grad[p] = pref * np.real(np.conj(g) * np.trace(A @ _embed_1q(_drmat(kind, params[p]), w)))
    return cost, grad


_I32 = np.eye(DIM, dtype=complex)


def _train(target, ops, npar, restarts, seed):
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    for _ in range(restarts):
        x0 = rng.uniform(-np.pi, np.pi, npar)
        r = minimize(lambda x: cost_grad(ops, x, target, _I32), x0, jac=True,
                     method="L-BFGS-B", options=dict(maxiter=4000, ftol=1e-15, gtol=1e-13))
        if r.fun < best[0]:
            best = (float(r.fun), r.x)
    return best


EPS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]


def _obs_err(ops, x):
    V = unitary(ops, x)
    return max(abs(F_from_unitary(V, e) - F_exact(e)) for e in EPS)


def success_rate(ops, npar, restarts, seed0):
    """Over `restarts` independent single restarts, return (best_delta, best_params,
    n_success) with success := final delta(U) < 1e-6."""
    best, nsucc = (np.inf, None), 0
    for r in range(restarts):
        x0 = np.random.default_rng(seed0 + r).uniform(-np.pi, np.pi, npar)
        res = minimize(lambda x: cost_grad(ops, x, U_TARGET, _I32), x0, jac=True,
                       method="L-BFGS-B", options=dict(maxiter=6000, ftol=1e-15, gtol=1e-13))
        if res.fun < 1e-6:
            nsucc += 1
        if res.fun < best[0]:
            best = (float(res.fun), res.x)
    return best[0], best[1], nsucc


PARAMS_FILE = "pqc_ring_L3_params.npy"


def main(restarts=20):
    print("=" * 84)
    print(" Step 5 (success) — RX-RY-RZ after EVERY CNOT, gadget connectivity, L=3")
    print("=" * 84)
    ops, npar = ansatz_percnot(GADGET_PAIRS, 3)
    ncx = len(GADGET_PAIRS) * 3
    print(f"  pairs = {GADGET_PAIRS};  CNOTs = {ncx};  params = {npar}")
    print(f"  training U = H_a CSWAP CSWAP H_a  ({restarts} restarts)...\n", flush=True)

    d, x, ns = success_rate(ops, npar, restarts, seed0=1000)
    print(f"  best delta(U) = 1 - fidelity = {d:.2e}")
    print(f"  robustness: delta(U) < 1e-6 in {ns}/{restarts} restarts ({100*ns/restarts:.0f}%)")

    V = unitary(ops, x)
    oe = _obs_err(ops, x)
    print(f"  |Tr(U^dag V)|/32 = {abs(np.vdot(U_TARGET, V))/DIM:.12f}   (=1 => exact)")
    print(f"  observable error  max|F - F_exact| over eps = {oe:.2e}")

    np.save(PARAMS_FILE, x)
    print(f"\n  saved best params -> {PARAMS_FILE}")
    print(f"  reconstruct with:  ops,_ = ansatz_percnot(GADGET_PAIRS, 3);  V = unitary(ops, params)")


def load_and_verify(path=PARAMS_FILE):
    """Reload saved params and re-check exact compilation."""
    x = np.load(path)
    ops, _ = ansatz_percnot(GADGET_PAIRS, 3)
    V = unitary(ops, x)
    return abs(np.vdot(U_TARGET, V)) / DIM, _obs_err(ops, x)


if __name__ == "__main__":
    main()
