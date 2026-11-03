"""
Independent reference trajectories for the BBPSSW mode of the interactive simulator.
=====================================================================================

The simulator iterates the scalar recurrence F_{n+1} = G_q(F_n).  This script produces data that
uses NONE of that closed form: the two-qubit state is iterated through the noisy four-qubit circuit
(CNOT A1->A2, D_q on A1 A2, CNOT B1->B2, D_q on B1 B2, keep 00/11 on A2 B2, twirl), by two routes
that share nothing beyond the circuit's definition.

  R1  verification/bbpssw_stage1/frozen_circuit_reference.py: bit-permuted CNOT matrices, D_q as an
      average over 16 Pauli conjugations, twirl as 24 explicit Clifford U (x) U* conjugations.
  R2  written here: CNOT by index relabeling, D_q literally as (1-q) rho + q (I/4) (x) Tr rho on the
      gate's two wires, twirl as the projection onto (F, (1-F)/3, (1-F)/3, (1-F)/3).

Also recorded, through R1 and R2: arbitrary Bell-diagonal inputs with and without the input twirl,
and through R1: the three Stage 1 negative controls.  Output: reference_bbpssw.json, consumed by
tests/referencebbpssw.test.ts.  Nothing else is modified.

Run (from anywhere; needs numpy):
    python interactive/type3-dynamics/verification/generate_reference_bbpssw.py
"""

import datetime
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STAGE1 = os.path.join(HERE, "bbpssw_stage1")
sys.path.insert(0, STAGE1)
sys.dont_write_bytecode = True
from frozen_circuit_reference import independent_circuit_round  # noqa: E402

PHI = [np.array(v, dtype=complex) / np.sqrt(2)
       for v in ([1, 0, 0, 1], [1, 0, 0, -1], [0, 1, 1, 0], [0, 1, -1, 0])]
PROJ = [np.outer(v, v.conj()) for v in PHI]


def r2_cnot(rho, c, t):
    def perm(i):
        b = [(i >> (3 - k)) & 1 for k in range(4)]
        b[t] ^= b[c]
        return sum(x << (3 - k) for k, x in enumerate(b))
    idx = [perm(i) for i in range(16)]
    out = np.zeros_like(rho)
    for i in range(16):
        for j in range(16):
            out[idx[i], idx[j]] = rho[i, j]
    return out


def r2_replace(rho, i, j, q):
    t = np.moveaxis(rho.reshape((2,) * 8), [i, j, i + 4, j + 4], [0, 1, 2, 3])
    tr = sum(t[a, b, a, b] for a in (0, 1) for b in (0, 1))
    new = np.zeros_like(t)
    for a in (0, 1):
        for b in (0, 1):
            new[a, b, a, b] = tr / 4
    new = np.moveaxis(new, [0, 1, 2, 3], [i, j, i + 4, j + 4]).reshape(16, 16)
    return (1 - q) * rho + q * new


def r2_isotropic(F):
    return [F, (1 - F) / 3, (1 - F) / 3, (1 - F) / 3]


def r2_round(p, q, pre=True, post=True):
    if pre:
        p = r2_isotropic(p[0])
    rho = sum(x * P for x, P in zip(p, PROJ))
    s = np.kron(rho, rho)
    s = r2_replace(r2_cnot(s, 0, 2), 0, 2, q)
    s = r2_replace(r2_cnot(s, 1, 3), 1, 3, q)
    m = s.reshape((2,) * 8)
    out = sum(m[:, :, a, a, :, :, a, a].reshape(4, 4) for a in (0, 1))
    P = float(np.trace(out).real)
    pops = [float(np.trace(X @ out).real) / P for X in PROJ]
    return (r2_isotropic(pops[0]) if post else pops), P


def r1_round(p, q, **kw):
    r = independent_circuit_round(p, q, **kw)
    return [float(x) for x in r["populations"]], float(r["success"])


def trajectory(route, eps, q, N):
    F0 = 1 - 3 * eps / 4
    p = r2_isotropic(F0)
    pops, succ = [p], []
    for _ in range(N):
        p, P = route(p, q)
        pops.append(p)
        succ.append(P)
    return pops, succ


def sha(name):
    with open(os.path.join(STAGE1, name), "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    N = 15
    cases = []
    for eps in (0, 0.001, 0.1, 0.3, 0.5, 2 / 3, 0.8, 0.95, 1):
        for q in (0, 1e-6, 0.001, 0.0327, 0.05, 0.08, 0.13, 0.18, 0.3, 0.6, 0.9, 1):
            r1 = trajectory(lambda p, x: r1_round(p, x, twirl=True, pre_twirl=True), eps, q, N)
            r2 = trajectory(lambda p, x: r2_round(p, x), eps, q, N)
            gap = max(np.abs(np.array(r1[0]) - np.array(r2[0])).max(),
                      np.abs(np.array(r1[1]) - np.array(r2[1])).max())
            cases.append({"eps": eps, "q": q, "N": N, "R1_circuit": r1[0], "R2_circuit": r2[0],
                          "R1_psucc": r1[1], "R2_psucc": r2[1], "max_abs_R1_minus_R2": float(gap)})

    rng = np.random.default_rng(20261103)
    arbitrary = []
    for k in range(30):
        p = [float(x) for x in rng.dirichlet(np.ones(4) * [0.3, 1.0, 5.0][k % 3])]
        for q in (0, 0.01, 0.08, 0.3, 1):
            arbitrary.append({
                "p": p, "q": q,
                "original_R1": r1_round(p, q, twirl=True, pre_twirl=True),
                "original_R2": r2_round(p, q),
                "no_input_twirl_R1": r1_round(p, q, twirl=True, pre_twirl=False),
                "no_input_twirl_R2": r2_round(p, q, pre=False),
                "no_twirls_R1": r1_round(p, q, twirl=False, pre_twirl=False),
                "no_twirls_R2": r2_round(p, q, pre=False, post=False),
            })

    p, q = (0.67, 0.11, 0.15, 0.07), 0.17
    controls = {
        "p": p, "q": q,
        "reverse_Bob_CNOT": r1_round(p, q, reverse_bob=True),
        "omit_Bob_depolarization": r1_round(p, q, bob_noise_override=0.0),
        "accept_unequal_measurements": r1_round(p, q, accepted_outcomes=((0, 1), (1, 0))),
    }

    meta = {
        "n_cnot": 2,
        "date": datetime.date.today().isoformat(),
        "numpy": np.__version__,
        "frozen_circuit_reference_sha256": sha("frozen_circuit_reference.py"),
        "frozen_circuit_spec_sha256": sha("frozen_circuit_spec.json"),
    }
    out = os.path.join(HERE, "reference_bbpssw.json")
    with open(out, "w") as fh:
        json.dump({"meta": meta, "cases": cases, "arbitrary": arbitrary, "negative_controls": controls}, fh)
    print(out, os.path.getsize(out), "bytes;", len(cases), "trajectories,", len(arbitrary), "arbitrary inputs;",
          "max |R1-R2| =", max(c["max_abs_R1_minus_R2"] for c in cases))


if __name__ == "__main__":
    main()
