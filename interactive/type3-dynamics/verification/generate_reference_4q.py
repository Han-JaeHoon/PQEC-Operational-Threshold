"""
Independent reference trajectories for the 4Q mode of the interactive simulator.
================================================================================

The simulator iterates the closed-form 4Q map on its invariant plane (a, b).  This script
produces trajectories that use NONE of that closed form: the full two-qubit density matrix
is iterated through the actual noisy four-qubit 5-CNOT circuit (+ H on q3, postselect
(m3, m4) = (0, 0)), by two routes that share nothing beyond the circuit's gate list.

  R1  the Bell-state-distillation package's dense 16x16 simulator:
      pqec_distill.noise.noisy_full_channel_dm(rho (x) rho, q, "replace") followed by
      pqec_distill.measurement.postselect_branch(., 0, 0), normalised.
  R2  PennyLane default.mixed, written here from the gate list only
      (pqec_distill.circuit.CNOT_SEQUENCE, Q3): qml.CNOT, the replacement channel as its own
      16 Pauli Kraus operators, qml.Hadamard.  The next state is rebuilt from all 16
      expectations <P (x) |00><00|_{q3 q4}>, i.e. without assuming Bell-diagonality.

Both routes record the per-round success probability P_succ = Tr(rho_tilde).  Unlike the
Type-3/4 parity visibility, P_succ stays >= 1/4 for every q (it is 1/4 at q = 1), so q = 1
is included here.

Nothing in either repository is modified.  The output JSON is consumed by
tests/reference4q.test.ts.

Run (from anywhere):
    python interactive/type3-dynamics/verification/generate_reference_4q.py \
        [--bsd PATH/TO/Bell-state-distillation]
Needs numpy + pennylane and a checkout of Bell-state-distillation (default: next to this
repository).
"""

import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "reference_4q.json")

# (eps, q, rounds).  4Q thresholds: q_abs(0.15) = 0.06066542, q_rel(0.15) = 0.06266623,
# q_abs(0.10) = 0.04329794, q_rel(0.10) = 0.04428694, q_rel(0.75) = 0.15865660,
# q_rel(0.90) = 0.15406397 (calculation notes, Tables).
CASES = [
    (0.15, 0.0606, 40),   # just below q_abs(0.15)
    (0.15, 0.0628, 40),   # just above q_rel(0.15)
    (0.10, 0.0438, 40),   # between q_abs(0.10) and q_rel(0.10)
    (0.10, 0.0615, 40),   # F decreases from round 0 (just above the F transition near 0.0612)
    (0.10, 0.0100, 40),   # small q
    (0.02, 0.0010, 30),   # very small eps and q
    (0.30, 0.1000, 40),   # moderate
    (0.50, 0.1500, 40),
    (0.60, 0.1300, 40),
    (0.75, 0.1500, 40),   # separable input, below q_rel(0.75): q_rel is defined there for 4Q
    (0.90, 0.1600, 40),   # separable input, above q_rel(0.90)
    (0.05, 0.1500, 40),   # strong noise on a good input
    (0.40, 0.3000, 40),   # high noise: decays towards I/4
    (0.00, 0.0200, 40),   # pure Bell input
    (0.20, 0.0000, 20),   # ideal gadget
    (0.30, 0.6000, 10),   # very strong noise
    (0.30, 1.0000, 5),    # q = 1: one round sends any input to I/4, P_succ = 1/4
    (1.00, 0.0500, 10),   # maximally mixed input
]

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
P1 = {"I": I2, "X": X, "Y": Y, "Z": Z}
LABELS = [a + b for a in "IXYZ" for b in "IXYZ"]
P2 = {l: np.kron(P1[l[0]], P1[l[1]]) for l in LABELS}
BELL = np.array([[1, 0, 0, 1], [1, 0, 0, -1], [0, 1, 1, 0], [0, 1, -1, 0]], dtype=complex) / np.sqrt(2)


def bell_pops(rho):
    R = BELL.conj() @ rho @ BELL.T
    pops = np.real(np.diag(R))
    off = float(np.linalg.norm(R - np.diag(np.diag(R))))
    return pops.tolist(), off


def rho_isotropic(eps):
    phi = BELL[0].reshape(4, 1)
    return (1 - eps) * (phi @ phi.conj().T) + eps * np.eye(4) / 4


def git(path, *args):
    return subprocess.check_output(["git", "-C", path, *args], text=True).strip()


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def make_r1():
    from pqec_distill.measurement import postselect_branch
    from pqec_distill.noise import noisy_full_channel_dm

    def step(rho, q):
        tau, prob = postselect_branch(noisy_full_channel_dm(np.kron(rho, rho), q, "replace"), 0, 0)
        return tau / prob, float(prob)

    return step


def make_r2(qml, cnots, q3):
    proj00 = np.zeros((4, 4), dtype=complex)
    proj00[0, 0] = 1.0
    obs = [np.kron(P2[l], proj00) for l in LABELS]
    dev = qml.device("default.mixed", wires=4)

    def kraus(q):
        # D_q(s) = (1-q) s + q I/4 (x) Tr s  =  sum_P K_P s K_P^dag over the 16 two-qubit Paulis
        return [np.sqrt(1 - q + q / 16) * P2["II"]] + [np.sqrt(q / 16) * P2[l] for l in LABELS[1:]]

    @qml.qnode(dev)
    def node(rho_in, q):
        qml.QubitDensityMatrix(rho_in, wires=[0, 1, 2, 3])
        for c, t in cnots:
            qml.CNOT(wires=[c, t])
            qml.QubitChannel(kraus(q), wires=[c, t])
        qml.Hadamard(q3)
        return [qml.expval(qml.Hermitian(o, wires=[0, 1, 2, 3])) for o in obs]

    def step(rho, q):
        c = np.array(node(np.kron(rho, rho), q), dtype=float)
        tau = sum(ci * P2[l] for ci, l in zip(c, LABELS)) / 4.0
        prob = float(np.real(np.trace(tau)))
        return tau / prob, prob

    return step


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bsd", default=os.path.join(REPO, "..", "Bell-state-distillation"))
    args = ap.parse_args()
    bsd = os.path.abspath(args.bsd)
    sys.path.insert(0, os.path.join(bsd, "src"))

    import pennylane as qml
    from pqec_distill.circuit import CNOT_SEQUENCE, Q3

    n_cx = len(CNOT_SEQUENCE)
    assert n_cx == 5, n_cx
    r1_step = make_r1()
    r2_step = make_r2(qml, CNOT_SEQUENCE, Q3)

    cases = []
    for eps, q, N in CASES:
        out = {}
        for name, stepf in (("R1", r1_step), ("R2", r2_step)):
            rho = rho_isotropic(eps).astype(complex)
            pops, psucc, offm, herm, trerr, mineig = [], [], 0.0, 0.0, 0.0, 1.0
            for n in range(N + 1):
                p, off = bell_pops(rho)
                pops.append(p)
                offm = max(offm, off)
                trerr = max(trerr, abs(float(np.real(np.trace(rho))) - 1.0))
                mineig = min(mineig, float(np.min(np.linalg.eigvalsh(0.5 * (rho + rho.conj().T)))))
                if n < N:
                    nxt, prob = stepf(rho, q)
                    psucc.append(prob)
                    herm = max(herm, float(np.linalg.norm(nxt - nxt.conj().T)))
                    rho = 0.5 * (nxt + nxt.conj().T)
            out[name] = dict(pops=pops, psucc=psucc, off=offm, herm=herm, tr=trerr, mineig=mineig)
        r1, r2 = out["R1"], out["R2"]
        d12 = float(np.max(np.abs(np.array(r1["pops"]) - np.array(r2["pops"]))))
        dp = float(np.max(np.abs(np.array(r1["psucc"]) - np.array(r2["psucc"]))))
        cases.append(dict(
            eps=eps, q=q, N=N,
            R1_dense=r1["pops"], R2_pennylane=r2["pops"],
            R1_psucc=r1["psucc"], R2_psucc=r2["psucc"],
            max_abs_R1_minus_R2=d12, max_abs_psucc_R1_minus_R2=dp,
            max_bell_offdiag_R1=r1["off"], max_bell_offdiag_R2=r2["off"],
            max_herm_residual_R2=r2["herm"], max_trace_err_R2=r2["tr"],
            min_eig_R1=r1["mineig"], min_eig_R2=r2["mineig"],
        ))
        print(f"eps={eps:.3f} q={q:.4f} N={N:3d}  |R1-R2|={d12:.2e}  |dPsucc|={dp:.1e}  "
              f"offdiag R1={r1['off']:.1e} R2={r2['off']:.1e}  min eig={min(r1['mineig'], r2['mineig']):+.1e}  "
              f"F1={r1['pops'][1][0]:.15f}")

    meta = dict(
        generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        repo_head=git(REPO, "rev-parse", "HEAD"),
        bsd_commit=git(bsd, "rev-parse", "HEAD"),
        circuit="Bell-state-distillation pqec_distill.circuit.CNOT_SEQUENCE + H(q3), "
                "postselect (m3, m4) = (0, 0)",
        circuit_file_sha256=sha256(os.path.join(bsd, "src", "pqec_distill", "circuit.py")),
        noise_file_sha256=sha256(os.path.join(bsd, "src", "pqec_distill", "noise.py")),
        n_cnot=n_cx,
        pennylane=qml.__version__, numpy=np.__version__,
        bell_order=["Phi+", "Phi-", "Psi+", "Psi-"],
        note="R1/R2 iterate the full 4x4 state through the noisy 5-CNOT circuit; "
             "neither uses the (a, b) closed form.",
    )
    with open(OUT, "w") as fh:
        json.dump(dict(meta=meta, cases=cases), fh, indent=1)
    print("wrote", os.path.relpath(OUT, REPO))


if __name__ == "__main__":
    main()
