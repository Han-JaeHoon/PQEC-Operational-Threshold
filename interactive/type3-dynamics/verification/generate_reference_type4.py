"""
Independent reference trajectories for the Type-4 mode of the interactive simulator.
====================================================================================

The simulator iterates the closed-form Type-4 map on the invariant plane (u, v).  This
script produces trajectories that use NONE of that closed form: the full two-qubit
density matrix is iterated through the actual noisy resynthesized 14-CNOT circuit
(pqec_resynth_noise.GATES), by two routes that share no code with each other beyond the
repository's gate list.  It mirrors generate_reference.py (Type 3), which is left as is.

  R1  dense 32x32 propagation -- iterated_noisy_pqec.iterate_effective_map("step4"),
      read from the `iterated-noisy-pqec` branch (circuit = pqec_resynth_noise.GATES,
      replacement channel implemented from its definition).
  R2  PennyLane default.mixed -- the repository's own Step-4 executor
      (pqec_resynth_noise._apply, global_depol_kraus after every CNOT).  The next state is
      rebuilt from all 16 parity-weighted Pauli correlators
      rho' = sum_P <Z_a (x) P> P / (4 <Z_a>), i.e. without assuming Bell-diagonality.
      The normalisation <Z_a> (the parity visibility) is recorded every round.
  R3  the repository's own one-round read-out pqec_resynth_noise.F_resynth, for F_1 only.

Nothing in the repository is modified: the branch file is extracted to a temporary
directory with `git show`.  The output JSON is consumed by tests/reference4.test.ts.

Run (from the repository root or anywhere):
    python interactive/type3-dynamics/verification/generate_reference_type4.py
Needs numpy + pennylane (requirements.txt of the repository) and the branch ref
`origin/iterated-noisy-pqec` (override with --ref).
"""

import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "reference_type4.json")

# (eps, q, rounds)
CASES = [
    (0.10, 0.0405, 40),   # rise, then decrease at n = 2 (Type-4 analogue of the Type-3 case)
    (0.15, 0.0590, 40),   # rise, then late decrease at n = 4
    (0.15, 0.0478, 40),   # just below q_abs(0.15) = 0.04781013
    (0.15, 0.0490, 40),   # just above q_rel(0.15) = 0.04898169
    (0.10, 0.0100, 40),   # small q
    (0.02, 0.0010, 30),   # very small eps and q
    (0.30, 0.0800, 40),   # moderate
    (0.50, 0.1300, 40),
    (0.60, 0.1500, 40),
    (0.05, 0.1500, 40),   # strong noise on a good input
    (0.40, 0.3000, 40),   # high noise: decays towards I/4
    (0.90, 0.0200, 40),   # separable input
    (0.00, 0.0200, 40),   # pure Bell input
    (0.20, 0.0000, 20),   # ideal gadget
    # Not q -> 1: the circuit's normalisation <Z_a> vanishes there (it underflows float64
    # for q >~ 0.97 and is exactly 0 at q = 1), so the full-circuit estimator is undefined.
    (0.30, 0.6000, 10),   # very strong noise, still well-conditioned
    (1.00, 0.0500, 10),   # maximally mixed input
]


def git(*args):
    return subprocess.check_output(["git", "-C", REPO, *args], text=True).strip()


def load_iterated_module(ref):
    src = subprocess.check_output(["git", "-C", REPO, "show", f"{ref}:iterated_noisy_pqec.py"])
    tmp = tempfile.mkdtemp(prefix="type4ref_")
    with open(os.path.join(tmp, "iterated_noisy_pqec.py"), "wb") as fh:
        fh.write(src)
    sys.path.insert(0, tmp)
    import iterated_noisy_pqec as inp  # noqa: E402
    return inp


def make_r2(qml):
    """PennyLane route: parity-weighted Pauli correlators of the repo's Step-4 circuit."""
    import pqec_resynth_noise as r4

    I2 = np.eye(2)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.diag([1.0, -1.0]).astype(complex)
    P1 = {"I": I2, "X": X, "Y": Y, "Z": Z}
    labels = [a + b for a in "IXYZ" for b in "IXYZ"]
    P2 = {l: np.kron(P1[l[0]], P1[l[1]]) for l in labels}

    dev = qml.device("default.mixed", wires=5)

    @qml.qnode(dev)
    def node(rho_AB, q):
        qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
        r4._apply(q)
        return [qml.expval(qml.PauliZ(0) @ qml.Hermitian(P2[l], wires=[1, 2])) for l in labels]

    def step(rho, q):
        c = np.array(node(np.kron(rho, rho), q), dtype=float)
        tau = sum(ci * P2[l] for ci, l in zip(c, labels)) / 4.0
        return tau / c[0], float(c[0])

    return step


BELL = np.array([[1, 0, 0, 1], [1, 0, 0, -1], [0, 1, 1, 0], [0, 1, -1, 0]], dtype=complex) / np.sqrt(2)


def bell_pops(rho):
    R = BELL.conj() @ rho @ BELL.T
    pops = np.real(np.diag(R))
    off = float(np.linalg.norm(R - np.diag(np.diag(R))))
    return pops.tolist(), off


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="origin/iterated-noisy-pqec")
    args = ap.parse_args()

    sys.path.insert(0, REPO)
    import pennylane as qml
    import pqec_resynth_noise as r4

    n_cx = sum(1 for g in r4.GATES if g[0] == "cx")
    assert n_cx == 14, n_cx
    inp = load_iterated_module(args.ref)
    r2_step = make_r2(qml)

    cases = []
    for eps, q, N in CASES:
        t = 1.0 - eps
        rho0 = inp.rho_isotropic(t)

        # R1: dense propagation, no early stop (tol = -1)
        res = inp.iterate_effective_map(rho0, q, "step4", tol=-1.0, max_iter=N)
        hist = sorted(res["history"], key=lambda h: h["n"])
        r1 = [[h["p_Phi_plus"], h["p_Phi_minus"], h["p_Psi_plus"], h["p_Psi_minus"]] for h in hist]
        r1_off = max(h["bell_offdiag"] for h in hist)
        assert [h["n"] for h in hist] == list(range(N + 1)), "R1 history incomplete"

        # R2: PennyLane full-state iteration
        rho = rho0.astype(complex)
        r2, vis, r2_off, r2_herm, r2_tr, r2_mineig = [], [], 0.0, 0.0, 0.0, 1.0
        for n in range(N + 1):
            pops, off = bell_pops(rho)
            r2.append(pops)
            r2_off = max(r2_off, off)
            r2_tr = max(r2_tr, abs(float(np.real(np.trace(rho))) - 1.0))
            r2_mineig = min(r2_mineig, float(np.min(np.linalg.eigvalsh(0.5 * (rho + rho.conj().T)))))
            if n < N:
                nxt, z = r2_step(rho, q)
                vis.append(z)
                r2_herm = max(r2_herm, float(np.linalg.norm(nxt - nxt.conj().T)))
                rho = 0.5 * (nxt + nxt.conj().T)

        f1_r3 = float(r4.F_resynth(eps, q))
        d12 = float(np.max(np.abs(np.array(r1) - np.array(r2))))
        cases.append(dict(
            eps=eps, q=q, N=N,
            R1_dense=r1, R2_pennylane=r2,
            R2_parity_visibility=vis,
            R3_F1_resynth=f1_r3,
            max_abs_R1_minus_R2=d12,
            max_bell_offdiag_R1=r1_off, max_bell_offdiag_R2=r2_off,
            max_herm_residual_R2=r2_herm, max_trace_err_R2=r2_tr, min_eig_R2=r2_mineig,
        ))
        print(f"eps={eps:.3f} q={q:.4f} N={N:3d}  |R1-R2|={d12:.2e}  "
              f"offdiag R1={r1_off:.1e} R2={r2_off:.1e}  min eig={r2_mineig:+.1e}  "
              f"F1: R1={r1[1][0]:.15f} R3={f1_r3:.15f}")

    with open(os.path.join(REPO, "pqec_resynth_noise.py"), "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    meta = dict(
        generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        repo_head=git("rev-parse", "HEAD"),
        iterated_branch_ref=args.ref,
        iterated_branch_commit=git("rev-parse", args.ref),
        circuit="pqec_resynth_noise.GATES (resynthesized 14-CNOT)",
        circuit_file_sha256=sha,
        n_cnot=n_cx,
        pennylane=qml.__version__, numpy=np.__version__,
        bell_order=["Phi+", "Phi-", "Psi+", "Psi-"],
        note="R1/R2 iterate the full 4x4 state through the noisy 14-CNOT circuit; "
             "neither uses the (u, v) closed form.",
    )
    with open(OUT, "w") as fh:
        json.dump(dict(meta=meta, cases=cases), fh, indent=1)
    print("wrote", os.path.relpath(OUT, REPO))


if __name__ == "__main__":
    main()
