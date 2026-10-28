"""
Independent reference trajectories for the Type-3 interactive simulator.
=======================================================================

The simulator iterates the closed-form Type-3 map on the invariant plane (u, v).  This
script produces trajectories that use NONE of that closed form: the full two-qubit
density matrix is iterated through the actual noisy 16-CNOT circuit, by two routes that
share no code with each other beyond the repository's gate list.

  R1  dense 32x32 propagation -- iterated_noisy_pqec.iterate_effective_map("step3"),
      read from the `iterated-noisy-pqec` branch (circuit = verify_analytic_decomposed._fred,
      replacement channel implemented from its definition).
  R2  PennyLane default.mixed -- the repository's own Step-3 executor
      (verify_analytic_decomposed._fred with global_depol_kraus after every CNOT).  The next
      state is rebuilt from all 16 parity-weighted Pauli correlators
      rho' = sum_P <Z_a (x) P> P / (4 <Z_a>), i.e. without assuming Bell-diagonality.
  R3  the one-round closed form already in the repository, pqec_cnot_threshold.F_dec,
      for F_1 only.

Nothing in the repository is modified: the branch file is extracted to a temporary
directory with `git show`.  The output JSON is consumed by tests/reference.test.ts.

Run (from the repository root or anywhere):
    python interactive/type3-dynamics/verification/generate_reference.py
Needs numpy + pennylane (requirements.txt of the repository) and the branch ref
`origin/iterated-noisy-pqec` (override with --ref).
"""

import argparse
import datetime
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "reference_type3.json")

# (eps, q, rounds)
CASES = [
    (0.10, 0.0327, 60),   # the case of interest: rise, then slow later-round decrease
    (0.10, 0.0100, 40),
    (0.10, 0.0330, 40),   # just above the one-round threshold at eps = 0.1
    (0.30, 0.0500, 40),
    (0.50, 0.1000, 40),
    (0.60, 0.1250, 40),
    (0.05, 0.1200, 40),
    (0.90, 0.0200, 40),   # separable input
    (0.40, 0.1400, 40),   # above q_SN: decays towards I/4
    (0.00, 0.0200, 40),   # pure Bell input
    (0.20, 0.0000, 20),   # ideal gadget
    (1.00, 0.0500, 10),   # maximally mixed input
]


def git(*args):
    return subprocess.check_output(["git", "-C", REPO, *args], text=True).strip()


def load_iterated_module(ref):
    src = subprocess.check_output(["git", "-C", REPO, "show", f"{ref}:iterated_noisy_pqec.py"])
    tmp = tempfile.mkdtemp(prefix="type3ref_")
    with open(os.path.join(tmp, "iterated_noisy_pqec.py"), "wb") as fh:
        fh.write(src)
    sys.path.insert(0, tmp)
    import iterated_noisy_pqec as inp  # noqa: E402
    return inp


def make_r2(qml):
    """PennyLane route: parity-weighted Pauli correlators of the repo's Step-3 circuit."""
    import verify_analytic_decomposed as v3

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
        qml.Hadamard(0)
        v3._fred(0, 1, 3, q)
        v3._fred(0, 2, 4, q)
        qml.Hadamard(0)
        return [qml.expval(qml.PauliZ(0) @ qml.Hermitian(P2[l], wires=[1, 2])) for l in labels]

    def step(rho, q):
        c = np.array(node(np.kron(rho, rho), q), dtype=float)
        tau = sum(ci * P2[l] for ci, l in zip(c, labels)) / 4.0
        return tau / c[0]

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
    from pqec_cnot_threshold import F_dec

    inp = load_iterated_module(args.ref)
    r2_step = make_r2(qml)

    cases = []
    for eps, q, N in CASES:
        t = 1.0 - eps
        rho0 = inp.rho_isotropic(t)

        # R1: dense propagation, no early stop (tol = -1)
        res = inp.iterate_effective_map(rho0, q, "step3", tol=-1.0, max_iter=N)
        hist = sorted(res["history"], key=lambda h: h["n"])
        r1 = [[h["p_Phi_plus"], h["p_Phi_minus"], h["p_Psi_plus"], h["p_Psi_minus"]] for h in hist]
        r1_off = max(h["bell_offdiag"] for h in hist)
        assert [h["n"] for h in hist] == list(range(N + 1)), "R1 history incomplete"

        # R2: PennyLane full-state iteration
        rho = rho0.astype(complex)
        r2, r2_off, r2_herm = [], 0.0, 0.0
        for n in range(N + 1):
            pops, off = bell_pops(rho)
            r2.append(pops)
            r2_off = max(r2_off, off)
            if n < N:
                nxt = r2_step(rho, q)
                r2_herm = max(r2_herm, float(np.linalg.norm(nxt - nxt.conj().T)))
                rho = 0.5 * (nxt + nxt.conj().T)

        d12 = float(np.max(np.abs(np.array(r1) - np.array(r2))))
        cases.append(dict(
            eps=eps, q=q, N=N,
            R1_dense=r1, R2_pennylane=r2,
            R3_F1_closed_form=float(F_dec(q, t)),
            max_abs_R1_minus_R2=d12,
            max_bell_offdiag_R1=r1_off, max_bell_offdiag_R2=r2_off,
            max_herm_residual_R2=r2_herm,
        ))
        print(f"eps={eps:.3f} q={q:.4f} N={N:3d}  |R1-R2|={d12:.2e}  "
              f"offdiag R1={r1_off:.1e} R2={r2_off:.1e}  "
              f"F1: R1={r1[1][0]:.15f} R3={F_dec(q, t):.15f}")

    meta = dict(
        generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        repo_head=git("rev-parse", "HEAD"),
        iterated_branch_ref=args.ref,
        iterated_branch_commit=git("rev-parse", args.ref),
        pennylane=qml.__version__, numpy=np.__version__,
        bell_order=["Phi+", "Phi-", "Psi+", "Psi-"],
        note="R1/R2 iterate the full 4x4 state through the noisy 16-CNOT circuit; "
             "neither uses the (u, v) closed form.",
    )
    with open(OUT, "w") as fh:
        json.dump(dict(meta=meta, cases=cases), fh, indent=1)
    print("wrote", os.path.relpath(OUT, REPO))


if __name__ == "__main__":
    main()
