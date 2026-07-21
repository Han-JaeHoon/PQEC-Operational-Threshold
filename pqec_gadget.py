"""
Step 2 -- The PQEC purification gadget (ideal / noiseless).
===========================================================

The SWAP-test "gadget" is the primitive of PQEC.  Given two identical noisy copies
rho (x) rho of an M-qubit register, an ancilla-controlled SWAP test followed by
reading the ancilla extracts the *purified* component

    P(rho) = rho^2 / Tr[rho^2]          (Raghoonanan & Byrnes, Eq. 7)

which concentrates weight on the dominant eigenvector of rho.  Here the register is
the 2-qubit Bell register (M = 2), so the controlled-SWAP is two parallel Fredkin
gates on an ancilla.  This module implements the gadget as a genuine mixed-state
circuit on top of the Step-1 input state rho_eps and verifies it.

Two equivalent read-outs of the same gadget (both used later, when the gadget is
made noisy in Step 3):

  * state extraction (Eq. 9):
        rho^2 = P_+ rho_+ - P_- rho_-  =  (ancilla |0> block) - (|1> block),
    so purify_once(rho) returns rho^2/Tr[rho^2].
  * observable / parity correlator (the paper's actual protocol): after the gadget
    the joint (ancilla, register) state is  1/2 ( I(x)rho + Z(x)rho^2 ), so
        <O>_purified = <Z(x)O> / <Z(x)I> = Tr(O rho^2)/Tr(rho^2).

On the isotropic input rho_eps the ideal eigenstate |Phi+> is the strictly dominant
eigenvector for every eps < 1 (eigenvalues 1-3eps/4 vs eps/4), so ideal PQEC drives
fidelity and entanglement back to 1 for all eps < 1 -- it even re-entangles a
*separable* input (2/3 <= eps < 1).  The only fixed point is the fully mixed rho=I/4
at eps = 1.

Run:  python pqec_gadget.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

from noisy_bell_state import (
    make_noisy_bell, rho_eps_analytic, PHI_PLUS, O_PHI_PLUS,
    fidelity_phi_plus, purity, X, Y, Z, I2,
)

np.set_printoptions(precision=4, suppress=True)

_YY = np.kron(Y, Y)


# ---------------------------------------------------------------------------
# Entanglement measure
# ---------------------------------------------------------------------------
def concurrence(rho):
    """Wootters concurrence of a 2-qubit state (0 separable ... 1 max entangled)."""
    R = rho @ _YY @ rho.conj() @ _YY
    ev = np.sort(np.sqrt(np.clip(np.real(np.linalg.eigvals(R)), 0, None)))[::-1]
    return max(0.0, ev[0] - ev[1] - ev[2] - ev[3])


# ---------------------------------------------------------------------------
# The genuine SWAP-test gadget (M = 2 register)
#   wire 0 = ancilla,  wires [1,2] = register A (kept),  [3,4] = register B
#   controlled-SWAP of the two 2-qubit registers = two parallel Fredkin gates
# ---------------------------------------------------------------------------
_gadget_dev = qml.device("default.mixed", wires=5)


@qml.qnode(_gadget_dev)
def _gadget_state(rho_AB):
    """Run the SWAP test on two 2-qubit copies; return the joint (ancilla, A) state."""
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])   # two noisy copies rho (x) rho
    qml.Hadamard(0)
    qml.ctrl(qml.SWAP, control=0)(wires=[1, 3])          # Fredkin on qubit pair 0
    qml.ctrl(qml.SWAP, control=0)(wires=[2, 4])          # Fredkin on qubit pair 1
    qml.Hadamard(0)
    return qml.density_matrix(wires=[0, 1, 2])           # ancilla + register A (8x8)


@qml.qnode(_gadget_dev)
def _gadget_obs(rho_AB, O):
    """Same gadget; measure the parity correlators <Z(x)O> and <Z(x)I> on (anc, A)."""
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
    qml.Hadamard(0)
    qml.ctrl(qml.SWAP, control=0)(wires=[1, 3])
    qml.ctrl(qml.SWAP, control=0)(wires=[2, 4])
    qml.Hadamard(0)
    return (qml.expval(qml.PauliZ(0) @ qml.Hermitian(O, wires=[1, 2])),
            qml.expval(qml.PauliZ(0)))


# ---------------------------------------------------------------------------
# Purification API
# ---------------------------------------------------------------------------
def purify_once(rho):
    """One ideal PQEC round on a 2-qubit state via the genuine gadget circuit.
    Returns rho^2/Tr[rho^2] extracted from the ancilla |0>/|1> blocks (Eq. 9).
    Already-pure states are fixed points (skip to avoid 0/0)."""
    if purity(rho) > 1 - 1e-12:
        return rho
    sigma = _gadget_state(np.kron(rho, rho))       # joint (anc, A), 8x8
    rho_sq = sigma[:4, :4] - sigma[4:, 4:]         # |0> block minus |1> block == rho^2
    return rho_sq / np.trace(rho_sq)


def purify_rounds(rho, ell):
    """Apply ell ideal rounds; ell rounds consume N = 2^ell copies -> rho^N/Tr[rho^N]."""
    for _ in range(ell):
        rho = purify_once(rho)
    return rho


def obs_purified(rho, O=O_PHI_PLUS):
    """Purified expectation value <O> = Tr(O rho^2)/Tr(rho^2), measured via the
    ancilla-parity correlator on the genuine gadget."""
    zO, zI = _gadget_obs(np.kron(rho, rho), O)
    return float(zO) / float(zI)


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
def verify(tol=1e-13, verbose=True):
    """Check that the gadget extracts rho^2/Tr[rho^2] (state) and Tr(O rho^2)/Tr(rho^2)
    (observable), on random states and on the Step-1 input rho_eps."""
    rng = np.random.default_rng(1)
    ok = True

    # (1) random 2-qubit states: circuit == analytic purification map
    worst_state = 0.0
    for _ in range(500):
        A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        rho = A @ A.conj().T
        rho = rho / np.trace(rho)
        circ = purify_once(rho)
        alg = (rho @ rho) / np.trace(rho @ rho)
        worst_state = max(worst_state, np.abs(circ - alg).max())
    ok = ok and worst_state < tol
    if verbose:
        print(f"  (1) purify_once == rho^2/Tr[rho^2]   (500 random states): "
              f"max err {worst_state:.2e}  {'PASS' if worst_state < tol else 'FAIL'}")

    # (2) observable protocol: <Z(x)O>/<Z(x)I> == Tr(O rho^2)/Tr(rho^2)
    worst_obs = 0.0
    for _ in range(200):
        A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        rho = A @ A.conj().T
        rho = rho / np.trace(rho)
        meas = obs_purified(rho, O_PHI_PLUS)
        r2 = rho @ rho
        analytic = float(np.real(np.trace(O_PHI_PLUS @ r2) / np.trace(r2)))
        worst_obs = max(worst_obs, abs(meas - analytic))
    ok = ok and worst_obs < tol
    if verbose:
        print(f"  (2) <Z(x)O>/<Z(x)I> == Tr(O rho^2)/Tr(rho^2)  (200 random): "
              f"max err {worst_obs:.2e}  {'PASS' if worst_obs < tol else 'FAIL'}")

    # (3) on rho_eps: one round equals the analytic purified isotropic state
    worst_eps = 0.0
    for eps in np.linspace(0.0, 0.99, 40):
        rho = make_noisy_bell(eps)
        circ = purify_once(rho)
        alg = (rho @ rho) / np.trace(rho @ rho)
        worst_eps = max(worst_eps, np.abs(circ - alg).max())
    ok = ok and worst_eps < tol
    if verbose:
        print(f"  (3) purify_once(rho_eps) == rho_eps^2/Tr   (40 eps values): "
              f"max err {worst_eps:.2e}  {'PASS' if worst_eps < tol else 'FAIL'}")
    return ok


# ===========================================================================
def main():
    print("=" * 78)
    print(" Step 2 -- ideal PQEC purification gadget on rho_eps")
    print("=" * 78)

    print("\n Verification:")
    all_ok = verify()

    # -- recovery demo: fidelity + concurrence over rounds --------------------
    print("\n" + "-" * 78)
    print(" Recovery of |Phi+> from rho_eps (ideal gadget):")
    print("-" * 78)
    print(f"  {'eps':>5} {'input F':>8} {'input C':>8} | "
          f"{'F(l=3)':>8} {'C(l=3)':>8}  note")
    for eps in [0.30, 0.50, 2 / 3, 0.80, 0.95]:
        rho0 = make_noisy_bell(eps)
        r3 = purify_rounds(rho0, 3)
        note = "input separable" if concurrence(rho0) < 1e-9 else "input entangled"
        print(f"  {eps:>5.2f} {fidelity_phi_plus(rho0):>8.4f} {concurrence(rho0):>8.4f} | "
              f"{fidelity_phi_plus(r3):>8.4f} {concurrence(r3):>8.4f}  {note}")
    print("  -> |Phi+> is the dominant eigenvector for every eps < 1, so ideal PQEC")
    print("     restores F, C -> 1 for all eps < 1 -- even re-entangling a separable")
    print("     input (eps >= 2/3).  Only rho=I/4 at eps=1 is a fixed point.")

    # -- Figure ---------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    rounds = range(7)

    ax = axes[0]
    for eps, c in [(0.30, "C0"), (0.50, "C1"), (2 / 3, "C2"), (0.80, "C3")]:
        F = []
        r = make_noisy_bell(eps)
        F.append(fidelity_phi_plus(r))
        for _ in range(6):
            r = purify_once(r)
            F.append(fidelity_phi_plus(r))
        ax.plot(rounds, F, "-o", ms=4, color=c, label=f"$\\varepsilon$={eps:.2f}")
    ax.axhline(1.0, color="k", ls=":", lw=1)
    ax.set_xlabel(r"purification rounds  $\ell$")
    ax.set_ylabel(r"fidelity  $\langle\Phi^+|\sigma_\ell|\Phi^+\rangle$")
    ax.set_title("(a) Fidelity recovery (ideal gadget)")
    ax.set_ylim(0.2, 1.03)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    for eps, c in [(0.30, "C0"), (0.50, "C1"), (2 / 3, "C2"), (0.80, "C3")]:
        C = []
        r = make_noisy_bell(eps)
        C.append(concurrence(r))
        for _ in range(6):
            r = purify_once(r)
            C.append(concurrence(r))
        ax.plot(rounds, C, "-s", ms=4, color=c, label=f"$\\varepsilon$={eps:.2f}")
    ax.axhline(1.0, color="k", ls=":", lw=1)
    ax.set_xlabel(r"purification rounds  $\ell$")
    ax.set_ylabel("concurrence")
    ax.set_title("(b) Entanglement recovery (even from separable inputs)")
    ax.set_ylim(0.0, 1.03)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig("pqec_gadget_recovery.png", dpi=140)
    print("\n  saved  pqec_gadget_recovery.png")

    print("\n" + "=" * 78)
    print(f"  STEP 2 {'OK' if all_ok else 'FAILED'}")
    print("=" * 78)


if __name__ == "__main__":
    main()
