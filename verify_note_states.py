"""
Verify the step-by-step states of the Step-3a note against the circuit.
=======================================================================

Note: "Effect of Three-Qubit Global Depolarizing Noise after Each CSWAP Gate in
One-Round PQEC" (Jaehun Han, July 22).  The note evolves the five-qubit state
    ordering (a, A1, A2, B1, B2),  wires (0, 1, 2, 3, 4)
through one PQEC round with a 3-qubit GLOBAL replacement depolarizing channel
    D_g(sigma) = (1-g) sigma + g (I_Q/8) (x) Tr_Q(sigma)
right after each controlled-SWAP, and gives a closed form after every step.

With  s = 1-g,  R = rho_A (x) rho_B,  S1 = SWAP_{A1B1},  S2 = SWAP_{A2B2},
SAB = S1 S2 (which commutes with R since the two copies are identical):

  sig0  = |+><+|_a (x) R                                              (after 1st H)
  sig1  = 1/2 [[R, R S1],[S1 R, S1 R S1]]                            (after CSWAP1) (1)
  sig1' = s sig1 + (1-s) I^5/32                                       (after depol1) (4)
  sig2  = (s/2)[[R, R SAB],[R SAB, R]] + (1-s) I^5/32                 (after CSWAP2) (6)
  sig2' = (s^2/2)[[R, R SAB],[R SAB, R]] + (1-s^2) I^5/32             (after depol2) (13)
  sigout= (s^2/2)[[R+R SAB, 0],[0, R-R SAB]] + (1-s^2) I^5/32         (after last H) (14)

Blocks are in the ancilla basis: [[X,Y],[Z,W]] = |0><0|X + |0><1|Y + |1><0|Z +
|1><1|W.  We build each closed form and compare it to the genuine mixed-state
circuit stopped after the corresponding gate (density_matrix on all 5 wires).

Run:  python verify_note_states.py
"""

import numpy as np
import pennylane as qml

from noisy_bell_state import make_noisy_bell
from pqec_gadget_noise import global_depol3_kraus

np.set_printoptions(precision=4, suppress=True)


# --- register operators, ordering (A1, A2, B1, B2) = qubits (0,1,2,3) --------
def _swap_qubits(i, j, n=4):
    dim = 2 ** n
    M = np.zeros((dim, dim))
    for x in range(dim):
        b = [(x >> (n - 1 - k)) & 1 for k in range(n)]
        b[i], b[j] = b[j], b[i]
        y = sum(bit << (n - 1 - k) for k, bit in enumerate(b))
        M[y, x] = 1
    return M


S1 = _swap_qubits(0, 2)          # SWAP A1,B1
S2 = _swap_qubits(1, 3)          # SWAP A2,B2
SAB = S1 @ S2                    # full register SWAP

# ancilla-block builder: |0><0|(x)X + |0><1|(x)Y + |1><0|(x)Z + |1><1|(x)W
_E00 = np.array([[1, 0], [0, 0]]); _E01 = np.array([[0, 1], [0, 0]])
_E10 = _E01.T;                     _E11 = np.array([[0, 0], [0, 1]])


def _ab(X, Y, Z, W):
    return np.kron(_E00, X) + np.kron(_E01, Y) + np.kron(_E10, Z) + np.kron(_E11, W)


_PLUS = np.array([[.5, .5], [.5, .5]])
_I32 = np.eye(32)

_dev = qml.device("default.mixed", wires=5)


def note_states(eps, g):
    """The six closed-form states of the note."""
    s = 1 - g
    rho = make_noisy_bell(eps)
    R = np.kron(rho, rho)                      # rho_A on (A1,A2), rho_B on (B1,B2)
    K = R @ SAB
    return R, [
        ("sig0  (after 1st H)", np.kron(_PLUS, R)),
        ("sig1  (after CSWAP1)", 0.5 * _ab(R, R @ S1, S1 @ R, S1 @ R @ S1)),
        ("sig1' (after depol1)", s * (0.5 * _ab(R, R @ S1, S1 @ R, S1 @ R @ S1))
                                 + (1 - s) * _I32 / 32),
        ("sig2  (after CSWAP2)", (s / 2) * _ab(R, K, K, R) + (1 - s) * _I32 / 32),
        ("sig2' (after depol2)", (s ** 2 / 2) * _ab(R, K, K, R)
                                 + (1 - s ** 2) * _I32 / 32),
        ("sigout(after last H)", (s ** 2 / 2) * _ab(R + K, 0 * R, 0 * R, R - K)
                                 + (1 - s ** 2) * _I32 / 32),
    ]


def circuit_state(R, g, nsteps):
    """The genuine circuit stopped after `nsteps` gates; full 5-qubit density matrix."""
    @qml.qnode(_dev)
    def c():
        qml.QubitDensityMatrix(R, wires=[1, 2, 3, 4])
        if nsteps >= 1: qml.Hadamard(0)
        if nsteps >= 2: qml.ctrl(qml.SWAP, control=0)(wires=[1, 3])   # CSWAP(a;A1,B1)
        if nsteps >= 3: qml.QubitChannel(global_depol3_kraus(g), wires=[0, 1, 3])
        if nsteps >= 4: qml.ctrl(qml.SWAP, control=0)(wires=[2, 4])   # CSWAP(a;A2,B2)
        if nsteps >= 5: qml.QubitChannel(global_depol3_kraus(g), wires=[0, 2, 4])
        if nsteps >= 6: qml.Hadamard(0)
        return qml.density_matrix(wires=[0, 1, 2, 3, 4])
    return c()


def worst_over_steps(eps, g):
    R, note = note_states(eps, g)
    return max(np.abs(circuit_state(R, g, k + 1) - M).max() for k, (_, M) in enumerate(note))


def main():
    print("=" * 74)
    print(" Step-3a note states  vs  circuit  (ordering a,A1,A2,B1,B2 = wires 0..4)")
    print("=" * 74)

    # commutation the note relies on
    rho = make_noisy_bell(0.30); R = np.kron(rho, rho)
    print(f"\n  [R, SAB] = 0 ?   max|R SAB - SAB R| = {np.abs(R @ SAB - SAB @ R).max():.2e}")

    eps, g = 0.30, 0.20
    print(f"\n  eps={eps}, g={g}:   max|circuit - note| after each step")
    R, note = note_states(eps, g)
    for k, (name, M) in enumerate(note, start=1):
        err = np.abs(circuit_state(R, g, k) - M).max()
        print(f"    step {k}:  {name:22s}  {err:.2e}")

    print("\n  robustness over (eps, g) incl. edge cases:")
    worst = 0.0
    for eps, g in [(0.0, 0.0), (0.2, 0.05), (0.4, 0.30), (0.6, 0.50),
                   (0.5, 0.90), (0.9, 0.15)]:
        w = worst_over_steps(eps, g)
        worst = max(worst, w)
        print(f"    eps={eps:.2f} g={g:.2f}:  worst over all 6 steps = {w:.2e}")

    print("\n" + "=" * 74)
    print(f"  NOTE STATES MATCH CIRCUIT AT EVERY STEP: "
          f"{'YES' if worst < 1e-12 else 'CHECK'}  (worst {worst:.1e})")
    print("=" * 74)


if __name__ == "__main__":
    main()
