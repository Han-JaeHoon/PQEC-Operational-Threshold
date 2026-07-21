"""
Draw the full decomposed, noisy PQEC gadget (Step 3b).
======================================================

Two figures:
  * circuit_decomposed_fredkin.png -- a SINGLE decomposed Fredkin building block
    (CNOT + Clifford+T Toffoli + CNOT) with per-gate noise, for readability.
  * circuit_decomposed_gadget.png  -- the FULL gadget: ancilla H, two decomposed
    Fredkins, H, measurement, with a depolarizing box after every native gate.

Noise shown at p1 = p2 = 0.05 so the channel boxes are visible.

Run:  python draw_decomposed.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

from noisy_bell_state import make_noisy_bell, O_PHI_PLUS
from pqec_decomposed_noise import _fredkin, _g1

P1 = P2 = 0.05

_dev5 = qml.device("default.mixed", wires=5)
_dev3 = qml.device("default.mixed", wires=3)


@qml.qnode(_dev5)
def _full_gadget(rho_AB, p1, p2):
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
    _g1(qml.Hadamard, 0, p1)
    _fredkin(0, 1, 3, p1, p2)
    _fredkin(0, 2, 4, p1, p2)
    _g1(qml.Hadamard, 0, p1)
    return qml.expval(qml.PauliZ(0))


@qml.qnode(_dev3)
def _one_fredkin(p1, p2):
    # control = 0, swap targets 1,2  (decomposed, with per-gate noise)
    _fredkin(0, 1, 2, p1, p2)
    return qml.expval(qml.PauliZ(0))


def main():
    # (1) single decomposed Fredkin building block
    fig, ax = qml.draw_mpl(_one_fredkin, decimals=None, style="pennylane")(P1, P2)
    ax.set_title("Decomposed Fredkin  CSWAP(0;1,2) = CNOT + Clifford+T Toffoli + CNOT\n"
                 "(Depol boxes = per-gate noise; 2-qubit 'QubitChannel' after CNOTs)",
                 fontsize=10)
    fig.savefig("circuit_decomposed_fredkin.png", dpi=150, bbox_inches="tight")
    print("saved  circuit_decomposed_fredkin.png")

    # (2) full gadget
    rho_AB = np.kron(make_noisy_bell(0.40), make_noisy_bell(0.40))
    fig, ax = qml.draw_mpl(_full_gadget, decimals=None, style="pennylane",
                           wire_order=[0, 1, 2, 3, 4])(rho_AB, P1, P2)
    ax.set_title("Full decomposed noisy PQEC gadget: anc H | Fredkin(0;1,3) | "
                 "Fredkin(0;2,4) | anc H | measure", fontsize=10)
    fig.savefig("circuit_decomposed_gadget.png", dpi=130, bbox_inches="tight")
    print("saved  circuit_decomposed_gadget.png")

    # text view of the full gadget
    print("\nText circuit (full gadget):\n")
    print(qml.draw(_full_gadget, decimals=None, max_length=200)(rho_AB, P1, P2))


if __name__ == "__main__":
    main()
