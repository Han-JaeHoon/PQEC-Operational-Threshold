"""
Draw the ideal PQEC purification gadget (M = 2 register).
=========================================================

    anc: |0> --H--*--*--H--(measure Z / parity)
                  |  |
    A1 : ~~~~~~~~[x]~|~~~~~~   register A (kept)
    A2 : ~~~~~~~~~~~[x]~~~~~~
    B1 : ~~~~~~~~[x]~|~~~~~~   register B (discarded)
    B2 : ~~~~~~~~~~~[x]~~~~~~

Two 2-qubit copies rho (x) rho enter on [A1,A2] and [B1,B2]; the ancilla-controlled
SWAP is two parallel Fredkin gates.  Saves circuit_pqec_gadget.png.

Run:  python draw_pqec_gadget.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

from noisy_bell_state import make_noisy_bell

_dev = qml.device("default.mixed", wires=5)


@qml.qnode(_dev)
def _gadget(rho_AB):
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
    qml.Hadamard(0)
    qml.ctrl(qml.SWAP, control=0)(wires=[1, 3])
    qml.ctrl(qml.SWAP, control=0)(wires=[2, 4])
    qml.Hadamard(0)
    return qml.expval(qml.PauliZ(0))


def main():
    rho = make_noisy_bell(0.40)
    rho_AB = np.kron(rho, rho)
    fig, ax = qml.draw_mpl(_gadget, decimals=None, style="pennylane",
                           wire_order=[0, 1, 2, 3, 4])(rho_AB)
    ax.set_title("PQEC SWAP-test gadget  (M=2): anc + two Fredkins on register A / B",
                 fontsize=11)
    fig.savefig("circuit_pqec_gadget.png", dpi=150, bbox_inches="tight")
    print("saved  circuit_pqec_gadget.png")
    print("\nText circuit:\n")
    print(qml.draw(_gadget, decimals=None)(rho_AB))


if __name__ == "__main__":
    main()
