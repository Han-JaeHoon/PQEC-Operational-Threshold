"""
Draw the Step 3a gadget: 3-qubit global depolarizing after each CSWAP.
=====================================================================

Wires: 0 = ancilla a, [1,2] = register A (kept), [3,4] = register B (discarded).
    a  : |+> --*--------*----------- H -- <Z>
              |        |
    A1 : -----x--[G]---|-------------      (G = 3-qubit global depol on {a,A1,B1})
    B1 : -----x--[G]---|-------------
    A2 : -------------x--[G']--------      (G'= 3-qubit global depol on {a,A2,B2})
    B2 : -------------x--[G']--------

Saves circuit_gadget_noise.png.  Run:  python draw_gadget_noise.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

from noisy_bell_state import make_noisy_bell
from pqec_gadget_noise import global_depol3_kraus

G = 0.10  # shown so the depolarizing boxes are visible

_dev = qml.device("default.mixed", wires=5)


@qml.qnode(_dev)
def _gadget(rho_AB, g):
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])   # two noisy copies rho (x) rho
    qml.Hadamard(0)
    qml.ctrl(qml.SWAP, control=0)(wires=[1, 3])          # CSWAP(a; A1,B1)
    qml.QubitChannel(global_depol3_kraus(g), wires=[0, 1, 3])   # 3q global depol on Q1
    qml.ctrl(qml.SWAP, control=0)(wires=[2, 4])          # CSWAP(a; A2,B2)
    qml.QubitChannel(global_depol3_kraus(g), wires=[0, 2, 4])   # 3q global depol on Q2
    qml.Hadamard(0)
    return qml.expval(qml.PauliZ(0))


def main():
    rho_AB = np.kron(make_noisy_bell(0.40), make_noisy_bell(0.40))
    fig, ax = qml.draw_mpl(_gadget, decimals=None, style="pennylane",
                           wire_order=[0, 1, 2, 3, 4])(rho_AB, G)
    ax.set_title("Step 3a: 3-qubit global depolarizing (QubitChannel) after each CSWAP\n"
                 "wires: 0=ancilla, [1,2]=register A (kept), [3,4]=register B",
                 fontsize=10)
    fig.savefig("circuit_gadget_noise.png", dpi=150, bbox_inches="tight")
    print("saved  circuit_gadget_noise.png")
    print("\nText circuit:\n")
    print(qml.draw(_gadget, decimals=None, max_length=160)(rho_AB, G))


if __name__ == "__main__":
    main()
