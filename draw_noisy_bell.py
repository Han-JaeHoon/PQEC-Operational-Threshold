"""
Draw the circuit that prepares  rho_eps = (1-eps)|Phi+><Phi+| + eps I/4.
========================================================================

    q0: |0> --H--*--[ Depol_eps ]--
                 |  [  (global)  ]
    q1: |0> -----X--[           ]--

H . CNOT builds |Phi+>; the two-wire box is the global 2-qubit depolarizing
channel of strength eps (a single qml.QubitChannel with 16 two-qubit-Pauli Kraus
operators).  Saves circuit_noisy_bell.png.

Run:  python draw_noisy_bell.py
"""

import matplotlib.pyplot as plt
import pennylane as qml

from noisy_bell_state import make_noisy_bell

EPS = 0.40   # value shown in the drawing (structure is identical for any eps)


def main():
    fig, ax = qml.draw_mpl(make_noisy_bell, decimals=2, style="pennylane")(EPS)
    ax.set_title(
        r"Noisy Bell prep:  $\rho_\varepsilon=(1-\varepsilon)|\Phi^+\rangle\langle\Phi^+|"
        r"+\varepsilon\, I/4$   ($\varepsilon=%.2f$)" % EPS,
        fontsize=11,
    )
    fig.savefig("circuit_noisy_bell.png", dpi=150, bbox_inches="tight")
    print("saved  circuit_noisy_bell.png")

    # Also print the text circuit for a quick terminal view.
    print("\nText circuit:\n")
    print(qml.draw(make_noisy_bell, decimals=2)(EPS))


if __name__ == "__main__":
    main()
