"""
Three circuit diagrams for the CNOT-only noise model (single-qubit gates ideal).
===============================================================================

  1. circuit_cswap_decomp.png        -- the CSWAP (Fredkin) decomposition alone
  2. circuit_swaptest_decomp.png     -- the full SWAP test using that decomposition
  3. circuit_swaptest_cnot_noise.png -- the same SWAP test with a two-qubit
                                        depolarizing channel after EACH CNOT
                                        (single-qubit gates left ideal)

Uses the textbook decomposition CSWAP(c;x,y) = CNOT(x->y) . Toffoli(c,y;x) .
CNOT(x->y) with the Clifford+T Toffoli (8 CNOTs per Fredkin).  Reuses _fredkin from
pqec_decomposed_noise.py with p1=0 (no single-qubit noise); p2 turns on the
two-qubit depolarizing after each CNOT.

Run:  python draw_cnot_noise.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

from noisy_bell_state import make_noisy_bell
from pqec_decomposed_noise import _fredkin

P2 = 0.10  # shown so the two-qubit depolarizing boxes are visible

_dev3 = qml.device("default.mixed", wires=3)
_dev5 = qml.device("default.mixed", wires=5)


# 1) single decomposed CSWAP (no noise) --------------------------------------
@qml.qnode(_dev3)
def _cswap_decomp():
    _fredkin(0, 1, 2, 0.0, 0.0)          # control 0; swap 1,2; ideal
    return qml.expval(qml.PauliZ(0))


_ALL = [0, 1, 2, 3, 4]


def _bar():
    qml.Barrier(wires=_ALL, only_visual=True)   # visual stage separator


# 2) full SWAP test with the decomposed CSWAP (no noise) ----------------------
# barriers delimit:  state prep | CSWAP1 decomp | CSWAP2 decomp | final Hadamard
@qml.qnode(_dev5)
def _swaptest_decomp(rho_AB):
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
    qml.Hadamard(0)                      # ancilla -> |+>  (state prep)
    _bar()
    _fredkin(0, 1, 3, 0.0, 0.0)          # CSWAP(a; A1,B1)
    _bar()
    _fredkin(0, 2, 4, 0.0, 0.0)          # CSWAP(a; A2,B2)
    _bar()
    qml.Hadamard(0)                      # final Hadamard
    return qml.expval(qml.PauliZ(0))


# 3) full SWAP test with 2-qubit depolarizing after each CNOT (CNOT-only) -----
@qml.qnode(_dev5)
def _swaptest_cnot_noise(rho_AB, p2):
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
    qml.Hadamard(0)                      # ancilla -> |+>  (state prep)
    _bar()
    _fredkin(0, 1, 3, 0.0, p2)           # single-qubit gates ideal (p1=0),
    _bar()
    _fredkin(0, 2, 4, 0.0, p2)           # 2-qubit depol after each CNOT
    _bar()
    qml.Hadamard(0)                      # final Hadamard
    return qml.expval(qml.PauliZ(0))


def _save(qnode, args, fname, title, wire_order=None):
    kw = {} if wire_order is None else {"wire_order": wire_order}
    fig, ax = qml.draw_mpl(qnode, decimals=None, style="pennylane", **kw)(*args)
    ax.set_title(title, fontsize=10)
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("saved ", fname)


def main():
    rho_AB = np.kron(make_noisy_bell(0.40), make_noisy_bell(0.40))

    _save(_cswap_decomp, (), "circuit_cswap_decomp.png",
          "CSWAP decomposition  CSWAP(0;1,2) = CNOT(1->2) . Toffoli(0,2;1) . CNOT(1->2)\n"
          "(Clifford+T Toffoli; 8 CNOTs, single-qubit gates H/T/T†)")

    _save(_swaptest_decomp, (rho_AB,), "circuit_swaptest_decomp.png",
          "SWAP test with decomposed CSWAP (ideal, no noise)\n"
          "wires: 0=ancilla, [1,2]=register A (kept), [3,4]=register B",
          wire_order=[0, 1, 2, 3, 4])

    _save(_swaptest_cnot_noise, (rho_AB, P2), "circuit_swaptest_cnot_noise.png",
          "SWAP test: two-qubit depolarizing (QubitChannel) after each CNOT, "
          "single-qubit gates ideal\n"
          "wires: 0=ancilla, [1,2]=register A (kept), [3,4]=register B",
          wire_order=[0, 1, 2, 3, 4])


if __name__ == "__main__":
    main()
