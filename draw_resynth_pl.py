"""Draw the PennyLane implementation of the Step-4 14-CNOT gadget (pqec_resynth_noise)."""
import matplotlib
matplotlib.use("Agg")
import pennylane as qml
from pqec_resynth_noise import _apply, N_CX

fig, ax = qml.draw_mpl(_apply, decimals=2, style="pennylane",
                       wire_order=range(5))(0.0)   # eps2=0: show gate structure
ax.set_title(
    f"Step 4 gadget, PennyLane implementation (pqec_resynth_noise.py) — {N_CX} CNOTs\n"
    "0 = ancilla a,  1,2 = register A (kept),  3,4 = register B (discarded)\n"
    "U3 = single-qubit rotation; in the threshold analysis a 2-qubit depol "
    "$\\varepsilon_2$ is inserted after each CNOT",
    fontsize=10)
fig.savefig("circuit_resynth_pennylane.png", dpi=140, bbox_inches="tight")
print("saved circuit_resynth_pennylane.png ; N_CX =", N_CX)
