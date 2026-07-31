"""Draw the Step-5 hardware-efficient PQC ansatz structure (default B=6).

Three stages: a full Rot layer, then per-CNOT Rot pairs (ancilla-centric CNOT
schedule), then a final full Rot layer -> 30 + 6B parameters. In the auxiliary
observable study the same circuit runs on default.mixed with a 2-qubit depolarizing
channel after each CNOT.

Run:  python draw_pqc_ansatz.py
"""
import numpy as np
import pennylane as qml
from pqc_common import ansatz_ops, N_WIRES, _PAIR_CYCLE

B = 6
ops, npar = ansatz_ops(B)
theta = np.random.default_rng(0).uniform(-np.pi, np.pi, npar)
dev = qml.device("default.qubit", wires=N_WIRES)


@qml.qnode(dev)
def circ(params):
    p = 0
    for w in range(N_WIRES):                       # stage 1: full Rot layer
        qml.Rot(params[p], params[p + 1], params[p + 2], wires=w); p += 3
    qml.Barrier(wires=range(N_WIRES))
    for k in range(B):                             # stage 2: per-CNOT Rot pairs
        c, t = _PAIR_CYCLE[k % len(_PAIR_CYCLE)]
        qml.CNOT(wires=[c, t])
        qml.Rot(params[p], params[p + 1], params[p + 2], wires=c); p += 3
        qml.Rot(params[p], params[p + 1], params[p + 2], wires=t); p += 3
    qml.Barrier(wires=range(N_WIRES))
    for w in range(N_WIRES):                       # stage 3: full Rot layer
        qml.Rot(params[p], params[p + 1], params[p + 2], wires=w); p += 3
    return qml.state()


def main():
    fig, ax = qml.draw_mpl(circ, decimals=None, style="pennylane",
                           wire_order=range(N_WIRES))(theta)
    ax.set_title(
        f"Step-5 PQC ansatz  (B={B} CNOTs, {npar} params = 30+6B)\n"
        "wires: 0=ancilla a, 1=A1, 2=A2 (kept), 3=B1, 4=B2 (discarded)\n"
        "stage 1: full Rot layer  |  stage 2: per-CNOT Rot pairs  |  stage 3: full Rot layer",
        fontsize=10)
    fig.savefig("circuit_pqc_ansatz.png", dpi=140, bbox_inches="tight")
    print(f"saved circuit_pqc_ansatz.png ; params = {npar}")


if __name__ == "__main__":
    main()
