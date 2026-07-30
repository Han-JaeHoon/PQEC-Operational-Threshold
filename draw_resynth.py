"""Draw the Step-4a optimized gadget: the CNOT-reduced circuit that Qiskit produces
from H_a . CSWAP(0;1,3) . CSWAP(0;2,4) . H_a, verified unitary-equivalent."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator

qc = QuantumCircuit(5)                 # 0=a, 1=A1, 2=A2, 3=B1, 4=B2
qc.h(0); qc.cswap(0, 1, 3); qc.cswap(0, 2, 4); qc.h(0)
U_ref = Operator(qc)

# baseline (no optimization) and best optimized
base = transpile(qc, basis_gates=["cx", "u"], optimization_level=0)
best = None
for lvl in (3, 2, 1):
    t = transpile(qc, basis_gates=["cx", "u"], optimization_level=lvl, seed_transpiler=7)
    if Operator(t).equiv(U_ref):
        best = t if best is None or t.count_ops().get("cx", 0) < best.count_ops().get("cx", 0) else best
cx0 = base.count_ops().get("cx", 0)
cx1 = best.count_ops().get("cx", 0)
print(f"baseline CX={cx0}  optimized CX={cx1}  unitary-equiv={Operator(best).equiv(U_ref)}")

fig = best.draw(output="mpl", fold=-1, idle_wires=True, scale=0.9)
fig.suptitle(f"Step 4a: optimized SWAP-test gadget — {cx1} CNOTs "
             f"(from {cx0}), unitary-equivalent to $H_a\\,CSWAP\\,CSWAP\\,H_a$\n"
             f"q0 = ancilla a,  q1,q2 = register A (kept),  q3,q4 = register B",
             fontsize=10)
fig.savefig("circuit_resynth_14cnot.png", dpi=140, bbox_inches="tight")
print("saved circuit_resynth_14cnot.png")
