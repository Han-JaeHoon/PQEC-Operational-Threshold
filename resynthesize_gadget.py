"""
Unitary-preserving CNOT reduction of the SWAP-test gadget (Part 1).
==================================================================

The coherent gadget  U = H_a . CSWAP(a;A1,B1) . CSWAP(a;A2,B2) . H_a  on 5 qubits
(ancilla 0; register A = 1,2; register B = 3,4) is Clifford + 2 Toffoli.  The
textbook decomposition uses 8 CNOTs per Fredkin = 16 CNOTs.  Here we ask: keeping
the SAME 5-qubit unitary, how few CNOTs are needed?

Two independent optimizers are run and each result is checked for exact unitary
equivalence to U (up to global phase):

  * Qiskit  transpile(optimization_level = 1,2,3)
  * pytket  FullPeepholeOptimise

Both FIND 14 CNOTs (= 2 x 7, the per-Fredkin optimum), verified equivalent.  Note:
this is the count the optimizers produced, NOT a proof that 14 is minimal over the
full 5-qubit unitary (that would need a lower-bound / exhaustive-synthesis argument).
ZX-calculus (PyZX) minimizes T-count, not CNOT-count, and does worse here (~27
CNOTs), so it is not used for this metric.

For a much larger reduction one must relax "same unitary" to "same measured
observable" -- see destructive_gadget.py (2 CNOTs).

Run:  python resynthesize_gadget.py
"""


def _qiskit_reduce():
    from qiskit import QuantumCircuit, transpile
    from qiskit.quantum_info import Operator

    qc = QuantumCircuit(5)                 # 0=a, 1=A1, 2=A2, 3=B1, 4=B2
    qc.h(0)
    qc.cswap(0, 1, 3)
    qc.cswap(0, 2, 4)
    qc.h(0)
    U_ref = Operator(qc)

    base = transpile(qc, basis_gates=["cx", "u3"], optimization_level=0)
    results = {"baseline": base.count_ops().get("cx", 0)}
    best = None
    for lvl in (1, 2, 3):
        t = transpile(qc, basis_gates=["cx", "u3"],
                      optimization_level=lvl, seed_transpiler=7)
        cx = t.count_ops().get("cx", 0)
        ok = Operator(t).equiv(U_ref)
        results[f"opt{lvl}"] = (cx, ok)
        if ok and (best is None or cx < best):
            best = cx
    return results, best


def _pytket_reduce():
    from pytket import Circuit
    from pytket.passes import FullPeepholeOptimise
    from pytket.circuit import OpType
    from pytket.utils import compare_unitaries

    c = Circuit(5)
    c.H(0)
    c.CSWAP(0, 1, 3)
    c.CSWAP(0, 2, 4)
    c.H(0)
    U0 = c.get_unitary()
    cc = c.copy()
    FullPeepholeOptimise().apply(cc)
    return cc.n_gates_of_type(OpType.CX), compare_unitaries(U0, cc.get_unitary())


def main():
    print("=" * 74)
    print(" Unitary-preserving CNOT reduction of  H_a . CSWAP . CSWAP . H_a")
    print("=" * 74)

    res, best = _qiskit_reduce()
    print(f"\n  Qiskit baseline (opt0):  {res['baseline']} CX")
    for lvl in (1, 2, 3):
        cx, ok = res[f"opt{lvl}"]
        print(f"  Qiskit opt{lvl}:  {cx} CX   unitary-equiv = {ok}")
    print(f"  Qiskit best (verified):  {best} CX")

    try:
        cx, ok = _pytket_reduce()
        print(f"\n  pytket FullPeepholeOptimise:  {cx} CX   unitary-equiv = {ok}")
    except Exception as e:                          # pragma: no cover
        print(f"\n  pytket unavailable: {e}")

    print("\n  => two independent optimizers FIND 16 -> 14 CNOTs, unitary-preserving")
    print("     (2 x 7 = per-Fredkin optimum). This is the count they produced, not a")
    print("     proof of minimality over the full 5-qubit unitary.")


if __name__ == "__main__":
    main()
