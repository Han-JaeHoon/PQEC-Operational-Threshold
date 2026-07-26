# Reducing the CNOT count of the SWAP-test gadget

The one-round gadget `H_a · CSWAP(a;A1,B1) · CSWAP(a;A2,B2) · H_a` (ancilla `a` =
wire 0, retained register `A=(A1,A2)` = wires 1,2, discarded `B=(B1,B2)` = wires
3,4) costs **16 CNOTs** with the textbook 8-CNOT Fredkin decomposition. We reduce
this in two senses of "same role":

- **(1) Same unitary** — keep the exact 5-qubit unitary, fewer CNOTs.
- **(2) Same measured observable** — reproduce only `F = Tr(Oρ²)/Tr(ρ²)`, which is
  all PQEC uses; this allows a much larger reduction.

All results are checked rigorously (unitary equivalence for (1); machine-precision
observable equivalence for (2)).

---

## Part 1 — Unitary-preserving reduction: 16 → 14 CNOTs

Script: [`resynthesize_gadget.py`](resynthesize_gadget.py). The gadget is
Clifford + 2 Toffoli. Two independent optimizers were run, and each output was
verified to implement the **same 5-qubit unitary** (up to global phase):

| optimizer | CNOTs | unitary-equivalent |
|-----------|:-----:|:------------------:|
| textbook baseline | 16 | — |
| Qiskit `optimization_level=1` | 16 | ✓ |
| Qiskit `optimization_level=2,3` | **14** | ✓ |
| pytket `FullPeepholeOptimise` | **14** | ✓ |
| PyZX `full_reduce` (ZX) | 27 | ✓ (T-count tool — worse for CNOTs) |

Both peephole optimizers converge to **14 CNOTs = 2 × 7**, i.e. the per-Fredkin
optimum (cf. the 7-CNOT Fredkin, Cruz–Murta arXiv:2305.18128); the shared ancilla
control gives no further CNOT saving at fixed unitary. PyZX minimizes T-count, not
CNOT-count, so it is not the right tool for this metric.

**Takeaway:** a modest, safe 12.5% CNOT reduction is available for free while
keeping the exact coherent gadget.

---

## Part 2 — Measurement-equivalent gadget: 2 CNOTs

Script: [`destructive_gadget.py`](destructive_gadget.py). PQEC only needs the
purified observable `F = Tr(Oρ²)/Tr(ρ²)`. This is exactly what the **destructive
SWAP test / virtual-distillation measurement** returns, with no ancilla and no
controlled-SWAP (Garcia-Escartín & Chamorro-Posada 2013; Huggins et al. 2021).

**Construction** (registers `A=(A1,A2)`, `B=(B1,B2)`, qubits (0,1,2,3), no ancilla):

```
SWAP_reg = SWAP_{A1B1} · SWAP_{A2B2}
V        = [H_{A1} CNOT_{A1→B1}] [H_{A2} CNOT_{A2→B2}]      (Bell-basis change; 2 CNOTs)
Π        = V SWAP_reg V†                                    (diagonalized SWAP)
O_den = Π                                          ⟨O_den⟩ = Tr(ρ²)
O_num = V [ ½(O_A SWAP_reg + SWAP_reg O_A) ] V†     ⟨O_num⟩ = Tr(Oρ²)
F     = ⟨O_num⟩ / ⟨O_den⟩
```

The experiment applies the (noisy) Bell change `V` and reads the fixed ideal-frame
observables `O_den`, `O_num`.

**Verification (ideal, `ε₂=0`).** `F_dest` equals both the controlled-SWAP gadget
and the exact trace `Tr(Oρ²)/Tr(ρ²)` to `2e-16`, for the fidelity projector
`O=|Φ⁺⟩⟨Φ⁺|` **and** a generic Pauli `O=ZZ`:

| `ε` | `Φ⁺`: dest / cSWAP / exact | `ZZ`: dest / exact |
|-----|---------------------------|--------------------|
| 0.2 | 0.98973 / 0.98973 / 0.98973 | 0.98630 / 0.98630 |
| 0.4 | 0.94231 / 0.94231 / 0.94231 | 0.92308 / 0.92308 |
| 0.6 | 0.81757 / 0.81757 / 0.81757 | 0.75676 / 0.75676 |

**CNOT-noise threshold** (two-qubit depol `ε₂` after each CNOT, single-qubit gates
ideal — same convention as the controlled-SWAP study):

| input `ε` | dest `ε₂*` (2 CNOT) | cSWAP `ε₂*` (16 CNOT) | ratio |
|-----------|:------------------:|:---------------------:|:-----:|
| 0.10 | 0.1456 | 0.0330 | 4.41 |
| 0.20 | 0.2285 | 0.0612 | 3.73 |
| 0.30 | 0.2802 | 0.0845 | 3.32 |
| 0.40 | 0.3132 | 0.1029 | 3.04 |
| 0.50 | 0.3333 | 0.1167 | 2.86 |
| 0.60 | 0.3435 | 0.1257 | 2.73 |

The destructive gadget tolerates **~3–4× higher per-CNOT noise**, because only
2 CNOTs carry noise instead of 16. Figure: `destructive_gadget.png`.

**Caveat.** The destructive gadget is **measurement-only**: it yields `⟨O⟩` but not
a coherent purified state to feed forward. It is the correct tool for the
observable/threshold analysis (what this project computes), but not for
interleaving a purified state inside an algorithm — that still needs the
controlled-SWAP gadget.

---

## Summary

| | CNOTs | keeps | best `ε₂*` @ `ε=0.4` |
|--|:-----:|-------|:--------------------:|
| textbook controlled-SWAP | 16 | coherent purified state | 0.103 |
| optimized controlled-SWAP (Part 1) | 14 | coherent purified state | (≳0.103) |
| destructive / VD (Part 2) | **2** | observable `⟨O⟩` only | **0.313** |

For the operational-threshold study — which needs only `⟨O⟩` — the destructive/VD
gadget is both the largest CNOT reduction and the most noise-tolerant, and it
matches the virtual-distillation framing of the original PQEC paper.

Optional dependencies for Part 1: `qiskit`, `pytket` (`resynthesize_gadget.py`).
