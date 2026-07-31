# Reducing the CNOT count of the SWAP-test gadget

The one-round gadget `H_a · CSWAP(a;A1,B1) · CSWAP(a;A2,B2) · H_a` (ancilla `a` =
wire 0, retained register `A=(A1,A2)` = wires 1,2, discarded `B=(B1,B2)` = wires
3,4) costs **16 CNOTs** with the textbook 8-CNOT Fredkin decomposition. We reduce
this in two senses of "same role":

- **Step 4 — Same unitary** — keep the exact 5-qubit unitary, fewer CNOTs.
- **Auxiliary — Same measured observable (destructive gadget)** — reproduce only
  `F = Tr(Oρ²)/Tr(ρ²)`, which is all PQEC uses; this allows a much larger reduction.
  (Not part of the main linear flow; kept as a low-CNOT reference baseline.)

All results are checked rigorously (unitary equivalence for Step 4; machine-precision
observable equivalence for the destructive gadget).

---

## Step 4 — Unitary-preserving reduction: 16 → 14 CNOTs

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

Both peephole optimizers **found** 14 CNOTs, which equals **2 × 7** — the known
per-Fredkin optimum (cf. the 7-CNOT Fredkin, Cruz–Murta arXiv:2305.18128). PyZX
minimizes T-count, not CNOT-count, so it is not the right tool for this metric.

**Scope of the claim.** "16 → 14" means *these optimizers produced a verified
14-CNOT circuit with the same unitary*. It is **not** a proof that 14 is the
minimum: minimality of the CNOT count over the full 5-qubit controlled
register-SWAP unitary would need a lower-bound argument or exhaustive synthesis,
which we have not done. The shared ancilla control gave no *further* reduction with
these tools.

**Takeaway:** a modest, safe 12.5% CNOT reduction is available for free while
keeping the exact coherent gadget.

**CNOT-noise threshold (computed, not assumed).** The 14-CNOT circuit is re-implemented
*identically* in PennyLane (`pqec_resynth_noise.py`, hard-coded gate list; Hilbert–
Schmidt overlap with the gadget `U` = `1.0000000000` up to global phase; `ε₂=0` read-out
matches `Tr(Oρ²)/Tr(ρ²)` for `Φ⁺` and `ZZ` to `1e-16`). Putting a 2-qubit depolarizing
`ε₂` after each of the 14 CNOTs (single-qubit gates ideal — same convention as Step 3 / the destructive gadget)
gives, for the purified fidelity `F` with `O=|Φ⁺⟩⟨Φ⁺|`:

| input `ε` | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 |
|-----------|------|------|------|------|------|------|
| textbook cSWAP (16) `ε₂*` | 0.033 | 0.061 | 0.085 | 0.103 | 0.117 | 0.126 |
| **Step 4 optimized (14) `ε₂*`** | **0.041** | **0.079** | **0.112** | **0.140** | **0.162** | **0.178** |
| ratio 14/16 | 1.25 | 1.29 | 1.32 | 1.36 | 1.39 | 1.42 |

So the exact-unitary reduction lifts the threshold **~1.25–1.42×** — *more* than the naive
`16/14 = 1.14` count ratio, i.e. the machine-found 14-CNOT layout also propagates noise a
little more favourably. It stays well below the measurement-only destructive gadget (2 CNOTs), as
expected. Figure: `pqec_resynth_threshold.png`.

---

## Auxiliary — Measurement-equivalent (destructive) gadget: 2 CNOTs

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

**Measurement cost — an important clarification.** The denominator `O_den = Π` is
**diagonal** in the post-`V` computational basis, so `Tr(ρ²)` comes from a *single*
Bell measurement. The numerator `O_num` is **not** diagonal there — its Pauli
expansion has **40 nonzero strings for `O=|Φ⁺⟩⟨Φ⁺|`** and **8 for `O=ZZ`** — so a
general numerator needs *several Pauli measurement settings* (each a 2-CNOT circuit
+ single-qubit rotations, combined classically). Hence **"2 CNOTs" is the per-setting
two-qubit-gate cost**, not a claim that one circuit yields every correlator. The
*mean* value `F` (and therefore the threshold) is unaffected by this — it only
changes the shot/setting budget. So this is best read as an **alternative
destructive measurement of the same VD estimator**, not a compilation of the
controlled-SWAP unitary.

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

The destructive gadget tolerates **~2.7–4.4× higher per-CNOT noise** (this is a
**mean-fidelity threshold** at fixed per-CNOT `ε₂`, not a total-resource claim —
the full sampling cost also depends on the number of Pauli settings for the
numerator). Figure: `destructive_gadget.png`.

**Exact closed form** (single-qubit gates ideal, `s=1−ε₂`, `t=1−ε`; verified vs
circuit to `~1e-16`):

```
F_dest(t, ε₂) = (1 + 6 s² t + 9 s² t²) / (4 (1 + 3 s² t²))       [F_dest(t,0) = (1+3t)²/(4(1+3t²))]
threshold     ε₂* = q_th(t) = 1 − 1/√(2 + 2t − 3t²)
```

**Caveat.** The destructive gadget is **measurement-only**: it yields `⟨O⟩` but not
a coherent purified state to feed forward. It is the correct tool for the
observable/threshold analysis (what this project computes), but not for
interleaving a purified state inside an algorithm — that still needs the
controlled-SWAP gadget.

---

## Summary

| | CNOTs / measurement circuit | keeps | `ε₂*` @ `ε=0.4` |
|--|:---------------------------:|-------|:---------------:|
| textbook controlled-SWAP | 16 | coherent purified state | 0.103 |
| optimized controlled-SWAP (Step 4) | 14 | coherent purified state | **0.140** (computed) |
| destructive / VD (auxiliary) | **2** (per setting) | observable `⟨O⟩` only | **0.313** |

The per-circuit two-qubit-gate count drops **7–8×** (16→2, or 14→2), and the
destructive/VD gadget has a **2.7–4.4× higher mean-fidelity threshold**. For an
M-qubit register the Bell-change costs `M` CNOTs vs `8M` (textbook) / `7M`
(optimized) / for the controlled-SWAP.

Caveats to keep the claim honest: (i) 14 is what the optimizers found, not a proven
minimum; (ii) fewer CNOTs *usually* helps but the threshold also depends on CNOT
topology and noise propagation (recomputed here, not assumed); (iii) the numerator
needs several Pauli measurement settings, so "2 CNOTs" is a per-setting cost and the
*total* sampling cost (shots × settings) is a separate question from the mean-bias
threshold. The destructive route is an alternative destructive **measurement** of
the same virtual-distillation estimator — matching the VD framing of the original
PQEC paper — not a re-compilation of the controlled-SWAP unitary.

Optional dependencies for Step 4: `qiskit`, `pytket` (`resynthesize_gadget.py`).
