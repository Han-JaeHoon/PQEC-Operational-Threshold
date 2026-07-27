# Step 5 — Variational (PQC) approximation of the gadget

Step 4 reduced the gadget's CNOT count in two *exact* senses (4a: same unitary,
16→14; 4b: same observable, 2 CNOTs). Step 5 asks the *approximate / learned*
version of the same question: **can a parameterized quantum circuit (PQC) be trained
to play the gadget's role with fewer CNOTs?** As in Step 4, "role" has a ladder of
meanings, and we do all three:

| | target the PQC reproduces | keeps | Step-4 analogue |
|--|--|--|--|
| **5a** | the full 5-qubit **unitary** `U = H_a·CSWAP·CSWAP·H_a` | everything | 4a |
| **5b** | `U` on the **physical input subspace** (ancilla `|0⟩`) — the coherent purified **state** | the fed-forward state | (between 4a/4b) |
| **5c** | only the purified **observable** `F = Tr(Oρ²)/Tr(ρ²)`, **noise-aware** | the threshold | 4b |

All conventions match the rest of the project: `ρ_ε = (1−ε)|Φ⁺⟩⟨Φ⁺| + ε I/4`
(`t = 1−ε`), kept register `A` = wires 1,2, observable `O = |Φ⁺⟩⟨Φ⁺|`, and a 2-qubit
global depolarizing channel `(1−ε₂)ρ + ε₂ I/4` after each CNOT (`s = 1−ε₂`).

Ansatz: a hardware-efficient PQC with a **CNOT budget `B`** as the knob
(`pqc_common.ansatz_ops`): a full `Rot` layer, then for each of `B` CNOTs (cycling an
ancilla-centric connectivity) a `Rot` on its two wires, then a final `Rot` layer
(`30 + 6B` parameters). The same ansatz is executed two ways — a fast pure-numpy
unitary builder (5a/5b) and a PennyLane `default.mixed` circuit with per-CNOT
depolarizing (5c) — verified to implement the identical map to machine precision
(`pqc_common._selftest`).

---

## The trainability barrier (why 5a and 5b are hard)

The natural cost for 5a is the **Hilbert–Schmidt test**
`C = 1 − |Tr(U†V)|²/d²` (`d = 32`); for 5b the same with `U→U₀` (the `32×16`
ancilla-`|0⟩` block of `U`) and `d→16`. Both are **global** cost functions and
suffer **barren plateaus** (McClean et al. 2018; Cerezo et al. 2021): random-init
gradients are tiny and the landscape is riddled with poor minima.

We differentiate the cost **exactly** (analytic backprop through the gate product,
`pqc_common._cost_grad`, matched to finite differences at `1e-11`) and run L-BFGS
with many restarts. We also implemented the **local Hilbert–Schmidt test (LHST)** of
Khatri et al. (*Quantum* **3**, 140, 2019), the standard plateau mitigation. Its
closed operator form (derived below) is

```
C_LHST = 3/4 − (1/(4 d n)) Re Tr(T M†),   M = V†U,   T = Σ_{j=1}^n Σ_{P∈{X,Y,Z}} P_j M P_j,
```

with `n = 5`, `d = 32`; `C_LHST = 0` iff `V = U` up to global phase, and its gradient
does **not** vanish exponentially.

**Derivation of the LHST form.** The Choi vector of `W = V†U` is
`|w⟩ = (1/√d) Σ_i (W|i⟩)_A ⊗ |i⟩_B`, i.e. amplitude matrix `W/√d`. The global HST
fidelity is the overlap of `|w⟩` with the product-of-Bell-pairs state
`⊗_j |Φ⁺⟩_{A_jB_j}`; LHST replaces it by the **average per-qubit** Bell overlap
`F_j = ⟨w| (|Φ⁺⟩⟨Φ⁺|_{A_jB_j} ⊗ I) |w⟩`. Using
`|Φ⁺⟩⟨Φ⁺| = ¼(II + XX − YY + ZZ)` and
`⟨w|(P_A⊗Q_B)|w⟩ = (1/d) Tr(P W Qᵀ W†)` (direct from the amplitude matrix), together
with `Xᵀ=X, Zᵀ=Z, Yᵀ=−Y`, the four terms give
`F_j = ¼ + (1/4d) Σ_{P∈{X,Y,Z}} Tr(P_j W P_j W†)`. Averaging over `j` and writing
`f_{P,j} = Tr(P_j W P_j W†) = Tr(P_j M P_j M†)` collapses the sum to
`Σ_{j,P} f_{P,j} = Tr(T M†)` with `T = Σ_{j,P} P_j M P_j`, giving the boxed `C_LHST`.
Its gradient uses the same prefix/suffix backprop as the global cost, with the single
extra object `T` (details in `pqc_common.cost_grad_lhst`, checked vs finite
differences at `1e-11`).

**Result: the barrier persists.** Best-of-16-restarts infidelity `δ = 1 − fidelity`
(so `δ = 0` is exact; `B = 16` CNOTs is already enough to represent `U` exactly, since
the textbook decomposition uses 16):

| `B` | 5a δ (full unitary) | 5b δ (coherent state) | 5a δ (LHST local cost) |
|----:|:-------------------:|:---------------------:|:----------------------:|
| 6 | 0.86 | 0.72 | 0.99 |
| 8 | 0.86 | 0.78 | 0.91 |
| 10 | 0.85 | 0.57 | 0.91 |
| 12 | 0.61 | 0.38 | 0.61 |
| 14 | 0.54 | 0.35 | 0.47 |
| 16 | **0.44** | **0.29** | 0.59 |

Even at `B = 16` the best infidelity is `0.44` (5a) / `0.29` (5b) — nowhere near `0`.
The LHST local cost does not rescue the global fidelity (its optimum trades global
faithfulness for local terms). Pushing further — a **richer ansatz** with a full
rotation layer after every CNOT, up to **375 parameters / 24 CNOTs** — still leaves
`δ ≈ 0.52`. So this is a genuine trainability/expressibility barrier, not an
under-optimization artifact. Figure: `pqc_compile_pareto.png`.

Conclusion: from-scratch variational compiling does **not** reproduce this gadget at
practical CNOT counts and does **not** beat the structured Step-4a decomposition.
Reproducing the coherent **state** (5b) is consistently **easier** than the full
**unitary** (5a) — lower infidelity at every budget — because it need only match the
16 physical columns; but it is still barrier-limited. The difficulty ordering
unitary (5a) > state (5b) > observable (5c) mirrors Step 4 exactly and points to
targeting the observable.

---

## Step 5c — noise-aware training of the observable (the useful route)

Relaxing all the way to the operational scalar `F(ε)` for `O ∈ {|Φ⁺⟩⟨Φ⁺|, ZZ}` gives
a **low-dimensional** target (a few numbers per `ε`) that trains easily even at a
**small** CNOT budget — no barrier. Two observables are used in the loss so the
circuit is forced to reproduce genuine purification behaviour, not curve-fit one
number.

The question 5c exists to answer: at a **fixed** ansatz with a **fixed** CNOT count,
does training that **includes** the CNOT depolarizing channel (noise-aware) give a
**higher operational threshold** than training noise-free and deploying? This
isolates *noise-aware compilation* from the trivial *fewer-CNOTs-less-noise* effect.

We compare, on `O = |Φ⁺⟩⟨Φ⁺|` (threshold = where `F_PQC > F_bare = (1+3t)/4`):
`θ_free` (trained at `ε₂ = 0`), `θ_aware` (trained at `ε₂ = 0.25`, init from `θ_free`),
the 16-CNOT textbook controlled-SWAP, and the 2-CNOT Step-4b destructive gadget.

**Results.**

*Trainability.* Unlike 5a/5b, the observable target trains to `loss ≈ 1e-5` — `θ_free`
reproduces `F_exact(ε)` to `~0.001` at only `B = 6` CNOTs. So a legitimate,
observable-equivalent gadget **is** learnable at a small CNOT count (whereas the
unitary/state was not).

*Positive result — fewer CNOTs lift the threshold.* Deployed under per-CNOT
depolarizing, `θ_free` (6 CNOTs) beats the 16-CNOT textbook at every input noise and
approaches the exact 2-CNOT Step-4b ceiling:

| input `ε` | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 |
|-----------|------|------|------|------|------|------|
| textbook cSWAP (16) `ε₂*` | 0.033 | 0.061 | 0.085 | 0.103 | 0.117 | 0.126 |
| **PQC `θ_free` (6) `ε₂*`** | **0.051** | **0.100** | **0.148** | **0.194** | **0.236** | **0.273** |
| Step 4b destructive (2) `ε₂*` | 0.146 | 0.229 | 0.280 | 0.313 | 0.333 | 0.343 |

*Negative result — noise-aware training does not help.* Refining `θ_free` with the
CNOT depolarizing **in** the loss gives **no threshold gain** (indeed a small loss),
and this is robust to how we pose it:

- **matching the ratio `F`** is *denominator-degenerate* — the optimizer drives
  `⟨Z_a⟩ → 0` and reproduces any `F`, giving unphysical "gadgets" (`F > 1`);
- **matching the absolute correlators** `⟨Z_a⟩, ⟨Z_a⊗O⟩` is *misaligned with the
  ratio* — a better absolute fit can still lower `F` — and at strong noise its global
  minimum is an **input-independent collapse** (a constant `⟨Z_a⟩` that ignores the
  input purity), which wins on loss regardless of warm-start.

Physically this is expected: the CNOT channel is **unital** (depolarizing contracts
every Bloch direction uniformly), so the `F`-bias it injects **cannot** be removed by
re-choosing single-qubit rotations. Any conceivable gain is limited to bias
cancellation in the ratio `F = ⟨Z_aO⟩/⟨Z_a⟩`, and we observe none. The threshold win
is entirely **structural** (fewer CNOTs), exactly as in Step 4. Figure:
`pqc_noise_aware.png`.

---

## Files

| File | Description |
|--|--|
| [`pqc_common.py`](pqc_common.py) | shared: target `U`, ansatz, fast numpy unitary + exact gradients, LHST cost, PennyLane noisy executor, read-out, references (self-test) |
| [`pqc_compile.py`](pqc_compile.py) | Step 5a/5b: variational compiling sweep (global + LHST costs), Pareto figure `pqc_compile_pareto.png` |
| [`pqc_noise_aware.py`](pqc_noise_aware.py) | Step 5c: noise-aware observable training, threshold comparison, figure `pqc_noise_aware.png` |

## Takeaway

The PQC study reproduces, in a *learned* setting, the central lesson of Step 4:
**the observable is the right thing to target.** Trying to learn the unitary (5a) or
the coherent state (5b) runs into the barren-plateau barrier and buys nothing over
the structured decomposition; learning the **observable** (5c) is easy, works at a
small CNOT count, and — because it uses far fewer noisy CNOTs — lifts the operational
threshold well above the 16-CNOT gadget, approaching the exact 2-CNOT Step-4b ceiling.
Noise-aware training adds only a marginal further gain, as expected for **unital**
depolarizing noise (uncompensable by unitary re-parameterization; any gain is bias
cancellation in the ratio `F = ⟨Z_aO⟩/⟨Z_a⟩`).
