# Step 5 — Variational (PQC) approximation of the gadget

Step 4 reduced the gadget's CNOT count in two *exact* senses (4a: same unitary,
16→14; 4b: same observable, 2 CNOTs). Step 5 asks the *approximate / learned*
version of the same question: **can a parameterized quantum circuit (PQC) be trained
to play the gadget's role with fewer CNOTs?** As in Step 4, "role" has a ladder of
meanings, and we do all three:

| | target the PQC reproduces | keeps | Step-4 analogue |
|--|--|--|--|
| **5a** | the full 5-qubit **unitary** `U = H_a·CSWAP·CSWAP·H_a` | everything | 4a |
| **5b** | `U` on the **ancilla-`|0⟩` input subspace** — the `32×16` **isometry** `U₀ = U E₀` (every physical input, since `ρ_ε` is full rank) | the fed-forward state | (between 4a/4b) |
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

## Why 5a and 5b fail: an expressibility / trainability squeeze

The natural cost for 5a is the **Hilbert–Schmidt test**
`C = 1 − |Tr(U†V)|²/d²` (`d = 32`); for 5b the same with `U→U₀` (the `32×16`
ancilla-`|0⟩` block of `U`) and `d→16`. Both are **global** cost functions, which are
prone to **barren plateaus** (McClean et al. 2018; Cerezo et al. 2021). But — as a
teacher–student test below shows — the plateau is only *half* the story: at shallow
depth the obstruction is instead the ansatz's limited **expressibility**.

We differentiate the cost **exactly** (analytic backprop through the gate product,
`pqc_common._cost_grad`, matched to finite differences at `1e-11`) and run L-BFGS
with many restarts. We also implemented the **local Hilbert–Schmidt test (LHST)** of
Khatri et al. (*Quantum* **3**, 140, 2019), the standard plateau mitigation. Its
closed operator form (derived below) is

```
C_LHST = 3/4 − (1/(4 d n)) Re Tr(T M†),   M = V†U,   T = Σ_{j=1}^n Σ_{P∈{X,Y,Z}} P_j M P_j,
```

with `n = 5`, `d = 32`; `C_LHST = 0` iff `V = U` up to global phase. LHST is designed
to *mitigate* (under suitable locality/depth conditions) the exponential gradient
decay of the global cost — it is not a guarantee of convergence for an arbitrary
ansatz, and below we find it does not rescue the global fidelity here.

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

**Result 1 — no budget reaches `U`.** Best-of-16-restarts infidelity `δ = 1 − fidelity`
(`δ = 0` is exact):

| `B` | 5a δ (full unitary) | 5b δ (anc-`|0⟩` isometry) | 5a δ (LHST local cost) |
|----:|:-------------------:|:-------------------------:|:----------------------:|
| 6 | 0.86 | 0.72 | 0.99 |
| 8 | 0.86 | 0.78 | 0.91 |
| 10 | 0.85 | 0.57 | 0.91 |
| 12 | 0.61 | 0.38 | 0.61 |
| 14 | 0.54 | 0.35 | 0.47 |
| 16 | **0.44** | **0.29** | 0.59 |

Even at `B = 16` the best `δ` is `0.44` (5a) / `0.29` (5b). LHST does not rescue the
global fidelity, and a richer ansatz (full rotation layer after every CNOT, up to
**375 params / 24 CNOTs**) still leaves `δ ≈ 0.52`. 5b is consistently **easier** than
5a — it need only match the 16 physical columns. Figure: `pqc_compile_pareto.png`.

**Result 2 — the cause is an expressibility/trainability *squeeze*, not a plain
plateau.** A **teacher–student** test (`test_inclass.py`) isolates the mechanism: we
recompile a *reachable* target `V(θ*)` drawn from the **same** ansatz with the **same**
optimizer, so `δ → 0` is achievable by construction and any residual `δ` is purely the
optimizer's failure:

| `B` | in-class `δ` (reachable targets, 30 restarts) | reading |
|----:|-----------------------------------------------|---------|
| 12 | `4e-16, 1e-15, 9e-16` — **all succeed** | optimizer fully capable, yet `U` unreached (`δ=0.61`) → **expressibility**: the fixed topology cannot hold `U` at 12 CNOTs |
| 16 | `9e-15, 0.53, 5e-15, 0.54` — **2 of 4** | transition |
| 20 | `0.81, 0.82, 0.80` — **all fail** | optimizer fails on *provably reachable* targets → genuine **barren-plateau / trainability** barrier |

So the two obstructions **trade off with depth**: where the ansatz is **trainable**
(`B ≤ 12`) it is **not expressive enough** to contain `U`; by the depth where it might
be (`B ≥ 16`) the optimizer **can no longer navigate** the landscape — it fails even on
targets that are reachable by construction. No budget both expresses and reaches `U`.
Calling this simply a "barren plateau" (as an earlier draft did) was imprecise: at low
depth it is an **expressibility/topology** limit; only at higher depth does the
**trainability** barrier bind.

Conclusion: from-scratch variational compiling with this **fixed-topology** ansatz does
**not** reproduce the gadget and does **not** beat the structured Step-4a
decomposition — but the obstruction is this ansatz-specific squeeze, **not** a claim
that the gadget is uncompilable in general (a topology-matched ansatz built to contain
the known 14-CNOT circuit would trivially succeed — that *is* Step 4a). The empirical
difficulty ordering unitary (5a) > isometry (5b) > observable (5c) mirrors Step 4 and
points to targeting the observable.

---

## Step 5c — noise-aware training of the observable (the useful route)

Relaxing all the way to the operational scalar `F(ε)` for `O ∈ {|Φ⁺⟩⟨Φ⁺|, ZZ}` gives
a **low-dimensional** target (a few numbers per `ε`) that trains easily even at a
**small** CNOT budget — no barrier. Two observables plus the denominator anchor are
used in the loss (see below); this pins down `F` for those observables on the
Bell-isotropic input family, but — as the out-of-sample check shows — does **not**
make the circuit a general-purpose gadget.

The question 5c exists to answer: at a **fixed** ansatz with a **fixed** CNOT count,
does training that **includes** the CNOT depolarizing channel (noise-aware) give a
**higher operational threshold** than training noise-free and deploying? This
isolates *noise-aware compilation* from the trivial *fewer-CNOTs-less-noise* effect.

We compare, on `O = |Φ⁺⟩⟨Φ⁺|` (threshold = where `F_PQC > F_bare = (1+3t)/4`):
`θ_free` (trained at `ε₂ = 0`), `θ_aware` (trained at `ε₂ = 0.25`, init from `θ_free`),
the 16-CNOT textbook controlled-SWAP, and the 2-CNOT Step-4b destructive gadget.

**Results.**

*Trainability.* Unlike 5a/5b, the observable target trains to `loss ≈ 1e-5` — `θ_free`
reproduces `F_exact(ε)` to `~0.001` at only `B = 6` CNOTs.

*Scope — a specialized estimator, not a general gadget.* An out-of-sample check
(`test_5c_oos.py`) makes the claim precise. Over a **dense unseen `ε` grid** (40
points in `[0.02, 0.72]`, well outside the 4 training points) `θ_free` still matches
`F_exact` to `0.0014` for `Φ⁺` and `0.0021` for `ZZ` — so it is **not** `ε`-grid
curve-fitting; it genuinely generalizes across the input family. But on **unseen
observables** it is far off (`XX`, `YY`: `~0.13`; `ZI`: `0.05`), whereas Step 4b is
exact for *all* `O`. So `θ_free` is best described as **a learned estimator
specialized to the Bell-isotropic input family and the observables `{Φ⁺, ZZ}`**, not a
general observable-equivalent gadget. (The threshold uses `O = Φ⁺` only, so this
scoping does not affect the threshold numbers below; and `F(ε₂)` is verified
**monotone** decreasing, so the threshold crossing is well-defined.)

*Positive result — fewer CNOTs lift the threshold.* Deployed under per-CNOT
depolarizing, `θ_free` (6 CNOTs) beats the 16-CNOT textbook at every input noise and
approaches the exact 2-CNOT Step-4b reference:

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

This is consistent with the physics: the CNOT channel is **unital** (depolarizing
contracts every Bloch direction uniformly), so the noise cannot be *inverted* by
unitary re-parameterization — any gain could only come from redistributing which
Pauli paths traverse which noisy CNOTs, or from bias cancellation in the ratio
`F = ⟨Z_aO⟩/⟨Z_a⟩`. **We did not observe such a gain** for the tested ansatz,
objectives (ratio-matching and correlator-matching), and noise strengths. This is an
empirical negative, **not** a proof of impossibility: other objectives we did not try
— e.g. a denominator-constrained margin objective `max_θ min_{ε, q≤q₀} [F_θ(ε,q) −
F_bare(ε)]` — remain open. What we *can* say is that on this problem the threshold win
came entirely from the **structural** change (fewer noisy CNOTs), exactly as in
Step 4. Figure: `pqc_noise_aware.png`.

---

## Files

| File | Description |
|--|--|
| [`pqc_common.py`](pqc_common.py) | shared: target `U`, ansatz, fast numpy unitary + exact gradients, LHST cost, PennyLane noisy executor, read-out, references (self-test) |
| [`pqc_compile.py`](pqc_compile.py) | Step 5a/5b: variational compiling sweep (global + LHST costs), Pareto figure `pqc_compile_pareto.png` |
| [`pqc_noise_aware.py`](pqc_noise_aware.py) | Step 5c: noise-aware observable training, threshold comparison, figure `pqc_noise_aware.png` |
| [`test_inclass.py`](test_inclass.py) | teacher–student test isolating expressibility vs trainability (in-class recompilation at `B=12/16/20`) |
| [`test_5c_oos.py`](test_5c_oos.py) | Step-5c out-of-sample check: unseen observables, dense `ε` grid, `F(ε₂)` monotonicity |
| [`draw_pqc_ansatz.py`](draw_pqc_ansatz.py) | draws the ansatz structure (`circuit_pqc_ansatz.png`) |

## Takeaway

The PQC study reproduces, in a *learned* setting, the central lesson of Step 4:
**the observable is the right thing to target.** Learning the unitary (5a) or the
ancilla-`|0⟩` isometry (5b) with a generic fixed-topology ansatz fails — to an
**expressibility/trainability squeeze** (not a plain plateau), and tied to that ansatz
rather than to the gadget itself — and buys nothing over the structured decomposition.
Learning the **observable** (5c) trains easily at a small CNOT count and, because it
uses far fewer noisy CNOTs, lifts the operational threshold well above the 16-CNOT
gadget, approaching the exact 2-CNOT Step-4b **reference**; but the learned circuit is
a **specialized estimator** (Bell-isotropic inputs, `{Φ⁺, ZZ}`), not a general gadget.
Noise-aware training added **no** further gain for the objectives we tried — an
empirical negative consistent with the noise being unital, **not** a proof of
impossibility.
