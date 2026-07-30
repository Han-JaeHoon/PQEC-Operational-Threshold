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

## 5a/5b with generic ansätze: an expressibility / trainability squeeze
<!-- (These generic ansätze fail; the gadget-matched ansatz that SUCCEEDS is in
     "Step 5a — solved with a gadget-matched ansatz" below.) -->

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

Conclusion (for **these** ansätze): the obstruction is ansatz-specific, **not** a claim
that the gadget is uncompilable — a better-designed ansatz *does* compile it (next
section). The empirical difficulty ordering unitary (5a) > isometry (5b) > observable
(5c) still mirrors Step 4.

---

## Step 5a — solved with a gadget-matched ansatz (`pqc_ring_ansatz.py`)

The squeeze above is an ansatz limitation. Fixing three things makes a **generic
RX-RY-RZ PQC compile `U` to machine precision**:

1. **connectivity matched to the gadget** — `GADGET_PAIRS = (0,1)(0,3)(1,3)(0,2)(0,4)(2,4)`
   (ancilla to both registers + the swap pairs). A linear chain `0-1-2-3-4` cannot
   express `U` even at 20 CNOTs (teacher–student: in-class solves to `1e-15`, `U` stays
   at `δ≈0.6`).
2. **a full RX,RY,RZ layer after *every* CNOT** (not just between CNOT blocks).
3. **enough depth**: with `ansatz_percnot(GADGET_PAIRS, L)` —

   | `L` | CNOTs | params | `δ(U)` |
   |----:|:-----:|:------:|:------:|
   | 1 | 6 | 105 | 0.86 |
   | 2 | 12 | 195 | 0.44 |
   | **3** | **18** | **285** | **≈3e-15 (exact)** |

At `L=3` the landscape is already in the barren-plateau onset (random in-class targets
start to fail), yet the structured target `U` is found in ~15% of restarts.

**Pruning to 14** (`pqc_ring_prune.py`). Greedily removing CNOTs from the exact 18-CNOT
solution — warm-started, keeping all rotation layers so parameters stay aligned —
reaches **14 CNOTs, still exact** (`δ=0`). This is the *same* count as the Step-4a
peephole optimum and the `2×7` per-Fredkin optimum. **13 is unreachable**: removing any
one of the 14 (with warm + 11 random restarts, `reach13.py`) leaves `δ ≥ 0.146`. So
learn-then-prune and peephole optimization independently converge to 14.

**The arrangement, not the count, sets the noise threshold** (`pqc_ring_threshold.py`).
With `ε₂` on each of the learned 14 CNOTs (single-qubit gates ideal), the learned
14-CNOT circuit tolerates **~1.5–1.9× more per-CNOT noise than Step 4a's 14-CNOT
circuit**, nearly reaching the 2-CNOT destructive gadget at high input noise — although
both are exact 14-CNOT realizations of the same `U`:

| input `ε` | textbook (16) | Step 4a (14) | learned (14) | dest (2) |
|-----------|:---:|:---:|:---:|:---:|
| 0.10 | 0.033 | 0.041 | **0.060** | 0.146 |
| 0.40 | 0.103 | 0.140 | **0.237** | 0.313 |
| 0.60 | 0.126 | 0.178 | **0.338** | 0.343 |

So **CNOT count alone does not set the operational threshold — the arrangement does**,
with ~2× of room at fixed count. The robustness is *emergent*: the circuit was trained
and pruned noise-free, so a noise-aware selection among 14-CNOT realizations could do
better still. (Caveats: one specific pruned solution; 14-as-floor is strong convergent
evidence within this ansatz family, not a universal minimality proof.)

---

## The relaxation ladder on one ansatz+prune footing (5a / 5b / 5c)

With the gadget-matched ansatz solved for 5a, we can run the **same** learn-then-prune
procedure for the weaker 5b (isometry) and 5c (observable) targets and read off the
minimum CNOT count each requires. This puts all three rungs of the Step-5 ladder on a
common footing — same ansatz, same greedy pruning — and answers *where* in the
unitary → state → observable relaxation the CNOT savings actually appear.

**5b — the ancilla-`|0⟩` isometry `U₀` (`pqc_ring_5b.py`, `pqc_ring_5b_from14.py`).**
The isometry cost is `1 − |Tr(U₀†V₀)|²/16²` on the `32×16` block. Compiling from
scratch is slightly easier than the full unitary at intermediate depth (`δ_iso ≈ 0.24`
vs `0.44` at `L=2`) but still needs `L=3` (18 CNOTs) for exact compilation. Pruning is
**path-dependent**: greedily pruning a *separately* isometry-trained 18-CNOT solution
gets stuck at 17 (an artifact of that basin), so the honest test starts from the
**14-CNOT full-`U` solution** — which already compiles the isometry exactly
(`δ_iso = 0`) — and prunes with the isometry cost. It **cannot go below 14** (best
13-CNOT try `δ_iso = 0.13`). So `min-CNOT(isometry) = 14 = min-CNOT(full unitary)`:
**relaxing from the full unitary to the coherent state / isometry saves no CNOTs.**

**5c — the purified observable via the ancilla-parity read-out (`pqc_ring_5c.py`).**
Here we keep only the operational read-out `F = ⟨Z_a⊗O⟩/⟨Z_a⟩`. We prune the 14-CNOT
circuit under an **expectation-value cost** that matches the anchored correlators
`⟨Z_a⟩ → Tr(ρ²)`, `⟨Z_a⊗O⟩ → Tr(Oρ²)` for `O ∈ {|Φ⁺⟩⟨Φ⁺|, ZZ}` over an `ε` grid
(anchoring the denominator avoids the `⟨Z_a⟩→0` degeneracy). The cost is differentiated
exactly (`obs_cost_grad`, `Bracket = Σ_t 2(e_t−τ_t)·ρ_t V† M_t`, matched to finite
differences at `5.6e-10`; the 14-CNOT start has obs cost `1.4e-30`). Greedy pruning
then peels the circuit all the way down:

```
14 → 13 → 12 → 11 → 10 → 9 → 8 → 7 → 6 → 5   (every step obs cost ~1e-11)
```

reaching a **5-CNOT observable floor** (remaining pairs `(0,1)(2,4)(0,4)(2,4)(0,3)`);
4 is unreachable (best try `δ = 8.4e-2`). So relaxing to the observable takes the count
from 14 down to 5 on the *same* ancilla-parity architecture — and the *structured*
destructive read-out (Step 4b, no ancilla) reaches **2**.

**The ladder (all on the gadget-matched ansatz + greedy prune, `O ∈ {Φ⁺, ZZ}`):**

| rung | target the circuit reproduces | min CNOTs | Step-4 analogue |
|--|--|:--:|:--:|
| **5a** | full 5-qubit unitary `U` | **14** (13 impossible) | 4a (14) |
| **5b** | ancilla-`|0⟩` isometry `U₀` | **14** (same as 5a) | between 4a/4b |
| **5c** | purified observable `F`, ancilla-parity read-out | **5** | 4b |
| — | purified observable `F`, *structured* destructive read-out | **2** (Step 4b, exact) | 4b |

**Reading.** The CNOT savings are concentrated entirely at the **observable
relaxation**. Going unitary → coherent-state (5a → 5b) buys nothing (14 → 14); only
discarding the coherent output and keeping just `F` (5b → 5c) collapses the count
(14 → 5), and switching from an ancilla-parity to a structured destructive read-out
takes it further (5 → 2). This mirrors Step 4 exactly (4a keeps the unitary at 14; 4b
keeps only the observable at 2) — the learned setting reproduces the same lesson:
**relax the requirement to the observable, not to the state.** The 5-vs-2 gap is
architectural (a generic ancilla-parity gadget vs the purpose-built virtual-distillation
measurement), not a limit of the training. (Caveats as before: specific pruned solutions
on one ansatz family; floors are strong convergent evidence, not universal minimality
proofs.)

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
| [`pqc_ring_ansatz.py`](pqc_ring_ansatz.py) | Step 5a **success**: RX-RY-RZ ansätze; the gadget-matched, per-CNOT-rotation ansatz compiling `U` exactly at 18 CNOTs; saves `pqc_ring_L3_params.npy` |
| [`pqc_ring_prune.py`](pqc_ring_prune.py) | greedy CNOT pruning of the 18-CNOT solution → 14 (exact); saves `pqc_ring_pruned_*.{npy,json}` |
| [`reach13.py`](reach13.py) | rigorous test that 13 CNOTs is unreachable (14 is the floor here) |
| [`pqc_ring_threshold.py`](pqc_ring_threshold.py) | CNOT-noise threshold of the learned 14-CNOT gadget vs Step 4a/textbook/4b (`pqc_ring_threshold.png`) |
| [`pqc_ring_5b.py`](pqc_ring_5b.py), [`pqc_ring_5b_from14.py`](pqc_ring_5b_from14.py) | 5b rung: isometry compile + prune (from-scratch and from the 14-CNOT full-`U` solution) → isometry floor **14** |
| [`pqc_ring_5c.py`](pqc_ring_5c.py) | 5c rung: prune the 14-CNOT circuit under the ancilla-parity **observable** cost (exact expectation gradient) → observable floor **5**; saves `pqc_ring_5c_{params.npy,.json}` |

## Takeaway

Two complementary lessons.

**(1) The gadget *can* be variationally compiled — the ansatz is what matters.** A
generic RX-RY-RZ PQC reproduces `U` to machine precision at **18 CNOTs** once it has
gadget-matched connectivity and a rotation layer after every CNOT; greedy pruning then
reaches **14 CNOTs** (= the Step-4a peephole optimum; 13 is unreachable). Poorly-designed
ansätze (linear chain, sparse rotations) fail — to an *expressibility/trainability
squeeze* diagnosed by the teacher–student test — but that is an ansatz limitation, not a
property of the gadget. And strikingly, **CNOT count alone does not fix the noise
threshold**: the learned 14-CNOT arrangement is ~1.5–1.9× more noise-robust than the
Qiskit 14-CNOT one, nearly matching the 2-CNOT destructive gadget — so *how* the CNOTs
are arranged matters as much as *how many*.

**(2) For the observable, target the observable (5c).** Learning only `F(ε)` trains
easily at a small CNOT count and, using far fewer noisy CNOTs, lifts the threshold well
above the 16-CNOT gadget toward the exact 2-CNOT Step-4b **reference** — but the learned
circuit is a **specialized estimator** (Bell-isotropic inputs, `{Φ⁺, ZZ}`), not a general
gadget, and noise-aware training added no further gain for the objectives tried (an
empirical negative consistent with unital noise, not a proof of impossibility).

**(3) The CNOT savings live at the observable relaxation, not the state relaxation.**
Running the same learn-then-prune down the whole ladder gives min-CNOT counts **5a = 14,
5b = 14, 5c = 5** (ancilla-parity), with the structured Step-4b destructive read-out at
**2**. Relaxing unitary → coherent-state costs nothing (14 → 14); only relaxing to the
observable collapses the count (14 → 5 → 2) — precisely the Step-4 lesson, now
reproduced by training and pruning on a single ansatz.
