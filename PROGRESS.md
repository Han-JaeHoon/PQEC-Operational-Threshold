# Progress log

Dated log of what has been built in the PQEC operational-threshold project.

## Project goal

Go beyond the ideal-gadget threshold analysis of Raghoonanan & Byrnes
(arXiv:2603.11568) — which applies noise to the data and then a *perfect*
purification — and find the threshold on the **noise of the PQEC operations
themselves** (controlled-SWAP etc.) below which purification still **recovers
entanglement**. Built up step by step, starting from the noisy input state.

---

## 2026-07-21 — Step 1: noisy input state `ρ_ε` + verification

The noisy input for the whole study is the isotropic (global-depolarizing) Bell
state

```
ρ_ε = (1 − ε) |Φ⁺⟩⟨Φ⁺| + ε · I/4.
```

**Built.**

- **`noisy_bell_state.py`** — prepares `ρ_ε` with a genuine PennyLane
  mixed-state circuit: `H · CNOT` to make `|Φ⁺⟩`, then a global 2-qubit
  depolarizing channel of strength `ε` implemented as one `QubitChannel` with
  the 16 two-qubit-Pauli Kraus operators
  (`K₀ = √(1−ε+ε/16) I`, `Kᵢ = √(ε/16) Pᵢ`). Includes `verify(eps)`.
- **`draw_noisy_bell.py`** — circuit diagram (`circuit_noisy_bell.png`).
- **`requirements.txt`**, **`.gitignore`** — environment and hygiene.

**Verified.** `verify(eps)` checks, for the circuit output vs the analytic
target:

| Check | Criterion |
|-------|-----------|
| analytic match | `max |ρ − ρ_ε| < 1e-12` |
| valid state | `Tr ρ = 1`, Hermitian, positive-semidefinite |
| Bell spectrum | `|Φ⁺⟩ → 1 − 3ε/4`, other three Bell states → `ε/4` |
| fidelity | `F = 1 − 3ε/4` |
| purity | `(1 − 3ε/4)² + 3(ε/4)²` |

**Result.** All checks pass for `ε ∈ {0, 0.1, 0.2, ⅓, 0.5, ⅔, 0.8, 1.0}` and for
500 random `ε ∈ [0,1]`, with `max |ρ − ρ_ε| ≈ 2.8e-16` — the circuit reproduces
the target to machine precision for arbitrary `ε`.

**Note.** `ρ_ε` is entangled iff `F > 1/2`, i.e. `ε < 2/3`; this is the bare-input
entanglement boundary that later gadget-noise thresholds will be compared against.

### Next

- Step 2: implement the SWAP-test purification gadget on this input.

---

## 2026-07-21 — Step 2: ideal PQEC purification gadget

The PQEC primitive: two identical noisy copies `ρ ⊗ ρ` enter, an ancilla-controlled
SWAP test (for the 2-qubit register, two parallel Fredkin gates) is applied, and
reading the ancilla extracts `P(ρ) = ρ²/Tr[ρ²]`, concentrating weight on the
dominant eigenvector.

**Built.**

- **`pqec_gadget.py`** — genuine 5-wire gadget (ancilla + register A `[1,2]` +
  register B `[3,4]`). Two read-outs:
  - `purify_once(ρ)` / `purify_rounds(ρ, ℓ)` — state extraction
    `ρ² = |0⟩ block − |1⟩ block` (Eq. 9), returning `ρ²/Tr[ρ²]`;
  - `obs_purified(ρ, O)` — the paper's protocol
    `⟨O⟩ = ⟨Z⊗O⟩/⟨Z⊗I⟩ = Tr(Oρ²)/Tr(ρ²)`.
  Plus `concurrence` and the `ρ_ε` recovery demo.
- **`draw_pqec_gadget.py`** — gadget circuit (`circuit_pqec_gadget.png`).

**Verified** (`verify()`), all to `~1e-16`:

| Check | Criterion |
|-------|-----------|
| state extraction | `purify_once(ρ) == ρ²/Tr[ρ²]` on 500 random states |
| observable protocol | `⟨Z⊗O⟩/⟨Z⊗I⟩ == Tr(Oρ²)/Tr(ρ²)` on 200 random states |
| on `ρ_ε` | `purify_once(ρ_ε) == ρ_ε²/Tr` over 40 values of `ε` |

**Result.** On `ρ_ε`, `|Φ⁺⟩` is the strictly dominant eigenvector for every
`ε < 1` (`1−3ε/4 > ε/4`), so the ideal gadget drives `F, C → 1` for **all `ε < 1`**
— even re-entangling a **separable** input (`2/3 ≤ ε < 1`, where input concurrence
is 0). Only `ρ = I/4` at `ε = 1` is a fixed point. Recovery curves in
`pqec_gadget_recovery.png`.

Concretely (ideal gadget, ℓ=3 rounds): `ε=0.30: C 0.55→1.00`, `ε=0.50: 0.25→1.00`,
`ε=2/3: 0.00→1.00`, `ε=0.80: 0.00→0.98`.

**Note for Step 3.** The `obs_purified` (ancilla-parity) read-out is the handle for
the noisy-gadget analysis: inserting depolarizing after the Fredkins/H and a readout
bit-flip, then measuring `⟨O⟩`, is exactly how the operational threshold will be
probed.

### Next

- Step 3: add noise to the gadget operations (Fredkin/H depolarizing, ancilla
  readout) and locate the operational threshold on the gadget error rate.

---

## 2026-07-21 — Step 3 (v1): 3-qubit global depolarizing on the Fredkins

First noisy-gadget model: right after each Fredkin `cSWAP(0;1,3)` and
`cSWAP(0;2,4)`, apply a **3-qubit global depolarizing channel** of strength `g_F`
to the three qubits it touched (ancilla included):

```
G_gF(σ) = (1 − g_F) σ + g_F · (I₈/8) ⊗ Tr_S(σ).
```

Ancilla H's and readout left ideal. Purified value read out with the parity
correlator `⟨O⟩ = ⟨Z_a⊗O⟩/⟨Z_a⊗I⟩`.

**Built.**

- **`pqec_gadget_noise.py`** — 3-qubit global-depol Kraus (`64` three-qubit
  Paulis), the noisy gadget circuit (`_gadget_obs_noisy`, `_gadget_state_noisy`),
  `obs_pqec_noisy(eps, g_F)`, the parity-weighted `effective_state_noisy`, and
  `no_qec`.
- **`verify_analytic_global_depol.py`** — checks the analytic derivation below
  against the circuit.
- **`plot_global_depol_benchmark.py`** — `global_depol_benchmark.png`.

**Key result — this model self-mitigates: signal loss, no threshold.**
Both correlators scale by exactly `(1−g_F)²` and the factor cancels in the ratio, so
the purified fidelity is **independent of `g_F`** for `0 ≤ g_F < 1`:

```
F_PQEC(p, g_F) = (1+3α²)² / (4(1+3α⁴)) = F_ideal-PQEC(p),   α = 1−4p/3  (α² = 1−ε).
```

The mechanism (Heisenberg picture): the measured observables `X_a⊗Φ_A`, `X_a⊗I_A`
are traceless on each noisy 3-qubit subsystem (because `Tr X_a = 0`, and back-
propagating through the second Fredkin keeps the ancilla part off-diagonal), so the
global-depol adjoint just multiplies each by `s = 1−g_F`. The error branch fully
randomizes the ancilla and carries no parity signal.

**Verified against the analytic derivation** (`verify_analytic_global_depol.py`),
all to `~1e-13` or better:

| Analytic | circuit error |
|----------|---------------|
| `A = (1−g_F)²(1+3α²)²/16` (numerator) | `2e-15` |
| `B = (1−g_F)²(1+3α⁴)/4` (denominator) | `1e-15` |
| `F_PQEC = (1+3α²)²/(4(1+3α⁴))`, `g_F`-independent | `1e-13` |
| `F_bare = 1−3ε/4 = 1−2p+4p²/3` | `7e-16` |
| `ΔF = 3α²(1−α²)(1+3α²)/(4(1+3α⁴)) > 0` for `0<p<3/4` | `6e-16` |
| sampling overhead `(B₀/B_g)² = (1−g_F)^{-4}` | exact (`1.52/4.16/16/10⁴`) |

Both parametrizations (global `ε`, local `p` with `α²=1−ε`) verified.

**Consequence.** No finite `g_F` fidelity threshold here; the only cost is a
sampling-overhead divergence `N_samp ∼ (1−g_F)^{-4}` (figure
`global_depol_benchmark.png`). This is a useful analytic **benchmark**, but the
global channel models Fredkin noise too symmetrically (it removes the whole parity
signal in the error branch, biasing nothing). A real operational threshold needs a
model that attenuates numerator and denominator **asymmetrically**.

**Step-by-step cross-check with the hand derivation.** A note (Han, July 22;
ordering `(a,A1,A2,B1,B2)`) evolves the 5-qubit state through the round and gives a
closed form after every gate: `σ₀=|+⟩⟨+|⊗R`, `σ₁=½[[R,RS₁],[S₁R,S₁RS₁]]`,
`σ₁'=sσ₁+(1−s)I⊗⁵/32`, `σ₂=(s/2)[[R,RS_AB],[RS_AB,R]]+(1−s)I⊗⁵/32`,
`σ₂'=(s²/2)[…]+(1−s²)I⊗⁵/32`, `σ_out=(s²/2)[[R+RS_AB,0],[0,R−RS_AB]]+(1−s²)I⊗⁵/32`
(`s=1−g`). `verify_note_states.py` reproduces **every** intermediate state on the
genuine circuit to `~1e-16` across several `(ε,g)` (incl. `g=0`, `g=0.9`); the note
and the implementation agree at every step, and the circuit diagram is
`circuit_gadget_noise.png` (`draw_gadget_noise.py`).

### Next

- Step 3 (v2): **independent local depolarizing** on the three Fredkin qubits
  (the derivation's §18). There the ancilla attenuation cancels but the data Bell
  correlator keeps an extra factor, giving `F_PQEC = ¼ + β²(F_ideal − ¼)` and a
  **finite `g_F` threshold**. Implement + verify.

---

## 2026-07-21 — Step 3b: realistic gate noise on a DECOMPOSED controlled-SWAP

Instead of an abstract 3-qubit channel, decompose each Fredkin into native 1- and
2-qubit gates and put realistic depolarizing noise on **each native gate**.

**Representative decompositions (literature).**

- **Textbook:** `CSWAP(q;a,b) = CNOT(b→a) · Toffoli(q,a;b) · CNOT(b→a)`, and the
  Toffoli as the standard Clifford+T circuit (Nielsen & Chuang, Fig. 4.9): 6 CNOT
  + 2 H + 7 T/T†. So **one Fredkin = 8 CNOTs** + single-qubit gates.
- **2-qubit-gate optimum:** 5 two-qubit gates (Smolin & DiVincenzo, PRA 53, 2855
  (1996)).
- **Recent connectivity-aware low-CNOT:** arXiv:2305.18128 (APL Quantum 1, 016105,
  2024) — lowest CNOT counts under all-to-all / linear connectivity.

We simulate the textbook Clifford+T version (all 2-qubit gates are CNOTs, so one
`p2` applies uniformly).

**Noise model.** After every CNOT, a 2-qubit depolarizing channel of strength
`p2` (`(1−p2)ρ + p2 I₄/4` — this equals the analytic `ε₂`); after every 1-qubit
gate (H, T, T†), a 1-qubit depolarizing `p1`; default `p1 = p2/10`. Read out with
the parity correlator `⟨O⟩ = ⟨Z_a⊗O⟩/⟨Z_a⊗I⟩`.

**Built.** `pqec_decomposed_noise.py` — noisy native gates, the Clifford+T Toffoli,
the decomposed Fredkin (with an `orient` option, below), `obs_pqec_decomposed`,
`threshold_p2`, scan + figure. `p1=p2=0` reproduces the ideal gadget to `6.7e-16`
(the CNOT+Toffoli decomposition is exact).

**Key result — a finite operational threshold appears.** Unlike the 3-qubit global
depolarizing (Step 3a, no threshold), native-gate noise hits the data qubits
asymmetrically, so it **biases** `⟨O⟩` and a finite `p2*` exists. With `p1=p2/10`,
`⟨O⟩` at `ε=0.40` falls from `0.942` and crosses the no-QEC baseline `1−3ε/4=0.70`
around `p2* ≈ 0.09`; `p2*` grows with input noise. `p2* ≈ 0.05–0.12` sits far above
realistic hardware 2-qubit errors (`~10⁻³–10⁻²`), so one PQEC round comfortably
tolerates realistic gate noise. (`ℓ=1`; multi-round would tighten it.)

**Exact analytic result (verified to ~1e-14).** With `u=1−ε₁`, `v=1−ε₂`,
`t=(1−4p/3)²`, `C=2u⁴v⁶+u⁶v⁵+3u⁶v⁶`, `D=1+(1+2u⁴)v⁴t²`:

```
B = (u⁹v¹⁰/4) D,   A = (u⁹v¹⁰/16)[D + t(1+t)C],   F_dec = A/B = ¼[1 + t(1+t)C/D].
```

`verify_analytic_decomposed.py` reproduces `A`, `B`, `F_dec`, the ideal limit and
the slopes on the genuine circuit (worst `3.9e-14`).

**Orientation matters (new finding).** The Toffoli target leg carries the H/T/T†
single-qubit gates, so it absorbs most of the single-qubit noise. The same Fredkin
unitary can put that target on either swapped qubit:

| orientation | Toffoli target on | single-qubit slope `K₁` | `ε₁` threshold @ p=0.4 |
|-------------|-------------------|-------------------------|------------------------|
| `retain`    | retained register A | **5/2** (= 2.500) | 0.140 |
| `discard`   | discarded register B | **2** (= 2.000) | 0.191 |

- The **denominator `B` and the CNOT slope `K₂ = 17/8`** (and the `ε₂` threshold)
  are **orientation-independent** — CNOT noise is symmetric across the swap.
- The analytic formula above corresponds to `orient="retain"`.
- **`orient="discard"` is the better layout**: putting the noisy target leg on the
  register you throw away shields the kept register, giving the milder `K₁=2` and a
  higher `ε₁` threshold. A protocol should choose this orientation.

**Both closed forms** (same `D`; verified to `~1e-14`):

```
C_retain (u,v) = 2u⁴v⁶ + u⁶v⁵ + 3u⁶v⁶            (P_R(u)=2u⁴+3u⁶,  P_R'(1)=26)
C_discard(u,v) =  u⁶v⁵ + (1+u²+2u⁶+u⁸) v⁶         (P_D(u)=1+u²+2u⁶+u⁸, P_D'(1)=22)
```

t-dependent slopes `F_dec ≈ F_ideal − K₁ε₁ − K₂ε₂` (verified against the circuit):

```
K₁_retain(t)  = 4t(1+t)(3t²+2)/(1+3t²)²        → 5/2  at t=1
K₁_discard(t) =  t(1+t)(9t²+7)/(1+3t²)²        → 2    at t=1
K₂(t)         =  t(1+t)(33t²+35)/(4(1+3t²)²)   → 17/8 at t=1   (both orientations)
```

Small-noise advantage condition: `K₁ε₁ + K₂ε₂ < Δ₀ ≈ 2p`; e.g. retain
`(5/2)ε₁+(17/8)ε₂ < 2p`, discard `2ε₁+(17/8)ε₂ < 2p`.

**Consistency notes (conventions).**
- `ε₂ = p₂`; `ε₁ = 4p₁/3` (analytic replacement-depolarizing vs PennyLane 1-qubit
  `DepolarizingChannel(p₁)`, contraction `u=1−4p₁/3`). Asserted in the verifier.
- Input `t`: `t = 1−ε` for the global Bell-depolarizing input `ρ_ε` (main study),
  or `t = (1−4p/3)²` for local depolarizing `p` per Bell qubit (verifier) — same
  isotropic family.
- Outer ancilla Hadamards: the verifier keeps them **ideal**, so its raw `A,B`
  carry `u⁹v¹⁰`; the main circuit also depolarizes the two outer H's, multiplying
  both `A` and `B` by a common `u²` (so `u¹¹v¹⁰`) that cancels in `F_dec`. Verified.

Figure: `pqec_decomposed_threshold.png` — (a) `⟨O⟩` vs 1q noise (orientation-split,
solid) vs 2q noise (coincident, dashed); (b) `p2*` vs input noise for both
orientations.

**Cross-check history.** An independent analytic derivation (GPT) first gave the
`retain` formula; an early review flagged its numerator as an error, but the
discrepancy was entirely the **orientation choice** — flipping the circuit
reproduces that formula to `1e-15`, and the `discard` closed form `C_discard` was
then confirmed the same way. Both analyses are correct for their respective layouts.

### Next

- entanglement (concurrence/negativity) threshold from the effective state, and
  multi-round accumulation of gate noise;
- coherent / biased-CNOT Fredkin error (Cruz–Murta) and the 7-CNOT optimized
  decomposition.
