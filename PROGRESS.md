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
