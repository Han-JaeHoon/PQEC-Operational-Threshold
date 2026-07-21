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
