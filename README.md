# PQEC Operational Threshold

Studying the **operational error threshold** of Purification Quantum Error
Correction (PQEC) for **entanglement distillation**.

The paper this builds on

> J. Raghoonanan & T. Byrnes, *Quantum Error Correction by Purification*,
> arXiv:2603.11568 (2026)

analyzes PQEC by applying a noise channel to the **data** and then a **perfect**
purification step. Real PQEC hardware — above all the 3-qubit controlled-SWAP
(Fredkin) at the heart of the SWAP test — is itself noisy. The goal of this
project is to go beyond the ideal-gadget analysis and find the threshold on the
**noise of the PQEC operations themselves** below which purification still
**recovers entanglement**.

This repository is being built up step by step. It starts from the noisy input
state and the tooling to certify it, before adding the (noisy) purification gadget.

## Status

| Step | Item | State |
|------|------|-------|
| 1 | Noisy input state `ρ_ε` — genuine preparation circuit + verification | **done** |
| 2 | Purification (SWAP-test) gadget — ideal, verified on `ρ_ε` | **done** |
| 3 | Noise on the gadget operations; operational threshold | planned |

## The noisy input state

The noisy input is the isotropic (global-depolarizing) Bell state

```
ρ_ε = (1 − ε) |Φ⁺⟩⟨Φ⁺| + ε · I/4,      |Φ⁺⟩ = (|00⟩ + |11⟩)/√2.
```

It is prepared by a genuine mixed-state circuit — `H · CNOT` to build `|Φ⁺⟩`,
then a **global 2-qubit depolarizing channel** of strength `ε` (a single
`QubitChannel` with the 16 two-qubit-Pauli Kraus operators):

![noisy Bell prep circuit](circuit_noisy_bell.png)

Key closed forms (all checked by the verification code):

- fidelity `F = ⟨Φ⁺|ρ_ε|Φ⁺⟩ = 1 − 3ε/4`
- Bell-basis spectrum: `|Φ⁺⟩ → 1 − 3ε/4`, the three other Bell states → `ε/4` each
- purity `Tr(ρ_ε²) = (1 − 3ε/4)² + 3(ε/4)²`
- `ρ_ε` is entangled iff `F > 1/2`, i.e. `ε < 2/3`

## The purification gadget (Step 2)

The PQEC primitive is the **SWAP-test gadget**: two identical noisy copies
`ρ ⊗ ρ` enter, an ancilla-controlled SWAP (for the 2-qubit register, two parallel
Fredkin gates) is applied, and reading the ancilla extracts the purified component

```
P(ρ) = ρ² / Tr[ρ²]        (concentrates weight on the dominant eigenvector)
```

![PQEC SWAP-test gadget](circuit_pqec_gadget.png)

The gadget is implemented as a genuine 5-wire circuit with two equivalent
read-outs (both used when the gadget is made noisy in Step 3):

- **state extraction** (Eq. 9): `ρ² = (ancilla |0⟩ block) − (|1⟩ block)`, so
  `purify_once(ρ)` returns `ρ²/Tr[ρ²]`;
- **observable / parity correlator** (the paper's actual protocol):
  `⟨O⟩_purified = ⟨Z⊗O⟩ / ⟨Z⊗I⟩ = Tr(Oρ²)/Tr(ρ²)`.

Both are verified to machine precision (`~1e-16`) on 500 random states and on `ρ_ε`.

On the isotropic input `ρ_ε`, `|Φ⁺⟩` is the strictly dominant eigenvector for
**every `ε < 1`** (eigenvalues `1−3ε/4` vs `ε/4`), so the ideal gadget restores
fidelity and entanglement to 1 for all `ε < 1` — it even **re-entangles a
separable input** (`2/3 ≤ ε < 1`); only `ρ = I/4` at `ε = 1` is a fixed point.

![recovery over rounds](pqec_gadget_recovery.png)

## Files

| File | Description |
|------|-------------|
| [`noisy_bell_state.py`](noisy_bell_state.py) | Prepares `ρ_ε` with a genuine circuit and `verify(eps)` — checks the analytic match, unit trace, Hermiticity, positive-semidefiniteness, the Bell spectrum, fidelity and purity |
| [`draw_noisy_bell.py`](draw_noisy_bell.py) | Draws the preparation circuit (`circuit_noisy_bell.png`) |
| [`pqec_gadget.py`](pqec_gadget.py) | Ideal SWAP-test gadget: `purify_once` / `purify_rounds` / `obs_purified`, verification, and the `ρ_ε` recovery demo |
| [`draw_pqec_gadget.py`](draw_pqec_gadget.py) | Draws the 5-wire gadget (`circuit_pqec_gadget.png`) |
| [`requirements.txt`](requirements.txt) | Dependencies (pinned minimums + tested versions) |

## Setup & run

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt

python noisy_bell_state.py         # build ρ_ε and verify it (sweep + 500 random ε)
python draw_noisy_bell.py          # regenerate the input circuit diagram
python pqec_gadget.py              # ideal gadget: verify + ρ_ε recovery demo
python draw_pqec_gadget.py         # regenerate the gadget circuit diagram
```

### Verification output (excerpt)

```
eps = 0.200 | F = 0.8500 (=1-3eps/4=0.8500) | purity = 0.7300 | max|rho-target| = 1.67e-16 | PASS
eps = 0.500 | F = 0.6250 (=1-3eps/4=0.6250) | purity = 0.4375 | max|rho-target| = 5.55e-17 | PASS
...
500 random eps in [0,1]: max |rho - target| = 2.78e-16
ALL CHECKS PASSED
```

The circuit reproduces `ρ_ε` to machine precision (`~1e-16`) for arbitrary `ε ∈ [0,1]`.
