# Branch `step-1-ideal-pqec` — Step 1: Ideal one-round PQEC baseline

- **Step:** 1 (ideal, noiseless SWAP-test purification gadget).
- **Input state:** `ρ_ε` from the Setup step.
- **Circuit:** genuine 5-wire gadget `U = H_a·CSWAP(a;A1,B1)·CSWAP(a;A2,B2)·H_a`
  (wire 0 = ancilla, register A = wires 1,2 [kept], register B = wires 3,4).
  Two read-outs: state extraction `ρ²/Tr[ρ²]` (Eq. 9) and the ancilla-parity observable
  `⟨O⟩ = ⟨Z⊗O⟩/⟨Z⊗I⟩ = Tr(Oρ²)/Tr(ρ²)`.
- **Noise model:** none (ideal baseline).
- **Core script:** `pqec_gadget.py` (`purify_once`/`purify_rounds`/`obs_purified`,
  verification, `ρ_ε` recovery demo); `draw_pqec_gadget.py` (circuit diagram).
- **Dependencies:** `noisy_bell_state.py`, `requirements.txt`.
- **Data files:** none.
- **Outputs:** `circuit_pqec_gadget.png`, `pqec_gadget_recovery.png`; verification to `~1e-16`.
- **Result:** on `ρ_ε`, `|Φ⁺⟩` is the strictly dominant eigenvector for every `ε<1`, so the
  ideal gadget drives `F,C→1` for all `ε<1` (even re-entangling a separable input).
- **Relation to other steps:** the ideal reference that Steps 2–5 make noisy; the
  ancilla-parity read-out here is the handle used for every operational threshold.

## Reproduce

```bash
git switch step-1-ideal-pqec
pip install -r requirements.txt
python pqec_gadget.py          # verify + ρ_ε recovery demo
python draw_pqec_gadget.py     # circuit_pqec_gadget.png
```

## git worktree

```bash
git worktree add ../pqec-step1 step-1-ideal-pqec
```

> Every step branch is a full copy of `main`; this guide marks the *core* files for Step 1.
> The integrated/stable branch is `main`.
