# Branch `setup-bell-input` — Setup: Bell-isotropic input preparation & verification

This branch is the entry point of the study: prepare and certify the noisy input state
that every later step consumes.

- **Step:** Setup (Bell-isotropic input `ρ_ε`).
- **Input state:** `ρ_ε = (1−ε)|Φ⁺⟩⟨Φ⁺| + ε·I/4` (isotropic / global-depolarizing Bell
  state), prepared by `H·CNOT` then a global 2-qubit depolarizing channel of strength `ε`
  (one `QubitChannel` with the 16 two-qubit-Pauli Kraus operators).
- **Circuit / noise model:** no gadget yet; this step only builds and verifies the input.
- **Core script:** `noisy_bell_state.py` (prepares `ρ_ε`, `verify(eps)`); `draw_noisy_bell.py`
  (circuit diagram).
- **Dependencies:** `requirements.txt` (PennyLane, NumPy, SciPy, Matplotlib).
- **Data files:** none required.
- **Outputs:** `circuit_noisy_bell.png`; `verify(eps)` PASS for a sweep + 500 random `ε`
  (`max|ρ−ρ_ε| ≈ 1e-16`).
- **Relation to other steps:** provides `rho_eps_analytic`, `O_PHI_PLUS`, `PHI_PLUS`
  imported by every subsequent step (Steps 1–5). `ρ_ε` is entangled iff `ε < 2/3`.

## Reproduce

```bash
git switch setup-bell-input
pip install -r requirements.txt
python noisy_bell_state.py     # build ρ_ε and verify (sweep + 500 random ε)
python draw_noisy_bell.py      # circuit_noisy_bell.png
```

## Working on several steps at once (git worktree)

```bash
git worktree add ../pqec-setup setup-bell-input
git worktree add ../pqec-step3 step-3-textbook-16cnot
git worktree add ../pqec-step5 step-5-learned-14cnot
```

> Note: every step branch is a full copy of `main` (all files present) so each is
> independently checkout-able and runnable; this guide marks the files that are *core* to
> this step. The integrated/stable branch with all results and comparison code is `main`.
