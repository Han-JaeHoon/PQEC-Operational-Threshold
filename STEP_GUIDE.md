# Branch `step-4-resynthesized-14cnot` — Step 4: Resynthesized 14-CNOT implementation with per-CNOT noise

- **Step:** 4 (unitary-preserving CNOT reduction).
- **Input state:** `ρ_ε` from the Setup step (`t = 1−ε`).
- **Circuit:** the *same* 5-qubit gadget unitary `U`, resynthesized from 16 → **14 CNOTs**
  by two independent peephole optimizers (Qiskit `opt2/3` and pytket `FullPeepholeOptimise`),
  each verified to implement the identical unitary (up to global phase); `14 = 2×7`.
- **Noise model:** 2-qubit depolarizing `(1−ε₂)ρ + ε₂ I/4` after each of the 14 CNOTs;
  single-qubit gates ideal (same convention as Step 3).
- **Core scripts:** `resynthesize_gadget.py` (16→14 reduction, needs `qiskit`/`pytket`);
  `pqec_resynth_noise.py` (pure-PennyLane replay, unitary-verified, + CNOT-noise threshold);
  `draw_resynth_pl.py` (PennyLane circuit diagram). Write-up: `SWAP_GADGET_OPTIMIZATION.md`.
- **Dependencies:** `noisy_bell_state.py`, `pqc_common.py`, `pqec_cnot_threshold.py`,
  `verify_analytic_decomposed.py`, `destructive_gadget.py` (imported for the comparison
  baseline), `requirements.txt`; optional `qiskit`, `pytket` for `resynthesize_gadget.py`.
- **Data files:** none (the 14-CNOT gate list is hard-coded in `pqec_resynth_noise.py`).
- **Outputs:** `pqec_resynth_threshold.png`, `circuit_resynth_pennylane.png`.
- **Key result:** the 14-CNOT circuit keeps the exact unitary and its threshold is
  **1.25–1.42×** the 16-CNOT textbook (`ε₂*` = 0.041 / 0.140 / 0.178 at `ε` = 0.1 / 0.4 / 0.6).
- **Relation to other steps:** same job as Step 3 with fewer noisy CNOTs; the coherent
  14-CNOT baseline that Step 5's *learned* 14-CNOT arrangement is compared against. (The
  2-CNOT destructive gadget in `destructive_gadget.py` is an **auxiliary** reference, not a
  main-flow step.)

## Reproduce

```bash
git switch step-4-resynthesized-14cnot
pip install -r requirements.txt
python pqec_resynth_noise.py    # unitary check + threshold table + pqec_resynth_threshold.png
python draw_resynth_pl.py       # circuit_resynth_pennylane.png
python resynthesize_gadget.py   # (optional) re-derive 16->14 (needs qiskit/pytket)
```

## git worktree

```bash
git worktree add ../pqec-step4 step-4-resynthesized-14cnot
```

> Every step branch is a full copy of `main`; this guide marks the *core* files for Step 4.
> The integrated/stable branch is `main`.
