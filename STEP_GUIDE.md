# Branch `step-3-textbook-16cnot` — Step 3: Textbook 16-CNOT implementation with per-CNOT noise

- **Step:** 3 (first model with a finite operational threshold).
- **Input state:** `ρ_ε` from the Setup step (`t = 1−ε`).
- **Circuit:** each Fredkin decomposed into native gates (textbook Clifford+T Toffoli,
  8 CNOTs/Fredkin → **16 CNOTs** per round).
- **Noise model:** a 2-qubit **replacement depolarizing** channel `(1−ε₂)ρ + ε₂ I/4` after
  **each CNOT**; single-qubit gates and ancilla Hadamards ideal (CNOT-only premise).
- **Core scripts:** `pqec_cnot_threshold.py` (closed forms + threshold `ε₂*` + figure);
  `verify_analytic_decomposed.py` (analytic `A/B/F_dec` and the `K₂` slope vs circuit,
  `~1e-14`); `draw_cnot_noise.py` (CSWAP decomposition, SWAP test, +CNOT-noise diagrams).
  Write-up: `CNOT_NOISE_ANALYSIS.md`.
- **Dependencies:** `noisy_bell_state.py`, `pqec_gadget.py`, `requirements.txt`.
- **Data files:** none.
- **Outputs:** `pqec_cnot_threshold.png`, `circuit_cswap_decomp.png`,
  `circuit_swaptest_decomp.png`, `circuit_swaptest_cnot_noise.png`.
- **Key result:** a finite threshold appears — `F_dec = ¼[1 + t(1+t)(v⁵+5v⁶)/(1+3v⁴t²)]`
  (`v=1−ε₂`); `ε₂*` = 0.033 / 0.061 / 0.085 / 0.103 / 0.117 / 0.126 at `ε` = 0.1…0.6.
- **Relation to other steps:** the **baseline 16-CNOT threshold** that Steps 4 (14-CNOT
  resynthesis) and 5 (learned 14-CNOT) are compared against, under the identical noise model.

## Reproduce

```bash
git switch step-3-textbook-16cnot
pip install -r requirements.txt
python verify_analytic_decomposed.py   # analytic A/B/F_dec + K2 slope vs circuit (~1e-14)
python pqec_cnot_threshold.py          # threshold ε₂* table + pqec_cnot_threshold.png
python draw_cnot_noise.py              # the three CNOT-only circuit diagrams
```

## git worktree

```bash
git worktree add ../pqec-step3 step-3-textbook-16cnot
```

> Every step branch is a full copy of `main`; this guide marks the *core* files for Step 3.
> The integrated/stable branch is `main`.
