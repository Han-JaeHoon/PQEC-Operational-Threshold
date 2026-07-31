# Branch `step-2-fredkin-noise` — Step 2: Fredkin-level global replacement depolarizing noise

- **Step:** 2 (first noisy-gadget model).
- **Input state:** `ρ_ε` from the Setup step.
- **Circuit:** the Step-1 gadget, with noise inserted right after each Fredkin (CSWAP).
- **Noise model:** a **3-qubit global depolarizing channel** of strength `g_F` on the three
  qubits each Fredkin touches (ancilla included):
  `G(σ) = (1−g_F)σ + g_F·(I₈/8)⊗Tr_S(σ)`. Ancilla H's and read-out ideal.
- **Core scripts:** `pqec_gadget_noise.py` (noisy gadget, `obs_pqec_noisy`, effective state);
  `verify_analytic_global_depol.py` (analytic formulas vs circuit, `~1e-13`);
  `verify_note_states.py` (step-by-step hand-derivation states vs circuit, `~1e-16`);
  `plot_global_depol_benchmark.py`; `draw_gadget_noise.py`.
- **Dependencies:** `noisy_bell_state.py`, `pqec_gadget.py`, `requirements.txt`.
- **Data files:** none.
- **Outputs:** `global_depol_benchmark.png`, `circuit_gadget_noise.png`.
- **Key result:** this model **self-mitigates** — both parity correlators scale by `(1−g_F)²`,
  which cancels in the ratio, so `F_PQEC` is **independent of `g_F`** (no threshold); the only
  cost is a sampling-overhead divergence `N_samp ∼ (1−g_F)^{-4}`.
- **Relation to other steps:** shows why a *symmetric* channel gives no threshold, motivating
  the **asymmetric** per-CNOT noise of Step 3.

## Reproduce

```bash
git switch step-2-fredkin-noise
pip install -r requirements.txt
python pqec_gadget_noise.py            # ⟨O⟩ vs g_F (self-mitigates)
python verify_analytic_global_depol.py # analytic vs circuit (~1e-13)
python verify_note_states.py           # note states vs circuit (~1e-16)
python plot_global_depol_benchmark.py  # global_depol_benchmark.png
python draw_gadget_noise.py            # circuit_gadget_noise.png
```

## git worktree

```bash
git worktree add ../pqec-step2 step-2-fredkin-noise
```

> Every step branch is a full copy of `main`; this guide marks the *core* files for Step 2.
> The integrated/stable branch is `main`.
