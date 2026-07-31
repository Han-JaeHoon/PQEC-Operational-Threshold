# Branch `step-5-learned-14cnot` — Step 5: Learned & pruned 14-CNOT implementation with per-CNOT noise

- **Step:** 5 (variationally-compiled gadget).
- **Input state:** `ρ_ε` from the Setup step (`t = 1−ε`).
- **Circuit:** a gadget-matched RX-RY-RZ PQC (connectivity `(0,1)(0,3)(1,3)(0,2)(0,4)(2,4)`,
  a full RX-RY-RZ layer after every CNOT) trained to compile `U` exactly at 18 CNOTs, then
  greedily pruned to **14 CNOTs** (still exact; 13 unreachable).
- **Noise model:** 2-qubit depolarizing `(1−ε₂)ρ + ε₂ I/4` after each of the learned 14
  CNOTs; single-qubit gates ideal (same convention as Steps 3–4).
- **Core scripts:** `pqc_common.py` (target `U`, ansatz, fast numpy unitary + exact
  gradients, LHST cost, PennyLane noisy executor, read-out); `pqc_ring_ansatz.py`
  (gadget-matched ansatz, compiles `U`); `pqc_ring_prune.py` (prune 18→14);
  `reach13.py` (13 unreachable); `pqc_ring_threshold.py` (CNOT-noise threshold vs
  Step 3/4/destructive); `draw_pqc_5abc.py` (circuit diagrams + analytic spec);
  `draw_pqc_ansatz.py`. Write-up: `PQC_APPROX.md`, spec `PQC_CIRCUITS_FOR_ANALYSIS.md`.
- **Dependencies:** `noisy_bell_state.py`, `pqec_cnot_threshold.py`, `pqec_resynth_noise.py`,
  `destructive_gadget.py` (comparison baselines), `requirements.txt`.
- **Data files:** `pqc_ring_L3_params.npy` (18-CNOT solution), `pqc_ring_pruned_params.npy`
  + `pqc_ring_pruned.json` (the learned 14-CNOT solution). *(Auxiliary relaxations also ship
  `pqc_ring_5b_params.npy`, `pqc_ring_5c_params.npy`, `pqc_ring_5c.json`.)*
- **Outputs:** `circuit_pqc_5a.png`, `circuit_pqc_5a_raw.png`, `pqc_ring_threshold.png`.
- **Key result:** learned 14-CNOT threshold is **~1.5–1.9×** the Step-4 14-CNOT circuit
  (`ε₂*` = 0.060 / 0.237 / 0.338 at `ε` = 0.1 / 0.4 / 0.6) — CNOT *arrangement*, not just
  count, sets the threshold.
- **Auxiliary (not part of the main flow):** `pqc_compile.py` (generic-ansatz squeeze),
  `pqc_ring_5b*.py` (isometry floor = 14), `pqc_ring_5c.py` + `pqc_noise_aware.py`
  (observable relaxation → 5 CNOTs; noise-aware training adds no gain), `test_inclass.py`,
  `test_5c_oos.py`. See `PQC_APPROX.md`.

## Reproduce

```bash
git switch step-5-learned-14cnot
pip install -r requirements.txt
python pqc_ring_ansatz.py      # (re)train the gadget-matched ansatz (saves L3 params)
python pqc_ring_prune.py       # prune 18 -> 14 (saves pruned params/json)
python pqc_ring_threshold.py   # learned 14-CNOT threshold + pqc_ring_threshold.png
python draw_pqc_5abc.py        # circuit diagrams + PQC_CIRCUITS_FOR_ANALYSIS.md
```
(The pruned solution is committed, so `pqc_ring_threshold.py` and `draw_pqc_5abc.py` run
without retraining.)

## git worktree

```bash
git worktree add ../pqec-step5 step-5-learned-14cnot
```

> Every step branch is a full copy of `main`; this guide marks the *core* files for Step 5.
> The integrated/stable branch is `main`.
