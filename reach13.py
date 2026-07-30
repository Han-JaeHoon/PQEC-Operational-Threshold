"""Rigorous test: can ANY single CNOT be removed from the 14-CNOT solution to give an
exact 13-CNOT compilation? For each of the 14 CNOTs, remove it and train with many
restarts (warm-start + 11 random). Report the best 13-CNOT delta found."""
import json, numpy as np
from scipy.optimize import minimize
from pqc_ring_prune import ansatz_masked, SEQ, NPAR
from pqc_ring_ansatz import cost_grad, _I32, _obs_err
from pqc_common import U_TARGET

info = json.load(open("pqc_ring_pruned.json"))
mask14 = info["mask"]
theta14 = np.load("pqc_ring_pruned_params.npy")
active = [i for i in range(len(SEQ)) if mask14[i]]
print(f"14-CNOT solution active slots: {active}\n", flush=True)

rng = np.random.default_rng(1)
best_overall = (np.inf, None)
for i in active:
    m = mask14.copy(); m[i] = False
    ops = ansatz_masked(m)
    starts = [theta14] + [rng.uniform(-np.pi, np.pi, NPAR) for _ in range(11)]
    best = np.inf
    for x0 in starts:
        r = minimize(lambda x: cost_grad(ops, x, U_TARGET, _I32), x0, jac=True,
                     method="L-BFGS-B", options=dict(maxiter=8000, ftol=1e-15, gtol=1e-13))
        best = min(best, float(r.fun))
    print(f"  remove slot {i:>2} {str(SEQ[i]):>7} -> best 13-CNOT delta = {best:.2e}", flush=True)
    if best < best_overall[0]:
        best_overall = (best, i)
print(f"\n  best achievable 13-CNOT delta over all removals = {best_overall[0]:.2e} "
      f"(slot {best_overall[1]})")
print("  => 13 CNOTs " + ("ACHIEVABLE" if best_overall[0] < 1e-6 else "NOT reached (14 is the floor for this ansatz)"))
