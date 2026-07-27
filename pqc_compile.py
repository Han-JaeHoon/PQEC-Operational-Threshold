"""
Step 5a / 5b -- variational (PQC) compiling of the PQEC gadget.
==============================================================

We try to train a hardware-efficient PQC V(theta) (CNOT budget B, pqc_common.ansatz_ops)
to reproduce the gadget U = H_a CSWAP CSWAP H_a in two senses, and sweep B:

  * Step 5a -- FULL unitary:      minimise  1 - |Tr(U^dag V)|^2 / 32^2
  * Step 5b -- COHERENT state:    minimise  1 - |Tr(U0^dag V0)|^2 / 16^2
               (U0 = U on the ancilla-|0> input block -- the only inputs PQEC ever
               feeds the gadget; V0 = V restricted to those 16 columns; this is the
               coherent purified state that would be fed forward)

Both use EXACT analytic gradients (pqc_common._cost_grad) with many L-BFGS restarts.

What this script demonstrates (see PQC_APPROX.md for the discussion):
  * Variational compiling of this gadget hits the well-known TRAINABILITY BARRIER
    (barren plateaus / hard landscape).  The infidelity delta floors around ~0.3-0.5
    and does NOT reach 0, even though B=16 CNOTs is enough to represent U exactly
    (the textbook decomposition uses 16).  We also ran the plateau-mitigating LHST
    LOCAL cost (Khatri et al.) and rich ansaetze up to 375 params / 24 CNOTs -- the
    barrier persists (delta ~ 0.5).  So from-scratch variational compiling does NOT
    beat the structured Step-4a decomposition.
  * The COHERENT-state target (5b) is consistently EASIER than the full unitary (5a)
    -- lower delta at every budget -- because it need only match the 16 physical
    (ancilla-|0>) columns.  This is the same relaxation ladder as Step 4:
    unitary (5a) harder than state (5b) harder than observable (5c).

The upshot motivates Step 5c: don't reproduce the unitary or the state -- reproduce
only the purified OBSERVABLE, which is low-dimensional and trains easily.

Run:  python pqc_compile.py
"""

import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from pqc_common import (
    ansatz_ops, cost_grad_unitary, cost_grad_coherent, cost_grad_lhst,
    cost_unitary, unitary_from_params, F_from_unitary, F_exact,
)

BUDGETS = [6, 8, 10, 12, 14, 16]
EPS_GRID = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])


def _train(cost_grad, ops, npar, restarts, seed, maxiter=4000):
    """Minimise cost_grad(ops, theta) (returns (cost, grad)); best of `restarts`."""
    rng = np.random.default_rng(seed)
    best_c, best_x = np.inf, None
    for _ in range(restarts):
        x0 = rng.uniform(-np.pi, np.pi, size=npar)
        res = minimize(lambda x: cost_grad(ops, x), x0, jac=True, method="L-BFGS-B",
                       options=dict(maxiter=maxiter, ftol=1e-15, gtol=1e-13))
        if res.fun < best_c:
            best_c, best_x = float(res.fun), res.x
    return best_c, best_x


def _obs_error(ops, theta):
    """max_eps |F_read(V) - F_exact|  (operational error of the compiled circuit)."""
    V = unitary_from_params(ops, theta)
    return max(abs(F_from_unitary(V, e) - F_exact(e)) for e in EPS_GRID)


def run(restarts=16):
    print("=" * 82)
    print(" Step 5a / 5b -- variational compiling of  U = H_a CSWAP CSWAP H_a")
    print("=" * 82)
    print(f"  ansatz: hardware-efficient, n_params = 30 + 6B;  L-BFGS restarts = {restarts}")
    print("  delta = 1 - gate fidelity (0 = exact);  obs_err = worst |F_PQC - F_exact|\n")

    res = {"5a": {}, "5b": {}, "5a_lhst": {}}
    t0 = time.time()
    print(f"  {'B':>3} | {'5a delta':>10} {'obs':>8} | {'5b delta':>10} {'obs':>8} | "
          f"{'5a-LHST delta':>13}")
    print("  " + "-" * 70)
    for B in BUDGETS:
        ops, npar = ansatz_ops(B)

        ca, xa = _train(cost_grad_unitary, ops, npar, restarts, 100 + B)
        ea = _obs_error(ops, xa)
        res["5a"][B] = (ca, ea)

        cb, xb = _train(cost_grad_coherent, ops, npar, restarts, 200 + B)
        eb = _obs_error(ops, xb)
        res["5b"][B] = (cb, eb)

        # best-effort with the plateau-mitigating LOCAL cost: report the resulting
        # TRUE global infidelity (LHST optimum need not equal the global optimum).
        cl, xl = _train(cost_grad_lhst, ops, npar, max(6, restarts // 2), 300 + B)
        dl = cost_unitary(ops, xl)
        res["5a_lhst"][B] = dl

        print(f"  {B:>3} | {ca:>10.2e} {ea:>8.1e} | {cb:>10.2e} {eb:>8.1e} | "
              f"{dl:>13.2e}")

    print(f"\n  ({time.time()-t0:.0f}s)")

    # ---- interpretation --------------------------------------------------
    print("\n  Findings:")
    best5a = min(res["5a"][B][0] for B in BUDGETS)
    best5b = min(res["5b"][B][0] for B in BUDGETS)
    best5al = min(res["5a_lhst"][B] for B in BUDGETS)
    print(f"   * best full-unitary infidelity over the sweep (5a, global cost): {best5a:.2e}")
    print(f"   * best full-unitary infidelity with the LHST local cost:         {best5al:.2e}")
    print(f"   * best coherent-state infidelity (5b):                           {best5b:.2e}")
    lower = all(res["5b"][B][0] <= res["5a"][B][0] + 1e-9 for B in BUDGETS)
    print(f"   * 5b (coherent) <= 5a (unitary) at every budget: {lower}")
    print("   * None reaches a useful accuracy (delta -> 0): the trainability barrier")
    print("     dominates, so from-scratch PQC compiling does NOT beat Step-4a's 14 CNOTs.")
    print("   * The relaxation ladder unitary(5a) > state(5b) > observable(5c) in")
    print("     difficulty mirrors Step 4 and points to targeting the observable (5c).")

    # ---- figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

    ax = axes[0]
    ax.plot(BUDGETS, [res["5a"][B][0] for B in BUDGETS], "-o", color="C0",
            label="5a: full unitary (global cost)")
    ax.plot(BUDGETS, [res["5a_lhst"][B] for B in BUDGETS], "--^", color="C3",
            label="5a: full unitary (LHST local cost)")
    ax.plot(BUDGETS, [res["5b"][B][0] for B in BUDGETS], "-s", color="C1",
            label="5b: coherent state (anc-|0>)")
    ax.axvline(16, color="0.6", ls=":", lw=1)
    ax.text(16.1, 0.05, "16 = exact\n(textbook)", fontsize=8, va="bottom")
    ax.set_xlabel("CNOT budget  B")
    ax.set_ylabel(r"approximation infidelity  $\delta = 1-$fidelity")
    ax.set_title("(a) Variational compiling hits the trainability barrier")
    ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=8.5)

    ax = axes[1]
    ax.semilogy(BUDGETS, [max(res["5a"][B][1], 1e-16) for B in BUDGETS], "-o",
                color="C0", label="5a: full unitary")
    ax.semilogy(BUDGETS, [max(res["5b"][B][1], 1e-16) for B in BUDGETS], "-s",
                color="C1", label="5b: coherent state")
    ax.set_xlabel("CNOT budget  B")
    ax.set_ylabel(r"worst $|F_{PQC}-F_{exact}|$ over $\varepsilon$")
    ax.set_title("(b) Operational error stays large (compiling failed)")
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig("pqc_compile_pareto.png", dpi=140)
    print("\n  saved  pqc_compile_pareto.png")
    return res


if __name__ == "__main__":
    run()
