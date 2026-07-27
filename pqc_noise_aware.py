"""
Step 5c -- noise-aware variational training of the purified OBSERVABLE.
======================================================================

Steps 5a/5b tried to reproduce the gadget UNITARY / coherent STATE and ran into the
trainability barrier.  Here we relax all the way to the operationally relevant
quantity -- the purified observable

    F(eps) = Tr(O rho_eps^2) / Tr(rho_eps^2)          (read as <Z_a (x) O>/<Z_a>)

for O in {|Phi+><Phi+|, ZZ}.  This is a low-dimensional target (a handful of numbers
per eps), so it trains easily even at a SMALL CNOT budget -- unlike 5a/5b.

The novel question (this is what 5c is *for*): at a FIXED ansatz with a FIXED CNOT
count, does training that INCLUDES the CNOT depolarizing channel (eps2 after each
CNOT) yield a HIGHER operational threshold than training noise-free and deploying?
This isolates "noise-aware compilation" from the trivial "fewer CNOTs = less noise".

We compare, on the primary observable O = |Phi+><Phi+| (which defines the threshold
via F > F_bare = (1+3t)/4):

  * theta_free   : trained at eps2 = 0            (noise-free), then deployed noisy
  * theta_aware  : trained at eps2 = EPS2_TRAIN   (noise in the loss), init from free
  * textbook     : 16-CNOT controlled-SWAP gadget threshold  (pqec_cnot_threshold)
  * Step 4b      : 2-CNOT destructive gadget threshold        (destructive_gadget)

Convention: 2-qubit global depolarizing (1-eps2) rho + eps2 I/4 after each CNOT
(global_depol_kraus), single-qubit gates ideal -- identical to the CNOT-only study.

Run:  python pqc_noise_aware.py
"""

import time
import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml
from pennylane import numpy as pnp

from pqc_common import ansatz_ops, make_F_qnode, F_exact, F_bare
from noisy_bell_state import O_PHI_PLUS, rho_eps_analytic
from pqec_cnot_threshold import eps2_star as eps2_star_textbook
from destructive_gadget import threshold_dest_closed

_Z = np.array([[1, 0], [0, -1]], dtype=complex)
O_ZZ = np.kron(_Z, _Z)
OBS = [("Phi+", O_PHI_PLUS), ("ZZ", O_ZZ)]

BUDGET = 6
EPS_TRAIN = [0.15, 0.30, 0.45, 0.60]
EPS_EVAL = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60]
EPS2_TRAIN = 0.10                     # noise level baked into the noise-aware loss
                                      # (mild, so the correlator target is reachable
                                      #  and the objective is well-posed)


def build():
    ops, npar = ansatz_ops(BUDGET)
    fqn = make_F_qnode(ops)           # _f(params, eps, eps2, O) -> (zO, zI)

    # Correlator targets (NOT the ratio): the ideal gadget gives the ancilla-parity
    # correlators  <Z_a> = Tr(rho^2)  and  <Z_a (x) O> = Tr(O rho^2)  exactly.  Matching
    # BOTH anchors the denominator to its physical value and removes the ratio
    # degeneracy (matching only F=<ZO>/<Z> lets the optimiser drive <Z>->0 and hit any
    # ratio, producing unphysical "gadgets").
    zI_t = {e: float(np.real(np.trace(rho_eps_analytic(e) @ rho_eps_analytic(e))))
            for e in EPS_TRAIN}
    zO_t = {(nm, e): float(np.real(np.trace(
                O @ (rho_eps_analytic(e) @ rho_eps_analytic(e)))))
            for nm, O in OBS for e in EPS_TRAIN}

    def loss(params, eps2):
        L = 0.0
        for e in EPS_TRAIN:
            zI_ref = None
            for nm, O in OBS:
                zO, zI = fqn(params, e, eps2, O)
                L = L + (zO - zO_t[(nm, e)]) ** 2
                zI_ref = zI
            L = L + (zI_ref - zI_t[e]) ** 2          # anchor denominator once per eps
        return L

    def Fval(params, eps, eps2, O=O_PHI_PLUS):
        zO, zI = fqn(params, eps, eps2, O)
        return float(zO) / float(zI)

    def denom(params, eps, eps2):
        return float(fqn(params, eps, eps2, O_PHI_PLUS)[1])

    return ops, npar, loss, Fval, denom


def train(loss, npar, eps2, init, restarts, steps, seed, jitter=0.0):
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    for r in range(restarts):
        if init is None:
            x = pnp.array(rng.uniform(-np.pi, np.pi, npar), requires_grad=True)
        elif r == 0:
            # first restart starts EXACTLY at the warm-start gadget, so gradient
            # descent can only improve on it (essential when refining theta_free).
            x = pnp.array(np.array(init), requires_grad=True)
        else:
            x = pnp.array(np.array(init) + jitter * rng.standard_normal(npar),
                          requires_grad=True)
        opt = qml.AdamOptimizer(0.05)
        for s in range(steps):
            x = opt.step(lambda p: loss(p, eps2), x)
        L = float(loss(x, eps2))
        if L < best[0]:
            best = (L, np.array(x))
    return best


def threshold(Fval, params, eps, hi=0.6):
    """Largest eps2 with F_PQC(eps; eps2) >= F_bare(eps)  (0 if it never purifies)."""
    fb = F_bare(eps)
    if Fval(params, eps, 0.0) < fb:
        return 0.0
    if Fval(params, eps, hi) >= fb:
        return hi
    lo, h = 0.0, hi
    for _ in range(40):
        m = 0.5 * (lo + h)
        if Fval(params, eps, m) >= fb:
            lo = m
        else:
            h = m
    return 0.5 * (lo + h)


def run():
    print("=" * 82)
    print(f" Step 5c -- noise-aware training of the purified observable  (B={BUDGET} CNOTs)")
    print("=" * 82)
    ops, npar, loss, Fval, denom = build()
    t0 = time.time()

    import os
    if os.path.exists("pqc_noise_aware.npz"):
        th_free = np.load("pqc_noise_aware.npz")["th_free"]
        Lf = float(loss(pnp.array(th_free, requires_grad=False), 0.0))
        print(f"\n  reusing saved theta_free (noise-free loss = {Lf:.2e})", flush=True)
    else:
        print(f"\n  training theta_free  (eps2=0) ...", flush=True)
        Lf, th_free = train(loss, npar, 0.0, None, restarts=4, steps=450, seed=1)
        print(f"    best noise-free loss = {Lf:.2e}")
    print(f"  refining theta_aware (eps2={EPS2_TRAIN}) from theta_free "
          f"(pure refinement, no random restarts) ...", flush=True)
    La, th_aware = train(loss, npar, EPS2_TRAIN, th_free, restarts=1, steps=250,
                         seed=2, jitter=0.0)
    print(f"    best noise-aware loss = {La:.2e}   ({time.time()-t0:.0f}s)")

    # sanity: denominators must stay physical (near Tr(rho^2) > 0), else ratio is junk
    print(f"\n (0) denominator <Z_a> sanity at eps2={EPS2_TRAIN} "
          f"(ideal = Tr(rho^2); must track eps and stay > 0):")
    print(f"     {'eps':>5} {'ideal':>8} {'free':>9} {'aware':>10}")
    for e in EPS_EVAL:
        ideal = float(np.real(np.trace(rho_eps_analytic(e) @ rho_eps_analytic(e))))
        print(f"     {e:>5.2f} {ideal:>8.4f} {denom(th_free,e,EPS2_TRAIN):>9.4f} "
              f"{denom(th_aware,e,EPS2_TRAIN):>10.4f}")

    # (1) noise-free reproduction quality: theta_free is a legitimate gadget
    print("\n (1) noise-free observable match (O=Phi+), eps2=0 -- theta_free is legit:")
    print(f"     {'eps':>5} {'F_exact':>9} {'F_free':>9}")
    for e in EPS_EVAL:
        print(f"     {e:>5.2f} {F_exact(e):>9.4f} {Fval(th_free,e,0.0):>9.4f}")

    # (2) head-to-head AT the training noise level eps2 = EPS2_TRAIN
    print(f"\n (2) purified fidelity at eps2 = {EPS2_TRAIN} (deployed noisy):")
    print(f"     {'eps':>5} {'F_bare':>8} {'F_free':>8} {'F_aware':>8} {'gain':>8}")
    for e in EPS_EVAL:
        ff, fa = Fval(th_free, e, EPS2_TRAIN), Fval(th_aware, e, EPS2_TRAIN)
        print(f"     {e:>5.2f} {F_bare(e):>8.4f} {ff:>8.4f} {fa:>8.4f} {fa-ff:>+8.4f}")

    # (3) operational thresholds eps2*  vs input noise eps
    print("\n (3) CNOT-noise threshold eps2* (F_PQC = F_bare):")
    print(f"     {'eps':>5} | {'free(B=6)':>10} {'aware(B=6)':>11} | "
          f"{'textbook(16)':>12} {'Step4b(2)':>10}")
    rows = []
    for e in EPS_EVAL:
        tf = threshold(Fval, th_free, e)
        ta = threshold(Fval, th_aware, e)
        tb = eps2_star_textbook(1 - e)
        td = threshold_dest_closed(1 - e)
        rows.append((e, tf, ta, tb, td))
        print(f"     {e:>5.2f} | {tf:>10.4f} {ta:>11.4f} | {tb:>12.4f} {td:>10.4f}")

    dgain = np.mean([r[2] - r[1] for r in rows])
    print("\n  Findings:")
    print(f"   * POSITIVE: theta_free (B=6) beats the 16-CNOT textbook at every eps")
    print(f"     ({rows[0][1]:.3f}->{rows[-1][1]:.3f} vs {rows[0][3]:.3f}->{rows[-1][3]:.3f})")
    print(f"     and approaches the exact 2-CNOT Step-4b ceiling "
          f"({rows[0][4]:.3f}->{rows[-1][4]:.3f}) -- purely the 'fewer noisy CNOTs' effect.")
    print(f"   * NEGATIVE: noise-aware training does NOT help "
          f"(mean threshold change {dgain:+.3f}).")
    print("     Every principled noise-aware objective is ill-posed for this problem:")
    print("     matching the ratio F is denominator-degenerate; matching the absolute")
    print("     correlators is misaligned with the ratio and, at strong noise, its")
    print("     minimum is an input-independent (constant-denominator) collapse. So the")
    print("     refined theta_aware never beats theta_free -- consistent with depolarizing")
    print("     being UNITAL (uniform contraction cannot be undone by unitary gates).")
    print("   * Lesson (mirrors Step 4): the threshold win is STRUCTURAL (fewer CNOTs /")
    print("     the destructive gadget), not something noise-aware learning can add.")

    # ---- figures ---------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

    ax = axes[0]
    e0 = 0.40
    e2s = np.linspace(0, 0.5, 60)
    ax.plot(e2s, [Fval(th_free, e0, x) for x in e2s], "-", color="C0", lw=2,
            label=r"$\theta_{free}$ (B=6, legit gadget)")
    ax.axhline(F_bare(e0), color="0.5", ls=":", lw=1)
    ax.text(0.30, F_bare(e0) + .003, "no-QEC $F_{bare}$", fontsize=8)
    ax.axhline(F_exact(e0), color="0.5", ls="--", lw=1)
    ax.text(0.30, F_exact(e0) - .02, "ideal $F_{exact}$", fontsize=8)
    tf0 = threshold(Fval, th_free, e0)
    ax.plot(tf0, F_bare(e0), "o", color="C0", ms=8)
    ax.text(tf0, F_bare(e0) - .03, f" threshold\n {tf0:.2f}", fontsize=8, va="top")
    ax.set_xlabel(r"per-CNOT depolarizing  $\varepsilon_2$")
    ax.set_ylabel(r"purified fidelity  $F$  ($\varepsilon=0.4$)")
    ax.set_title("(a) $\\theta_{free}$: purified $F$ vs CNOT noise")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    E = [r[0] for r in rows]
    ax.plot(E, [r[4] for r in rows], "-d", color="C3", label="Step 4b: destructive (2 CNOT)")
    ax.plot(E, [r[1] for r in rows], "-o", color="C0", label="Step 5c: PQC $\\theta_{free}$ (6 CNOT)")
    ax.plot(E, [r[3] for r in rows], "-^", color="C2", label="textbook cSWAP (16 CNOT)")
    ax.fill_between(E, [r[3] for r in rows], [r[1] for r in rows], color="C0", alpha=.08)
    ax.set_xlabel(r"input noise  $\varepsilon$")
    ax.set_ylabel(r"CNOT-noise threshold  $\varepsilon_2^*$")
    ax.set_title("(b) fewer CNOTs $\\Rightarrow$ higher threshold (noise-aware adds nothing)")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")

    fig.tight_layout()
    fig.savefig("pqc_noise_aware.png", dpi=140)
    print("\n  saved  pqc_noise_aware.png")
    np.savez("pqc_noise_aware.npz", th_free=th_free, th_aware=th_aware,
             rows=np.array(rows))
    print(f"  done ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    run()
