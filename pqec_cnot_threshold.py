"""
CNOT-only noise threshold for decomposed-Fredkin PQEC (single-qubit gates ideal).
=================================================================================

Model: the textbook decomposition (8 CNOTs per Fredkin, two Fredkins = 16 CNOTs);
single-qubit gates are IDEAL and a two-qubit depolarizing channel of strength
eps2 acts after EACH CNOT.  With u = 1-eps1 = 1 (ideal single-qubit gates), the
two orientation numerators coincide,

    C(1, v) = v^5 + 5 v^6,     v = 1 - eps2,

so the purified Bell fidelity is orientation-independent:

    F_dec(eps2, t) = 1/4 [ 1 + t(1+t)(v^5 + 5 v^6) / (1 + 3 v^4 t^2) ],
    F_bare(t)      = (1 + 3 t) / 4,     t = 1 - eps  (global Bell-depol input).

PQEC helps iff F_dec > F_bare, i.e.

    (1 + t)(v^5 + 5 v^6) > 3 (1 + 3 v^4 t^2).

The CNOT-noise threshold eps2*(t) is the root of the equality.  Above eps2* the
noisy CNOTs inject more error than one PQEC round removes.

Run:  python pqec_cnot_threshold.py
"""

import numpy as np
import matplotlib.pyplot as plt

from verify_analytic_decomposed import circuit_AB   # genuine circuit (e1=0 => CNOT-only)

np.set_printoptions(precision=4, suppress=True)


# --- analytic (single-qubit ideal, u=1) ------------------------------------
def F_dec(eps2, t):
    v = 1 - eps2
    return 0.25 * (1 + t * (1 + t) * (v**5 + 5 * v**6) / (1 + 3 * v**4 * t**2))


def F_bare(t):
    return (1 + 3 * t) / 4


def eps2_star(t):
    """Root of F_dec = F_bare in eps2 (0 if PQEC never helps)."""
    f = lambda e: F_dec(e, t) - F_bare(t)
    if f(1e-9) <= 0:
        return 0.0
    lo, hi = 0.0, 0.999
    for _ in range(80):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if f(m) > 0 else (lo, m)
    return 0.5 * (lo + hi)


def _eps_to_local_p(eps):
    """Global-input eps -> local depolarizing p with t=(1-4p/3)^2 = 1-eps (for the
    circuit, which prepares the isotropic state via local depolarizing)."""
    return 0.75 * (1 - np.sqrt(1 - eps))


# ===========================================================================
def main():
    print("=" * 74)
    print(" CNOT-only noise threshold (single-qubit gates ideal)")
    print("=" * 74)

    # (0) analytic vs genuine circuit, and orientation-independence
    print("\n (0) analytic F_dec vs circuit (e1=0), and retain==discard:")
    worst = worst_or = 0.0
    for eps in [0.2, 0.4, 0.6]:
        t = 1 - eps
        p = _eps_to_local_p(eps)
        for e2 in [0.05, 0.12, 0.20]:
            fr = (lambda z: z[0] / z[1])(circuit_AB(p, 0.0, e2, "retain"))
            fd = (lambda z: z[0] / z[1])(circuit_AB(p, 0.0, e2, "discard"))
            worst = max(worst, abs(fr - F_dec(e2, t)))
            worst_or = max(worst_or, abs(fr - fd))
    print(f"     max|circuit - analytic| = {worst:.2e}   "
          f"max|retain - discard| = {worst_or:.2e}")

    # (1) threshold table
    print("\n (1) CNOT-noise threshold eps2* vs input noise eps (t = 1-eps):")
    print(f"     {'eps':>5} {'t':>5} {'F_bare':>8} {'eps2*':>8} {'16*eps2* (budget)':>18}")
    for eps in [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60]:
        t = 1 - eps
        es = eps2_star(t)
        print(f"     {eps:>5.2f} {t:>5.2f} {F_bare(t):>8.4f} {es:>8.4f} {16*es:>18.3f}")
    print("     (16 = total CNOT count; 16*eps2* is a rough per-round CNOT budget.)")

    # (2) small-eps2 behaviour
    print("\n (2) near eps2=0:  d F_dec / d eps2 |_0  (= -K2, orientation-independent)")
    for eps in [0.2, 0.4, 0.6]:
        t = 1 - eps
        h = 1e-6
        slope = -(F_dec(h, t) - F_dec(0, t)) / h
        print(f"     eps={eps:.1f} (t={t:.1f}):  K2 = {slope:.4f}  "
              f"(=t(1+t)(33t^2+35)/(4(1+3t^2)^2))")

    # ---- Figure -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))

    ax = axes[0]
    e2s = np.linspace(0, 0.25, 120)
    for eps, c in [(0.20, "C0"), (0.40, "C1"), (0.60, "C2")]:
        t = 1 - eps
        ax.plot(e2s, [F_dec(e, t) for e in e2s], "-", color=c, lw=2,
                label=f"$\\varepsilon$={eps}")
        ax.axhline(F_bare(t), color=c, ls=":", lw=1)
        es = eps2_star(t)
        ax.plot(es, F_bare(t), "o", color=c, ms=7)
    ax.text(0.252, F_bare(0.20), " no-QEC", va="center", fontsize=8)
    ax.set_xlabel(r"per-CNOT two-qubit depolarizing  $\varepsilon_2$")
    ax.set_ylabel(r"purified fidelity  $F_{dec}$")
    ax.set_title("(a) CNOT-only: $F_{dec}$ vs $\\varepsilon_2$ (dots = threshold)")
    ax.set_xlim(0, 0.27)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    es_in = np.linspace(0.02, 0.66, 60)
    ax.plot(es_in, [eps2_star(1 - e) for e in es_in], "-", color="C3", lw=2)
    ax.axhspan(0, 0.01, color="0.85", alpha=0.6)
    ax.text(0.35, 0.006, "realistic hardware CNOT error $\\sim10^{-2}$", fontsize=8,
            va="center")
    ax.set_xlabel(r"input noise  $\varepsilon$")
    ax.set_ylabel(r"CNOT-noise threshold  $\varepsilon_2^*$")
    ax.set_title("(b) $\\varepsilon_2^*$ grows with input noise")
    ax.set_ylim(0, None)

    fig.tight_layout()
    fig.savefig("pqec_cnot_threshold.png", dpi=140)
    print("\n  saved  pqec_cnot_threshold.png")


if __name__ == "__main__":
    main()
