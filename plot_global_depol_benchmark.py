"""
Figure: 3-qubit global depolarizing on the Fredkins -- signal loss, no threshold.
================================================================================

(a) The purified fidelity F_PQEC(eps, g_F) is flat in g_F for every eps (the
    (1-g_F)^2 visibility cancels in the ratio) and stays above the bare fidelity
    F_bare = 1-3eps/4 for all 0 < eps < 3/4.
(b) The only cost is sampling: the parity signal B ∝ (1-g_F)^2, so the shot
    overhead diverges as N_samp ~ (1-g_F)^{-4}.

Saves global_depol_benchmark.png.  Run:  python plot_global_depol_benchmark.py
"""

import numpy as np
import matplotlib.pyplot as plt

from pqec_gadget_noise import obs_pqec_noisy, no_qec


def main():
    gs = np.linspace(0.0, 0.95, 40)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    ax = axes[0]
    for eps, c in [(0.20, "C0"), (0.40, "C1"), (0.60, "C2")]:
        F = [obs_pqec_noisy(eps, g) for g in gs]
        ax.plot(gs, F, "-", color=c, lw=2, label=f"$F_{{PQEC}}$, $\\varepsilon$={eps}")
        ax.axhline(no_qec(eps), color=c, ls="--", lw=1)
    ax.text(0.02, no_qec(0.20) + 0.01, "dashed = $F_{bare}=1-3\\varepsilon/4$",
            fontsize=8)
    ax.set_xlabel(r"Fredkin global-depol strength  $g_F$")
    ax.set_ylabel(r"purified fidelity  $F_{PQEC}$")
    ax.set_title("(a) $F_{PQEC}$ is independent of $g_F$ (and $> F_{bare}$)")
    ax.set_ylim(0.3, 1.02)
    ax.legend(frameon=False, fontsize=9, loc="center right")

    ax = axes[1]
    ax.plot(gs, 1 / (1 - gs) ** 4, "k-", lw=2)
    for g in [0.3, 0.5, 0.9]:
        ax.plot(g, 1 / (1 - g) ** 4, "o", color="C3")
        ax.annotate(f"{1/(1-g)**4:.0f}x", (g, 1 / (1 - g) ** 4),
                    textcoords="offset points", xytext=(-28, 2), fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel(r"Fredkin global-depol strength  $g_F$")
    ax.set_ylabel(r"sampling overhead  $N_{samp}(g_F)/N_{samp}(0)$")
    ax.set_title(r"(b) Cost is sampling only:  $\sim (1-g_F)^{-4}$")

    fig.tight_layout()
    fig.savefig("global_depol_benchmark.png", dpi=140)
    print("saved  global_depol_benchmark.png")


if __name__ == "__main__":
    main()
