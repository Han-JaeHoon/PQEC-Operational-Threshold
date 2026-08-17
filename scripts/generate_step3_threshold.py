"""
Standalone Step 3 threshold plot for the calculation note.
==========================================================

Plots the Step-3 one-round operational threshold q_th^(3) against the input noise
epsilon, using the supplied values verbatim (no recomputation). Caption-free; clean
matplotlib style on a white background for direct paper insertion.

Output:  figures/step3_threshold_vs_epsilon.png   (raster, as requested)
         figures/step3_threshold_vs_epsilon.pdf   (vector, same plot; suggested name)

Run:  python scripts/generate_step3_threshold.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
OUT_PNG = os.path.join(FIG, "step3_threshold_vs_epsilon.png")
OUT_PDF = os.path.join(FIG, "step3_threshold_vs_epsilon.pdf")

# supplied data (used verbatim)
EPS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
Q_TH = [0.032991, 0.061194, 0.084499, 0.102949, 0.116657, 0.125715]

DPI = 300


def main():
    fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=DPI)
    ax.plot(EPS, Q_TH, "-o", color="C0", lw=1.9, ms=7,
            label="Step 3: textbook 16 CNOT")

    ax.set_xlabel(r"input noise  $\varepsilon$", fontsize=13)
    ax.set_ylabel(r"CNOT-noise threshold  $q_{\mathrm{th}}^{(3)}$", fontsize=13)
    ax.set_title("Operational threshold vs input noise", fontsize=14)
    ax.grid(True, ls=":", lw=0.6, color="0.7")
    ax.legend(frameon=False, fontsize=11, loc="upper left")
    ax.set_xlim(0.08, 0.62)
    ax.set_ylim(0.0, 0.135)

    fig.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")   # vector
    plt.close(fig)

    from PIL import Image
    w, h = Image.open(OUT_PNG).size
    print(f"wrote {os.path.relpath(OUT_PNG, ROOT)}  ({w}x{h}px, {DPI} dpi)")
    print(f"wrote {os.path.relpath(OUT_PDF, ROOT)}  (vector PDF)")
    print(f"  data points: {list(zip(EPS, Q_TH))}")


if __name__ == "__main__":
    main()
