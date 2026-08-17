"""
Standalone Step 4 figures for the calculation note.
===================================================

(1) The resynthesized 14-CNOT PQEC circuit, rendered directly from the authoritative
    gate list pqec_resynth_noise.GATES (Qiskit transpile, basis {u, cx},
    optimization_level 3, seed_transpiler 7) -- so the drawing matches the simulated
    circuit exactly (14 CNOTs; Hilbert-Schmidt overlap with U_PQEC = 1, i.e. equal up to
    global phase).
(2) A Step 3 vs Step 4 operational-threshold comparison plot, using supplied values
    verbatim (no recomputation).

Caption-free (small in-figure title only). White background, clean paper style.

Outputs:
  figures/step4_resynthesized_14cnot_circuit.png   (standalone circuit)
  figures/step4_vs_step3_thresholds.png            (comparison graph)

Run:  python scripts/generate_step4_standalone.py
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from generate_paper_figures import (
    _assign_columns, _y, WIRES5, DPI, COL_DX, BOX_W, BOX_H,
    FS_GATE, FS_WIRE, FS_NOTE, LW,
)
import pqec_resynth_noise as r4

FIG = os.path.join(ROOT, "figures")
OUT_CIRCUIT = os.path.join(FIG, "step4_resynthesized_14cnot_circuit.png")
OUT_GRAPH = os.path.join(FIG, "step4_vs_step3_thresholds.png")

TITLE = "Step 4: resynthesized 14-CNOT circuit"
NOTE = (r"Qiskit transpile (basis $\{u,\mathrm{cx}\}$, optimization_level 3, "
        r"seed 7);  $U_4=e^{i\phi_4}U_{\mathrm{PQEC}}$.  A two-qubit replacement "
        r"depolarizing channel $D_q$ is applied after every CNOT.")


# ---------------------------------------------------------------------------
# (1) circuit  (rendered from the actual GATES list)
# ---------------------------------------------------------------------------
def _circuit_items():
    items, ncx, seq = [], 0, []
    for name, wires, params in r4.GATES:
        if name == "u":
            items.append({"type": "box", "wire": int(wires[0]), "label": r"$U_3$"})
        elif name == "cx":
            c, t = int(wires[0]), int(wires[1])
            items.append({"type": "cnot", "ctrl": c, "targ": t})
            ncx += 1
            seq.append((c, t))
    return items, ncx, seq


def render_circuit(items, col_scale=0.92):
    nw = len(WIRES5)
    ncols = _assign_columns(items, nw)
    dx = COL_DX * col_scale
    x0 = 0.0

    def cx(col):
        return x0 + (col + 0.5) * dx

    xL, xR = -0.35 * dx, ncols * dx + 0.15 * dx
    fig, ax = plt.subplots(figsize=(max(6.0, 1.1 * (xR - xL) + 1.6), 0.95 * nw + 1.6))
    ax.set_axis_off()
    ax.patch.set_visible(False)

    for r in range(nw):
        y = _y(r, nw)
        ax.plot([xL, xR], [y, y], color="k", lw=1.0, zorder=1)
        ax.text(xL - 0.18 * dx, y, WIRES5[r], ha="right", va="center", fontsize=FS_WIRE)

    for it in items:
        x = cx(it["col"])
        if it["type"] == "box":
            y = _y(it["wire"], nw)
            ax.add_patch(Rectangle((x - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H,
                         facecolor="white", edgecolor="k", lw=LW, zorder=3))
            ax.text(x, y, it["label"], ha="center", va="center", fontsize=FS_GATE, zorder=4)
        else:
            yc, yt = _y(it["ctrl"], nw), _y(it["targ"], nw)
            ax.plot([x, x], [yc, yt], color="k", lw=LW, zorder=2)
            ax.add_patch(Circle((x, yc), 0.085, color="k", zorder=4))
            ax.add_patch(Circle((x, yt), 0.185, facecolor="white", edgecolor="k",
                         lw=LW, zorder=4))
            ax.plot([x - 0.185, x + 0.185], [yt, yt], color="k", lw=LW, zorder=5)
            ax.plot([x, x], [yt - 0.185, yt + 0.185], color="k", lw=LW, zorder=5)

    ax.set_xlim(xL - 0.6 * dx, xR + 0.6 * dx)
    ax.set_ylim(_y(nw - 1, nw) - 1.25, _y(0, nw) + 1.25)
    ax.set_aspect("equal")

    # center title/note over the drawn content
    fig.canvas.draw()
    rr = fig.canvas.get_renderer()
    xs = []
    for art in list(ax.lines) + list(ax.patches) + list(ax.texts):
        try:
            bb = art.get_window_extent(renderer=rr)
            if bb.width > 0:
                xs += [bb.x0, bb.x1]
        except Exception:
            pass
    inv = ax.transData.inverted()
    xc = (inv.transform((min(xs), 0))[0] + inv.transform((max(xs), 0))[0]) / 2 if xs \
        else (xL + xR) / 2

    ax.text(xc, _y(0, nw) + 1.02, TITLE, ha="center", va="bottom", fontsize=FS_WIRE)
    ax.text(xc, _y(nw - 1, nw) - 0.78, NOTE, ha="center", va="top", fontsize=FS_NOTE,
            color="0.1")

    fig.savefig(OUT_CIRCUIT, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------------
# (2) threshold comparison  (supplied values, verbatim)
# ---------------------------------------------------------------------------
EPS = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60]
Q3 = [0.032991, 0.061194, 0.084499, 0.102949, 0.116657, 0.125715]
Q4 = [0.041263, 0.078817, 0.111911, 0.139861, 0.162075, 0.178027]


def render_graph():
    fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=DPI)
    ax.plot(EPS, Q3, "-^", color="C2", lw=1.9, ms=7, label="Step 3: textbook 16-CNOT")
    ax.plot(EPS, Q4, "-o", color="C0", lw=1.9, ms=7, label="Step 4: resynthesized 14-CNOT")
    ax.set_xlabel(r"input noise  $\varepsilon$", fontsize=13)
    ax.set_ylabel(r"CNOT-noise threshold  $q_{\mathrm{th}}$", fontsize=13)
    ax.set_title("Operational threshold comparison", fontsize=14)
    ax.grid(True, ls=":", lw=0.6, color="0.7")
    ax.legend(frameon=False, fontsize=11, loc="upper left")
    ax.set_xlim(0.08, 0.62)
    ax.set_ylim(0.0, 0.19)
    fig.savefig(OUT_GRAPH, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    items, ncx, seq = _circuit_items()
    assert ncx == 14, f"Step 4 must have 14 CNOTs, got {ncx}"
    render_circuit(items)
    render_graph()

    from PIL import Image
    wc, hc = Image.open(OUT_CIRCUIT).size
    wg, hg = Image.open(OUT_GRAPH).size
    print(f"wrote {os.path.relpath(OUT_CIRCUIT, ROOT)}  ({wc}x{hc}px, {DPI} dpi)")
    print(f"wrote {os.path.relpath(OUT_GRAPH, ROOT)}  ({wg}x{hg}px, {DPI} dpi)")
    print(f"  circuit CNOT count = {ncx} (expect 14); overlap with U = {r4.verify_unitary():.12f}")
    print(f"  circuit CNOT sequence = {seq}")
    print(f"  graph: Step 3 = {Q3}")
    print(f"  graph: Step 4 = {Q4}")


if __name__ == "__main__":
    main()
