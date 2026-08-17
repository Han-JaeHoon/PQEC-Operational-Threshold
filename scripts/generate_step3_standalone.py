"""
Standalone Step 3 circuit figure for the calculation note (vector PDF).
======================================================================

Renders ONLY the textbook 16-CNOT decomposition used in Step 3 as a clean, caption-free
standalone figure (small title on top, one small D_q note at the bottom). The gate
sequence is the SAME object Step 3 simulates: it is captured from the actual
verify_analytic_decomposed._fred replay (via the helpers in generate_paper_figures),
so the circuit here is identical to figures/step3_textbook_16cnot_circuit.png -- only the
output format (vector PDF) and file name differ.

Per Controlled-SWAP(a; x, y) the textbook decomposition is (y = swap/Toffoli target):
  C1=CNOT(y->x); H_y; C2=CNOT(x->y); Tdg_y; C3=CNOT(a->y); T_y; C4=CNOT(x->y); Tdg_y;
  C5=CNOT(a->y); T_y,T_x; C6=CNOT(a->x); H_y; T_a,Tdg_x; C7=CNOT(a->x); C8=CNOT(y->x)
i.e. 8 CNOTs per Controlled-SWAP, single-qubit gates ideal. The first acts on (a,A1,B1),
the second on (a,A2,B2); the ancilla carries an H before and after.

Output:  figures/step3_16cnot_circuit.pdf   (vector; primary deliverable)
         figures/step3_16cnot_circuit.png   (raster preview, same drawing)

Run:  python scripts/generate_step3_standalone.py
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

# Reuse the EXACT Step-3 gate sequence + geometry from the main figure script,
# without modifying it.
from generate_paper_figures import (
    _tape_ops, _step3_qfunc, _ops_to_items, _assign_columns, _y, WIRES5,
    DPI, COL_DX, BOX_W, BOX_H, FS_GATE, FS_WIRE, FS_NOTE, LW,
)

FIG = os.path.join(ROOT, "figures")
OUT_PDF = os.path.join(FIG, "step3_16cnot_circuit.pdf")
OUT_PNG = os.path.join(FIG, "step3_16cnot_circuit.png")

TITLE = "Step 3: textbook two-Controlled-SWAP decomposition (16 CNOTs)"
NOTE = (r"A two-qubit replacement depolarizing channel $D_q$ is applied after "
        r"every CNOT.")


def build():
    ops = _tape_ops(_step3_qfunc)
    items, ncx, cx_seq = _ops_to_items(ops)
    assert ncx == 16, f"Step 3 must have 16 CNOTs, got {ncx}"
    return items, ncx, cx_seq


def render(items, col_scale=0.92):
    nw = len(WIRES5)
    ncols = _assign_columns(items, nw)
    dx = COL_DX * col_scale
    x0 = 0.0

    def cx(col):
        return x0 + (col + 0.5) * dx

    xL, xR = -0.35 * dx, ncols * dx + 0.15 * dx

    fig_w = max(6.0, 1.1 * (xR - xL) + 1.6)
    fig_h = 0.95 * nw + 1.6
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off()
    ax.patch.set_visible(False)

    # wires + labels
    for r in range(nw):
        y = _y(r, nw)
        ax.plot([xL, xR], [y, y], color="k", lw=1.0, zorder=1)
        ax.text(xL - 0.18 * dx, y, WIRES5[r], ha="right", va="center", fontsize=FS_WIRE)

    def box(col, row, label):
        x, y = cx(col), _y(row, nw)
        ax.add_patch(Rectangle((x - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H,
                     facecolor="white", edgecolor="k", lw=LW, zorder=3))
        ax.text(x, y, label, ha="center", va="center", fontsize=FS_GATE, zorder=4)

    for it in items:
        col = it["col"]
        x = cx(col)
        if it["type"] == "box":
            box(col, it["wire"], it["label"])
        else:  # cnot
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

    # center the small title/note over the drawn circuit content
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    xs = []
    for art in list(ax.lines) + list(ax.patches) + list(ax.texts):
        try:
            bb = art.get_window_extent(renderer=r)
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
    return fig


def main():
    items, ncx, seq = build()
    fig = render(items)
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")   # vector
    fig.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight", facecolor="white")  # preview
    plt.close(fig)

    from PIL import Image
    w, h = Image.open(OUT_PNG).size
    print(f"wrote {os.path.relpath(OUT_PDF, ROOT)}  (vector PDF)")
    print(f"wrote {os.path.relpath(OUT_PNG, ROOT)}  ({w}x{h}px preview, {DPI} dpi)")
    print(f"  CNOT count = {ncx} (expect 16)")
    print(f"  CNOT sequence = {seq}")


if __name__ == "__main__":
    main()
