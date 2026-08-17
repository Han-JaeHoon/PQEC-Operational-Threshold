"""
Two-row rendering of the Step 5 learned 14-CNOT circuit (for the paper body).
=============================================================================

Renders the SAME circuit as figures/step5_learned_14cnot_circuit.png, but broken across
two rows so it fits a portrait-page \\textwidth without being unreadably wide. The CNOT
topology and L_k numbering are taken from the actual saved pruning data
(pqc_ring_pruned.json + pqc_ring_pruned_params.npy via draw_pqc_5abc.merged_circuit),
and asserted against the required sequence before drawing.

Row split (no gate duplicated across rows):
  Row 1:  L0 C1 L1 C2 L2 C3 L3 C4 L4 C5 L5 C6 L6 C7 L7
  Row 2:  (from L7) C8 L8 C9 L9 C10 L10 C11 L11 C12 L12 C13 L13 C14 L14

Output (new file; does NOT touch the existing one-row figure):
  figures/step5_learned_14cnot_circuit_two_rows.png

Run:  python scripts/generate_step5_two_rows.py
"""
import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import draw_pqc_5abc as d5

FIG = os.path.join(ROOT, "figures")
OUT = os.path.join(FIG, "step5_learned_14cnot_circuit_two_rows.png")

# ---- style (matched to the other paper figures) ----------------------------
DPI = 300
DX = 1.0            # column spacing
DY = 0.92           # wire spacing within a row
ROW_GAP = 1.9       # vertical gap between the two 5-wire rows
BOX_W, BOX_H = 0.60, 0.64
FS_WIRE = 13
FS_L = 13
FS_C = 10
FS_NOTE = 11
FS_TITLE = 15.5
LW = 1.5

WIRES = [r"$a$", r"$A_1$", r"$A_2$", r"$B_1$", r"$B_2$"]   # rows 0..4
NW = 5

REQUIRED = [(0, 1), (1, 3), (0, 4), (2, 4), (0, 3), (1, 3), (0, 4),
            (2, 4), (0, 1), (0, 3), (1, 3), (0, 2), (0, 4), (2, 4)]


def _load_cnots():
    mask = json.load(open(os.path.join(ROOT, "pqc_ring_pruned.json")))["mask"]
    params = np.load(os.path.join(ROOT, "pqc_ring_pruned_params.npy"))
    _, cnots = d5.merged_circuit(mask, params)
    return [(int(c), int(t)) for c, t in cnots]


def _row_items(kind):
    """Build the ordered item list for a row. Uses the real cnots (asserted == REQUIRED)."""
    cn = _load_cnots()
    assert cn == REQUIRED, "saved CNOT sequence does not match the required C1..C14"
    if kind == "top":
        # L0 C1 L1 ... C7 L7
        items = [{"type": "L", "k": 0}]
        for i in range(0, 7):                 # C1..C7 -> cn[0..6], L1..L7
            c, t = cn[i]
            items.append({"type": "cnot", "ctrl": c, "targ": t, "ck": i + 1})
            items.append({"type": "L", "k": i + 1})
        return items
    else:
        # (from L7) C8 L8 ... C14 L14
        items = []
        for i in range(7, 14):                # C8..C14 -> cn[7..13], L8..L14
            c, t = cn[i]
            items.append({"type": "cnot", "ctrl": c, "targ": t, "ck": i + 1})
            items.append({"type": "L", "k": i + 1})
        return items


def _draw_row(ax, items, y_top, x0):
    """Draw one 5-wire row; return (xL, xR, y_mid)."""
    ncols = len(items)

    def yw(i):
        return y_top - i * DY

    def cx(col):
        return x0 + (col + 0.5) * DX

    xL, xR = x0 - 0.30 * DX, x0 + ncols * DX + 0.10 * DX
    for i in range(NW):
        ax.plot([xL, xR], [yw(i), yw(i)], color="k", lw=1.0, zorder=1)
        ax.text(xL - 0.16 * DX, yw(i), WIRES[i], ha="right", va="center", fontsize=FS_WIRE)

    for col, it in enumerate(items):
        x = cx(col)
        if it["type"] == "L":
            ytop, ybot = yw(0), yw(NW - 1)
            ax.add_patch(Rectangle((x - BOX_W / 2, ybot - BOX_H / 2), BOX_W,
                         (ytop - ybot) + BOX_H, facecolor="white", edgecolor="k",
                         lw=LW, zorder=3))
            ax.text(x, (ytop + ybot) / 2, rf"$L_{{{it['k']}}}$", ha="center",
                    va="center", fontsize=FS_L, zorder=4)
        else:
            yc, yt = yw(it["ctrl"]), yw(it["targ"])
            ax.plot([x, x], [yc, yt], color="k", lw=LW, zorder=2)
            ax.add_patch(Circle((x, yc), 0.09, color="k", zorder=4))
            ax.add_patch(Circle((x, yt), 0.19, facecolor="white", edgecolor="k",
                         lw=LW, zorder=4))
            ax.plot([x - 0.19, x + 0.19], [yt, yt], color="k", lw=LW, zorder=5)
            ax.plot([x, x], [yt - 0.19, yt + 0.19], color="k", lw=LW, zorder=5)
            ax.text(x, yw(NW - 1) - 0.62, rf"$C_{{{it['ck']}}}$", ha="center",
                    va="top", fontsize=FS_C, color="0.25")
    return xL, xR, (yw(0) + yw(NW - 1)) / 2


def main():
    top = _row_items("top")
    bot = _row_items("bot")

    # verification
    cn = _load_cnots()
    n_cnot = sum(1 for it in top + bot if it["type"] == "cnot")
    seq = [(it["ctrl"], it["targ"]) for it in (top + bot) if it["type"] == "cnot"]
    Ls = [it["k"] for it in (top + bot) if it["type"] == "L"]
    assert n_cnot == 14, f"expected 14 CNOTs, got {n_cnot}"
    assert seq == REQUIRED == cn, "CNOT sequence mismatch"
    assert Ls == list(range(15)), f"L layers must be L0..L14 once each, got {Ls}"

    fig, ax = plt.subplots(figsize=(9.6, 6.2), dpi=DPI)
    ax.set_axis_off()

    y_top1 = 0.0
    row_h = (NW - 1) * DY
    y_top2 = y_top1 - row_h - ROW_GAP
    x0 = 0.0

    xL1, xR1, ymid1 = _draw_row(ax, top, y_top1, x0)
    xL2, xR2, ymid2 = _draw_row(ax, bot, y_top2, x0)

    # continuation markers (make it clear it is ONE circuit continuing).
    # Row 1 (right, at wire level -> clear space beyond L7):
    ax.annotate("", xy=(xR1 + 1.05 * DX, ymid1), xytext=(xR1 + 0.15 * DX, ymid1),
                arrowprops=dict(arrowstyle="-|>", color="0.15", lw=1.6))
    ax.text(xR1 + 1.15 * DX, ymid1, "continued\nbelow", ha="left", va="center",
            fontsize=FS_NOTE, color="0.15")
    # Row 2 (in the left-margin corner of the gap, clear of wire labels and C_k labels):
    ax.text(x0 - 0.45 * DX, y_top2 + 0.52, r"from $L_7$", ha="right",
            va="bottom", fontsize=FS_NOTE, color="0.15")
    ax.annotate("", xy=(x0 + 0.18 * DX, y_top2 + 0.06), xytext=(x0 - 0.35 * DX, y_top2 + 0.52),
                arrowprops=dict(arrowstyle="-|>", color="0.15", lw=1.6))

    # title + single-line note
    xmid = (min(xL1, xL2) + max(xR1, xR2)) / 2
    ax.text(xmid, y_top1 + 1.35, "Step 5: learned and pruned 14-CNOT circuit",
            ha="center", va="bottom", fontsize=FS_TITLE)
    ax.text(xmid, y_top2 - row_h - 1.35,
            r"$L_k=\bigotimes_j U_{(k,j)}$;  $D_q$ is applied after every CNOT.",
            ha="center", va="top", fontsize=FS_NOTE, color="0.1")

    ax.set_aspect("equal")
    # limits (leave room for arrows/labels); tight bbox crops the rest
    ax.set_xlim(min(xL1, xL2) - 2.3 * DX, max(xR1, xR2) + 2.7 * DX)
    ax.set_ylim(y_top2 - row_h - 2.1, y_top1 + 2.1)

    fig.savefig(OUT, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    from PIL import Image
    w, h = Image.open(OUT).size
    print(f"wrote {os.path.relpath(OUT, ROOT)}  ({w}x{h}px, {DPI} dpi, aspect {w/h:.2f}:1)")
    print(f"  CNOT count = {n_cnot} (expect 14)")
    print(f"  CNOT sequence == required == saved: {seq == REQUIRED == cn}")
    print(f"  sequence = {seq}")
    print(f"  L layers = {Ls}  (L0..L14 once each: {Ls == list(range(15))})")
    print(f"  row1 = L0..L7 + C1..C7 ; row2 = C8..C14 + L8..L14 (no duplication)")


if __name__ == "__main__":
    main()
