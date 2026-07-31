"""
Generate the paper figures for the Operational-Threshold write-up.
==================================================================

This script ONLY renders figures. It does not change any circuit, formula, noise model,
or numerical result: the Step 3 / Step 4 / Step 5 circuit drawings are built *directly*
from the actual simulation objects, so the picture is identical to what is simulated:

  * Step 3 (textbook 16-CNOT)      -- replays verify_analytic_decomposed._fred (the exact
                                       8-CNOT/Fredkin decomposition) and renders the
                                       captured operation tape (16 CNOTs).
  * Step 4 (resynthesized 14-CNOT) -- renders pqec_resynth_noise.GATES verbatim (14 cx).
  * Step 5 (learned 14-CNOT)       -- renders draw_pqc_5abc.merged_circuit() from the
                                       saved pruning mask + params (14 CNOTs; merged
                                       local-SU(2) layers L_k).

All figures share one renderer (wire spacing, font, gate/label style). Output: 300-dpi
white-background PNGs under figures/ with tight bounding boxes.

Run:  python scripts/generate_paper_figures.py
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pennylane as qml

FIG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")
os.makedirs(FIG, exist_ok=True)

# ---- unified style ---------------------------------------------------------
DPI = 300
WIRE_DY = 1.0          # vertical spacing between wires
COL_DX = 1.0           # horizontal spacing between columns
BOX_W, BOX_H = 0.62, 0.62
FS_GATE = 11
FS_WIRE = 12
FS_NOTE = 10.5
LW = 1.4
NOISE_FILL = "0.86"    # light gray for noise-channel boxes (minimal color)


# ---------------------------------------------------------------------------
# Column packing + renderer
# ---------------------------------------------------------------------------
def _span(item, nw):
    tp = item["type"]
    if tp == "box":
        return [item["wire"]]
    if tp == "cnot":
        lo, hi = sorted((item["ctrl"], item["targ"]))
        return list(range(lo, hi + 1))
    if tp == "cswap":
        ws = [item["ctrl"], item["t1"], item["t2"]]
        return list(range(min(ws), max(ws) + 1))
    if tp == "noise":
        ws = item["wires"]
        return list(range(min(ws), max(ws) + 1))
    if tp == "noisebox":
        ws = item["wires"]
        return list(range(min(ws), max(ws) + 1))
    if tp == "full":
        return list(range(nw))
    raise ValueError(tp)


def _assign_columns(items, nw):
    free = [0] * nw
    for it in items:
        sp = _span(it, nw)
        col = max(free[w] for w in sp)
        it["col"] = col
        for w in sp:
            free[w] = col + 1
    ncols = max((it["col"] for it in items), default=0) + 1
    return ncols


def _y(row, nw):
    # row 0 (ancilla / first wire) drawn at the TOP
    return (nw - 1 - row) * WIRE_DY


def draw_circuit(items, wire_labels, outname, title=None, notes=None,
                 left_labels=None, right_labels=None, brace_groups=None,
                 col_scale=1.0):
    """items: list of dicts (see _span). Renders and saves figures/<outname>."""
    nw = len(wire_labels)
    ncols = _assign_columns(items, nw)
    dx = COL_DX * col_scale

    x0 = 0.0
    xr = ncols * dx
    fig_w = max(5.0, 1.15 * (xr) + 2.4)
    fig_h = 0.95 * nw + 1.4
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=DPI)
    ax.set_axis_off()

    def cx(col):
        return x0 + (col + 0.5) * dx

    # wires
    xL, xR = -0.35 * dx, xr + 0.15 * dx
    for r in range(nw):
        y = _y(r, nw)
        ax.plot([xL, xR], [y, y], color="k", lw=1.0, zorder=1)
        ax.text(xL - 0.18 * dx, y, wire_labels[r], ha="right", va="center",
                fontsize=FS_WIRE)
        if left_labels and left_labels[r]:
            ax.text(xL - 0.02 * dx, y + 0.28, left_labels[r], ha="right",
                    va="bottom", fontsize=FS_NOTE, color="0.25")

    def box(col, row, label, fill="white", edge="k", dashed=False):
        x, y = cx(col), _y(row, nw)
        ax.add_patch(Rectangle((x - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H,
                     facecolor=fill, edgecolor=edge, lw=LW, zorder=3,
                     linestyle="--" if dashed else "-"))
        ax.text(x, y, label, ha="center", va="center", fontsize=FS_GATE, zorder=4)

    for it in items:
        tp, col = it["type"], it["col"]
        x = cx(col)
        if tp == "box":
            box(col, it["wire"], it["label"])
        elif tp == "full":
            ytop, ybot = _y(0, nw), _y(nw - 1, nw)
            ax.add_patch(Rectangle((x - BOX_W / 2, ybot - BOX_H / 2),
                         BOX_W, (ytop - ybot) + BOX_H,
                         facecolor="white", edgecolor="k", lw=LW, zorder=3))
            ax.text(x, (ytop + ybot) / 2, it["label"], ha="center", va="center",
                    fontsize=FS_GATE + 1, zorder=4)
        elif tp == "cnot":
            yc, yt = _y(it["ctrl"], nw), _y(it["targ"], nw)
            ax.plot([x, x], [yc, yt], color="k", lw=LW, zorder=2)
            ax.add_patch(Circle((x, yc), 0.085, color="k", zorder=4))
            ax.add_patch(Circle((x, yt), 0.185, facecolor="white", edgecolor="k",
                         lw=LW, zorder=4))
            ax.plot([x - 0.185, x + 0.185], [yt, yt], color="k", lw=LW, zorder=5)
            ax.plot([x, x], [yt - 0.185, yt + 0.185], color="k", lw=LW, zorder=5)
        elif tp == "cswap":
            yc = _y(it["ctrl"], nw)
            y1, y2 = _y(it["t1"], nw), _y(it["t2"], nw)
            ax.plot([x, x], [min(yc, y1, y2), max(yc, y1, y2)], color="k", lw=LW, zorder=2)
            ax.add_patch(Circle((x, yc), 0.085, color="k", zorder=4))
            for yy in (y1, y2):
                d = 0.16
                ax.plot([x - d, x + d], [yy - d, yy + d], color="k", lw=LW, zorder=5)
                ax.plot([x - d, x + d], [yy + d, yy - d], color="k", lw=LW, zorder=5)
        elif tp == "noisebox":
            # ONE joint multi-qubit channel: a single dashed light-gray rectangle
            # covering all its wires together (not per-wire boxes).
            ws = sorted(it["wires"])
            yhi, ylo = _y(ws[0], nw), _y(ws[-1], nw)
            bw = it.get("box_w", BOX_W)
            ax.add_patch(Rectangle((x - bw / 2, ylo - BOX_H / 2), bw,
                         (yhi - ylo) + BOX_H, facecolor=NOISE_FILL, edgecolor="k",
                         lw=LW, linestyle="--", zorder=3))
            ax.text(x, yhi + BOX_H / 2 + 0.34, it["label"], ha="center", va="bottom",
                    fontsize=FS_GATE, zorder=4)
        elif tp == "noise":
            ws = sorted(it["wires"])
            ys = [_y(w, nw) for w in ws]
            ax.plot([x, x], [min(ys), max(ys)], color="0.35", lw=1.0,
                    ls="--", zorder=2)
            for w in ws:
                box(col, w, it.get("wire_label", "D"), fill=NOISE_FILL, dashed=True)
            ax.text(x, max(ys) + 0.52, it["label"], ha="center", va="bottom",
                    fontsize=FS_NOTE, color="0.15")

    # right end-state labels
    if right_labels:
        for r in range(nw):
            if right_labels[r]:
                ax.text(xR + 0.10 * dx, _y(r, nw), right_labels[r], ha="left",
                        va="center", fontsize=FS_NOTE, color="0.25")

    # brace groups: list of (row_lo, row_hi, text) drawn on the far left
    if brace_groups:
        xb = xL - 1.05 * dx
        for lo, hi, text in brace_groups:
            ytop, ybot = _y(lo, nw) + 0.32, _y(hi, nw) - 0.32
            ax.plot([xb, xb], [ybot, ytop], color="0.3", lw=1.2)
            ax.plot([xb, xb + 0.10 * dx], [ytop, ytop], color="0.3", lw=1.2)
            ax.plot([xb, xb + 0.10 * dx], [ybot, ybot], color="0.3", lw=1.2)
            ax.text(xb - 0.10 * dx, (ytop + ybot) / 2, text, ha="right", va="center",
                    fontsize=FS_NOTE, color="0.2", rotation=90)

    ax.set_xlim(xL - 1.7 * dx, xR + 1.9 * dx)
    nlines = len(notes) if notes else 0
    ax.set_ylim(_y(nw - 1, nw) - 1.0 - 0.44 * max(1, nlines), _y(0, nw) + 1.35)
    ax.set_aspect("equal")
    ax.patch.set_visible(False)

    # Center title + notes over the ACTUAL drawn circuit content (not the wire span),
    # so they are visually centered left-right in the saved (tight-bbox) image.
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
    xc = (xL + xR) / 2
    if xs:
        x0d = inv.transform((min(xs), 0))[0]
        x1d = inv.transform((max(xs), 0))[0]
        xc = 0.5 * (x0d + x1d)

    if title:
        ax.text(xc, _y(0, nw) + 1.08, title, ha="center", va="bottom",
                fontsize=FS_WIRE + 1)
    if notes:
        ybase = _y(nw - 1, nw) - 0.72
        for i, nline in enumerate(notes):
            ax.text(xc, ybase - i * 0.44, nline, ha="center", va="top",
                    fontsize=FS_NOTE, color="0.1")

    out = os.path.join(FIG, outname)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    from PIL import Image
    w, h = Image.open(out).size
    print(f"  wrote {outname}  ({w}x{h}px, {DPI} dpi)")
    return out


# ===========================================================================
# Figure 1 — Setup: Bell-isotropic input
# ===========================================================================
def fig_setup():
    items = [
        {"type": "box", "wire": 0, "label": "H"},
        {"type": "cnot", "ctrl": 0, "targ": 1},
        # ONE joint two-qubit channel acting on q0 and q1 together (not two 1-qubit boxes)
        {"type": "noisebox", "wires": [0, 1], "box_w": 0.92,
         "label": r"$D_\varepsilon^{(q_0,q_1)}$"},
    ]
    draw_circuit(
        items, [r"$q_0$", r"$q_1$"], "setup_bell_input.png",
        title="Setup: Bell-isotropic input preparation",
        left_labels=[r"$|0\rangle$", r"$|0\rangle$"],
        right_labels=[r"$\rho_t$", ""],
        notes=[r"$D_\varepsilon^{(q_0,q_1)}$: one joint two-qubit replacement depolarizing"
               r" channel,  $D_\varepsilon^{(q_0,q_1)}(\rho)=(1-\varepsilon)\rho"
               r"+\varepsilon\,(I_4/4)\,\mathrm{Tr}(\rho)$",
               r"$\Rightarrow\ \rho_t=(1-\varepsilon)\,\Phi+\varepsilon\,I_4/4,"
               r"\quad t=1-\varepsilon.$"],
        col_scale=1.35,
    )


# ===========================================================================
# Figure 2 — Step 1: ideal PQEC
# ===========================================================================
WIRES5 = [r"$a$", r"$A_1$", r"$A_2$", r"$B_1$", r"$B_2$"]  # rows 0..4


def fig_step1():
    items = [
        {"type": "box", "wire": 0, "label": "H"},
        {"type": "cswap", "ctrl": 0, "t1": 1, "t2": 3},
        {"type": "cswap", "ctrl": 0, "t1": 2, "t2": 4},
        {"type": "box", "wire": 0, "label": "H"},
    ]
    draw_circuit(
        items, WIRES5, "step1_ideal_pqec_circuit.png",
        title=r"Step 1: ideal PQEC   $U_{\mathrm{PQEC}}=H_a\,\mathrm{CSWAP}_{a;A_1B_1}\,"
              r"\mathrm{CSWAP}_{a;A_2B_2}\,H_a$",
        left_labels=[r"$|0\rangle$", "", "", "", ""],
        brace_groups=[(1, 2, r"copy 1: $\rho_t$"), (3, 4, r"copy 2: $\rho_t$")],
        notes=[r"Register $A=(A_1,A_2)$ and $B=(B_1,B_2)$ each hold one Bell-state"
               r" copy $\rho_t$; the ancilla controls the swaps $A_1\!\leftrightarrow\!B_1$"
               r" and $A_2\!\leftrightarrow\!B_2$."],
        col_scale=1.5,
    )


# ===========================================================================
# Figure 3 — Step 2: Fredkin-level global noise
# ===========================================================================
def fig_step2():
    items = [
        {"type": "box", "wire": 0, "label": "H"},
        {"type": "cswap", "ctrl": 0, "t1": 1, "t2": 3},
        {"type": "noise", "wires": [0, 1, 3], "label": r"$G_g^{(a,A_1,B_1)}$"},
        {"type": "cswap", "ctrl": 0, "t1": 2, "t2": 4},
        {"type": "noise", "wires": [0, 2, 4], "label": r"$G_g^{(a,A_2,B_2)}$"},
        {"type": "box", "wire": 0, "label": "H"},
    ]
    draw_circuit(
        items, WIRES5, "step2_fredkin_noise_circuit.png",
        title="Step 2: Fredkin-level global replacement depolarizing noise",
        brace_groups=[(1, 2, r"copy 1: $\rho_t$"), (3, 4, r"copy 2: $\rho_t$")],
        notes=[r"$G_g^{(ijk)}$: three-qubit replacement depolarizing channel after each"
               r" Fredkin — the first acts on exactly $(a,A_1,B_1)$, the second on"
               r" exactly $(a,A_2,B_2)$."],
        col_scale=1.5,
    )


# ===========================================================================
# Figure 4 — Step 3: textbook 16-CNOT (rendered from the actual _fred replay)
# ===========================================================================
def _tape_ops(qfunc):
    with qml.queuing.AnnotatedQueue() as q:
        qfunc()
    return qml.tape.QuantumScript.from_queue(q).operations


def _step3_qfunc():
    import verify_analytic_decomposed as v3
    qml.Hadamard(0)
    v3._fred(0, 1, 3, 0.0)     # CSWAP(a;A1,B1), textbook 8-CNOT decomposition
    v3._fred(0, 2, 4, 0.0)     # CSWAP(a;A2,B2)
    qml.Hadamard(0)


_LABELMAP = {"Hadamard": "H", "T": "T", "Adjoint(T)": r"$T^\dagger$"}


def _ops_to_items(ops):
    items, ncx, cx_seq = [], 0, []
    for op in ops:
        nm = op.name
        if nm == "CNOT":
            c, t = int(op.wires[0]), int(op.wires[1])
            items.append({"type": "cnot", "ctrl": c, "targ": t})
            ncx += 1
            cx_seq.append((c, t))
        elif nm in ("U3",):
            items.append({"type": "box", "wire": int(op.wires[0]), "label": r"$U_3$"})
        else:
            lab = _LABELMAP.get(nm, nm)
            items.append({"type": "box", "wire": int(op.wires[0]), "label": lab})
    return items, ncx, cx_seq


def fig_step3():
    ops = _tape_ops(_step3_qfunc)
    items, ncx, cx_seq = _ops_to_items(ops)
    assert ncx == 16, f"Step 3 must have 16 CNOTs, got {ncx}"
    draw_circuit(
        items, WIRES5, "step3_textbook_16cnot_circuit.png",
        title=r"Step 3: textbook two-Fredkin decomposition (16 CNOTs)",
        notes=[r"Textbook Clifford+$T$ Toffoli (6 CNOTs) $\Rightarrow$ 8 CNOTs per Fredkin,"
               r" 16 total; single-qubit gates ideal.",
               r"A two-qubit replacement depolarizing channel $D_q$ is applied after"
               r" every CNOT (16 in total)."],
        col_scale=0.92,
    )
    return ncx, cx_seq


# ===========================================================================
# Figure 5 — Step 4: resynthesized 14-CNOT (rendered from GATES verbatim)
# ===========================================================================
def fig_step4():
    import pqec_resynth_noise as r4
    items, ncx, cx_seq = [], 0, []
    for name, wires, params in r4.GATES:
        if name == "u":
            items.append({"type": "box", "wire": int(wires[0]), "label": r"$U_3$"})
        elif name == "cx":
            c, t = int(wires[0]), int(wires[1])
            items.append({"type": "cnot", "ctrl": c, "targ": t})
            ncx += 1
            cx_seq.append((c, t))
    assert ncx == 14, f"Step 4 must have 14 CNOTs, got {ncx}"
    draw_circuit(
        items, WIRES5, "step4_resynthesized_14cnot_circuit.png",
        title=r"Step 4: resynthesized 14-CNOT circuit ($U_4=e^{i\phi_4}U_{\mathrm{PQEC}}$)",
        notes=[r"Exact-unitary resynthesis (Qiskit, basis $\{u,\mathrm{cx}\}$): 14 CNOTs,"
               r" single-qubit gates $U_3$; same wire order and convention as Fig. 4.",
               r"A two-qubit replacement depolarizing channel $D_q$ is applied after"
               r" every CNOT (14 in total)."],
        col_scale=0.92,
    )
    return ncx, cx_seq


# ===========================================================================
# Figure 6 — Step 5: learned 14-CNOT (merged local-SU(2) layers)
# ===========================================================================
def fig_step5():
    import json
    import draw_pqc_5abc as d5
    mask = json.load(open(os.path.join(os.path.dirname(FIG), "pqc_ring_pruned.json")))["mask"]
    params = np.load(os.path.join(os.path.dirname(FIG), "pqc_ring_pruned_params.npy"))
    blocks, cnots = d5.merged_circuit(mask, params)
    ncx = len(cnots)
    assert ncx == 14, f"Step 5 must have 14 CNOTs, got {ncx}"
    assert len(blocks) == 15, f"expected 15 local layers, got {len(blocks)}"

    items = [{"type": "full", "label": r"$L_{0}$"}]
    for k, (c, t) in enumerate(cnots, start=1):
        items.append({"type": "cnot", "ctrl": int(c), "targ": int(t)})
        items.append({"type": "full", "label": rf"$L_{{{k}}}$"})
    draw_circuit(
        items, WIRES5, "step5_learned_14cnot_circuit.png",
        title=r"Step 5: learned 14-CNOT circuit   "
              r"$V_5=L_{14}C_{14}L_{13}\cdots L_1 C_1 L_0$",
        notes=[r"$L_k=\bigotimes_j U_{(k,j)}$ is a merged single-qubit (local SU(2)) layer;"
               r" $C_1,\dots,C_{14}$ are the learned CNOTs (exact full unitary).",
               r"A two-qubit replacement depolarizing channel $D_q$ is applied after"
               r" every CNOT (14 in total)."],
        col_scale=1.0,
    )
    return ncx, [(int(c), int(t)) for c, t in cnots]


# ===========================================================================
# Figure 7 — threshold comparison (values supplied, not recomputed)
# ===========================================================================
THRESHOLD_TABLE = {
    "eps":    [0.10, 0.20, 0.30, 0.40, 0.50, 0.60],
    "Step 3": [0.0330, 0.0610, 0.0850, 0.1030, 0.1170, 0.1260],
    "Step 4": [0.0413, 0.0788, 0.1119, 0.1399, 0.1621, 0.1780],
    "Step 5": [0.0602, 0.1203, 0.1797, 0.2371, 0.2908, 0.3384],
}


def fig_threshold():
    t = THRESHOLD_TABLE
    fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=DPI)
    styles = [("Step 3: textbook 16 CNOT", "-^", "C2"),
              ("Step 4: resynthesized 14 CNOT", "-o", "C0"),
              ("Step 5: learned 14 CNOT", "-s", "C3")]
    for (label, mk, col), key in zip(styles, ["Step 3", "Step 4", "Step 5"]):
        ax.plot(t["eps"], t[key], mk, color=col, lw=1.8, ms=6, label=label)
    ax.set_xlabel(r"input noise  $\varepsilon$", fontsize=FS_WIRE)
    ax.set_ylabel(r"CNOT-noise threshold  $q_{\mathrm{th}}$", fontsize=FS_WIRE)
    ax.set_title("Operational threshold vs input noise", fontsize=FS_WIRE + 1)
    ax.grid(True, ls=":", lw=0.6, color="0.7")
    ax.legend(frameon=False, fontsize=FS_NOTE, loc="upper left")
    ax.set_xlim(0.08, 0.62)
    ax.set_ylim(0, 0.37)
    out = os.path.join(FIG, "threshold_comparison_steps3_5.png")
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    from PIL import Image
    w, h = Image.open(out).size
    print(f"  wrote threshold_comparison_steps3_5.png  ({w}x{h}px, {DPI} dpi)")


def main():
    print("Generating paper figures into figures/ ...")
    fig_setup()
    fig_step1()
    fig_step2()
    n3, seq3 = fig_step3()
    n4, seq4 = fig_step4()
    n5, seq5 = fig_step5()
    fig_threshold()
    print("\nVerification (rendered == simulated gate data):")
    print(f"  Step 3 CNOT count = {n3}  (expect 16)")
    print(f"  Step 4 CNOT count = {n4}  (expect 14);  sequence = {seq4}")
    print(f"  Step 5 CNOT count = {n5}  (expect 14);  sequence = {seq5}")
    print("\nDone.")


if __name__ == "__main__":
    main()
