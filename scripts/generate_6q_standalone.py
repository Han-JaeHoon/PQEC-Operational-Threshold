"""
Standalone 6Q distributed (LOCC) circuit figures for the calculation note.
==========================================================================

Renders the two six-qubit distributed virtual-distillation circuits in the same
visual language as figures/step3_16cnot_circuit.* and
figures/step4_resynthesized_14cnot_circuit.png: black line art, white gate
boxes, filled control dots, small in-figure title, one or two note lines, no
caption.

  6Q-L16 (Type 3)  each party runs the textbook Clifford+T Controlled-SWAP
                   used in Step 3 (8 CNOTs), replayed from the SAME object
                   Step 3 simulates: verify_analytic_decomposed._fred.
                   16 CNOTs in total, 8 per party.

  6Q-L14 (Type 4)  each party runs the seven-CNOT block G7 extracted from the
                   SAME authoritative Step-4 gate list pqec_resynth_noise.GATES
                   (Qiskit transpile, basis {u, cx}, optimization_level 3,
                   seed 7).  14 CNOTs in total, 7 per party.

Wire convention follows the rest of this repository: A and B label the two
COPIES, the subscript labels the party.  Party 1 holds (A1, B1) and party 2
holds (A2, B2), so the local Controlled-SWAPs are

    party 1:  CSWAP(a; A1, B1)        party 2:  CSWAP(b; A2, B2)

exactly the wire triples {0,1,3} and {0,2,4} that the five-qubit Step 3 and
Step 4 circuits act on.  Replacing the single shared ancilla by one ancilla per
party is the only change; no two-qubit gate crosses the party cut.

Both blocks are drawn in parallel because they act on disjoint wires -- that is
the column packer's doing, not a layout choice.

Outputs:
  figures/6q_type3_l16_circuit.pdf  /  .png
  figures/6q_type4_l14_circuit.pdf  /  .png

Run:  python scripts/generate_6q_standalone.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import pennylane as qml

from generate_paper_figures import (
    _assign_columns, _y, DPI, COL_DX, BOX_W, BOX_H,
    FS_GATE, FS_WIRE, FS_NOTE, LW,
)
import verify_analytic_decomposed as v3
import pqec_resynth_noise as r4

FIG = os.path.join(ROOT, "figures")

# ---------------------------------------------------------------------------
# Wires.  Rows 0..2 are party 1 (ancilla a, data A1 B1); rows 3..5 are party 2
# (data A2 B2, ancilla b).  Each party's block is contiguous, so the party cut
# is a single horizontal line and no vertical gate line crosses it.
# ---------------------------------------------------------------------------
WIRES6 = [r"$a$", r"$A_1$", r"$B_1$", r"$A_2$", r"$B_2$", r"$b$"]
P1 = (0, 1, 2)          # a, A1, B1
P2 = (5, 3, 4)          # b, A2, B2
PARTY_CUT = 2           # dashed line between row 2 and row 3

_LABELMAP = {"Hadamard": "H", "T": "T", "Adjoint(T)": r"$T^\dagger$"}


# ===========================================================================
# 6Q-L16  (Type 3) -- replayed from verify_analytic_decomposed._fred
# ===========================================================================
def _l16_qfunc():
    qml.Hadamard(P1[0])
    qml.Hadamard(P2[0])
    v3._fred(P1[0], P1[1], P1[2], 0.0)     # party 1: CSWAP(a; A1, B1)
    v3._fred(P2[0], P2[1], P2[2], 0.0)     # party 2: CSWAP(b; A2, B2)


def _tape_ops(qfunc):
    with qml.queuing.AnnotatedQueue() as q:
        qfunc()
    return qml.tape.QuantumScript.from_queue(q).operations


def build_l16():
    items, ncx, cx_seq = [], 0, []
    for op in _tape_ops(_l16_qfunc):
        if op.name == "CNOT":
            c, t = int(op.wires[0]), int(op.wires[1])
            items.append({"type": "cnot", "ctrl": c, "targ": t})
            ncx += 1
            cx_seq.append((c, t))
        else:
            items.append({"type": "box", "wire": int(op.wires[0]),
                          "label": _LABELMAP.get(op.name, op.name)})
    assert ncx == 16, f"6Q-L16 must have 16 CNOTs, got {ncx}"
    for c, t in cx_seq:                     # locality is checked, not assumed
        assert {c, t} <= set(P1) or {c, t} <= set(P2), \
            f"CNOT ({c},{t}) crosses the party cut"
    assert sum(1 for c, t in cx_seq if {c, t} <= set(P1)) == 8
    return items, ncx, cx_seq


# ===========================================================================
# 6Q-L14  (Type 4) -- G7 extracted from pqec_resynth_noise.GATES
# ===========================================================================
# Documented split of the 36-operation Step-4 list (five-qubit wires
# [a, A1, A2, B1, B2] = [0, 1, 2, 3, 4]); every claim below is asserted.
H_IDX = [0, 33]                                        # the two ancilla H gates
A_IDX = [1] + list(range(3, 19))                       # block on wires {0,1,3}
B_IDX = [2] + list(range(19, 33)) + [34, 35]           # block on wires {0,2,4}
HALF_PI, PI = 1.5707963267948966, 3.141592653589793


def _g7_blocks():
    """Split GATES into the two seven-CNOT blocks and check they are the same
    gate sequence up to the names of the three wires each acts on."""
    G = r4.GATES
    assert len(G) == 36, f"Step-4 list must have 36 operations, got {len(G)}"
    assert sorted(H_IDX + A_IDX + B_IDX) == list(range(36)), "index split is not a partition"

    for i in H_IDX:                                     # the displaced Hadamards
        name, wires, params = G[i]
        assert name == "u" and tuple(int(w) for w in wires) == (0,)
        assert max(abs(p - t) for p, t in zip(params, (HALF_PI, 0.0, PI))) < 1e-12, \
            f"operation {i} is not a Hadamard on the ancilla"

    blocks = {}
    for tag, idx, support, local in (("A", A_IDX, {0, 1, 3}, {0: 0, 1: 1, 3: 2}),
                                     ("B", B_IDX, {0, 2, 4}, {0: 0, 2: 1, 4: 2})):
        ops, ncx = [], 0
        for i in idx:
            name, wires, params = G[i]
            ws = tuple(int(w) for w in wires)
            assert set(ws) <= support, f"operation {i} leaves the {tag} block support"
            ops.append((name, tuple(local[w] for w in ws),
                        None if params is None else tuple(params)))
            ncx += (name == "cx")
        assert ncx == 7, f"block {tag} must have 7 CNOTs, got {ncx}"
        blocks[tag] = ops

    assert len(blocks["A"]) == len(blocks["B"]) == 17
    for (nA, wA, pA), (nB, wB, pB) in zip(blocks["A"], blocks["B"]):
        assert nA == nB and wA == wB, "the two blocks differ in structure"
        if pA:
            assert max(abs(x - y) for x, y in zip(pA, pB)) == 0.0, \
                "the two blocks differ in a rotation angle"
    return blocks["A"]                       # the common block, on local wires 0,1,2


def build_l14():
    g7 = _g7_blocks()
    items, ncx, cx_seq = [], 0, []
    items.append({"type": "box", "wire": P1[0], "label": "H"})
    items.append({"type": "box", "wire": P2[0], "label": "H"})
    for party in (P1, P2):
        for name, lw, _params in g7:
            ws = tuple(party[w] for w in lw)
            if name == "cx":
                items.append({"type": "cnot", "ctrl": ws[0], "targ": ws[1]})
                ncx += 1
                cx_seq.append(ws)
            else:
                items.append({"type": "box", "wire": ws[0], "label": r"$U_3$"})
    assert ncx == 14, f"6Q-L14 must have 14 CNOTs, got {ncx}"
    for c, t in cx_seq:
        assert {c, t} <= set(P1) or {c, t} <= set(P2), \
            f"CNOT ({c},{t}) crosses the party cut"
    assert sum(1 for c, t in cx_seq if {c, t} <= set(P1)) == 7
    return items, ncx, cx_seq


# ===========================================================================
# read-out: the two ancilla settings
# ===========================================================================
def add_readout(items):
    """S-dagger (dashed: YY setting only) then H, on both ancillas."""
    for w in (P1[0], P2[0]):
        items.append({"type": "box", "wire": w, "label": r"$S^\dagger$", "dashed": True})
    for w in (P1[0], P2[0]):
        items.append({"type": "box", "wire": w, "label": "H"})
    return items


# ===========================================================================
# renderer -- same geometry and weights as the Step 3 / Step 4 standalones
# ===========================================================================
def render(items, title, notes, out_stem, col_scale=0.92):
    nw = len(WIRES6)
    ncols = _assign_columns(items, nw)
    dx = COL_DX * col_scale
    x0 = 0.0

    def cx(col):
        return x0 + (col + 0.5) * dx

    xL, xR = -0.35 * dx, ncols * dx + 0.15 * dx
    fig, ax = plt.subplots(figsize=(max(6.0, 1.1 * (xR - xL) + 2.0), 0.95 * nw + 1.8),
                           dpi=DPI)
    ax.set_axis_off()
    ax.patch.set_visible(False)

    # party cut
    ycut = 0.5 * (_y(PARTY_CUT, nw) + _y(PARTY_CUT + 1, nw))
    ax.plot([xL - 1.25 * dx, xR], [ycut, ycut], color="0.62", lw=1.0, ls=(0, (6, 4)),
            zorder=0)

    # wires
    for r in range(nw):
        y = _y(r, nw)
        ax.plot([xL, xR], [y, y], color="k", lw=1.0, zorder=1)
        ax.text(xL - 0.18 * dx, y, WIRES6[r], ha="right", va="center", fontsize=FS_WIRE)

    # party brackets
    xb = xL - 0.98 * dx
    for lo, hi, text in ((0, 2, "party 1"), (3, 5, "party 2")):
        ytop, ybot = _y(lo, nw) + 0.34, _y(hi, nw) - 0.34
        ax.plot([xb, xb], [ybot, ytop], color="0.3", lw=1.2)
        ax.plot([xb, xb + 0.10 * dx], [ytop, ytop], color="0.3", lw=1.2)
        ax.plot([xb, xb + 0.10 * dx], [ybot, ybot], color="0.3", lw=1.2)
        ax.text(xb - 0.12 * dx, (ytop + ybot) / 2, text, ha="center", va="center",
                fontsize=FS_NOTE, color="0.2", rotation=90)

    for it in items:
        x = cx(it["col"])
        if it["type"] == "box":
            y = _y(it["wire"], nw)
            ax.add_patch(Rectangle((x - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H,
                         facecolor="white", edgecolor="k", lw=LW, zorder=3,
                         linestyle="--" if it.get("dashed") else "-"))
            ax.text(x, y, it["label"], ha="center", va="center", fontsize=FS_GATE,
                    zorder=4)
        else:
            yc, yt = _y(it["ctrl"], nw), _y(it["targ"], nw)
            ax.plot([x, x], [yc, yt], color="k", lw=LW, zorder=2)
            ax.add_patch(Circle((x, yc), 0.085, color="k", zorder=4))
            ax.add_patch(Circle((x, yt), 0.185, facecolor="white", edgecolor="k",
                         lw=LW, zorder=4))
            ax.plot([x - 0.185, x + 0.185], [yt, yt], color="k", lw=LW, zorder=5)
            ax.plot([x, x], [yt - 0.185, yt + 0.185], color="k", lw=LW, zorder=5)

    ax.set_xlim(xL - 1.6 * dx, xR + 0.6 * dx)
    ax.set_ylim(_y(nw - 1, nw) - 0.95 - 0.44 * len(notes), _y(0, nw) + 1.3)
    ax.set_aspect("equal")

    # centre title and notes over the drawn content
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

    ax.text(xc, _y(0, nw) + 1.02, title, ha="center", va="bottom", fontsize=FS_WIRE)
    ybase = _y(nw - 1, nw) - 0.78
    for i, line in enumerate(notes):
        ax.text(xc, ybase - i * 0.44, line, ha="center", va="top", fontsize=FS_NOTE,
                color="0.1")

    outs = []
    for ext in ("pdf", "png"):
        out = os.path.join(FIG, f"{out_stem}.{ext}")
        fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
        outs.append(out)
    plt.close(fig)
    return outs


READOUT_NOTE = [
    r"Dashed $S^\dagger$: present only in the $YY$ setting.  The estimator "
    r"$\Delta(O)=\langle X_aX_b\otimes O\rangle-\langle Y_aY_b\otimes O\rangle$"
    r" requires both settings.",
]
NOISE_NOTE = (r"no two-qubit gate crosses the party cut.  A two-qubit replacement "
              r"depolarizing channel $D_q$ is applied after every CNOT.")


def main():
    i16, n16, s16 = build_l16()
    render(add_readout(i16),
           "6Q-L16 (Type 3): two local textbook Controlled-SWAPs (16 CNOTs)",
           [r"Each party runs the Step-3 textbook Clifford+$T$ Controlled-SWAP "
            r"(8 CNOTs) on its own ancilla;",
            NOISE_NOTE] + READOUT_NOTE,
           "6q_type3_l16_circuit")

    i14, n14, s14 = build_l14()
    render(add_readout(i14),
           "6Q-L14 (Type 4): two local resynthesized Fredkin blocks (14 CNOTs)",
           [r"Each party runs the seven-CNOT block $G_7=e^{-i\pi/8}F$ extracted from "
            r"the Step-4 gate list (Qiskit, basis $\{u,\mathrm{cx}\}$);",
            NOISE_NOTE] + READOUT_NOTE,
           "6q_type4_l14_circuit")

    from PIL import Image
    for stem, ncx, seq in (("6q_type3_l16_circuit", n16, s16),
                           ("6q_type4_l14_circuit", n14, s14)):
        w, h = Image.open(os.path.join(FIG, stem + ".png")).size
        print(f"wrote figures/{stem}.pdf  and  figures/{stem}.png  ({w}x{h}px, {DPI} dpi)")
        print(f"  CNOTs = {ncx}  (party 1 | party 2 = "
              f"{sum(1 for c, t in seq if {c, t} <= set(P1))} | "
              f"{sum(1 for c, t in seq if {c, t} <= set(P2))}),  cross-party = 0")
        print(f"  sequence = {seq}")


if __name__ == "__main__":
    main()
