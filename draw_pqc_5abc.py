"""
Draw the learned Step-5 circuit (14 CNOTs) and the auxiliary observable circuit
(5 CNOTs) and emit an analytic circuit spec.
=====================================================================================

The ansatz (`pqc_ring_prune.ansatz_masked`) always carries a full RX-RY-RZ layer
before the first CNOT and one after every CNOT slot; pruning only deletes CNOTs, so a
pruned solution has several *consecutive* rotation layers with no CNOT between them.
For a clean, analysis-ready picture we compose each maximal run of single-qubit
rotations (between two CNOTs, or before the first / after the last) into a single net
SU(2) per wire and render it as one `Rot(phi, theta, omega)` (= RZ.RY.RZ) gate. The
drawn circuit is verified to equal the stored solution to machine precision (up to a
global phase).

Outputs:
  circuit_pqc_5a.png   -- the 14-CNOT circuit (compiles U exactly; the isometry uses the SAME circuit)
  circuit_pqc_5c.png   -- the auxiliary 5-CNOT circuit (reproduces the ancilla-parity observable)
  PQC_CIRCUITS_FOR_ANALYSIS.md -- exact ordered gate lists + merged per-block SU(2)
                                  matrices, for an independent analytic check.

Run:  python draw_pqc_5abc.py
"""
import io
import json

import numpy as np
import pennylane as qml
import matplotlib.pyplot as plt

from pqc_ring_prune import ansatz_masked, SEQ
from pqc_ring_ansatz import unitary, _rmat
from pqc_common import DIM, U_TARGET

WIRE_LABELS = {0: "a (anc)", 1: "A1", 2: "A2", 3: "B1", 4: "B2"}


# ---------------------------------------------------------------------------
# merge the ansatz into [block] (CNOT [block])*  with one net SU(2) per wire
# ---------------------------------------------------------------------------
def merged_circuit(mask, params):
    """Return (blocks, cnots):
       blocks -- list of dicts {wire: 2x2 net unitary}, length = #CNOT + 1
       cnots  -- list of (control, target), length = #CNOT
    The k-th block is applied before the k-th CNOT; the last block is applied at the end.
    """
    ops = ansatz_masked(mask)
    blocks, cnots = [], []
    cur = {w: np.eye(2, dtype=complex) for w in range(5)}
    for op in ops:
        if op[0] == "g":
            _, kind, w, p = op
            cur[w] = _rmat(kind, params[p]) @ cur[w]        # left-apply, same as unitary()
        else:                                               # CNOT closes the current block
            blocks.append(cur)
            cnots.append((op[1], op[2]))
            cur = {w: np.eye(2, dtype=complex) for w in range(5)}
    blocks.append(cur)                                      # trailing rotation block
    return blocks, cnots


def _zyz(U):
    """(phi, theta, omega) with U = e^{i.} RZ(omega) RY(theta) RZ(phi) (qml.Rot)."""
    det = U[0, 0] * U[1, 1] - U[0, 1] * U[1, 0]
    V = U / np.sqrt(det)                                    # to SU(2)
    theta = 2.0 * np.arctan2(abs(V[1, 0]), abs(V[0, 0]))
    a00, a10 = np.angle(V[0, 0]), np.angle(V[1, 0])
    omega = a10 - a00
    phi = -a10 - a00
    return float(phi), float(theta), float(omega)


def _build_qnode(blocks, cnots):
    def circuit():
        for k, cn in enumerate(cnots):
            for w in range(5):
                qml.Rot(*_zyz(blocks[k][w]), wires=w)
            qml.CNOT(wires=list(cn))
        for w in range(5):                                  # trailing block
            qml.Rot(*_zyz(blocks[-1][w]), wires=w)
        return qml.state()
    return circuit


def _verify(blocks, cnots, mask, params):
    """|Tr(V_drawn^dag V_ansatz)|/32 -- 1.0 means the drawn circuit == the solution."""
    Vans = unitary(ansatz_masked(mask), params)
    Vdrawn = qml.matrix(_build_qnode(blocks, cnots), wire_order=range(5))()
    return abs(np.vdot(Vdrawn, Vans)) / DIM


def draw(mask, params, fname, title):
    blocks, cnots = merged_circuit(mask, params)
    ov = _verify(blocks, cnots, mask, params)
    fig, ax = qml.draw_mpl(_build_qnode(blocks, cnots), wire_order=range(5),
                           show_all_wires=True, style="pennylane")()
    for i, lab in WIRE_LABELS.items():
        ax.text(-1.6, i, lab, ha="right", va="center", fontsize=11)
    ax.set_title(f"{title}\n{len(cnots)} CNOTs, {len(cnots)+1} rotation blocks "
                 f"(each Rot = RZ·RY·RZ);  |Tr(V†V_sol)|/32 = {ov:.10f}", fontsize=12)
    fig.savefig(fname, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fname}  (drawn==solution overlap {ov:.12f})")
    return blocks, cnots, ov


# ---------------------------------------------------------------------------
# analytic spec (exact primitive gate list + merged SU(2) blocks) for GPT
# ---------------------------------------------------------------------------
def _build_raw_qnode(mask, params, barriers=False):
    """QNode of the EXACT ansatz: RX/RY/RZ per wire + CNOTs, no merging.

    With `barriers=True` a visual-only barrier is placed on both sides of every
    CNOT, so the drawing reads as the alternating structure it really is:
    rotation layer L_0 | CNOT C_1 | L_1 | C_2 | ...  The barriers carry
    `only_visual=True` and therefore change nothing about the circuit.
    """
    ops = ansatz_masked(mask)
    gate = {"rx": qml.RX, "ry": qml.RY, "rz": qml.RZ}

    def circuit():
        prev = None
        for op in ops:
            kind_now = "g" if op[0] == "g" else "cx"
            # one barrier at each rotation-layer <-> CNOT boundary, never two in a row
            if barriers and prev is not None and kind_now != prev:
                qml.Barrier(wires=range(5), only_visual=True)
            if op[0] == "g":
                _, kind, w, p = op
                gate[kind](params[p], wires=w)
            else:
                qml.CNOT(wires=[op[1], op[2]])
            prev = kind_now
        return qml.state()
    return circuit


def draw_raw(mask, params, fname, title, max_length=30, dpi=110):
    """Draw the exact RX-RY-RZ + CNOT ansatz (all rotation layers kept).

    The primitive form is ~120 columns wide, so it is wrapped into several rows
    with `max_length`.  PennyLane returns one (fig, ax) per row; we render each
    row and stack them into a single image, which keeps the drawer's own
    rendering (`style="pennylane"`, matching the Step-3 / Step-4 figures).
    """
    ncx = sum(mask)
    rows = qml.draw_mpl(_build_raw_qnode(mask, params, barriers=True),
                        wire_order=range(5), show_all_wires=True,
                        style="pennylane", max_length=max_length)()
    if not isinstance(rows, list):
        rows = [rows]

    imgs = []
    for k, (f, a) in enumerate(rows):
        for i, lab in WIRE_LABELS.items():
            a.text(-2.4, i, lab, ha="right", va="center", fontsize=11)
        a.text(-2.4, -0.9, f"part {k + 1} of {len(rows)}", ha="right", va="center",
               fontsize=10, style="italic", color="0.35")
        buf = io.BytesIO()
        f.savefig(buf, format="png", dpi=dpi, bbox_inches="tight",
                  facecolor="white")
        plt.close(f)
        buf.seek(0)
        imgs.append(plt.imread(buf))

    gap = 24                                        # px between rows
    W = max(im.shape[1] for im in imgs)
    head = 300                                      # px reserved for the title
    H = head + sum(im.shape[0] for im in imgs) + gap * (len(imgs) - 1)
    fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi, facecolor="white")
    y = H - head
    for im in imgs:
        h, w = im.shape[0], im.shape[1]
        ax = fig.add_axes([0.0, (y - h) / H, w / W, h / H])
        ax.imshow(im)
        ax.axis("off")
        y -= h + gap
    fig.text(0.5, 1 - 0.28 * head / H,
             f"{title}\n{ncx} CNOTs; exact ansatz = initial RX-RY-RZ layer + a full "
             f"RX-RY-RZ layer after every CNOT slot (19 rotation layers, 285 params).\n"
             f"Barriers are visual only: they separate each rotation layer from the "
             f"CNOT that follows it.  Read the parts top to bottom.",
             ha="center", va="center", fontsize=22, linespacing=1.6)
    fig.savefig(fname, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {fname}  ({len(rows)} rows)")


def _primitive_lines(mask, params):
    """Exact ordered primitive gates with numeric angles."""
    lines, layer = [], 0
    ops = ansatz_masked(mask)
    for op in ops:
        if op[0] == "g":
            _, kind, w, p = op
            lines.append(f"    {kind.upper()}(theta={params[p]:+.10f}, wire={w})")
        else:
            lines.append(f"  CNOT(control={op[1]}, target={op[2]})")
    return ops, lines


def _fmt2x2(U):
    def c(z):
        return f"{z.real:+.6f}{z.imag:+.6f}j"
    return (f"[[{c(U[0,0])}, {c(U[0,1])}],\n"
            f"       [{c(U[1,0])}, {c(U[1,1])}]]")


def write_spec(specs):
    L = []
    L.append("# Step-5 learned circuits — exact spec for analytic analysis\n")
    L.append("Convention: 5 qubits, wire 0 = ancilla `a`, wires 1,2 = kept register A "
             "(A1,A2), wires 3,4 = discarded register B (B1,B2).  Single-qubit gates "
             "are ideal; the noise model puts a 2-qubit depolarizing channel "
             "`(1-e2)rho + e2 I/4` after **each CNOT**.  Rotation convention "
             "`RX(t)=exp(-i t X/2)` etc.; `Rot(phi,theta,omega) = RZ(omega) RY(theta) "
             "RZ(phi)`.  Input `rho_eps = (1-eps)|Phi+><Phi+| + eps I/4` on register "
             "(A,B) pairs (A1B1 and A2B2 each carry one noisy Bell copy).\n")
    L.append("Two representations are given per circuit: (1) the **exact primitive gate "
             "list** as compiled (RX/RY/RZ per wire, CNOTs), fully reproducible; (2) the "
             "**merged form** — consecutive single-qubit layers composed into one net "
             "SU(2) per wire between CNOTs — which is the minimal-depth circuit and the "
             "convenient object for a pen-and-paper conjugation of the observables "
             "through the CNOTs.  Both implement the identical unitary.\n")
    for s in specs:
        L.append("\n" + "=" * 78)
        L.append(f"## {s['name']}")
        L.append("=" * 78 + "\n")
        L.append(s["desc"] + "\n")
        L.append(f"- CNOT sequence (in order): {s['cnots']}")
        L.append(f"- CNOT count: {len(s['cnots'])}    rotation blocks: {len(s['cnots'])+1}")
        L.append(f"- verification: |Tr(U_target^dag V)|/32 = {s['tgt_overlap']}")
        if s.get("obs_note"):
            L.append(f"- {s['obs_note']}")
        L.append(f"- drawn(merged)==solution overlap = {s['drawn_overlap']:.12f}\n")

        L.append("### (1) Exact primitive gate list\n```")
        L.extend(s["prim_lines"])
        L.append("```\n")

        L.append("### (2) Merged form: net SU(2) per wire per block, alternating with CNOTs\n")
        blocks, cnots = s["blocks"], s["cnots"]
        for k in range(len(cnots) + 1):
            L.append(f"**Block {k}** (before "
                     + (f"CNOT {cnots[k]}" if k < len(cnots) else "measurement") + "):")
            L.append("```")
            for w in range(5):
                U = blocks[k][w]
                phi, th, om = _zyz(U)
                L.append(f"  wire {w} ({WIRE_LABELS[w]:>7}): "
                         f"Rot(phi={phi:+.6f}, theta={th:+.6f}, omega={om:+.6f})")
                L.append(f"      U = {_fmt2x2(U)}")
            L.append("```")
            if k < len(cnots):
                L.append(f"then **CNOT(control={cnots[k][0]}, target={cnots[k][1]})** "
                         f"followed by depolarizing `(1-e2)rho+e2 I/4` on wires "
                         f"{list(cnots[k])}.\n")
    open("PQC_CIRCUITS_FOR_ANALYSIS.md", "w").write("\n".join(L))
    print("  wrote PQC_CIRCUITS_FOR_ANALYSIS.md")


def main():
    print("Drawing Step-5 learned circuits + writing analytic spec\n")

    # ---- Step 5 (unitary; the isometry uses the same circuit) : the 14-CNOT circuit ----
    mask_a = json.load(open("pqc_ring_pruned.json"))["mask"]
    par_a = np.load("pqc_ring_pruned_params.npy")
    b_a, c_a, ov_a = draw(mask_a, par_a, "circuit_pqc_5a.png",
                          "Step 5 — learned 14-CNOT circuit "
                          "(compiles U = H_a·CSWAP·CSWAP·H_a exactly)")
    draw_raw(mask_a, par_a, "circuit_pqc_5a_raw.png",
             "Step 5 — learned 14-CNOT circuit, primitive RX-RY-RZ + CNOT form")
    tgt_a = abs(np.vdot(U_TARGET, unitary(ansatz_masked(mask_a), par_a))) / DIM
    prim_ops_a, prim_a = _primitive_lines(mask_a, par_a)

    # ---- Auxiliary (observable relaxation) : the 5-CNOT circuit ----
    mask_c = json.load(open("pqc_ring_5c.json"))["mask"]
    par_c = np.load("pqc_ring_5c_params.npy")
    b_c, c_c, ov_c = draw(mask_c, par_c, "circuit_pqc_5c.png",
                          "Auxiliary (observable relaxation) — learned 5-CNOT circuit "
                          "(reproduces the purified observable F=<Z_a⊗O>/<Z_a>)")
    tgt_c = abs(np.vdot(U_TARGET, unitary(ansatz_masked(mask_c), par_c))) / DIM
    prim_ops_c, prim_c = _primitive_lines(mask_c, par_c)

    write_spec([
        dict(name="Step 5 — 14-CNOT full-unitary circuit",
             desc="This single circuit compiles the full 5-qubit gadget unitary "
                  "U = H_a . CSWAP(0;1,3) . CSWAP(0;2,4) . H_a to machine precision, and "
                  "therefore also realizes the ancilla-|0> isometry exactly. It is "
                  "the greedy-pruning floor of the gadget-matched ansatz (13 CNOTs "
                  "unreachable). Under per-CNOT depolarizing noise it defines the "
                  "operational threshold analysed in pqc_ring_threshold.py.",
             cnots=c_a, blocks=b_a, prim_lines=prim_a,
             tgt_overlap=f"{tgt_a:.12f}  (=1 => exact compilation of U)",
             drawn_overlap=ov_a),
        dict(name="Auxiliary (observable relaxation) — 5-CNOT observable circuit",
             desc="Pruned from the 14-CNOT circuit under the OBSERVABLE cost: it need "
                  "only reproduce the ancilla-parity correlators <Z_a> -> Tr(rho^2) and "
                  "<Z_a (x) O> -> Tr(O rho^2) for O in {|Phi+><Phi+|, ZZ} over an eps "
                  "grid (NOT the full unitary). This is the observable rung of the relaxation "
                  "ladder; 4 CNOTs is unreachable.",
             cnots=c_c, blocks=b_c, prim_lines=prim_c,
             tgt_overlap=f"{tgt_c:.12f}  (NOT ~1: it is not the full unitary, by design)",
             obs_note="observable read-out F = <Z_a (x) O> / <Z_a>, with the two register "
                      "copies measured jointly; O acts on kept register A = wires (1,2).",
             drawn_overlap=ov_c),
    ])
    print("\nDone.")


if __name__ == "__main__":
    main()
