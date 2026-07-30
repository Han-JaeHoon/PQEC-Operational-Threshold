"""
CNOT-noise threshold of the LEARNED 14-CNOT gadget (from pruning).
=================================================================

Take the pruned 14-CNOT exact compilation (pqc_ring_pruned_*.{npy,json}), put a 2-qubit
depolarizing eps2 after each of its 14 CNOTs (single-qubit rotations ideal -- same
convention as Step 3b/4a/4b), and find the threshold eps2*(eps) where the purified
fidelity F = <Z_a (x) O>/<Z_a> drops to F_bare = (1+3t)/4.

Compare four 14/16/2-CNOT gadgets under the identical noise convention:
  * textbook controlled-SWAP  (16)  -- pqec_cnot_threshold.eps2_star
  * Step 4a, Qiskit-optimized (14)  -- pqec_resynth_noise.threshold_resynth
  * Step 5 learned + pruned   (14)  -- computed here
  * destructive Step 4b       ( 2)  -- destructive_gadget.threshold_dest_closed

Run:  python pqc_ring_threshold.py
"""
import json
import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

from pqc_ring_prune import ansatz_masked, SEQ
from pqc_common import U_TARGET, DIM, F_exact
from noisy_bell_state import rho_eps_analytic, O_PHI_PLUS, global_depol_kraus
from pqec_cnot_threshold import eps2_star as eps2_star_textbook
from pqec_resynth_noise import threshold_resynth as threshold_step4a
from destructive_gadget import threshold_dest_closed

MASK = json.load(open("pqc_ring_pruned.json"))["mask"]
THETA = np.load("pqc_ring_pruned_params.npy")
OPS = ansatz_masked(MASK)
N_CX = sum(1 for op in OPS if op[0] == "cnot")

_dev = qml.device("default.mixed", wires=5)


def _apply(eps2):
    for op in OPS:
        if op[0] == "g":
            kind, w, p = op[1], op[2], op[3]
            (qml.RX if kind == "rx" else qml.RY if kind == "ry" else qml.RZ)(THETA[p], wires=w)
        else:
            qml.CNOT(wires=[op[1], op[2]])
            if eps2 > 0:
                qml.QubitChannel(global_depol_kraus(eps2), wires=[op[1], op[2]])


@qml.qnode(_dev)
def _readout(rho_AB, eps2, O):
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
    _apply(eps2)
    return (qml.expval(qml.PauliZ(0) @ qml.Hermitian(O, wires=[1, 2])),
            qml.expval(qml.PauliZ(0)))


def F_learned(eps, eps2, O=O_PHI_PLUS):
    r = rho_eps_analytic(eps)
    zO, zI = _readout(np.kron(r, r), eps2, O)
    return float(zO) / float(zI)


def F_bare(eps):
    return (1 + 3 * (1 - eps)) / 4


def threshold_learned(eps, hi=0.6):
    fb = F_bare(eps)
    if F_learned(eps, 0.0) < fb:
        return 0.0
    if F_learned(eps, hi) >= fb:
        return hi
    lo, h = 0.0, hi
    for _ in range(50):
        m = 0.5 * (lo + h)
        lo, h = (m, h) if F_learned(eps, m) >= fb else (lo, m)
    return 0.5 * (lo + h)


def main():
    print("=" * 80)
    print(f" CNOT-noise threshold of the LEARNED pruned {N_CX}-CNOT gadget")
    print("=" * 80)

    # (0) exactness at eps2=0
    err = max(abs(F_learned(e, 0.0) - F_exact(e)) for e in [0.1, 0.2, 0.4, 0.6])
    print(f"\n (0) eps2=0 read-out == F_exact:  max err = {err:.2e}  (14-CNOT circuit is exact)")

    # (1) monotonicity
    up = -1.0
    for eps in [0.2, 0.4, 0.6]:
        Fs = [F_learned(eps, x) for x in np.linspace(0, 0.5, 30)]
        up = max(up, max(Fs[i+1]-Fs[i] for i in range(len(Fs)-1)))
    print(f" (1) F(eps2) max upward step = {up:+.4f}  (<=0 monotone)")

    # (2) threshold comparison
    print("\n (2) CNOT-noise threshold eps2* (single-qubit gates ideal):")
    print(f"     {'eps':>5} | {'textbook(16)':>12} {'Step4a(14)':>11} {'learned(14)':>12} "
          f"{'dest(2)':>9}")
    EPS = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60]
    rows = []
    for e in EPS:
        t16 = eps2_star_textbook(1 - e)
        t4a = threshold_step4a(e)
        tL = threshold_learned(e)
        t2 = threshold_dest_closed(1 - e)
        rows.append((e, t16, t4a, tL, t2))
        print(f"     {e:>5.2f} | {t16:>12.4f} {t4a:>11.4f} {tL:>12.4f} {t2:>9.4f}")

    print("\n  Both 14-CNOT gadgets (Step 4a Qiskit-optimized vs Step 5 learned+pruned)")
    print("  keep the exact unitary; comparing their thresholds shows whether the")
    print("  specific 14-CNOT layout matters for noise robustness.")

    # figure
    fig, ax = plt.subplots(figsize=(7, 5))
    E = [r[0] for r in rows]
    ax.plot(E, [r[4] for r in rows], "-d", color="C3", label="Step 4b destructive (2)")
    ax.plot(E, [r[3] for r in rows], "-o", color="C1", label=f"Step 5 learned+pruned ({N_CX})")
    ax.plot(E, [r[2] for r in rows], "-s", color="C0", label="Step 4a Qiskit-opt (14)")
    ax.plot(E, [r[1] for r in rows], "-^", color="C2", label="textbook cSWAP (16)")
    ax.set_xlabel(r"input noise  $\varepsilon$")
    ax.set_ylabel(r"CNOT-noise threshold  $\varepsilon_2^*$")
    ax.set_title("Learned 14-CNOT gadget: CNOT-noise threshold")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig("pqc_ring_threshold.png", dpi=140)
    print("\n  saved  pqc_ring_threshold.png")


if __name__ == "__main__":
    main()
