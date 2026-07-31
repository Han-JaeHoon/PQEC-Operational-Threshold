"""
Step 4 — CNOT-noise threshold of the optimized 14-CNOT gadget (PennyLane).
==========================================================================

Step 4 reduced the SWAP-test gadget  H_a . CSWAP(0;1,3) . CSWAP(0;2,4) . H_a  from 16
to 14 CNOTs while keeping the *exact* unitary (resynthesize_gadget.py, via Qiskit).
That established unitary equivalence but NOT the operational threshold: fewer CNOTs
usually help, but the threshold also depends on where the CNOTs sit and how errors
propagate, which the peephole optimizer did not optimize for.  Here we settle it.

The Qiskit-optimized 14-CNOT gate list (basis {u, cx}, optimization_level 3, seed 7)
is hard-coded below as `GATES`, so this analysis is pure PennyLane and reproducible
without Qiskit.  We (1) verify the PennyLane replay is unitarily IDENTICAL to the
gadget U (Hilbert-Schmidt overlap = 1 up to global phase), then (2) put a 2-qubit
depolarizing channel of strength eps2 after EACH of the 14 CNOTs (single-qubit gates
ideal -- the same CNOT-only convention as Step 3 / the destructive gadget) and find the
threshold eps2*(eps) where the purified fidelity  F = <Z_a (x) O> / <Z_a>  drops to
F_bare = (1+3t)/4.

We compare three gadgets under the identical noise convention:
  * textbook controlled-SWAP  (16 CNOTs)  -- eps2* from pqec_cnot_threshold
  * optimized Step 4          (14 CNOTs)  -- computed here
  * destructive gadget         ( 2 CNOTs)  -- threshold_dest_closed  (auxiliary)

Run:  python pqec_resynth_noise.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

from noisy_bell_state import rho_eps_analytic, O_PHI_PLUS, global_depol_kraus
from pqc_common import U_TARGET, DIM
from pqec_cnot_threshold import eps2_star as eps2_star_textbook
from destructive_gadget import threshold_dest_closed

np.set_printoptions(precision=4, suppress=True)

# Qiskit-optimized 14-CNOT gadget (basis {u, cx}, opt_level 3, seed_transpiler 7);
# global phase pi/4 (irrelevant to observables). qml.U3(theta,phi,lambda) == qiskit u.
GATES = [
    ('u', [0], [1.5707963267948966, 0.0, 3.141592653589793]),
    ('u', [1], [1.5707963267948966, -1.5707963267948966, 1.5707963267948966]),
    ('u', [2], [1.5707963267948966, -1.5707963267948966, 1.5707963267948966]),
    ('u', [3], [1.5707963267948968, -3.141592653589793, -2.356194490192345]),
    ('cx', [1, 3], []),
    ('u', [1], [1.5707963267948966, 1.5707963267948966, 1.5707963267948966]),
    ('u', [3], [0.7853981633974475, -2.356194490192345, -1.5707963267948966]),
    ('cx', [0, 3], []),
    ('u', [3], [0.0, 0.0, 0.7853981633974483]),
    ('cx', [1, 3], []),
    ('u', [1], [0.0, 0.0, 0.7853981633974483]),
    ('u', [3], [0.0, 0.0, -0.7853981633974483]),
    ('cx', [0, 3], []),
    ('cx', [0, 1], []),
    ('u', [0], [0.0, 0.0, 0.7853981633974483]),
    ('u', [1], [0.0, 0.0, -0.7853981633974483]),
    ('cx', [0, 1], []),
    ('u', [3], [1.5707963267948966, 0.0, -2.3561944901923453]),
    ('cx', [3, 1], []),
    ('u', [4], [1.5707963267948968, -3.141592653589793, -2.356194490192345]),
    ('cx', [2, 4], []),
    ('u', [2], [1.5707963267948966, 1.5707963267948966, 1.5707963267948966]),
    ('u', [4], [0.7853981633974475, -2.356194490192345, -1.5707963267948966]),
    ('cx', [0, 4], []),
    ('u', [4], [0.0, 0.0, 0.7853981633974483]),
    ('cx', [2, 4], []),
    ('u', [2], [0.0, 0.0, 0.7853981633974483]),
    ('u', [4], [0.0, 0.0, -0.7853981633974483]),
    ('cx', [0, 4], []),
    ('cx', [0, 2], []),
    ('u', [0], [0.0, 0.0, 0.7853981633974483]),
    ('u', [2], [0.0, 0.0, -0.7853981633974483]),
    ('cx', [0, 2], []),
    ('u', [0], [1.5707963267948966, 0.0, 3.141592653589793]),
    ('u', [4], [1.5707963267948966, 0.0, -2.3561944901923453]),
    ('cx', [4, 2], []),
]
N_CX = sum(1 for g in GATES if g[0] == "cx")


def _apply(eps2):
    """Replay the 14-CNOT gadget; 2-qubit depol(eps2) after each CNOT."""
    for name, q, p in GATES:
        if name == "u":
            qml.U3(p[0], p[1], p[2], wires=q[0])
        else:
            qml.CNOT(wires=q)
            if eps2 > 0:
                qml.QubitChannel(global_depol_kraus(eps2), wires=q)


_devq = qml.device("default.qubit", wires=5)
_devm = qml.device("default.mixed", wires=5)


@qml.qnode(_devq)
def _unitary_state():
    _apply(0.0)
    return qml.state()


def verify_unitary():
    """Hilbert-Schmidt overlap of the PennyLane replay with the gadget U (=1 up to
    global phase iff identical)."""
    V = qml.matrix(_apply, wire_order=range(5))(0.0)
    return abs(np.vdot(U_TARGET, V)) / DIM        # |Tr(U^dag V)|/32 == 1 iff equal


@qml.qnode(_devm)
def _readout(rho_AB, eps2, O):
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
    _apply(eps2)
    return (qml.expval(qml.PauliZ(0) @ qml.Hermitian(O, wires=[1, 2])),
            qml.expval(qml.PauliZ(0)))


def F_resynth(eps, eps2, O=O_PHI_PLUS):
    r = rho_eps_analytic(eps)
    zO, zI = _readout(np.kron(r, r), eps2, O)
    return float(zO) / float(zI)


def F_bare(eps):
    return (1 + 3 * (1 - eps)) / 4


def F_exact(eps, O=O_PHI_PLUS):
    r = rho_eps_analytic(eps)
    r2 = r @ r
    return float(np.real(np.trace(O @ r2) / np.trace(r2)))


def threshold_resynth(eps, hi=0.6):
    fb = F_bare(eps)
    if F_resynth(eps, 0.0) < fb:
        return 0.0
    if F_resynth(eps, hi) >= fb:
        return hi
    lo, h = 0.0, hi
    for _ in range(50):
        m = 0.5 * (lo + h)
        if F_resynth(eps, m) >= fb:
            lo = m
        else:
            h = m
    return 0.5 * (lo + h)


# ===========================================================================
def main():
    print("=" * 80)
    print(f" Step 4 -- CNOT-noise threshold of the optimized {N_CX}-CNOT gadget (PennyLane)")
    print("=" * 80)

    # (0) the PennyLane replay is unitarily identical to U
    ov = verify_unitary()
    print(f"\n (0) |Tr(U^dag V)|/32 = {ov:.12f}   (=1 => identical unitary up to phase)")

    # (1) at eps2=0 the read-out reproduces the exact purified observable
    O_ZZ = np.kron(np.diag([1., -1]), np.diag([1., -1]))
    err0 = max(abs(F_resynth(e, 0.0) - F_exact(e)) for e in [0.1, 0.2, 0.4, 0.6])
    errZZ = max(abs(F_resynth(e, 0.0, O_ZZ) - F_exact(e, O_ZZ)) for e in [0.2, 0.4, 0.6])
    print(f" (1) eps2=0 read-out == F_exact:  Phi+ err {err0:.2e},  ZZ err {errZZ:.2e}")

    # (2) monotonicity of F(eps2) (assumption behind the threshold bisection)
    up = -1.0
    for eps in [0.2, 0.4, 0.6]:
        Fs = [F_resynth(eps, x) for x in np.linspace(0, 0.5, 40)]
        up = max(up, max(Fs[i + 1] - Fs[i] for i in range(len(Fs) - 1)))
    print(f" (2) F(eps2) max upward step over eps in (0.2,0.4,0.6): {up:+.5f}  (<=0 monotone)")

    # (3) threshold table: Step 4 (14) vs textbook (16) vs destructive (2)
    print("\n (3) CNOT-noise threshold eps2* (single-qubit gates ideal):")
    print(f"     {'eps':>5} | {'textbook(16)':>12} {'Step4(14)':>11} {'dest(2)':>10} | "
          f"{'14/16':>6}")
    EPS = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60]
    rows = []
    for e in EPS:
        t16 = eps2_star_textbook(1 - e)
        t14 = threshold_resynth(e)
        t2 = threshold_dest_closed(1 - e)
        rows.append((e, t16, t14, t2))
        ratio = t14 / t16 if t16 > 0 else float("nan")
        print(f"     {e:>5.2f} | {t16:>12.4f} {t14:>11.4f} {t2:>10.4f} | {ratio:>6.2f}")

    print("\n  Reading it:")
    print("   * Step 4 (14) keeps the EXACT gadget unitary, so it is a coherent")
    print("     purified-state gadget (unlike the measurement-only destructive gadget).")
    print("   * Whether 14 CNOTs beat 16 for the threshold is now computed, not assumed.")

    # ---- figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

    ax = axes[0]
    e0 = 0.40
    e2s = np.linspace(0, 0.35, 60)
    ax.plot(e2s, [F_resynth(e0, x) for x in e2s], "-", color="C0", lw=2,
            label="Step 4 (14 CNOT)")
    ax.axhline(F_bare(e0), color="0.5", ls=":", lw=1)
    ax.text(0.0, F_bare(e0) + .003, "no-QEC $F_{bare}$", fontsize=8)
    ax.axhline(F_exact(e0), color="0.5", ls="--", lw=1)
    ax.text(0.0, F_exact(e0) + .003, "ideal $F_{exact}$", fontsize=8)
    tr = threshold_resynth(e0)
    ax.plot(tr, F_bare(e0), "o", color="C0", ms=8)
    ax.text(tr, F_bare(e0) - .03, f" $\\varepsilon_2^*$={tr:.3f}", fontsize=8, va="top")
    ax.set_xlabel(r"per-CNOT depolarizing  $\varepsilon_2$")
    ax.set_ylabel(r"purified fidelity  $F$  ($\varepsilon=0.4$)")
    ax.set_title(f"(a) optimized {N_CX}-CNOT gadget: $F$ vs CNOT noise")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    E = [r[0] for r in rows]
    ax.plot(E, [r[3] for r in rows], "-d", color="C3", label="destructive (2, auxiliary)")
    ax.plot(E, [r[2] for r in rows], "-o", color="C0", label=f"Step 4: optimized ({N_CX})")
    ax.plot(E, [r[1] for r in rows], "-^", color="C2", label="textbook cSWAP (16)")
    ax.set_xlabel(r"input noise  $\varepsilon$")
    ax.set_ylabel(r"CNOT-noise threshold  $\varepsilon_2^*$")
    ax.set_title("(b) threshold vs input noise (identical noise convention)")
    ax.legend(frameon=False, fontsize=9, loc="upper left")

    fig.tight_layout()
    fig.savefig("pqec_resynth_threshold.png", dpi=140)
    print("\n  saved  pqec_resynth_threshold.png")


if __name__ == "__main__":
    main()
