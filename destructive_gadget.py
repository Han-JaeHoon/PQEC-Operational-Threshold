"""
Destructive / virtual-distillation SWAP gadget: same PQEC observable, 2 CNOTs.
=============================================================================

The controlled-SWAP gadget (H_a . CSWAP . CSWAP . H_a) needs, for a 2-qubit Bell
register, two decomposed Fredkins = 16 CNOTs.  But the quantity PQEC actually uses
is the purified observable

    F = <O>_purified = Tr(O rho^2) / Tr(rho^2)  (for two identical copies rho).

That can be obtained *without* an ancilla or any controlled-SWAP, via the
destructive SWAP test / virtual-distillation measurement
(Garcia-Escartin & Chamorro-Posada 2013; Huggins et al. 2021): apply a Bell-basis
change between corresponding qubits of the two copies and read out.  For an
M-qubit register this uses only **M CNOTs** (here M = 2) instead of 16.

Construction (registers A=(A1,A2), B=(B1,B2), qubits (0,1,2,3), no ancilla):

  SWAP_reg = SWAP_{A1B1} SWAP_{A2B2}
  V        = [H_{A1} CNOT_{A1->B1}] [H_{A2} CNOT_{A2->B2}]        (Bell-basis change)
  Pi       = V SWAP_reg V^dagger        (diagonalized SWAP; a Z-basis observable)

  denominator observable  O_den = Pi                       -> <O_den> = Tr(rho^2)
  numerator observable    O_num = V [ 1/2 (O_A SWAP_reg + SWAP_reg O_A) ] V^dagger
                                                            -> <O_num> = Tr(O rho^2)

The experiment applies the (noisy) Bell-change V and measures the FIXED
ideal-frame observables O_den, O_num; F = <O_num>/<O_den>.

Results (verified in main()):
  * Ideal (no gate noise): F equals the controlled-SWAP gadget and Tr(Orho^2)/Tr(rho^2)
    to ~1e-15, for the fidelity projector and for a generic Pauli observable.
  * With a two-qubit depolarizing channel eps2 after each CNOT (single-qubit gates
    ideal, same convention as the controlled-SWAP study), the CNOT-noise threshold
    eps2* is ~3-4x higher than the 16-CNOT controlled-SWAP gadget -- because only
    2 CNOTs carry noise.

Caveat: the destructive gadget is measurement-only; it yields <O> but not a
coherent purified state to feed forward.  It is the right tool for the
observable/threshold analysis, not for interleaving a purified state in an algorithm.

Run:  python destructive_gadget.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

from noisy_bell_state import make_noisy_bell, O_PHI_PLUS, global_depol_kraus
from pqec_gadget import obs_purified                     # controlled-SWAP gadget (ideal)
from pqec_cnot_threshold import F_dec as F_cswap, F_bare, eps2_star as eps2_star_cswap

np.set_printoptions(precision=4, suppress=True)

# ---------------------------------------------------------------------------
# 4-qubit operators, ordering (A1,A2,B1,B2) = qubits (0,1,2,3), qubit 0 = MSB
# ---------------------------------------------------------------------------
_I2 = np.eye(2, dtype=complex)
_H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)


def _embed1(g, i, n=4):
    ops = [_I2] * n
    ops[i] = g
    M = ops[0]
    for o in ops[1:]:
        M = np.kron(M, o)
    return M


def _cnot4(c, t, n=4):
    dim = 2 ** n
    M = np.zeros((dim, dim), dtype=complex)
    for x in range(dim):
        b = [(x >> (n - 1 - k)) & 1 for k in range(n)]
        if b[c]:
            b[t] ^= 1
        M[sum(v << (n - 1 - k) for k, v in enumerate(b)), x] = 1
    return M


def _swap4(i, j, n=4):
    dim = 2 ** n
    M = np.zeros((dim, dim), dtype=complex)
    for x in range(dim):
        b = [(x >> (n - 1 - k)) & 1 for k in range(n)]
        b[i], b[j] = b[j], b[i]
        M[sum(v << (n - 1 - k) for k, v in enumerate(b)), x] = 1
    return M


SWAP_REG = _swap4(0, 2) @ _swap4(1, 3)                       # SWAP_{A1B1} SWAP_{A2B2}
V_BELL = _embed1(_H, 0) @ _cnot4(0, 2) @ _embed1(_H, 1) @ _cnot4(1, 3)   # Bell-change
PI = V_BELL @ SWAP_REG @ V_BELL.conj().T                     # diagonalized SWAP


def _observables(O2):
    """Fixed (ideal-frame) numerator and denominator observables for a 2-qubit
    register observable O2 acting on A=(A1,A2)."""
    O_A = np.kron(O2, np.eye(4, dtype=complex))
    O_num = V_BELL @ (0.5 * (O_A @ SWAP_REG + SWAP_REG @ O_A)) @ V_BELL.conj().T
    return O_num, PI


# ---------------------------------------------------------------------------
# Genuine destructive gadget circuit (no ancilla; 2 CNOTs; noise after each CNOT)
# ---------------------------------------------------------------------------
_dev = qml.device("default.mixed", wires=4)


def _dest_state(rho_AB, eps2):
    @qml.qnode(_dev)
    def c():
        qml.QubitDensityMatrix(rho_AB, wires=[0, 1, 2, 3])
        qml.CNOT(wires=[0, 2])
        if eps2 > 0:
            qml.QubitChannel(global_depol_kraus(eps2), wires=[0, 2])
        qml.Hadamard(0)
        qml.CNOT(wires=[1, 3])
        if eps2 > 0:
            qml.QubitChannel(global_depol_kraus(eps2), wires=[1, 3])
        qml.Hadamard(1)
        return qml.density_matrix(wires=[0, 1, 2, 3])
    return c()


def F_dest(eps, eps2, O2=O_PHI_PLUS):
    """Purified observable from the destructive gadget with 2-qubit depol eps2 after
    each CNOT (single-qubit gates ideal)."""
    rho = make_noisy_bell(eps)
    s = _dest_state(np.kron(rho, rho), eps2)
    O_num, O_den = _observables(O2)
    return float((np.trace(O_num @ s) / np.trace(O_den @ s)).real)


def threshold_dest(eps, hi=0.98):
    t = 1 - eps
    fb = F_bare(t)
    f = lambda e: F_dest(eps, e) - fb
    if f(1e-9) <= 0:
        return 0.0
    lo, hi_ = 0.0, hi
    for _ in range(60):
        m = 0.5 * (lo + hi_)
        lo, hi_ = (m, hi_) if f(m) > 0 else (lo, m)
    return 0.5 * (lo + hi_)


# ===========================================================================
def main():
    print("=" * 78)
    print(" Destructive / VD gadget (2 CNOTs) vs controlled-SWAP gadget (16 CNOTs)")
    print("=" * 78)

    # (0) IDEAL equivalence: destructive == controlled-SWAP == exact trace
    print("\n (0) ideal equivalence  F = Tr(O rho^2)/Tr(rho^2)  (eps2 = 0):")
    ZZ = np.kron(_Z, _Z)
    worst = 0.0
    for eps in [0.2, 0.4, 0.6]:
        rho = make_noisy_bell(eps)
        # fidelity projector O = |Phi+><Phi+|
        fd = F_dest(eps, 0.0, O_PHI_PLUS)
        fc = obs_purified(rho)
        fe = float((np.trace(O_PHI_PLUS @ rho @ rho) / np.trace(rho @ rho)).real)
        # generic Pauli O = ZZ
        gd = F_dest(eps, 0.0, ZZ)
        ge = float((np.trace(ZZ @ rho @ rho) / np.trace(rho @ rho)).real)
        worst = max(worst, abs(fd - fc), abs(fd - fe), abs(gd - ge))
        print(f"   eps={eps}:  Phi+  dest={fd:.8f} cSWAP={fc:.8f} exact={fe:.8f} | "
              f"ZZ dest={gd:.8f} exact={ge:.8f}")
    print(f"   max deviation = {worst:.2e}   "
          f"{'PASS' if worst < 1e-12 else 'FAIL'}")

    # (1) CNOT-noise threshold comparison
    print("\n (1) CNOT-noise threshold eps2*  (single-qubit gates ideal):")
    print(f"   {'eps':>5} {'dest (2 CNOT)':>14} {'cSWAP (16 CNOT)':>16} {'ratio':>7}")
    rows = []
    for eps in [0.10, 0.20, 0.30, 0.40, 0.50, 0.60]:
        d = threshold_dest(eps)
        c = eps2_star_cswap(1 - eps)
        rows.append((eps, d, c))
        print(f"   {eps:>5.2f} {d:>14.4f} {c:>16.4f} {d / c:>7.2f}")
    print("   -> the destructive gadget tolerates ~3-4x higher per-CNOT noise")
    print("      (2 noisy CNOTs instead of 16).")

    # ---- Figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))

    ax = axes[0]
    e2s = np.linspace(0, 0.45, 46)
    for eps, c in [(0.2, "C0"), (0.4, "C1"), (0.6, "C2")]:
        ax.plot(e2s, [F_dest(eps, e) for e in e2s], "-", color=c, lw=2,
                label=f"dest, $\\varepsilon$={eps}")
        ax.plot(e2s, [F_cswap(e, 1 - eps) for e in e2s], "--", color=c, lw=1.3,
                label=f"cSWAP, $\\varepsilon$={eps}")
        ax.axhline(F_bare(1 - eps), color=c, ls=":", lw=0.8)
    ax.set_xlabel(r"per-CNOT depolarizing  $\varepsilon_2$")
    ax.set_ylabel(r"$F_{PQEC}$")
    ax.set_title("(a) Destructive (solid) vs controlled-SWAP (dashed)")
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="lower left")

    ax = axes[1]
    es = np.linspace(0.05, 0.65, 25)
    ax.plot(es, [threshold_dest(e) for e in es], "-o", color="C0", ms=3,
            label="destructive (2 CNOT)")
    ax.plot(es, [eps2_star_cswap(1 - e) for e in es], "-s", color="C3", ms=3,
            label="controlled-SWAP (16 CNOT)")
    ax.axhspan(0, 0.01, color="0.85", alpha=0.6)
    ax.text(0.33, 0.006, "hardware CNOT error $\\sim10^{-2}$", fontsize=8, va="center")
    ax.set_xlabel(r"input noise  $\varepsilon$")
    ax.set_ylabel(r"CNOT-noise threshold  $\varepsilon_2^*$")
    ax.set_title("(b) Threshold: fewer noisy CNOTs -> higher tolerance")
    ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig("destructive_gadget.png", dpi=140)
    print("\n  saved  destructive_gadget.png")


if __name__ == "__main__":
    main()
