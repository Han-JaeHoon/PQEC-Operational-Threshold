"""
Step 3b -- Realistic gate noise on a DECOMPOSED controlled-SWAP.
===============================================================

Decompose each Fredkin into native 1- and 2-qubit gates and put realistic
depolarizing noise on *each native gate*.

Representative decomposition (textbook):
    CSWAP(c; x, y) = CNOT(x->y) . Toffoli(c, y; x) . CNOT(x->y),
    Toffoli = the standard Clifford+T circuit (Nielsen & Chuang, Fig. 4.9):
              6 CNOTs + 2 H + 7 T/T^dagger.
So one Fredkin = 8 CNOTs (2 outer + 6 in the Toffoli) + 9 single-qubit gates.
(5-two-qubit-gate optimum: Smolin & DiVincenzo, PRA 53, 2855 (1996); a 7-CNOT
connectivity-aware version: Cruz & Murta, arXiv:2305.18128 / APL Quantum 2024.)

Noise model (realistic):
  * after every CNOT:            2-qubit depolarizing p2  ((1-p2)rho + p2 I4/4)
  * after every 1-qubit gate:    1-qubit depolarizing p1  (PennyLane convention)
  * default p1 = p2/10 (typical hardware 2q:1q error ratio ~ 10:1)

Read out with the parity correlator  <O>_PQEC = <Z_a (x) O> / <Z_a (x) I>.

ORIENTATION (matters!).  The Toffoli TARGET leg carries all the H/T/T^dagger
single-qubit gates, so it absorbs far more single-qubit noise than the CONTROL
leg.  The same Fredkin unitary can put that target on either swapped qubit:
  * orient="discard": Toffoli target on the DISCARDED register B -> the kept
    register A only sees CNOTs -> milder single-qubit slope (K1 = 2).
  * orient="retain":  Toffoli target on the RETAINED register A -> A absorbs the
    single-qubit noise -> steeper slope (K1 = 5/2).
CNOT noise is symmetric across the swap, so it (and its slope K2 = 17/8) is
orientation-independent.  See verify_analytic_decomposed.py for the exact
analytic match of the "retain" orientation.

Unlike the 3-qubit GLOBAL depolarizing benchmark (which cancels in the ratio and
gives no threshold), native-gate noise hits the data qubits asymmetrically, so a
finite operational threshold p2* appears.

Run:  python pqec_decomposed_noise.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

from noisy_bell_state import make_noisy_bell, O_PHI_PLUS, global_depol_kraus
from pqec_gadget import obs_purified

np.set_printoptions(precision=4, suppress=True)


# ---------------------------------------------------------------------------
# Noisy native gates
# ---------------------------------------------------------------------------
def _cnot(c, t, p2):
    qml.CNOT(wires=[c, t])
    if p2 > 0:
        qml.QubitChannel(global_depol_kraus(p2), wires=[c, t])   # 2-qubit depol


def _g1(gate, w, p1, adjoint=False):
    (qml.adjoint(gate)(wires=w) if adjoint else gate(wires=w))
    if p1 > 0:
        qml.DepolarizingChannel(p1, wires=w)                     # 1-qubit depol


# ---------------------------------------------------------------------------
# Toffoli (controls c1,c2; target t) -- Clifford+T, 6 CNOTs  (N&C Fig. 4.9)
# ---------------------------------------------------------------------------
def _toffoli(c1, c2, t, p1, p2):
    _g1(qml.Hadamard, t, p1)
    _cnot(c2, t, p2); _g1(qml.T, t, p1, adjoint=True)
    _cnot(c1, t, p2); _g1(qml.T, t, p1)
    _cnot(c2, t, p2); _g1(qml.T, t, p1, adjoint=True)
    _cnot(c1, t, p2); _g1(qml.T, t, p1); _g1(qml.T, c2, p1)
    _cnot(c1, c2, p2); _g1(qml.Hadamard, t, p1)
    _g1(qml.T, c1, p1); _g1(qml.T, c2, p1, adjoint=True)
    _cnot(c1, c2, p2)


# ---------------------------------------------------------------------------
# Decomposed Fredkin (controlled-SWAP):  CNOT(b,a) . Toffoli(q,a;b) . CNOT(b,a)
# swaps qubits a,b with control q; the Toffoli TARGET is b.
# ---------------------------------------------------------------------------
def _fredkin(q, a, b, p1, p2):
    _cnot(b, a, p2)
    _toffoli(q, a, b, p1, p2)
    _cnot(b, a, p2)


# ---------------------------------------------------------------------------
# Gadget with two decomposed, noisy Fredkins
#   wire 0 = ancilla, [1,2] = register A (kept), [3,4] = register B (discarded)
#   orient="discard": Toffoli target on B (wires 3,4)  -> kept register protected
#   orient="retain":  Toffoli target on A (wires 1,2)
# ---------------------------------------------------------------------------
_dev = qml.device("default.mixed", wires=5)


@qml.qnode(_dev)
def _gadget_decomposed(rho_AB, p1, p2, O, orient):
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
    _g1(qml.Hadamard, 0, p1)
    if orient == "discard":
        _fredkin(0, 1, 3, p1, p2)      # target = wire 3 (register B)
        _fredkin(0, 2, 4, p1, p2)      # target = wire 4 (register B)
    else:                              # "retain"
        _fredkin(0, 3, 1, p1, p2)      # target = wire 1 (register A)
        _fredkin(0, 4, 2, p1, p2)      # target = wire 2 (register A)
    _g1(qml.Hadamard, 0, p1)
    return (qml.expval(qml.PauliZ(0) @ qml.Hermitian(O, wires=[1, 2])),
            qml.expval(qml.PauliZ(0)))


def obs_pqec_decomposed(eps, p2, p1=None, orient="discard", O=O_PHI_PLUS):
    """Purified observable with realistic native-gate noise (2q depol p2, 1q depol
    p1; default p1 = p2/10) on the decomposed Fredkins."""
    if p1 is None:
        p1 = p2 / 10
    rho = make_noisy_bell(eps)
    zO, zI = _gadget_decomposed(np.kron(rho, rho), p1, p2, O, orient)
    return float(zO) / float(zI)


def no_qec(eps):
    return 1 - 3 * eps / 4


def threshold_p2(eps, p1_ratio=0.1, orient="discard", hi=0.20, n=41):
    """Find p2* where obs_pqec_decomposed crosses no_qec(eps) (PQEC stops helping)."""
    p2s = np.linspace(0, hi, n)
    base = no_qec(eps)
    prev_p, prev_v = 0.0, obs_pqec_decomposed(eps, 0.0, orient=orient)
    for p2 in p2s[1:]:
        v = obs_pqec_decomposed(eps, p2, p1=p1_ratio * p2, orient=orient)
        if v < base:
            return prev_p + (p2 - prev_p) * (prev_v - base) / (prev_v - v)
        prev_p, prev_v = p2, v
    return None


# ===========================================================================
def main():
    print("=" * 78)
    print(" Step 3b -- realistic gate noise on the decomposed controlled-SWAP")
    print("=" * 78)

    # -- sanity: p1=p2=0 reproduces the ideal gadget (both orientations) -----
    worst = 0.0
    for orient in ("discard", "retain"):
        for eps in [0.1, 0.3, 0.5, 0.8]:
            a = obs_pqec_decomposed(eps, 0.0, 0.0, orient)
            b = obs_purified(make_noisy_bell(eps))
            worst = max(worst, abs(a - b))
    print(f"\n  sanity: decomposed(p1=p2=0) == ideal gadget, max err = {worst:.2e}  "
          f"{'PASS' if worst < 1e-12 else 'FAIL'}")

    # -- threshold vs input noise, BOTH orientations (p1 = p2/10) ------------
    print("\n" + "-" * 78)
    print(" Operational threshold p2* (p1 = p2/10), by Fredkin orientation")
    print("-" * 78)
    print(f"  {'eps':>5} {'p2* (target on discard)':>24} {'p2* (target on retain)':>24}")
    for e in [0.2, 0.3, 0.4, 0.5, 0.6]:
        pd = threshold_p2(e, orient="discard")
        pr = threshold_p2(e, orient="retain")
        print(f"  {e:>5.2f} {pd:>24.4f} {pr:>24.4f}")
    print("  -> 'discard' (Toffoli target on the discarded register) tolerates a bit"
          " more\n     gate noise: the kept register is shielded from the single-qubit"
          " gates.")

    # -- isolate the orientation effect: single-qubit vs CNOT noise ----------
    print("\n" + "-" * 78)
    print(" Orientation effect is a SINGLE-QUBIT-noise effect (eps=0.40):")
    print("-" * 78)
    print(f"  {'':>14}{'<O> vs p1 (p2=0)':>26}{'<O> vs p2 (p1=0)':>26}")
    print(f"  {'noise':>6}{'discard':>10}{'retain':>10}{'  ':>6}{'discard':>10}{'retain':>10}")
    for x in [0.0, 0.05, 0.10, 0.20]:
        od1 = obs_pqec_decomposed(0.40, 0.0, p1=x, orient="discard")
        or1 = obs_pqec_decomposed(0.40, 0.0, p1=x, orient="retain")
        od2 = obs_pqec_decomposed(0.40, x, p1=0.0, orient="discard")
        or2 = obs_pqec_decomposed(0.40, x, p1=0.0, orient="retain")
        print(f"  {x:>6.2f}{od1:>10.4f}{or1:>10.4f}{'  ':>6}{od2:>10.4f}{or2:>10.4f}")
    print("  -> single-qubit noise (p1): orientations DIFFER (K1 = 2 vs 5/2).")
    print("     CNOT noise (p2):        orientations COINCIDE (K2 = 17/8).")

    # -- Figure -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))

    ax = axes[0]
    xs = np.linspace(0, 0.20, 21)
    for orient, c in [("discard", "C0"), ("retain", "C3")]:
        v1 = [obs_pqec_decomposed(0.40, 0.0, p1=x, orient=orient) for x in xs]
        ax.plot(xs, v1, "-", color=c, lw=2, label=f"1q noise, {orient}")
        v2 = [obs_pqec_decomposed(0.40, x, p1=0.0, orient=orient) for x in xs]
        ax.plot(xs, v2, "--", color=c, lw=1.5, label=f"2q noise, {orient}")
    ax.axhline(no_qec(0.40), color="k", ls=":", lw=1)
    ax.text(0.201, no_qec(0.40), " no-QEC", va="center", fontsize=8)
    ax.set_xlabel(r"gate noise strength")
    ax.set_ylabel(r"$\langle O\rangle_{PQEC}$  ($\varepsilon=0.4$)")
    ax.set_title("(a) Orientation splits 1q noise (solid),\nnot 2q noise (dashed)")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    es = np.linspace(0.1, 0.65, 23)
    for orient, c in [("discard", "C0"), ("retain", "C3")]:
        ps = [threshold_p2(e, orient=orient, hi=0.20) for e in es]
        es2 = [e for e, p in zip(es, ps) if p is not None]
        ps2 = [p for p in ps if p is not None]
        ax.plot(es2, ps2, "-o", color=c, ms=3, label=f"target on {orient}")
    ax.set_xlabel(r"input noise  $\varepsilon$")
    ax.set_ylabel(r"gate-error threshold  $p_2^*$   ($p_1=p_2/10$)")
    ax.set_title("(b) Threshold $p_2^*$ vs input noise")
    ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig("pqec_decomposed_threshold.png", dpi=140)
    print("\n  saved  pqec_decomposed_threshold.png")


if __name__ == "__main__":
    main()
