"""
Step 3 (v1) -- Noise on the Fredkin gates of the PQEC gadget.
============================================================

Noise model requested: immediately after each controlled-SWAP (Fredkin), apply a
3-qubit GLOBAL depolarizing channel of strength g_F to the three qubits that
Fredkin acted on:

    E_gF(sigma) = (1 - g_F) sigma + g_F * I_8 / 8        (d = 8, three qubits).

Gadget wiring (see draw_pqec_gadget.py):
    wire 0 = ancilla (control),  [1,2] = register A (kept),  [3,4] = register B.
    Fredkin 1 = cSWAP(0; 1,3)  acts on qubits {0,1,3}
    Fredkin 2 = cSWAP(0; 2,4)  acts on qubits {0,2,4}
so the two depolarizing insertions are on wires [0,1,3] and [0,2,4].

Only the Fredkins are noisy here (the ancilla H's and the readout are left ideal);
those are added later.  The purified value is read out with the same parity
correlator as the ideal gadget,

    <O>_PQEC = <Z_anc (x) O> / <Z_anc (x) I>,

which reduces to Tr(O rho^2)/Tr(rho^2) when g_F = 0.

3-qubit global depolarizing Kraus operators (64 three-qubit Paulis P_i):
    K_0 = sqrt(1 - g_F + g_F/64) I,     K_i = sqrt(g_F/64) P_i  (i = 1..63),
which is trace preserving:  sum K_i^dag K_i = I.

Run:  python pqec_gadget_noise.py
"""

import numpy as np
import pennylane as qml

from noisy_bell_state import make_noisy_bell, O_PHI_PLUS, I2, X, Y, Z
from pqec_gadget import obs_purified as obs_purified_ideal

np.set_printoptions(precision=4, suppress=True)

# ---------------------------------------------------------------------------
# 3-qubit global depolarizing channel
# ---------------------------------------------------------------------------
_PAULIS_1Q = (I2, X, Y, Z)
_PAULIS_3Q = [np.kron(np.kron(a, b), c)
              for a in _PAULIS_1Q for b in _PAULIS_1Q for c in _PAULIS_1Q]  # 64


def global_depol3_kraus(g):
    """Kraus operators of the 3-qubit global depolarizing channel of strength g:
        E(sigma) = (1-g) sigma + g I_8/8.
    """
    K = [np.sqrt(1 - g + g / 64) * _PAULIS_3Q[0]]
    K += [np.sqrt(g / 64) * P for P in _PAULIS_3Q[1:]]
    return K


# ---------------------------------------------------------------------------
# Noisy gadget: 3-qubit depolarizing after each Fredkin
#   wire 0 = ancilla, [1,2] = register A (kept), [3,4] = register B (discarded)
# ---------------------------------------------------------------------------
_dev = qml.device("default.mixed", wires=5)


@qml.qnode(_dev)
def _gadget_obs_noisy(rho_AB, g_F, O):
    """SWAP-test gadget with a 3-qubit global depolarizing (strength g_F) right
    after each Fredkin.  Returns the correlators (<Z(x)O>, <Z(x)I>) on (anc, A)."""
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])   # two noisy copies rho (x) rho
    qml.Hadamard(0)
    qml.ctrl(qml.SWAP, control=0)(wires=[1, 3])          # Fredkin 1 on {0,1,3}
    if g_F > 0:
        qml.QubitChannel(global_depol3_kraus(g_F), wires=[0, 1, 3])
    qml.ctrl(qml.SWAP, control=0)(wires=[2, 4])          # Fredkin 2 on {0,2,4}
    if g_F > 0:
        qml.QubitChannel(global_depol3_kraus(g_F), wires=[0, 2, 4])
    qml.Hadamard(0)
    return (qml.expval(qml.PauliZ(0) @ qml.Hermitian(O, wires=[1, 2])),
            qml.expval(qml.PauliZ(0)))


@qml.qnode(_dev)
def _gadget_state_noisy(rho_AB, g_F):
    """Same noisy gadget; return the joint (ancilla, register A) state (8x8)."""
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
    qml.Hadamard(0)
    qml.ctrl(qml.SWAP, control=0)(wires=[1, 3])
    if g_F > 0:
        qml.QubitChannel(global_depol3_kraus(g_F), wires=[0, 1, 3])
    qml.ctrl(qml.SWAP, control=0)(wires=[2, 4])
    if g_F > 0:
        qml.QubitChannel(global_depol3_kraus(g_F), wires=[0, 2, 4])
    qml.Hadamard(0)
    return qml.density_matrix(wires=[0, 1, 2])


# ---------------------------------------------------------------------------
# Read-outs
# ---------------------------------------------------------------------------
def obs_pqec_noisy(eps, g_F, O=O_PHI_PLUS):
    """Purified observable <O> = <Z(x)O>/<Z(x)I> with Fredkin depolarizing g_F,
    on two copies of the input rho_eps."""
    rho = make_noisy_bell(eps)
    zO, zI = _gadget_obs_noisy(np.kron(rho, rho), g_F, O)
    return float(zO) / float(zI)


def effective_state_noisy(eps, g_F):
    """Parity-weighted effective (purified) state = (|0> block - |1> block)/trace.
    NOTE: under noise this quasi-state need not be positive-semidefinite; it is the
    object whose expectation values obs_pqec_noisy measures (like the virtual state
    in virtual distillation)."""
    rho = make_noisy_bell(eps)
    sigma = _gadget_state_noisy(np.kron(rho, rho), g_F)
    m = sigma[:4, :4] - sigma[4:, 4:]
    return m / np.trace(m)


def no_qec(eps):
    """<O_fid> with no purification = Tr(O_fid rho_eps) = 1 - 3 eps/4."""
    return 1 - 3 * eps / 4


# ===========================================================================
def main():
    print("=" * 74)
    print(" Step 3 (v1) -- 3-qubit global depolarizing g_F after each Fredkin")
    print("=" * 74)

    # -- sanity: g_F = 0 reproduces the ideal gadget ------------------------
    print("\n Sanity check (g_F = 0 must equal the ideal gadget):")
    worst = 0.0
    for eps in [0.1, 0.3, 0.5, 0.8]:
        a = obs_pqec_noisy(eps, 0.0)
        b = obs_purified_ideal(make_noisy_bell(eps))
        worst = max(worst, abs(a - b))
    print(f"   max |obs_noisy(g_F=0) - obs_ideal| over eps = {worst:.2e}  "
          f"{'PASS' if worst < 1e-13 else 'FAIL'}")

    # -- trace-preserving check for the channel -----------------------------
    K = global_depol3_kraus(0.37)
    tp = np.max(np.abs(sum(k.conj().T @ k for k in K) - np.eye(8)))
    print(f"   3-qubit channel trace-preserving:  |sum K^dag K - I| = {tp:.2e}  "
          f"{'PASS' if tp < 1e-13 else 'FAIL'}")

    # -- effect of g_F on the purified observable ---------------------------
    eps = 0.40
    print("\n" + "-" * 74)
    print(f" Purified observable vs Fredkin noise  (eps = {eps}, "
          f"no-QEC <O> = {no_qec(eps):.4f})")
    print("-" * 74)
    print(f"  {'g_F':>6} {'<O>_PQEC':>10} {'vs no-QEC':>11}   verdict")
    for g_F in [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50]:
        v = obs_pqec_noisy(eps, g_F)
        verdict = "PQEC helps" if v > no_qec(eps) else "PQEC hurts"
        print(f"  {g_F:>6.2f} {v:>10.4f} {v - no_qec(eps):>+11.4f}   {verdict}")
    print("\n  -> g_F = 0 matches the ideal one-round value; increasing g_F degrades")
    print("     the purified observable toward (and below) the no-QEC baseline.")


if __name__ == "__main__":
    main()
