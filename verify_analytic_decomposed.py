"""
Check the CNOT-only analytic decomposed-Fredkin result against the circuit.
==========================================================================

Model (CNOT-only): the textbook 8-CNOT Fredkin decomposition, **single-qubit gates
ideal**, and a two-qubit *replacement* depolarizing channel of strength e2 after
EACH CNOT.  With v = 1-e2 and the input Bell-correlation t:

    C(v)   = v^5 + 5 v^6
    D(v,t) = 1 + 3 v^4 t^2
    B = v^10/4 · D,   A = v^10/16 · [ D + t(1+t) C ],
    F_dec = A/B = (1/4)[ 1 + t(1+t) C / D ]
    F_bare = (1+3t)/4,   ideal (v=1): (1+3t)^2 / (4(1+3t^2))
    CNOT small-noise slope  K2(t) = t(1+t)(33 t^2 + 35) / (4(1+3 t^2)^2) -> 17/8 at t=1.

INPUT t.  t = (1-4p/3)^2 for local depolarizing p per Bell qubit (used here), or
equivalently t = 1-eps for a global Bell depolarizing input rho_eps (same isotropic
family).

CONVENTION.  The analytic e2 equals p2 of the circuit (both 2-qubit global depol).

Run:  python verify_analytic_decomposed.py
"""

import numpy as np
import pennylane as qml

from noisy_bell_state import global_depol_kraus, O_PHI_PLUS
from pqec_gadget import obs_purified

np.set_printoptions(precision=6, suppress=True)


# --- isotropic input via local depolarizing p per Bell qubit ---------------
_dev2 = qml.device("default.mixed", wires=2)


@qml.qnode(_dev2)
def _bell_local_p(p):
    qml.Hadamard(0)
    qml.CNOT(wires=[0, 1])
    qml.DepolarizingChannel(p, wires=0)   # Bloch x (1-4p/3) => t=(1-4p/3)^2
    qml.DepolarizingChannel(p, wires=1)
    return qml.density_matrix(wires=[0, 1])


# --- faithful decomposed gadget: 2-qubit depol e2 after each CNOT, 1q ideal --
_dev5 = qml.device("default.mixed", wires=5)


def _c2(c, t, e2):
    qml.CNOT(wires=[c, t])
    if e2 > 0:
        qml.QubitChannel(global_depol_kraus(e2), wires=[c, t])


def _tof(c1, c2, t, e2):        # Clifford+T Toffoli (6 CNOTs), single-qubit gates ideal
    qml.Hadamard(t)
    _c2(c2, t, e2); qml.adjoint(qml.T)(wires=t)
    _c2(c1, t, e2); qml.T(wires=t)
    _c2(c2, t, e2); qml.adjoint(qml.T)(wires=t)
    _c2(c1, t, e2); qml.T(wires=t); qml.T(wires=c2)
    _c2(c1, c2, e2); qml.Hadamard(t)
    qml.T(wires=c1); qml.adjoint(qml.T)(wires=c2)
    _c2(c1, c2, e2)


def _fred(q, a, b, e2):         # swaps a,b with control q; Toffoli target = b
    _c2(b, a, e2); _tof(q, a, b, e2); _c2(b, a, e2)


@qml.qnode(_dev5)
def _gadget_f(rho_AB, e2, O):
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
    qml.Hadamard(0)                       # prepare |+> (ideal)
    _fred(0, 1, 3, e2)                    # CSWAP(a; A1,B1)
    _fred(0, 2, 4, e2)                    # CSWAP(a; A2,B2)
    qml.Hadamard(0)                       # final H (ideal)
    return (qml.expval(qml.PauliZ(0) @ qml.Hermitian(O, wires=[1, 2])),
            qml.expval(qml.PauliZ(0)))


def circuit_AB(p, e2):
    rho = _bell_local_p(p)
    zO, zI = _gadget_f(np.kron(rho, rho), e2, O_PHI_PLUS)
    return float(zO), float(zI)


# --- analytic (CNOT-only, single-qubit ideal) ------------------------------
def t_of(p):
    return (1 - 4 * p / 3) ** 2


def C_v(v):
    return v**5 + 5 * v**6


def D_vt(v, t):
    return 1 + 3 * v**4 * t**2


def A_ana(p, e2):
    v, t = 1 - e2, t_of(p)
    return v**10 / 16 * (D_vt(v, t) + t * (1 + t) * C_v(v))


def B_ana(p, e2):
    v, t = 1 - e2, t_of(p)
    return v**10 / 4 * D_vt(v, t)


def Fdec_ana(p, e2):
    v, t = 1 - e2, t_of(p)
    return 0.25 * (1 + t * (1 + t) * C_v(v) / D_vt(v, t))


def Fbare(p):
    return (1 + 3 * t_of(p)) / 4


def K2_slope(t):
    return t * (1 + t) * (33 * t**2 + 35) / (4 * (1 + 3 * t**2) ** 2)  # ->17/8 at t=1


def main():
    print("=" * 78)
    print(" CNOT-only analytic decomposed-Fredkin result  vs  circuit")
    print("=" * 78)

    # (1) exact A, B, F_dec vs circuit
    errA = errB = errF = 0.0
    for p in [0.0, 0.05, 0.15, 0.30, 0.5]:
        for e2 in [0.0, 0.05, 0.15]:
            zO, zI = circuit_AB(p, e2)
            errA = max(errA, abs(zO - A_ana(p, e2)))
            errB = max(errB, abs(zI - B_ana(p, e2)))
            errF = max(errF, abs(zO / zI - Fdec_ana(p, e2)))
    print(f"\n (1) numerator  A = v^10/16 [D + t(1+t)C]   : max err = {errA:.2e}")
    print(f" (2) denominator B = v^10/4 D               : max err = {errB:.2e}")
    print(f" (3) F_dec = A/B = 1/4[1 + t(1+t)C/D]        : max err = {errF:.2e}")

    # (4) ideal limit
    errI = 0.0
    for p in [0.05, 0.2, 0.4, 0.6]:
        errI = max(errI, abs(Fdec_ana(p, 0) - obs_purified(_bell_local_p(p))))
    print(f" (4) ideal limit F_dec(p,0) = ideal PQEC     : max err = {errI:.2e}")

    # (5) CNOT small-noise slope K2 (circuit vs analytic)
    print("\n (5) CNOT small-noise slope  K2 = -dF/de2 at e2=0  (circuit vs analytic):")
    h = 1e-6
    print(f"     {'t':>5} {'K2 circuit':>12} {'K2 analytic':>12}")
    for p in [0.0, 0.15, 0.30]:
        t = t_of(p)
        k2 = -((circuit_AB(p, h)[0] / circuit_AB(p, h)[1]) - Fdec_ana(p, 0)) / h
        print(f"     {t:>5.2f} {k2:>12.4f} {K2_slope(t):>12.4f}")
    print("     -> K2(1) = 17/8.")

    worst = max(errA, errB, errF, errI)
    print("\n" + "=" * 78)
    print(f"  ANALYTIC MATCHES CIRCUIT: {'YES' if worst < 1e-11 else 'CHECK'}  "
          f"(worst {worst:.1e})")
    print("  (Single-qubit gates ideal, so orientation is irrelevant and the outer")
    print("   Hadamards carry no noise; the common v^10 cancels in F_dec.)")
    print("=" * 78)


if __name__ == "__main__":
    main()
