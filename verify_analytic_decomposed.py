"""
Check the analytic decomposed-Fredkin result against the circuit.
=================================================================

Analytic model (gate-local depolarizing on the 8-CNOT textbook Fredkin):
  * 1-qubit gates:  E1(e1)(rho) = (1-e1) rho + e1 I2/2   (u = 1-e1)
  * CNOTs:          E2(e2)(rho) = (1-e2) rho + e2 I4/4   (v = 1-e2)
  * input:          local depolarizing p per Bell qubit, t = (1-4p/3)^2

Analytic result (with C = 2u^4v^6+u^6v^5+3u^6v^6,  D = 1+(1+2u^4)v^4 t^2):
  B = (u^9 v^10 / 4)  D
  A = (u^9 v^10 /16) [ D + t(1+t) C ]
  F_dec = A/B = (1/4)[ 1 + t(1+t) C / D ]
  F_bare = (1+3t)/4,   ideal (u=v=1): (1+3t)^2/(4(1+3t^2))
  slopes  K1 = 5/2 (single-qubit), K2 = 17/8 (CNOT) at t->1.

ORIENTATION.  The result above corresponds to the Toffoli TARGET on the RETAINED
register ("retain"): the kept register absorbs the H/T single-qubit gates.  The
other orientation ("discard", target on the discarded register) shares the same
denominator D but has a different numerator,
    C_D(u,v) = u^6 v^5 + (1+u^2+2u^6+u^8) v^6,
giving the milder single-qubit slope K1 = 2 (vs 5/2).  The t-dependent slopes are
    K1_retain(t)  = 4 t(1+t)(3t^2+2)/(1+3t^2)^2      -> 5/2  at t=1,
    K1_discard(t) =   t(1+t)(9t^2+7)/(1+3t^2)^2      -> 2    at t=1,
    K2(t)         =   t(1+t)(33t^2+35)/(4(1+3t^2)^2) -> 17/8 at t=1  (both).
We verify both orientations against the circuit.

OUTER HADAMARDS.  This verifier uses IDEAL outer gadget Hadamards, so its raw
(A,B) carry the common factor u^9 v^10.  The main circuit (pqec_decomposed_noise.py)
also depolarizes the two outer ancilla H's, which multiplies both A and B by a
further u^2 (common factor u^11 v^10) -- this cancels in F_dec, so thresholds are
unchanged.

INPUT t.  t = (1-4p/3)^2 for local depolarizing p per Bell qubit (used here), or
equivalently t = 1-eps for a global Bell depolarizing input rho_eps -- the same
isotropic family.

CONVENTION.  The analytic e2 equals p2 of pqec_decomposed_noise.py (both 2-qubit
global depol); the analytic e1 relates to that file's PennyLane 1-qubit p1 by
e1 = 4 p1 / 3 (asserted below).

Run:  python verify_analytic_decomposed.py
"""

import numpy as np
import pennylane as qml

from noisy_bell_state import global_depol_kraus, O_PHI_PLUS
from pqec_gadget import obs_purified

np.set_printoptions(precision=6, suppress=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def depol1_kraus(e1):
    """Exact E1: (1-e1) rho + e1 I2/2."""
    return [np.sqrt(1 - e1 + e1 / 4) * I2,
            np.sqrt(e1 / 4) * X, np.sqrt(e1 / 4) * Y, np.sqrt(e1 / 4) * Z]


_dev2 = qml.device("default.mixed", wires=2)


@qml.qnode(_dev2)
def _bell_local_p(p):
    qml.Hadamard(0)
    qml.CNOT(wires=[0, 1])
    qml.DepolarizingChannel(p, wires=0)   # Bloch x (1-4p/3) => t=(1-4p/3)^2
    qml.DepolarizingChannel(p, wires=1)
    return qml.density_matrix(wires=[0, 1])


# --- faithful decomposed gadget (E1/E2 channels; no noise on the gadget H's) --
_dev5 = qml.device("default.mixed", wires=5)


def _c1(gate, w, e1, adj=False):
    (qml.adjoint(gate)(wires=w) if adj else gate(wires=w))
    if e1 > 0:
        qml.QubitChannel(depol1_kraus(e1), wires=w)


def _c2(c, t, e2):
    qml.CNOT(wires=[c, t])
    if e2 > 0:
        qml.QubitChannel(global_depol_kraus(e2), wires=[c, t])


def _tof(c1, c2, t, e1, e2):
    _c1(qml.Hadamard, t, e1)
    _c2(c2, t, e2); _c1(qml.T, t, e1, adj=True)
    _c2(c1, t, e2); _c1(qml.T, t, e1)
    _c2(c2, t, e2); _c1(qml.T, t, e1, adj=True)
    _c2(c1, t, e2); _c1(qml.T, t, e1); _c1(qml.T, c2, e1)
    _c2(c1, c2, e2); _c1(qml.Hadamard, t, e1)
    _c1(qml.T, c1, e1); _c1(qml.T, c2, e1, adj=True)
    _c2(c1, c2, e2)


def _fred(q, a, b, e1, e2):     # swaps a,b; Toffoli target = b
    _c2(b, a, e2); _tof(q, a, b, e1, e2); _c2(b, a, e2)


@qml.qnode(_dev5)
def _gadget_f(rho_AB, e1, e2, O, orient):
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
    qml.Hadamard(0)                       # prepare |+> (no noise)
    if orient == "retain":                # Toffoli target on retained register A
        _fred(0, 3, 1, e1, e2); _fred(0, 4, 2, e1, e2)
    else:                                 # "discard": target on discarded register B
        _fred(0, 1, 3, e1, e2); _fred(0, 2, 4, e1, e2)
    qml.Hadamard(0)                       # final H (no noise)
    return (qml.expval(qml.PauliZ(0) @ qml.Hermitian(O, wires=[1, 2])),
            qml.expval(qml.PauliZ(0)))


def circuit_AB(p, e1, e2, orient="retain"):
    rho = _bell_local_p(p)
    zO, zI = _gadget_f(np.kron(rho, rho), e1, e2, O_PHI_PLUS, orient)
    return float(zO), float(zI)


# --- analytic (matches the "retain" orientation) ---------------------------
def t_of(p):
    return (1 - 4 * p / 3) ** 2


def C_uv(u, v):
    return 2 * u**4 * v**6 + u**6 * v**5 + 3 * u**6 * v**6


def D_uvt(u, v, t):
    return 1 + (1 + 2 * u**4) * v**4 * t**2


def A_ana(p, e1, e2):
    u, v, t = 1 - e1, 1 - e2, t_of(p)
    return u**9 * v**10 / 16 * (D_uvt(u, v, t) + t * (1 + t) * C_uv(u, v))


def B_ana(p, e1, e2):
    u, v, t = 1 - e1, 1 - e2, t_of(p)
    return u**9 * v**10 / 4 * D_uvt(u, v, t)


def Fdec_ana(p, e1, e2):
    u, v, t = 1 - e1, 1 - e2, t_of(p)
    return 0.25 * (1 + t * (1 + t) * C_uv(u, v) / D_uvt(u, v, t))


def Fbare(p):
    return (1 + 3 * t_of(p)) / 4


# --- discard orientation: closed form C_D (denominator D is the same) --------
def C_D_uv(u, v):
    """v^6 coefficient P_D(u) = 1+u^2+2u^6+u^8 (P_D(1)=5, P_D'(1)=22); v^5 term u^6."""
    return u**6 * v**5 + (1 + u**2 + 2 * u**6 + u**8) * v**6


def Fdec_D_ana(p, e1, e2):
    u, v, t = 1 - e1, 1 - e2, t_of(p)
    return 0.25 * (1 + t * (1 + t) * C_D_uv(u, v) / D_uvt(u, v, t))


# --- small-noise slopes  F_dec ~= F_ideal - K1 e1 - K2 e2 -------------------
def K1_retain(t):
    return 4 * t * (1 + t) * (3 * t**2 + 2) / (1 + 3 * t**2) ** 2      # ->5/2 at t=1


def K1_discard(t):
    return t * (1 + t) * (9 * t**2 + 7) / (1 + 3 * t**2) ** 2          # ->2 at t=1


def K2_slope(t):
    return t * (1 + t) * (33 * t**2 + 35) / (4 * (1 + 3 * t**2) ** 2)  # ->17/8 at t=1


def _thr(p, which, orient):
    """Circuit threshold: solve F_dec(circuit) = F_bare for e1 (e2=0) or e2 (e1=0)."""
    fb = Fbare(p)

    def f(e):
        zO, zI = circuit_AB(p, e, 0.0, orient) if which == "e1" \
            else circuit_AB(p, 0.0, e, orient)
        return zO / zI - fb
    lo, hi = 0.0, 0.9
    if f(1e-9) <= 0:
        return 0.0
    for _ in range(40):
        m = 0.5 * (lo + hi)
        if f(m) > 0:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def main():
    print("=" * 78)
    print(" Analytic decomposed-Fredkin result  vs  circuit")
    print("=" * 78)

    # (0) convention: PennyLane DepolarizingChannel(p1) == E1(4 p1/3)
    rho = np.array([[0.6, 0.3], [0.3, 0.4]], complex)
    p1 = 0.03
    lhs = sum(k @ rho @ k.conj().T for k in depol1_kraus(4 * p1 / 3))
    rhs = (1 - p1) * rho + p1 / 3 * (X @ rho @ X + Y @ rho @ Y + Z @ rho @ Z)
    print(f"\n (0) E1(4 p1/3) == PennyLane DepolarizingChannel(p1):  "
          f"diff = {np.max(np.abs(lhs - rhs)):.2e}")

    # (1) exact A, B, F_dec vs the "retain" orientation (the analytic circuit)
    errA = errB = errF = 0.0
    for p in [0.0, 0.05, 0.15, 0.30, 0.5]:
        for e1 in [0.0, 0.03, 0.10]:
            for e2 in [0.0, 0.05, 0.15]:
                zO, zI = circuit_AB(p, e1, e2, "retain")
                errA = max(errA, abs(zO - A_ana(p, e1, e2)))
                errB = max(errB, abs(zI - B_ana(p, e1, e2)))
                errF = max(errF, abs(zO / zI - Fdec_ana(p, e1, e2)))
    print(f"\n (1) numerator  A  (Sec.8, retain orientation) : max err = {errA:.2e}")
    print(f" (2) denominator B (Sec.7)                     : max err = {errB:.2e}")
    print(f" (3) F_dec = A/B  (Sec.9)                      : max err = {errF:.2e}")

    # (4) ideal limit
    errI = 0.0
    for p in [0.05, 0.2, 0.4, 0.6]:
        errI = max(errI, abs(Fdec_ana(p, 0, 0) - obs_purified(_bell_local_p(p))))
    print(f" (4) ideal limit F_dec(p,0,0) = ideal PQEC     : max err = {errI:.2e}")

    # (5) denominator B is orientation-INDEPENDENT
    errBo = 0.0
    for p in [0.1, 0.3]:
        for e1 in [0.05, 0.15]:
            for e2 in [0.05, 0.15]:
                _, bR = circuit_AB(p, e1, e2, "retain")
                _, bD = circuit_AB(p, e1, e2, "discard")
                errBo = max(errBo, abs(bR - bD))
    print(f" (5) denominator B is orientation-independent  : max diff = {errBo:.2e}")

    # (5b) discard-orientation closed form C_D(u,v) = u^6 v^5 + (1+u^2+2u^6+u^8) v^6
    errFD = 0.0
    for p in [0.0, 0.05, 0.15, 0.30, 0.5]:
        for e1 in [0.0, 0.05, 0.15, 0.25]:
            for e2 in [0.0, 0.05, 0.15]:
                zO, zI = circuit_AB(p, e1, e2, "discard")
                errFD = max(errFD, abs(zO / zI - Fdec_D_ana(p, e1, e2)))
    print(f"(5b) discard closed form F_D (C_D, Sec.4) vs circuit : max err = {errFD:.2e}")

    # (6) slopes: K1 orientation-DEPENDENT, K2 orientation-independent.
    #     Verify the t-dependent formulas K1_retain/K1_discard/K2_slope on the circuit.
    print("\n (6) small-noise slopes  K = -dF/de at e=0  (circuit vs analytic K(t)):")
    h = 1e-6
    print(f"     {'t':>5} {'K1 retain':>20} {'K1 discard':>20} {'K2':>18}")
    for p in [0.0, 0.15, 0.30]:
        t = t_of(p)
        kr = -((circuit_AB(p, h, 0, "retain")[0] / circuit_AB(p, h, 0, "retain")[1])
               - Fdec_ana(p, 0, 0)) / h
        kd = -((circuit_AB(p, h, 0, "discard")[0] / circuit_AB(p, h, 0, "discard")[1])
               - Fdec_D_ana(p, 0, 0)) / h
        k2 = -((circuit_AB(p, 0, h, "retain")[0] / circuit_AB(p, 0, h, "retain")[1])
               - Fdec_ana(p, 0, 0)) / h
        print(f"     {t:>5.2f} {kr:>9.4f}/{K1_retain(t):<9.4f} "
              f"{kd:>9.4f}/{K1_discard(t):<9.4f} {k2:>8.4f}/{K2_slope(t):<8.4f}")
    print("     (circuit / analytic)  ->  K1_retain(1)=5/2, K1_discard(1)=2, K2(1)=17/8")

    # (7) orientation-dependent threshold table (circuit)
    print("\n (7) single-noise thresholds by orientation (circuit):")
    print(f"     {'p':>5} {'e1_th retain':>13}{'e1_th discard':>14}   "
          f"{'e2_th retain':>13}{'e2_th discard':>14}")
    for p in [0.05, 0.10, 0.20, 0.30, 0.40]:
        print(f"     {p:>5.2f} {_thr(p,'e1','retain'):>13.5f}"
              f"{_thr(p,'e1','discard'):>14.5f}   "
              f"{_thr(p,'e2','retain'):>13.5f}{_thr(p,'e2','discard'):>14.5f}")

    worst = max(errA, errB, errF, errI, errFD)
    print("\n" + "=" * 78)
    print(f"  ANALYTIC MATCHES CIRCUIT (both orientations): "
          f"{'YES' if worst < 1e-11 else 'CHECK'}  (worst {worst:.1e})")
    print("  B and K2 orientation-independent; K1, C(u,v), e1-threshold are not.")
    print("  (Verifier uses ideal outer gadget Hadamards; the main circuit's outer-H")
    print("   noise multiplies A and B by a common u^2 that cancels in F_dec.)")
    print("=" * 78)


if __name__ == "__main__":
    main()
