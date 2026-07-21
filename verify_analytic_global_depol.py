"""
Check the analytic global-depolarizing result against the circuit implementation.
================================================================================

Analytic claims (Fredkin-only 3-qubit GLOBAL depolarizing g_F), with
alpha = 1 - 4p/3 for a local-depolarizing-p Bell input, equivalently alpha^2 = 1-eps
for the global-depolarizing-eps input make_noisy_bell(eps):

    A(g_F) = <Z_a (x) Phi_A>  = (1-g_F)^2 (1+3 alpha^2)^2 / 16          (numerator)
    B(g_F) = <Z_a (x) I_A>    = (1-g_F)^2 (1+3 alpha^4) / 4            (denominator)
    F_PQEC = A/B = (1+3 alpha^2)^2 / (4 (1+3 alpha^4))     (independent of g_F<1)
    F_bare = (1+3 alpha^2)/4 = 1 - 3 eps/4 = 1 - 2p + 4p^2/3
    dF     = F_PQEC - F_bare = 3 alpha^2 (1-alpha^2)(1+3 alpha^2) / (4(1+3 alpha^4))
    N_samp ~ 1/B^2  ->  overhead vs ideal Fredkin  ~ (1-g_F)^{-4}

Run:  python verify_analytic_global_depol.py
"""

import numpy as np
import pennylane as qml

from noisy_bell_state import make_noisy_bell, O_PHI_PLUS, fidelity_phi_plus
from pqec_gadget import obs_purified
from pqec_gadget_noise import _gadget_obs_noisy, obs_pqec_noisy, no_qec


# --- analytic formulas (in alpha^2) ---------------------------------------
def A_ana(a2, g):
    return (1 - g) ** 2 * (1 + 3 * a2) ** 2 / 16


def B_ana(a2, g):
    return (1 - g) ** 2 * (1 + 3 * a2 ** 2) / 4


def Fpqec_ana(a2):
    return (1 + 3 * a2) ** 2 / (4 * (1 + 3 * a2 ** 2))


def Fbare_ana(a2):
    return (1 + 3 * a2) / 4


def dF_ana(a2):
    return 3 * a2 * (1 - a2) * (1 + 3 * a2) / (4 * (1 + 3 * a2 ** 2))


# --- local-depolarizing Bell input, exactly as the derivation sets it up ----
_dev2 = qml.device("default.mixed", wires=2)


@qml.qnode(_dev2)
def _bell_local_depol(p):
    qml.Hadamard(0)
    qml.CNOT(wires=[0, 1])
    qml.DepolarizingChannel(p, wires=0)
    qml.DepolarizingChannel(p, wires=1)
    return qml.density_matrix(wires=[0, 1])


def main():
    print("=" * 78)
    print(" Analytic global-depolarizing result  vs  circuit")
    print("=" * 78)

    eps_grid = np.linspace(0.0, 0.99, 34)
    g_grid = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 0.95]

    # (1) correlators A, B and their (1-g)^2 scaling; (2) ratio g-independence
    errA = errB = errF = 0.0
    for eps in eps_grid:
        a2 = 1 - eps
        rho = make_noisy_bell(eps)
        rr = np.kron(rho, rho)
        for g in g_grid:
            zO, zI = _gadget_obs_noisy(rr, g, O_PHI_PLUS)
            errA = max(errA, abs(float(zO) - A_ana(a2, g)))
            errB = max(errB, abs(float(zI) - B_ana(a2, g)))
            errF = max(errF, abs(float(zO) / float(zI) - Fpqec_ana(a2)))
    print(f"\n (1) numerator  A = (1-g_F)^2 (1+3a^2)^2/16 : max err = {errA:.2e}")
    print(f" (2) denominator B = (1-g_F)^2 (1+3a^4)/4   : max err = {errB:.2e}")
    print(f" (3) ratio F_PQEC = (1+3a^2)^2/(4(1+3a^4)), g_F-independent : "
          f"max err = {errF:.2e}")

    # (4) F_bare and dF
    errFb = errdF = 0.0
    for eps in eps_grid:
        a2 = 1 - eps
        rho = make_noisy_bell(eps)
        Fb = fidelity_phi_plus(rho)
        errFb = max(errFb, abs(Fb - Fbare_ana(a2)), abs(Fb - no_qec(eps)))
        Fp = obs_purified(rho)                          # ideal one-round PQEC
        errdF = max(errdF, abs((Fp - Fb) - dF_ana(a2)))
    print(f" (4) F_bare = 1-3eps/4 = (1+3a^2)/4         : max err = {errFb:.2e}")
    print(f" (5) dF = F_PQEC - F_bare (boxed formula)   : max err = {errdF:.2e}")

    # (6) direct check with the derivation's own local-depolarizing-p input
    print("\n (6) local-depolarizing-p input (the derivation's own parametrization):")
    print(f"     {'p':>5} {'eps=1-a^2':>10} {'F_bare(p)':>10} {'1-2p+4p^2/3':>12} "
          f"{'F_PQEC circ':>12} {'analytic':>10}")
    errP = 0.0
    for p in [0.0, 0.1, 0.2, 0.3, 0.5, 0.75]:
        a2 = (1 - 4 * p / 3) ** 2
        rho = _bell_local_depol(p)
        # purified fidelity via the noisy gadget (pick a nonzero g_F to prove independence)
        zO, zI = _gadget_obs_noisy(np.kron(rho, rho), 0.37, O_PHI_PLUS)
        Fp = float(zO) / float(zI)
        Fb = fidelity_phi_plus(rho)
        errP = max(errP, abs(Fp - Fpqec_ana(a2)), abs(Fb - (1 - 2 * p + 4 * p ** 2 / 3)))
        print(f"     {p:>5.2f} {1-a2:>10.4f} {Fb:>10.4f} {1-2*p+4*p**2/3:>12.4f} "
              f"{Fp:>12.4f} {Fpqec_ana(a2):>10.4f}")
    print(f"     max err (local-p form, g_F=0.37) = {errP:.2e}")

    # (7) sampling overhead ~ (1-g_F)^-4
    print("\n (7) sampling overhead  N_samp(g_F)/N_samp(0) ~ 1/(1-g_F)^4  "
          "(from B ∝ (1-g_F)^2):")
    eps = 0.40
    rr = np.kron(make_noisy_bell(eps), make_noisy_bell(eps))
    _, B0 = _gadget_obs_noisy(rr, 0.0, O_PHI_PLUS)
    for g in [0.1, 0.3, 0.5, 0.9]:
        _, Bg = _gadget_obs_noisy(rr, g, O_PHI_PLUS)
        print(f"     g_F={g:.2f}:  (B0/Bg)^2 = {(float(B0)/float(Bg))**2:>8.2f}   "
              f"1/(1-g)^4 = {1/(1-g)**4:>8.2f}")

    # ratios of expectation values are division-amplified near eps->1, so use a
    # tolerance appropriate for a ratio (1e-12) rather than the 1e-15 of a raw matrix
    worst = max(errA, errB, errF, errFb, errdF, errP)
    ok = worst < 1e-12
    print("\n" + "=" * 78)
    print(f"  ANALYTIC MATCHES CIRCUIT: {'YES' if ok else 'MISMATCH'}  "
          f"(worst deviation {worst:.1e})")
    print("=" * 78)


if __name__ == "__main__":
    main()
