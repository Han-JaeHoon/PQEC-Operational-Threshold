"""
Preparing the noisy Bell state  rho_eps = (1-eps)|Phi+><Phi+| + eps I/4.
========================================================================

This is the isotropic (global-depolarizing) Bell state used as the noisy input
throughout the operational-threshold study.  We build it with a *genuine* circuit
on PennyLane's mixed-state simulator -- prepare |Phi+> with H . CNOT, then apply a
global 2-qubit depolarizing channel of strength eps -- and verify that the result
matches the target density matrix.

Global depolarizing channel on d = 4:

    E(rho) = (1-eps) rho + eps I/4
           = (1-eps + eps/16) rho + (eps/16) sum_{P != I} P rho P,

using that (1/16) sum over all 16 two-qubit Paulis P of  P rho P  =  I/4 * Tr(rho).
So its Kraus operators are

    K_0   = sqrt(1 - eps + eps/16) I,
    K_i   = sqrt(eps/16) P_i     (P_i = the 15 non-identity two-qubit Paulis).

Reference spectrum of rho_eps:
    |Phi+>  eigenvalue  1 - 3 eps/4      (fidelity  F = <Phi+|rho_eps|Phi+>)
    the 3 other Bell states  eigenvalue  eps/4  each.

Run:  python noisy_bell_state.py
"""

import numpy as np
import pennylane as qml

np.set_printoptions(precision=4, suppress=True)

# ---------------------------------------------------------------------------
# Fixed operators
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
_PAULIS_1Q = (I2, X, Y, Z)
_PAULIS_2Q = [np.kron(a, b) for a in _PAULIS_1Q for b in _PAULIS_1Q]  # 16 of them

PHI_PLUS = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)   # (|00> + |11>)/sqrt2
O_PHI_PLUS = np.outer(PHI_PLUS, PHI_PLUS.conj())               # |Phi+><Phi+|

# The four Bell states, for the spectral check.
_BELL = {
    "Phi+": np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2),
    "Phi-": np.array([1, 0, 0, -1], dtype=complex) / np.sqrt(2),
    "Psi+": np.array([0, 1, 1, 0], dtype=complex) / np.sqrt(2),
    "Psi-": np.array([0, 1, -1, 0], dtype=complex) / np.sqrt(2),
}


# ---------------------------------------------------------------------------
# Analytic target and its Kraus-form global depolarizing channel
# ---------------------------------------------------------------------------
def rho_eps_analytic(eps):
    """Target state  rho_eps = (1-eps)|Phi+><Phi+| + eps I/4."""
    return (1 - eps) * O_PHI_PLUS + eps * np.eye(4, dtype=complex) / 4


def global_depol_kraus(eps):
    """Kraus operators of the global 2-qubit depolarizing channel of strength eps:
        E(rho) = (1-eps) rho + eps I/4.
    """
    K = [np.sqrt(1 - eps + eps / 16) * _PAULIS_2Q[0]]
    K += [np.sqrt(eps / 16) * P for P in _PAULIS_2Q[1:]]
    return K


# ---------------------------------------------------------------------------
# The circuit that PREPARES rho_eps
# ---------------------------------------------------------------------------
_dev = qml.device("default.mixed", wires=2)


@qml.qnode(_dev)
def make_noisy_bell(eps):
    """Prepare |Phi+> then apply global depolarizing of strength eps; return rho."""
    qml.Hadamard(wires=0)
    qml.CNOT(wires=[0, 1])
    if eps > 0:
        qml.QubitChannel(global_depol_kraus(eps), wires=[0, 1])
    return qml.density_matrix(wires=[0, 1])


# ---------------------------------------------------------------------------
# Figures of merit
# ---------------------------------------------------------------------------
def fidelity_phi_plus(rho):
    """F = <Phi+| rho |Phi+>."""
    return float(np.real(np.vdot(PHI_PLUS, rho @ PHI_PLUS)))


def purity(rho):
    return float(np.real(np.trace(rho @ rho)))


# ---------------------------------------------------------------------------
# Verification: is the noisy rho built correctly?
# ---------------------------------------------------------------------------
def verify(eps, tol=1e-12, verbose=True):
    """Check the circuit output against every property of rho_eps.
    Returns (all_ok, details)."""
    rho = make_noisy_bell(eps)
    target = rho_eps_analytic(eps)

    checks = {}
    # 1. matches the analytic target density matrix
    checks["matches_analytic"] = np.max(np.abs(rho - target)) < tol
    # 2. valid density matrix: unit trace, Hermitian, positive semidefinite
    checks["trace_one"] = abs(np.trace(rho) - 1) < tol
    checks["hermitian"] = np.max(np.abs(rho - rho.conj().T)) < tol
    evals = np.linalg.eigvalsh(rho)
    checks["psd"] = evals.min() > -tol
    # 3. correct spectrum: Phi+ -> 1-3eps/4, other Bell states -> eps/4
    bell_evals = {name: float(np.real(np.vdot(v, rho @ v))) for name, v in _BELL.items()}
    checks["eig_phi_plus"] = abs(bell_evals["Phi+"] - (1 - 3 * eps / 4)) < tol
    checks["eig_others"] = all(
        abs(bell_evals[n] - eps / 4) < tol for n in ("Phi-", "Psi+", "Psi-")
    )
    # 4. fidelity and purity match closed forms
    F = fidelity_phi_plus(rho)
    checks["fidelity"] = abs(F - (1 - 3 * eps / 4)) < tol
    P = purity(rho)
    checks["purity"] = abs(P - ((1 - 3 * eps / 4) ** 2 + 3 * (eps / 4) ** 2)) < tol

    all_ok = all(checks.values())
    details = dict(rho=rho, target=target, F=F, purity=P,
                   bell_evals=bell_evals, checks=checks,
                   max_abs_err=float(np.max(np.abs(rho - target))))

    if verbose:
        status = "PASS" if all_ok else "FAIL"
        print(f"  eps = {eps:.3f} | F = {F:.4f} (=1-3eps/4={1-3*eps/4:.4f}) | "
              f"purity = {P:.4f} | max|rho-target| = {details['max_abs_err']:.2e} "
              f"| {status}")
        if not all_ok:
            for k, ok in checks.items():
                if not ok:
                    print(f"      FAILED: {k}")
    return all_ok, details


# ===========================================================================
def main():
    print("=" * 74)
    print(" Noisy Bell state  rho_eps = (1-eps)|Phi+><Phi+| + eps I/4")
    print("=" * 74)

    # Show one state explicitly.
    eps0 = 0.40
    rho0 = make_noisy_bell(eps0)
    print(f"\n  Circuit output at eps = {eps0} (real part):\n{np.real(rho0)}")
    print(f"\n  Analytic target        (real part):\n{np.real(rho_eps_analytic(eps0))}")
    print(f"\n  Bell-basis eigenvalues at eps = {eps0}:")
    for name, v in _BELL.items():
        print(f"    {name}: {float(np.real(np.vdot(v, rho0 @ v))):.4f}")

    # Sweep and verify.
    print("\n" + "-" * 74)
    print(" Verification across eps:")
    print("-" * 74)
    all_ok = True
    for eps in [0.0, 0.1, 0.2, 1 / 3, 0.5, 2 / 3, 0.8, 1.0]:
        ok, _ = verify(eps)
        all_ok = all_ok and ok

    # Random-eps stress test.
    rng = np.random.default_rng(0)
    worst = 0.0
    for eps in rng.uniform(0, 1, size=500):
        ok, d = verify(eps, verbose=False)
        all_ok = all_ok and ok
        worst = max(worst, d["max_abs_err"])
    print("-" * 74)
    print(f"  500 random eps in [0,1]: max |rho - target| = {worst:.2e}")

    print("\n" + "=" * 74)
    print(f"  ALL CHECKS {'PASSED' if all_ok else 'FAILED'}")
    print("  (note: rho_eps is entangled iff F > 1/2, i.e. eps < 2/3.)")
    print("=" * 74)


if __name__ == "__main__":
    main()
