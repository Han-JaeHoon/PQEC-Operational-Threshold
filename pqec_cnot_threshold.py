"""
CNOT-only noise threshold for decomposed-Fredkin PQEC (single-qubit gates ideal).
=================================================================================

Model: the textbook decomposition (8 CNOTs per Fredkin, two Fredkins = 16 CNOTs);
single-qubit gates are IDEAL and a two-qubit depolarizing channel of strength
eps2 acts after EACH CNOT.  With v = 1 - eps2 and numerator

    C(v) = v^5 + 5 v^6,

the purified Bell fidelity is

    F_dec(eps2, t) = 1/4 [ 1 + t(1+t)(v^5 + 5 v^6) / (1 + 3 v^4 t^2) ],
    F_bare(t)      = (1 + 3 t) / 4,     t = 1 - eps  (global Bell-depol input).

PQEC helps iff F_dec > F_bare, i.e.

    (1 + t)(v^5 + 5 v^6) > 3 (1 + 3 v^4 t^2).

The CNOT-noise threshold eps2*(t) is the root of the equality.  Above eps2* the
noisy CNOTs inject more error than one PQEC round removes.

Run:  python pqec_cnot_threshold.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

from verify_analytic_decomposed import circuit_AB, _fred, _bell_local_p  # genuine circuit

np.set_printoptions(precision=4, suppress=True)

_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)
_dev5 = qml.device("default.mixed", wires=5)


# --- analytic (single-qubit ideal, u=1;  s = 1-eps2) -----------------------
def F_dec(eps2, t):
    """Purified Bell fidelity  F = 1/4 [1 + t(1+t)(v^5+5v^6)/(1+3 v^4 t^2)],  v=1-eps2
    (equivalently s^5(1+5s)... with s=v).  See CNOT_NOISE_ANALYSIS.md."""
    v = 1 - eps2
    return 0.25 * (1 + t * (1 + t) * (v**5 + 5 * v**6) / (1 + 3 * v**4 * t**2))


def F_bare(t):
    return (1 + 3 * t) / 4


# parity denominator, Bell-projector numerator, and the anisotropic effective-state
# Bell-diagonal correlators (rho_eff = 1/4[II + c_perp(XX-YY) + c_z ZZ]).
def Q_denom(t, eps2):
    s = 1 - eps2
    return s**10 / 4 * (1 + 3 * s**4 * t**2)


def N_num(t, eps2):
    s = 1 - eps2
    return s**10 / 16 * (1 + 3 * s**4 * t**2 + s**5 * (1 + 5 * s) * t * (1 + t))


def c_perp(t, eps2):
    s = 1 - eps2
    return 2 * s**6 * t * (1 + t) / (1 + 3 * s**4 * t**2)


def c_z(t, eps2):
    s = 1 - eps2
    return s**5 * (1 + s) * t * (1 + t) / (1 + 3 * s**4 * t**2)


def eff_correlators_circuit(t, eps2):
    """Effective-state correlators (c_perp from XX, c_z from ZZ) from the genuine
    CNOT-only circuit (single-qubit gates ideal)."""
    p = _eps_to_local_p(1 - t)                       # t = 1-eps  ->  local p
    rho = np.kron(_bell_local_p(p), _bell_local_p(p))

    @qml.qnode(_dev5)
    def run(O):
        qml.QubitDensityMatrix(rho, wires=[1, 2, 3, 4])
        qml.Hadamard(0)
        _fred(0, 1, 3, eps2)                          # Toffoli target on discarded B
        _fred(0, 2, 4, eps2)
        qml.Hadamard(0)
        return qml.expval(qml.PauliZ(0) @ O) if O is not None \
            else qml.expval(qml.PauliZ(0))

    B = float(run(None))
    cx = float(run(qml.Hermitian(np.kron(_X, _X), wires=[1, 2]))) / B
    cz = float(run(qml.Hermitian(np.kron(_Z, _Z), wires=[1, 2]))) / B
    return cx, cz


def eps2_star(t):
    """Root of F_dec = F_bare in eps2 (0 if PQEC never helps)."""
    f = lambda e: F_dec(e, t) - F_bare(t)
    if f(1e-9) <= 0:
        return 0.0
    lo, hi = 0.0, 0.999
    for _ in range(80):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if f(m) > 0 else (lo, m)
    return 0.5 * (lo + hi)


def _eps_to_local_p(eps):
    """Global-input eps -> local depolarizing p with t=(1-4p/3)^2 = 1-eps (for the
    circuit, which prepares the isotropic state via local depolarizing)."""
    return 0.75 * (1 - np.sqrt(1 - eps))


# ===========================================================================
def main():
    print("=" * 74)
    print(" CNOT-only noise threshold (single-qubit gates ideal)")
    print("=" * 74)

    # (0) analytic vs genuine circuit
    print("\n (0) analytic F_dec vs circuit (single-qubit gates ideal):")
    worst = 0.0
    for eps in [0.2, 0.4, 0.6]:
        t = 1 - eps
        p = _eps_to_local_p(eps)
        for e2 in [0.05, 0.12, 0.20]:
            zO, zI = circuit_AB(p, e2)
            worst = max(worst, abs(zO / zI - F_dec(e2, t)))
    print(f"     max|circuit - analytic| = {worst:.2e}")

    # (1) threshold table
    print("\n (1) CNOT-noise threshold eps2* vs input noise eps (t = 1-eps):")
    print(f"     {'eps':>5} {'t':>5} {'F_bare':>8} {'eps2*':>8} {'16*eps2* (budget)':>18}")
    for eps in [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60]:
        t = 1 - eps
        es = eps2_star(t)
        print(f"     {eps:>5.2f} {t:>5.2f} {F_bare(t):>8.4f} {es:>8.4f} {16*es:>18.3f}")
    print("     (16 = total CNOT count; 16*eps2* is a rough per-round CNOT budget.)")

    # (2) small-eps2 behaviour
    print("\n (2) near eps2=0:  d F_dec / d eps2 |_0  (= -K2, the CNOT slope)")
    for eps in [0.2, 0.4, 0.6]:
        t = 1 - eps
        h = 1e-6
        slope = -(F_dec(h, t) - F_dec(0, t)) / h
        print(f"     eps={eps:.1f} (t={t:.1f}):  K2 = {slope:.4f}  "
              f"(=t(1+t)(33t^2+35)/(4(1+3t^2)^2))")

    # (3) effective state is ANISOTROPIC Bell-diagonal under CNOT noise
    print("\n (3) effective state rho_eff = 1/4[II + c_perp(XX-YY) + c_z ZZ]:")
    print(f"     {'eps':>4}{'eps2':>6} | {'c_perp circ':>12}{'(ana)':>9} "
          f"{'c_z circ':>10}{'(ana)':>9} {'c_z-c_perp':>11}")
    for eps, e2 in [(0.4, 0.0), (0.4, 0.12), (0.4, 0.25)]:
        t = 1 - eps
        cx, cz = eff_correlators_circuit(t, e2)
        print(f"     {eps:>4.1f}{e2:>6.2f} | {cx:>12.5f}{c_perp(t,e2):>9.5f} "
              f"{cz:>10.5f}{c_z(t,e2):>9.5f} {cz-cx:>+11.5f}")
    print("     -> q>0 breaks isotropy (c_z > c_perp): Z-correlation is better preserved;")
    print("        q=0 gives c_z = c_perp (isotropic Bell-diagonal).")

    # ---- Figure -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))

    ax = axes[0]
    e2s = np.linspace(0, 0.25, 120)
    for eps, c in [(0.20, "C0"), (0.40, "C1"), (0.60, "C2")]:
        t = 1 - eps
        ax.plot(e2s, [F_dec(e, t) for e in e2s], "-", color=c, lw=2,
                label=f"$\\varepsilon$={eps}")
        ax.axhline(F_bare(t), color=c, ls=":", lw=1)
        es = eps2_star(t)
        ax.plot(es, F_bare(t), "o", color=c, ms=7)
    ax.text(0.252, F_bare(0.20), " no-QEC", va="center", fontsize=8)
    ax.set_xlabel(r"per-CNOT two-qubit depolarizing  $\varepsilon_2$")
    ax.set_ylabel(r"purified fidelity  $F_{dec}$")
    ax.set_title("(a) CNOT-only: $F_{dec}$ vs $\\varepsilon_2$ (dots = threshold)")
    ax.set_xlim(0, 0.27)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    es_in = np.linspace(0.02, 0.66, 60)
    ax.plot(es_in, [eps2_star(1 - e) for e in es_in], "-", color="C3", lw=2)
    ax.axhspan(0, 0.01, color="0.85", alpha=0.6)
    ax.text(0.35, 0.006, "realistic hardware CNOT error $\\sim10^{-2}$", fontsize=8,
            va="center")
    ax.set_xlabel(r"input noise  $\varepsilon$")
    ax.set_ylabel(r"CNOT-noise threshold  $\varepsilon_2^*$")
    ax.set_title("(b) $\\varepsilon_2^*$ grows with input noise")
    ax.set_ylim(0, None)

    fig.tight_layout()
    fig.savefig("pqec_cnot_threshold.png", dpi=140)
    print("\n  saved  pqec_cnot_threshold.png")


if __name__ == "__main__":
    main()
