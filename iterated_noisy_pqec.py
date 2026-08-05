"""
Iterated noisy PQEC: the effective nonlinear round map and its fixed points.
============================================================================

One-round PQEC takes two identical copies rho (x) rho plus a SWAP-test ancilla, runs the
(noisy) gadget, and reconstructs the parity-weighted operator on the retained register:

    sigma_in  = |0><0|_a (x) rho_A (x) rho_B ,        rho_A = rho_B = rho
    sigma_out = E_q(sigma_in)                          (full circuit + CNOT noise)
    tau_A     = Tr_{a,B} [ (Z_a (x) I_A (x) I_B) sigma_out ]
    P_q(rho)  = tau_A / Tr(tau_A)                      <-- NONLINEAR (normalisation)

This module iterates rho_{n+1} = P_q(rho_n) for the three circuits already verified in
this repository and analyses the limiting behaviour.

SOURCE OF TRUTH.  The circuits are NOT re-derived here.  Their gate sequences are
captured directly from the verified implementations:

  step3 : verify_analytic_decomposed._fred / _tof / _c2   (textbook 16-CNOT)
  step4 : pqec_resynth_noise.GATES                        (resynthesised 14-CNOT)
  step5 : pqc_ring_prune.ansatz_masked + pqc_ring_pruned_params.npy (learned 14-CNOT)

QUBIT ORDERING (identical in all three, as in the repo):
  wire 0 = ancilla a,  wires 1,2 = retained register A,  wires 3,4 = discarded B.
  The input density matrix kron(rho, rho) is placed on wires [1,2,3,4], and the read-out
  uses Z_a AFTER the final Hadamard -- i.e. exactly the tau_A above (no X_a convention).

NOISE.  After every CNOT on wires (i,j) the two-qubit *replacement* depolarizing channel

    D_q^(ij)(rho) = (1-q) rho + q [ I_ij/4 (x) Tr_ij(rho) ]

is applied.  It is implemented here straight from that definition (partial trace +
identity replacement), and checked against the repository's Kraus implementation
(noisy_bell_state.global_depol_kraus) in `selftest()` / the test-suite.

PHYSICAL CAVEAT.  tau_A is a parity-WEIGHTED (post-selection-like) operator, not the
unconditional physical output state of the round.  Iterating it assumes the reconstructed
effective state can be re-prepared as the two identical inputs of the next round; it is
an *effective nonlinear map*, not a shot-by-shot experimental protocol.  Every iterate is
therefore checked for hermiticity, trace, positivity and denominator size, and violations
are reported rather than silently normalised away.

Run:  python iterated_noisy_pqec.py      # self-test + one-round validation
"""

import functools
import json
import os

import numpy as np
import pennylane as qml

ROOT = os.path.dirname(os.path.abspath(__file__))

N_WIRES = 5
DIM = 2 ** N_WIRES
ANCILLA = 0
KEEP = (1, 2)          # retained register A
DISCARD = (0, 3, 4)    # ancilla + discarded register B

# --- single-qubit / two-qubit constants ------------------------------------
I2 = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)
_PAULI1 = {"I": I2, "X": _X, "Y": _Y, "Z": _Z}

PHI_PLUS = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
PHI = np.outer(PHI_PLUS, PHI_PLUS.conj())                     # |Phi+><Phi+|

# Bell basis (Phi+, Phi-, Psi+, Psi-)
_BELL = {
    "Phi+": np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2),
    "Phi-": np.array([1, 0, 0, -1], dtype=complex) / np.sqrt(2),
    "Psi+": np.array([0, 1, 1, 0], dtype=complex) / np.sqrt(2),
    "Psi-": np.array([0, 1, -1, 0], dtype=complex) / np.sqrt(2),
}
BELL_LABELS = ["Phi+", "Phi-", "Psi+", "Psi-"]
BELL_U = np.stack([_BELL[k] for k in BELL_LABELS])            # rows = Bell vectors

PAULI_LABELS = [a + b for a in "IXYZ" for b in "IXYZ"]
PAULI_2Q = {a + b: np.kron(_PAULI1[a], _PAULI1[b]) for a in "IXYZ" for b in "IXYZ"}


# ===========================================================================
# input state
# ===========================================================================
def rho_isotropic(t):
    """rho_t = 1/4 [ II + t (XX - YY + ZZ) ]   (= (1-eps) Phi + eps I/4 with t = 1-eps)."""
    return 0.25 * (PAULI_2Q["II"] + t * (PAULI_2Q["XX"] - PAULI_2Q["YY"] + PAULI_2Q["ZZ"]))


def F_bare(t):
    """Bell fidelity of the bare input rho_t."""
    return (1 + 3 * t) / 4


# ===========================================================================
# dense 32x32 utilities
# ===========================================================================
def _embed(op, wires):
    """Embed a 1- or 2-qubit matrix `op` acting on `wires` into the full 32x32 space."""
    wires = tuple(int(w) for w in wires)
    k = len(wires)
    full = np.zeros((DIM, DIM), dtype=complex)
    # build via tensor reshaping: op (x) I on the complementary wires, then permute
    rest = [w for w in range(N_WIRES) if w not in wires]
    M = np.kron(op, np.eye(2 ** len(rest), dtype=complex))     # order: wires + rest
    T = M.reshape([2] * N_WIRES + [2] * N_WIRES)
    perm = list(wires) + rest
    inv = np.argsort(perm)
    T = np.transpose(T, list(inv) + [N_WIRES + p for p in inv])
    full = T.reshape(DIM, DIM)
    return full


def _ptrace_keep(rho, keep):
    """Partial trace of a 32x32 operator, keeping `keep` (ordered) wires."""
    keep = list(keep)
    rest = [w for w in range(N_WIRES) if w not in keep]
    perm = keep + rest
    T = rho.reshape([2] * N_WIRES + [2] * N_WIRES)
    T = np.transpose(T, perm + [N_WIRES + p for p in perm])
    dk, dt = 2 ** len(keep), 2 ** len(rest)
    T = T.reshape(dk, dt, dk, dt)
    return np.einsum("abcb->ac", T)


def replacement_depol(rho, wires, q):
    """D_q^(ij)(rho) = (1-q) rho + q [ I_ij/4 (x) Tr_ij(rho) ]  -- straight from the
    definition (partial trace + identity replacement), on the full 32x32 operator."""
    if q == 0.0:
        return rho
    wires = list(wires)
    rest = [w for w in range(N_WIRES) if w not in wires]
    perm = wires + rest
    T = rho.reshape([2] * N_WIRES + [2] * N_WIRES)
    T = np.transpose(T, perm + [N_WIRES + p for p in perm])
    T = T.reshape(4, 8, 4, 8)
    traced = np.einsum("abad->bd", T)                       # Tr_ij(rho), 8x8
    repl = np.eye(4, dtype=complex)[:, None, :, None] * traced[None, :, None, :] / 4.0
    out = (1 - q) * T + q * repl
    out = out.reshape([2] * (2 * N_WIRES))
    inv = np.argsort(perm)
    out = np.transpose(out, list(inv) + [N_WIRES + p for p in inv])
    return out.reshape(DIM, DIM)


# ===========================================================================
# circuit definitions -- captured from the VERIFIED repository implementations
# ===========================================================================
def _tape_ops(qfunc, *args):
    with qml.queuing.AnnotatedQueue() as qu:
        qfunc(*args)
    return qml.tape.QuantumScript.from_queue(qu).operations


def _qfunc_step3(q):
    """Textbook 16-CNOT decomposition -- verify_analytic_decomposed._fred."""
    import verify_analytic_decomposed as v3
    qml.Hadamard(0)
    v3._fred(0, 1, 3, q)            # CSWAP(a; A1,B1)
    v3._fred(0, 2, 4, q)            # CSWAP(a; A2,B2)
    qml.Hadamard(0)


def _qfunc_step4(q):
    """Resynthesised 14-CNOT circuit -- pqec_resynth_noise.GATES."""
    import pqec_resynth_noise as r4
    r4._apply(q)


def _qfunc_step5(q):
    """Learned & pruned 14-CNOT circuit -- pqc_ring_threshold._apply."""
    import pqc_ring_threshold as r5
    r5._apply(q)


CIRCUITS = {
    "step3": ("textbook 16-CNOT", _qfunc_step3),
    "step4": ("resynthesized 14-CNOT", _qfunc_step4),
    "step5": ("learned 14-CNOT", _qfunc_step5),
}


@functools.lru_cache(maxsize=None)
def build_program(circuit, q):
    """Return (ops, n_cnot) where ops is a list of ('U', M32) / ('D', (i,j)).

    Consecutive unitaries are merged into a single 32x32 matrix (exact, just matrix
    products) so one round is a short sequence of conjugations and channels.
    """
    if circuit not in CIRCUITS:
        raise KeyError(f"unknown circuit {circuit!r}; choose from {list(CIRCUITS)}")
    raw = _tape_ops(CIRCUITS[circuit][1], q)
    ops, acc, n_cnot = [], np.eye(DIM, dtype=complex), 0
    for op in raw:
        if op.name == "QubitChannel":
            # the repository applies global_depol_kraus(q) here; we implement the
            # channel from its definition (checked to agree in selftest()).
            ops.append(("U", acc))
            acc = np.eye(DIM, dtype=complex)
            ops.append(("D", tuple(int(w) for w in op.wires)))
        else:
            if op.name == "CNOT":
                n_cnot += 1
            acc = _embed(qml.matrix(op), op.wires) @ acc
    ops.append(("U", acc))
    return tuple(ops), n_cnot


def n_cnots(circuit):
    return build_program(circuit, 0.0)[1]


# ===========================================================================
# the effective one-round map
# ===========================================================================
_Z_ANC = _embed(_Z, (ANCILLA,))
_KET0 = np.array([[1, 0], [0, 0]], dtype=complex)


def one_round_tau(rho, q, circuit):
    """Parity-weighted effective operator tau_A = Tr_{a,B}[(Z_a (x) I) sigma_out] (4x4)."""
    ops, _ = build_program(circuit, float(q))
    sigma = np.kron(_KET0, np.kron(rho, rho))          # a | A(1,2) | B(3,4)
    for kind, payload in ops:
        if kind == "U":
            sigma = payload @ sigma @ payload.conj().T
        else:
            sigma = replacement_depol(sigma, payload, float(q))
    return _ptrace_keep(_Z_ANC @ sigma, KEEP)


def one_round_effective_map(rho, q, circuit):
    """rho -> (rho_next, info).  rho_next = tau_A / Tr(tau_A); `info` carries the
    diagnostics that decide whether the iterate is still a physical state.

    NUMERICAL NOTE (important).  tau_A is Hermitian in exact arithmetic (Z_a is supported
    entirely on the traced-out ancilla, so partial-trace cyclicity applies).  The map is
    however QUADRATIC in rho (it consumes rho (x) rho), so a floating-point anti-Hermitian
    residue i*A enters through both slots and is amplified by exactly x2 per round -- an
    unstable direction of the map *lifted* to general matrices, lying outside the physical
    Hermitian manifold.  Left alone it reaches O(1) after ~50 rounds and destroys the
    iterate.  We therefore project each iterate back onto the Hermitian manifold and
    RECORD the discarded anti-Hermitian norm (`herm_err`) so the artifact stays visible;
    it must remain at machine-precision level for the run to be trustworthy.
    """
    tau = one_round_tau(rho, q, circuit)
    herm = float(np.linalg.norm(tau - tau.conj().T))
    tau = 0.5 * (tau + tau.conj().T)                    # exact-arithmetic identity
    Q = float(np.real(np.trace(tau)))
    if Q == 0.0 or not np.isfinite(Q):
        info = dict(Q=Q, herm_err=herm, min_eig=np.nan, trace_err=np.nan, ok=False)
        return None, info
    rho_next = tau / Q
    eig = np.linalg.eigvalsh(rho_next)
    info = dict(Q=Q, herm_err=herm, min_eig=float(eig.min()),
                trace_err=float(abs(np.trace(rho_next) - 1.0)),
                eigs=eig, ok=True)
    return rho_next, info


# ===========================================================================
# state diagnostics
# ===========================================================================
def bell_populations(rho):
    """(p_Phi+, p_Phi-, p_Psi+, p_Psi-) and the off-diagonal norm in the Bell basis."""
    R = BELL_U.conj() @ rho @ BELL_U.T                  # rho in the Bell basis
    pops = np.real(np.diag(R))
    off = R - np.diag(np.diag(R))
    return pops, float(np.linalg.norm(off))


def pauli_coeffs(rho):
    """{label: Tr(P rho)} for all 16 two-qubit Paulis (real parts)."""
    return {lab: float(np.real(np.trace(P @ rho))) for lab, P in PAULI_2Q.items()}


def werner_deviation(rho):
    """||rho - rho_W(F)||_F with rho_W the isotropic state of the same Bell fidelity."""
    F = float(np.real(np.trace(PHI @ rho)))
    t = (4 * F - 1) / 3.0
    return float(np.linalg.norm(rho - rho_isotropic(t))), F, t


def state_record(rho, Q=None):
    """All per-iteration quantities requested for the analysis."""
    rho_h = 0.5 * (rho + rho.conj().T)
    eig = np.linalg.eigvalsh(rho_h)
    pops, off = bell_populations(rho)
    pc = pauli_coeffs(rho)
    dev, F, t_eff = werner_deviation(rho)
    rec = dict(
        F=F, purity=float(np.real(np.trace(rho @ rho))),
        eig_min=float(eig.min()), eig_max=float(eig.max()),
        eigs=eig.tolist(),
        herm_err=float(np.linalg.norm(rho - rho.conj().T)),
        trace_err=float(abs(np.trace(rho) - 1.0)),
        p_Phi_plus=float(pops[0]), p_Phi_minus=float(pops[1]),
        p_Psi_plus=float(pops[2]), p_Psi_minus=float(pops[3]),
        bell_offdiag=off,
        delta_iso=dev, t_eff=t_eff,
        XX=pc["XX"], YY=pc["YY"], ZZ=pc["ZZ"],
        iso_XXplusYY=abs(pc["XX"] + pc["YY"]), iso_XXminusZZ=abs(pc["XX"] - pc["ZZ"]),
    )
    rec.update({f"P_{k}": v for k, v in pc.items()})
    if Q is not None:
        rec["Q"] = Q
    return rec


# ===========================================================================
# iteration
# ===========================================================================
def iterate_effective_map(rho0, q, circuit, tol=1e-12, max_iter=1000, record_every=1):
    """Iterate rho_{n+1} = P_q(rho_n).

    Returns dict with the trajectory records, convergence flags, the final state, the
    fixed-point residual, and the 2-cycle distance ||rho_{n+2}-rho_n||_F.
    """
    rho = np.array(rho0, dtype=complex)
    hist, prev, prev2 = [], None, None
    status, n_used = "max_iter", 0
    diff = np.nan
    diff2 = np.nan
    max_herm, min_eig_seen, min_Q = 0.0, np.inf, np.inf

    rec0 = state_record(rho)
    rec0.update(n=0, diff=np.nan, diff2=np.nan, Q=np.nan)
    hist.append(rec0)

    for n in range(1, max_iter + 1):
        nxt, info = one_round_effective_map(rho, q, circuit)
        n_used = n
        max_herm = max(max_herm, info["herm_err"])
        min_Q = min(min_Q, info["Q"])
        if nxt is None or not np.all(np.isfinite(nxt)):
            status = "denominator_vanished"
            break
        min_eig_seen = min(min_eig_seen, info["min_eig"])
        diff = float(np.linalg.norm(nxt - rho))
        diff2 = float(np.linalg.norm(nxt - prev)) if prev is not None else np.nan
        prev2, prev = prev, rho
        rho = nxt
        if n % record_every == 0 or n <= 20:
            rec = state_record(rho, Q=info["Q"])
            rec.update(n=n, diff=diff, diff2=diff2)
            hist.append(rec)
        if diff < tol:
            status = "converged"
            break

    # fixed-point residual, recomputed on the final state
    res_rho, res_info = one_round_effective_map(rho, q, circuit)
    residual = float(np.linalg.norm(res_rho - rho)) if res_rho is not None else np.nan
    # 2-cycle check: ||P^2(rho) - rho||
    if res_rho is not None:
        r2, _ = one_round_effective_map(res_rho, q, circuit)
        cyc2 = float(np.linalg.norm(r2 - rho)) if r2 is not None else np.nan
    else:
        cyc2 = np.nan

    final = state_record(rho, Q=res_info.get("Q", np.nan))
    final.update(n=n_used, diff=diff, diff2=diff2)
    return dict(circuit=circuit, q=float(q), status=status, n_iter=n_used,
                rho=rho, final=final, residual=residual, cycle2=cyc2, history=hist,
                max_herm_err=max_herm,
                min_eig_seen=(None if min_eig_seen is np.inf else float(min_eig_seen)),
                min_Q_seen=(None if min_Q is np.inf else float(min_Q)))


# ===========================================================================
# local stability (Jacobian of the nonlinear map on trace-one Hermitian coords)
# ===========================================================================
def _rho_to_coords(rho):
    """Real coordinates r_P = Tr(P rho) for the 15 traceless Paulis (trace fixed to 1)."""
    return np.array([np.real(np.trace(PAULI_2Q[l] @ rho)) for l in PAULI_LABELS[1:]])


def _coords_to_rho(r):
    rho = PAULI_2Q["II"].astype(complex).copy()
    for c, l in zip(r, PAULI_LABELS[1:]):
        rho = rho + c * PAULI_2Q[l]
    return rho / 4.0


def coherence_amplification(rho):
    """Structural instability factor of the purification map.

    For rho -> rho^2/Tr(rho^2) a coherence |i><j| between eigenvectors with eigenvalues
    lambda_i, lambda_j is mapped to (lambda_i+lambda_j)/Tr(rho^2) times itself.  Since
    Tr(rho^2) = sum lambda^2 <= lambda_1 < lambda_1 + lambda_2, this factor EXCEEDS 1 for
    every mixed state: the mixed fixed point is a saddle, attracting inside the
    Bell-diagonal (population) manifold and expanding along eigenbasis coherences.
    Returns the largest such factor (lambda_1+lambda_2)/Tr(rho^2).
    """
    lam = np.sort(np.linalg.eigvalsh(0.5 * (rho + rho.conj().T)))[::-1]
    S = float(np.sum(lam ** 2))
    return float((lam[0] + lam[1]) / S) if S > 0 else np.nan


def solve_fixed_point(q, circuit, rho0=None, tol=1e-13):
    """Solve rho = P_q(rho) directly (Newton/hybr on the 15 Pauli coordinates).

    Root-finding does not care about linear stability, so this returns the fixed point
    even where plain iteration escapes it (see coherence_amplification).
    """
    from scipy.optimize import root
    if rho0 is None:
        rho0 = rho_isotropic(0.9)
        for _ in range(8):                       # a few iterations as a warm start
            nxt, _ = one_round_effective_map(rho0, q, circuit)
            if nxt is None:
                break
            rho0 = nxt

    def resid(r):
        nxt, _ = one_round_effective_map(_coords_to_rho(r), q, circuit)
        if nxt is None:
            return np.full(15, 1e3)
        return _rho_to_coords(nxt) - r

    sol = root(resid, _rho_to_coords(rho0), method="hybr", tol=tol)
    r = sol.x
    # Newton refinement with the central-difference Jacobian (quadratic convergence;
    # root-finding is insensitive to the saddle instability of the iteration).
    for _ in range(40):
        g = resid(r)
        if np.linalg.norm(g) < 1e-15:
            break
        Jm, _, _ = jacobian(_coords_to_rho(r), q, circuit, h=1e-6)
        # J has an exact eigenvalue 1 (marginal direction), so J - I is singular:
        # use the pseudo-inverse, which steps only in the complementary directions.
        step = np.linalg.pinv(Jm - np.eye(15), rcond=1e-8) @ g
        if not np.all(np.isfinite(step)):
            break
        r_new = r - step
        if np.linalg.norm(resid(r_new)) > np.linalg.norm(g):
            break
        r = r_new
    rho_star = _coords_to_rho(r)
    nxt, info = one_round_effective_map(rho_star, q, circuit)
    residual = float(np.linalg.norm(nxt - rho_star)) if nxt is not None else np.nan
    return rho_star, residual, bool(sol.success), info


def jacobian(rho_star, q, circuit, h=1e-6):
    """Finite-difference Jacobian J = dP_q/dr at rho_star (15x15) and its spectral radius."""
    r0 = _rho_to_coords(rho_star)
    f0, _ = one_round_effective_map(_coords_to_rho(r0), q, circuit)
    base = _rho_to_coords(f0)
    J = np.zeros((15, 15))
    for k in range(15):
        rp = r0.copy(); rp[k] += h
        rm = r0.copy(); rm[k] -= h
        fp, _ = one_round_effective_map(_coords_to_rho(rp), q, circuit)
        fm, _ = one_round_effective_map(_coords_to_rho(rm), q, circuit)
        J[:, k] = (_rho_to_coords(fp) - _rho_to_coords(fm)) / (2 * h)
    ev = np.linalg.eigvals(J)
    return J, float(np.max(np.abs(ev))), base


# ===========================================================================
# validation against the repository's verified one-round results
# ===========================================================================
def _pennylane_reference(rho, q, circuit):
    """(<Z_a>, <Z_a (x) Phi_A>) from the repository's own PennyLane executors."""
    from noisy_bell_state import O_PHI_PLUS
    dev = qml.device("default.mixed", wires=5)
    qfunc = CIRCUITS[circuit][1]

    @qml.qnode(dev)
    def node(rho_AB):
        qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])
        qfunc(q)
        return (qml.expval(qml.PauliZ(0) @ qml.Hermitian(O_PHI_PLUS, wires=[1, 2])),
                qml.expval(qml.PauliZ(0)))

    zO, zI = node(np.kron(rho, rho))
    return float(zI), float(zO)


def validate_one_round(ts=(0.95, 0.9, 0.8, 0.7, 0.5), qs=(0.0, 1e-3, 1e-2, 5e-2, 0.15),
                       verbose=True):
    """Compare the dense-matrix map against (a) the repo's PennyLane executors and
    (b) the Step-3 / Step-4 closed forms.  Returns a dict of max abs errors."""
    from pqec_cnot_threshold import Q_denom, N_num, F_dec

    err = {c: dict(Q=0.0, N=0.0, F=0.0) for c in CIRCUITS}
    err_ana = {"step3": dict(Q=0.0, N=0.0, F=0.0), "step4": dict(Q=0.0, N=0.0, F=0.0)}
    err_ideal = 0.0

    for t in ts:
        rho = rho_isotropic(t)
        for q in qs:
            for c in CIRCUITS:
                tau = one_round_tau(rho, q, c)
                Q = float(np.real(np.trace(tau)))
                N = float(np.real(np.trace(PHI @ tau)))
                zI, zO = _pennylane_reference(rho, q, c)
                err[c]["Q"] = max(err[c]["Q"], abs(Q - zI))
                err[c]["N"] = max(err[c]["N"], abs(N - zO))
                err[c]["F"] = max(err[c]["F"], abs(N / Q - zO / zI))

            # Step 3 closed form (pqec_cnot_threshold)
            tau3 = one_round_tau(rho, q, "step3")
            Q3 = float(np.real(np.trace(tau3)))
            N3 = float(np.real(np.trace(PHI @ tau3)))
            err_ana["step3"]["Q"] = max(err_ana["step3"]["Q"], abs(Q3 - Q_denom(t, q)))
            err_ana["step3"]["N"] = max(err_ana["step3"]["N"], abs(N3 - N_num(t, q)))
            err_ana["step3"]["F"] = max(err_ana["step3"]["F"], abs(N3 / Q3 - F_dec(q, t)))

            # Step 4 closed form (write-up):  s = 1-q
            s = 1 - q
            Q4a = s ** 10 / 4 * (1 + 3 * s ** 2 * t ** 2)
            N4a = s ** 10 / 16 * (1 + 3 * s ** 2 * t ** 2 + s ** 3 * (1 + 5 * s) * t * (1 + t))
            tau4 = one_round_tau(rho, q, "step4")
            Q4 = float(np.real(np.trace(tau4)))
            N4 = float(np.real(np.trace(PHI @ tau4)))
            err_ana["step4"]["Q"] = max(err_ana["step4"]["Q"], abs(Q4 - Q4a))
            err_ana["step4"]["N"] = max(err_ana["step4"]["N"], abs(N4 - N4a))
            err_ana["step4"]["F"] = max(err_ana["step4"]["F"], abs(N4 / Q4 - N4a / Q4a))

        # ideal map must be exactly rho^2 / Tr(rho^2)
        for c in CIRCUITS:
            tau0 = one_round_tau(rho, 0.0, c)
            err_ideal = max(err_ideal, float(np.linalg.norm(tau0 - rho @ rho)))

    if verbose:
        print("  one-round validation (dense map vs repository executors)")
        for c in CIRCUITS:
            print(f"    {c}: max|dQ| = {err[c]['Q']:.2e}, max|dN| = {err[c]['N']:.2e}, "
                  f"max|dF| = {err[c]['F']:.2e}")
        print("  one-round validation (dense map vs closed forms)")
        for c in err_ana:
            print(f"    {c}: max|dQ| = {err_ana[c]['Q']:.2e}, max|dN| = {err_ana[c]['N']:.2e}, "
                  f"max|dF| = {err_ana[c]['F']:.2e}")
        print(f"  ideal map  tau_A == rho^2 :  max||.|| = {err_ideal:.2e}")
    return dict(vs_pennylane=err, vs_analytic=err_ana, ideal=err_ideal)


def selftest(verbose=True):
    """Channel definition vs the repository Kraus channel, plus basic sanity."""
    from noisy_bell_state import global_depol_kraus
    rng = np.random.default_rng(0)
    A = rng.normal(size=(DIM, DIM)) + 1j * rng.normal(size=(DIM, DIM))
    rho = A @ A.conj().T
    rho /= np.trace(rho)
    worst = 0.0
    for q in (0.0, 1e-3, 0.05, 0.3, 1.0):
        for wires in [(0, 1), (1, 3), (2, 4), (0, 4), (3, 1)]:
            mine = replacement_depol(rho, wires, q)
            K = [_embed(k, wires) for k in global_depol_kraus(q)]
            kra = sum(k @ rho @ k.conj().T for k in K)
            worst = max(worst, float(np.linalg.norm(mine - kra)))
    if verbose:
        print(f"  replacement channel (definition) vs global_depol_kraus: "
              f"max||.|| = {worst:.2e}")
    for c in CIRCUITS:
        assert n_cnots(c) == (16 if c == "step3" else 14), c
    if verbose:
        print(f"  CNOT counts: " + ", ".join(f"{c}={n_cnots(c)}" for c in CIRCUITS))
    return worst


if __name__ == "__main__":
    print("=" * 78)
    print(" Iterated noisy PQEC -- self-test and one-round validation")
    print("=" * 78)
    selftest()
    validate_one_round()
