"""
Shared infrastructure for Step 5 -- variational (PQC) approximation of the gadget.
=================================================================================

The one-round PQEC gadget is the fixed 5-qubit unitary

    U = H_a . CSWAP(a;A1,B1) . CSWAP(a;A2,B2) . H_a
        (wire 0 = ancilla a, register A = wires 1,2 [kept], register B = wires 3,4)

Step 5 (and two auxiliary relaxations) ask whether a *parameterized* quantum circuit
(PQC) with FEWER CNOTs can play the same role, in three increasingly relaxed senses
(mirroring Step 4):

  * Step 5 (main)   -- reproduce the FULL 5-qubit unitary U    (all 32 columns)
  * isometry (aux)  -- reproduce U on the PHYSICAL input subspace (ancilla |0>, 16 cols)
                       i.e. keep the coherent purified state that gets fed forward
  * observable (aux)-- noise-aware training: keep only the purified OBSERVABLE
                       F = Tr(O rho^2)/Tr(rho^2), with CNOT depolarizing IN the loss

This module holds everything shared by those three scripts:

  * the target unitary U and its ancilla-|0> restriction U0 (32x16 isometry),
  * a hardware-efficient ansatz whose CNOT budget B is a knob,
  * a FAST pure-numpy executor (reshape-based gate application) that returns the
    ansatz unitary V(theta) -- used for the unitary/coherent-state costs,
  * a PennyLane default.mixed executor of the SAME ansatz with a 2-qubit
    depolarizing channel of strength eps2 after each CNOT -- used for the noisy
    observable read-out (auxiliary observable study),
  * the ancilla-parity read-out F = <Z_a (x) O>/<Z_a> from any 5-qubit unitary,
  * exact references: F_exact(eps) (ideal purified fidelity) and F_bare(eps).

All three notions of "same role" reuse the Setup/2/3 conventions exactly:
rho_eps = (1-eps)|Phi+><Phi+| + eps I/4 (t = 1-eps), the kept register is A = wires
1,2, the observable is O = |Phi+><Phi+|, and the CNOT channel is the replacement
2-qubit global depolarizing  (1-eps2) rho + eps2 I/4  (same as global_depol_kraus).

The qml.Rot convention (RZ(w) RY(t) RZ(p)) is matched bit-for-bit by _rot_np so the
numpy unitary and the PennyLane circuit implement the identical map (checked in the
self-test below).

Run:  python pqc_common.py     # self-tests (numpy<->PennyLane, U read-out, U0)
"""

import numpy as np
import pennylane as qml

from noisy_bell_state import rho_eps_analytic, O_PHI_PLUS, PHI_PLUS

# ---------------------------------------------------------------------------
# Fixed single-qubit matrices
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)
_PAULI = (I2, _X, _Y, _Z)

N_WIRES = 5
DIM = 2 ** N_WIRES


# ---------------------------------------------------------------------------
# Target unitary  U = H_a CSWAP(a;1,3) CSWAP(a;2,4) H_a  (wire 0 = ancilla)
# ---------------------------------------------------------------------------
def _target_qfunc():
    qml.Hadamard(0)
    qml.ctrl(qml.SWAP, control=0)(wires=[1, 3])
    qml.ctrl(qml.SWAP, control=0)(wires=[2, 4])
    qml.Hadamard(0)


U_TARGET = qml.matrix(_target_qfunc, wire_order=range(N_WIRES))()

# ancilla-|0> input isometry:  columns of I_32 whose wire-0 bit is 0.
# wire 0 is the most significant qubit, so those are basis indices 0..15.
_ANC0_COLS = [i for i in range(DIM) if ((i >> (N_WIRES - 1)) & 1) == 0]
U0_TARGET = U_TARGET[:, _ANC0_COLS]                      # 32 x 16 isometry


# ---------------------------------------------------------------------------
# Hardware-efficient ansatz  (CNOT budget B is the knob)
# ---------------------------------------------------------------------------
# Ancilla-centric connectivity: the gadget entangles the ancilla with each pair
# and swaps (1,3),(2,4).  We cycle a schedule that offers exactly those links.
_PAIR_CYCLE = [(0, 1), (0, 3), (0, 2), (0, 4), (1, 3), (2, 4), (1, 2), (3, 4)]


def ansatz_ops(budget):
    """Return (ops, n_params) for a CNOT budget `budget`.

    ops is a flat list of gate specs executed left-to-right:
        ('rot', wire, p)   -- qml.Rot(params[p:p+3], wire)
        ('cnot', c, t)     -- CNOT(c->t)  (noise is added by the executor)
    Layout: full Rot layer, then for each CNOT a Rot on its two wires, then a
    final full Rot layer.  n_params = 30 + 6*budget.
    """
    ops, p = [], 0
    for w in range(N_WIRES):                    # initial full rotation layer
        ops.append(("rot", w, p)); p += 3
    for k in range(budget):
        c, t = _PAIR_CYCLE[k % len(_PAIR_CYCLE)]
        ops.append(("cnot", c, t))
        ops.append(("rot", c, p)); p += 3
        ops.append(("rot", t, p)); p += 3
    for w in range(N_WIRES):                    # final full rotation layer
        ops.append(("rot", w, p)); p += 3
    return ops, p


# ---------------------------------------------------------------------------
# FAST numpy executor: build the ansatz unitary V(theta)  (no PennyLane)
# ---------------------------------------------------------------------------
def _rot_np(phi, theta, omega):
    """qml.Rot convention: RZ(omega) RY(theta) RZ(phi)."""
    cz1 = np.exp(-0.5j * phi); cz3 = np.exp(-0.5j * omega)
    rz1 = np.array([[cz1, 0], [0, np.conj(cz1)]], dtype=complex)
    rz3 = np.array([[cz3, 0], [0, np.conj(cz3)]], dtype=complex)
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    ry = np.array([[c, -s], [s, c]], dtype=complex)
    return rz3 @ ry @ rz1


def _cnot_perm(c, t):
    """Row permutation implementing left-multiplication by CNOT(c->t) on 5 qubits."""
    perm = np.arange(DIM)
    cb, tb = N_WIRES - 1 - c, N_WIRES - 1 - t     # bit positions (wire 0 = MSB)
    for i in range(DIM):
        if (i >> cb) & 1:
            perm[i] = i ^ (1 << tb)
    return perm


_CNOT_PERM_CACHE = {}


def _apply_1q_left(V, g, w):
    """Return (G_w (x) I) @ V for single-qubit g on wire w. V is (DIM, cols)."""
    cols = V.shape[1]
    Vt = V.reshape([2] * N_WIRES + [cols])
    Vt = np.tensordot(g, Vt, axes=([1], [w]))     # new axis goes to front
    Vt = np.moveaxis(Vt, 0, w)
    return Vt.reshape(DIM, cols)


def _apply_1q_right(M, g, w):
    """Return M @ (g_w (x) I) for single-qubit g on wire w. M is (rows, DIM)."""
    rows = M.shape[0]
    Mt = M.reshape([rows] + [2] * N_WIRES)
    Mt = np.tensordot(Mt, g, axes=([1 + w], [0]))  # contract col-axis w; new axis last
    Mt = np.moveaxis(Mt, -1, 1 + w)
    return Mt.reshape(rows, DIM)


def _embed_1q(g, w):
    """Full 32x32 embedding of single-qubit g on wire w."""
    return _apply_1q_left(np.eye(DIM, dtype=complex), g, w)


# --- Rot derivatives (qml.Rot = RZ(omega) RY(theta) RZ(phi)) -----------------
def _rz(a):
    return np.array([[np.exp(-0.5j * a), 0], [0, np.exp(0.5j * a)]], dtype=complex)


def _ry(b):
    c, s = np.cos(b / 2), np.sin(b / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def _drot(phi, theta, omega):
    """Return (dRot/dphi, dRot/dtheta, dRot/domega) as 2x2 matrices."""
    RZp, RY, RZo = _rz(phi), _ry(theta), _rz(omega)
    dRZp = (-0.5j) * (_Z @ RZp)
    dRZo = (-0.5j) * (_Z @ RZo)
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    dRY = 0.5 * np.array([[-s, -c], [c, -s]], dtype=complex)
    return (RZo @ RY @ dRZp, RZo @ dRY @ RZp, dRZo @ RY @ RZp)


def unitary_from_params(ops, params, cols_isometry=None):
    """Build V(theta) (or V restricted to given input columns) as a dense matrix.

    cols_isometry: if None, start from I_32 (full 32x32 unitary). Otherwise start
    from that 32xk isometry (e.g. the ancilla-|0> columns) and return 32xk -- this
    is exactly V0 = V @ E used by the isometry (auxiliary) cost.
    """
    params = np.asarray(params, dtype=float)
    V = np.eye(DIM, dtype=complex) if cols_isometry is None \
        else np.asarray(cols_isometry, dtype=complex)
    for op in ops:
        if op[0] == "rot":
            _, w, p = op
            V = _apply_1q_left(V, _rot_np(params[p], params[p + 1], params[p + 2]), w)
        else:
            _, c, t = op
            key = (c, t)
            if key not in _CNOT_PERM_CACHE:
                _CNOT_PERM_CACHE[key] = _cnot_perm(c, t)
            V = V[_CNOT_PERM_CACHE[key], :]
    return V


# ---------------------------------------------------------------------------
# Read-out F = <Z_a (x) O> / <Z_a>  from an arbitrary 5-qubit unitary W
# ---------------------------------------------------------------------------
def _rho_in(eps):
    """|0><0|_a (x) rho_eps(1,2) (x) rho_eps(3,4)  as a 32x32 density matrix."""
    r = rho_eps_analytic(eps)
    ket0 = np.array([[1, 0], [0, 0]], dtype=complex)
    return np.kron(ket0, np.kron(r, r))


def _obs_ZO(O):
    """Operator  Z_0 (x) O_{12} (x) I_{34}  as 32x32."""
    return np.kron(_Z, np.kron(O, np.eye(4, dtype=complex)))


_OBS_Z = np.kron(_Z, np.eye(16, dtype=complex))          # Z_0 (x) I


def F_from_unitary(W, eps, O=O_PHI_PLUS):
    """Purified observable read out from a *noiseless* 5-qubit unitary W."""
    rho_out = W @ _rho_in(eps) @ W.conj().T
    zO = np.real(np.trace(_obs_ZO(O) @ rho_out))
    zI = np.real(np.trace(_OBS_Z @ rho_out))
    return float(zO / zI)


# ---------------------------------------------------------------------------
# Exact references
# ---------------------------------------------------------------------------
def F_exact(eps, O=O_PHI_PLUS):
    """Ideal purified observable  Tr(O rho_eps^2)/Tr(rho_eps^2)."""
    r = rho_eps_analytic(eps)
    r2 = r @ r
    return float(np.real(np.trace(O @ r2) / np.trace(r2)))


def F_bare(eps):
    """Bare (no-purification) Bell fidelity of rho_eps:  (1+3t)/4,  t = 1-eps."""
    return (1 + 3 * (1 - eps)) / 4


# ---------------------------------------------------------------------------
# PennyLane default.mixed executor of the SAME ansatz, with CNOT depol eps2
# ---------------------------------------------------------------------------
from noisy_bell_state import global_depol_kraus     # noqa: E402

_dev5 = qml.device("default.mixed", wires=N_WIRES)


def apply_ansatz_qml(ops, params, eps2):
    """Apply the ansatz on the current PennyLane tape; depol(eps2) after each CNOT."""
    for op in ops:
        if op[0] == "rot":
            _, w, p = op
            qml.Rot(params[p], params[p + 1], params[p + 2], wires=w)
        else:
            _, c, t = op
            qml.CNOT(wires=[c, t])
            if eps2 > 0:
                qml.QubitChannel(global_depol_kraus(eps2), wires=[c, t])


def make_F_qnode(ops):
    """Return f(params, eps, eps2, O) -> (<Z_a (x) O>, <Z_a>) for this ansatz."""
    @qml.qnode(_dev5)
    def _f(params, eps, eps2, O):
        r = rho_eps_analytic(eps)
        qml.QubitDensityMatrix(np.kron(r, r), wires=[1, 2, 3, 4])
        apply_ansatz_qml(ops, params, eps2)
        return (qml.expval(qml.PauliZ(0) @ qml.Hermitian(O, wires=[1, 2])),
                qml.expval(qml.PauliZ(0)))
    return _f


# ---------------------------------------------------------------------------
# Costs (Hilbert-Schmidt test)  -- 0 iff exact reproduction (up to global phase)
# ---------------------------------------------------------------------------
def cost_unitary(ops, params):
    """Step 5 cost:  1 - |Tr(U^dag V)|^2 / 32^2."""
    V = unitary_from_params(ops, params)
    ov = np.vdot(U_TARGET, V)                     # Tr(U^dag V)
    return 1.0 - (abs(ov) ** 2) / DIM ** 2


# isometry E (32x16) embedding ancilla-|0> inputs: columns of I_32 at _ANC0_COLS
_ANC0_ISO = np.eye(DIM, dtype=complex)[:, _ANC0_COLS]


def cost_coherent(ops, params):
    """Isometry (auxiliary) cost:  1 - |Tr(U0^dag V0)|^2 / 16^2   (ancilla-|0> input block)."""
    V0 = unitary_from_params(ops, params, cols_isometry=_ANC0_ISO)
    ov = np.vdot(U0_TARGET, V0)
    return 1.0 - (abs(ov) ** 2) / (len(_ANC0_COLS) ** 2)


# ---------------------------------------------------------------------------
# Exact analytic (cost, gradient) via backprop through the gate product.
# The Hilbert-Schmidt cost 1-|Tr(T^dag V E)|^2/dc^2 is a GLOBAL cost with barren
# plateaus; finite-difference gradients stall, so we differentiate it exactly.
# target T (32 x dc), input isometry E (32 x dc):  T=U,E=I (unitary);  T=U0,E=anc0 (isometry).
# ---------------------------------------------------------------------------
def _apply_left(op, M, params):
    if op[0] == "rot":
        _, w, p = op
        return _apply_1q_left(M, _rot_np(params[p], params[p + 1], params[p + 2]), w)
    _, c, t = op
    key = (c, t)
    if key not in _CNOT_PERM_CACHE:
        _CNOT_PERM_CACHE[key] = _cnot_perm(c, t)
    return M[_CNOT_PERM_CACHE[key], :]


def _apply_right(op, M, params):
    if op[0] == "rot":
        _, w, p = op
        return _apply_1q_right(M, _rot_np(params[p], params[p + 1], params[p + 2]), w)
    _, c, t = op
    key = (c, t)
    if key not in _CNOT_PERM_CACHE:
        _CNOT_PERM_CACHE[key] = _cnot_perm(c, t)
    return M[:, _CNOT_PERM_CACHE[key]]


def _cost_grad(ops, params, target, E):
    """Return (cost, grad) for  1 - |Tr(target^dag V E)|^2 / dc^2,  V = ansatz(params)."""
    params = np.asarray(params, dtype=float)
    m = len(ops)
    dc = target.shape[1]

    # prefixes P[i] = g_{i-1}..g_0   (i gates), P[0]=I, P[m]=V
    P = [np.eye(DIM, dtype=complex)]
    for op in ops:
        P.append(_apply_left(op, P[-1], params))
    V = P[m]
    Veff = V @ E
    g = np.vdot(target, Veff)                       # Tr(target^dag V E)
    cost = 1.0 - (abs(g) ** 2) / dc ** 2

    # suffixes Suf[i] = g_{m-1}..g_i,  Suf[m]=I
    Suf = [None] * (m + 1)
    Suf[m] = np.eye(DIM, dtype=complex)
    for i in range(m - 1, -1, -1):
        Suf[i] = _apply_right(ops[i], Suf[i + 1], params)

    Tdag = target.conj().T                          # dc x 32
    pref = -2.0 / dc ** 2
    grad = np.zeros_like(params)
    for i, op in enumerate(ops):
        if op[0] != "rot":
            continue
        _, w, p = op
        # A_i = (P[i] E) Tdag Suf[i+1]   (32x32);  dg/dparam = Tr(A_i dG_i)
        A = (P[i] @ E) @ Tdag @ Suf[i + 1]
        for k, dg2 in enumerate(_drot(params[p], params[p + 1], params[p + 2])):
            dgval = np.trace(A @ _embed_1q(dg2, w))
            grad[p + k] = pref * np.real(np.conj(g) * dgval)
    return cost, grad


def cost_grad_unitary(ops, params):
    """(cost, grad) for Step 5 (full unitary; E=I, target=U)."""
    return _cost_grad(ops, params, U_TARGET, np.eye(DIM, dtype=complex))


# --- Local Hilbert-Schmidt-Test (LHST) cost for Step 5 ---------------------
# The global cost 1-|Tr(U^dag V)|^2/d^2 is prone to barren plateaus; LHST (Khatri et al.,
# Quantum 3, 140 (2019)) replaces the global Bell overlap of the Choi state of
# W = V^dag U by an average of per-qubit Bell overlaps.  Closed operator form
# (derived in PQC_APPROX.md):
#     C_LHST = 3/4 - (1/(4 d n)) Re Tr(T M^dag),   M = V^dag U,
#     T = sum_{j=1..n} sum_{P in {X,Y,Z}}  P_j M P_j.
# C_LHST = 0 iff W = I up to phase (exact compilation), and its gradient does not
# vanish exponentially, so it trains where the global cost stalls.
_PAULI_EMB = [[_embed_1q(P, w) for P in (_X, _Y, _Z)] for w in range(N_WIRES)]
_U_DAG = U_TARGET.conj().T


def cost_grad_lhst(ops, params):
    """(cost, grad) for the LHST (local) compiling cost of the full unitary U."""
    params = np.asarray(params, dtype=float)
    m = len(ops)
    d, n = DIM, N_WIRES

    P = [np.eye(DIM, dtype=complex)]
    for op in ops:
        P.append(_apply_left(op, P[-1], params))
    V = P[m]
    M = V.conj().T @ U_TARGET                       # W = V^dag U
    # twirl sum  T = sum_{j,P} P_j M P_j
    T = np.zeros((DIM, DIM), dtype=complex)
    for w in range(n):
        for Pj in _PAULI_EMB[w]:
            T += Pj @ M @ Pj
    cost = 0.75 - np.real(np.trace(T @ M.conj().T)) / (4.0 * d * n)

    Suf = [None] * (m + 1)
    Suf[m] = np.eye(DIM, dtype=complex)
    for i in range(m - 1, -1, -1):
        Suf[i] = _apply_right(ops[i], Suf[i + 1], params)

    TUd = T @ _U_DAG
    coef = -2.0 / (4.0 * d * n)
    grad = np.zeros_like(params)
    for i, op in enumerate(ops):
        if op[0] != "rot":
            continue
        _, w, p = op
        Msum = P[i] @ TUd @ Suf[i + 1]
        for k, dg2 in enumerate(_drot(params[p], params[p + 1], params[p + 2])):
            grad[p + k] = coef * np.real(np.trace(Msum @ _embed_1q(dg2, w)))
    return cost, grad


def cost_grad_coherent(ops, params):
    """(cost, grad) for the isometry, auxiliary (ancilla-|0> block; E=anc0 isometry, target=U0)."""
    return _cost_grad(ops, params, U0_TARGET, _ANC0_ISO)


# ===========================================================================
def _selftest():
    print("=" * 74)
    print(" pqc_common self-test")
    print("=" * 74)
    rng = np.random.default_rng(0)

    # (1) target U reproduces the exact purified observable for several eps
    err = 0.0
    for eps in [0.1, 0.2, 0.35, 0.5, 0.6]:
        err = max(err, abs(F_from_unitary(U_TARGET, eps) - F_exact(eps)))
    print(f"\n (1) read-out from U  ==  F_exact          : max err {err:.2e}")

    # (2) numpy unitary  ==  PennyLane unitary  for the SAME random params
    ops, npar = ansatz_ops(6)
    th = rng.uniform(-np.pi, np.pi, size=npar)
    V_np = unitary_from_params(ops, th)
    V_ql = qml.matrix(lambda: apply_ansatz_qml(ops, th, 0.0), wire_order=range(N_WIRES))()
    # compare up to nothing (should be identical, same gate order/convention)
    print(f" (2) numpy V  ==  PennyLane V (B=6)        : max err "
          f"{np.max(np.abs(V_np - V_ql)):.2e}   (npar={npar})")

    # (3) V0 from isometry  ==  full V restricted to ancilla-|0> columns
    V0a = unitary_from_params(ops, th, cols_isometry=_ANC0_ISO)
    V0b = V_np[:, _ANC0_COLS]
    print(f" (3) V0 (isometry) == V[:,anc0]           : max err "
          f"{np.max(np.abs(V0a - V0b)):.2e}")

    # (4) costs vanish at the exact gadget: fit is trivially the target itself.
    #     (sanity) cost_unitary(U) via a params-free check using U directly:
    ovU = np.vdot(U_TARGET, U_TARGET)
    print(f" (4) |Tr(U^dag U)|^2/32^2 = {abs(ovU)**2/DIM**2:.6f}  (=1 -> cost 0)")

    # (5) noisy read-out at eps2=0 matches the noiseless read-out
    f = make_F_qnode(ops)
    zO, zI = f(th, 0.3, 0.0, O_PHI_PLUS)
    F_qml = float(zO) / float(zI)
    F_np = F_from_unitary(V_np, 0.3)
    print(f" (5) noisy qnode(eps2=0) == numpy read-out: err "
          f"{abs(F_qml - F_np):.2e}")

    # (6) exact references
    print("\n (6) references (O = |Phi+><Phi+|):")
    print(f"     {'eps':>5} {'F_bare':>8} {'F_exact':>8}")
    for eps in [0.1, 0.2, 0.4, 0.6]:
        print(f"     {eps:>5.2f} {F_bare(eps):>8.4f} {F_exact(eps):>8.4f}")
    print("\n  self-test done.")


if __name__ == "__main__":
    _selftest()
