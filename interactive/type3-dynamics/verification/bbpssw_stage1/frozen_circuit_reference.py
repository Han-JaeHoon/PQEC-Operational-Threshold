"""Independent, frozen, 4-qubit circuit-level BBPSSW reference.

No imports from the analytical recurrence module.  The depolarizing channel
is independently realized by averaging conjugations by all 16 Pauli products
on each two-qubit gate support (not by partial trace/replacement helper).
The twirling stage is realized by an explicit uniform 24-element Clifford
U (x) U* ensemble rather than by editing the Bell probabilities.
"""

from __future__ import annotations
import numpy as np
from functools import lru_cache

WIRE_ORDER = ('A1', 'B1', 'A2', 'B2')
BXOR_GATES = ((0,2), (1,3))
NOISE_AFTER_EACH = True
MEASUREMENT_WIRES = (2,3)
ACCEPTED_READOUTS = ((0,0),(1,1))
BELL_ORDER = ('Phi+', 'Phi-', 'Psi+', 'Psi-')

I = np.eye(2,dtype=complex)
X = np.array([[0,1],[1,0]],dtype=complex)
Y = np.array([[0,-1j],[1j,0]],dtype=complex)
Z = np.diag([1,-1]).astype(complex)
H = np.array([[1,1],[1,-1]],dtype=complex)/np.sqrt(2)
S = np.diag([1,1j]).astype(complex)
PAULIS = (I,X,Y,Z)


def kron_all(parts):
    result=parts[0]
    for p in parts[1:]:
        result=np.kron(result,p)
    return result


def bell_projectors():
    vectors=(np.array([1,0,0,1],complex),np.array([1,0,0,-1],complex),
             np.array([0,1,1,0],complex),np.array([0,1,-1,0],complex))
    return tuple(np.outer(v,v.conj())/2 for v in vectors)


def bell_density(probabilities):
    return sum(float(p)*P for p,P in zip(probabilities,bell_projectors()))


def cnot_permutation(control, target):
    n=4
    permutation=np.zeros((1<<n,1<<n),dtype=complex)
    for col in range(1<<n):
        bits=[(col>>(n-1-i))&1 for i in range(n)]
        if bits[control]:
            bits[target] ^= 1
        row=sum(b << (n-1-i) for i,b in enumerate(bits))
        permutation[row,col]=1
    return permutation


@lru_cache(None)
def pauli_operators_for_pair(wire0,wire1):
    ops=[]
    for op0 in PAULIS:
        for op1 in PAULIS:
            ops.append(kron_all([op0 if i==wire0 else op1 if i==wire1 else I
                                 for i in range(4)]))
    return tuple(ops)


def noisy_cnot_channel(sigma, noise_q, controls):
    control, target=controls
    U=cnot_permutation(control,target)
    out=U@sigma@U.conj().T
    if noise_q == 0.0:
        return out
    total=np.zeros_like(out)
    for P in pauli_operators_for_pair(control,target):
        total += P @ out @ P.conj().T
    return (1-noise_q)*out+(noise_q/16.0)*total


def _global_phase_key(U):
    vec=U.flatten()
    z=next(z for z in vec if abs(z)>1e-9)
    W=U/(z/abs(z))
    return tuple(np.round(W.real.flatten(),11))+tuple(np.round(W.imag.flatten(),11))


@lru_cache(None)
def clifford_twirl_ops():
    # Generate projective single-qubit Clifford group through H and S.
    queue=[I]
    seen={_global_phase_key(I)}
    for m in queue:
        for generator in (H,S):
            new=generator@m
            key=_global_phase_key(new)
            if key not in seen:
                queue.append(new)
                seen.add(key)
    if len(queue)!=24:
        raise AssertionError('Clifford group is not 24 elements: '+str(len(queue)))
    return tuple(np.kron(U,U.conj()) for U in queue)


def physical_twirl(sigma):
    total=np.zeros_like(sigma)
    for V in clifford_twirl_ops():
        total += V@sigma@V.conj().T
    return total/len(clifford_twirl_ops())


def postselect_success(sigma, outcomes=ACCEPTED_READOUTS):
    # Trace over Z-measured target qubits, keeping 00 and 11 incoherently.
    source=np.zeros((4,4),dtype=complex)
    for a2,b2 in outcomes:
        for src_row in range(4):
            a1,b1=divmod(src_row,2)
            r=8*a1+4*b1+2*a2+b2
            for src_col in range(4):
                aa1,bb1=divmod(src_col,2)
                c=8*aa1+4*bb1+2*a2+b2
                source[src_row,src_col]+=sigma[r,c]
    return source


def independent_circuit_round(p,q,*,twirl=True,pre_twirl=False,bob_noise_override=None,
                              reverse_bob=False,accepted_outcomes=ACCEPTED_READOUTS):
    if not (0 <= q <= 1): raise ValueError('q out of range')
    if bob_noise_override is not None and not 0 <= bob_noise_override <=1:
        raise ValueError('bob noise out of range')
    initial=bell_density(p)
    if pre_twirl:
        initial=physical_twirl(initial)
    state=np.kron(initial,initial)
    # A1 -> A2, noise on (A1,A2)
    state=noisy_cnot_channel(state,q,BXOR_GATES[0])
    # B1 -> B2, noise on (B1,B2)
    bob_pair=(3,1) if reverse_bob else BXOR_GATES[1]
    state=noisy_cnot_channel(state,q if bob_noise_override is None else bob_noise_override,bob_pair)
    unnorm=postselect_success(state,accepted_outcomes)
    success=float(np.trace(unnorm).real)
    conditional=unnorm/success
    final=physical_twirl(conditional) if twirl else conditional
    pops=np.array([float(np.trace(P@final).real) for P in bell_projectors()])
    return {'success':success,'populations':pops,'conditional':conditional,
            'final':final,'full_state':state}
