"""Additional reproducible consistency tests for BBPSSW circuit/reference."""
import numpy as np
import pytest
from bbpssw_map import (
    isotropic_input, first_round_isotropic, round_twirl, round_untwirled, round_original_bbpssw, trajectory
)
from frozen_circuit_reference import (
    independent_circuit_round, noisy_cnot_channel, cnot_permutation,
    bell_projectors, physical_twirl, clifford_twirl_ops
)


def partial_replacement_independent(rho, pair):
    keep=[w for w in range(4) if w not in pair]
    def bit(j,w): return (j>>(3-w))&1
    def repl(j,bits,ws):
        for k,w in enumerate(ws):
            mask=1<<(3-w)
            if (bits>>k)&1: j|=mask
            else: j&=~mask
        return j
    s=np.zeros_like(rho)
    for i in range(16):
        for j in range(16):
            if any(bit(i,w)!=bit(j,w) for w in pair): continue
            for k in range(4):
                s[i,j] += rho[repl(i,k,pair),repl(j,k,pair)]/4
    return s


def test_frozen_gates_have_correct_controls_and_are_unitary():
    for c,t in ((0,2),(1,3)):
        U=cnot_permutation(c,t)
        assert np.max(abs(U.conj().T@U-np.eye(16)))==0
        assert np.max(abs(U@U-np.eye(16)))==0


@pytest.mark.parametrize('pair',[(0,2),(1,3)])
@pytest.mark.parametrize('q',[0,0.07,0.5,1])
def test_pauli_conjugation_implements_replacement(pair,q):
    rng=np.random.default_rng(16)
    for _ in range(3):
        M=rng.normal(size=(16,16))+1j*rng.normal(size=(16,16))
        rho=M@M.conj().T;rho=rho/np.trace(rho)
        U=cnot_permutation(*pair)
        ideal=U@rho@U.conj().T
        wanted=(1-q)*ideal+q*partial_replacement_independent(ideal,pair)
        got=noisy_cnot_channel(rho,q,pair)
        assert np.max(abs(got-wanted))<1e-14


def test_twirl_uniform_clifford_reproduces_isotropic_density():
    rng=np.random.default_rng(17)
    assert len(clifford_twirl_ops())==24
    Phi=bell_projectors()[0]
    for _ in range(5):
        M=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        rho=M@M.conj().T;rho=rho/np.trace(rho)
        F=np.trace(Phi@rho).real
        expected=F*Phi+(1-F)*(np.eye(4)-Phi)/3
        assert np.max(abs(physical_twirl(rho)-expected))<2e-14


@pytest.mark.parametrize('epsilon',[0,0.00001,0.01,0.1,0.4,2/3,0.9,1])
@pytest.mark.parametrize('q',[0,1e-7,0.01,0.1,0.4,0.8,1])
def test_analytical_map_matches_frozen_circuit(epsilon,q):
    p=isotropic_input(epsilon)
    expected,P=round_twirl(p,q)
    got=independent_circuit_round(p,q,twirl=True)
    assert np.max(abs(np.asarray(expected)-got['populations']))<2e-14
    assert abs(P-got['success'])<2e-14
    f_scalar,P_scalar=first_round_isotropic(epsilon,q)
    assert abs(f_scalar-expected[0])<1e-14
    assert abs(P_scalar-P)<1e-14


def test_twirl_placement_important_for_non_isotropic_input():
    p=(0.67,0.03,0.19,0.11)
    f_nopre=round_twirl(p,0.10)[0][0]
    F=p[0]
    f_pretwirl=round_original_bbpssw(p,0.10)[0][0]
    assert abs(f_nopre-f_pretwirl)>1e-3
    phys=independent_circuit_round(p,0.10,twirl=True,pre_twirl=True)
    assert abs(f_pretwirl-phys["populations"][0])<1e-14


def test_twirl_preserves_phi_fidelity_and_success():
    p=(0.37,0.19,0.23,0.21)
    for q in (0,0.005,0.15,0.8,1):
        pre,p0=round_untwirled(p,q)
        post,p1=round_twirl(p,q)
        assert abs(pre[0]-post[0])<1e-14
        assert p0==p1
        assert max(abs(post[j]-post[k]) for j,k in ((1,2),(2,3)))<1e-14


def test_q1_is_maximally_mixed_and_success_half():
    for p in [(1,0,0,0),(.45,.03,.15,.37),(.25,)*4]:
        pops,success=round_twirl(p,1)
        assert np.allclose(pops,[.25]*4,atol=1e-15)
        assert success==.5


def test_ideal_agrees_with_original_bennett_recurrence():
    for F in np.linspace(0.25,1,17):
        p=(F,(1-F)/3,(1-F)/3,(1-F)/3)
        out,s=round_twirl(p,0)
        G=(F*F+((1-F)/3)**2)/(F*F+2*F*(1-F)/3+5*((1-F)/3)**2)
        assert abs(out[0]-G)<1e-14


def test_optimal_q_first_gain_isotropic_formula():
    for epsilon in (0.01,0.1,0.2,0.4,0.6):
        q0=1-np.sqrt(3/(3+2*epsilon-3*epsilon**2))
        F=1-3*epsilon/4
        fout,ps=first_round_isotropic(epsilon,q0)
        assert abs(fout-F)<1e-13
        assert first_round_isotropic(epsilon,q0*0.8)[0]>F
        assert first_round_isotropic(epsilon,q0*1.2)[0]<F


def test_repeated_input_remains_isotropic():
    for e in (0,0.05,0.4,0.8,1):
        for q in (0,0.05,0.12,0.4,1):
            for item in trajectory(e,q,25):
                a,b,c,d=item['populations']
                assert abs(a+b+c+d-1)<1e-13
                assert min(a,b,c,d)>-1e-13
                assert b==c==d

@pytest.mark.parametrize('q',[0,0.01,0.13,0.42,1])
def test_original_full_bbpssw_arbitrary_bell_diagonal(q):
    p=(0.62,0.03,0.25,0.10)
    analytical,success=round_original_bbpssw(p,q)
    frozen=independent_circuit_round(p,q,twirl=True,pre_twirl=True)
    assert np.max(abs(np.asarray(analytical)-frozen['populations']))<1e-14
    assert abs(success-frozen['success'])<1e-14
