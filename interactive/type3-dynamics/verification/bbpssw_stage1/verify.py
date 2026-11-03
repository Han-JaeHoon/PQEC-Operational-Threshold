"""Numerical study for the analytical BBPSSW noisy map and independent frozen circuit."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from bbpssw_map import isotropic_input, round_twirl, round_untwirled, round_original_bbpssw, first_round_isotropic
from frozen_circuit_reference import (
    independent_circuit_round, bell_projectors, physical_twirl,
    clifford_twirl_ops, WIRE_ORDER, BXOR_GATES, MEASUREMENT_WIRES,
    ACCEPTED_READOUTS, BELL_ORDER,
)

OUTPUT=Path(__file__).parent
RNG=np.random.default_rng(20261103)


def emax(a,b):
    return float(np.max(np.abs(np.asarray(a)-np.asarray(b))))


def test_families():
    cases=[]
    for j in range(4):
        vec=np.eye(4)[j]
        cases.append(('bell_vertex_%d'%j,tuple(vec)))
    for i in range(4):
        for j in range(i+1,4):
            cases.append(('bell_edge_%d_%d'%(i,j),tuple((np.eye(4)[i]*0.3+np.eye(4)[j]*0.7))))
    cases += [('uniform',(0.25,)*4),('asymmetric',(0.67,0.11,0.15,0.07)),
              ('asymmetric2',(0.17,0.43,0.09,0.31)),
              ('near_pure',(1-3e-10,1e-10,1e-10,1e-10))]
    for epsilon in [0,1e-6,0.001,0.01,0.05,0.1,0.3,0.4,0.5387,2/3,0.8,1]:
        cases.append(('epsilon_'+str(epsilon),isotropic_input(epsilon)))
    for concentration in [0.2,1.0,5.0]:
        for k in range(40):
            cases.append(('dirichlet_%.1f_%02d'%(concentration,k),
                          tuple(RNG.dirichlet(np.ones(4)*concentration))))
    return cases


def run():
    qs=[0,1e-8,1e-5,0.001,0.01,0.03,0.08,0.12,0.17,0.3,0.5,0.8,1]
    cases=test_families()
    metrics={'case_count':len(cases),'q_grid':qs,'tested_one_round_cases':0,
             'max_error_pretwirl_pop':0.0,'max_error_posttwirl_pop':0.0,
             'max_error_original_bbpssw_preposttwirl':0.0,'max_error_original_success_prob':0.0,
             'max_error_success_prob':0.0,'max_wernertwirl_offdiagonal':0.0,
             'min_output_eigenvalue':1.0,'max_trace_error':0.0,
             'max_bell_population_sum_error':0.0,
             'max_q0_original_bennett_fidelity_error':0.0,
             'max_error_closed_eps_formula':0.0,
             'repeat_case_count':0,'max_repeat_fidelity_error':0.0,
             'max_repeat_success_error':0.0,
             'max_clifford_group_unitarity_error':0.0,
             'noise_negative_controls':{},'frozen_circuit':{
              'wire_order':WIRE_ORDER,'cnot_control_target':BXOR_GATES,
              'noise_after_each_cnot':'two-qubit Pauli-average realization of replacement depolarization',
              'measurement_wires':MEASUREMENT_WIRES,
              'accepted_readouts':ACCEPTED_READOUTS,
              'twirl':'original BBPSSW: pre-twirl each input Bell pair and post-twirl retained Bell pair by 24 Clifford U tensor U*'},
             'seed':20261103}
    ref=metrics['frozen_circuit']
    metrics['frozen_circuit_sha256']=hashlib.sha256(json.dumps(ref,sort_keys=True).encode()).hexdigest()
    for V in clifford_twirl_ops():
        metrics['max_clifford_group_unitarity_error']=max(
            metrics['max_clifford_group_unitarity_error'],emax(V.conj().T@V,np.eye(4)))
    samples=[]
    for label,p in cases:
        for q in qs:
            pre,pa=round_untwirled(p,q)
            post,pb=round_twirl(p,q)
            ref_pre=independent_circuit_round(p,q,twirl=False)
            ref_post=independent_circuit_round(p,q,twirl=True)
            full_map,full_ps=round_original_bbpssw(p,q)
            ref_original=independent_circuit_round(p,q,twirl=True,pre_twirl=True)
            metrics['tested_one_round_cases']+=1
            metrics['max_error_pretwirl_pop']=max(metrics['max_error_pretwirl_pop'],emax(pre,ref_pre['populations']))
            metrics['max_error_posttwirl_pop']=max(metrics['max_error_posttwirl_pop'],emax(post,ref_post['populations']))
            metrics['max_error_original_bbpssw_preposttwirl']=max(metrics['max_error_original_bbpssw_preposttwirl'],emax(full_map,ref_original['populations']))
            metrics['max_error_original_success_prob']=max(metrics['max_error_original_success_prob'],abs(full_ps-ref_original['success']))
            metrics['max_error_success_prob']=max(metrics['max_error_success_prob'],abs(pa-ref_pre['success']),abs(pb-ref_post['success']))
            matrix=ref_post['final']
            metrics['min_output_eigenvalue']=min(metrics['min_output_eigenvalue'],float(np.linalg.eigvalsh(matrix).min()))
            metrics['max_trace_error']=max(metrics['max_trace_error'],abs(np.trace(matrix)-1))
            metrics['max_bell_population_sum_error']=max(metrics['max_bell_population_sum_error'],abs(sum(ref_post['populations'])-1))
            target=np.zeros((4,4),dtype=complex)
            for pj,P in zip(post,bell_projectors()):
                target += pj*P
            metrics['max_wernertwirl_offdiagonal']=max(metrics['max_wernertwirl_offdiagonal'],emax(target,matrix))
            if label.startswith('epsilon_'):
                eps=4*p[1]
                fe,ps=first_round_isotropic(eps,q)
                metrics['max_error_closed_eps_formula']=max(metrics['max_error_closed_eps_formula'],abs(fe-post[0]),abs(ps-pa))
                if q==0:
                    F=p[0]; e=(1-F)/3
                    original=(F*F+e*e)/(F*F+2*F*e+5*e*e)
                    metrics['max_q0_original_bennett_fidelity_error']=max(metrics['max_q0_original_bennett_fidelity_error'],abs(post[0]-original))
                samples.append({'epsilon':eps,'q':q,'f_in':p[0], 'f_out':post[0],
                                'p_success':pa,'max_p_diff':emax(post,ref_post['populations'])})
    # Compare repeated trajectories, feeding independent reference states back to circuit.
    for eps in [0,0.01,0.1,0.35,0.5387,2/3,0.8,1]:
        for q in [0,0.001,0.05,0.12,0.17,0.5,1]:
            p=isotropic_input(eps)
            for n in range(1,13):
                nxt,p_success=round_twirl(p,q)
                chk=independent_circuit_round(p,q,twirl=True)
                metrics['repeat_case_count']+=1
                metrics['max_repeat_fidelity_error']=max(metrics['max_repeat_fidelity_error'],emax(nxt,chk['populations']))
                metrics['max_repeat_success_error']=max(metrics['max_repeat_success_error'],abs(p_success-chk['success']))
                # Circuit reference result is used as the next-round input
                p=tuple(chk['populations'])
    p=(0.67,0.11,0.15,0.07);q=0.17
    expected,_=round_twirl(p,q)
    controls={
        'reverse_Bob_CNOT':independent_circuit_round(p,q,reverse_bob=True),
        'omit_Bob_depolarization':independent_circuit_round(p,q,bob_noise_override=0.0),
        'accept_unequal_measurements':independent_circuit_round(p,q,accepted_outcomes=((0,1),(1,0)))
    }
    for label,out in controls.items():
        metrics['noise_negative_controls'][label]={
            'population_difference':emax(out['populations'],expected),
            'success_probability':out['success']}
    assert len(clifford_twirl_ops())==24
    assert metrics['max_error_pretwirl_pop']<1e-12,metrics
    assert metrics['max_error_posttwirl_pop']<1e-12,metrics
    assert metrics['max_error_original_bbpssw_preposttwirl']<1e-12,metrics
    assert metrics['max_error_original_success_prob']<1e-12,metrics
    assert metrics['max_error_success_prob']<1e-12,metrics
    assert metrics['max_wernertwirl_offdiagonal']<1e-12,metrics
    assert metrics['min_output_eigenvalue']>-1e-12,metrics
    assert metrics['max_trace_error']<1e-12,metrics
    assert metrics['max_q0_original_bennett_fidelity_error']<1e-12,metrics
    assert metrics['max_error_closed_eps_formula']<1e-12,metrics
    assert metrics['max_repeat_fidelity_error']<1e-12,metrics
    assert metrics['max_repeat_success_error']<1e-12,metrics
    assert all(x['population_difference']>1e-3 for x in metrics['noise_negative_controls'].values()),metrics
    (OUTPUT/'verification_results.json').write_text(json.dumps(metrics,indent=2,ensure_ascii=True))
    with (OUTPUT/'isotropic_one_round_samples.csv').open('w',newline='') as fh:
        dw=csv.DictWriter(fh,fieldnames=list(samples[0]));dw.writeheader();dw.writerows(samples)
    print(json.dumps(metrics,indent=2))


if __name__=='__main__':
    run()
