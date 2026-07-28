"""Out-of-sample validation of theta_free: unseen observables, dense eps grid,
and monotonicity of F(eps2) (assumption behind the threshold bisection).

Requires the trained parameters: run  python pqc_noise_aware.py  first (it writes
pqc_noise_aware.npz, which is git-ignored as a generated artifact)."""
import numpy as np
from pqc_common import ansatz_ops, make_F_qnode
from noisy_bell_state import rho_eps_analytic, O_PHI_PLUS
Z=np.array([[1,0],[0,-1]],complex); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]],complex)
ops,npar=ansatz_ops(6)
th=np.load("pqc_noise_aware.npz")["th_free"]
fqn=make_F_qnode(ops)
def Fv(eps,eps2,O):
    zO,zI=fqn(th,eps,eps2,O); return float(zO)/float(zI)
def Fex(eps,O):
    r=rho_eps_analytic(eps); r2=r@r; return float(np.real(np.trace(O@r2)/np.trace(r2)))
OBS={'Phi+':O_PHI_PLUS,'ZZ':np.kron(Z,Z),'XX':np.kron(X,X),'YY':np.kron(Y,Y),
     'ZI':np.kron(Z,np.eye(2)),'IZ':np.kron(np.eye(2),Z)}
TRAINED={'Phi+','ZZ'}
eps_dense=np.linspace(0.02,0.72,40)
print("noise-free out-of-sample (dense eps in [0.02,0.72], 40 pts):")
for nm,O in OBS.items():
    err=max(abs(Fv(e,0.0,O)-Fex(e,O)) for e in eps_dense)
    tag="(trained)" if nm in TRAINED else "(UNSEEN) "
    print(f"  O={nm:4s} {tag}: max|F-Fexact| = {err:.4f}",flush=True)
print("monotonicity of F(eps2) at fixed eps (max upward step; <=0 => monotone decreasing):")
for eps in [0.2,0.4,0.6]:
    q=np.linspace(0,0.5,60); Fs=[Fv(eps,x,O_PHI_PLUS) for x in q]
    up=max(Fs[i+1]-Fs[i] for i in range(len(Fs)-1))
    print(f"  eps={eps}: max upward step = {up:+.5f}",flush=True)
