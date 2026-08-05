"""
Sweep: iterated noisy PQEC for Step 3 / Step 4 / Step 5.
========================================================

For every (circuit, q, initial epsilon) this
  * iterates rho_{n+1} = P_q(rho_n) and stores the full per-iteration diagnostics,
  * solves the fixed-point equation rho = P_q(rho) directly (root-finding, which is
    insensitive to the saddle instability of the iteration),
  * computes the local Jacobian spectral radius and the structural coherence
    amplification factor (lambda_1+lambda_2)/Tr(rho^2),
and writes machine-readable results + the analysis figures.

Run:  python scripts/run_iterated_pqec.py
"""
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import iterated_noisy_pqec as ip

RES = os.path.join(ROOT, "results", "iterated_pqec")
FIGS = os.path.join(ROOT, "figures", "iterated_pqec")
os.makedirs(RES, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)

Q_GRID = [0.0, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 2e-2, 5e-2]
EPS_GRID = [0.05, 0.10, 0.20, 0.30, 0.50]
TOL = 1e-12
MAX_ITER = 1000
CIRCUITS = ["step3", "step4", "step5"]
LABEL = {"step3": "Step 3: textbook 16-CNOT",
         "step4": "Step 4: resynthesized 14-CNOT",
         "step5": "Step 5: learned 14-CNOT"}
COLOR = {"step3": "C2", "step4": "C0", "step5": "C3"}
MARK = {"step3": "^", "step4": "o", "step5": "s"}

TRAJ_COLS = ["circuit", "q", "eps", "n", "F", "purity", "Q", "diff", "diff2",
             "eig0", "eig1", "eig2", "eig3", "eig_min",
             "p_Phi_plus", "p_Phi_minus", "p_Psi_plus", "p_Psi_minus", "bell_offdiag",
             "delta_iso", "iso_XXplusYY", "iso_XXminusZZ",
             "herm_err", "trace_err"] + [f"P_{l}" for l in ip.PAULI_LABELS]


def _git_hash():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                       text=True).strip()
    except Exception:
        return "unknown"


def run_sweep():
    traj_rows = {c: [] for c in CIRCUITS}
    fixed = []
    rho_star_store = {}
    t0 = time.time()

    for c in CIRCUITS:
        print(f"\n=== {LABEL[c]}  ({ip.n_cnots(c)} CNOTs) ===", flush=True)
        for q in Q_GRID:
            finals = []
            for eps in EPS_GRID:
                rho0 = ip.rho_isotropic(1 - eps)
                out = ip.iterate_effective_map(rho0, q, c, tol=TOL, max_iter=MAX_ITER,
                                               record_every=10)
                for rec in out["history"]:
                    row = dict(circuit=c, q=q, eps=eps, n=rec["n"], F=rec["F"],
                               purity=rec["purity"], Q=rec.get("Q", np.nan),
                               diff=rec.get("diff", np.nan), diff2=rec.get("diff2", np.nan),
                               eig_min=rec["eig_min"],
                               p_Phi_plus=rec["p_Phi_plus"], p_Phi_minus=rec["p_Phi_minus"],
                               p_Psi_plus=rec["p_Psi_plus"], p_Psi_minus=rec["p_Psi_minus"],
                               bell_offdiag=rec["bell_offdiag"], delta_iso=rec["delta_iso"],
                               iso_XXplusYY=rec["iso_XXplusYY"],
                               iso_XXminusZZ=rec["iso_XXminusZZ"],
                               herm_err=rec["herm_err"], trace_err=rec["trace_err"])
                    for k, e in enumerate(rec["eigs"]):
                        row[f"eig{k}"] = e
                    for l in ip.PAULI_LABELS:
                        row[f"P_{l}"] = rec[f"P_{l}"]
                    traj_rows[c].append(row)
                finals.append((eps, out))

            # --- initial-state dependence: spread of the iterated end states ---
            rhos = [o["rho"] for _, o in finals]
            spread = max(float(np.linalg.norm(rhos[i] - rhos[j]))
                         for i in range(len(rhos)) for j in range(i + 1, len(rhos)))

            # --- direct fixed-point solve + stability ---
            rs, resid, ok, info = ip.solve_fixed_point(q, c)
            J, sr, _ = ip.jacobian(rs, q, c)
            coh = ip.coherence_amplification(rs)
            rec = ip.state_record(rs, Q=info.get("Q", np.nan))
            rho_star_store[f"{c}_q{q:.8g}"] = rs

            statuses = [o["status"] for _, o in finals]
            fixed.append(dict(
                circuit=c, q=q,
                F_star=rec["F"], purity_star=rec["purity"], Q_star=rec.get("Q", np.nan),
                eig_min_star=rec["eig_min"], residual=resid, solver_ok=ok,
                spectral_radius=sr, coherence_factor=coh,
                delta_iso_star=rec["delta_iso"], bell_offdiag_star=rec["bell_offdiag"],
                XX=rec["XX"], YY=rec["YY"], ZZ=rec["ZZ"],
                iter_spread_over_eps=spread,
                n_converged=sum(s == "converged" for s in statuses),
                n_runs=len(statuses),
                iter_F_mean=float(np.mean([o["final"]["F"] for _, o in finals])),
                iter_F_min=float(np.min([o["final"]["F"] for _, o in finals])),
                iter_F_max=float(np.max([o["final"]["F"] for _, o in finals])),
                max_herm_err=max(o["max_herm_err"] for _, o in finals),
                min_eig_seen=min(o["min_eig_seen"] for _, o in finals
                                 if o["min_eig_seen"] is not None),
                min_Q_seen=min(o["min_Q_seen"] for _, o in finals
                               if o["min_Q_seen"] is not None),
                cycle2=float(np.mean([o["cycle2"] for _, o in finals])),
                mean_n_iter=float(np.mean([o["n_iter"] for _, o in finals])),
            ))
            print(f"  q={q:<8g} F*={rec['F']:.10f}  P*={rec['purity']:.8f}  "
                  f"Q*={rec.get('Q', np.nan):.4e}  resid={resid:.1e}  "
                  f"rho(J)={sr:.5f}  coh={coh:.5f}  conv={fixed[-1]['n_converged']}/"
                  f"{len(statuses)}", flush=True)

    print(f"\nsweep done in {time.time()-t0:.0f}s")
    return traj_rows, fixed, rho_star_store


def write_results(traj_rows, fixed, rho_star_store):
    import csv
    for c in CIRCUITS:
        path = os.path.join(RES, f"{c}_results.csv")
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=TRAJ_COLS, extrasaction="ignore")
            w.writeheader()
            for r in traj_rows[c]:
                w.writerow(r)
        print(f"  wrote {os.path.relpath(path, ROOT)}  ({len(traj_rows[c])} rows)")

    fpath = os.path.join(RES, "fixed_points.csv")
    keys = list(fixed[0].keys())
    with open(fpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(fixed)
    print(f"  wrote {os.path.relpath(fpath, ROOT)}  ({len(fixed)} rows)")

    np.savez(os.path.join(RES, "fixed_points.npz"), **rho_star_store)
    print(f"  wrote {os.path.relpath(os.path.join(RES,'fixed_points.npz'), ROOT)}")

    meta = dict(
        git_commit=_git_hash(),
        branch="iterated-noisy-pqec",
        circuits={c: dict(identifier=c, label=LABEL[c], n_cnot=ip.n_cnots(c),
                          source={"step3": "verify_analytic_decomposed._fred (textbook "
                                           "8-CNOT/Controlled-SWAP decomposition)",
                                  "step4": "pqec_resynth_noise.GATES (Qiskit transpile, "
                                           "basis {u,cx}, opt_level 3, seed 7)",
                                  "step5": "pqc_ring_prune.ansatz_masked + "
                                           "pqc_ring_pruned_params.npy"}[c])
                  for c in CIRCUITS},
        q_grid=Q_GRID, initial_epsilon_grid=EPS_GRID,
        convergence_tolerance=TOL, max_iterations=MAX_ITER,
        convergence_criterion="||rho_{n+1}-rho_n||_F < tol",
        noise_channel="D_q^(ij)(rho) = (1-q) rho + q [ I_ij/4 (x) Tr_ij(rho) ]  "
                      "(two-qubit replacement depolarizing, applied after every CNOT; "
                      "all single-qubit gates ideal)",
        qubit_ordering="wire 0 = ancilla a; wires 1,2 = retained register A; "
                       "wires 3,4 = discarded register B",
        effective_map="tau_A = Tr_{a,B}[(Z_a (x) I) sigma_out];  "
                      "rho_{n+1} = tau_A / Tr(tau_A)   (Z_a taken AFTER the final H, "
                      "matching the repository read-out convention)",
        input_state="rho_t = 1/4[II + t(XX - YY + ZZ)], t = 1 - epsilon",
        hermitian_projection="each iterate is projected onto the Hermitian manifold "
                             "(exact-arithmetic identity); the discarded anti-Hermitian "
                             "norm is recorded as herm_err and stays ~1e-16",
        timestamp=datetime.now(timezone.utc).isoformat(),
        backend=f"numpy {np.__version__} dense 32x32 density matrix, dtype=complex128",
        python=platform.python_version(),
        validation="one-round Q, N_Phi, F reproduce the repository PennyLane executors "
                   "and the Step-3/Step-4 closed forms to <= 1.5e-15",
    )
    mpath = os.path.join(RES, "metadata.json")
    json.dump(meta, open(mpath, "w"), indent=2)
    print(f"  wrote {os.path.relpath(mpath, ROOT)}")


# ===========================================================================
# figures
# ===========================================================================
def make_figures(traj_rows, fixed):
    fx = {c: [f for f in fixed if f["circuit"] == c] for c in CIRCUITS}
    NOISE = (r"noise: $D_q$ after every CNOT, $D_q(\rho)=(1-q)\rho+q\,"
             r"[I_{ij}/4\otimes\mathrm{Tr}_{ij}\rho]$")

    # 1. F_n vs n
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), dpi=200, sharey=True)
    for ax, c in zip(axes, CIRCUITS):
        for q, ls in zip([1e-3, 1e-2, 5e-2], ["-", "--", ":"]):
            for eps, col in zip([0.05, 0.20, 0.50], ["C0", "C1", "C3"]):
                rows = [r for r in traj_rows[c] if r["q"] == q and r["eps"] == eps]
                rows.sort(key=lambda r: r["n"])
                ax.plot([r["n"] for r in rows], [r["F"] for r in rows], ls, color=col,
                        lw=1.3, label=f"q={q:g}, ε₀={eps}" if c == "step3" else None)
        ax.set_xscale("symlog", linthresh=1)
        ax.set_xlabel("iteration  $n$")
        ax.set_title(LABEL[c], fontsize=10)
        ax.grid(True, ls=":", lw=.5)
    axes[0].set_ylabel(r"Bell fidelity  $F_n=\mathrm{Tr}(\Phi\rho_n)$")
    axes[0].legend(fontsize=7, frameon=False, ncol=1)
    fig.suptitle("Iterated noisy PQEC: fidelity vs iteration      " + NOISE, fontsize=10)
    fig.savefig(os.path.join(FIGS, "F_vs_iteration.png"), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    # 2. F*(q)
    fig, ax = plt.subplots(figsize=(6.6, 4.6), dpi=200)
    for c in CIRCUITS:
        ax.plot([f["q"] for f in fx[c]], [f["F_star"] for f in fx[c]],
                "-" + MARK[c], color=COLOR[c], lw=1.8, ms=6, label=LABEL[c])
    ax.set_xscale("symlog", linthresh=1e-5)
    ax.set_xlabel(r"per-CNOT noise  $q$"); ax.set_ylabel(r"$F_*(q)=\mathrm{Tr}(\Phi\rho_*)$")
    ax.set_title("Fixed-point Bell fidelity\n" + NOISE, fontsize=10)
    ax.grid(True, ls=":", lw=.6); ax.legend(frameon=False, fontsize=9)
    fig.savefig(os.path.join(FIGS, "Fstar_vs_q.png"), bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # 3. 1-F*(q) log-log + fit
    fig, ax = plt.subplots(figsize=(6.6, 4.6), dpi=200)
    fits = {}
    for c in CIRCUITS:
        qs = np.array([f["q"] for f in fx[c]]); Fs = np.array([f["F_star"] for f in fx[c]])
        m = (qs > 0) & (1 - Fs > 1e-14)
        ax.loglog(qs[m], 1 - Fs[m], MARK[c], color=COLOR[c], ms=6, label=LABEL[c])
        sel = m & (qs <= 1e-3)
        if sel.sum() >= 3:
            a, b = np.polyfit(np.log(qs[sel]), np.log(1 - Fs[sel]), 1)
            fits[c] = (a, float(np.exp(b)))
            ax.loglog(qs[m], np.exp(b) * qs[m] ** a, "-", color=COLOR[c], lw=1.1, alpha=.7)
    ax.set_xlabel(r"per-CNOT noise  $q$"); ax.set_ylabel(r"$1-F_*(q)$")
    ttl = "  ".join(f"{c}: α={fits[c][0]:.3f}, A={fits[c][1]:.3g}" for c in fits)
    ax.set_title("Small-$q$ scaling of the fixed-point infidelity\n" + ttl, fontsize=9)
    ax.grid(True, which="both", ls=":", lw=.5); ax.legend(frameon=False, fontsize=9)
    fig.savefig(os.path.join(FIGS, "one_minus_Fstar_loglog.png"), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    # 4. purity
    fig, ax = plt.subplots(figsize=(6.6, 4.6), dpi=200)
    for c in CIRCUITS:
        ax.plot([f["q"] for f in fx[c]], [f["purity_star"] for f in fx[c]],
                "-" + MARK[c], color=COLOR[c], lw=1.8, ms=6, label=LABEL[c])
    ax.set_xscale("symlog", linthresh=1e-5)
    ax.set_xlabel(r"per-CNOT noise  $q$"); ax.set_ylabel(r"$\mathrm{Tr}(\rho_*^2)$")
    ax.set_title("Fixed-point purity\n" + NOISE, fontsize=10)
    ax.grid(True, ls=":", lw=.6); ax.legend(frameon=False, fontsize=9)
    fig.savefig(os.path.join(FIGS, "purity_vs_q.png"), bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # 5. Pauli coefficients
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), dpi=200, sharex=True)
    for ax, key in zip(axes, ["XX", "YY", "ZZ"]):
        for c in CIRCUITS:
            ax.plot([f["q"] for f in fx[c]], [f[key] for f in fx[c]],
                    "-" + MARK[c], color=COLOR[c], lw=1.6, ms=5, label=LABEL[c])
        ax.set_xscale("symlog", linthresh=1e-5); ax.set_xlabel(r"$q$")
        ax.set_ylabel(rf"$\langle {key}\rangle_*$"); ax.grid(True, ls=":", lw=.5)
        ax.set_title(f"fixed-point {key}", fontsize=10)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Fixed-point Pauli correlators      " + NOISE, fontsize=10)
    fig.savefig(os.path.join(FIGS, "pauli_vs_q.png"), bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # 6. Werner deviation
    fig, ax = plt.subplots(figsize=(6.6, 4.6), dpi=200)
    for c in CIRCUITS:
        ax.plot([f["q"] for f in fx[c]], [max(f["delta_iso_star"], 1e-18) for f in fx[c]],
                "-" + MARK[c], color=COLOR[c], lw=1.8, ms=6, label=LABEL[c])
    ax.set_xscale("symlog", linthresh=1e-5); ax.set_yscale("log")
    ax.set_xlabel(r"per-CNOT noise  $q$")
    ax.set_ylabel(r"$\|\rho_*-\rho_W(F_*)\|_F$")
    ax.set_title("Deviation of the fixed point from the isotropic (Werner) family\n" + NOISE,
                 fontsize=10)
    ax.grid(True, which="both", ls=":", lw=.5); ax.legend(frameon=False, fontsize=9)
    fig.savefig(os.path.join(FIGS, "werner_deviation_vs_q.png"), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    # 7. Q*
    fig, ax = plt.subplots(figsize=(6.6, 4.6), dpi=200)
    for c in CIRCUITS:
        ax.plot([f["q"] for f in fx[c]], [f["Q_star"] for f in fx[c]],
                "-" + MARK[c], color=COLOR[c], lw=1.8, ms=6, label=LABEL[c])
    ax.set_xscale("symlog", linthresh=1e-5)
    ax.set_xlabel(r"per-CNOT noise  $q$")
    ax.set_ylabel(r"$Q_*=\mathrm{Tr}(\tau_A)=\langle Z_a\rangle$")
    ax.set_title("Parity visibility at the fixed point (sampling overhead $\\propto Q_*^{-2}$)\n"
                 + NOISE, fontsize=10)
    ax.grid(True, ls=":", lw=.6); ax.legend(frameon=False, fontsize=9)
    fig.savefig(os.path.join(FIGS, "Qstar_vs_q.png"), bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # 8. stability (extra, central to the interpretation)
    fig, ax = plt.subplots(figsize=(6.6, 4.6), dpi=200)
    for c in CIRCUITS:
        ax.plot([f["q"] for f in fx[c]], [f["spectral_radius"] for f in fx[c]],
                "-" + MARK[c], color=COLOR[c], lw=1.8, ms=6, label=LABEL[c] + r"  $\rho(J)$")
        ax.plot([f["q"] for f in fx[c]], [f["coherence_factor"] for f in fx[c]],
                ":", color=COLOR[c], lw=1.2,
                label=LABEL[c] + r"  $(\lambda_1+\lambda_2)/\mathrm{Tr}\rho^2$")
    ax.axhline(1.0, color="k", lw=1, ls="--")
    ax.set_xscale("symlog", linthresh=1e-5)
    ax.set_xlabel(r"per-CNOT noise  $q$"); ax.set_ylabel("linear growth factor")
    ax.set_title("Local stability of the fixed point (>1 ⇒ saddle)", fontsize=10)
    ax.grid(True, ls=":", lw=.6); ax.legend(frameon=False, fontsize=7)
    fig.savefig(os.path.join(FIGS, "stability_vs_q.png"), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    print(f"  wrote 8 figures -> {os.path.relpath(FIGS, ROOT)}")
    return fits


def main():
    print("=" * 78)
    print(" Iterated noisy PQEC sweep")
    print("=" * 78)
    print("\n[validation] one-round map vs repository executors / closed forms")
    ip.selftest()
    ip.validate_one_round()

    traj, fixed, store = run_sweep()
    print("\n[results]")
    write_results(traj, fixed, store)
    print("\n[figures]")
    fits = make_figures(traj, fixed)
    print("\n[small-q fits]  1-F*(q) ~ A q^alpha  (q <= 1e-3)")
    for c, (a, A) in fits.items():
        print(f"    {c}: alpha = {a:.4f},  A = {A:.6g}")


if __name__ == "__main__":
    main()
