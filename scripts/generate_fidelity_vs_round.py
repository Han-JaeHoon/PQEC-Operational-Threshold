#!/usr/bin/env python3
"""Bell fidelity vs PQEC round for Steps 3, 4 and 5 at ONE common condition.

Figure for Section 2 ("Main Results") of the summary document:

    x : PQEC round n        y : Bell fidelity F_n = <Phi+|rho_n|Phi+>
    same initial Bell-isotropic state rho_{t0} and the same per-CNOT noise q
    for all three circuits.

Condition
---------
    rho_0 = rho_{t0} = 1/4 [ II + t0 (XX - YY + ZZ) ],   t0 = 0.9  (eps_0 = 0.1)
    q     = 0.01   (two-qubit replacement depolarizing after EVERY CNOT)
    N     = 5000 rounds

Model (unchanged from the verified repository implementations)
-------------------------------------------------------------
The round map is `iterated_noisy_pqec.one_round_effective_map`:

    sigma_in = |0><0|_a (x) rho (x) rho
    tau_A    = Tr_{a,B}[ (Z_a (x) I) E_q(sigma_in) ]
    rho_{n+1}= tau_A / Tr(tau_A)

with the gate sequences captured from the verified Step-3/4/5 circuits and the
replacement depolarizing channel applied after every CNOT.  Nothing about the
noise model, the noise locations or q is altered by this script.

Two model choices, both taken from the Aug-2026 notes:

* Step 5 uses the EXACT-COMPILATION-CALIBRATED circuit (`step5cal`).  The learned
  circuit realises U_PQEC only to ||V - U||_F ~ 2e-7, which leaves a q-independent
  residual that survives even at q = 0 (where tau_A must equal rho^2 exactly).
  `iterated_noisy_pqec.compilation_correction` post-composes the round with
  C_corr = U Vbar^dag so the q = 0 map is exactly rho -> rho^2/Tr(rho^2).  The
  noise model is untouched.

* Steps 3 and 4 are iterated on the Bell-diagonal sector S_BD = span{II,XX,YY,ZZ},
  which their noisy round map leaves EXACTLY invariant (notes 01/02 Sec. 2; the
  transverse leakage measured per round here stays at ~1e-16).  In floating point
  that invariance is broken at the 1e-16 level, and because the Bell-directed fixed
  state is a full-state saddle (rho_sp(J_full) > 1) the roundoff residue is amplified
  and eventually produces an ARTIFICIAL escape -- Step 4 escapes near n ~ 1500 in an
  unprojected double-precision run.  Projecting each iterate back onto S_BD removes
  that artificial seed (note 04 Sec. 6).  The discarded transverse norm is recorded
  every round so the artifact stays visible.

Both raw (unprojected / uncalibrated) runs and the closed-form (u,v) recursions of
notes 01/02 are computed as controls and written to the CSV.

Outputs
-------
    figures/fidelity_vs_round.pdf                (+ .png)
    results/fidelity_vs_round/fidelity_vs_round.csv
    results/fidelity_vs_round/metadata.json

Run:  python scripts/generate_fidelity_vs_round.py
"""
import csv
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import iterated_noisy_pqec as ip

# ---------------------------------------------------------------------------
# condition
# ---------------------------------------------------------------------------
T0 = 0.9                       # Bell-isotropic input rho_{t0};  eps_0 = 1 - t0 = 0.1
Q = 0.01                       # per-CNOT replacement depolarizing strength
N_ROUNDS = 5000

FIG_DIR = os.path.join(ROOT, "figures")
RES_DIR = os.path.join(ROOT, "results", "fidelity_vs_round")

# Okabe-Ito blue / vermillion / bluish-green.  Validated with the dataviz
# palette checker (light surface, all-pairs): lightness band, chroma floor,
# CVD separation (worst dE 11.0 deutan), normal-vision floor (18.7) and
# contrast vs surface all PASS.  Line style is a second, colour-free encoding.
STYLE = {
    "step3": dict(color="#0072B2", ls="-",   label="Step 3 — textbook 16-CNOT"),
    "step4": dict(color="#D55E00", ls="--",  label="Step 4 — resynthesized 14-CNOT"),
    "step5": dict(color="#009E73", ls="-.",  label="Step 5 — learned 14-CNOT"),
}

# ---------------------------------------------------------------------------
# Bell-diagonal sector
# ---------------------------------------------------------------------------
_BD_LABELS = ("II", "XX", "YY", "ZZ")
_BD_OPS = [ip.PAULI_2Q[l] for l in _BD_LABELS]


def project_bell_diagonal(rho):
    """Orthogonal projection onto S_BD = span{II,XX,YY,ZZ}, and the discarded norm."""
    out = np.zeros_like(rho)
    for P in _BD_OPS:
        out = out + float(np.real(np.trace(P @ rho))) * P / 4.0
    return out, float(np.linalg.norm(rho - out))


def bell_fidelity(rho):
    """F = <Phi+|rho|Phi+> = Tr(Phi rho) -- the SAME definition for every circuit."""
    return float(np.real(np.trace(ip.PHI @ rho)))


# ---------------------------------------------------------------------------
# trajectories
# ---------------------------------------------------------------------------
def run_dense(circuit, project, n_rounds=N_ROUNDS, t0=T0, q=Q):
    """Iterate the full 32x32 five-qubit round map, recording one row per round."""
    rho = ip.rho_isotropic(t0)
    rows = []

    def record(n, Q_val, transverse, herm):
        eig = np.linalg.eigvalsh(0.5 * (rho + rho.conj().T))
        rows.append(dict(
            n=n, F=bell_fidelity(rho),
            purity=float(np.real(np.trace(rho @ rho))),
            Q=Q_val, transverse_norm=transverse, herm_err=herm,
            eig_min=float(eig.min()),
            trace_err=float(abs(np.trace(rho).real - 1.0)),
        ))

    _, tr0 = project_bell_diagonal(rho)
    record(0, float("nan"), tr0, 0.0)

    for n in range(1, n_rounds + 1):
        nxt, info = ip.one_round_effective_map(rho, q, circuit)
        if nxt is None:
            raise RuntimeError(f"{circuit}: denominator vanished at n = {n}")
        rho = nxt
        _, transverse = project_bell_diagonal(rho)
        if project:
            rho, _ = project_bell_diagonal(rho)
        record(n, info["Q"], transverse, info["herm_err"])
    return rows


def run_uv(circuit, n_rounds=N_ROUNDS, t0=T0, q=Q):
    """Closed-form (u,v) recursion on the invariant plane y = -x (notes 01/02).

        rho(u,v) = 1/4 [ II + u (XX - YY) + v ZZ ],   F = (1 + 2u + v)/4

    Step 3:  u' = 2 s^6 u (1+v) / D3,   v' = s^5 (1+s)(v + u^2) / D3,
             D3 = 1 + s^4 (2u^2 + v^2)
    Step 4:  u' = 2 s^4 u (1+v) / D4,   v' = s^3 (1+s)(v + u^2) / D4,
             D4 = 1 + s^2 (2u^2 + v^2)
    """
    s = 1.0 - q
    if circuit == "step3":
        a, b, c = s ** 6, s ** 5 * (1 + s), s ** 4
    elif circuit == "step4":
        a, b, c = s ** 4, s ** 3 * (1 + s), s ** 2
    else:
        raise KeyError(circuit)
    u, v = t0, t0
    rows = [dict(n=0, F=(1 + 2 * u + v) / 4.0, u=u, v=v)]
    for n in range(1, n_rounds + 1):
        D = 1.0 + c * (2 * u * u + v * v)
        u, v = 2 * a * u * (1 + v) / D, b * (v + u * u) / D
        rows.append(dict(n=n, F=(1 + 2 * u + v) / 4.0, u=u, v=v))
    return rows


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
def make_figure(series):
    plt.rcParams.update({
        "font.size": 9, "axes.labelsize": 9.5, "axes.titlesize": 9.5,
        "legend.fontsize": 8.5, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "mathtext.fontset": "dejavuserif", "font.family": "serif",
    })
    fig, (axL, axR) = plt.subplots(
        1, 2, figsize=(9.6, 3.9), dpi=300,
        gridspec_kw=dict(width_ratios=[1.0, 1.2], wspace=0.30))

    for ax in (axL, axR):
        ax.grid(True, ls=":", lw=0.45, color="0.78")
        ax.set_axisbelow(True)
        ax.set_xlabel(r"PQEC round  $n$")
        ax.set_ylabel(r"Bell fidelity  $F_n=\langle\Phi^+|\rho_n|\Phi^+\rangle$")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    ZOOM = (0.920, 0.9925)          # y-window of panel (a), marked on panel (b)

    # ---- (a) early rounds: the convergence ordering F*(5) > F*(4) > F*(3) ----
    for key in ("step3", "step4", "step5"):
        n, F = series[key]
        m = n <= 10
        axL.plot(n[m], F[m], lw=1.7, ms=4.2, marker="o",
                 markeredgecolor="white", markeredgewidth=0.6, **STYLE[key])
    axL.set_xlim(-0.4, 15.6)
    axL.set_ylim(*ZOOM)
    axL.set_xticks([0, 2, 4, 6, 8, 10])
    axL.set_title("(a)  early rounds — Bell-sector convergence",
                  loc="left", fontweight="bold")
    for key in ("step3", "step4", "step5"):
        n, F = series[key]
        axL.annotate(f"{STYLE[key]['label'].split(' — ')[0]}   {F[10]:.4f}",
                     xy=(10.55, F[10]), color=STYLE[key]["color"],
                     fontsize=8.2, fontweight="bold", ha="left", va="center")

    # ---- (b) full range on a log axis: the Step-5 plateau and its escape ----
    axR.axhspan(*ZOOM, color="0.90", zorder=0)
    axR.annotate("range of panel (a)", xy=(4300, ZOOM[0] + 0.005), fontsize=7.5,
                 color="0.45", ha="right", va="bottom")
    for key in ("step3", "step4", "step5"):
        n, F = series[key]
        m = n >= 1
        axR.plot(n[m], F[m], lw=1.8, **STYLE[key])
    axR.set_xscale("log")
    axR.set_xlim(1, N_ROUNDS)
    axR.set_ylim(0.38, 1.075)
    axR.set_title("(b)  long-time behaviour", loc="left", fontweight="bold")
    axR.legend(frameon=False, loc="center left", bbox_to_anchor=(0.015, 0.34),
               handlelength=2.6, labelspacing=0.45)

    axR.annotate("Steps 3 / 4 — flat through  $n=5000$", xy=(1.8, 1.030),
                 fontsize=8.4, color="0.2", ha="left", va="center")
    axR.annotate("Step 5 leaves\nthe plateau", xy=(690, 0.905), xytext=(33, 0.70),
                 fontsize=8.4, color="0.2", ha="center",
                 arrowprops=dict(arrowstyle="->", lw=0.9, color="0.4",
                                 connectionstyle="arc3,rad=-0.2"))
    axR.annotate("separable product state\n" r"$F_\infty\simeq0.4080$",
                 xy=(3600, 0.4085), xytext=(23, 0.437), fontsize=8.4, color="0.2",
                 ha="left", va="center",
                 arrowprops=dict(arrowstyle="->", lw=0.9, color="0.4",
                                 connectionstyle="arc3,rad=-0.12"))

    fig.suptitle(
        r"Repeated noisy PQEC:  $\rho_0=\rho_{t_0}$ with $t_0=%.1f$,  "
        r"per-CNOT depolarizing $q=%.2f$ — identical for all three circuits" % (T0, Q),
        fontsize=9.5, y=1.02)

    for ext in ("pdf", "png"):
        path = os.path.join(FIG_DIR, f"fidelity_vs_round.{ext}")
        fig.savefig(path, bbox_inches="tight", facecolor="white")
        print(f"  wrote {os.path.relpath(path, ROOT)}")
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(RES_DIR, exist_ok=True)
    t_start = time.time()

    print("=" * 78)
    print(f" Bell fidelity vs PQEC round   t0 = {T0}, q = {Q}, N = {N_ROUNDS}")
    print("=" * 78)

    print("\n[validation] one-round map vs the repository executors / closed forms")
    ip.validate_one_round(verbose=True)

    print("\n[validation] phase-aligned compilation residual  ||Vbar - U_PQEC||")
    residuals = {}
    for c in ("step3", "step4", "step5"):
        mx, fr = ip.compilation_residual(c)
        residuals[c] = dict(max_abs=mx, frobenius=fr)
        print(f"    {c:9s} max|dV| = {mx:.3e}   ||dV||_F = {fr:.3e}")

    print("\n[validation] ideal limit  max||tau_A - rho^2||  at q = 0  (5 random states)")
    rng = np.random.default_rng(0)
    states = []
    for _ in range(5):
        A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        r = A @ A.conj().T
        states.append(r / np.trace(r).real)
    ideal_err = {}
    for c in ("step3", "step4", "step5", "step5cal"):
        e = max(float(np.linalg.norm(ip.one_round_tau(r, 0.0, c) - r @ r)) for r in states)
        ideal_err[c] = e
        print(f"    {c:9s} {e:.3e}")

    # ---- trajectories -----------------------------------------------------
    print(f"\n[run] {N_ROUNDS} rounds per series")
    runs = {}
    for name, circuit, project in (
            ("step3", "step3", True),
            ("step4", "step4", True),
            ("step5", "step5cal", False),
            ("step3_raw", "step3", False),
            ("step4_raw", "step4", False),
            ("step5_uncalibrated", "step5", False)):
        t0 = time.time()
        runs[name] = run_dense(circuit, project)
        print(f"    {name:20s} circuit={circuit:9s} project_S_BD={str(project):5s} "
              f"({time.time() - t0:.0f}s)")
    for name in ("step3", "step4"):
        runs[name + "_uv_closed_form"] = run_uv(name)
    print("    step3_uv_closed_form / step4_uv_closed_form  (notes 01/02 recursion)")

    # ---- consistency checks ----------------------------------------------
    print("\n[check] Bell-diagonal leakage per round  max_n ||Pi_perp rho_n||")
    leak = {}
    for name in ("step3", "step4", "step5", "step3_raw", "step4_raw"):
        leak[name] = max(r["transverse_norm"] for r in runs[name])
        print(f"    {name:12s} {leak[name]:.3e}")

    print("\n[check] projected dense run vs the closed-form (u,v) recursion  max_n |dF|")
    uv_err = {}
    for name in ("step3", "step4"):
        d = max(abs(a["F"] - b["F"])
                for a, b in zip(runs[name], runs[name + "_uv_closed_form"]))
        uv_err[name] = d
        print(f"    {name:12s} {d:.3e}")

    print("\n[check] physicality over the plotted runs")
    phys = {}
    for name in ("step3", "step4", "step5"):
        rs = runs[name]
        phys[name] = dict(min_eig=min(r["eig_min"] for r in rs),
                          max_trace_err=max(r["trace_err"] for r in rs),
                          max_herm_err=max(r["herm_err"] for r in rs),
                          min_Q=min(r["Q"] for r in rs[1:]))
        print(f"    {name:8s} min eig = {phys[name]['min_eig']:+.2e}   "
              f"max |Tr-1| = {phys[name]['max_trace_err']:.1e}   "
              f"max herm err = {phys[name]['max_herm_err']:.1e}   "
              f"min Q = {phys[name]['min_Q']:.3f}")

    # ---- the numbers quoted in the summary document -----------------------
    marks = [0, 1, 2, 3, 5, 10, 100, 500, 1000, 2000, 5000]
    print("\n[table] F_n at the quoted rounds")
    print("    " + "n".rjust(6) + "".join(f"{k:>14s}" for k in ("step3", "step4", "step5")))
    F = {k: {r["n"]: r["F"] for r in runs[k]} for k in ("step3", "step4", "step5")}
    for n in marks:
        print(f"    {n:>6d}" + "".join(f"{F[k][n]:>14.6f}" for k in ("step3", "step4", "step5")))

    # ---- write data -------------------------------------------------------
    csv_path = os.path.join(RES_DIR, "fidelity_vs_round.csv")
    cols = ["series", "n", "F", "purity", "Q", "transverse_norm", "herm_err",
            "eig_min", "trace_err", "u", "v"]
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for name, rows in runs.items():
            for r in rows:
                w.writerow(dict(series=name, **r))
    print(f"\n[data] wrote {os.path.relpath(csv_path, ROOT)} "
          f"({sum(len(v) for v in runs.values())} rows, {len(runs)} series)")

    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                         text=True).strip()
    except Exception:
        commit = "unknown"
    meta = dict(
        script="scripts/generate_fidelity_vs_round.py",
        git_commit=commit,
        generated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        condition=dict(t0=T0, eps0=1 - T0, q=Q, n_rounds=N_ROUNDS,
                       rho0="1/4[II + t0(XX - YY + ZZ)]",
                       fidelity="F_n = <Phi+|rho_n|Phi+> = Tr(Phi rho_n)"),
        model=dict(
            round_map="iterated_noisy_pqec.one_round_effective_map",
            noise="two-qubit replacement depolarizing D_q after every CNOT; "
                  "single-qubit gates ideal",
            step3="iterated_noisy_pqec circuit 'step3' (verify_analytic_decomposed._fred), "
                  "dense 32x32 map, projected onto S_BD every round",
            step4="iterated_noisy_pqec circuit 'step4' (pqec_resynth_noise.GATES), "
                  "dense 32x32 map, projected onto S_BD every round",
            step5="iterated_noisy_pqec circuit 'step5cal' (pqc_ring_prune.ansatz_masked + "
                  "pqc_ring_pruned_params.npy, post-composed with the exact-compilation "
                  "correction C_corr = U_PQEC Vbar^dag), dense 32x32 map, no projection",
            controls="step3_raw / step4_raw (no S_BD projection), step5_uncalibrated "
                     "(no C_corr), step3_uv_closed_form / step4_uv_closed_form "
                     "((u,v) recursion of notes 01/02)"),
        n_cnots={c: ip.n_cnots(c) for c in ("step3", "step4", "step5", "step5cal")},
        validation=dict(compilation_residual=residuals, ideal_limit_q0=ideal_err,
                        max_bell_diagonal_leakage=leak,
                        dense_vs_closed_form_max_dF=uv_err, physicality=phys),
        F_at_rounds={k: {str(n): F[k][n] for n in marks} for k in F},
        palette={k: STYLE[k]["color"] for k in STYLE},
        outputs=["figures/fidelity_vs_round.pdf", "figures/fidelity_vs_round.png",
                 "results/fidelity_vs_round/fidelity_vs_round.csv",
                 "results/fidelity_vs_round/metadata.json"],
    )
    meta_path = os.path.join(RES_DIR, "metadata.json")
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"[data] wrote {os.path.relpath(meta_path, ROOT)}")

    # ---- figure -----------------------------------------------------------
    print("\n[figure]")
    series = {k: (np.array([r["n"] for r in runs[k]]),
                  np.array([r["F"] for r in runs[k]])) for k in ("step3", "step4", "step5")}
    make_figure(series)

    print(f"\ndone in {time.time() - t_start:.0f}s")


if __name__ == "__main__":
    main()
