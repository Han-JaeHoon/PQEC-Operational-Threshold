#!/usr/bin/env python3
"""Tests for the Bell-fidelity-vs-round figure (scripts/generate_fidelity_vs_round.py).

Covers the two model choices the figure relies on:
  * the exact-compilation calibration of the learned Step-5 circuit, and
  * the exact Bell-diagonal invariance that licenses the S_BD projection in Steps 3/4,
cross-checked against the closed-form (u,v) recursion of the Aug-2026 notes.

Run:  python tests/test_fidelity_vs_round.py
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import iterated_noisy_pqec as ip
import generate_fidelity_vs_round as g

N_SHORT = 200            # enough for both maps to sit on their fixed states


def test_calibration_makes_the_q0_map_exactly_rho_squared():
    """step5 floors at ~1e-8 (its learned-unitary residual); step5cal removes it."""
    rng = np.random.default_rng(1)
    raw = cal = 0.0
    for _ in range(5):
        A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        r = A @ A.conj().T
        r /= np.trace(r).real
        raw = max(raw, np.linalg.norm(ip.one_round_tau(r, 0.0, "step5") - r @ r))
        cal = max(cal, np.linalg.norm(ip.one_round_tau(r, 0.0, "step5cal") - r @ r))
    assert raw > 1e-10, raw
    assert cal < 1e-13, cal


def test_compilation_residual_matches_the_note():
    """||Vbar - U_PQEC||: 2.34e-8 (max) / 2.04e-7 (Frobenius) in note 03 Sec. 1."""
    mx, fr = ip.compilation_residual("step5")
    assert abs(mx - 2.34e-8) < 5e-10, mx
    assert abs(fr - 2.04e-7) < 5e-9, fr
    for c in ("step3", "step4"):
        assert ip.compilation_residual(c)[1] < 1e-13, c


def test_calibration_leaves_the_noise_model_alone():
    """Same CNOT count, and the one-round read-out moves only by the O(1e-7) residual."""
    assert ip.n_cnots("step5cal") == ip.n_cnots("step5") == 14
    rho = ip.rho_isotropic(0.9)
    for q in (0.0, 0.01, 0.05):
        d = np.linalg.norm(ip.one_round_tau(rho, q, "step5cal")
                           - ip.one_round_tau(rho, q, "step5"))
        assert d < 1e-8, (q, d)


def test_step3_step4_keep_the_bell_diagonal_sector():
    """One-round leakage out of S_BD from a Bell-diagonal input is at roundoff.

    `transverse_norm` is recorded BEFORE the projection, so on a projected run it is
    exactly the per-round leakage from a Bell-diagonal state -- the quantity that
    licenses the projection.
    """
    for c in ("step3", "step4"):
        rows = g.run_dense(c, project=True, n_rounds=N_SHORT)
        assert max(r["transverse_norm"] for r in rows) < 1e-14, c
    # Step 5, by contrast, leaks at O(q) -- S_BD is genuinely not invariant there.
    rows5 = g.run_dense("step5cal", project=False, n_rounds=50)
    assert max(r["transverse_norm"] for r in rows5) > 1e-4


def test_unprojected_step4_amplifies_its_roundoff_seed():
    """The artifact the projection removes: without it the roundoff residue grows.

    The Bell-directed fixed state is a full-state saddle, so a 1e-16 transverse residue
    is amplified round after round and eventually drives an ARTIFICIAL escape.
    """
    rows = g.run_dense("step4", project=False, n_rounds=N_SHORT)
    leak = [r["transverse_norm"] for r in rows]
    assert leak[1] < 1e-14, leak[1]
    assert leak[N_SHORT] > 50 * leak[1], (leak[1], leak[N_SHORT])
    assert leak[N_SHORT] > leak[50] > leak[10]


def test_projected_dense_map_matches_the_closed_form_uv_recursion():
    for c in ("step3", "step4"):
        dense = g.run_dense(c, project=True, n_rounds=N_SHORT)
        uv = g.run_uv(c, n_rounds=N_SHORT)
        d = max(abs(a["F"] - b["F"]) for a, b in zip(dense, uv))
        assert d < 1e-14, (c, d)


def test_fixed_point_fidelities_match_the_notes():
    """F*(3) = 0.978346, F*(4) = 0.982351 at q = 0.01 (notes 01/02 Sec. 7/8)."""
    for c, target in (("step3", 0.978345804077), ("step4", 0.982350951536)):
        F = g.run_dense(c, project=True, n_rounds=N_SHORT)[-1]["F"]
        assert abs(F - target) < 1e-11, (c, F)


def test_step5_plateau_then_escape():
    """F rises to the 0.98751 plateau, holds it past n = 100, and has left it by n = 1000."""
    rows = g.run_dense("step5cal", project=False, n_rounds=1000)
    F = {r["n"]: r["F"] for r in rows}
    assert abs(F[1] - 0.985465) < 1e-5, F[1]
    assert abs(F[10] - 0.987507) < 1e-5, F[10]
    assert abs(F[100] - 0.987507) < 1e-5, F[100]
    assert abs(F[1000] - 0.5442) < 1e-3, F[1000]          # note 03 Sec. 7: F1000 ~ 0.5442
    assert F[1000] < F[100]


def test_iterates_stay_physical():
    for c, project in (("step3", True), ("step4", True), ("step5cal", False)):
        for r in g.run_dense(c, project, n_rounds=N_SHORT):
            assert r["eig_min"] > -1e-12, (c, r["n"], r["eig_min"])
            assert r["trace_err"] < 1e-12, (c, r["n"])
            assert r["herm_err"] < 1e-12, (c, r["n"])


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} tests passed")
    sys.exit(1 if failed else 0)
