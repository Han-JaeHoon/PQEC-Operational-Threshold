"""
Tests for the iterated noisy-PQEC effective map.

Run:  python tests/test_iterated_noisy_pqec.py     (or: pytest tests/)
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import iterated_noisy_pqec as ip

TS = (0.95, 0.9, 0.8, 0.7, 0.5)
CIRCUITS = ("step3", "step4", "step5")

# Step 5's stored learned parameters realise U only to the training precision
# (delta_U ~ 3e-15 in fidelity => ~2e-7 in amplitude), which floors every
# "exactness" statement about that circuit at ~1e-8.
STEP5_FLOOR = 1e-7


def test_channel_matches_repository_kraus():
    """D_q from the definition == global_depol_kraus (the repository's channel)."""
    worst = ip.selftest(verbose=False)
    assert worst < 1e-14, worst


def test_cnot_counts():
    assert ip.n_cnots("step3") == 16
    assert ip.n_cnots("step4") == 14
    assert ip.n_cnots("step5") == 14


def test_ideal_map_is_rho_squared():
    """q = 0: tau_A must equal rho^2 exactly (Steps 3/4) / to the Step-5 floor."""
    for c in CIRCUITS:
        tol = STEP5_FLOOR if c == "step5" else 1e-13
        for t in TS:
            rho = ip.rho_isotropic(t)
            tau = ip.one_round_tau(rho, 0.0, c)
            assert np.linalg.norm(tau - rho @ rho) < tol, (c, t)


def test_one_round_matches_repository_and_closed_forms():
    """Q, N_Phi and F must reproduce the repo executors and the Step-3/4 closed forms."""
    errs = ip.validate_one_round(verbose=False)
    for c in CIRCUITS:
        for k in ("Q", "N", "F"):
            assert errs["vs_pennylane"][c][k] < 1e-12, (c, k, errs["vs_pennylane"][c][k])
    for c in ("step3", "step4"):
        for k in ("Q", "N", "F"):
            assert errs["vs_analytic"][c][k] < 1e-12, (c, k, errs["vs_analytic"][c][k])


def test_iterates_are_hermitian_unit_trace_and_psd():
    """Every iterate must stay Hermitian, unit-trace and (numerically) PSD, and the
    denominator Q must stay bounded away from zero."""
    for c in CIRCUITS:
        for q in (0.0, 1e-3, 1e-2, 5e-2):
            rho = ip.rho_isotropic(0.9)
            for _ in range(30):
                rho, info = ip.one_round_effective_map(rho, q, c)
                assert rho is not None, (c, q)
                assert info["herm_err"] < 1e-12, (c, q, info["herm_err"])
                assert abs(np.trace(rho) - 1) < 1e-12, (c, q)
                assert info["min_eig"] > -1e-12, (c, q, info["min_eig"])
                assert info["Q"] > 1e-3, (c, q, info["Q"])


def test_fixed_point_residual():
    """Steps 3/4 have exact fixed points; Step 5 only an approximate one (documented)."""
    for c in CIRCUITS:
        for q in (1e-3, 1e-2, 5e-2):
            rs, resid, _, _ = ip.solve_fixed_point(q, c)
            limit = 1e-3 if c == "step5" else 1e-12
            assert resid < limit, (c, q, resid)
            assert np.linalg.eigvalsh(rs).min() > -1e-12, (c, q)


def test_q0_fixed_point_is_the_bell_state():
    """At q = 0 the limiting state must be the pure target Bell state."""
    for c in CIRCUITS:
        out = ip.iterate_effective_map(ip.rho_isotropic(0.9), 0.0, c,
                                       tol=1e-13, max_iter=30)
        assert out["final"]["F"] > 1 - 1e-8, (c, out["final"]["F"])
        assert out["final"]["purity"] > 1 - 1e-8, c


def test_regression_fixed_point_fidelities():
    """Regression values of F_*(q) (see results/iterated_pqec/fixed_points.csv)."""
    expected = {
        ("step3", 1e-3): 0.997871123138, ("step3", 1e-2): 0.978345804077,
        ("step3", 5e-2): 0.881236792358,
        ("step4", 1e-3): 0.998248557358, ("step4", 1e-2): 0.982350951536,
        ("step4", 5e-2): 0.908137508560,
        ("step5", 1e-3): 0.998750074013, ("step5", 1e-2): 0.987507065567,
        ("step5", 5e-2): 0.937304324697,
    }
    for (c, q), F in expected.items():
        rs, _, _, _ = ip.solve_fixed_point(q, c)
        got = float(np.real(np.trace(ip.PHI @ rs)))
        assert abs(got - F) < 1e-8, (c, q, got, F)


def test_ordering_step5_gt_step4_gt_step3():
    """The fixed-point fidelity must follow the one-round threshold ordering."""
    for q in (1e-3, 1e-2, 5e-2):
        F = {}
        for c in CIRCUITS:
            rs, _, _, _ = ip.solve_fixed_point(q, c)
            F[c] = float(np.real(np.trace(ip.PHI @ rs)))
        assert F["step5"] > F["step4"] > F["step3"], (q, F)


def test_fixed_point_is_a_saddle():
    """Structural instability: a mixed fixed point always has coherence amplification
    (lambda_1+lambda_2)/Tr(rho^2) > 1, and the numerical Jacobian agrees (rho(J) > 1)."""
    for c in CIRCUITS:
        for q in (1e-2, 5e-2):
            rs, _, _, _ = ip.solve_fixed_point(q, c)
            coh = ip.coherence_amplification(rs)
            _, sr, _ = ip.jacobian(rs, q, c)
            assert coh > 1.0, (c, q, coh)
            assert sr > 1.0, (c, q, sr)


def test_initial_state_independence_step3_step4():
    """All initial isotropic states flow to the same fixed point (Steps 3/4)."""
    for c in ("step3", "step4"):
        for q in (1e-3, 1e-2):
            finals = []
            for eps in (0.05, 0.1, 0.2, 0.3, 0.5):
                out = ip.iterate_effective_map(ip.rho_isotropic(1 - eps), q, c,
                                               tol=1e-12, max_iter=200)
                finals.append(out["rho"])
            spread = max(np.linalg.norm(finals[i] - finals[j])
                         for i in range(len(finals)) for j in range(i + 1, len(finals)))
            assert spread < 1e-10, (c, q, spread)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    npass = 0
    for fn in fns:
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
            npass += 1
        except AssertionError as e:
            print(f"  FAIL  {fn.__name__}: {e}")
    print(f"\n{npass}/{len(fns)} tests passed")
    sys.exit(0 if npass == len(fns) else 1)
