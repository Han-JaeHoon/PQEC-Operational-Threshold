/**
 * Type-4 (resynthesized 14-CNOT, pqec_resynth_noise.GATES) noisy PQEC — exact map on the
 * invariant Bell-diagonal plane.
 *
 * Source of the map: the exact Bell-diagonal recurrence of the Type-4 circuit
 * (Bell-state-distillation, src/pqec_distill/analytic_exact_map.py, k = 2; the calculation
 * notes "Analytic Calculation of the Absolute / Relative Bell-Error-Suppression Threshold for
 * the Type-4 ..." Eqs. 12-28 / 21-39).  It agrees with the frozen circuit-only reference
 * (results/data/circuit_reference/circuit_reference.csv, 957 Type-4 rows) to 2.6e-15 and is
 * checked here against verification/reference_type4.json.  With r = 1 - q:
 *
 *   D_n     = 1 + r^2 (2 u_n^2 + v_n^2)
 *   u_{n+1} = 2 r^4 u_n (1 + v_n) / D_n
 *   v_{n+1} = r^3 (1 + r) (u_n^2 + v_n) / D_n
 *
 * on the same plane as Type 3, rho = 1/4 [II + u XX - u YY + v ZZ], u_0 = v_0 = 1 - eps
 * (y = -x is preserved by the general map: y' = -x' whenever y = -x).  Bell populations are
 * the Type-3 ones, p_Phi+- = (1 +- 2u + v)/4, p_Psi+ = p_Psi- = (1 - v)/4.
 *
 * Fidelity drift.  As for Type 3, F_{n+1} - F_n is also evaluated from the map itself:
 *
 *   F_{n+1} - F_n = N_n / (4 D_n),
 *   N_n = 2u (2 r^4 (1+v) - D) + r^3 (1+r)(u^2 + v) - v D,
 *
 * with its own forward error bound; the operation count per term is no larger than in the
 * Type-3 evaluation, and the bound is checked against exact rational arithmetic in
 * tests/type4map.test.ts.  This is a finite-round numerical tool, not a proof.
 */

import { EPS64, ETA64, Step, type Trajectory } from "./type3map";

/** One exact Type-4 round on the invariant plane. Returns [u', v', D]. */
export function type4Step(u: number, v: number, q: number): [number, number, number] {
  const r = 1 - q;
  const r2 = r * r;
  const r3 = r2 * r;
  const r4 = r2 * r2;
  const D = 1 + r2 * (2 * u * u + v * v);
  return [(2 * r4 * u * (1 + v)) / D, (r3 * (1 + r) * (u * u + v)) / D, D];
}

/**
 * Fill `t` with the Type-4 trajectory for (eps, q) over N displayed rounds (one extra state
 * is computed so that dF_N is defined).  Reuses the buffers; allocates nothing.
 */
export function computeTrajectory4(t: Trajectory, eps: number, q: number, N: number): Trajectory {
  if (N + 2 > t.u.length) throw new RangeError(`N=${N} exceeds buffer capacity`);
  const r = 1 - q;
  const r2 = r * r;
  const r3 = r2 * r;
  const r4 = r2 * r2;
  const c3 = r3 * (1 + r);

  t.eps = eps;
  t.q = q;
  t.N = N;

  let u = 1 - eps;
  let v = 1 - eps;
  for (let n = 0; n <= N + 1; n++) {
    t.u[n] = u;
    t.v[n] = v;
    const psi = (1 - v) / 4;
    t.pPhiP[n] = (1 + 2 * u + v) / 4;
    t.pPhiM[n] = (1 - 2 * u + v) / 4;
    t.pPsiP[n] = psi;
    t.pPsiM[n] = psi;
    if (n === N + 1) break;

    const uu = u * u;
    const D = 1 + r2 * (2 * uu + v * v);
    const a = 2 * r4 * (1 + v); // u' = a u / D
    const b = c3 * (uu + v); //    v' = b / D
    const un = (a * u) / D;
    const vn = b / D;

    // Route 2: drift numerator, N/(4D) = F_{n+1} - F_n, with a forward error bound
    // (sum of magnitudes of the cancelling terms, times a small multiple of the roundoff),
    // plus an absolute term for products that underflow (at most 7 of them reach N, each off
    // by <= 2^-1075 and amplified by < 8) and for the final quotient.
    const t1 = 2 * u * (a - D);
    const t2 = b;
    const t3 = v * D;
    const num = t1 + t2 - t3;
    const mag = Math.abs(2 * u) * (Math.abs(a) + D) + Math.abs(t2) + Math.abs(t3);
    t.dFmap[n] = num / (4 * D);
    t.dFmapErr[n] = (16 * EPS64 * mag + 64 * ETA64) / (4 * D) + ETA64;

    u = un;
    v = vn;
  }

  for (let n = 0; n <= N; n++) {
    const f0 = t.pPhiP[n];
    const f1 = t.pPhiP[n + 1];
    const d = f1 - f0;
    t.dF[n] = d;
    // Same structure as Type 3: each F carries a few ulps from (1 + 2u + v)/4 and from the
    // rounded u', v' of the preceding step; the subtraction itself is exact when close.
    const err = 8 * EPS64 * (Math.abs(f0) + Math.abs(f1));
    t.dFerr[n] = err;
    const dm = t.dFmap[n];
    const em = t.dFmapErr[n];
    let s = Step.Unresolved;
    if (Math.abs(d) > err && Math.abs(dm) > em && Math.sign(d) === Math.sign(dm)) {
      s = d > 0 ? Step.Increase : Step.Decrease;
    }
    t.step[n] = s;
  }
  return t;
}
