/**
 * 4Q (four-qubit, 5-CNOT, postselected) noisy Bell-state distillation — exact map on its
 * invariant Bell-diagonal plane.
 *
 * Source of the map: the exact Bell-diagonal recurrence of the 4Q circuit
 * (Bell-state-distillation, src/pqec_distill/analytic_exact_map.py, analytic_map_4q; the
 * calculation notes "Analytic Calculation of the Absolute Bell-Error-Suppression Threshold
 * for the 4Q ..." Eqs. 15-29 and "Analytic Derivation of the Relative ..." Eqs. 7-22).  It
 * agrees with the frozen circuit-only reference (results/data/circuit_reference/
 * circuit_reference.csv, 957 4Q rows) and is checked here against
 * verification/reference_4q.json.  On rho = 1/4 (II + xXX + yYY + zZZ), r = 1 - q:
 *
 *   D  = 1 + r^5 (x^2 + y^2) + r^3 z^2,       P_succ = D / 4
 *   x' = r^3 [(1+r) x - 2 r^2 y z] / D
 *   y' = r^4 [(1+r) y - 2 r   x z] / D
 *   z' = r^4 [(1+r) z - 2 r   x y] / D
 *
 * The 4Q plane is NOT the Type-3/4 plane: the Bell-isotropic input (x, y, z) = (t, -t, t)
 * keeps y = -z (not y = -x).  With a = <XX> = x and b = <ZZ> = z (stored in Trajectory.u, .v,
 * so the inspector labels <XX>, <ZZ> stay correct):
 *
 *   D_n     = 1 + r^5 a^2 + (r^5 + r^3) b^2
 *   a_{n+1} = r^3 [(1+r) a + 2 r^2 b^2] / D_n
 *   b_{n+1} = r^4 b [(1+r) + 2 r a] / D_n
 *
 *   p_Phi+ = (1 + a + 2b)/4,  p_Psi+ = (1 + a - 2b)/4,  p_Phi- = p_Psi- = (1 - a)/4.
 *
 * Fidelity drift, evaluated a second way from the map itself (route 2):
 *
 *   F_{n+1} - F_n = N_n / (4 D_n),   N_n = A + 2B - (a + 2b) D,
 *   A = r^3 (1+r) a + 2 r^5 b^2,  B = r^4 b (1 + r + 2 r a).
 *
 * Error bounds.  Counted operation by operation (first order, u = 2^-53, r = fl(1 - q) itself
 * one rounding): D carries <= 13u relative, A and 2B <= 12u, (a+2b)D <= 15u, so the computed
 * N is within 17u * mag; the quotient adds the 13u of D and one rounding, <= 31u * mag/(4D).
 * a' and b' are each within 26u relative, so F_{n+1} = (1 + a' + 2b')/4 is within 28u of
 * itself (a, b >= 0 on the whole trajectory for 0 <= eps <= 1).  Both bounds use 32u; the
 * absolute 2^-1074 terms cover products that underflow when q -> 1 drives a, b to 0.
 * Checked against exact rational arithmetic in tests/q4map.test.ts.
 * This is a finite-round numerical tool, not a proof.
 */

import { EPS64, ETA64, Step, type Trajectory } from "./type3map";

/** One exact 4Q round on the invariant plane (a, b) = (<XX>, <ZZ>). Returns [a', b', D]. */
export function q4Step(a: number, b: number, q: number): [number, number, number] {
  const r = 1 - q;
  const r2 = r * r;
  const r3 = r2 * r;
  const r4 = r3 * r;
  const r5 = r4 * r;
  const D = 1 + r5 * a * a + (r5 + r3) * b * b;
  return [(r3 * ((1 + r) * a + 2 * r2 * b * b)) / D, (r4 * b * (1 + r + 2 * r * a)) / D, D];
}

export function bellPopulations4Q(a: number, b: number): [number, number, number, number] {
  const m = (1 - a) / 4;
  return [(1 + a + 2 * b) / 4, m, (1 + a - 2 * b) / 4, m];
}

/**
 * Fill `t` with the 4Q trajectory for (eps, q) over N displayed rounds (one extra state is
 * computed so that dF_N is defined).  t.u holds a_n = <XX>, t.v holds b_n = <ZZ>.
 * Reuses the buffers; allocates nothing.
 */
export function computeTrajectory4Q(t: Trajectory, eps: number, q: number, N: number): Trajectory {
  if (N + 2 > t.u.length) throw new RangeError(`N=${N} exceeds buffer capacity`);
  const r = 1 - q;
  const r2 = r * r;
  const r3 = r2 * r;
  const r4 = r3 * r;
  const r5 = r4 * r;
  const c1 = 1 + r;
  const c3 = r3 * c1;
  const c35 = r5 + r3;
  const r5x2 = 2 * r5;

  t.eps = eps;
  t.q = q;
  t.N = N;

  let a = 1 - eps;
  let b = 1 - eps;
  for (let n = 0; n <= N + 1; n++) {
    t.u[n] = a;
    t.v[n] = b;
    const m = (1 - a) / 4;
    t.pPhiP[n] = (1 + a + 2 * b) / 4;
    t.pPhiM[n] = m;
    t.pPsiP[n] = (1 + a - 2 * b) / 4;
    t.pPsiM[n] = m;
    if (n === N + 1) break;

    const bb = b * b;
    const D = 1 + r5 * (a * a) + c35 * bb;
    const A = c3 * a + r5x2 * bb; //            a' = A / D
    const B = r4 * b * (c1 + 2 * r * a); //      b' = B / D
    const an = A / D;
    const bn = B / D;

    // Route 2: drift numerator, N/(4D) = F_{n+1} - F_n, with its forward error bound
    // (see the header for the operation count) plus an absolute term for underflow.
    const s = a + 2 * b;
    const sD = s * D;
    const num = A + 2 * B - sD;
    const mag = Math.abs(c3 * a) + r5x2 * bb + 2 * Math.abs(r4 * b) * (c1 + 2 * r * Math.abs(a)) +
      (Math.abs(a) + 2 * Math.abs(b)) * D;
    t.dFmap[n] = num / (4 * D);
    t.dFmapErr[n] = (32 * EPS64 * mag + 128 * ETA64) / (4 * D) + ETA64;

    a = an;
    b = bn;
  }

  for (let n = 0; n <= N; n++) {
    const f0 = t.pPhiP[n];
    const f1 = t.pPhiP[n + 1];
    const d = f1 - f0;
    t.dF[n] = d;
    // F_{n+1} carries the rounding of a', b' (<= 26u each) and of (1 + a' + 2b')/4; F_n
    // only the latter; the subtraction is exact when the two are close (Sterbenz).
    const err = 32 * EPS64 * (Math.abs(f0) + Math.abs(f1));
    t.dFerr[n] = err;
    const dm = t.dFmap[n];
    const em = t.dFmapErr[n];
    let st = Step.Unresolved;
    if (Math.abs(d) > err && Math.abs(dm) > em && Math.sign(d) === Math.sign(dm)) {
      st = d > 0 ? Step.Increase : Step.Decrease;
    }
    t.step[n] = st;
  }
  return t;
}
