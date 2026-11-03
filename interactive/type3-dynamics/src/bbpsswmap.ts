/**
 * Bennett et al. (1996) BBPSSW recurrence purification, Protocol A — exact scalar recurrence
 * under the per-CNOT replacement depolarizing channel.
 *
 * Source of the map: verification/bbpssw_stage1 (analytic derivation and frozen circuit), checked
 * here against verification/reference_bbpssw.json.  One round on two identical pairs, wires
 * [A1, B1, A2, B2]: input twirl, CNOT(A1 -> A2), D_q on (A1, A2), CNOT(B1 -> B2), D_q on (B1, B2),
 * Z-measure (A2, B2), keep (A1, B1) on 00 or 11, output twirl.  Twirls are ideal U (x) U* averages.
 *
 * For Bell populations p = (a, b, c, d) in the order (Phi+, Phi-, Psi+, Psi-) and s = (1 - q)^2:
 *
 *   w_q    = s (a^2 + b^2, 2ab, c^2 + d^2, 2cd) + (1 - s)/8 (1, 1, 1, 1)
 *   P_succ = s [(a + b)^2 + (c + d)^2] + (1 - s)/2
 *   before the output twirl: p' = w_q / P_succ
 *
 * The twirls make the state isotropic, p_n = (F_n, h_n, h_n, h_n) with h_n = (1 - F_n)/3, so the
 * protocol is one scalar recurrence, F_0 = 1 - 3 eps/4:
 *
 *   F_{n+1} = [s (F^2 + h^2) + (1 - s)/8] / P_succ,   P_succ = s [(F + h)^2 + 4 h^2] + (1 - s)/2
 *
 * F_n is conditional on every earlier round succeeding; P_succ^(n) belongs to the single round
 * n -> n+1, and nothing here accumulates it into a yield.
 *
 * Fidelity drift, evaluated a second way from a factored numerator (route 2):
 *
 *   F_{n+1} - F_n = w [8 s h (2F - 1) - 3 (1 - s)] / (8 P_succ),   w = (4F - 1)/3,
 *
 * with 1 - s = q (2 - q).  Identity: 1 - 4h = w, 1 - 6h = 2F - 1.
 *
 * Error bounds.  Counted operation by operation (first order, u = 2^-53): the factored numerator
 * is within 12u of `mag`, P_succ within 13u, so route 2 is within 26u mag / (8 P_succ).  F_{n+1}
 * is within 25u of itself, plus one rounding for the subtraction when F_{n+1}, F_n are not close.
 * Both bounds use 32u.  P_succ >= 1/2 for every q and F, so the quotient never blows up, and F >= 1/4
 * keeps every term far from underflow apart from the ETA64 safety term.  These cover one evaluation
 * from the stored F_n, as for the other maps, and are checked against exact rational arithmetic in
 * tests/bbpsswmap.test.ts.  This is a finite-round numerical tool, not a proof.
 */

import { EPS64, ETA64, Step, type Trajectory } from "./type3map";

export type BellPops = readonly [number, number, number, number];

/** One full round from fidelity F: returns [F', P_succ]. */
export function bbpsswStep(F: number, q: number): [number, number] {
  const r = 1 - q;
  const s = r * r;
  const h = (1 - F) / 3;
  const fh = F + h;
  const P = s * (fh * fh + 4 * h * h) + (1 - s) / 2;
  return [(s * (F * F + h * h) + (1 - s) / 8) / P, P];
}

/** P_succ^(n) of the round n -> n+1, from F_n. */
export function bbpsswPSucc(F: number, q: number): number {
  return bbpsswStep(F, q)[1];
}

/** Noisy BXOR + postselection on a Bell-diagonal input, with NO twirl on either side. */
export function bbpsswRoundUntwirled(p: BellPops, q: number): [BellPops, number] {
  const [a, b, c, d] = p;
  const r = 1 - q;
  const s = r * r;
  const bg = (1 - s) / 8;
  const P = s * ((a + b) * (a + b) + (c + d) * (c + d)) + (1 - s) / 2;
  return [[
    (s * (a * a + b * b) + bg) / P,
    (s * 2 * a * b + bg) / P,
    (s * (c * c + d * d) + bg) / P,
    (s * 2 * c * d + bg) / P,
  ], P];
}

/** The original protocol on a Bell-diagonal input: determined by its fidelity p[0] alone. */
export function bbpsswRoundOriginal(p: BellPops, q: number): [BellPops, number] {
  const [F, P] = bbpsswStep(p[0], q);
  const h = (1 - F) / 3;
  return [[F, h, h, h], P];
}

/**
 * Fill `t` with the BBPSSW trajectory for (eps, q) over N displayed rounds (one extra state is
 * computed so that dF_N is defined).  t.u and t.v both hold the correlator (4F - 1)/3, which is
 * <XX> = -<YY> = <ZZ> of the isotropic state.  Reuses the buffers; allocates nothing.
 */
export function computeTrajectoryBBPSSW(t: Trajectory, eps: number, q: number, N: number): Trajectory {
  if (N + 2 > t.u.length) throw new RangeError(`N=${N} exceeds buffer capacity`);
  const r = 1 - q;
  const s = r * r;
  const c = 1 - s;
  const one8 = c / 8;
  const one2 = c / 2;
  const qc = 3 * q * (2 - q); // 3 (1 - s), without the cancellation in 1 - s

  t.eps = eps;
  t.q = q;
  t.N = N;

  let F = 1 - (3 * eps) / 4;
  for (let n = 0; n <= N + 1; n++) {
    const h = (1 - F) / 3;
    const w = (4 * F - 1) / 3;
    t.u[n] = w;
    t.v[n] = w;
    t.pPhiP[n] = F;
    t.pPhiM[n] = h;
    t.pPsiP[n] = h;
    t.pPsiM[n] = h;
    if (n === N + 1) break;

    const fh = F + h;
    const P = s * (fh * fh + 4 * h * h) + one2;
    const Fn = (s * (F * F + h * h) + one8) / P;

    // Route 2: factored drift numerator, with its forward error bound.
    const t1 = 8 * s * h * (2 * F - 1);
    const den = 8 * P;
    t.dFmap[n] = (w * (t1 - qc)) / den;
    t.dFmapErr[n] = (32 * EPS64 * Math.abs(w) * (Math.abs(t1) + qc)) / den + ETA64;

    F = Fn;
  }

  for (let n = 0; n <= N; n++) {
    const f0 = t.pPhiP[n];
    const f1 = t.pPhiP[n + 1];
    const d = f1 - f0;
    t.dF[n] = d;
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
