/**
 * Type-3 (textbook 16-CNOT) noisy PQEC — exact map on the invariant Bell-diagonal plane.
 *
 * Source of the map: the exact Bell-diagonal recurrence of the Type-3 circuit
 * (RESEARCH_SUMMARY.md §3.3; branch `iterated-noisy-pqec`, iterated_noisy_pqec.py, which
 * verifies it against the dense 32x32 circuit simulation).  With r = 1 - q:
 *
 *   D_n     = 1 + r^4 (2 u_n^2 + v_n^2)
 *   u_{n+1} = 2 r^6 u_n (1 + v_n) / D_n
 *   v_{n+1} = r^5 (1 + r) (u_n^2 + v_n) / D_n
 *
 * on rho = 1/4 [II + u XX - u YY + v ZZ], starting from u_0 = v_0 = 1 - eps.
 * q is the strength of the two-qubit replacement depolarizing channel applied after every
 * CNOT.  Bell populations:
 *
 *   p_Phi+ = (1 + 2u + v)/4,  p_Phi- = (1 - 2u + v)/4,  p_Psi+ = p_Psi- = (1 - v)/4.
 *
 * Fidelity drift.  Besides the direct difference F_{n+1} - F_n, the change is evaluated a
 * second way, from the map itself without forming F_{n+1}:
 *
 *   F_{n+1} - F_n = N_n / (4 D_n),
 *   N_n = 2u (2 r^6 (1+v) - D) + r^5 (1+r)(u^2 + v) - v D.
 *
 * Each route carries a floating-point error bound; a sign is reported only when both routes
 * clear their bounds and agree.  This is a finite-round numerical tool, not a proof.
 */

/** Unit roundoff of IEEE-754 binary64. */
export const EPS64 = 2 ** -53;

/** Classification of one fidelity step F_n -> F_{n+1}. */
export const enum Step {
  /** F_{n+1} > F_n, resolved by both routes. */
  Increase = 1,
  /** F_{n+1} < F_n, resolved by both routes. */
  Decrease = -1,
  /** |change| is at or below the floating-point error bound of at least one route,
   *  or the two routes disagree in sign: no sign is claimed. */
  Unresolved = 0,
}

export interface Trajectory {
  eps: number;
  q: number;
  /** Number of displayed rounds; states are stored for n = 0..N+1 so that dF_N exists. */
  N: number;
  u: Float64Array;
  v: Float64Array;
  pPhiP: Float64Array; // = F_n
  pPhiM: Float64Array;
  pPsiP: Float64Array;
  pPsiM: Float64Array;
  /** dF_n = F_{n+1} - F_n (direct difference), n = 0..N. */
  dF: Float64Array;
  /** Same quantity from the drift numerator N_n/(4 D_n). */
  dFmap: Float64Array;
  /** Rounding-error bounds for the two routes. */
  dFerr: Float64Array;
  dFmapErr: Float64Array;
  /** Step classification (Int8: +1 / -1 / 0). */
  step: Int8Array;
}

/** Allocate a trajectory buffer able to hold up to `maxN` displayed rounds. */
export function allocTrajectory(maxN: number): Trajectory {
  const n = maxN + 2;
  return {
    eps: NaN, q: NaN, N: 0,
    u: new Float64Array(n), v: new Float64Array(n),
    pPhiP: new Float64Array(n), pPhiM: new Float64Array(n),
    pPsiP: new Float64Array(n), pPsiM: new Float64Array(n),
    dF: new Float64Array(n - 1), dFmap: new Float64Array(n - 1),
    dFerr: new Float64Array(n - 1), dFmapErr: new Float64Array(n - 1),
    step: new Int8Array(n - 1),
  };
}

/** One exact Type-3 round on the invariant plane. Returns [u', v', D]. */
export function type3Step(u: number, v: number, q: number): [number, number, number] {
  const r = 1 - q;
  const r2 = r * r;
  const r4 = r2 * r2;
  const r5 = r4 * r;
  const r6 = r5 * r;
  const D = 1 + r4 * (2 * u * u + v * v);
  return [(2 * r6 * u * (1 + v)) / D, (r5 * (1 + r) * (u * u + v)) / D, D];
}

export function bellPopulations(u: number, v: number): [number, number, number, number] {
  const psi = (1 - v) / 4;
  return [(1 + 2 * u + v) / 4, (1 - 2 * u + v) / 4, psi, psi];
}

/**
 * Fill `t` with the trajectory for (eps, q) over N displayed rounds (one extra state is
 * computed so that dF_N is defined).  Reuses the buffers; allocates nothing.
 */
export function computeTrajectory(t: Trajectory, eps: number, q: number, N: number): Trajectory {
  if (N + 2 > t.u.length) throw new RangeError(`N=${N} exceeds buffer capacity`);
  const r = 1 - q;
  const r2 = r * r;
  const r4 = r2 * r2;
  const r5 = r4 * r;
  const r6 = r5 * r;
  const c5 = r5 * (1 + r);

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
    const D = 1 + r4 * (2 * uu + v * v);
    const a = 2 * r6 * (1 + v); // u' = a u / D
    const b = c5 * (uu + v); //    v' = b / D
    const un = (a * u) / D;
    const vn = b / D;

    // Route 2: drift numerator, N/(4D) = F_{n+1} - F_n, with a forward error bound
    // (sum of magnitudes of the cancelling terms, times a small multiple of the roundoff).
    const t1 = 2 * u * (a - D);
    const t2 = b;
    const t3 = v * D;
    const num = t1 + t2 - t3;
    const mag = Math.abs(2 * u) * (Math.abs(a) + D) + Math.abs(t2) + Math.abs(t3);
    t.dFmap[n] = num / (4 * D);
    t.dFmapErr[n] = (16 * EPS64 * mag) / (4 * D);

    u = un;
    v = vn;
  }

  for (let n = 0; n <= N; n++) {
    const f0 = t.pPhiP[n];
    const f1 = t.pPhiP[n + 1];
    const d = f1 - f0;
    t.dF[n] = d;
    // Each F carries a few ulps from (1 + 2u + v)/4; the subtraction is exact (Sterbenz)
    // when F_{n+1}, F_n are close, so the bound is set by the operands themselves.
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

export interface Summary {
  F0: number;
  FN: number;
  maxF: number;
  maxFRound: number;
  minF: number;
  minFRound: number;
  /** First n in [0, N-1] with a resolved decrease F_{n+1} < F_n, or -1. */
  firstDecrease: number;
  /** Number of steps n in [0, N-1] whose sign could not be resolved. */
  unresolved: number;
  /** First unresolved step in [0, N-1], or -1. */
  firstUnresolved: number;
  finalPops: [number, number, number, number];
}

/** Summary over the displayed rounds n = 0..N (steps n = 0..N-1). */
export function summarize(t: Trajectory): Summary {
  const N = t.N;
  let maxF = -Infinity, minF = Infinity, maxR = 0, minR = 0;
  for (let n = 0; n <= N; n++) {
    const f = t.pPhiP[n];
    if (f > maxF) { maxF = f; maxR = n; }
    if (f < minF) { minF = f; minR = n; }
  }
  let firstDecrease = -1, unresolved = 0, firstUnresolved = -1;
  for (let n = 0; n < N; n++) {
    const s = t.step[n];
    if (s === Step.Decrease && firstDecrease < 0) firstDecrease = n;
    if (s === Step.Unresolved) {
      unresolved++;
      if (firstUnresolved < 0) firstUnresolved = n;
    }
  }
  return {
    F0: t.pPhiP[0], FN: t.pPhiP[N],
    maxF, maxFRound: maxR, minF, minFRound: minR,
    firstDecrease, unresolved, firstUnresolved,
    finalPops: [t.pPhiP[N], t.pPhiM[N], t.pPsiP[N], t.pPsiM[N]],
  };
}
