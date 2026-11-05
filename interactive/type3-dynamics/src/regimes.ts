/**
 * Global attractor regimes of the Type-3, Type-4 and 4Q maps: the protocol constants q_bi, q_SN
 * and q_m, and the limit of an orbit for given (eps, q).
 *
 * With s = 1 - q and the Bell-isotropic input:
 *
 *   q < q_bi            every eps in (0,1) converges to the same purified fixed point P_B
 *   q_bi <= q <= q_SN   bistable: small eps goes to P_B, eps near 1 goes to the low attractor
 *   q > q_SN            every eps in (0,1) converges to the same separable point Z(q), F <= 1/2
 *
 * The constants are proved analytically (see THRESHOLDS for each polynomial).  They are not
 * functions of eps, and they are never re-estimated from a simulation.  The tests check each
 * one against bisection of its polynomial.
 *
 * Z(q): for q < q_m the boundary fixed point B* on the invariant line (transverse = 0),
 *   Type 3/4:  v* = sqrt((s^(A-1) (1+s) - 1) / s^B),  F = (1 + v*)/4    (A,B) = (6,4), (4,2)
 *   4Q:        x* = sqrt((s^4 + s^3 - 1) / s^5),      F = (1 + x*)/4
 * and for q >= q_m the mixed point, F = 1/4.  Only 4Q has q_SN < q_m, so only 4Q has Z = B*.
 *
 * Orbit limit.  Convergence can take far more rounds than the charts show, and the transverse
 * component (u for Type 3/4, z for 4Q) can shrink a long way before it regrows, so the limit
 * orbit carries the transverse component as its logarithm tau:
 *
 *   Type 3/4:  D = 1 + s^B (2 w + L^2),  L' = s^(A-1) (1+s) (L + w) / D,  g = 2 s^A (1 + L) / D
 *   4Q:        D = 1 + s^5 L^2 + (s^3 + s^5) w,  L' = (s^3 (1+s) L + 2 s^5 w) / D,
 *              g = (s^4 (1+s) + 2 s^5 L) / D
 *
 * where L is the line coordinate (v or x), w = exp(2 tau) and tau' = tau + log g.  This is the
 * map of type3map.ts / type4map.ts / q4map.ts rewritten, not a new one: w may underflow to 0,
 * which is harmless next to 1 and L, but tau stays finite, so the orbit cannot be trapped on the
 * invariant line by underflow.
 *
 * A small |change| alone does not mean convergence, because on a plateau the orbit barely moves
 * while a tiny transverse component grows.  A tiny component that is still growing therefore
 * blocks the verdict, and rounds run out as "not converged".  A finite-round numerical tool,
 * not a proof.
 */

export type RegimeProtocol = "type3" | "type4" | "q4";
export type Regime = "unique-purified" | "bistable" | "unique-unpurified";

export interface ThresholdRoot {
  /** Proved value of q. */
  q: number;
  /** Polynomial in s = 1 - q, highest degree first. */
  poly: readonly number[];
  /** Bracket in s containing the relevant root. */
  bracket: readonly [number, number];
}

export interface ProtocolThresholds {
  qBi: ThresholdRoot;
  qSN: ThresholdRoot;
  qM: ThresholdRoot;
}

export const THRESHOLDS: Record<RegimeProtocol, ProtocolThresholds> = {
  q4: {
    // s^5 - 4s^4 - 6s^3 + s + 4
    qBi: { q: 0.1758332653, poly: [1, -4, -6, 0, 1, 4], bracket: [0.824, 0.825] },
    // 41s^8 + 32s^7 + 16s^6 + 16s^5 + 6s^4 - 32s^3 - 16s + 1 (largest root in (0,1))
    qSN: { q: 0.180669726, poly: [41, 32, 16, 16, 6, -32, 0, -16, 1], bracket: [0.8193, 0.8194] },
    // s^4 + s^3 - 1
    qM: { q: 0.1808274866, poly: [1, 1, 0, 0, -1], bracket: [0.81, 0.83] },
  },
  type4: {
    // 4s^4 + 4s^3 - s^2 + 2s - 5
    qBi: { q: 0.1789326522, poly: [4, 4, -1, 2, -5], bracket: [0.821, 0.822] },
    // s^8 + 2s^7 + 9s^6 + 12s^5 + 5s^4 - 2s^3 - 4s^2 - 6s - 1 (largest root in (0,1))
    qSN: { q: 0.1894166559, poly: [1, 2, 9, 12, 5, -2, -4, -6, -1], bracket: [0.8105, 0.8106] },
    // s^4 + s^3 - 1
    qM: { q: 0.1808274866, poly: [1, 1, 0, 0, -1], bracket: [0.81, 0.83] },
  },
  type3: {
    // 4s^6 + 4s^5 - s^4 + 2s^3 - s^2 - 4
    qBi: { q: 0.1182907629, poly: [4, 4, -1, 2, -1, 0, -4], bracket: [0.8816, 0.8818] },
    // s^10 + 2s^9 + 9s^8 + 12s^7 + 5s^6 - 2s^5 + s^4 - 5s^2 - 6s - 1 (largest root in (0,1))
    qSN: { q: 0.130579268, poly: [1, 2, 9, 12, 5, -2, 1, 0, -5, -6, -1], bracket: [0.8694, 0.8695] },
    // s^6 + s^5 - 1
    qM: { q: 0.1187285384, poly: [1, 1, 0, 0, 0, 0, -1], bracket: [0.875, 0.89] },
  },
};

/** Horner evaluation, coefficients highest degree first. */
export function evalPoly(poly: readonly number[], s: number): number {
  let out = 0;
  for (const c of poly) out = out * s + c;
  return out;
}

/** q = 1 - s* for the root s* of `root.poly` in `root.bracket`, by bisection. */
export function solveRoot(root: ThresholdRoot): number {
  let [lo, hi] = root.bracket;
  const fLo = evalPoly(root.poly, lo);
  if (fLo * evalPoly(root.poly, hi) > 0) throw new RangeError("bracket has no sign change");
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    if (mid <= lo || mid >= hi) break; // interval is one ulp wide
    if (evalPoly(root.poly, mid) * fLo > 0) lo = mid; else hi = mid;
  }
  return 1 - 0.5 * (lo + hi);
}

export function regime(protocol: RegimeProtocol, q: number): Regime {
  const t = THRESHOLDS[protocol];
  if (q < t.qBi.q) return "unique-purified";
  return q <= t.qSN.q ? "bistable" : "unique-unpurified";
}

// ---------------------------------------------------------------- exact side of a threshold
/** A finite float64 as an exact dyadic rational num / 2^k. */
function dyadic(x: number): { num: bigint; k: number } {
  if (x === 0) return { num: 0n, k: 0 };
  const dv = new DataView(new ArrayBuffer(8));
  dv.setFloat64(0, x);
  const bits = dv.getBigUint64(0);
  const exp = Number((bits >> 52n) & 0x7ffn);
  const frac = bits & 0xfffffffffffffn;
  const mant = exp === 0 ? frac : frac | (1n << 52n);
  const e = (exp === 0 ? 1 : exp) - 1075;
  const signed = bits >> 63n ? -mant : mant;
  return e >= 0 ? { num: signed << BigInt(e), k: 0 } : { num: signed, k: -e };
}

/** Exact sign of the polynomial at the rational s = num / 2^k (integer coefficients). */
function polySignAt(poly: readonly number[], num: bigint, k: number): number {
  // Horner on the numerator of poly(s) * 2^(k * degree)
  const K = BigInt(k);
  let acc = 0n;
  let scale = 0n;
  for (const c of poly) {
    acc = acc * num + (BigInt(c) << scale);
    scale += K;
  }
  return acc > 0n ? 1 : acc < 0n ? -1 : 0;
}

/**
 * Exact sign of q - q_root for this float64 q: 1 above, -1 below, 0 exactly on the root.  It uses
 * the polynomial and its bracket, not the rounded decimal `root.q`.  A float64 q is a dyadic
 * rational and these polynomials have no dyadic root in their brackets, so 0 does not occur for
 * the protocol constants (verification/verify_epsc_polys.py).
 */
export function sideOfRoot(root: ThresholdRoot, q: number): number {
  const { num: qn, k } = dyadic(q);
  const sNum = (1n << BigInt(k)) - qn; // s = 1 - q exactly
  const [loS, hiS] = root.bracket;
  const lo = dyadic(loS), hi = dyadic(hiS);
  const cmp = (a: { num: bigint; k: number }) => {
    const K = Math.max(a.k, k);
    return (sNum << BigInt(K - k)) - (a.num << BigInt(K - a.k));
  };
  if (cmp(lo) <= 0n) return 1; //  s <= bracket low end: s below the root, q above it
  if (cmp(hi) >= 0n) return -1;
  const here = polySignAt(root.poly, sNum, k);
  if (here === 0) return 0;
  return here === polySignAt(root.poly, hi.num, hi.k) ? -1 : 1; // same sign as hi end: s above root
}

export type ExactRegime = Regime | "on-threshold";

export interface RegimeReport {
  regime: ExactRegime;
  /** q < q_m exactly (so Z(q) = B*), otherwise Z(q) is the mixed point. */
  belowQm: boolean;
  /** q agrees with this constant up to the digits it is quoted to. */
  near: "qBi" | "qSN" | null;
  /** 0 < q < 1: the range of the theorems. */
  inRange: boolean;
}

/** |q - constant| below this counts as "at the threshold" for display. */
export const BOUNDARY_TOL = 1e-9;

/** Regime of this exact float64 q, from the sign of the polynomials. */
export function regimeReport(protocol: RegimeProtocol, q: number): RegimeReport {
  const t = THRESHOLDS[protocol];
  const bi = sideOfRoot(t.qBi, q);
  const sn = sideOfRoot(t.qSN, q);
  const regime: ExactRegime = bi === 0 || sn === 0 ? "on-threshold"
    : bi < 0 ? "unique-purified" : sn < 0 ? "bistable" : "unique-unpurified";
  return {
    regime,
    belowQm: sideOfRoot(t.qM, q) < 0,
    near: Math.abs(q - t.qBi.q) < BOUNDARY_TOL ? "qBi" : Math.abs(q - t.qSN.q) < BOUNDARY_TOL ? "qSN" : null,
    inRange: q > 0 && q < 1,
  };
}

/** (A, B) of the Type-3/4 map. */
const AB = { type3: [6, 4], type4: [4, 2] } as const;

export interface ZPoint {
  kind: "boundary" | "mixed";
  /** Line coordinate of the point: v* (Type 3/4) or x* (4Q); 0 at the mixed point. */
  line: number;
  F: number;
}

/** The separable attractor Z(q): B* for q < q_m, otherwise the mixed point. */
export function zPoint(protocol: RegimeProtocol, q: number): ZPoint {
  if (sideOfRoot(THRESHOLDS[protocol].qM, q) >= 0) return { kind: "mixed", line: 0, F: 0.25 };
  const s = 1 - q;
  let line: number;
  if (protocol === "q4") {
    line = Math.sqrt(Math.max(0, s ** 4 + s ** 3 - 1) / s ** 5);
  } else {
    const [A, B] = AB[protocol];
    line = Math.sqrt(Math.max(0, s ** (A - 1) * (1 + s) - 1) / s ** B);
  }
  return { kind: "boundary", line, F: (1 + line) / 4 };
}

/** Rounds available to the limit orbit before it is reported as not converged. */
export const LIMIT_MAX_ROUNDS = 200_000;

/** Below this the transverse component is "tiny": roughly exp(-30) = 1e-13 and exp(-40) = 4e-18. */
const TAU_PLATEAU = -30;
const TAU_GONE = -40;
const TOL_LINE = 1e-14;
const TOL_TAU = 1e-13;
/** w = exp(2 tau) is exactly 0 from here down. */
const TAU_W_ZERO = -400;
/** tau below this: a plain float64 orbit would be at or below the subnormal range. */
const TAU_PLAIN_FLOOR = Math.log(2) * -1022;

export interface OrbitLimit {
  /** Fidelity at the last round computed. */
  F: number;
  /** Line coordinate v (Type 3/4) or x (4Q). */
  line: number;
  /** Transverse coordinate u (Type 3/4) or z (4Q); may read 0 once below float64 range. */
  transverse: number;
  /** log10 of the transverse coordinate (finite even when `transverse` underflows). */
  log10Transverse: number;
  rounds: number;
  converged: boolean;
  /** Smallest log10 of the transverse coordinate along the orbit. */
  minLog10Transverse: number;
  /** The orbit went below the range of a plain float64 transverse component. */
  belowFloat64: boolean;
}

/** Limit of the orbit from the isotropic input u0 = v0 = 1 - eps (x0 = z0 = 1 - eps for 4Q). */
export function orbitLimit(
  protocol: RegimeProtocol, eps: number, q: number, maxRounds = LIMIT_MAX_ROUNDS,
): OrbitLimit {
  return iterateLogOrbit(protocol, 1 - eps, Math.log(1 - eps), q, maxRounds);
}

/** Iterate (L, tau) from an arbitrary start until convergence or `maxRounds`. */
export function iterateLogOrbit(
  protocol: RegimeProtocol, L0: number, tau0: number, q: number, maxRounds = LIMIT_MAX_ROUNDS,
): OrbitLimit {
  const s = 1 - q;
  const s2 = s * s;
  const s3 = s2 * s;
  const s4 = s2 * s2;
  const s5 = s4 * s;
  const is4Q = protocol === "q4";
  const [A, B] = is4Q ? [0, 0] : AB[protocol];
  const sB = is4Q ? 0 : s ** B;
  const sA = is4Q ? 0 : s ** A;
  const line1 = is4Q ? s3 * (1 + s) : s ** (A - 1) * (1 + s);
  const g4 = s4 * (1 + s);

  let L = L0;
  let tau = tau0;
  let minTau = tau;
  let rounds = 0;
  let converged = false;
  while (rounds < maxRounds) {
    const w = tau < TAU_W_ZERO ? 0 : Math.exp(2 * tau);
    let Ln: number, g: number;
    if (is4Q) {
      const D = 1 + s5 * L * L + (s3 + s5) * w;
      Ln = (line1 * L + 2 * s5 * w) / D;
      g = (g4 + 2 * s5 * L) / D;
    } else {
      const D = 1 + sB * (2 * w + L * L);
      Ln = (line1 * (L + w)) / D;
      g = (2 * sA * (1 + L)) / D;
    }
    const tauN = tau + Math.log(g);
    const dL = Ln - L;
    const dTau = tauN - tau; // NaN once tau is -Infinity (eps = 1, or g = 0 at q = 1)
    L = Ln;
    tau = tauN;
    rounds++;
    if (tau < minTau) minTau = tau;

    const lineDone = L > 1e-6 ? Math.abs(dL) < TOL_LINE : dL <= 0 && L < TOL_LINE;
    const gone = tau === -Infinity || (dTau < 0 && tau < TAU_GONE);
    const blocked = dTau > 0 && tau < TAU_PLATEAU;
    if (lineDone && !blocked && (gone || Math.abs(dTau) < TOL_TAU)) {
      converged = true;
      break;
    }
  }
  const transverse = Math.exp(tau);
  return {
    F: is4Q ? (1 + L + 2 * transverse) / 4 : (1 + 2 * transverse + L) / 4,
    line: L,
    transverse,
    log10Transverse: tau / Math.LN10,
    rounds,
    converged,
    minLog10Transverse: minTau / Math.LN10,
    belowFloat64: minTau < TAU_PLAIN_FLOOR,
  };
}
