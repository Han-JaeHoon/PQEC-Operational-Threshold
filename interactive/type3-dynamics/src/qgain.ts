/**
 * Asymptotic fidelity-gain threshold q_gain(eps) of the Type-3, Type-4 and 4Q maps.
 *
 * Source: "Exact Algebraic Determination of the Asymptotic Fidelity-Gain Threshold" (Nov 6, 2026).
 *
 *   q_gain(eps) = sup { q : F_inf(eps, q) > F_0(eps) },   F_0 = 1 - 3 eps / 4.
 *
 * Theorem 7.2: for 0 < eps < 2/3, q_gain(eps) = q_F(eps), the unique q in (0, q_ent) with
 * F_B(q) = F_0(eps), where F_B is the fidelity of the Bell-connected fixed point P_B (Eq. 12) and
 * F_B(q_ent) = 1/2.  F_B is strictly decreasing (Lemma 7.1), and the gain set is exactly
 * (0, q_F): no recurrence, no eps_c(q) and no table is needed.
 *
 * Computation.  F_B(q) is enclosed from its closed form with outward rounding (epsc.ts,
 * bellFixedPoint).  Bisection in q on the certified sign of F_B(q) - F_0 starts from
 * [0, upper end of the q_ent isolating interval]: F_B(0) = 1 > F_0, and F_B < 1/2 < F_0 at the
 * upper end, which is checked.  Every end point it keeps has a certified sign, so lo < q_gain < hi.
 * It stops when the midpoint is a neighbour of an end point or its sign is not certified.
 *
 * q_ent is the root of R_P(s) (Eqs. 22-24) in the isolating interval of Section 3; which side of it
 * a float q lies on is decided by the exact polynomial sign (regimes.sideOfRoot).
 */
import { bellFixedPoint, type Interval } from "./epsc";
import { sideOfRoot, solveRoot, type RegimeProtocol, type ThresholdRoot } from "./regimes";

const ENT_POLY: Record<RegimeProtocol, { poly: number[]; bracket: [number, number] }> = {
  // 112 s^8 - 16 s^7 - 19 s^6 + 10 s^5 - 3 s^4 - 49 s^2 + 14 s - 1
  type3: { poly: [112, -16, -19, 10, -3, 0, -49, 14, -1], bracket: [0.8708, 0.8709] },
  // 112 s^6 - 16 s^5 - 19 s^4 + 10 s^3 - 52 s^2 + 14 s - 1
  type4: { poly: [112, -16, -19, 10, -52, 14, -1], bracket: [0.815, 0.816] },
  // 271 s^11 + 176 s^10 - 128 s^9 - 76 s^8 + 10 s^7 - 184 s^6 + 76 s^5 + 80 s^4 - s^3 + 8 s^2 - 36 s - 4
  q4: { poly: [271, 176, -128, -76, 10, -184, 76, 80, -1, 8, -36, -4], bracket: [0.82, 0.821] },
};

/** q_ent with F_B(q_ent) = 1/2, as a ThresholdRoot: `q` is for display, `poly` decides sides. */
export const Q_ENT: Record<RegimeProtocol, ThresholdRoot> = Object.fromEntries(
  (Object.keys(ENT_POLY) as RegimeProtocol[]).map((p) => {
    const r = { ...ENT_POLY[p], q: NaN };
    return [p, { ...r, q: solveRoot(r) }];
  }),
) as unknown as Record<RegimeProtocol, ThresholdRoot>;

// ---------------------------------------------------------------- eps range of the theorem
export type EpsCoverage = "zero" | "covered" | "not-covered";

/** Exactly: 0 < eps < 2/3 is covered.  3 eps < 2 is decided without rounding. */
export function epsCoverage(eps: number): EpsCoverage {
  if (!(eps > 0)) return "zero";
  // eps = n / 2^k exactly; 3 eps < 2  <=>  3 n < 2^(k+1)
  let k = 0;
  let x = eps;
  while (!Number.isInteger(x)) { x *= 2; k++; }
  return 3n * BigInt(x) < 2n << BigInt(k) ? "covered" : "not-covered";
}

// ---------------------------------------------------------------- F_B and the bracket
const ULP = 2 ** -52;
const up = (x: number) => x + (Math.abs(x) * ULP + 2 ** -1074);
const dn = (x: number) => x - (Math.abs(x) * ULP + 2 ** -1074);

/** Certified enclosure of F_0 = 1 - 3 eps / 4. */
export const f0Enclosure = (eps: number): Interval => [dn(1 - up(0.75 * eps)), up(1 - dn(0.75 * eps))];

/** Certified enclosure of F_B(q), or null if P_B cannot be enclosed at this q. */
export function fbEnclosure(protocol: RegimeProtocol, q: number): Interval | null {
  const pb = bellFixedPoint(protocol, q);
  return pb ? pb.F : null;
}

/** Sign of F_B(q) - F_0: 1 above, -1 below, 0 not certified. */
function gainSign(protocol: RegimeProtocol, q: number, F0: Interval): number {
  const F = fbEnclosure(protocol, q);
  if (!F) return 0;
  if (F[0] > F0[1]) return 1;
  if (F[1] < F0[0]) return -1;
  return 0;
}

export interface QGainBracket {
  protocol: RegimeProtocol;
  eps: number;
  /** Certified: F_B(lo) > F_0 > F_B(hi), hence lo < q_gain(eps) < hi. */
  lo: number;
  hi: number;
  steps: number;
}

const cache = new Map<string, QGainBracket | null>();
const CACHE_MAX = 5000;

/**
 * Certified bracket of q_gain(eps) for 0 < eps < 2/3, else null.  Synchronous (about 60
 * closed-form evaluations) and cached per (protocol, eps).
 */
export function qGainBracket(protocol: RegimeProtocol, eps: number): QGainBracket | null {
  const key = `${protocol}|${eps}`;
  const hit = cache.get(key);
  if (hit !== undefined) return hit;
  let out: QGainBracket | null = null;
  if (epsCoverage(eps) === "covered") {
    const F0 = f0Enclosure(eps);
    let lo = 0; // F_B(0) = 1 > F_0 (Lemma 7.1)
    let hi = 1 - Q_ENT[protocol].bracket[0]; // q_ent < hi < q_SN
    let steps = 0;
    if (gainSign(protocol, hi, F0) === -1) {
      for (;;) {
        const mid = lo + (hi - lo) / 2;
        if (mid <= lo || mid >= hi) break;
        const sg = gainSign(protocol, mid, F0);
        steps++;
        if (sg === 0) break; // F_B(mid) and F_0 not separated: resolution limit
        if (sg > 0) lo = mid; else hi = mid;
      }
      out = { protocol, eps, lo, hi, steps };
    }
  }
  if (cache.size >= CACHE_MAX) cache.delete(cache.keys().next().value as string);
  cache.set(key, out);
  return out;
}

// ---------------------------------------------------------------- verdict for the current q
export type GainVerdict = "gain" | "no-gain" | "uncertain";

/**
 * By Theorem 7.2 the gain set is (0, q_gain): q <= lo gains, q >= hi does not, in between nothing
 * is claimed (this covers q = q_gain, which a float q cannot be shown to equal).  Requires 0 < q.
 */
export function gainVerdict(b: Pick<QGainBracket, "lo" | "hi">, q: number): GainVerdict {
  if (q <= b.lo) return "gain";
  if (q >= b.hi) return "no-gain";
  return "uncertain";
}

/** Side of q_ent for this float q, exactly: -1 below, 1 above. */
export const sideOfQEnt = (protocol: RegimeProtocol, q: number) => sideOfRoot(Q_ENT[protocol], q);
