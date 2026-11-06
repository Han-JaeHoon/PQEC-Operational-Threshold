/**
 * Asymptotic target-fidelity threshold q_target(F_target) of the Type-3, Type-4 and 4Q maps.
 *
 * Sources: "Exact Algebraic Calculation of the Asymptotic Target-Fidelity Noise Threshold"
 * (Nov 6, 2026, analytic note, Theorems 1.1 and 5.1) and "The target-fidelity threshold
 * q_target(F_target)" (Nov 2026, Theorem 2).
 *
 *   q_target(F_target) = sup { q : F_inf(eps, q) >= F_target }
 *
 * For 0 < eps < 2/3 and 1/2 < F_target < 1 it is the unique q_* in (0, q_ent) with
 * F_B(q_*) = F_target, F_B the fidelity of the Bell-connected fixed point P_B.  It does not depend
 * on eps, and the supremum is attained (the inequality is >=).  Equivalently
 * q_target(F) = q_gain(4 (1 - F) / 3) (Eq. 59).
 *
 * Computation.  The same certified bisection as q_gain (qgain.fbCrossingBracket): the sign of
 * F_B(q) - F_target from the outward-rounded closed form of P_B (larger root v_B of H_s, positive
 * u_B), so only the physical Bell-connected branch is ever evaluated; saddle, negative-u and
 * extraneous resultant roots never enter.  F_target is taken as the decimal the user typed: its
 * float64 value x is within half an ulp of it, so [x - ulp, x + ulp] encloses it.
 *
 * Algebraic cross-check.  The note's normalized polynomials P_4, P_2 (Eqs. 22-24, degree 14 and
 * 10 in s) and P_4Q (Eqs. 33-34, degree 11), at lambda = 4 F_target - 1, must change sign across
 * the bracket in s = 1 - q.  The signs are evaluated exactly in dyadic rational arithmetic.
 */
import type { Interval } from "./epsc";
import { fbCrossingBracket } from "./qgain";
import { dyadic, type RegimeProtocol } from "./regimes";

// ---------------------------------------------------------------- range of the theorem
export type TargetCoverage = "valid" | "too-low" | "too-high" | "not-a-number";

/** 1/2 < F_target < 1, compared on the float64 value. */
export function targetCoverage(F: number): TargetCoverage {
  if (!Number.isFinite(F)) return "not-a-number";
  if (!(F > 0.5)) return "too-low";
  if (!(F < 1)) return "too-high";
  return "valid";
}

const ULP = 2 ** -52;
/** An enclosure of the decimal that parsed to x (correct rounding: within half an ulp of x). */
export const targetEnclosure = (x: number): Interval => [
  x - (Math.abs(x) * ULP + 2 ** -1074), x + (Math.abs(x) * ULP + 2 ** -1074),
];

// ---------------------------------------------------------------- the bracket
export interface QTargetBracket {
  protocol: RegimeProtocol;
  F: number;
  /** Certified: F_B(lo) > F_target > F_B(hi), hence lo < q_target < hi. */
  lo: number;
  hi: number;
  steps: number;
  /** P_P(s, 4F - 1) changes sign between s = 1 - hi and s = 1 - lo (exact). */
  polyCheck: boolean;
}

const cache = new Map<string, QTargetBracket | null>();
const CACHE_MAX = 2000;

/**
 * Certified bracket of q_target(F) for 1/2 < F < 1, else null.  Synchronous (about 60 closed-form
 * evaluations plus two exact polynomial signs) and cached per (protocol, F).  Independent of eps.
 */
export function qTargetBracket(protocol: RegimeProtocol, F: number): QTargetBracket | null {
  const key = `${protocol}|${F}`;
  const hit = cache.get(key);
  if (hit !== undefined) return hit;
  let out: QTargetBracket | null = null;
  if (targetCoverage(F) === "valid") {
    const b = fbCrossingBracket(protocol, targetEnclosure(F));
    if (b) out = { protocol, F, ...b, polyCheck: polySignChange(protocol, F, b.lo, b.hi) };
  }
  if (cache.size >= CACHE_MAX) cache.delete(cache.keys().next().value as string);
  cache.set(key, out);
  return out;
}

// ---------------------------------------------------------------- verdict for the current q
export type TargetVerdict = "reach" | "miss" | "near";

/**
 * Inside the theorem's range: q <= lo gives F_inf = F_B(q) > F_target, q >= hi gives
 * F_inf < F_target, and in between nothing is claimed ("near": this contains q = q_target, where
 * F_inf = F_target and the target is met).  Requires 0 < q.
 */
export function targetVerdict(b: Pick<QTargetBracket, "lo" | "hi">, q: number): TargetVerdict {
  if (q <= b.lo) return "reach";
  if (q >= b.hi) return "miss";
  return "near";
}

// ---------------------------------------------------------------- exact normalized polynomials
/** Dyadic rational n / 2^k. */
export interface Dy { n: bigint; k: number }
const dy = (x: number): Dy => { const d = dyadic(x); return { n: d.num, k: d.k }; };
const int = (n: number): Dy => ({ n: BigInt(n), k: 0 });
const align = (a: Dy, b: Dy): [bigint, bigint, number] => {
  const k = Math.max(a.k, b.k);
  return [a.n << BigInt(k - a.k), b.n << BigInt(k - b.k), k];
};
export const dAdd = (a: Dy, b: Dy): Dy => { const [x, y, k] = align(a, b); return { n: x + y, k }; };
export const dSub = (a: Dy, b: Dy): Dy => { const [x, y, k] = align(a, b); return { n: x - y, k }; };
export const dMul = (a: Dy, b: Dy): Dy => ({ n: a.n * b.n, k: a.k + b.k });
export const dPow = (a: Dy, e: number): Dy => { let o = int(1); for (let i = 0; i < e; i++) o = dMul(o, a); return o; };
const sgn = (a: Dy) => (a.n > 0n ? 1 : a.n < 0n ? -1 : 0);

/** A, ..., G of Eq. (22) (Type 3: k = 4, Type 4: k = 2) or Eq. (33) (4Q). */
export function normalizedCoeffs(protocol: RegimeProtocol, s: Dy, L: Dy) {
  const one = int(1);
  const ops = dAdd(one, s);
  if (protocol === "q4") {
    const s2 = dPow(s, 2), s3 = dPow(s, 3), s4 = dPow(s, 4), s5 = dPow(s, 5);
    return {
      A: dMul(dMul(int(2), s3), dAdd(one, dMul(int(2), s2))),
      B: dMul(int(-1), dMul(s, dAdd(one, dMul(int(3), s4)))),
      C: dMul(int(2), dSub(dSub(one, s4), s5)),
      D: dMul(int(-3), s2),
      E: dSub(dMul(int(2), dSub(one, s2)), dMul(dMul(int(2), s2), L)),
      G: dMul(s2, dMul(L, L)),
    };
  }
  const k = protocol === "type3" ? 4 : 2;
  const h = dPow(s, k);
  return {
    A: dMul(h, dAdd(dMul(int(5), s), one)),
    B: dMul(dMul(int(2), h), dSub(dSub(dSub(s, one), dPow(s, 2)), dPow(s, 3))),
    C: dMul(ops, dSub(one, dMul(int(2), dPow(s, k + 2)))),
    D: dSub(one, dMul(int(7), s)),
    E: dSub(dMul(int(4), dSub(one, s)), dMul(dMul(int(2), ops), L)),
    G: dMul(ops, dMul(L, L)),
  };
}

/** N = (AG - CD)^2 - (AE - BD)(BG - CE), Eq. (23) / numerator of Eq. (34). */
export function resultantN(protocol: RegimeProtocol, s: Dy, L: Dy): Dy {
  const { A, B, C, D, E, G } = normalizedCoeffs(protocol, s, L);
  const t = dSub(dMul(A, G), dMul(C, D));
  return dSub(dMul(t, t), dMul(dSub(dMul(A, E), dMul(B, D)), dSub(dMul(B, G), dMul(C, E))));
}

/**
 * Exact sign of P_P(s, lambda): P_k = -N / (1+s)^2 (Eq. 24), P_4Q = N / s^3 (Eq. 34); the divisors
 * are positive for 0 < s, so only the sign of N is needed.
 */
export function targetPolySign(protocol: RegimeProtocol, s: Dy, L: Dy): number {
  const n = sgn(resultantN(protocol, s, L));
  return protocol === "q4" ? n : -n;
}

/** lambda = 4 x - 1, exactly. */
export const lambdaOf = (x: number): Dy => dSub(dMul(int(4), dy(x)), int(1));
/** s = 1 - q, exactly. */
export const sOf = (q: number): Dy => dSub(int(1), dy(q));

function polySignChange(protocol: RegimeProtocol, F: number, lo: number, hi: number): boolean {
  const L = lambdaOf(F);
  const a = targetPolySign(protocol, sOf(hi), L);
  const b = targetPolySign(protocol, sOf(lo), L);
  return a !== 0 && b !== 0 && a !== b;
}
