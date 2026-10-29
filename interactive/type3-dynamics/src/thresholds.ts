/**
 * Analytic Type-3 suppression thresholds for the Bell-isotropic input (calculation notes
 * "Analytic Calculation of the Absolute-Suppression Threshold" and "... Relative
 * Bell-Error-Suppression Threshold", J. Han, October 2026).  With t = 1 - eps and r = 1 - q:
 *
 *   q_abs(eps) = 1 - r*,   3(1+t) r^6 - (1+t) r^5 - 3 t^2 r^4 - 1 = 0          (abs note, Eq. 86)
 *   q_rel(eps) = 1 - r*,   (1+t)(2+t) r^6 - t(1+t) r^5 - 3 t^2 r^4 - 1 = 0     (rel note, Eq. 225)
 *
 * q_abs: largest q for which every non-target Bell population decreases absolutely at
 * every round; q_rel: largest q for which every ratio p_i / p_Phi+ decreases at every round.
 * Both notes prove first-round = all-round and that the limiting component is Phi-.
 *
 * Validity: entangled input, 0 < eps < 2/3 (1/3 < t < 1).  There f(0) = -1 < 0,
 * f(1) = (1-t)(1+3t) > 0, and Descartes' rule gives a unique positive root, so plain
 * bisection on (0, 1) is exact to the last bit.  Outside that range NaN is returned.
 */

export const EPS_MAX = 2 / 3;

/** Unique root r* in (0, 1) of a polynomial with f(0) < 0 < f(1), by bisection. */
function bisectUnit(f: (r: number) => number): number {
  let lo = 0, hi = 1;
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    if (mid <= lo || mid >= hi) break; // interval is one ulp wide
    if (f(mid) > 0) hi = mid; else lo = mid;
  }
  return 0.5 * (lo + hi);
}

function inRange(eps: number): boolean {
  return eps > 0 && eps < EPS_MAX;
}

/** Absolute-suppression boundary polynomial f_abs(r; t). */
export function fAbs(r: number, t: number): number {
  const r4 = r * r * r * r;
  return r4 * (r * (3 * (1 + t) * r - (1 + t)) - 3 * t * t) - 1;
}

/** Relative-suppression boundary polynomial f_rel(r; t) (= F_Phi of the rel note). */
export function fRel(r: number, t: number): number {
  const r4 = r * r * r * r;
  return r4 * (r * ((1 + t) * (2 + t) * r - t * (1 + t)) - 3 * t * t) - 1;
}

/** q_abs^(3)(eps); NaN outside 0 < eps < 2/3. */
export function qAbs(eps: number): number {
  if (!inRange(eps)) return NaN;
  const t = 1 - eps;
  return 1 - bisectUnit((r) => fAbs(r, t));
}

/** q_rel^(3)(eps); NaN outside 0 < eps < 2/3. */
export function qRel(eps: number): number {
  if (!inRange(eps)) return NaN;
  const t = 1 - eps;
  return 1 - bisectUnit((r) => fRel(r, t));
}
