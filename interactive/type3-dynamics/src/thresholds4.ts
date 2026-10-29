/**
 * Analytic Type-4 suppression thresholds for the Bell-isotropic input (calculation notes
 * "Analytic Calculation of the Absolute Bell-Error-Suppression Threshold for the Type-4 ..."
 * and "... Relative ...", J. Han, October 2026).  With t = 1 - eps and r = 1 - q:
 *
 *   q_abs(eps) = 1 - r*,   3(1+t) r^4 - (1+t) r^3 - 3 t^2 r^2 - 1 = 0          (abs note, Eq. 191)
 *   q_rel(eps) = 1 - r*,   (1+t)(2+t) r^4 - t(1+t) r^3 - 3 t^2 r^2 - 1 = 0     (rel note, Eq. 225)
 *
 * Same definitions as the Type-3 thresholds (./thresholds.ts) with every power of r lowered
 * by two.  Both notes prove first-round = all-round and that the limiting component is Phi-.
 *
 * Validity: entangled input, 0 < eps < 2/3 (1/3 < t < 1).  There f(0) = -1 < 0,
 * f(1) = (1-t)(1+3t) > 0, and Descartes' rule gives a unique positive root, so plain
 * bisection on (0, 1) is exact to the last bit.  Outside that range NaN is returned.
 */

import { EPS_MAX } from "./thresholds";

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

/** Type-4 absolute-suppression boundary polynomial f_abs(r; t). */
export function fAbs4(r: number, t: number): number {
  const r2 = r * r;
  return r2 * (r * (3 * (1 + t) * r - (1 + t)) - 3 * t * t) - 1;
}

/** Type-4 relative-suppression boundary polynomial f_rel(r; t) (= F_Phi of the rel note). */
export function fRel4(r: number, t: number): number {
  const r2 = r * r;
  return r2 * (r * ((1 + t) * (2 + t) * r - t * (1 + t)) - 3 * t * t) - 1;
}

/** q_abs^(4)(eps); NaN outside 0 < eps < 2/3. */
export function qAbs4(eps: number): number {
  if (!inRange(eps)) return NaN;
  const t = 1 - eps;
  return 1 - bisectUnit((r) => fAbs4(r, t));
}

/** q_rel^(4)(eps); NaN outside 0 < eps < 2/3. */
export function qRel4(eps: number): number {
  if (!inRange(eps)) return NaN;
  const t = 1 - eps;
  return 1 - bisectUnit((r) => fRel4(r, t));
}
