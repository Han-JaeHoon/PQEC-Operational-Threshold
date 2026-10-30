/**
 * Analytic 4Q suppression thresholds for the Bell-isotropic input (calculation notes
 * "Analytic Calculation of the Absolute Bell-Error-Suppression Threshold for the 4Q Noisy
 * Bell-State Distillation Protocol" and "Analytic Derivation of the Relative Bell-Error
 * Suppression Threshold for the Noisy Four-Qubit ...", J. Han, October 2026).
 * With t = 1 - eps and r = 1 - q:
 *
 *   q_abs(eps) = 1 - r*,  2(1 + t - t^2) r^5 + r^4 - (1 + t^2) r^3 - 1 = 0   (abs note, Eqs. 65, 142)
 *   q_rel(eps) = 1 - r*,  1 - r^4 + t(t+1) r^3 + (2t^2 - 3t - 1) r^5 = 0     (rel note, Eqs. 45, 92)
 *
 * Both notes prove first-round = all-round; the limiting component is Psi+ for both (not
 * Phi- as for Type 3 / Type 4).
 *
 * Validity differs between the two:
 *   q_abs: 0 < eps < 2/3 (the proof uses eta = t > 1/3).  f(0) = -1 < 0,
 *          f(1) = (1-t)(1+3t) > 0, one sign change (Descartes): unique root.
 *   q_rel: 0 < eps < 1 (the whole interior; rel note Sec. 14).  f(0) = 1 > 0,
 *          f(1) = (3t+1)(t-1) < 0, f rises then falls (Sec. 14): unique root.
 * Plain bisection on (0, 1) is then exact to the last bit.  Outside its range NaN is returned.
 */

import { EPS_MAX } from "./thresholds";

/** q_rel^(4Q) is defined on the whole interior 0 < eps < 1. */
export const EPS_MAX_REL_4Q = 1;

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

/** 4Q absolute-suppression boundary polynomial f_abs(r; t) (> 0 <=> Psi+ suppressed at round 1). */
export function fAbs4Q(r: number, t: number): number {
  const r3 = r * r * r;
  return r3 * (r * (2 * (1 + t - t * t) * r + 1) - (1 + t * t)) - 1;
}

/** 4Q relative-suppression boundary polynomial f_rel(r; t) (< 0 <=> Psi+ ratio suppressed at round 1). */
export function fRel4Q(r: number, t: number): number {
  const r3 = r * r * r;
  return 1 + r3 * (r * ((2 * t * t - 3 * t - 1) * r - 1) + t * (t + 1));
}

/** q_abs^(4Q)(eps); NaN outside 0 < eps < 2/3. */
export function qAbs4Q(eps: number): number {
  if (!(eps > 0 && eps < EPS_MAX)) return NaN;
  const t = 1 - eps;
  return 1 - bisectUnit((r) => fAbs4Q(r, t));
}

/** q_rel^(4Q)(eps); NaN outside 0 < eps < 1. */
export function qRel4Q(eps: number): number {
  if (!(eps > 0 && eps < EPS_MAX_REL_4Q)) return NaN;
  const t = 1 - eps;
  return 1 - bisectUnit((r) => -fRel4Q(r, t));
}
