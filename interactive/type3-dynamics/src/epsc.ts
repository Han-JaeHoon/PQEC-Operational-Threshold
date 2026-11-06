/**
 * Certified brackets of the separating noise eps_c(q) in the bistable window q_bi < q < q_SN
 * of the Type-3, Type-4 and 4Q maps.
 *
 * Source: "Analytic Calculation Notes for the Unique Basin Threshold eps_c(q)" (Nov 5, 2026),
 * Sec. 1-2, and "The separating noise eps_c(q) in the bistable window" (Nov 5, 2026), Sec. 1-3.
 * With s = 1 - q, in the coordinates t = 2u/(1+v), m = (1-v)/(2(1+v)), one round is
 * T(t, m) = (Delta/Sigma, L3/Sigma), where
 *
 *   Sigma = 2 [4 m^2 A + 4 m B + t^2 C + E],   L3 = 4 m^2 E + 4 m B + t^2 F + G,
 *   Delta = t (H + I m)
 *
 *            A = G              B        C                E               F            H             I
 *   Type 3/4 1+e1-e2-e3         1-e1     2e1+e2+e3        1+e1+e2+e3      2e1-e2-e3    16 e3         0
 *   4Q       1-s^3-s^4+s^5      1-s^5    3s^5+s^3         1+s^3+s^4+s^5   s^3-s^5      4s^4(3s+1)    8s^4(1-s)
 *
 * with e1 = s^B', e2 = s^(A'-1), e3 = s^A', (A', B') = (6, 4) for Type 3 and (4, 2) for Type 4.
 * Every coefficient is >= 0 on 0 < s < 1, which makes one-sided evaluation immediate.
 *
 * The certificate is Theorem 2.1 of the calculation notes, evaluated at a test input:
 *
 *   T^n p_eps strictly above S (t_n > t_S and m_n < m_S)  ==>  eps < eps_c, orbit -> P_B
 *   T^n p_eps strictly below S (t_n < t_S and m_n > m_S)  ==>  eps > eps_c, orbit -> Z(q)
 *
 * Exact real-root isolation of t_n(eps) = t_S is out of reach at the depths needed, because the
 * degree in eps grows like 2^n.  So the sign is decided at chosen inputs, and a bisection between
 * the closed-form bounds L_0, U_0 (Eqs. 39-42) shrinks the bracket.  The bisection only chooses
 * which inputs to test; each endpoint of the result carries its own certificate.
 *
 * Enclosure.  Each float64 result is pushed past its rounding error, so the enclosure holds for
 * the exact map.  T preserves the order (Fact 1/2; Eqs. 12-13), so a box [l, u] maps inside
 * [T(l), T(u)]: only the two corners are evaluated, and the box grows with the true dynamics
 * instead of with interval dependency.  This uses the proved order property, nothing numerical.
 *
 * A bracket may stay wide (close to q_bi or q_SN convergence is very slow).  It is still rigorous;
 * it is reported as partial and never tightened by estimation.  No table or interpolation is used.
 */

import type { RegimeProtocol } from "./regimes";

// ---------------------------------------------------------------- directed rounding
const ULP = 2 ** -52;
const TINY = 2 ** -1074;
/** Next float above / below at least; a rounded result pushed this way bounds the exact value. */
const up = (x: number) => x + (Math.abs(x) * ULP + TINY);
const dn = (x: number) => x - (Math.abs(x) * ULP + TINY);

export type Interval = readonly [number, number];

const iv = (x: number): Interval => [dn(x), up(x)];
const add = (a: Interval, b: Interval): Interval => [dn(a[0] + b[0]), up(a[1] + b[1])];
const sub = (a: Interval, b: Interval): Interval => [dn(a[0] - b[1]), up(a[1] - b[0])];
function mul(a: Interval, b: Interval): Interval {
  const c = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]];
  return [dn(Math.min(...c)), up(Math.max(...c))];
}
function div(a: Interval, b: Interval): Interval {
  if (!(b[0] > 0)) return [NaN, NaN]; // only positive denominators occur here
  const c = [a[0] / b[0], a[0] / b[1], a[1] / b[0], a[1] / b[1]];
  return [dn(Math.min(...c)), up(Math.max(...c))];
}
function sqrt(a: Interval): Interval {
  if (!(a[0] > 0)) return [NaN, NaN];
  // ECMAScript does not require Math.sqrt to be correctly rounded: allow one more unit each way
  return [dn(dn(Math.sqrt(a[0]))), up(up(Math.sqrt(a[1])))];
}
function pow(a: Interval, k: number): Interval {
  let out: Interval = [1, 1];
  for (let i = 0; i < k; i++) out = mul(out, a);
  return out;
}
const scale = (k: number, a: Interval): Interval => mul([k, k], a);
/** Coefficients are >= 0 on 0 < s < 1; a rounding-sized dip below 0 is removed. */
const nonneg = (a: Interval): Interval => [Math.max(0, a[0]), Math.max(0, a[1])];

// ---------------------------------------------------------------- the map
export interface MapCoeffs {
  A: Interval; B: Interval; C: Interval; E: Interval; F: Interval; H: Interval; I: Interval;
}

const EXPONENTS = { type3: [6, 4], type4: [4, 2] } as const;

export function mapCoeffs(protocol: RegimeProtocol, q: number): MapCoeffs {
  const s = iv(1 - q);
  const one: Interval = [1, 1];
  if (protocol === "q4") {
    const s3 = pow(s, 3), s4 = pow(s, 4), s5 = pow(s, 5);
    return {
      A: nonneg(add(sub(sub(one, s3), s4), s5)),
      B: nonneg(sub(one, s5)),
      C: add(scale(3, s5), s3),
      E: add(add(add(one, s3), s4), s5),
      F: nonneg(sub(s3, s5)),
      H: scale(4, mul(s4, add(scale(3, s), one))),
      I: nonneg(scale(8, mul(s4, sub(one, s)))),
    };
  }
  const [a, b] = EXPONENTS[protocol];
  const e1 = pow(s, b), e2 = pow(s, a - 1), e3 = pow(s, a);
  return {
    A: nonneg(sub(sub(add(one, e1), e2), e3)),
    B: nonneg(sub(one, e1)),
    C: add(add(scale(2, e1), e2), e3),
    E: add(add(add(one, e1), e2), e3),
    F: nonneg(sub(sub(scale(2, e1), e2), e3)),
    H: scale(16, e3),
    I: [0, 0],
  };
}

/**
 * Enclosure of T at the single point (t, m), t >= 0, m >= 0: [t'lo, t'hi, m'lo, m'hi].
 * Exported for the tests; the certifier below inlines the two one-sided parts it needs.
 */
export function stepEnclosure(c: MapCoeffs, t: number, m: number): [number, number, number, number] {
  const m2H = up(up(m * m) * 4), m2L = dn(dn(m * m) * 4);
  const m4H = up(4 * m), m4L = dn(4 * m);
  const t2H = up(t * t), t2L = dn(t * t);
  const sumH = up(up(up(up(m2H * c.A[1]) + up(m4H * c.B[1])) + up(t2H * c.C[1])) + c.E[1]);
  const sumL = dn(dn(dn(dn(m2L * c.A[0]) + dn(m4L * c.B[0])) + dn(t2L * c.C[0])) + c.E[0]);
  const SH = up(2 * sumH), SL = dn(2 * sumL);
  const DH = up(t * up(c.H[1] + up(c.I[1] * m))), DL = dn(t * dn(c.H[0] + dn(c.I[0] * m)));
  const LH = up(up(up(up(m2H * c.E[1]) + up(m4H * c.B[1])) + up(t2H * c.F[1])) + c.A[1]);
  const LL = dn(dn(dn(dn(m2L * c.E[0]) + dn(m4L * c.B[0])) + dn(t2L * c.F[0])) + c.A[0]);
  return [dn(DL / SH), up(DH / SL), dn(LL / SH), up(LH / SL)];
}

// ---------------------------------------------------------------- fixed points
/** Common form of calculation notes Eq. (1): (alpha, beta, gamma, delta, eta, b). */
function commonForm(protocol: RegimeProtocol, q: number) {
  const s = iv(1 - q);
  const sp = (k: number) => pow(s, k);
  const onePlusS = add([1, 1], s);
  if (protocol === "q4") {
    return {
      alpha: scale(2, sp(5)), beta: mul(sp(4), onePlusS), gamma: scale(2, sp(5)),
      delta: mul(sp(3), onePlusS), eta: add(sp(3), sp(5)), b: sp(5),
    };
  }
  const [a, b] = EXPONENTS[protocol];
  const line = mul(sp(a - 1), onePlusS);
  return {
    alpha: scale(2, sp(a)), beta: scale(2, sp(a)), gamma: line, delta: line,
    eta: scale(2, sp(b)), b: sp(b),
  };
}

export interface FixedPoint {
  t: Interval; m: Interval;
  /** Plane coordinates: transverse u (z for 4Q) and longitudinal v (x for 4Q). */
  u: Interval; v: Interval;
  /** Fidelity, F = (1 + t) / (2 (1 + 2m)) = (1 + v + 2u)/4. */
  F: Interval;
}

function fixedPointAt(p: ReturnType<typeof commonForm>, v: Interval): FixedPoint {
  const d = sub(add(mul(p.alpha, v), p.beta), p.delta);
  const u = sqrt(div(mul(v, d), p.gamma)); // Eq. 21
  const onePlusV = add([1, 1], v);
  return {
    u, v,
    t: div(scale(2, u), onePlusV), // Eq. 22
    m: div(sub([1, 1], v), scale(2, onePlusV)),
    F: mul(add(add([1, 1], v), scale(2, u)), [0.25, 0.25]),
  };
}

/**
 * Rigorous enclosure of the Bell-connected fixed point P_B(q): the larger root v_B of H_s and the
 * positive square root for u (Nov 6 note, Eq. 12).  The same branch runs over 0 < q < q_SN; the
 * caller decides whether q is in that range.  null if the enclosure fails (discriminant or u^2 not
 * certified positive, e.g. right at the saddle-node).
 */
export function bellFixedPoint(protocol: RegimeProtocol, q: number): FixedPoint | null {
  const p = commonForm(protocol, q);
  const h2 = add(mul(p.gamma, p.b), mul(p.eta, p.alpha));
  const h1 = sub(mul(p.eta, sub(p.beta, p.delta)), mul(p.gamma, p.alpha));
  const h0 = mul(p.gamma, sub([1, 1], p.beta));
  const root = sqrt(sub(mul(h1, h1), scale(4, mul(h2, h0))));
  if (!(root[0] > 0) || !(h2[0] > 0)) return null;
  // (-h1 + root) / (2 h2); when h1 > 0 the two terms cancel, and then h0 < 0 and the equal form
  // 2(-h0) / (h1 + root) is used instead
  const vB = h1[0] > 0
    ? div(scale(2, [-h0[1], -h0[0]]), add(h1, root))
    : div(add([-h1[1], -h1[0]], root), scale(2, h2));
  const PB = fixedPointAt(p, vB);
  return [PB.u[0], PB.u[1], PB.F[0], PB.F[1]].every(Number.isFinite) && vB[0] > 0 ? PB : null;
}

/**
 * Rigorous enclosures of the saddle S and the sink P_B (Eqs. 16-22).  The caller is responsible
 * for q being inside the bistable window; null if the enclosures do not separate the two roots.
 */
export function bistableFixedPoints(
  protocol: RegimeProtocol, q: number,
): { S: FixedPoint; PB: FixedPoint } | null {
  const p = commonForm(protocol, q);
  const h2 = add(mul(p.gamma, p.b), mul(p.eta, p.alpha));
  const h1 = sub(mul(p.eta, sub(p.beta, p.delta)), mul(p.gamma, p.alpha));
  const h0 = mul(p.gamma, sub([1, 1], p.beta));
  const root = sqrt(sub(mul(h1, h1), scale(4, mul(h2, h0))));
  const minusH1: Interval = [-h1[1], -h1[0]];
  if (!(root[0] > 0) || !(minusH1[0] > 0)) return null;
  // conjugate form for the small root: the two terms do not cancel
  const vS = div(scale(2, h0), add(minusH1, root));
  const vB = div(add(minusH1, root), scale(2, h2));
  const S = fixedPointAt(p, vS);
  const PB = fixedPointAt(p, vB);
  const ok = (x: FixedPoint) => [x.t[0], x.t[1], x.m[0], x.m[1]].every(Number.isFinite);
  return ok(S) && ok(PB) ? { S, PB } : null;
}

// ---------------------------------------------------------------- side of one input
export const SIDE_PB = 1;
export const SIDE_Z = -1;
export const SIDE_UNDECIDED = 0;
export type Side = 1 | -1 | 0;

/** The enclosure is abandoned once this wide: S lies inside it and it only widens from there. */
const BOX_LIMIT = 1e-3;

export interface Classification {
  side: Side;
  /** Round at which the side was certified, or at which the search stopped. */
  rounds: number;
}

/**
 * Certify on which side of S the orbit of the isotropic input ends up, for this exact float eps
 * in (0, 1).  side = SIDE_UNDECIDED means only that the enclosure could not tell.
 */
export function classifyInput(
  c: MapCoeffs, S: Pick<FixedPoint, "t" | "m">, eps: number, maxRounds: number,
): Classification {
  const [cAl, cAh] = c.A, [cBl, cBh] = c.B, [cCl, cCh] = c.C, [cEl, cEh] = c.E;
  const [cFl, cFh] = c.F, [cHl, cHh] = c.H, [cIl, cIh] = c.I;
  const [tSl, tSh] = S.t, [mSl, mSh] = S.m;

  // Eq. 8: t = 2(1-eps)/(2-eps), m = eps/(2(2-eps)).
  const dL = dn(2 - eps), dH = up(2 - eps);
  // lower corner l = (tl, mh); upper corner u = (th, ml)
  let tl = dn(dn(2 * dn(1 - eps)) / dH);
  let th = up(up(2 * up(1 - eps)) / dL);
  let ml = dn(eps / up(2 * dH));
  let mh = up(eps / dn(2 * dL));

  for (let n = 0; n <= maxRounds; n++) {
    if (tl > tSh && mh < mSl) return { side: SIDE_PB, rounds: n };
    if (th < tSl && ml > mSh) return { side: SIDE_Z, rounds: n };
    if (th - tl > BOX_LIMIT || mh - ml > BOX_LIMIT) return { side: SIDE_UNDECIDED, rounds: n };
    if (n === maxRounds) break;

    // T(l): lower bound on t', upper bound on m'.  T(u): upper bound on t', lower bound on m'.
    // All terms are >= 0, so each bound is one chain of same-direction roundings.
    const a2h = up(mh * mh) * 4, a2l = dn(mh * mh) * 4; // m^2 at corner l
    const a1h = up(4 * mh), a1l = dn(4 * mh);
    const atH = up(tl * tl), atL = dn(tl * tl);
    const SlH = up(2 * up(up(up(up(up(a2h) * cAh) + up(a1h * cBh)) + up(atH * cCh)) + cEh));
    const SlL = dn(2 * dn(dn(dn(dn(dn(a2l) * cAl) + dn(a1l * cBl)) + dn(atL * cCl)) + cEl));
    const DlL = dn(tl * dn(cHl + dn(cIl * mh)));
    const LlH = up(up(up(up(up(a2h) * cEh) + up(a1h * cBh)) + up(atH * cFh)) + cAh);

    const b2h = up(ml * ml) * 4, b2l = dn(ml * ml) * 4; // m^2 at corner u
    const b1h = up(4 * ml), b1l = dn(4 * ml);
    const btH = up(th * th), btL = dn(th * th);
    const SuH = up(2 * up(up(up(up(up(b2h) * cAh) + up(b1h * cBh)) + up(btH * cCh)) + cEh));
    const SuL = dn(2 * dn(dn(dn(dn(dn(b2l) * cAl) + dn(b1l * cBl)) + dn(btL * cCl)) + cEl));
    const DuH = up(th * up(cHh + up(cIh * ml)));
    const LuL = dn(dn(dn(dn(dn(b2l) * cEl) + dn(b1l * cBl)) + dn(btL * cFl)) + cAl);

    // T(R) is inside R: clamping cannot lose the exact image.
    const tl2 = Math.max(0, dn(DlL / SlH));
    const th2 = Math.min(1, up(DuH / SuL));
    const mh2 = Math.min(0.5, up(LlH / SlL));
    const ml2 = Math.max(0, dn(LuL / SuH));
    tl = tl2; th = th2; ml = ml2; mh = mh2;
  }
  return { side: SIDE_UNDECIDED, rounds: maxRounds };
}

// ---------------------------------------------------------------- bracket
/** Bracket width at which the threshold counts as determined (the research tables: 2e-9). */
export const EPSC_TARGET_WIDTH = 2e-9;
/** Rounds allowed per tested input, and over one whole bracket computation. */
export const EPSC_MAX_ROUNDS = 2_000_000;
export const EPSC_TOTAL_ROUNDS = 40_000_000;

export interface EpscBracket {
  protocol: RegimeProtocol;
  q: number;
  /** Certified: lo is in the P_B basin and hi is in the Z basin, so lo < eps_c(q) < hi. */
  lo: number;
  hi: number;
  /** Rounds at which the sides of lo and of hi were certified. */
  roundsLo: number;
  roundsHi: number;
  /** Rounds over the whole computation. */
  totalRounds: number;
  /** False while the computation is still running. */
  done: boolean;
  /** Why the bisection stopped once done. */
  stop: "" | "resolution" | "budget" | "failed";
}

export const isDetermined = (b: Pick<EpscBracket, "lo" | "hi">) => b.hi - b.lo <= EPSC_TARGET_WIDTH;

/**
 * Yields a valid certified bracket after each tested input and returns the final one.  Yields
 * nothing and returns null if the fixed points cannot be enclosed.  Bracket widths never grow.
 */
export function* epscBrackets(
  protocol: RegimeProtocol, q: number,
  maxRounds = EPSC_MAX_ROUNDS, totalBudget = EPSC_TOTAL_ROUNDS,
): Generator<EpscBracket, EpscBracket | null, void> {
  const fp = bistableFixedPoints(protocol, q);
  if (!fp) return null;
  const c = mapCoeffs(protocol, q);
  const { t: tS, m: mS } = fp.S;

  // Eqs. 39-42, with a margin so that the n = 0 test is strict.
  const a0 = div(scale(2, sub([1, 1], tS)), sub([2, 2], tS));
  const b0 = div(scale(4, mS), add([1, 1], scale(2, mS)));
  let lo = dn(dn(Math.min(a0[0], b0[0])));
  let hi = up(up(Math.max(a0[1], b0[1])));
  if (!(lo > 0 && hi < 1)) return null;
  const first = classifyInput(c, fp.S, lo, 0);
  const last = classifyInput(c, fp.S, hi, 0);
  if (first.side !== SIDE_PB || last.side !== SIDE_Z) return null;

  let roundsLo = 0, roundsHi = 0, total = 0;
  const snap = (done: boolean, stop: EpscBracket["stop"]): EpscBracket => ({
    protocol, q, lo, hi, roundsLo, roundsHi, totalRounds: total, done, stop,
  });
  yield snap(false, "");

  let stop: EpscBracket["stop"] = "resolution";
  for (;;) {
    if (total >= totalBudget) { stop = "budget"; break; }
    const w = hi - lo;
    const mid = lo + w / 2;
    if (mid <= lo || mid >= hi) break;
    let r = classifyInput(c, fp.S, mid, maxRounds);
    total += r.rounds;
    let at = mid;
    if (r.side === SIDE_UNDECIDED) {
      // At the resolution limit the middle may be undecidable while a nearby input is not.
      const budgetHit = r.rounds >= maxRounds;
      if (budgetHit) { stop = "budget"; break; }
      let found = false;
      for (const x of [lo + w / 4, lo + (3 * w) / 4]) {
        r = classifyInput(c, fp.S, x, maxRounds);
        total += r.rounds;
        if (r.side !== SIDE_UNDECIDED) { at = x; found = true; break; }
      }
      if (!found) break;
    }
    if (r.side === SIDE_PB) { lo = at; roundsLo = r.rounds; } else { hi = at; roundsHi = r.rounds; }
    yield snap(false, "");
  }
  return snap(true, stop);
}

/** Run the whole computation synchronously (tests, tools). */
export function computeEpsc(protocol: RegimeProtocol, q: number): EpscBracket | null {
  const g = epscBrackets(protocol, q);
  for (;;) {
    const r = g.next();
    if (r.done) return r.value;
  }
}

// ---------------------------------------------------------------- verdict for one eps
export type BasinVerdict = "P_B" | "Z" | "undetermined";

/**
 * eps <= lo is in the P_B basin and eps >= hi in the Z basin: the endpoints themselves carry the
 * certificate, and the basins are the intervals (0, eps_c) and (eps_c, 1).  Between them nothing
 * is claimed, even if eps looks equal to eps_c in floating point.  Requires 0 < eps < 1.
 */
export function basinVerdict(b: Pick<EpscBracket, "lo" | "hi">, eps: number): BasinVerdict {
  if (eps <= b.lo) return "P_B";
  if (eps >= b.hi) return "Z";
  return "undetermined";
}

/** Initial Bell fidelity of (1-eps) |Phi+><Phi+| + eps I/4. */
export const initialFidelity = (eps: number) => 1 - 0.75 * eps;
