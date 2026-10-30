import { describe, expect, it } from "vitest";
import { allocTrajectory } from "../src/type3map";
import { computeTrajectory4Q } from "../src/q4map";
import { fAbs4Q, fRel4Q, qAbs4Q, qRel4Q } from "../src/thresholds4q";
import { qAbs, qRel } from "../src/thresholds";
import { qAbs4, qRel4 } from "../src/thresholds4";

// Table 1 of the 4Q q_abs note.
const TABLE_ABS: [number, number][] = [
  [0.02, 0.009707058830], [0.05, 0.023231092240], [0.10, 0.043297940301],
  [0.15, 0.060665417272], [0.20, 0.075697361355], [0.30, 0.099844658982],
  [0.40, 0.117399578335], [0.50, 0.129400238963], [0.60, 0.136455003648],
  [0.65, 0.138225184013],
];
// Sec. 17 of the 4Q q_rel note.
const TABLE_REL: [number, number][] = [
  [0.01, 0.004938087562], [0.02, 0.009754651976], [0.05, 0.023507999519],
  [0.10, 0.044286942985], [0.15, 0.062666225863], [0.20, 0.078916032920],
  [0.30, 0.105877025226], [0.40, 0.126513294424], [0.50, 0.141710916436],
  [0.60, 0.152003414201], [2 / 3, 0.156264485085], [0.75, 0.158656604525],
  [0.90, 0.154063967757],
];

const EPS_ABS = Array.from({ length: 65 }, (_, i) => 0.01 * (i + 1)); //  0.01 … 0.65
const EPS_REL = Array.from({ length: 99 }, (_, i) => 0.01 * (i + 1)); //  0.01 … 0.99

function run(eps: number, q: number, N: number) {
  return computeTrajectory4Q(allocTrajectory(N), eps, q, N);
}

describe("Test T1-4Q — analytic 4Q values reproduce the calculation notes", () => {
  it("matches both tables to the printed 12 digits", () => {
    for (const [eps, v] of TABLE_ABS) expect(Math.abs(qAbs4Q(eps) - v)).toBeLessThan(1e-12);
    for (const [eps, v] of TABLE_REL) expect(Math.abs(qRel4Q(eps) - v)).toBeLessThan(1e-12);
  });
  it("matches the limits, the q_rel maximum and the small-eps expansion quoted in the text", () => {
    expect(Math.abs(qAbs4Q(2 / 3 - 1e-12) - 0.138558310224)).toBeLessThan(1e-10); // abs Eq. 171
    expect(Math.abs(qRel4Q(1 - 1e-12) - 0.143325116145)).toBeLessThan(1e-10); //     rel Sec. 20.2
    // rel Eqs. 113-118: interior maximum at eps = (3 + r^2)/(2 + 4 r^2)
    const ePk = 0.767499473760;
    expect(Math.abs(qRel4Q(ePk) - 0.158732684144)).toBeLessThan(1e-11);
    expect(qRel4Q(ePk)).toBeGreaterThan(qRel4Q(ePk - 1e-3));
    expect(qRel4Q(ePk)).toBeGreaterThan(qRel4Q(ePk + 1e-3));
    const r = 1 - qRel4Q(ePk);
    expect(Math.abs((3 + r * r) / (2 + 4 * r * r) - ePk)).toBeLessThan(1e-11);
    // rel Eq. 122: q_rel = eps/2 - 5 eps^2/8 + 19 eps^3/32 + O(eps^4)
    for (const e of [1e-3, 3e-3, 1e-2]) {
      expect(Math.abs(qRel4Q(e) - (e / 2 - (5 * e * e) / 8 + (19 * e ** 3) / 32))).toBeLessThan(10 * e ** 4);
    }
    expect(Math.abs(qAbs4Q(0.15) - 0.060665417272)).toBeLessThan(1e-12); // abs Eq. 167
  });
  it("returns roots of the boundary polynomials, also in the q-form of rel Eq. 108", () => {
    for (const eps of EPS_ABS) expect(Math.abs(fAbs4Q(1 - qAbs4Q(eps), 1 - eps))).toBeLessThan(1e-14);
    for (const e of EPS_REL) {
      const q = qRel4Q(e);
      expect(Math.abs(fRel4Q(1 - q, 1 - e))).toBeLessThan(1e-14);
      const poly = (2 + e - 2 * e * e) * q ** 5 + (10 * e * e - 5 * e - 11) * q ** 4 +
        (-21 * e * e + 13 * e + 22) * q ** 3 + (23 * e * e - 19 * e - 20) * q ** 2 +
        (-13 * e * e + 14 * e + 8) * q + 3 * e * e - 4 * e;
      expect(Math.abs(poly)).toBeLessThan(1e-14);
    }
  });
});

describe("Test T2-4Q — structure", () => {
  it("q_abs is NaN outside 0 < eps < 2/3; q_rel is defined on the whole interior 0 < eps < 1", () => {
    for (const eps of [0, -0.1, 2 / 3, 0.7, 0.9, 1]) expect(qAbs4Q(eps)).toBeNaN();
    for (const eps of [0, -0.1, 1, 1.1]) expect(qRel4Q(eps)).toBeNaN();
    for (const eps of [2 / 3, 0.7, 0.9, 0.999]) expect(Number.isFinite(qRel4Q(eps))).toBe(true);
  });
  it("0 < q_abs < q_rel; q_abs increasing in eps with r* > 6/7 (abs Eqs. 84, 165)", () => {
    let pa = 0;
    for (const eps of EPS_ABS) {
      const a = qAbs4Q(eps), r = qRel4Q(eps);
      expect(a).toBeGreaterThan(0);
      expect(a).toBeLessThan(r);
      expect(a).toBeGreaterThan(pa);
      expect(1 - a).toBeGreaterThan(6 / 7);
      pa = a;
    }
  });
  it("q_rel is not monotone on (0, 1): it rises to its maximum and then falls (rel Sec. 17)", () => {
    expect(qRel4Q(0.9)).toBeLessThan(qRel4Q(0.75));
    expect(qRel4Q(0.75)).toBeGreaterThan(qRel4Q(0.6));
  });
  it("is not the Type-3 / Type-4 threshold: above Type 3 throughout, crossing Type 4 once", () => {
    for (const eps of EPS_ABS) {
      expect(qAbs4Q(eps)).toBeGreaterThan(qAbs(eps));
      expect(qRel4Q(eps)).toBeGreaterThan(qRel(eps));
      // 4Q vs Type 4: larger for small eps, smaller beyond eps ~ 0.451 (abs) / 0.483 (rel)
      expect(qAbs4Q(eps) > qAbs4(eps)).toBe(eps < 0.451);
      expect(qRel4Q(eps) > qRel4(eps)).toBe(eps < 0.483);
    }
  });
});

describe("Test T3-4Q — thresholds agree with the exact 4Q map (independent of the notes' proof)", () => {
  const N = 200;
  const DQ = 1e-6;
  const TOL = 1e-15; // once converged, successive values agree to float64 resolution
  const decreasing = (x: (n: number) => number) => {
    for (let n = 0; n < 3; n++) expect(x(n + 1)).toBeLessThan(x(n));
    for (let n = 0; n < N; n++) expect(x(n + 1) - x(n)).toBeLessThan(TOL);
  };
  it("q_abs: below -> every non-target population decreases each round; above -> Psi+ rises at round 1", () => {
    for (const eps of EPS_ABS) {
      const qa = qAbs4Q(eps);
      const lo = run(eps, qa - DQ, N);
      decreasing((n) => lo.pPhiM[n]);
      decreasing((n) => lo.pPsiP[n]);
      const hi = run(eps, qa + DQ, 1);
      expect(hi.pPsiP[1]).toBeGreaterThan(hi.pPsiP[0]);
      expect(hi.pPhiM[1]).toBeLessThan(hi.pPhiM[0]); // Psi+ is the limiting component
    }
  });
  it("q_rel: below -> every ratio p_i/p_Phi+ decreases each round; above -> Psi+ ratio rises (also eps > 2/3)", () => {
    for (const eps of EPS_REL) {
      const qr = qRel4Q(eps);
      const lo = run(eps, qr - DQ, N);
      decreasing((n) => lo.pPhiM[n] / lo.pPhiP[n]);
      decreasing((n) => lo.pPsiP[n] / lo.pPhiP[n]);
      const hi = run(eps, qr + DQ, 1);
      expect(hi.pPsiP[1] / hi.pPhiP[1]).toBeGreaterThan(hi.pPsiP[0] / hi.pPhiP[0]);
      expect(hi.pPhiM[1] / hi.pPhiP[1]).toBeLessThan(hi.pPhiM[0] / hi.pPhiP[0]);
    }
  });
});
