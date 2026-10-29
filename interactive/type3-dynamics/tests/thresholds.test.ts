import { describe, expect, it } from "vitest";
import { allocTrajectory, computeTrajectory } from "../src/type3map";
import { fAbs, fRel, qAbs, qRel } from "../src/thresholds";

// Table 1 of the q_rel note (q_abs column also matches Table 1 of the q_abs note).
const TABLE: [number, number, number][] = [
  [0.02, 0.005659445995, 0.005675467891],
  [0.05, 0.013937335868, 0.014034550412],
  [0.10, 0.027141672179, 0.027510478464],
  [0.15, 0.039562631537, 0.040346258209],
  [0.20, 0.051155686268, 0.052465899458],
  [0.30, 0.071710593708, 0.074289724421],
  [0.40, 0.088567511457, 0.092530582050],
  [0.50, 0.101555797133, 0.106859191408],
  [0.60, 0.110528200570, 0.117024018264],
  [0.65, 0.113454759292, 0.120470156477],
];

const EPS_DENSE = Array.from({ length: 65 }, (_, i) => 0.01 * (i + 1)); // 0.01 … 0.65

function run(eps: number, q: number, N: number) {
  return computeTrajectory(allocTrajectory(N), eps, q, N);
}

describe("Test T1 — analytic values reproduce the calculation notes", () => {
  it("matches Table 1 to the printed 12 digits", () => {
    for (const [eps, a, r] of TABLE) {
      expect(Math.abs(qAbs(eps) - a)).toBeLessThan(1e-12);
      expect(Math.abs(qRel(eps) - r)).toBeLessThan(1e-12);
    }
  });
  it("matches the eps -> 2/3 limits", () => {
    expect(Math.abs(qAbs(2 / 3 - 1e-12) - 0.114192361209)).toBeLessThan(1e-10);
    expect(Math.abs(qRel(2 / 3 - 1e-12) - 0.121368048928)).toBeLessThan(1e-10);
  });
  it("returns roots of the boundary polynomials", () => {
    for (const eps of EPS_DENSE) {
      const t = 1 - eps;
      expect(Math.abs(fAbs(1 - qAbs(eps), t))).toBeLessThan(1e-14);
      expect(Math.abs(fRel(1 - qRel(eps), t))).toBeLessThan(1e-14);
    }
  });
});

describe("Test T2 — structure", () => {
  it("is NaN outside the entangled regime 0 < eps < 2/3", () => {
    for (const eps of [0, -0.1, 2 / 3, 0.7, 1]) {
      expect(qAbs(eps)).toBeNaN();
      expect(qRel(eps)).toBeNaN();
    }
  });
  it("satisfies 0 < q_abs < q_rel, both increasing in eps", () => {
    let pa = 0, pr = 0;
    for (const eps of EPS_DENSE) {
      const a = qAbs(eps), r = qRel(eps);
      expect(a).toBeGreaterThan(0);
      expect(a).toBeLessThan(r);
      expect(a).toBeGreaterThan(pa);
      expect(r).toBeGreaterThan(pr);
      pa = a; pr = r;
    }
  });
});

describe("Test T3 — thresholds agree with the exact map (independent of the notes' proof)", () => {
  const N = 200;
  const DQ = 1e-6;
  const TOL = 1e-15; // once converged, successive values agree to float64 resolution
  /** Strict decrease for the first rounds, then no rise beyond rounding, over all N. */
  const decreasing = (x: (n: number) => number) => {
    for (let n = 0; n < 3; n++) expect(x(n + 1)).toBeLessThan(x(n));
    for (let n = 0; n < N; n++) expect(x(n + 1) - x(n)).toBeLessThan(TOL);
  };
  it("q_abs: below -> every non-target population decreases each round; above -> Phi- rises at round 1", () => {
    for (const eps of EPS_DENSE) {
      const qa = qAbs(eps);
      const lo = run(eps, qa - DQ, N);
      decreasing((n) => lo.pPhiM[n]);
      decreasing((n) => lo.pPsiP[n]);
      const hi = run(eps, qa + DQ, 1);
      expect(hi.pPhiM[1]).toBeGreaterThan(hi.pPhiM[0]);
    }
  });
  it("q_rel: below -> every ratio p_i/p_Phi+ decreases each round; above -> Phi- ratio rises at round 1", () => {
    for (const eps of EPS_DENSE) {
      const qr = qRel(eps);
      const lo = run(eps, qr - DQ, N);
      decreasing((n) => lo.pPhiM[n] / lo.pPhiP[n]);
      decreasing((n) => lo.pPsiP[n] / lo.pPhiP[n]);
      const hi = run(eps, qr + DQ, 1);
      expect(hi.pPhiM[1] / hi.pPhiP[1]).toBeGreaterThan(hi.pPhiM[0] / hi.pPhiP[0]);
    }
  });
});
