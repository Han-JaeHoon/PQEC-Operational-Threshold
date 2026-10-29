import { describe, expect, it } from "vitest";
import { allocTrajectory } from "../src/type3map";
import { computeTrajectory4 } from "../src/type4map";
import { fAbs4, fRel4, qAbs4, qRel4 } from "../src/thresholds4";
import { qAbs, qRel } from "../src/thresholds";

// Table 1 of the Type-4 q_rel note (q_abs column also matches Table 2 of the Type-4 q_abs note).
const TABLE: [number, number, number][] = [
  [0.02, 0.006632618362, 0.006654687458],
  [0.05, 0.016447044318, 0.016583423004],
  [0.10, 0.032408164382, 0.032942058135],
  [0.15, 0.047810129835, 0.048981686135],
  [0.20, 0.062575889256, 0.064599667015],
  [0.30, 0.089875575401, 0.094125003670],
  [0.40, 0.113641041067, 0.120578868012],
  [0.50, 0.133188123130, 0.142985185946],
  [0.60, 0.147822494728, 0.160376167384],
  [0.65, 0.153074008453, 0.166893343399],
];

const EPS_DENSE = Array.from({ length: 65 }, (_, i) => 0.01 * (i + 1)); // 0.01 … 0.65

function run(eps: number, q: number, N: number) {
  return computeTrajectory4(allocTrajectory(N), eps, q, N);
}

describe("Test T1-4 — analytic Type-4 values reproduce the calculation notes", () => {
  it("matches the tables to the printed 12 digits", () => {
    for (const [eps, a, r] of TABLE) {
      expect(Math.abs(qAbs4(eps) - a)).toBeLessThan(1e-12);
      expect(Math.abs(qRel4(eps) - r)).toBeLessThan(1e-12);
    }
  });
  it("matches the eps -> 2/3 limits and the eps = 0.15 values quoted in the text", () => {
    expect(Math.abs(qAbs4(2 / 3 - 1e-12) - 0.154495318610)).toBeLessThan(1e-10);
    expect(Math.abs(qRel4(2 / 3 - 1e-12) - 0.168713457068)).toBeLessThan(1e-10);
    expect(Math.abs(qRel4(0.15) - 0.0489816861346)).toBeLessThan(1e-13);
    expect(Math.abs(qRel4(0.15) - qAbs4(0.15) - 1.1715563e-3)).toBeLessThan(1e-10);
  });
  it("returns roots of the boundary polynomials", () => {
    for (const eps of EPS_DENSE) {
      const t = 1 - eps;
      expect(Math.abs(fAbs4(1 - qAbs4(eps), t))).toBeLessThan(1e-14);
      expect(Math.abs(fRel4(1 - qRel4(eps), t))).toBeLessThan(1e-14);
    }
  });
});

describe("Test T2-4 — structure", () => {
  it("is NaN outside the entangled regime 0 < eps < 2/3", () => {
    for (const eps of [0, -0.1, 2 / 3, 0.7, 1]) {
      expect(qAbs4(eps)).toBeNaN();
      expect(qRel4(eps)).toBeNaN();
    }
  });
  it("satisfies 0 < q_abs < q_rel, both increasing in eps; r* above the notes' bounds 5/6, 4/5", () => {
    let pa = 0, pr = 0;
    for (const eps of EPS_DENSE) {
      const a = qAbs4(eps), r = qRel4(eps);
      expect(a).toBeGreaterThan(0);
      expect(a).toBeLessThan(r);
      expect(a).toBeGreaterThan(pa);
      expect(r).toBeGreaterThan(pr);
      expect(1 - a).toBeGreaterThan(5 / 6);
      expect(1 - r).toBeGreaterThan(4 / 5);
      pa = a; pr = r;
    }
  });
  it("is not the Type-3 threshold (each protocol has its own)", () => {
    for (const eps of EPS_DENSE) {
      expect(qAbs4(eps)).toBeGreaterThan(qAbs(eps));
      expect(qRel4(eps)).toBeGreaterThan(qRel(eps));
    }
  });
});

describe("Test T3-4 — thresholds agree with the exact Type-4 map (independent of the notes' proof)", () => {
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
      const qa = qAbs4(eps);
      const lo = run(eps, qa - DQ, N);
      decreasing((n) => lo.pPhiM[n]);
      decreasing((n) => lo.pPsiP[n]);
      const hi = run(eps, qa + DQ, 1);
      expect(hi.pPhiM[1]).toBeGreaterThan(hi.pPhiM[0]);
      expect(hi.pPsiP[1]).toBeLessThan(hi.pPsiP[0]); // Phi- is the limiting component
    }
  });
  it("q_rel: below -> every ratio p_i/p_Phi+ decreases each round; above -> Phi- ratio rises at round 1", () => {
    for (const eps of EPS_DENSE) {
      const qr = qRel4(eps);
      const lo = run(eps, qr - DQ, N);
      decreasing((n) => lo.pPhiM[n] / lo.pPhiP[n]);
      decreasing((n) => lo.pPsiP[n] / lo.pPhiP[n]);
      const hi = run(eps, qr + DQ, 1);
      expect(hi.pPhiM[1] / hi.pPhiP[1]).toBeGreaterThan(hi.pPhiM[0] / hi.pPhiP[0]);
      expect(hi.pPsiP[1] / hi.pPhiP[1]).toBeLessThan(hi.pPsiP[0] / hi.pPhiP[0]);
    }
  });
});
