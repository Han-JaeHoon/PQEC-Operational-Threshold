import { describe, expect, it } from "vitest";
import {
  allocTrajectory, bellPopulations, computeTrajectory, Step, summarize, type3Step,
} from "../src/type3map";

const EPS_GRID = Array.from({ length: 101 }, (_, i) => i / 100);
const Q_GRID = [0, 1e-4, 0.001, 0.01, 0.0327, 0.05, 0.1, 0.13, 0.2, 0.5, 0.9, 1];

function run(eps: number, q: number, N: number) {
  return computeTrajectory(allocTrajectory(N), eps, q, N);
}

describe("Test A — initial state", () => {
  it("u0 = v0 = 1 - eps and F0 = 1 - 3 eps / 4 for every eps", () => {
    for (const eps of [...EPS_GRID, 0.0327, 2 / 3, 1e-9]) {
      const t = run(eps, 0.05, 1);
      expect(t.u[0]).toBe(1 - eps);
      expect(t.v[0]).toBe(1 - eps);
      expect(t.pPhiP[0]).toBeCloseTo(1 - (3 * eps) / 4, 15);
      expect(t.pPhiM[0]).toBeCloseTo(eps / 4, 15);
      expect(t.pPsiP[0]).toBeCloseTo(eps / 4, 15);
      expect(t.pPsiM[0]).toBeCloseTo(eps / 4, 15);
    }
  });
});

describe("Tests B, C, D — normalisation, positivity, invariant plane", () => {
  it("holds for every round on an (eps, q) grid, N = 500", () => {
    let worstSum = 0, worstNeg = 0;
    for (const eps of EPS_GRID.filter((_, i) => i % 5 === 0)) {
      for (const q of Q_GRID) {
        const t = run(eps, q, 500);
        for (let n = 0; n <= 501; n++) {
          const s = t.pPhiP[n] + t.pPhiM[n] + t.pPsiP[n] + t.pPsiM[n];
          worstSum = Math.max(worstSum, Math.abs(s - 1));
          worstNeg = Math.min(worstNeg, t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]);
          expect(t.pPsiP[n]).toBe(t.pPsiM[n]); // D: exact equality
          expect(Number.isFinite(t.u[n]) && Number.isFinite(t.v[n])).toBe(true);
        }
      }
    }
    expect(worstSum).toBeLessThan(1e-15); // B
    expect(worstNeg).toBeGreaterThanOrEqual(-1e-15); // C
  });
});

describe("Test E — ideal case q = 0", () => {
  it("reproduces rho -> rho^2 / Tr rho^2 (Bell populations p -> p^2 / sum p^2)", () => {
    for (const eps of [0.05, 0.1, 0.3, 0.5, 2 / 3, 0.8, 0.95]) {
      const t = run(eps, 0, 40);
      let p = [1 - (3 * eps) / 4, eps / 4, eps / 4, eps / 4];
      for (let n = 0; n <= 40; n++) {
        const got = [t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]];
        for (let k = 0; k < 4; k++) expect(Math.abs(got[k] - p[k])).toBeLessThan(1e-14);
        const s = p.reduce((a, x) => a + x * x, 0);
        p = p.map((x) => (x * x) / s);
      }
    }
  });
  it("one ideal round matches F1 = (1+3t)^2 / (4(1+3t^2)) (Step 1 closed form)", () => {
    for (const eps of EPS_GRID) {
      const t = 1 - eps;
      const F1 = ((1 + 3 * t) ** 2) / (4 * (1 + 3 * t * t));
      expect(Math.abs(run(eps, 0, 1).pPhiP[1] - F1)).toBeLessThan(1e-15);
    }
  });
});

describe("Test F — maximally mixed input", () => {
  it("eps = 1 keeps F_n = 1/4 exactly for every q", () => {
    for (const q of Q_GRID) {
      const t = run(1, q, 500);
      for (let n = 0; n <= 501; n++) expect(t.pPhiP[n]).toBe(0.25);
    }
  });
});

describe("one-round closed form already in the repository (pqec_cnot_threshold.F_dec)", () => {
  it("F1 = 1/4 [1 + t(1+t)(v^5 + 5v^6)/(1 + 3 v^4 t^2)], v = 1 - q", () => {
    let worst = 0;
    for (const eps of EPS_GRID) for (const q of Q_GRID) {
      const t = 1 - eps, v = 1 - q;
      const F1 = 0.25 * (1 + (t * (1 + t) * (v ** 5 + 5 * v ** 6)) / (1 + 3 * v ** 4 * t * t));
      worst = Math.max(worst, Math.abs(run(eps, q, 1).pPhiP[1] - F1));
    }
    expect(worst).toBeLessThan(2e-15);
  });
});

describe("fixed point vs. the closed form of note 01 (independent formula)", () => {
  // A = r^3 + r^2 - r + 1, Delta = r^4 A^2 + (5r+1)(1+r)(2r^6 - 1),
  // v* = [r^2 A + sqrt(Delta)] / [r^2 (5r+1)],  u* = sqrt(v*(r - 1 + 2 r v*)/(1+r))
  const fixed = (q: number) => {
    const r = 1 - q;
    const A = r ** 3 + r ** 2 - r + 1;
    const Dl = r ** 4 * A * A + (5 * r + 1) * (1 + r) * (2 * r ** 6 - 1);
    const v = (r * r * A + Math.sqrt(Dl)) / (r * r * (5 * r + 1));
    const u = Math.sqrt((v * (r - 1 + 2 * r * v)) / (1 + r));
    return { u, v, F: (1 + 2 * u + v) / 4 };
  };
  it("iterates converge to (u*, v*) and match the tabulated values", () => {
    const table: [number, number, number, number][] = [
      [0.01, 0.969470, 0.974443, 0.978346],
      [0.05, 0.833597, 0.857754, 0.881237],
      [0.10, 0.602256, 0.647161, 0.712918],
      [0.12, 0.450744, 0.501193, 0.600670],
    ];
    for (const [q, us, vs, Fs] of table) {
      const fp = fixed(q);
      expect(Math.abs(fp.u - us)).toBeLessThan(6e-7);
      expect(Math.abs(fp.v - vs)).toBeLessThan(6e-7);
      expect(Math.abs(fp.F - Fs)).toBeLessThan(6e-7);
      const t = run(0.1, q, 500);
      expect(Math.abs(t.u[500] - fp.u)).toBeLessThan(1e-12);
      expect(Math.abs(t.v[500] - fp.v)).toBeLessThan(1e-12);
      const [u1, v1] = type3Step(fp.u, fp.v, q);
      expect(Math.abs(u1 - fp.u) + Math.abs(v1 - fp.v)).toBeLessThan(1e-14);
    }
  });
});

describe("fidelity-drift cross-check", () => {
  it("drift numerator N/(4D) agrees with F_{n+1} - F_n within the combined error bounds", () => {
    for (const eps of EPS_GRID.filter((_, i) => i % 4 === 0)) for (const q of Q_GRID) {
      const t = run(eps, q, 300);
      for (let n = 0; n <= 300; n++) {
        const tol = t.dFerr[n] + t.dFmapErr[n];
        expect(Math.abs(t.dF[n] - t.dFmap[n])).toBeLessThanOrEqual(tol);
      }
    }
  });
  it("never claims a sign for a step at the fixed point (dF = 0 exactly)", () => {
    const t = run(1, 0.05, 50);
    for (let n = 0; n <= 50; n++) expect(t.step[n]).toBe(Step.Unresolved);
  });
});

describe("Test G (part) — eps = 0.1, q = 0.0327: rise, then later-round decrease", () => {
  it("first step increases, a later step is a resolved decrease", () => {
    const t = run(0.1, 0.0327, 30);
    expect(t.step[0]).toBe(Step.Increase);
    const s = summarize(t);
    expect(s.firstDecrease).toBeGreaterThan(0);
    expect(s.maxF).toBeGreaterThan(s.F0);
    // the decrease is many orders of magnitude above the rounding bound
    const n = s.firstDecrease;
    expect(-t.dF[n]).toBeGreaterThan(1e6 * t.dFerr[n]);
  });
});

describe("summary wording inputs", () => {
  it("reports no decrease for an ideal purification run", () => {
    const s = summarize(run(0.3, 0, 30));
    expect(s.firstDecrease).toBe(-1);
  });
  it("bellPopulations matches the trajectory buffers", () => {
    const t = run(0.37, 0.021, 10);
    for (let n = 0; n <= 10; n++) {
      expect(bellPopulations(t.u[n], t.v[n])).toEqual([t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]]);
    }
  });
});
