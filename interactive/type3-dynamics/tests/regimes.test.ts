import { describe, expect, it } from "vitest";
import { allocTrajectory, computeTrajectory, type Trajectory } from "../src/type3map";
import { computeTrajectory4 } from "../src/type4map";
import { computeTrajectory4Q } from "../src/q4map";
import {
  evalPoly, iterateLogOrbit, orbitLimit, regime, solveRoot, THRESHOLDS, zPoint, type RegimeProtocol,
} from "../src/regimes";

const PROTOS: RegimeProtocol[] = ["type3", "type4", "q4"];
const COMPUTE = { type3: computeTrajectory, type4: computeTrajectory4, q4: computeTrajectory4Q };

function log(label: string, x: unknown) {
  console.log(`${label}: ${JSON.stringify(x)}`);
}

describe("Test 1: constants against bisection of their polynomials", () => {
  for (const p of PROTOS) {
    for (const key of ["qBi", "qSN", "qM"] as const) {
      it(`${p} ${key}`, () => {
        const root = THRESHOLDS[p][key];
        const q = solveRoot(root);
        log(`${p} ${key} bisected - constant`, q - root.q);
        expect(Math.abs(q - root.q)).toBeLessThan(1e-9);
      });
    }
  }

  it("the q_SN bracket holds the largest root of its polynomial in (0,1)", () => {
    for (const p of PROTOS) {
      const r = THRESHOLDS[p].qSN;
      // no sign change between the bracket and s = 1
      let sign = Math.sign(evalPoly(r.poly, r.bracket[1]));
      for (let s = r.bracket[1]; s <= 1; s += 1e-4) {
        expect(Math.sign(evalPoly(r.poly, s))).toBe(sign);
        sign = Math.sign(evalPoly(r.poly, s));
      }
    }
  });
});

describe("Test 2: q = q_bi - 0.002, limit independent of eps", () => {
  const EXPECT = { type3: 0.6267, type4: 0.5546, q4: 0.5829 };
  for (const p of PROTOS) {
    it(p, () => {
      const q = THRESHOLDS[p].qBi.q - 0.002;
      expect(regime(p, q)).toBe("unique-purified");
      const res = [0.5, 0.9, 0.99].map((e) => orbitLimit(p, e, q));
      log(`${p} F, rounds`, res.map((r) => [r.F, r.rounds]));
      for (const r of res) {
        expect(r.converged).toBe(true);
        expect(r.F).toBeCloseTo(EXPECT[p], 3);
        expect(Math.abs(r.F - res[0].F)).toBeLessThan(1e-10);
      }
    });
  }
});

describe("Test 3: q = q_SN + 0.002, limit is Z(q)", () => {
  for (const p of PROTOS) {
    it(p, () => {
      const q = THRESHOLDS[p].qSN.q + 0.002;
      expect(regime(p, q)).toBe("unique-unpurified");
      const Z = zPoint(p, q);
      const res = [0.1, 0.5, 0.9].map((e) => orbitLimit(p, e, q));
      log(`${p} Z.F, F, rounds`, [Z.F, res.map((r) => [r.F, r.rounds])]);
      expect(Z).toEqual({ kind: "mixed", line: 0, F: 0.25 });
      for (const r of res) {
        expect(r.converged).toBe(true);
        expect(Math.abs(r.F - Z.F)).toBeLessThan(1e-10);
      }
    });
  }

  it("4Q on q_SN < q < q_m: Z is the boundary point B*, with z = 0", () => {
    const t = THRESHOLDS.q4;
    const q = 0.5 * (t.qSN.q + t.qM.q);
    const Z = zPoint("q4", q);
    expect(Z.kind).toBe("boundary");
    expect(Z.F).toBeGreaterThan(0.25);
    for (const e of [0.1, 0.5, 0.9]) {
      const r = orbitLimit("q4", e, q);
      log(`4Q window eps=${e} F, Z.F, rounds`, [r.F, Z.F, r.rounds]);
      expect(r.converged).toBe(true);
      expect(Math.abs(r.F - Z.F)).toBeLessThan(1e-9);
    }
  });

  it("Type 3/4 have q_SN > q_m, so Z is the mixed point on all of q > q_SN", () => {
    for (const p of ["type3", "type4"] as const) {
      expect(THRESHOLDS[p].qSN.q).toBeGreaterThan(THRESHOLDS[p].qM.q);
      expect(zPoint(p, THRESHOLDS[p].qSN.q + 1e-9).kind).toBe("mixed");
    }
  });
});

describe("Test 4: bistability between q_bi and q_SN", () => {
  for (const p of PROTOS) {
    it(`${p}, midpoint`, () => {
      const t = THRESHOLDS[p];
      const q = 0.5 * (t.qBi.q + t.qSN.q);
      expect(regime(p, q)).toBe("bistable");
      const hi = orbitLimit(p, 0.1, q);
      const lo = orbitLimit(p, 0.999, q);
      log(`${p} q, F(eps=.1), F(eps=.999), Z.F`, [q, hi.F, lo.F, zPoint(p, q).F]);
      expect(hi.converged && lo.converged).toBe(true);
      expect(hi.F).toBeGreaterThan(0.5);
      expect(lo.F).toBeLessThan(0.5);
      expect(Math.abs(lo.F - zPoint(p, q).F)).toBeLessThan(1e-10);
    });
  }

  it("just below q_SN the upper branch is itself below 1/2, still above the low attractor", () => {
    for (const p of PROTOS) {
      const q = THRESHOLDS[p].qSN.q - 1e-4;
      const hi = orbitLimit(p, 0.1, q);
      log(`${p} q_SN - 1e-4: F(eps=.1), Z.F`, [hi.F, zPoint(p, q).F]);
      expect(hi.converged).toBe(true);
      expect(hi.F).toBeLessThan(0.5);
      expect(hi.F).toBeGreaterThan(zPoint(p, q).F + 0.1);
    }
  });
});

describe("Test 5: the constants do not depend on eps", () => {
  it("regime() takes no eps, and the limit stays unique across eps at fixed q", () => {
    expect(regime.length).toBe(2);
    for (const p of PROTOS) {
      const before = JSON.stringify(THRESHOLDS[p]);
      for (const e of [0.01, 0.5, 0.99]) orbitLimit(p, e, 0.05);
      expect(JSON.stringify(THRESHOLDS[p])).toBe(before);
    }
  });
});

describe("regime boundaries", () => {
  it("q_bi and q_SN themselves are bistable", () => {
    for (const p of PROTOS) {
      const t = THRESHOLDS[p];
      expect(regime(p, t.qBi.q)).toBe("bistable");
      expect(regime(p, t.qSN.q)).toBe("bistable");
      expect(regime(p, t.qBi.q - 1e-12)).toBe("unique-purified");
      expect(regime(p, t.qSN.q + 1e-12)).toBe("unique-unpurified");
    }
  });
});

describe("the limit orbit is the same map as the displayed one", () => {
  const trans = (t: Trajectory, p: RegimeProtocol, n: number) => (p === "q4" ? t.v[n] : t.u[n]);
  const line = (t: Trajectory, p: RegimeProtocol, n: number) => (p === "q4" ? t.u[n] : t.v[n]);

  it("agrees with computeTrajectory* round by round", () => {
    for (const p of PROTOS) {
      for (const [e, q] of [[0.1, 0.03], [0.5, 0.1], [0.9, 0.15], [0.99, 0.2], [0.3, 0.5]]) {
        const t = COMPUTE[p](allocTrajectory(60), e, q, 60);
        for (const n of [1, 2, 5, 17, 60]) {
          const r = orbitLimit(p, e, q, n);
          if (r.rounds !== n) continue; // converged earlier
          expect(Math.abs(r.line - line(t, p, n))).toBeLessThan(1e-12);
          expect(Math.abs(r.transverse - trans(t, p, n))).toBeLessThan(1e-12);
          expect(Math.abs(r.F - t.pPhiP[n])).toBeLessThan(1e-12);
        }
      }
    }
  });
});

describe("numerics: underflow and non-convergence", () => {
  const nearestOne = 1 - 2 ** -53;

  it("eps next to 1 below q_bi still reaches P_B", () => {
    for (const p of PROTOS) {
      for (const gap of [0.002, 2e-4]) {
        const q = THRESHOLDS[p].qBi.q - gap;
        const r = orbitLimit(p, nearestOne, q);
        log(`${p} q_bi-${gap} eps=1-2^-53 F, rounds, min log10 transverse`,
          [r.F, r.rounds, r.minLog10Transverse]);
        expect(r.converged).toBe(true);
        expect(Math.abs(r.F - orbitLimit(p, 0.5, q).F)).toBeLessThan(1e-8);
      }
    }
  });

  it("a transverse part below float64 range regrows from its logarithm (plain float64: stuck at 0)", () => {
    for (const p of PROTOS) {
      const q = THRESHOLDS[p].qBi.q - 0.01;
      const start = 0.9;
      const tau0 = -2000; // 1e-869: exactly 0 as a plain float64
      expect(Math.exp(tau0)).toBe(0);
      const r = iterateLogOrbit(p, start, tau0, q);
      log(`${p} tau0=-2000 F, rounds`, [r.F, r.rounds]);
      expect(r.belowFloat64).toBe(true);
      expect(r.converged).toBe(true);
      expect(Math.abs(r.F - orbitLimit(p, 0.5, q).F)).toBeLessThan(1e-9);
    }
  });

  it("a plateau is not mistaken for convergence: too few rounds report not converged", () => {
    for (const p of PROTOS) {
      const q = THRESHOLDS[p].qBi.q - 0.002;
      const r = orbitLimit(p, 0.99, q, 100);
      expect(r.converged).toBe(false);
      expect(r.rounds).toBe(100);
    }
  });

  it("reports not converged at the round limit next to q_bi, rather than a wrong plateau", () => {
    const q = THRESHOLDS.type3.qBi.q - 1e-9;
    const r = orbitLimit("type3", 0.9999999, q);
    log("type3 q_bi-1e-9 converged, rounds, F", [r.converged, r.rounds, r.F]);
    expect(r.converged).toBe(false);
    expect(r.rounds).toBe(200_000);
  });

  it("handles the excluded ends without NaN", () => {
    for (const p of PROTOS) {
      for (const [e, q] of [[1, 0.05], [0, 0.05], [0.5, 1], [0.5, 0], [1, 1]]) {
        const r = orbitLimit(p, e, q);
        expect(Number.isFinite(r.F)).toBe(true);
        expect(r.converged).toBe(true);
      }
      expect(orbitLimit(p, 1, 0.05).F).toBe(0.25);
      expect(orbitLimit(p, 0.5, 1).F).toBe(0.25);
    }
  });

  it("the worst-case cost stays within a frame budget", () => {
    const t0 = performance.now();
    orbitLimit("type3", 0.9999999, THRESHOLDS.type3.qBi.q - 1e-9);
    const ms = performance.now() - t0;
    log("full-length limit orbit, ms", ms);
    expect(ms).toBeLessThan(100);
  });
});
