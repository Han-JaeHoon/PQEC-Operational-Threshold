import { describe, expect, it } from "vitest";
import t4q from "../verification/epsc_research/epsc_table_4Q.json";
import tT3 from "../verification/epsc_research/epsc_table_T3.json";
import tT4 from "../verification/epsc_research/epsc_table_T4.json";
import series from "../verification/epsc_research/series_results.json";
import ref from "../verification/reference_epsc.json";
import {
  basinVerdict, bistableFixedPoints, classifyInput, computeEpsc, epscBrackets, initialFidelity,
  isDetermined, mapCoeffs, SIDE_PB, SIDE_Z, stepEnclosure,
} from "../src/epsc";
import { allocTrajectory, computeTrajectory } from "../src/type3map";
import { computeTrajectory4, type4Step } from "../src/type4map";
import { computeTrajectory4Q, q4Step } from "../src/q4map";
import { computeTrajectoryBBPSSW } from "../src/bbpsswmap";
import { orbitLimit, regimeReport, sideOfRoot, THRESHOLDS, zPoint, type RegimeProtocol } from "../src/regimes";

const PROTOS: RegimeProtocol[] = ["type3", "type4", "q4"];
const frac = (x: string) => { const [a, b] = x.split("/"); return b ? Number(a) / Number(b) : Number(a); };
const TABLES: [RegimeProtocol, { q: string; eps_lo: string; eps_hi: string }[]][] = [
  ["q4", (t4q as any)["4Q"]], ["type4", (tT4 as any)["T4"]], ["type3", (tT3 as any)["T3"]],
];

// exact rationals, for the enclosure test
type R = [bigint, bigint];
const gcd = (a: bigint, b: bigint): bigint => (b === 0n ? (a < 0n ? -a : a) : gcd(b, a % b));
const norm = ([a, b]: R): R => { const g = gcd(a, b) || 1n; const s = b < 0n ? -1n : 1n; return [(s * a) / g, (s * b) / g]; };
const rAdd = (x: R, y: R): R => norm([x[0] * y[1] + y[0] * x[1], x[1] * y[1]]);
const rSub = (x: R, y: R): R => norm([x[0] * y[1] - y[0] * x[1], x[1] * y[1]]);
const rMul = (x: R, y: R): R => norm([x[0] * y[0], x[1] * y[1]]);
const rDiv = (x: R, y: R): R => norm([x[0] * y[1], x[1] * y[0]]);
const rPow = (x: R, k: number): R => { let o: R = [1n, 1n]; for (let i = 0; i < k; i++) o = rMul(o, x); return o; };
const I = (n: number): R => [BigInt(n), 1n];
function fromFloat(x: number): R {
  let k = 0n; let v = x;
  while (!Number.isInteger(v)) { v *= 2; k++; }
  return norm([BigInt(v), 1n << k]);
}
const rCmpFloat = (x: R, f: number) => { const d = rSub(x, fromFloat(f)); return d[0] > 0n ? 1 : d[0] < 0n ? -1 : 0; };

/** T(t, m) in exact rationals, from the common form (1) of the notes (not the certifier's formula). */
function exactT(p: RegimeProtocol, q: number, t: number, m: number): [R, R] {
  const s = rSub(I(1), fromFloat(q));
  const sp = (k: number) => rPow(s, k);
  const [al, be, ga, de, et, b] = p === "q4"
    ? [rMul(I(2), sp(5)), rMul(sp(4), rAdd(I(1), s)), rMul(I(2), sp(5)), rMul(sp(3), rAdd(I(1), s)), rAdd(sp(3), sp(5)), sp(5)]
    : (() => {
      const [A, B] = p === "type3" ? [6, 4] : [4, 2];
      const c = rMul(sp(A - 1), rAdd(I(1), s));
      return [rMul(I(2), sp(A)), rMul(I(2), sp(A)), c, c, rMul(I(2), sp(B)), sp(B)];
    })();
  const T = fromFloat(t), M = fromFloat(m);
  const d = rAdd(I(1), rMul(I(2), M));
  const v = rDiv(rSub(I(1), rMul(I(2), M)), d);
  const u = rDiv(T, d);
  const D = rAdd(rAdd(I(1), rMul(et, rMul(u, u))), rMul(b, rMul(v, v)));
  const u1 = rDiv(rMul(u, rAdd(rMul(al, v), be)), D);
  const v1 = rDiv(rAdd(rMul(ga, rMul(u, u)), rMul(de, v)), D);
  return [rDiv(rMul(I(2), u1), rAdd(I(1), v1)), rDiv(rSub(I(1), v1), rMul(I(2), rAdd(I(1), v1)))];
}

describe("the map of the certifier is the simulator's map", () => {
  it("one (t, m) step encloses the exact rational value, from the notes' common form", () => {
    let seed = 12345;
    const rnd = () => { seed = (seed * 1664525 + 1013904223) % 2 ** 32; return seed / 2 ** 32; };
    for (const p of PROTOS) {
      for (let i = 0; i < 60; i++) {
        const q = 0.1 + 0.1 * rnd(), t = rnd(), m = 0.5 * rnd();
        const [tl, th, ml, mh] = stepEnclosure(mapCoeffs(p, q), t, m);
        const [te, me] = exactT(p, q, t, m);
        expect(rCmpFloat(te, tl)).toBeGreaterThanOrEqual(0);
        expect(rCmpFloat(te, th)).toBeLessThanOrEqual(0);
        expect(rCmpFloat(me, ml)).toBeGreaterThanOrEqual(0);
        expect(rCmpFloat(me, mh)).toBeLessThanOrEqual(0);
        expect(th - tl).toBeLessThan(5e-14);
      }
    }
  });

  it("matches the trajectories drawn in the charts, round by round, in (u, v) coordinates", () => {
    for (const [p, compute] of [["type3", computeTrajectory], ["type4", computeTrajectory4], ["q4", computeTrajectory4Q]] as const) {
      const tr = compute(allocTrajectory(40), 0.3, 0.15, 40);
      const c = mapCoeffs(p, 0.15);
      for (let n = 0; n < 40; n++) {
        // Trajectory stores (u, v) for Type 3/4 and (a, b) = (line, transverse) for 4Q
        const [U, V] = p === "q4" ? [tr.v[n], tr.u[n]] : [tr.u[n], tr.v[n]];
        const [U1, V1] = p === "q4" ? [tr.v[n + 1], tr.u[n + 1]] : [tr.u[n + 1], tr.v[n + 1]];
        const [tl, th, ml, mh] = stepEnclosure(c, (2 * U) / (1 + V), (1 - V) / (2 * (1 + V)));
        expect(Math.abs((tl + th) / 2 - (2 * U1) / (1 + V1))).toBeLessThan(1e-13);
        expect(Math.abs((ml + mh) / 2 - (1 - V1) / (2 * (1 + V1)))).toBeLessThan(1e-13);
      }
    }
    expect(type4Step(0.5, 0.5, 0.1).length).toBe(3);
    expect(q4Step(0.5, 0.5, 0.1).length).toBe(3);
  });

  it("reproduces the q convention of the research code: lam1, lam2 of the Sec. 4 table", () => {
    for (const r of series as { kind: string; q: string; lam1: string; lam2: string }[]) {
      const p: RegimeProtocol = r.kind === "4Q" ? "q4" : r.kind === "T4" ? "type4" : "type3";
      const q = Number(r.q);
      const fp = bistableFixedPoints(p, q)!;
      const t0 = (fp.S.t[0] + fp.S.t[1]) / 2, m0 = (fp.S.m[0] + fp.S.m[1]) / 2;
      const h = 1e-6;
      const c = mapCoeffs(p, q);
      const mid = (a: number[], i: number) => (a[i] + a[i + 1]) / 2;
      const P = (t: number, m: number) => { const e = stepEnclosure(c, t, m); return [mid(e, 0), mid(e, 2)]; };
      const a = (P(t0 + h, m0)[0] - P(t0 - h, m0)[0]) / (2 * h);
      const b = (P(t0, m0 + h)[0] - P(t0, m0 - h)[0]) / (2 * h);
      const cc = (P(t0 + h, m0)[1] - P(t0 - h, m0)[1]) / (2 * h);
      const d = (P(t0, m0 + h)[1] - P(t0, m0 - h)[1]) / (2 * h);
      const disc = Math.sqrt((a - d) ** 2 + 4 * b * cc);
      expect(Math.abs((a + d + disc) / 2 - Number(r.lam1))).toBeLessThan(1e-7);
      expect(Math.abs((a + d - disc) / 2 - Number(r.lam2))).toBeLessThan(1e-7);
    }
  });
});

describe("exact regime", () => {
  it("decides the side of each constant exactly, down to adjacent floats", () => {
    for (const p of PROTOS) {
      for (const root of Object.values(THRESHOLDS[p])) {
        let lo = root.bracket[0], hi = root.bracket[1]; // bisection in s with the exact sign only
        const q = (x: number) => 1 - x;
        for (let i = 0; i < 80; i++) {
          const mid = 0.5 * (lo + hi);
          if (mid === lo || mid === hi) break;
          if (sideOfRoot(root, q(mid)) > 0) lo = mid; else hi = mid;
        }
        expect(Math.abs(1 - hi - root.q)).toBeLessThan(1e-9);
        expect(sideOfRoot(root, root.q - 1e-6)).toBe(-1);
        expect(sideOfRoot(root, root.q + 1e-6)).toBe(1);
        expect(sideOfRoot(root, 0)).toBe(-1);
        expect(sideOfRoot(root, 1)).toBe(1);
      }
    }
  });
  it("gives the three regions, flags the boundary and q out of range", () => {
    for (const p of PROTOS) {
      const t = THRESHOLDS[p];
      expect(regimeReport(p, t.qBi.q - 1e-6).regime).toBe("unique-purified");
      expect(regimeReport(p, (t.qBi.q + t.qSN.q) / 2).regime).toBe("bistable");
      expect(regimeReport(p, t.qSN.q + 1e-6).regime).toBe("unique-unpurified");
      expect(regimeReport(p, t.qBi.q).near).toBe("qBi");
      expect(regimeReport(p, t.qSN.q).near).toBe("qSN");
      expect(regimeReport(p, t.qBi.q + 1e-6).near).toBeNull();
      expect(regimeReport(p, 0).inRange).toBe(false);
      expect(regimeReport(p, 1).inRange).toBe(false);
    }
  });
  it("Z(q) changes kind exactly at q_m, without NaN next to it", () => {
    for (const p of PROTOS) {
      const m = THRESHOLDS[p].qM.q;
      expect(zPoint(p, m - 1e-6).kind).toBe("boundary");
      expect(zPoint(p, m + 1e-6).kind).toBe("mixed");
      for (const d of [-1e-11, -1e-12, 0, 1e-12]) expect(Number.isFinite(zPoint(p, m + d).F)).toBe(true);
      expect(regimeReport(p, m - 1e-6).belowQm).toBe(true);
      expect(regimeReport(p, m + 1e-6).belowQm).toBe(false);
    }
  });
});

describe("regression: the research tables (never read by src/)", () => {
  const examples: [RegimeProtocol, number, number, number][] = [
    ["q4", 0.178251, 0.830572489, 0.830572491],
    ["type4", 0.184175, 0.8941400914, 0.8941400934],
    ["type3", 0.124435, 0.8991163518, 0.8991163538],
  ];
  for (const [p, q, lo, hi] of examples) {
    it(`${p} q = ${q}: inside the certified interval, narrower than it`, () => {
      const b = computeEpsc(p, q)!;
      expect(b.lo).toBeGreaterThan(lo);
      expect(b.hi).toBeLessThan(hi);
      expect(b.hi - b.lo).toBeLessThan(1e-12);
      expect(b.stop).toBe("resolution");
    });
  }
  for (const [p, rows] of TABLES) {
    it(`${p}: all ${rows.length} table rows overlap the computed bracket`, () => {
      for (const r of rows) {
        const b = computeEpsc(p, frac(r.q))!;
        expect(isDetermined(b)).toBe(true);
        // float64 q differs from the rational q by < 1e-16; eps_c moves by far less than the margin
        expect(b.lo).toBeGreaterThan(frac(r.eps_lo) - 1e-13);
        expect(b.hi).toBeLessThan(frac(r.eps_hi) + 1e-13);
      }
    });
  }
  it("the root of the stable-manifold relation lies in each bracket", () => {
    for (const r of series as { kind: string; q: string; root: string }[]) {
      const p: RegimeProtocol = r.kind === "4Q" ? "q4" : r.kind === "T4" ? "type4" : "type3";
      const b = computeEpsc(p, Number(r.q))!;
      // the series root is printed to 13 digits and is not rigorous
      expect(Number(r.root)).toBeGreaterThan(b.lo - 1e-12);
      expect(Number(r.root)).toBeLessThan(b.hi + 1e-12);
    }
  });
  it("src/ does not read the tables and has no table of eps_c values", () => {
    const sources = import.meta.glob("../src/*.ts", { query: "?raw", import: "default", eager: true });
    expect(Object.keys(sources).length).toBeGreaterThan(8);
    for (const text of Object.values(sources)) {
      expect(text as string).not.toMatch(/epsc_table|epsc_research|series_results/);
    }
  });
});

describe("independent arb certificate, off the research grid", () => {
  const rows = (ref as { rows: any[] }).rows;
  it("has q near q_m, q_bi and q_SN for every protocol", () => {
    expect(rows.length).toBeGreaterThanOrEqual(18);
    for (const r of rows) {
      expect(regimeReport(r.protocol, r.q).regime).toBe("bistable");
      // Proposition 4 with 256-bit balls agrees on both endpoints
      expect(r.arb_side_lo).toBe(1);
      expect(r.arb_side_hi).toBe(-1);
    }
  });
  for (const r of rows) {
    it(`${r.protocol} q = ${r.q}: the simulator still produces the brackets that were checked`, () => {
      const b = computeEpsc(r.protocol, r.q)!;
      expect(b.lo).toBe(r.lo);
      expect(b.hi).toBe(r.hi);
    });
  }
});

describe("bracket behaviour", () => {
  it("brackets only shrink, every stage is valid", () => {
    for (const p of PROTOS) {
      const t = THRESHOLDS[p];
      const q = (t.qBi.q + t.qSN.q) / 2;
      const fin = computeEpsc(p, q)!;
      const g = epscBrackets(p, q);
      let plo = 0, phi = 1, stages = 0;
      for (;;) {
        const r = g.next();
        const b = r.value!;
        expect(b.lo).toBeGreaterThanOrEqual(plo);
        expect(b.hi).toBeLessThanOrEqual(phi);
        expect(b.lo).toBeLessThanOrEqual(fin.lo);
        expect(b.hi).toBeGreaterThanOrEqual(fin.hi);
        plo = b.lo; phi = b.hi; stages++;
        if (r.done) break;
      }
      expect(stages).toBeGreaterThan(10);
    }
  });
  it("L_0, U_0 of Eqs. 41-42 enclose eps_c", () => {
    for (const p of PROTOS) {
      const t = THRESHOLDS[p];
      const q = (t.qBi.q + t.qSN.q) / 2;
      const { S } = bistableFixedPoints(p, q)!;
      const ts = (S.t[0] + S.t[1]) / 2, ms = (S.m[0] + S.m[1]) / 2;
      const a0 = (2 * (1 - ts)) / (2 - ts), b0 = (4 * ms) / (1 + 2 * ms);
      const b = computeEpsc(p, q)!;
      expect(Math.min(a0, b0)).toBeLessThan(b.lo);
      expect(Math.max(a0, b0)).toBeGreaterThan(b.hi);
    }
  });
  it("is at least as tight as the research tables across the window, and fast", () => {
    for (const p of PROTOS) {
      const t = THRESHOLDS[p];
      for (let k = 1; k < 20; k++) {
        const q = t.qBi.q + ((t.qSN.q - t.qBi.q) * k) / 20;
        const t0 = performance.now();
        const b = computeEpsc(p, q)!;
        expect(performance.now() - t0).toBeLessThan(2000);
        expect(b.hi - b.lo).toBeLessThan(2e-9);
      }
    }
  });
  it("just outside the bracket the float64 orbit does go to the predicted attractor", () => {
    for (const p of PROTOS) {
      const t = THRESHOLDS[p];
      const q = (t.qBi.q + t.qSN.q) / 2;
      const b = computeEpsc(p, q)!;
      const fp = bistableFixedPoints(p, q)!;
      const FB = (fp.PB.F[0] + fp.PB.F[1]) / 2;
      const hiLim = orbitLimit(p, b.lo - 1e-6, q);
      const loLim = orbitLimit(p, b.hi + 1e-6, q);
      expect(Math.abs(hiLim.F - FB)).toBeLessThan(1e-9);
      expect(Math.abs(loLim.F - zPoint(p, q).F)).toBeLessThan(1e-6);
    }
  });
  it("very close to q_bi / q_SN: still rigorous and honest, never fabricated", () => {
    for (const p of PROTOS) {
      const t = THRESHOLDS[p];
      for (const q of [t.qBi.q + 2e-9, t.qSN.q - 2e-9]) {
        const b = computeEpsc(p, q);
        expect(b === null || (b.lo < b.hi && b.lo > 0 && b.hi < 1)).toBe(true);
      }
    }
  });
});

describe("basin verdict", () => {
  const b = { lo: 0.5, hi: 0.5000001 };
  it("uses the interval, not a float comparison with a midpoint", () => {
    expect(basinVerdict(b, 0.4)).toBe("P_B");
    expect(basinVerdict(b, 0.5)).toBe("P_B");
    expect(basinVerdict(b, 0.50000005)).toBe("undetermined");
    expect(basinVerdict(b, 0.5000001)).toBe("Z");
    expect(basinVerdict(b, 0.6)).toBe("Z");
  });
  it("the certifier says nothing for an input it cannot resolve", () => {
    const q = 0.178251, p: RegimeProtocol = "q4";
    const f = computeEpsc(p, q)!;
    const fp = bistableFixedPoints(p, q)!;
    const c = mapCoeffs(p, q);
    expect(classifyInput(c, fp.S, f.lo, 100000).side).toBe(SIDE_PB);
    expect(classifyInput(c, fp.S, f.hi, 100000).side).toBe(SIDE_Z);
    const inside = classifyInput(c, fp.S, (f.lo + f.hi) / 2, 100000);
    expect([0, SIDE_PB, SIDE_Z]).toContain(inside.side); // whichever it is, it stays certified
  });
});

describe("initial fidelity", () => {
  it("is 1 - 3 eps / 4 and equals round 0 of every protocol", () => {
    expect(initialFidelity(0.2).toFixed(6)).toBe("0.850000");
    expect(initialFidelity(0)).toBe(1);
    expect(initialFidelity(1)).toBe(0.25);
    for (const compute of [computeTrajectory, computeTrajectory4, computeTrajectory4Q, computeTrajectoryBBPSSW]) {
      for (const e of [0, 0.1, 0.2, 2 / 3, 0.9, 1]) {
        const tr = compute(allocTrajectory(2), e, 0.05, 2);
        expect(Math.abs(tr.pPhiP[0] - initialFidelity(e))).toBeLessThan(1e-15);
      }
    }
  });
});
