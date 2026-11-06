import { describe, expect, it } from "vitest";
import { bellFixedPoint, computeEpsc } from "../src/epsc";
import { gainReport } from "../src/gainPanel";
import {
  epsCoverage, f0Enclosure, fbEnclosure, gainVerdict, Q_ENT, qGainBracket, sideOfQEnt,
} from "../src/qgain";
import { type3Step } from "../src/type3map";
import { type4Step } from "../src/type4map";
import { q4Step } from "../src/q4map";
import { orbitLimit, regimeReport, THRESHOLDS, type RegimeProtocol } from "../src/regimes";

const PROTOS: RegimeProtocol[] = ["type3", "type4", "q4"];
const mid = (x: readonly [number, number]) => (x[0] + x[1]) / 2;

// exact rationals for the resultant test
type R = [bigint, bigint];
const gcd = (a: bigint, b: bigint): bigint => (b === 0n ? (a < 0n ? -a : a) : gcd(b, a % b));
const norm = ([a, b]: R): R => { const g = gcd(a, b) || 1n; const s = b < 0n ? -1n : 1n; return [(s * a) / g, (s * b) / g]; };
const add = (x: R, y: R): R => norm([x[0] * y[1] + y[0] * x[1], x[1] * y[1]]);
const sub = (x: R, y: R): R => norm([x[0] * y[1] - y[0] * x[1], x[1] * y[1]]);
const mul = (x: R, y: R): R => norm([x[0] * y[0], x[1] * y[1]]);
const pw = (x: R, k: number): R => { let o: R = [1n, 1n]; for (let i = 0; i < k; i++) o = mul(o, x); return o; };
const I = (n: number): R => [BigInt(n), 1n];
function fromFloat(x: number): R { let k = 0n; let v = x; while (!Number.isInteger(v)) { v *= 2; k++; } return norm([BigInt(v), 1n << k]); }

/** Coefficient table of Section 1, typed in from the note. */
function coeffs(p: RegimeProtocol, s: R) {
  const sp = (k: number) => pw(s, k);
  const ops = add(I(1), s);
  if (p === "q4") return { a: mul(I(2), sp(5)), be: mul(sp(4), ops), g: mul(I(2), sp(5)), d: mul(sp(3), ops), et: add(sp(3), sp(5)), b: sp(5) };
  const [A, B] = p === "type3" ? [6, 4] : [4, 2];
  const c = mul(sp(A - 1), ops);
  return { a: mul(I(2), sp(A)), be: mul(I(2), sp(A)), g: c, d: c, et: mul(I(2), sp(B)), b: sp(B) };
}

/** Sign of R_P(s, eps), Eq. (18), in exact rational arithmetic. */
function resultantSign(p: RegimeProtocol, s: R, eps: R): number {
  const { a, be, g, d, et, b } = coeffs(p, s);
  const th = sub(I(1), eps);
  const h2 = add(mul(g, b), mul(et, a));
  const h1 = sub(mul(et, sub(be, d)), mul(g, a));
  const h0 = mul(g, sub(I(1), be));
  const e2 = sub(g, mul(I(4), a));
  const e1 = sub(mul(mul(I(-6), g), th), mul(I(4), sub(be, d)));
  const e0 = mul(mul(I(9), g), mul(th, th));
  const t1 = sub(mul(h2, e0), mul(h0, e2));
  const r = sub(mul(t1, t1), mul(sub(mul(h2, e1), mul(h1, e2)), sub(mul(h1, e0), mul(h0, e1))));
  return r[0] > 0n ? 1 : r[0] < 0n ? -1 : 0;
}

const STEP = { type3: type3Step, type4: type4Step, q4: q4Step } as const;
const EPS_GRID = Array.from({ length: 133 }, (_, i) => 0.005 * (i + 1)); // 0.005 ... 0.665

describe("the note's fixed point is the simulator's fixed point", () => {
  it("P_B from Eq. (12) is fixed by the simulator's own one-round map, for q < q_bi and in the window", () => {
    for (const p of PROTOS) {
      const t = THRESHOLDS[p];
      for (const q of [0.001, 0.02, 0.05, t.qBi.q * 0.9, (t.qBi.q + t.qSN.q) / 2, t.qSN.q - 1e-4]) {
        const pb = bellFixedPoint(p, q)!;
        const [u, v] = [mid(pb.u), mid(pb.v)];
        // the simulator stores 4Q as (a, b) = (line x, transverse z) and Type 3/4 as (u, v)
        const [x1, y1] = p === "q4" ? STEP.q4(v, u, q) : STEP[p](u, v, q);
        const [u1, v1] = p === "q4" ? [y1, x1] : [x1, y1];
        // float64 evaluation of one round at the enclosure midpoint
        expect(Math.abs(u1 - u)).toBeLessThan(1e-12);
        expect(Math.abs(v1 - v)).toBeLessThan(1e-12);
        expect(pb.F[1] - pb.F[0]).toBeLessThan(1e-12); // enclosure width; largest next to q_SN
      }
    }
  });
  it("is the Bell-connected branch: the float64 orbit from a good input ends there for q < q_bi", () => {
    for (const p of PROTOS) {
      const q = THRESHOLDS[p].qBi.q * 0.8;
      expect(Math.abs(orbitLimit(p, 0.1, q).F - mid(fbEnclosure(p, q)!))).toBeLessThan(1e-9);
    }
  });
  it("F_B(0) = 1 (Lemma 7.1)", () => {
    for (const p of PROTOS) expect(Math.abs(mid(fbEnclosure(p, 0)!) - 1)).toBeLessThan(1e-14);
  });
});

describe("q_ent", () => {
  it("F_B(q_ent) = 1/2 and q_bi < q_ent < q_SN (Eqs. 25, 58)", () => {
    for (const p of PROTOS) {
      const qe = Q_ENT[p].q;
      const F = fbEnclosure(p, qe)!;
      expect(Math.abs(mid(F) - 0.5)).toBeLessThan(1e-12);
      expect(qe).toBeGreaterThan(THRESHOLDS[p].qBi.q);
      expect(qe).toBeLessThan(THRESHOLDS[p].qSN.q);
      expect(sideOfQEnt(p, qe - 1e-9)).toBe(-1);
      expect(sideOfQEnt(p, qe + 1e-9)).toBe(1);
    }
  });
  it("agrees with the earlier values in RESEARCH_SUMMARY.md", () => {
    expect(Math.abs(Q_ENT.type3.q - 0.129137470776)).toBeLessThan(1e-11);
    expect(Math.abs(Q_ENT.type4.q - 0.184735805365)).toBeLessThan(1e-11);
  });
});

describe("q_gain(eps) = q_F(eps)", () => {
  it("residual: F_B(q_gain) = F_0 to the enclosure width, and the end points are certified", () => {
    for (const p of PROTOS) {
      for (const eps of EPS_GRID) {
        const b = qGainBracket(p, eps)!;
        const F0 = f0Enclosure(eps);
        expect(fbEnclosure(p, b.lo)![0]).toBeGreaterThan(F0[1]);
        expect(fbEnclosure(p, b.hi)![1]).toBeLessThan(F0[0]);
        expect(Math.abs(mid(fbEnclosure(p, (b.lo + b.hi) / 2)!) - (1 - 0.75 * eps))).toBeLessThan(1e-12);
        expect(b.hi - b.lo).toBeLessThan(5e-11); // measured max 1.1e-11 (4Q)
        expect(b.hi).toBeLessThan(Q_ENT[p].q + 1e-12);
        expect(b.lo).toBeGreaterThan(0);
      }
    }
  });
  it("the bracket contains a root of the exact resultant R_P(s, eps), Eq. (18)", () => {
    for (const p of PROTOS) {
      for (const eps of [0.01, 0.1, 0.25, 0.4, 0.55, 0.66]) {
        const b = qGainBracket(p, eps)!;
        const e = fromFloat(eps);
        const a = resultantSign(p, sub(I(1), fromFloat(b.lo)), e);
        const c = resultantSign(p, sub(I(1), fromFloat(b.hi)), e);
        expect(a * c).toBe(-1);
      }
    }
  });
  it("is strictly increasing on 0 < eps < 2/3 (Eq. 83): brackets are disjoint and ordered", () => {
    for (const p of PROTOS) {
      let prev = -Infinity;
      for (const eps of EPS_GRID) {
        const b = qGainBracket(p, eps)!;
        expect(b.lo).toBeGreaterThan(prev);
        prev = b.hi;
      }
    }
  });
  it("limits: q_gain -> 0 as eps -> 0+, q_gain -> q_ent as eps -> 2/3- (Eq. 84)", () => {
    for (const p of PROTOS) {
      expect(qGainBracket(p, 1e-12)!.hi).toBeLessThan(1e-11);
      const b = qGainBracket(p, 0.6666666666666666)!; // the largest float below 2/3
      expect(Math.abs((b.lo + b.hi) / 2 - Q_ENT[p].q)).toBeLessThan(1e-12);
    }
  });
  it("verdict: gain below, no gain above, uncertain only inside the bracket", () => {
    for (const p of PROTOS) {
      for (const eps of [0.05, 0.3, 0.6]) {
        const b = qGainBracket(p, eps)!;
        expect(gainVerdict(b, b.lo * 0.99)).toBe("gain");
        expect(gainVerdict(b, b.lo)).toBe("gain");
        expect(gainVerdict(b, b.hi)).toBe("no-gain");
        expect(gainVerdict(b, (b.hi + Q_ENT[p].q) / 2)).toBe("no-gain");
        expect(gainVerdict(b, (b.lo + b.hi) / 2)).toBe("uncertain");
        // the verdict agrees with the direct certified comparison of F_B(q) and F_0
        for (const q of [b.lo * 0.5, b.hi * 1.01]) {
          const F = fbEnclosure(p, q)!, F0 = f0Enclosure(eps);
          expect(F[0] > F0[1]).toBe(gainVerdict(b, q) === "gain");
        }
      }
    }
  });
  it("Eq. (70): for q_bi < q <= q_ent the basin threshold eps_c(q) exceeds 2/3 (independent eps_c code)", () => {
    for (const p of PROTOS) {
      const { qBi } = THRESHOLDS[p];
      for (let k = 1; k <= 4; k++) {
        const q = qBi.q + ((Q_ENT[p].q - qBi.q) * k) / 5;
        expect(computeEpsc(p, q)!.lo).toBeGreaterThan(2 / 3);
      }
    }
  });
});

describe("theorem range in eps", () => {
  it("decides 0 < eps < 2/3 exactly, not by display digits", () => {
    expect(epsCoverage(0)).toBe("zero");
    expect(epsCoverage(1e-300)).toBe("covered");
    expect(epsCoverage(0.6666666666666666)).toBe("covered"); //   the float just below 2/3
    expect(epsCoverage(0.6666666666666667)).toBe("not-covered"); // the float just above 2/3
    expect(epsCoverage(0.7)).toBe("not-covered");
    expect(epsCoverage(1)).toBe("not-covered");
    for (const p of PROTOS) {
      expect(qGainBracket(p, 0)).toBeNull();
      expect(qGainBracket(p, 0.7)).toBeNull();
    }
  });
});

describe("gain panel", () => {
  const text = (r: ReturnType<typeof gainReport>) => r.rows.map((x) => `${x.k}: ${x.v}`).join("\n");
  it("keeps the same rows in every state", () => {
    for (const p of PROTOS) {
      const counts = [
        [0.1, 0.01], [0.1, 0.1], [0.1, 0.3], [0, 0.05], [0.7, 0.05], [1, 0.05], [0.5, 0], [0.5, 1],
        [0.6, (THRESHOLDS[p].qBi.q + THRESHOLDS[p].qSN.q) / 2],
      ].map(([e, q]) => gainReport(p, e, q).rows.length);
      expect(new Set(counts).size).toBe(1);
    }
  });
  it("'No' between q_gain and q_ent says P_B, not Z; above q_SN it says Z", () => {
    for (const p of PROTOS) {
      const eps = 0.3;
      const b = qGainBracket(p, eps)!;
      const between = gainReport(p, eps, (b.hi + Q_ENT[p].q) / 2);
      expect(between.verdict).toBe("no-gain");
      expect(text(between)).toContain("P<sub>B</sub>로 수렴하지만");
      expect(regimeReport(p, (b.hi + Q_ENT[p].q) / 2).regime).not.toBe("unique-unpurified");
      const above = gainReport(p, eps, THRESHOLDS[p].qSN.q + 0.01);
      expect(above.verdict).toBe("no-gain");
      expect(text(above)).toContain("Z(q)로 수렴");
      expect(gainReport(p, eps, b.lo / 2).verdict).toBe("gain");
    }
  });
  it("outside the theorem: eps = 0, eps >= 2/3, q = 0, q = 1", () => {
    for (const p of PROTOS) {
      expect(gainReport(p, 0, 0.05).verdict).toBe("not-covered");
      expect(text(gainReport(p, 0.8, 0.05))).toContain("Not covered by the current analytic result");
      expect(gainReport(p, 0.3, 0).verdict).toBe("out-of-range");
      expect(gainReport(p, 0.3, 1).verdict).toBe("out-of-range");
      expect(text(gainReport(p, 0.6666666666666666, 0.05))).toContain("정확히는 2/3보다 작음");
      expect(text(gainReport(p, 0.6666666666666667, 0.05))).toContain("정확히는 2/3보다 큼");
    }
  });
});
