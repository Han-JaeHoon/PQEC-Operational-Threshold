import { describe, expect, it } from "vitest";
import { bellFixedPoint, bistableFixedPoints } from "../src/epsc";
import { fbEnclosure, Q_ENT, qGainBracket, sideOfQEnt } from "../src/qgain";
import {
  dAdd, dMul, dPow, dSub, lambdaOf, normalizedCoeffs, qTargetBracket, resultantN, sOf,
  targetCoverage, targetPolySign, targetVerdict, type Dy,
} from "../src/qtarget";
import { targetReport, type TargetInput } from "../src/targetPanel";
import { regimeReport, THRESHOLDS, type RegimeProtocol } from "../src/regimes";
import { type3Step } from "../src/type3map";
import { type4Step } from "../src/type4map";
import { q4Step } from "../src/q4map";
import ref from "../verification/reference_qtarget.json";

const PROTOS: RegimeProtocol[] = ["type3", "type4", "q4"];
const STEP = { type3: type3Step, type4: type4Step, q4: q4Step } as const;
const mid = (x: readonly [number, number]) => (x[0] + x[1]) / 2;
const F_GRID = Array.from({ length: 49 }, (_, i) => 0.51 + 0.01 * i); // 0.51 ... 0.99

describe("q_target: regression values (q_target note, Section 6 table)", () => {
  // lower ends of certified intervals of width < 1e-12, quoted to 12 digits
  const TABLE: Record<string, [number, number, number]> = {
    "0.6": [0.171847661011, 0.167755801706, 0.120088280938],
    "0.7": [0.154361656727, 0.139524943834, 0.102834729206],
    "0.8": [0.124701910389, 0.101210084626, 0.077190723143],
    "0.9": [0.077335738856, 0.054164752388, 0.042902284650],
    "0.95": [0.043670518431, 0.027863906822, 0.022505129046],
    "0.99": [0.009720149273, 0.005687168622, 0.004665499548],
    "0.999": [0.000997132852, 0.000571160055, 0.000470185871],
    "0.9999": [0.000099971258, 0.000057140175, 0.000047054801],
  };
  it("every table entry agrees with the certified bracket to the 12 quoted digits", () => {
    for (const [F, vals] of Object.entries(TABLE)) {
      (["q4", "type4", "type3"] as const).forEach((p, i) => {
        const b = qTargetBracket(p, Number(F))!;
        expect(b.hi - b.lo).toBeLessThan(1e-12);
        expect(Math.abs(mid([b.lo, b.hi]) - vals[i])).toBeLessThan(1.5e-12);
      });
    }
  });
  it("F_target = 0.90 (task statement)", () => {
    expect(qTargetBracket("type3", 0.9)!.lo).toBeCloseTo(0.042902284650, 11);
    expect(qTargetBracket("type4", 0.9)!.lo).toBeCloseTo(0.054164752388, 11);
    expect(qTargetBracket("q4", 0.9)!.lo).toBeCloseTo(0.077335738856, 11);
  });
});

describe("q_target: independent arb reference (verification/reference_qtarget.json)", () => {
  it("each 1e-30-wide arb interval (Sturm-checked on P_P) lies in the certified float64 bracket", () => {
    expect(ref.rows.length).toBe(39);
    for (const r of ref.rows) {
      const b = qTargetBracket(r.protocol as RegimeProtocol, Number(r.F))!;
      const x = Number(r.lo30);
      const ulp = Math.abs(x) * 2 ** -52;
      expect(b.lo).toBeLessThanOrEqual(x + ulp);
      expect(b.hi).toBeGreaterThanOrEqual(x - ulp);
      expect(r.sturm_roots).toBe(1);
    }
  });
});

describe("q_target: algebraic consistency", () => {
  it("F_B(lo) > F_target > F_B(hi) (certified) and the residual F_B(q_target) - F_target is tiny", () => {
    for (const p of PROTOS) {
      for (const F of [...F_GRID, 0.9999, 0.99999]) {
        const b = qTargetBracket(p, F)!;
        expect(fbEnclosure(p, b.lo)![0]).toBeGreaterThan(F);
        expect(fbEnclosure(p, b.hi)![1]).toBeLessThan(F);
        expect(Math.abs(mid(fbEnclosure(p, mid([b.lo, b.hi]))!) - F)).toBeLessThan(1e-12);
        expect(b.polyCheck).toBe(true);
      }
    }
  });

  it("the root is the physical Bell-connected branch: u_B = (lambda - v_B)/2 > 0, larger root, fixed by the simulator's map", () => {
    for (const p of PROTOS) {
      for (const F of [0.51, 0.6, 0.75, 0.9, 0.99]) {
        const b = qTargetBracket(p, F)!;
        const q = mid([b.lo, b.hi]);
        const pb = bellFixedPoint(p, q)!;
        const [u, v] = [mid(pb.u), mid(pb.v)];
        expect(u).toBeGreaterThan(0);
        expect(v).toBeGreaterThan(0);
        expect(v).toBeLessThanOrEqual(1);
        expect(Math.abs(u - (4 * F - 1 - v) / 2)).toBeLessThan(1e-11);
        // the simulator stores 4Q as (a, b) = (line x, transverse z) and Type 3/4 as (u, v)
        const [x1, y1] = p === "q4" ? STEP.q4(v, u, q) : STEP[p](u, v, q);
        const [u1, v1] = p === "q4" ? [y1, x1] : [x1, y1];
        expect(Math.abs(u1 - u) + Math.abs(v1 - v)).toBeLessThan(1e-12);
        if (regimeReport(p, q).regime === "bistable") {
          const fp = bistableFixedPoints(p, q)!;
          expect(fp.S.v[1]).toBeLessThan(pb.v[0]); // the saddle is the smaller root
          expect(fp.S.F[1]).toBeLessThan(F);
        }
      }
    }
  });

  it("0 < q_target < q_ent, and q_target -> q_ent as F_target -> 1/2", () => {
    for (const p of PROTOS) {
      for (const F of F_GRID) {
        const b = qTargetBracket(p, F)!;
        expect(b.lo).toBeGreaterThan(0);
        expect(sideOfQEnt(p, b.hi)).toBe(-1);
      }
      const near = qTargetBracket(p, 0.5 + 1e-9)!;
      expect(Q_ENT[p].q - near.hi).toBeGreaterThan(0);
      expect(Q_ENT[p].q - near.lo).toBeLessThan(1e-7);
    }
  });

  it("strictly decreasing in F_target: brackets are disjoint and ordered", () => {
    for (const p of PROTOS) {
      let prev = qTargetBracket(p, F_GRID[0])!;
      for (const F of F_GRID.slice(1)) {
        const b = qTargetBracket(p, F)!;
        expect(b.hi).toBeLessThan(prev.lo);
        prev = b;
      }
    }
  });

  it("F_target -> 1: q_target -> 0 with the analytic slope 1/c_P (Section 6.3) and the series of the q_target note", () => {
    const series: Record<RegimeProtocol, [number, number, number]> = {
      q4: [1, -23 / 8, 63 / 8], type4: [4 / 7, -92 / 343, -4952 / 16807], type3: [8 / 17, -1976 / 4913, -235888 / 1419857],
    };
    for (const p of PROTOS) {
      for (const x of [1e-2, 1e-3, 1e-4, 1e-6]) {
        const b = qTargetBracket(p, 1 - x)!;
        const [c1, c2, c3] = series[p];
        expect(Math.abs(mid([b.lo, b.hi]) - (c1 * x + c2 * x * x + c3 * x ** 3))).toBeLessThan(50 * x ** 4 + 1e-14);
      }
    }
  });

  it("Eq. (59): q_target(F) = q_gain(4 (1 - F) / 3), and the order against q_gain(eps)", () => {
    for (const p of PROTOS) {
      for (const F of [0.55, 0.7, 0.85, 0.95]) {
        const t = qTargetBracket(p, F)!;
        const g = qGainBracket(p, (4 * (1 - F)) / 3)!;
        expect(Math.abs(mid([t.lo, t.hi]) - mid([g.lo, g.hi]))).toBeLessThan(1e-12);
      }
      const g = qGainBracket(p, 0.1)!; // F_0 = 0.925
      expect(qTargetBracket(p, 0.95)!.hi).toBeLessThan(g.lo); // F_target > F_0
      expect(qTargetBracket(p, 0.9)!.lo).toBeGreaterThan(g.hi); // F_target < F_0
    }
  });
});

describe("normalized polynomials P_4, P_2, P_4Q (analytic note Eqs. 22-24, 33-34)", () => {
  const I = (n: number): Dy => ({ n: BigInt(n), k: 0 });
  const eq = (a: Dy, b: Dy) => dSub(a, b).n === 0n;
  /** Raw resultant of Eq. (19) from the coefficient table of Section 1. */
  function rawResultant(p: RegimeProtocol, s: Dy, L: Dy): Dy {
    const sp = (k: number) => dPow(s, k);
    const ops = dAdd(I(1), s);
    let al: Dy, be: Dy, ga: Dy, de: Dy, et: Dy, b: Dy;
    if (p === "q4") {
      al = dMul(I(2), sp(5)); be = dMul(sp(4), ops); ga = dMul(I(2), sp(5)); de = dMul(sp(3), ops);
      et = dAdd(sp(3), sp(5)); b = sp(5);
    } else {
      const [A, B] = p === "type3" ? [6, 4] : [4, 2];
      al = be = dMul(I(2), sp(A)); ga = de = dMul(sp(A - 1), ops); et = dMul(I(2), sp(B)); b = sp(B);
    }
    const h2 = dAdd(dMul(ga, b), dMul(et, al));
    const h1 = dSub(dMul(et, dSub(be, de)), dMul(ga, al));
    const h0 = dMul(ga, dSub(I(1), be));
    const e2 = dSub(ga, dMul(I(4), al));
    const e1 = dSub(dMul(dMul(I(-2), ga), L), dMul(I(4), dSub(be, de)));
    const e0 = dMul(ga, dMul(L, L));
    const t = dSub(dMul(h2, e0), dMul(h0, e2));
    return dSub(dMul(t, t), dMul(dSub(dMul(h2, e1), dMul(h1, e2)), dSub(dMul(h1, e0), dMul(h0, e1))));
  }

  it("the raw resultant (Eq. 19) is the normalized one times s^{4(k+1)} (Type 3/4) or 4 s^16 (4Q), exactly", () => {
    const pts = [[0.8125, 1.5], [0.9, 2.2], [0.873, 1.0], [0.99, 2.97], [0.5, 2.5]];
    for (const p of PROTOS) {
      for (const [sv, Lv] of pts) {
        const s = sOf(1 - sv), L = lambdaOf((Lv + 1) / 4);
        const fac = p === "q4" ? dMul(I(4), dPow(s, 16)) : dPow(s, 4 * ((p === "type3" ? 4 : 2) + 1));
        expect(eq(rawResultant(p, s, L), dMul(fac, resultantN(p, s, L)))).toBe(true);
      }
    }
  });

  it("s = 1: P_k(1, L) = -12 (L-3)(L+1)^2(3L-1), P_4Q(1, L) = 12 (...) (Eqs. 51-52)", () => {
    for (const p of PROTOS) {
      for (const Lv of [1.25, 2, 2.5, 2.75]) {
        const L = lambdaOf((Lv + 1) / 4);
        const fac = dMul(dMul(dSub(L, I(3)), dPow(dAdd(L, I(1)), 2)), dSub(dMul(I(3), L), I(1)));
        // P = -N / (1+s)^2 at s = 1 is -N/4; P_4Q = N
        const N = resultantN(p, I(1), L);
        expect(eq(N, p === "q4" ? dMul(I(12), fac) : dMul(I(48), fac))).toBe(true);
      }
    }
  });

  it("lambda = 1: P(s, 1) = (3s^4+1) R_T3, (3s^2+1) R_T4, -R_4Q (Eq. 46), with the q_ent polynomials", () => {
    const horner = (poly: readonly number[], s: Dy) => poly.reduce<Dy>((acc, c) => dAdd(dMul(acc, s), I(c)), I(0));
    for (const p of PROTOS) {
      for (const sv of [0.3, 0.8125, 0.87, 0.95]) {
        const s = sOf(1 - sv);
        const R = horner(Q_ENT[p].poly, s);
        const P = p === "type3" ? dMul(dAdd(dMul(I(3), dPow(s, 4)), I(1)), R)
          : p === "type4" ? dMul(dAdd(dMul(I(3), dPow(s, 2)), I(1)), R) : dMul(I(-1), R);
        const N = resultantN(p, s, I(1));
        const expectN = p === "q4" ? dMul(P, dPow(s, 3)) : dMul(I(-1), dMul(P, dPow(dAdd(I(1), s), 2)));
        expect(eq(N, expectN)).toBe(true);
      }
    }
  });

  it("normalizedCoeffs reproduce the note's c_4 / d_4 leading coefficient", () => {
    // c_4 = -h^2 (5s+1)^2 is the L^4 coefficient of P_k; via N: AG term (A (1+s) L^2)^2 = A^2 (1+s)^2 L^4
    const s = sOf(1 - 0.9);
    const { A } = normalizedCoeffs("type3", s, I(0));
    expect(eq(A, dMul(dPow(s, 4), dAdd(dMul(I(5), s), I(1))))).toBe(true);
  });

  it("targetPolySign changes sign across every certified bracket (exactly one crossing direction)", () => {
    for (const p of PROTOS) {
      for (const F of [0.55, 0.8, 0.95]) {
        const b = qTargetBracket(p, F)!;
        const L = lambdaOf(F);
        expect(targetPolySign(p, sOf(b.hi), L) * targetPolySign(p, sOf(b.lo), L)).toBe(-1);
      }
    }
  });
});

describe("validation and verdicts", () => {
  it("1/2 < F_target < 1 only", () => {
    expect(targetCoverage(0.5)).toBe("too-low");
    expect(targetCoverage(1)).toBe("too-high");
    expect(targetCoverage(0.3)).toBe("too-low");
    expect(targetCoverage(NaN)).toBe("not-a-number");
    expect(targetCoverage(0.5000001)).toBe("valid");
    expect(targetCoverage(0.9999)).toBe("valid");
    expect(qTargetBracket("type3", 1)).toBeNull();
    expect(qTargetBracket("type3", 0.5)).toBeNull();
  });

  it("targetVerdict: <= lo reaches, >= hi misses, inside is near", () => {
    const b = qTargetBracket("type3", 0.9)!;
    expect(targetVerdict(b, b.lo)).toBe("reach");
    expect(targetVerdict(b, b.hi)).toBe("miss");
    expect(targetVerdict(b, mid([b.lo, b.hi]))).toBe("near");
  });

  const base: TargetInput = { protocol: "type3", eps: 0.1, q: 0.03, F: 0.9, attractor: "P_B", n: 30, Fn: 0.9 };
  const answer = (o: Partial<TargetInput>) => targetReport({ ...base, ...o }).answer;

  it("entangled input: yes below, no above, near inside; independent of eps", () => {
    const b = qTargetBracket("type3", 0.9)!;
    for (const eps of [0.01, 0.1, 0.4, 0.66]) {
      expect(answer({ eps, q: 0.03 })).toBe("yes");
      expect(answer({ eps, q: 0.05, attractor: "P_B" })).toBe("no");
      expect(answer({ eps, q: mid([b.lo, b.hi]) })).toBe("near");
      const r = targetReport({ ...base, eps });
      expect([r.lo, r.hi]).toEqual([b.lo, b.hi]);
    }
  });

  it("separable input: through the attractor below q_target, always no above", () => {
    const q = 0.1195; // Type 3, F = 0.6: q_bi < q < q_target(0.6) = 0.12009
    expect(regimeReport("type3", q).regime).toBe("bistable");
    expect(answer({ eps: 0.9, q, F: 0.6, attractor: "P_B" })).toBe("yes");
    expect(answer({ eps: 0.999, q, F: 0.6, attractor: "Z" })).toBe("no");
    expect(answer({ eps: 0.95, q, F: 0.6, attractor: "undetermined" })).toBe("not-guaranteed");
    expect(answer({ eps: 0.9, q: 0.125, F: 0.6, attractor: "P_B" })).toBe("no");
    expect(answer({ eps: 0, q: 0.03 })).toBe("not-covered");
    expect(answer({ q: 0 })).toBe("unresolved");
  });

  it("finite round and asymptotic verdicts are separate: q > q_target with F_n >= F_target", () => {
    const r = targetReport({ ...base, eps: 0.01, q: 0.05, F: 0.95, n: 0, Fn: 0.9925 });
    expect(r.answer).toBe("no");
    const fin = r.rows[r.rows.length - 1];
    expect(fin.v).toContain("Yes");
  });

  it("F_target = F_0: q_target and q_gain are the same threshold; strict gain is not claimed", () => {
    const r = targetReport({ ...base, eps: 0.1, F: 0.925 });
    expect(r.relation).toBe("equal");
    const g = qGainBracket("type3", 0.1)!;
    expect(Math.abs(mid([r.lo!, r.hi!]) - mid([g.lo, g.hi]))).toBeLessThan(1e-12);
    const near = targetReport({ ...base, eps: 0.1, F: 0.925, q: mid([g.lo, g.hi]) });
    expect(near.answer).toBe("near");
    expect(near.rows.find((x) => x.k.startsWith("Asymptotically"))!.v).toContain("strict gain");
  });

  it("q_target lies below q_bi exactly when F_target > F_B(q_bi) (q_target note, Thm 2(d))", () => {
    for (const p of PROTOS) {
      const FBbi = mid(fbEnclosure(p, THRESHOLDS[p].qBi.q)!);
      expect(qTargetBracket(p, FBbi + 1e-6)!.hi).toBeLessThan(THRESHOLDS[p].qBi.q);
      expect(qTargetBracket(p, FBbi - 1e-6)!.lo).toBeGreaterThan(THRESHOLDS[p].qBi.q);
    }
  });
});
