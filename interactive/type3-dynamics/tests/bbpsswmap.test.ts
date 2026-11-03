import { describe, expect, it } from "vitest";
import { allocTrajectory, Step, summarize, type Trajectory } from "../src/type3map";
import {
  bbpsswPSucc, bbpsswRoundOriginal, bbpsswRoundUntwirled, bbpsswStep, computeTrajectoryBBPSSW,
} from "../src/bbpsswmap";

const EPS = [0, 1e-6, 0.01, 0.1, 0.3, 0.5, 2 / 3, 0.8, 0.95, 1];
const QS = [0, 1e-8, 1e-5, 0.001, 0.01, 0.0327, 0.05, 0.08, 0.13, 0.18, 0.3, 0.6, 0.9, 1];

function run(eps: number, q: number, N: number) {
  return computeTrajectoryBBPSSW(allocTrajectory(N), eps, q, N);
}

describe("B, C, D: initial state, normalization, twirling invariance", () => {
  it("F_0 = 1 - 3 eps/4", () => {
    for (const e of EPS) expect(run(e, 0.1, 1).pPhiP[0]).toBe(1 - (3 * e) / 4);
  });

  it("populations sum to 1 and the three non-target populations stay equal", () => {
    for (const e of EPS) for (const q of QS) {
      const t = run(e, q, 60);
      for (let n = 0; n <= 60; n++) {
        expect(Math.abs(t.pPhiP[n] + t.pPhiM[n] + t.pPsiP[n] + t.pPsiM[n] - 1)).toBeLessThan(1e-14);
        expect(t.pPhiM[n]).toBe(t.pPsiP[n]);
        expect(t.pPsiP[n]).toBe(t.pPsiM[n]);
        expect(t.pPhiP[n]).toBeGreaterThanOrEqual(0.25 - 1e-14);
        expect(t.pPhiP[n]).toBeLessThanOrEqual(1 + 1e-14);
      }
    }
  });
});

describe("E, F, G: limits", () => {
  it("E: q = 0 is Bennett's recurrence (F^2 + h^2) / (F^2 + 2Fh + 5h^2)", () => {
    for (let F = 0.25; F <= 1; F += 0.03) {
      const h = (1 - F) / 3;
      const [G, P] = bbpsswStep(F, 0);
      expect(Math.abs(G - (F * F + h * h) / (F * F + 2 * F * h + 5 * h * h))).toBeLessThan(1e-15);
      expect(Math.abs(P - (F * F + 2 * F * h + 5 * h * h))).toBeLessThan(1e-15);
    }
  });

  it("F: q = 1 gives F' = 1/4 and P_succ = 1/2, exactly", () => {
    for (const e of EPS) {
      const t = run(e, 1, 5);
      for (let n = 1; n <= 5; n++) expect(t.pPhiP[n]).toBe(0.25);
      for (let n = 0; n <= 5; n++) expect(bbpsswPSucc(t.pPhiP[n], 1)).toBe(0.5);
    }
  });

  it("G: eps = 1 stays at F = 1/4 for every q", () => {
    for (const q of QS) {
      const t = run(1, q, 40);
      for (let n = 0; n <= 40; n++) expect(t.pPhiP[n]).toBe(0.25);
    }
  });

  it("P_succ lies in [1/2, 1]", () => {
    for (const e of EPS) for (const q of QS) {
      const P = bbpsswPSucc(1 - (3 * e) / 4, q);
      expect(P).toBeGreaterThanOrEqual(0.5 - 1e-15);
      expect(P).toBeLessThanOrEqual(1 + 1e-15);
    }
  });

  it("the trajectory's F agrees with bbpsswStep and the closed eps form of round 1", () => {
    for (const e of EPS) for (const q of QS) {
      const t = run(e, q, 25);
      let F = t.pPhiP[0];
      for (let n = 1; n <= 25; n++) {
        F = bbpsswStep(F, q)[0];
        expect(t.pPhiP[n]).toBe(F);
      }
      const s = (1 - q) ** 2;
      const A = 1 - 1.5 * e + (5 * e * e) / 8;
      const P0 = 1 - e + (e * e) / 2;
      expect(Math.abs(t.pPhiP[1] - (s * A + (1 - s) / 8) / (s * P0 + (1 - s) / 2))).toBeLessThan(1e-14);
    }
  });
});

describe("I: input twirl is part of the protocol", () => {
  const p = [0.67, 0.03, 0.19, 0.11] as const;

  it("the original map depends on the input through its fidelity only", () => {
    const [a] = bbpsswRoundOriginal(p, 0.1);
    const h = (1 - p[0]) / 3;
    const [b] = bbpsswRoundOriginal([p[0], h, h, h], 0.1);
    expect(a).toEqual(b);
  });

  it("omitting the input twirl gives a different fidelity on a non-isotropic input", () => {
    const original = bbpsswRoundOriginal(p, 0.1)[0][0];
    const noInput = bbpsswRoundUntwirled(p, 0.1)[0][0];
    expect(Math.abs(original - noInput)).toBeGreaterThan(1e-3);
  });

  it("on an isotropic input the two coincide, and the scalar recurrence is that common value", () => {
    for (const e of EPS) for (const q of QS) {
      const F = 1 - (3 * e) / 4;
      const h = (1 - F) / 3;
      const [u, Pu] = bbpsswRoundUntwirled([F, h, h, h], q);
      const [G, P] = bbpsswStep(F, q);
      expect(Math.abs(u[0] - G)).toBeLessThan(1e-14);
      expect(Math.abs(Pu - P)).toBeLessThan(1e-14);
    }
  });
});

describe("first-round gain sign", () => {
  it("F_1 - F_0 changes sign at 1 - sqrt(3/(3 + 2 eps - 3 eps^2)), used here only as a check", () => {
    for (const e of [0.01, 0.1, 0.2, 0.4, 0.6]) {
      const q0 = 1 - Math.sqrt(3 / (3 + 2 * e - 3 * e * e));
      expect(run(e, q0 * 0.8, 1).dF[0]).toBeGreaterThan(0);
      expect(run(e, q0 * 1.2, 1).dF[0]).toBeLessThan(0);
    }
  });
});

// Exact arithmetic on dyadic rationals m * 2^e: 24 F_{n+1} - 24 F_n is a ratio of polynomials in
// F and r = 1 - q, computed EXACTLY from the stored float F_n and the exact r.
type Dy = { m: bigint; e: number };

function dy(x: number): Dy {
  if (x === 0) return { m: 0n, e: 0 };
  const buf = new DataView(new ArrayBuffer(8));
  buf.setFloat64(0, x);
  const hi = buf.getUint32(0), lo = buf.getUint32(4);
  const sign = hi >>> 31 ? -1n : 1n;
  const exp = (hi >>> 20) & 0x7ff;
  let mant = (BigInt(hi & 0xfffff) << 32n) | BigInt(lo);
  let e: number;
  if (exp === 0) e = -1074; else { mant |= 1n << 52n; e = exp - 1075; }
  return { m: sign * mant, e };
}
const add = (a: Dy, b: Dy): Dy => {
  if (a.e > b.e) return { m: (a.m << BigInt(a.e - b.e)) + b.m, e: b.e };
  return { m: a.m + (b.m << BigInt(b.e - a.e)), e: a.e };
};
const neg = (a: Dy): Dy => ({ m: -a.m, e: a.e });
const sub = (a: Dy, b: Dy) => add(a, neg(b));
const mul = (a: Dy, b: Dy): Dy => ({ m: a.m * b.m, e: a.e + b.e });
const absd = (a: Dy): Dy => ({ m: a.m < 0n ? -a.m : a.m, e: a.e });
const le = (a: Dy, b: Dy) => sub(b, a).m >= 0n;
const int = (k: number) => dy(k);

/** dF = Num / Den exactly, from the direct definition F' - F (not the factored route 2). */
function exactDrift(F: number, q: number): { N: Dy; D: Dy } {
  const f = dy(F);
  const r = sub(int(1), dy(q));
  const s = mul(r, r);
  const c = sub(int(1), s);
  // with H = 3h = 1 - F, multiply through by 72:
  const H = sub(int(1), f);
  const FF = mul(f, f), HH = mul(H, H);
  // 72 w_q^Phi+ = 8 s (9F^2 + H^2) + 9 c
  const num = add(mul(mul(int(8), s), add(mul(int(9), FF), HH)), mul(int(9), c));
  // 72 P_succ = 8 s (9F^2 + 6FH + 5H^2) + 36 c
  const fh = mul(f, H);
  const den = add(
    mul(mul(int(8), s), add(add(mul(int(9), FF), mul(int(6), fh)), mul(int(5), HH))),
    mul(int(36), c),
  );
  return { N: sub(num, mul(f, den)), D: den };
}

function checkBounds(t: Trajectory) {
  let bad1 = 0, bad2 = 0, wrong = 0, resolved = 0;
  for (let n = 0; n <= t.N; n++) {
    const { N, D } = exactDrift(t.pPhiP[n], t.q);
    if (!le(absd(sub(mul(dy(t.dF[n]), D), N)), mul(dy(t.dFerr[n]), D))) bad1++;
    if (!le(absd(sub(mul(dy(t.dFmap[n]), D), N)), mul(dy(t.dFmapErr[n]), D))) bad2++;
    const st = t.step[n];
    if (st !== Step.Unresolved) {
      resolved++;
      const sgn = N.m > 0n ? 1 : N.m < 0n ? -1 : 0;
      if (sgn !== st) wrong++;
    }
  }
  return { bad1, bad2, wrong, resolved };
}

describe("numerical sign: dF error bounds against exact rational arithmetic", () => {
  it("both bounds hold and every resolved sign is exactly right", () => {
    let resolved = 0, total = 0;
    for (const e of [...EPS, 0.0443, 0.2, 0.4, 0.6, 0.65, 0.99, 0.999999]) {
      for (const q of [...QS, 0.02, 0.0433, 0.06, 0.07, 0.1, 0.2, 0.5, 0.99]) {
        const t = run(e, q, 200);
        const r = checkBounds(t);
        expect(r.bad1).toBe(0);
        expect(r.bad2).toBe(0);
        expect(r.wrong).toBe(0);
        resolved += r.resolved;
        total += t.N + 1;
      }
    }
    console.log(`resolved ${resolved} of ${total} steps`);
    expect(resolved).toBeGreaterThan(total / 10);
  });

  it("the bounds are not vacuous: a large drift is resolved", () => {
    const t = run(0.5, 0.01, 5);
    expect(t.step[0]).toBe(Step.Increase);
    expect(run(0.3, 0.3, 5).step[0]).toBe(Step.Decrease);
  });

  it("a fixed point is left unresolved rather than given a sign", () => {
    expect(run(1, 0.1, 3).step[0]).toBe(Step.Unresolved);
    expect(run(0, 0, 3).step[0]).toBe(Step.Unresolved);
  });

  it("summarize reports gain inputs and the first decrease", () => {
    const s = summarize(run(0.3, 0.2, 30));
    expect(s.firstDecrease).toBe(0);
    expect(s.FN).toBeLessThan(s.F0);
    expect(s.maxFRound).toBe(0);
  });
});
