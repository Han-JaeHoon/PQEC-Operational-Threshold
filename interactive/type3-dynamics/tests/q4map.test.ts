import { describe, expect, it } from "vitest";
import {
  allocTrajectory, computeTrajectory, Step, summarize, type Trajectory,
} from "../src/type3map";
import { computeTrajectory4 } from "../src/type4map";
import { bellPopulations4Q, computeTrajectory4Q, q4Step } from "../src/q4map";

const EPS_GRID = Array.from({ length: 101 }, (_, i) => i / 100);
const Q_GRID = [0, 1e-4, 0.001, 0.01, 0.0438, 0.0606, 0.0628, 0.1, 0.15, 0.2, 0.5, 0.9, 1];

function run(eps: number, q: number, N: number) {
  return computeTrajectory4Q(allocTrajectory(N), eps, q, N);
}

describe("4Q — Test A: initial state is the same rho_0(eps) as Type 3 / Type 4", () => {
  it("a0 = b0 = 1 - eps, F0 = 1 - 3 eps / 4, all non-target populations eps / 4", () => {
    for (const eps of [...EPS_GRID, 0.0327, 2 / 3, 1e-9]) {
      const t = run(eps, 0.05, 1);
      const t3 = computeTrajectory(allocTrajectory(1), eps, 0.05, 1);
      expect(t.u[0]).toBe(1 - eps);
      expect(t.v[0]).toBe(1 - eps);
      expect(t.pPhiP[0]).toBeCloseTo(1 - (3 * eps) / 4, 15);
      for (const p of [t.pPhiM[0], t.pPsiP[0], t.pPsiM[0]]) expect(p).toBeCloseTo(eps / 4, 15);
      // same state as Type 3; (1 + a + 2b) vs (1 + 2u + v) may round differently by an ulp
      expect(Math.abs(t.pPhiP[0] - t3.pPhiP[0])).toBeLessThan(3e-16);
    }
  });
});

describe("4Q — Tests B, C, D: normalisation, positivity, invariant plane", () => {
  it("holds for every round on an (eps, q) grid, N = 500", () => {
    let worstSum = 0, worstNeg = 0;
    for (const eps of EPS_GRID.filter((_, i) => i % 5 === 0)) {
      for (const q of Q_GRID) {
        const t = run(eps, q, 500);
        for (let n = 0; n <= 501; n++) {
          const s = t.pPhiP[n] + t.pPhiM[n] + t.pPsiP[n] + t.pPsiM[n];
          worstSum = Math.max(worstSum, Math.abs(s - 1));
          worstNeg = Math.min(worstNeg, t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]);
          expect(t.pPhiM[n]).toBe(t.pPsiM[n]); // D: exact equality on the 4Q plane
          // a, b >= 0 along the whole trajectory (assumed by the F_{n+1} error bound)
          expect(t.u[n]).toBeGreaterThanOrEqual(0);
          expect(t.v[n]).toBeGreaterThanOrEqual(0);
        }
      }
    }
    expect(worstSum).toBeLessThan(1e-15); // B
    expect(worstNeg).toBeGreaterThanOrEqual(-1e-15); // C
  });

  it("the general (x, y, z) map keeps y = -z and reduces to the (a, b) step", () => {
    // Pauli form of the 4Q map (Bell-state-distillation analytic_exact_map.analytic_map_4q)
    const pauli = (x: number, y: number, z: number, q: number) => {
      const r = 1 - q;
      const n = 1 + r ** 5 * (x * x + y * y) + r ** 3 * z * z;
      return [
        r ** 3 * ((1 + r) * x - 2 * r ** 2 * y * z) / n,
        r ** 4 * ((1 + r) * y - 2 * r * x * z) / n,
        r ** 4 * ((1 + r) * z - 2 * r * x * y) / n, n,
      ];
    };
    for (const a of [-0.9, -0.3, 0, 0.2, 0.7, 1]) for (const b of [-1, -0.4, 0, 0.5, 0.95]) for (const q of Q_GRID) {
      const [xp, yp, zp, n] = pauli(a, -b, b, q);
      const [ap, bp, D] = q4Step(a, b, q);
      expect(Math.abs(yp + zp)).toBeLessThan(1e-16); // two orders of the same product
      expect(Math.abs(xp - ap)).toBeLessThan(1e-15);
      expect(Math.abs(zp - bp)).toBeLessThan(1e-15);
      expect(Math.abs(n - D)).toBeLessThan(1e-15);
    }
  });

  it("this is not the Type-3/4 plane: Phi- = Psi- while Psi+ differs", () => {
    const t = run(0.2, 0.03, 5);
    for (let n = 1; n <= 5; n++) expect(Math.abs(t.pPsiP[n] - t.pPsiM[n])).toBeGreaterThan(1e-4);
  });
});

describe("4Q — Test E: ideal case q = 0", () => {
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
  it("one ideal round matches F1 = (1+3t)^2 / (4(1+3t^2)) (same as Type 3 / Type 4)", () => {
    for (const eps of EPS_GRID) {
      const t = 1 - eps;
      const F1 = ((1 + 3 * t) ** 2) / (4 * (1 + 3 * t * t));
      expect(Math.abs(run(eps, 0, 1).pPhiP[1] - F1)).toBeLessThan(1e-15);
    }
  });
});

describe("4Q — Test F: maximally mixed input and full noise", () => {
  it("eps = 1 keeps F_n = 1/4 exactly for every q", () => {
    for (const q of Q_GRID) {
      const t = run(1, q, 500);
      for (let n = 0; n <= 501; n++) expect(t.pPhiP[n]).toBe(0.25);
    }
  });
  it("q = 1 sends every input to I/4 in one round, with P_succ = D/4 = 1/4", () => {
    for (const eps of EPS_GRID) {
      const t = run(eps, 1, 3);
      for (let n = 1; n <= 4; n++) for (const p of [t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]]) expect(p).toBe(0.25);
      expect(q4Step(1 - eps, 1 - eps, 1)[2]).toBe(1);
    }
  });
});

describe("4Q is a different map from Type 3 and Type 4", () => {
  it("differs for 0 < q < 1 on an entangled input, and coincides at q = 0", () => {
    for (const q of [0.001, 0.0438, 0.1, 0.5]) {
      const f = run(0.1, q, 1).pPhiP[1];
      expect(f).not.toBe(computeTrajectory(allocTrajectory(1), 0.1, q, 1).pPhiP[1]);
      expect(f).not.toBe(computeTrajectory4(allocTrajectory(1), 0.1, q, 1).pPhiP[1]);
    }
    const a = run(0.2, 0, 20), b = computeTrajectory(allocTrajectory(20), 0.2, 0, 20);
    for (let n = 0; n <= 20; n++) expect(Math.abs(a.pPhiP[n] - b.pPhiP[n])).toBeLessThan(1e-15);
  });
  it("first-round closed forms of the q_abs note (Eqs. 48, 50, 53)", () => {
    for (const eps of EPS_GRID) for (const q of Q_GRID) {
      const e = 1 - eps, r = 1 - q, D0 = 1 + r ** 3 * (1 + 2 * r * r) * e * e;
      const tr = run(eps, q, 1);
      // two evaluation orders of the same rational expression (values <= 1): a few ulps apart
      expect(Math.abs(tr.u[1] - (r ** 3 * e * (1 + r + 2 * r * r * e)) / D0)).toBeLessThan(1e-15);
      expect(Math.abs(tr.v[1] - (r ** 4 * (1 + r) * e + 2 * r ** 5 * e * e) / D0)).toBeLessThan(1e-15);
    }
  });
  it("bellPopulations4Q matches the trajectory buffers", () => {
    const t = run(0.37, 0.021, 10);
    for (let n = 0; n <= 10; n++) {
      expect(bellPopulations4Q(t.u[n], t.v[n])).toEqual([t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]]);
    }
  });
});

// ---------------------------------------------------------------------------------------
// Exact check of the dF error bounds (as in type4map.test.ts).  Every float is a dyadic
// rational m * 2^e; N_n = A + 2B - (a + 2b) D and 4D are polynomials in a, b, r, computed
// EXACTLY in BigInt arithmetic from the float inputs (a_n, b_n) and the exact r = 1 - q.
// Then F_{n+1} - F_n = N/(4D) exactly, and each float route must satisfy
// |route * 4D - N| <= bound * 4D, compared exactly.
// ---------------------------------------------------------------------------------------
type Dy = { m: bigint; e: number }; // m * 2^e

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
const ONE = dy(1), TWO = dy(2), FOUR = dy(4);

function exactDrift4Q(a: number, b: number, q: number): { N: Dy; D4: Dy } {
  const A_ = dy(a), B_ = dy(b), r = sub(ONE, dy(q));
  const r2 = mul(r, r), r3 = mul(r2, r), r4 = mul(r3, r), r5 = mul(r4, r);
  const bb = mul(B_, B_);
  const D = add(add(ONE, mul(r5, mul(A_, A_))), mul(add(r5, r3), bb));
  const A = add(mul(mul(r3, add(ONE, r)), A_), mul(TWO, mul(r5, bb)));
  const B = mul(mul(r4, B_), add(add(ONE, r), mul(TWO, mul(r, A_))));
  const N = sub(add(A, mul(TWO, B)), mul(add(A_, mul(TWO, B_)), D));
  return { N, D4: mul(FOUR, D) };
}

function checkBounds(t: Trajectory) {
  let bad1 = 0, bad2 = 0, resolvedWrong = 0;
  for (let n = 0; n <= t.N; n++) {
    const { N, D4 } = exactDrift4Q(t.u[n], t.v[n], t.q);
    if (!le(absd(sub(mul(dy(t.dF[n]), D4), N)), mul(dy(t.dFerr[n]), D4))) bad1++;
    if (!le(absd(sub(mul(dy(t.dFmap[n]), D4), N)), mul(dy(t.dFmapErr[n]), D4))) bad2++;
    const s = t.step[n];
    if (s !== Step.Unresolved && Math.sign(Number(N.m)) !== s) resolvedWrong++;
  }
  return { bad1, bad2, resolvedWrong };
}

describe("4Q — numerical sign: dF error bounds against exact rational arithmetic", () => {
  // N = 500 and q up to 0.99 drive a, b into float64 underflow, where the absolute
  // (2^-1074) part of the route-2 bound is what keeps it valid.
  const cases: [number, number][] = [];
  for (const eps of [0, 0.02, 0.1, 0.15, 0.3, 0.5, 2 / 3, 0.75, 0.9, 1]) {
    for (const q of [0, 1e-4, 0.01, 0.0433, 0.0438, 0.0443, 0.0606, 0.0612, 0.0628, 0.1, 0.16, 0.2, 0.5, 0.9, 0.99, 1]) {
      cases.push([eps, q]);
    }
  }
  it("both routes stay within their bounds and every resolved sign is the exact sign", () => {
    let bad = 0;
    for (const [eps, q] of cases) {
      const r = checkBounds(run(eps, q, 500));
      const k = r.bad1 + r.bad2 + r.resolvedWrong;
      bad += k;
      if (k) console.error("4Q bound violated", eps, q, r);
    }
    expect(bad).toBe(0);
  });
  it("near the fixed point steps become unresolved instead of getting a spurious sign", () => {
    const t = run(0.1, 0.01, 400);
    let unres = 0;
    for (let n = 0; n < 400; n++) if (t.step[n] === Step.Unresolved) unres++;
    expect(unres).toBeGreaterThan(0);
  });
});

describe("4Q — drift routes agree within the combined bounds", () => {
  it("on an (eps, q) grid, N = 300", () => {
    for (const eps of EPS_GRID.filter((_, i) => i % 4 === 0)) for (const q of Q_GRID) {
      const t = run(eps, q, 300);
      for (let n = 0; n <= 300; n++) {
        expect(Math.abs(t.dF[n] - t.dFmap[n])).toBeLessThanOrEqual(t.dFerr[n] + t.dFmapErr[n]);
      }
    }
  });
});

describe("4Q — fidelity monotonicity on representative points", () => {
  it("eps = 0.1, q = 0.0438 (between q_abs and q_rel): no decrease over 500 rounds", () => {
    const s = summarize(run(0.1, 0.0438, 500));
    expect(s.firstDecrease).toBe(-1);
  });
  it("eps = 0.1, q = 0.0615: F decreases from the first step, far above the bound", () => {
    const t = run(0.1, 0.0615, 30);
    expect(t.step[0]).toBe(Step.Decrease);
    expect(-t.dF[0]).toBeGreaterThan(1e6 * t.dFerr[0]);
  });
});
