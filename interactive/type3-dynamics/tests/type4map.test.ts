import { describe, expect, it } from "vitest";
import {
  allocTrajectory, bellPopulations, computeTrajectory, EPS64, Step, summarize, type Trajectory,
} from "../src/type3map";
import { computeTrajectory4, type4Step } from "../src/type4map";

const EPS_GRID = Array.from({ length: 101 }, (_, i) => i / 100);
const Q_GRID = [0, 1e-4, 0.001, 0.01, 0.0405, 0.05, 0.1, 0.15, 0.2, 0.5, 0.9, 1];

function run4(eps: number, q: number, N: number) {
  return computeTrajectory4(allocTrajectory(N), eps, q, N);
}
function run3(eps: number, q: number, N: number) {
  return computeTrajectory(allocTrajectory(N), eps, q, N);
}

describe("Type 4 — Test A: initial state is the same rho_0(eps) as Type 3", () => {
  it("u0 = v0 = 1 - eps, F0 = 1 - 3 eps / 4, and round 0 is bitwise identical to Type 3", () => {
    for (const eps of [...EPS_GRID, 0.0327, 2 / 3, 1e-9]) {
      const t4 = run4(eps, 0.05, 1);
      const t3 = run3(eps, 0.05, 1);
      expect(t4.u[0]).toBe(1 - eps);
      expect(t4.v[0]).toBe(1 - eps);
      expect(t4.pPhiP[0]).toBeCloseTo(1 - (3 * eps) / 4, 15);
      expect(t4.pPhiM[0]).toBeCloseTo(eps / 4, 15);
      expect(t4.pPsiP[0]).toBeCloseTo(eps / 4, 15);
      expect([t4.pPhiP[0], t4.pPhiM[0], t4.pPsiP[0], t4.pPsiM[0]])
        .toEqual([t3.pPhiP[0], t3.pPhiM[0], t3.pPsiP[0], t3.pPsiM[0]]);
    }
  });
});

describe("Type 4 — Tests B, C, D: normalisation, positivity, invariant plane", () => {
  it("holds for every round on an (eps, q) grid, N = 500", () => {
    let worstSum = 0, worstNeg = 0;
    for (const eps of EPS_GRID.filter((_, i) => i % 5 === 0)) {
      for (const q of Q_GRID) {
        const t = run4(eps, q, 500);
        for (let n = 0; n <= 501; n++) {
          const s = t.pPhiP[n] + t.pPhiM[n] + t.pPsiP[n] + t.pPsiM[n];
          worstSum = Math.max(worstSum, Math.abs(s - 1));
          worstNeg = Math.min(worstNeg, t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]);
          expect(t.pPsiP[n]).toBe(t.pPsiM[n]); // D: exact equality on the plane
          expect(Number.isFinite(t.u[n]) && Number.isFinite(t.v[n])).toBe(true);
        }
      }
    }
    expect(worstSum).toBeLessThan(1e-15); // B
    expect(worstNeg).toBeGreaterThanOrEqual(-1e-15); // C
  });

  it("the general (x, y, z) map keeps y = -x and reduces to the (u, v) step", () => {
    // Pauli form of the Type-4 map (Bell-state-distillation analytic_exact_map, k = 2)
    const pauli = (x: number, y: number, z: number, q: number) => {
      const r = 1 - q;
      const n = 1 + r ** 2 * (x * x + y * y + z * z);
      return [2 * r ** 4 * (x - y * z) / n, 2 * r ** 4 * (y - x * z) / n, r ** 3 * (1 + r) * (z - x * y) / n, n];
    };
    for (const u of [-0.9, -0.3, 0, 0.2, 0.7, 1]) for (const v of [-1, -0.4, 0, 0.5, 0.95]) for (const q of Q_GRID) {
      const [xp, yp, zp, n] = pauli(u, -u, v, q);
      const [up, vp, D] = type4Step(u, v, q);
      expect(yp + xp).toBe(0); // exact; (toBe(-xp) would distinguish +0 from -0)
      expect(Math.abs(xp - up)).toBeLessThan(1e-15);
      expect(Math.abs(zp - vp)).toBeLessThan(1e-15);
      expect(Math.abs(n - D)).toBeLessThan(1e-15);
    }
  });
});

describe("Type 4 — Test E: ideal case q = 0", () => {
  it("reproduces rho -> rho^2 / Tr rho^2 (Bell populations p -> p^2 / sum p^2)", () => {
    for (const eps of [0.05, 0.1, 0.3, 0.5, 2 / 3, 0.8, 0.95]) {
      const t = run4(eps, 0, 40);
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
      expect(Math.abs(run4(eps, 0, 1).pPhiP[1] - F1)).toBeLessThan(1e-15);
    }
  });
});

describe("Type 4 — Test F: maximally mixed input", () => {
  it("eps = 1 keeps F_n = 1/4 exactly for every q", () => {
    for (const q of Q_GRID) {
      const t = run4(1, q, 500);
      for (let n = 0; n <= 501; n++) expect(t.pPhiP[n]).toBe(0.25);
    }
  });
});

describe("Type 4 is a different map from Type 3", () => {
  it("differs from Type 3 for 0 < q < 1 on an entangled input, and coincides at q = 0", () => {
    for (const q of [0.001, 0.0405, 0.1, 0.5]) {
      expect(run4(0.1, q, 1).pPhiP[1]).not.toBe(run3(0.1, q, 1).pPhiP[1]);
    }
    const a = run4(0.2, 0, 20), b = run3(0.2, 0, 20);
    for (let n = 0; n <= 20; n++) expect(Math.abs(a.pPhiP[n] - b.pPhiP[n])).toBeLessThan(1e-15);
  });
  it("first-round closed form u1, v1 of the notes (Eqs. 39-40)", () => {
    for (const eps of EPS_GRID) for (const q of Q_GRID) {
      const t = 1 - eps, r = 1 - q, D0 = 1 + 3 * r * r * t * t;
      const tr = run4(eps, q, 1);
      // two evaluation orders of the same rational expression: a couple of ulps apart at most
      expect(Math.abs(tr.u[1] - (2 * r ** 4 * t * (1 + t)) / D0)).toBeLessThan(5e-16);
      expect(Math.abs(tr.v[1] - (r ** 3 * (1 + r) * t * (1 + t)) / D0)).toBeLessThan(5e-16);
    }
  });
});

// ---------------------------------------------------------------------------------------
// Exact check of the dF error bounds.  Every float is a dyadic rational m * 2^e; the drift
// numerator N_n = 2u(a - D) + b - vD and the denominator 4D are polynomials in u, v, r, so
// they are computed EXACTLY in BigInt arithmetic from the float inputs (u_n, v_n) and the
// exact r = 1 - q.  Then F_{n+1} - F_n = N/(4D) exactly, and each float route must satisfy
// |route - N/(4D)| <= its bound, i.e. |route * 4D - N| <= bound * 4D, compared exactly.
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

/** Exact N and 4D for one step from float (u, v) with exact r = 1 - q; kPow = [p_D, p_a, p_b]. */
function exactDrift(u: number, v: number, q: number, k: number): { N: Dy; D4: Dy } {
  const U = dy(u), V = dy(v), r = sub(ONE, dy(q));
  const pw = (x: Dy, n: number) => { let o = ONE; for (let i = 0; i < n; i++) o = mul(o, x); return o; };
  // Type 3: k = 4 (D: r^4, a: r^6, b: r^5); Type 4: k = 2 (D: r^2, a: r^4, b: r^3)
  const D = add(ONE, mul(pw(r, k), add(mul(TWO, mul(U, U)), mul(V, V))));
  const a = mul(TWO, mul(pw(r, k + 2), add(ONE, V)));
  const b = mul(mul(pw(r, k + 1), add(ONE, r)), add(mul(U, U), V));
  const N = sub(add(mul(mul(TWO, U), sub(a, D)), b), mul(V, D));
  return { N, D4: mul(FOUR, D) };
}

function checkBounds(t: Trajectory, k: number) {
  let worst1 = 0, worst2 = 0, resolvedWrong = 0;
  for (let n = 0; n <= t.N; n++) {
    const { N, D4 } = exactDrift(t.u[n], t.v[n], t.q, k);
    const r1 = absd(sub(mul(dy(t.dF[n]), D4), N));
    const r2 = absd(sub(mul(dy(t.dFmap[n]), D4), N));
    const ok1 = le(r1, mul(dy(t.dFerr[n]), D4));
    const ok2 = le(r2, mul(dy(t.dFmapErr[n]), D4));
    if (!ok1) worst1++;
    if (!ok2) worst2++;
    // a resolved sign must be the sign of the exact value
    const s = t.step[n];
    if (s !== Step.Unresolved && Math.sign(Number(N.m)) !== s) resolvedWrong++;
  }
  return { worst1, worst2, resolvedWrong };
}

describe("Type 4 — numerical sign: dF error bounds against exact rational arithmetic", () => {
  // N = 500 and q up to 0.99 drive the state into float64 underflow (u, v < 1e-308), where
  // the absolute (2^-1074) part of the route-2 bound is what keeps it valid.
  const cases: [number, number, number][] = [];
  for (const eps of [0, 0.02, 0.1, 0.15, 0.3, 0.5, 2 / 3, 0.9, 1]) {
    for (const q of [0, 1e-4, 0.0101, 0.0405, 0.04781, 0.04898, 0.059, 0.1, 0.2, 0.5, 0.9, 0.95, 0.99]) {
      cases.push([eps, q, 500]);
    }
  }
  it("both routes stay within their bounds and every resolved sign is the exact sign", () => {
    let bad = 0;
    for (const [eps, q, N] of cases) {
      const r = checkBounds(run4(eps, q, N), 2);
      bad += r.worst1 + r.worst2 + r.resolvedWrong;
      if (r.worst1 + r.worst2 + r.resolvedWrong) console.error("Type 4 bound violated", eps, q, r);
    }
    expect(bad).toBe(0);
  });
  it("the same exact check passes for the unchanged Type-3 engine (regression, no change to it)", () => {
    let bad = 0;
    for (const [eps, q, N] of cases) {
      const r = checkBounds(run3(eps, q, N), 4);
      bad += r.worst1 + r.worst2 + r.resolvedWrong;
      if (r.worst1 + r.worst2 + r.resolvedWrong) console.error("Type 3 bound violated", eps, q, r);
    }
    expect(bad).toBe(0);
  });
  it("near the fixed point steps become unresolved instead of getting a spurious sign", () => {
    const t = run4(0.1, 0.0101, 400);
    let unres = 0;
    for (let n = 0; n < 400; n++) if (t.step[n] === Step.Unresolved) unres++;
    expect(unres).toBeGreaterThan(0);
    expect(EPS64).toBe(2 ** -53);
  });
});

describe("Type 4 — drift routes agree within the combined bounds (as for Type 3)", () => {
  it("on an (eps, q) grid, N = 300", () => {
    for (const eps of EPS_GRID.filter((_, i) => i % 4 === 0)) for (const q of Q_GRID) {
      const t = run4(eps, q, 300);
      for (let n = 0; n <= 300; n++) {
        expect(Math.abs(t.dF[n] - t.dFmap[n])).toBeLessThanOrEqual(t.dFerr[n] + t.dFmapErr[n]);
      }
    }
  });
});

describe("Type 4 — rise, then later-round decrease (eps = 0.1, q = 0.0405)", () => {
  it("first step increases, a later step is a resolved decrease far above the bound", () => {
    const t = run4(0.1, 0.0405, 30);
    expect(t.step[0]).toBe(Step.Increase);
    const s = summarize(t);
    expect(s.firstDecrease).toBe(2);
    expect(-t.dF[2]).toBeGreaterThan(1e6 * t.dFerr[2]);
  });
  it("bellPopulations matches the trajectory buffers", () => {
    const t = run4(0.37, 0.021, 10);
    for (let n = 0; n <= 10; n++) {
      expect(bellPopulations(t.u[n], t.v[n])).toEqual([t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]]);
    }
  });
});
