/**
 * Test G — comparison against independent reference trajectories.
 *
 * verification/reference_type3.json is produced by verification/generate_reference.py from
 * two full-state routes that never use the (u, v) closed form:
 *   R1  dense 32x32 propagation through the noisy 16-CNOT circuit (iterated_noisy_pqec.py)
 *   R2  PennyLane default.mixed with the repository's Step-3 executor; next state rebuilt
 *       from all 16 parity-weighted Pauli correlators
 * plus R3, the repository's one-round closed form F_dec.
 */
import { describe, expect, it } from "vitest";
import ref from "../verification/reference_type3.json";
import { allocTrajectory, computeTrajectory } from "../src/type3map";

type Case = {
  eps: number; q: number; N: number;
  R1_dense: number[][]; R2_pennylane: number[][];
  R3_F1_closed_form: number;
};

const cases = (ref as { cases: Case[] }).cases;
const TOL = 1e-12;

describe("Test G — simulator vs. independent full-circuit references", () => {
  it("has reference cases", () => expect(cases.length).toBeGreaterThanOrEqual(10));

  for (const c of cases) {
    it(`eps=${c.eps}, q=${c.q}: Bell populations for n = 0..${c.N}`, () => {
      const t = computeTrajectory(allocTrajectory(c.N), c.eps, c.q, c.N);
      let w1 = 0, w2 = 0;
      for (let n = 0; n <= c.N; n++) {
        const js = [t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]];
        for (let k = 0; k < 4; k++) {
          w1 = Math.max(w1, Math.abs(js[k] - c.R1_dense[n][k]));
          w2 = Math.max(w2, Math.abs(js[k] - c.R2_pennylane[n][k]));
        }
      }
      expect(w1).toBeLessThan(TOL);
      expect(w2).toBeLessThan(TOL);
      expect(Math.abs(t.pPhiP[1] - c.R3_F1_closed_form)).toBeLessThan(TOL);
    });
  }

  it("eps=0.1, q=0.0327: the references themselves show the rise and the later decrease", () => {
    const c = cases.find((x) => x.eps === 0.1 && x.q === 0.0327)!;
    const F = c.R1_dense.map((p) => p[0]);
    expect(F[1]).toBeGreaterThan(F[0]);
    const firstDec = F.findIndex((f, n) => n + 1 < F.length && F[n + 1] < f - 1e-12);
    expect(firstDec).toBeGreaterThan(0);
    const js = computeTrajectory(allocTrajectory(c.N), c.eps, c.q, c.N);
    let jsFirst = -1;
    for (let n = 0; n < c.N; n++) if (js.step[n] === -1) { jsFirst = n; break; }
    expect(jsFirst).toBe(firstDec);
  });
});
