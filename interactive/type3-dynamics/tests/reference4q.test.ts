/**
 * Test G4Q — 4Q mode against independent reference trajectories.
 *
 * verification/reference_4q.json is produced by verification/generate_reference_4q.py from
 * two full-state routes through the noisy four-qubit 5-CNOT circuit that never use the (a, b)
 * closed form:
 *   R1  Bell-state-distillation's dense 16x16 simulator (noisy_full_channel_dm + postselect_branch)
 *   R2  PennyLane default.mixed written from the gate list only; next state rebuilt from all
 *       16 expectations <P (x) |00><00|_{q3 q4}>
 * Both record the per-round success probability P_succ.
 */
import { describe, expect, it } from "vitest";
import ref from "../verification/reference_4q.json";
import { allocTrajectory, Step } from "../src/type3map";
import { computeTrajectory4Q, q4Step } from "../src/q4map";

type Case = {
  eps: number; q: number; N: number;
  R1_dense: number[][]; R2_pennylane: number[][];
  R1_psucc: number[]; R2_psucc: number[];
  max_abs_R1_minus_R2: number;
  max_bell_offdiag_R1: number; max_bell_offdiag_R2: number;
  max_trace_err_R2: number; min_eig_R1: number; min_eig_R2: number;
};

const meta = (ref as { meta: { n_cnot: number; circuit_file_sha256: string } }).meta;
const cases = (ref as { cases: Case[] }).cases;
const TOL = 1e-12;

describe("Test G4Q — 4Q simulator vs. independent full-circuit references", () => {
  it("references come from the frozen 5-CNOT circuit", () => {
    expect(cases.length).toBeGreaterThanOrEqual(15);
    expect(meta.n_cnot).toBe(5);
    // Bell-state-distillation src/pqec_distill/circuit.py at 13c138b
    expect(meta.circuit_file_sha256.startsWith("e475a54da16a3d1c")).toBe(true);
  });

  for (const c of cases) {
    it(`eps=${c.eps}, q=${c.q}: Bell populations and P_succ for n = 0..${c.N}`, () => {
      // the two circuit routes agree with each other to their own conditioning; never looser than 1e-12
      const tol = Math.max(TOL, 4 * c.max_abs_R1_minus_R2);
      const t = computeTrajectory4Q(allocTrajectory(c.N), c.eps, c.q, c.N);
      let w1 = 0, w2 = 0;
      for (let n = 0; n <= c.N; n++) {
        const js = [t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]];
        for (let k = 0; k < 4; k++) {
          w1 = Math.max(w1, Math.abs(js[k] - c.R1_dense[n][k]));
          w2 = Math.max(w2, Math.abs(js[k] - c.R2_pennylane[n][k]));
        }
        // the references are physical states: unit trace, Bell-diagonal, Phi- = Psi- (4Q plane)
        expect(Math.abs(c.R1_dense[n].reduce((a, x) => a + x, 0) - 1)).toBeLessThan(1e-13);
        expect(Math.abs(c.R1_dense[n][1] - c.R1_dense[n][3])).toBeLessThan(tol);
        expect(Math.abs(c.R2_pennylane[n][1] - c.R2_pennylane[n][3])).toBeLessThan(tol);
      }
      expect(w1).toBeLessThan(tol);
      expect(w2).toBeLessThan(tol);
      expect(c.max_bell_offdiag_R1).toBeLessThan(1e-9);
      expect(c.max_bell_offdiag_R2).toBeLessThan(1e-9);
      expect(c.max_trace_err_R2).toBeLessThan(1e-13);
      expect(Math.min(c.min_eig_R1, c.min_eig_R2)).toBeGreaterThan(-1e-13);

      // normalisation: the physical success probability equals D_n / 4 every round, >= 1/4
      for (let n = 0; n < c.N; n++) {
        const P = q4Step(t.u[n], t.v[n], c.q)[2] / 4;
        expect(Math.abs(c.R1_psucc[n] - P)).toBeLessThan(1e-13);
        expect(Math.abs(c.R2_psucc[n] - P)).toBeLessThan(1e-13);
        expect(P).toBeGreaterThanOrEqual(0.25);
      }
    });
  }

  it("q = 1 is well defined for 4Q: one round sends the input to I/4 with P_succ = 1/4", () => {
    const c = cases.find((x) => x.q === 1)!;
    for (let n = 1; n <= c.N; n++) for (const p of c.R1_dense[n]) expect(Math.abs(p - 0.25)).toBeLessThan(1e-15);
    for (const P of c.R1_psucc) expect(Math.abs(P - 0.25)).toBeLessThan(1e-15);
  });

  it("eps=0.1, q=0.0615: the references and the simulator both decrease from round 0", () => {
    const c = cases.find((x) => x.eps === 0.1 && x.q === 0.0615)!;
    const F = c.R1_dense.map((p) => p[0]);
    expect(F[1]).toBeLessThan(F[0]);
    const t = computeTrajectory4Q(allocTrajectory(c.N), c.eps, c.q, c.N);
    expect(t.step[0]).toBe(Step.Decrease);
  });
});
