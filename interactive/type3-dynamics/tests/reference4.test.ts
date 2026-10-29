/**
 * Test G4 — Type-4 mode against independent reference trajectories.
 *
 * verification/reference_type4.json is produced by verification/generate_reference_type4.py
 * from two full-state routes that never use the (u, v) closed form:
 *   R1  dense 32x32 propagation through the noisy 14-CNOT circuit (iterated_noisy_pqec.py, "step4")
 *   R2  PennyLane default.mixed with the repository's Step-4 executor pqec_resynth_noise._apply;
 *       next state rebuilt from all 16 parity-weighted Pauli correlators
 * plus R3, the repository's one-round read-out pqec_resynth_noise.F_resynth.
 */
import { describe, expect, it } from "vitest";
import ref from "../verification/reference_type4.json";
import { allocTrajectory } from "../src/type3map";
import { computeTrajectory4 } from "../src/type4map";

type Case = {
  eps: number; q: number; N: number;
  R1_dense: number[][]; R2_pennylane: number[][];
  R2_parity_visibility: number[];
  R3_F1_resynth: number;
  max_abs_R1_minus_R2: number;
  max_bell_offdiag_R1: number; max_bell_offdiag_R2: number;
  max_trace_err_R2: number; min_eig_R2: number;
};

const meta = (ref as { meta: { n_cnot: number; circuit_file_sha256: string } }).meta;
const cases = (ref as { cases: Case[] }).cases;
const TOL = 1e-12;

describe("Test G4 — Type-4 simulator vs. independent full-circuit references", () => {
  it("references come from the frozen 14-CNOT circuit", () => {
    expect(cases.length).toBeGreaterThanOrEqual(15);
    expect(meta.n_cnot).toBe(14);
    // pqec_resynth_noise.py, also recorded by Bell-state-distillation/scripts/verify_type4_*_circuit.py
    expect(meta.circuit_file_sha256.startsWith("8d66cea64f494443")).toBe(true);
  });

  for (const c of cases) {
    it(`eps=${c.eps}, q=${c.q}: Bell populations, normalisation and F_1 for n = 0..${c.N}`, () => {
      // the two circuit routes agree with each other to their own conditioning; never looser than 1e-12
      const tol = Math.max(TOL, 4 * c.max_abs_R1_minus_R2);
      const t = computeTrajectory4(allocTrajectory(c.N), c.eps, c.q, c.N);
      let w1 = 0, w2 = 0;
      for (let n = 0; n <= c.N; n++) {
        const js = [t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]];
        for (let k = 0; k < 4; k++) {
          w1 = Math.max(w1, Math.abs(js[k] - c.R1_dense[n][k]));
          w2 = Math.max(w2, Math.abs(js[k] - c.R2_pennylane[n][k]));
        }
        // the references are physical states: unit trace, Bell-diagonal, Psi+ = Psi-
        expect(Math.abs(c.R1_dense[n].reduce((a, x) => a + x, 0) - 1)).toBeLessThan(1e-13);
        expect(Math.abs(c.R1_dense[n][2] - c.R1_dense[n][3])).toBeLessThan(tol);
      }
      expect(w1).toBeLessThan(tol);
      expect(w2).toBeLessThan(tol);
      expect(Math.abs(t.pPhiP[1] - c.R3_F1_resynth)).toBeLessThan(tol);
      expect(c.max_bell_offdiag_R1).toBeLessThan(1e-9);
      expect(c.max_bell_offdiag_R2).toBeLessThan(1e-9);
      expect(c.max_trace_err_R2).toBeLessThan(1e-13);
      expect(c.min_eig_R2).toBeGreaterThan(-1e-13);

      // normalisation: the circuit's parity visibility <Z_a> equals D_n r^10 / 4 every round
      const r = 1 - c.q;
      for (let n = 0; n < c.N; n++) {
        const u = t.u[n], v = t.v[n];
        const D = 1 + r * r * (2 * u * u + v * v);
        const z = c.R2_parity_visibility[n];
        expect(Math.abs(z - (D * r ** 10) / 4)).toBeLessThan(1e-12 * Math.max(1, z));
      }
    });
  }

  it("eps=0.1, q=0.0405: the references themselves show the rise and the later decrease", () => {
    const c = cases.find((x) => x.eps === 0.1 && x.q === 0.0405)!;
    const F = c.R1_dense.map((p) => p[0]);
    expect(F[1]).toBeGreaterThan(F[0]);
    const firstDec = F.findIndex((f, n) => n + 1 < F.length && F[n + 1] < f - 1e-12);
    expect(firstDec).toBe(2);
    const js = computeTrajectory4(allocTrajectory(c.N), c.eps, c.q, c.N);
    let jsFirst = -1;
    for (let n = 0; n < c.N; n++) if (js.step[n] === -1) { jsFirst = n; break; }
    expect(jsFirst).toBe(firstDec);
  });
});
