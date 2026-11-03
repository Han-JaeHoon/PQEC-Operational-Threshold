/**
 * Test GB — BBPSSW mode against independent circuit-level references.
 *
 * verification/reference_bbpssw.json is produced by verification/generate_reference_bbpssw.py from
 * two routes through the noisy four-qubit circuit that never use the scalar recurrence:
 *   R1  verification/bbpssw_stage1/frozen_circuit_reference.py (16 Pauli conjugations, 24 Cliffords)
 *   R2  written in the generator (index-relabeled CNOT, literal partial-trace replacement channel)
 */
import { describe, expect, it } from "vitest";
import ref from "../verification/reference_bbpssw.json";
import { allocTrajectory } from "../src/type3map";
import {
  bbpsswPSucc, bbpsswRoundOriginal, bbpsswRoundUntwirled, computeTrajectoryBBPSSW, type BellPops,
} from "../src/bbpsswmap";

type Case = {
  eps: number; q: number; N: number;
  R1_circuit: number[][]; R2_circuit: number[][];
  R1_psucc: number[]; R2_psucc: number[];
  max_abs_R1_minus_R2: number;
};
type Round = [number[], number];
type Arbitrary = {
  p: number[]; q: number;
  original_R1: Round; original_R2: Round;
  no_input_twirl_R1: Round; no_input_twirl_R2: Round;
  no_twirls_R1: Round; no_twirls_R2: Round;
};

const data = ref as unknown as {
  meta: { n_cnot: number; frozen_circuit_reference_sha256: string };
  cases: Case[];
  arbitrary: Arbitrary[];
  negative_controls: Record<string, unknown> & { p: number[]; q: number };
};
const TOL = 1e-12;
const close = (a: readonly number[], b: readonly number[]) =>
  Math.max(...a.map((x, i) => Math.abs(x - b[i])));

describe("Test GB — BBPSSW simulator vs. independent circuit references", () => {
  it("references come from the frozen 2-CNOT circuit", () => {
    expect(data.meta.n_cnot).toBe(2);
    expect(data.meta.frozen_circuit_reference_sha256.startsWith("782d157d8bc3ccb3")).toBe(true);
    expect(data.cases.length).toBe(108);
  });

  it("H: repeated rounds: populations and P_succ on 108 (eps, q) grids", () => {
    let worst = 0, worstP = 0, worstRoutes = 0;
    for (const c of data.cases) {
      const t = computeTrajectoryBBPSSW(allocTrajectory(c.N), c.eps, c.q, c.N);
      worstRoutes = Math.max(worstRoutes, c.max_abs_R1_minus_R2);
      for (let n = 0; n <= c.N; n++) {
        const got = [t.pPhiP[n], t.pPhiM[n], t.pPsiP[n], t.pPsiM[n]];
        worst = Math.max(worst, close(got, c.R1_circuit[n]), close(got, c.R2_circuit[n]));
        if (n < c.N) {
          const P = bbpsswPSucc(t.pPhiP[n], c.q);
          worstP = Math.max(worstP, Math.abs(P - c.R1_psucc[n]), Math.abs(P - c.R2_psucc[n]));
        }
      }
    }
    console.log(`max error: populations ${worst}, P_succ ${worstP}; |R1-R2| ${worstRoutes}`);
    expect(worst).toBeLessThan(TOL);
    expect(worstP).toBeLessThan(TOL);
    expect(worstRoutes).toBeLessThan(TOL);
  });

  it("I: arbitrary inputs: the original map matches the circuit with both twirls", () => {
    for (const a of data.arbitrary) {
      const [p, P] = bbpsswRoundOriginal(a.p as unknown as BellPops, a.q);
      for (const [pops, succ] of [a.original_R1, a.original_R2]) {
        expect(close(p, pops)).toBeLessThan(TOL);
        expect(Math.abs(P - succ)).toBeLessThan(TOL);
      }
    }
  });

  it("I: arbitrary inputs: the no-twirl map matches the circuit without twirls", () => {
    for (const a of data.arbitrary) {
      const [p, P] = bbpsswRoundUntwirled(a.p as unknown as BellPops, a.q);
      for (const [pops, succ] of [a.no_twirls_R1, a.no_twirls_R2]) {
        expect(close(p, pops)).toBeLessThan(TOL);
        expect(Math.abs(P - succ)).toBeLessThan(TOL);
      }
    }
  });

  it("I: the circuit itself separates original from no-input-twirl, and the scalar map is not the latter", () => {
    let separated = 0;
    for (const a of data.arbitrary) {
      if (a.q === 1) continue;
      const original = bbpsswRoundOriginal(a.p as unknown as BellPops, a.q)[0][0];
      const noInput = a.no_input_twirl_R1[0][0];
      expect(Math.abs(noInput - a.no_input_twirl_R2[0][0])).toBeLessThan(TOL);
      if (Math.abs(original - noInput) > 1e-3) separated++;
    }
    expect(separated).toBeGreaterThan(100);
  });

  describe("J: negative controls are detected", () => {
    const { p, q } = data.negative_controls;
    const [expected] = bbpsswRoundUntwirled(p as unknown as BellPops, q);
    const h = (1 - expected[0]) / 3;
    const target = [expected[0], h, h, h];
    for (const k of ["reverse_Bob_CNOT", "omit_Bob_depolarization", "accept_unequal_measurements"]) {
      it(k, () => {
        const [pops] = data.negative_controls[k] as Round;
        const diff = close(target, pops);
        console.log(`${k}: population difference ${diff}`);
        expect(diff).toBeGreaterThan(1e-3);
      });
    }
  });
});
