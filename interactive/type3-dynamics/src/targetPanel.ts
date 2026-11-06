/**
 * Rows of the "Asymptotic Target-Fidelity Threshold" panel.  Pure.  The same rows appear in every
 * state, so the panel keeps its height while F_target, eps and q move.
 *
 * Two different questions are kept apart:
 *   asymptotic   F_inf >= F_target    decided by the theorems (q_target, attractor, eps_c);
 *   finite round F_n >= F_target      a float64 comparison at the inspected round, nothing more.
 * For q > q_target a high-fidelity input can still satisfy F_n >= F_target for some rounds.
 *
 * Asymptotic rules used (Nov 6 analytic note Thm 1.1/5.1; "The target-fidelity threshold" Thm 2):
 *   0 < eps < 2/3:  q <= q_target  =>  F_inf = F_B(q) >= F_target;  q > q_target  =>  F_inf < F_target.
 *   0 < eps < 1, q > q_target: no input reaches the target (greatest fixed point, Thm 2(c)).
 *   2/3 <= eps < 1, q <= q_target: only through the attractor of this input (P_B: yes, Z: no).
 * A float q is never claimed to equal q_target: inside the certified bracket the verdict is
 * "At/Near threshold" (q = q_target itself meets the target, since the condition is >=).
 */
import { fmtOutward, type Attractor, type AttractorRow } from "./attractor";
import { epsCoverage, f0Enclosure, fbEnclosure, Q_ENT, qGainBracket, sideOfQEnt } from "./qgain";
import { qTargetBracket, targetEnclosure, targetVerdict } from "./qtarget";
import { regimeReport, zPoint, type RegimeProtocol } from "./regimes";
import { EPS64 } from "./type3map";

export type TargetAnswer = "yes" | "no" | "near" | "unresolved" | "not-guaranteed" | "not-covered";
/** F_target against F_0, from certified enclosures; "equal" when they cannot be separated. */
export type GainRelation = "above" | "below" | "equal" | "none";

export interface TargetReport {
  answer: TargetAnswer;
  lo: number | null;
  hi: number | null;
  relation: GainRelation;
  rows: AttractorRow[];
}

export interface TargetInput {
  protocol: RegimeProtocol;
  eps: number;
  q: number;
  /** Last valid F_target (1/2 < F < 1). */
  F: number;
  /** Verdict of the attractor panel for this (protocol, q, eps). */
  attractor: Attractor;
  /** Inspected round and its float64 fidelity. */
  n: number;
  Fn: number;
}

const num = (s: string) => `<span class="num">${s}</span>`;
const unres = (s: string) => `<span class="st-unres">${s}</span>`;
const EMPTY = "—";
const TWO = "two-line";
const POLY: Record<RegimeProtocol, string> = { type3: "P<sub>4</sub>", type4: "P<sub>2</sub>", q4: "P<sub>4Q</sub>" };

export function targetReport(inp: TargetInput): TargetReport {
  const { protocol, eps, q, F, attractor } = inp;
  const rows: AttractorRow[] = [];
  const add = (k: string, v: string, cls?: string) => rows.push({ k, v, cls, txt: true });
  const b = qTargetBracket(protocol, F);
  const cov = epsCoverage(eps);
  const epsOpen = eps > 0 && eps < 1;
  const reg = regimeReport(protocol, q);

  // 1-6. the threshold itself (independent of eps and q)
  add("Target fidelity F<sub>target</sub>", num(String(F)));
  if (b) {
    const w = b.hi - b.lo;
    const d = Math.min(15, Math.max(4, Math.ceil(-Math.log10(w)) + 1));
    add("q<sub>target</sub>(F<sub>target</sub>)",
      `≈ ${num(((b.lo + b.hi) / 2).toFixed(Math.min(14, Math.max(3, Math.floor(-Math.log10(w))))))} (width ${num(w.toExponential(1))})`, TWO);
    add("Certified interval", num(`[${fmtOutward(b.lo, d, "down")}, ${fmtOutward(b.hi, d, "up")}]`), TWO);
    add("Algebraic check", b.polyCheck
      ? `${POLY[protocol]}(s, 4F−1) 부호 변화 (exact) ✓`
      : unres(`${POLY[protocol]}(s, 4F−1) 부호 변화 확인 안 됨`), TWO);
  } else {
    add("q<sub>target</sub>(F<sub>target</sub>)", unres("계산 실패: 인증 구간 없음"), TWO);
    add("Certified interval", EMPTY, TWO);
    add("Algebraic check", EMPTY, TWO);
  }
  add("q<sub>ent</sub> (F<sub>B</sub> = 1/2)", num(Q_ENT[protocol].q.toFixed(8)));

  // relation to q_gain(eps) = F_B^{-1}(F_0): a consistency check of the two brackets
  let relation: GainRelation = "none";
  const g = qGainBracket(protocol, eps);
  let relText = EMPTY;
  if (g && b) {
    const T = targetEnclosure(F), F0 = f0Enclosure(eps);
    relation = T[0] > F0[1] ? "above" : T[1] < F0[0] ? "below" : "equal";
    // the brackets must not contradict the order (they may overlap when F_target and F_0 are close)
    const ok = relation === "above" ? b.lo < g.hi : relation === "below" ? b.hi > g.lo : b.lo < g.hi && g.lo < b.hi;
    const claim = {
      above: "F<sub>target</sub> &gt; F<sub>0</sub> ⇒ q<sub>target</sub> &lt; q<sub>gain</sub>",
      below: "F<sub>target</sub> &lt; F<sub>0</sub> ⇒ q<sub>target</sub> &gt; q<sub>gain</sub>",
      equal: "F<sub>target</sub> = F<sub>0</sub> (float64) ⇒ q<sub>target</sub> = q<sub>gain</sub>, 같은 threshold",
    }[relation];
    relText = `${claim} ${ok ? "✓" : unres("· 구간 불일치")}`;
  } else if (b) {
    relText = "q<sub>gain</sub>은 0 &lt; ε &lt; 2/3에서만 정의";
  }
  add("vs q<sub>gain</sub>(ε)", relText, TWO);

  // 7-8. current q and what the limit is
  add("Current CNOT noise q", num(String(q)), "sep");
  const pbByTheorem = cov === "covered" && reg.inRange && sideOfQEnt(protocol, q) < 0; // Thm 6.2 / 1.1
  const pb = reg.inRange && (pbByTheorem || attractor === "P_B");
  const FB = pb ? fbEnclosure(protocol, q) : null;
  let finf = EMPTY;
  if (FB) finf = `F<sub>B</sub>(q) = ${num(((FB[0] + FB[1]) / 2).toFixed(6))} (P<sub>B</sub>)`;
  else if (reg.inRange && attractor === "Z") finf = `F(Z) = ${num(zPoint(protocol, q).F.toFixed(6))} ≤ 1/2 (Z(q))`;
  else if (reg.inRange && epsOpen && (attractor === "undetermined" || attractor === "pending")) finf = unres("basin 미확정");
  add("Predicted F<sub>∞</sub>", finf, TWO);

  // 9-10. asymptotic verdict and the theorem used
  let answer: TargetAnswer;
  let line: string;
  let scope: string;
  if (cov === "covered") scope = "Thm 1.1 / 5.1: 0 &lt; ε &lt; 2/3, 1/2 &lt; F<sub>target</sub> &lt; 1 · q<sub>target</sub>는 ε와 무관";
  else if (epsOpen) scope = unres("Separable (ε ≥ 2/3): q &gt; q<sub>target</sub>이면 No, 아래는 attractor로 판정");
  else scope = unres(`ε = ${eps}: 정리 범위(0 &lt; ε &lt; 1) 밖`);

  const tv = b ? targetVerdict(b, q) : null;
  const yes = `<b class="st-inc">Yes</b> · q &lt; q<sub>target</sub>: F<sub>∞</sub> = F<sub>B</sub>(q) &gt; F<sub>target</sub>`;
  if (!reg.inRange) {
    answer = "unresolved";
    line = unres(`q = ${q}: 정리 범위(0 &lt; q &lt; 1) 밖, 판정 없음`);
  } else if (!b || !tv) {
    answer = "unresolved";
    line = unres("Unresolved: q<sub>target</sub> 인증 구간 없음");
  } else if (!epsOpen) {
    answer = "not-covered";
    line = unres("Not covered: 이 입력에는 정리를 적용하지 않음");
  } else if (tv === "miss") {
    answer = "no";
    line = `<b class="st-dec">No</b> · q &gt; q<sub>target</sub>: F<sub>∞</sub> &lt; F<sub>target</sub>`;
  } else if (tv === "near") {
    answer = "near";
    line = unres(relation === "equal"
      ? "At/Near threshold: q가 인증 구간 안. q = q<sub>target</sub> = q<sub>gain</sub>이면 target(≥)은 달성, strict gain(&gt;)은 없음"
      : "At/Near threshold: q가 인증 구간 안 (q = q<sub>target</sub>이면 F<sub>∞</sub> = F<sub>target</sub>, 달성)");
  } else if (cov === "covered") {
    answer = "yes";
    line = yes;
  } else if (attractor === "P_B") {
    answer = "yes";
    line = yes;
  } else if (attractor === "Z") {
    answer = "no";
    line = `<b class="st-dec">No</b> · 이 입력은 Z(q)로 수렴: F<sub>∞</sub> ≤ 1/2 &lt; F<sub>target</sub>`;
  } else {
    answer = "not-guaranteed";
    line = unres("Not guaranteed for this input: basin 미확정 (ε<sub>c</sub> 구간 안 또는 계산 중)");
  }
  add("Asymptotically: F<sub>∞</sub> ≥ F<sub>target</sub>", line, `sep ${TWO}`);
  add("Theorem applicability", scope, TWO);

  // 11. the finite-round comparison, kept apart
  const T = targetEnclosure(F);
  const tol = 8 * EPS64 * (Math.abs(inp.Fn) + 1);
  const fin = inp.Fn - T[1] > tol
    ? `<b class="st-inc">Yes</b> · F<sub>${inp.n}</sub> = ${num(inp.Fn.toFixed(6))}`
    : T[0] - inp.Fn > tol
      ? `<b class="st-dec">No</b> · F<sub>${inp.n}</sub> = ${num(inp.Fn.toFixed(6))}`
      : unres(`F<sub>${inp.n}</sub> ≈ F<sub>target</sub> (float64 해상도 안)`);
  add(`At current round n = ${inp.n}: F<sub>n</sub> ≥ F<sub>target</sub>`, `${fin} <span class="hint-inline">(float64, 점근 판정과 별개)</span>`, `sep ${TWO}`);

  return { answer, lo: b?.lo ?? null, hi: b?.hi ?? null, relation, rows };
}
