/**
 * Rows of the "Asymptotic Fidelity Gain" panel (Nov 6 note, Theorem 7.2).  Pure.  The same rows
 * appear in every state, so the panel keeps its height while eps and q move.
 *
 * Fidelity gain is about F_inf versus F_0 only.  "No" does not mean convergence to Z: for
 * q_gain < q <= q_ent the orbit still converges to P_B (Theorem 6.2), only with F_B(q) < F_0.
 */
import { fmtOutward, type AttractorRow } from "./attractor";
import { initialFidelity } from "./epsc";
import { epsCoverage, fbEnclosure, gainVerdict, Q_ENT, qGainBracket, sideOfQEnt, type GainVerdict } from "./qgain";
import { regimeReport, type RegimeProtocol } from "./regimes";

export interface GainReport {
  verdict: GainVerdict | "not-covered" | "out-of-range";
  /** Certified bracket of q_gain, for the q slider tick. */
  lo: number | null;
  hi: number | null;
  rows: AttractorRow[];
}

const num = (s: string) => `<span class="num">${s}</span>`;
const EMPTY = "—";
const TWO = "two-line";

export function gainReport(protocol: RegimeProtocol, eps: number, q: number): GainReport {
  const rows: AttractorRow[] = [];
  const add = (k: string, v: string, cls?: string) => rows.push({ k, v, cls, txt: true });
  const cov = epsCoverage(eps);
  const F0 = initialFidelity(eps);
  const b = qGainBracket(protocol, eps);
  const near23 = Math.abs(eps - 2 / 3) < 1e-9;

  // 1. range of the theorem
  if (cov === "covered") {
    add("Theorem range", near23
      ? `0 &lt; ε &lt; 2/3 <span class="st-unres">(ε는 표시 자릿수로는 2/3과 같지만 정확히는 2/3보다 작음)</span>`
      : "0 &lt; ε &lt; 2/3 (Theorem 7.2)", TWO);
  } else if (cov === "zero") {
    add("Theorem range", `<span class="st-unres">ε = ${eps}: 열린 구간 밖. 극한 q<sub>gain</sub> → 0만 알려짐</span>`, TWO);
  } else {
    add("Theorem range", `<span class="st-unres">Not covered by the current analytic result (ε ≥ 2/3${near23 ? ", 표시 자릿수로는 2/3과 같지만 정확히는 2/3보다 큼" : ""})</span>`, TWO);
  }

  // 2-3. q_gain and its certified interval
  if (b) {
    const w = b.hi - b.lo;
    const d = Math.min(15, Math.max(4, Math.ceil(-Math.log10(w)) + 1));
    add("q<sub>gain</sub>(ε) = q<sub>F</sub>(ε)", `≈ ${num(((b.lo + b.hi) / 2).toFixed(Math.min(14, Math.max(3, Math.floor(-Math.log10(w))))))} (width ${num(w.toExponential(1))})`, TWO);
    add("Certified interval", num(`[${fmtOutward(b.lo, d, "down")}, ${fmtOutward(b.hi, d, "up")}]`), TWO);
  } else {
    add("q<sub>gain</sub>(ε) = q<sub>F</sub>(ε)", cov === "covered" ? `<span class="st-unres">계산 실패</span>` : EMPTY, TWO);
    add("Certified interval", EMPTY, TWO);
  }
  add("q<sub>ent</sub> (F<sub>B</sub> = 1/2)", num(Q_ENT[protocol].q.toFixed(8)));

  // 4-6. inputs and F_B(q)
  add("Initial Fidelity", num(F0.toFixed(6)), "sep");
  add("Current gate noise q", num(String(q)));
  const reg = regimeReport(protocol, q);
  const hasPB = reg.inRange && (reg.regime === "unique-purified" || reg.regime === "bistable");
  const FB = hasPB ? fbEnclosure(protocol, q) : null;
  add("F<sub>B</sub>(q) at P<sub>B</sub>", FB
    ? num(((FB[0] + FB[1]) / 2).toFixed(6))
    : reg.inRange && reg.regime === "unique-unpurified" ? "P<sub>B</sub> 없음 (q &gt; q<sub>SN</sub>)" : EMPTY);

  // 7-8. verdict and what it does not mean
  let verdict: GainReport["verdict"];
  let line: string, note: string;
  if (cov !== "covered") {
    verdict = "not-covered";
    line = `<span class="st-unres">${cov === "zero" ? "판정 없음" : "Not covered by the current analytic result"}</span>`;
    note = EMPTY;
  } else if (!reg.inRange) {
    verdict = "out-of-range";
    line = `<span class="st-unres">q = ${q}: 정리 범위(0 &lt; q &lt; 1) 밖, 판정 없음</span>`;
    note = EMPTY;
  } else if (!b) {
    verdict = "uncertain";
    line = `<span class="st-unres">판정 불가</span>`;
    note = EMPTY;
  } else {
    verdict = gainVerdict(b, q);
    const belowEnt = sideOfQEnt(protocol, q) < 0;
    if (verdict === "gain") {
      line = `<b class="st-inc">Yes</b> · q &lt; q<sub>gain</sub>: F<sub>∞</sub> = F<sub>B</sub>(q) &gt; F<sub>0</sub>`;
      note = "P<sub>B</sub>로 수렴 (Thm 6.2)";
    } else if (verdict === "no-gain") {
      line = `<b class="st-dec">No</b> · q &gt; q<sub>gain</sub>: F<sub>∞</sub> &lt; F<sub>0</sub>`;
      note = belowEnt
        ? "P<sub>B</sub>로 수렴하지만 F<sub>B</sub>(q) &lt; F<sub>0</sub> (Z 아님, Thm 6.2)"
        : reg.regime === "bistable"
          ? "P<sub>B</sub> 또는 Z (ε<sub>c</sub>는 위 패널), 어느 쪽이든 F<sub>∞</sub> &lt; 1/2 &lt; F<sub>0</sub>"
          : reg.regime === "unique-unpurified"
            ? "Z(q)로 수렴, F ≤ 1/2 &lt; F<sub>0</sub>"
            : "F<sub>∞</sub> &lt; 1/2 &lt; F<sub>0</sub>";
    } else {
      line = `<span class="st-unres">판정 불확실: q가 q<sub>gain</sub>의 인증 구간 안 (q = q<sub>gain</sub>이면 F<sub>∞</sub> = F<sub>0</sub>)</span>`;
      note = EMPTY;
    }
  }
  add("Asymptotic fidelity gain", line, `sep ${TWO}`);
  add("Attractor", note, TWO);
  return { verdict, lo: b?.lo ?? null, hi: b?.hi ?? null, rows };
}
