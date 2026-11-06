/**
 * What the theorems say about (protocol, q, eps), as data for the "Asymptotic attractor" panel.
 * Pure: everything shown comes from the exact regime, from certified brackets and enclosures, or
 * is labelled as an uncertified float64 orbit.
 *
 *   q < q_bi            P_B for every 0 < eps < 1
 *   q_bi < q < q_SN     P_B below eps_c(q), Z(q) above it
 *   q > q_SN            Z(q) for every 0 < eps < 1
 *
 * Outside these (q or eps at the ends of their range) no theorem is applied.
 */
import { basinVerdict, bistableFixedPoints, type BasinVerdict, type EpscBracket, initialFidelity, isDetermined } from "./epsc";
import type { EpscStatus } from "./epscClient";
import { regimeReport, type ExactRegime, type OrbitLimit, type RegimeProtocol, type RegimeReport, zPoint } from "./regimes";

export type Attractor = BasinVerdict | "out-of-range" | "pending";

export interface AttractorRow { k: string; v: string; cls?: string; txt?: boolean }

export interface AttractorReport {
  regime: ExactRegime;
  report: RegimeReport;
  /** The verdict for this eps. */
  attractor: Attractor;
  /** Certified bracket, for the eps slider. */
  bracket: EpscBracket | null;
  rows: AttractorRow[];
}

const num = (s: string) => `<span class="num">${s}</span>`;

/** Fixed-point text rounded away from the value, so a printed interval still contains it. */
export function fmtOutward(x: number, d: number, dir: "down" | "up"): string {
  const s = x.toFixed(d);
  const v = Number(s);
  const unit = 10 ** -d;
  if (dir === "down") return v < x ? s : (v - unit).toFixed(d);
  return v > x ? s : (v + unit).toFixed(d);
}

function intervalText(lo: number, hi: number): string {
  const w = hi - lo;
  const d = Math.min(15, Math.max(4, Math.ceil(-Math.log10(w)) + 1));
  return `[${fmtOutward(lo, d, "down")}, ${fmtOutward(hi, d, "up")}]`;
}

/** Text of Z(q) in the plane coordinates of the protocol. */
export function zText(protocol: RegimeProtocol, coords: readonly [string, string], q: number): string {
  const Z = zPoint(protocol, q);
  if (Z.kind === "mixed") return `혼합점 o (I/4), F = ${num(Z.F.toFixed(6))}`;
  // the invariant line is v (Type 3/4) or a (4Q); the transverse coordinate vanishes on it
  const [line, trans] = protocol === "q4" ? [coords[0], coords[1]] : [coords[1], coords[0]];
  return `B* (${line}* = ${num(Z.line.toFixed(6))}, ${trans} = 0), F = ${num(Z.F.toFixed(6))}`;
}

export function attractorReport(
  protocol: RegimeProtocol, coords: readonly [string, string], q: number, eps: number,
  status: EpscStatus, orbit: OrbitLimit,
  /** The 500-round chart orbit has a transverse component that is exactly 0. */
  chartTrapped = false,
): AttractorReport {
  const report = regimeReport(protocol, q);
  const { regime } = report;
  const rows: AttractorRow[] = [];
  const add = (k: string, v: string, cls?: string, txt = true) => rows.push({ k, v, cls, txt });
  const F0 = initialFidelity(eps);
  const epsOk = eps > 0 && eps < 1;
  const Z = zPoint(protocol, q);
  const fp = regime === "bistable" ? bistableFixedPoints(protocol, q) : null;
  const FB = fp ? (fp.PB.F[0] + fp.PB.F[1]) / 2 : NaN;
  const FS = fp ? (fp.S.F[0] + fp.S.F[1]) / 2 : NaN;
  let bracket: EpscBracket | null = null;
  let attractor: Attractor = "out-of-range";

  const labels: Record<ExactRegime, string> = {
    "unique-purified": "Region A: q &lt; q<sub>bi</sub>",
    bistable: "Bistable regime: q<sub>bi</sub> &lt; q &lt; q<sub>SN</sub>",
    "unique-unpurified": "Region C: q &gt; q<sub>SN</sub>",
    "on-threshold": "q = threshold",
  };
  if (!report.inRange) {
    add("Regime", `q = ${q}: 정리 범위(0 &lt; q &lt; 1) 밖, 판정 없음`, "st-unres");
    return { regime, report, attractor, bracket, rows };
  }
  add("Regime", `<b>${labels[regime]}</b>`);
  if (report.near) {
    const name = report.near === "qBi" ? "q<sub>bi</sub>" : "q<sub>SN</sub>";
    add("경계", `<span class="st-unres">q ≈ ${name} (차이 &lt; 1e-9). 위 영역은 exact polynomial 부호로 판정한 쪽이고, 등호에는 수렴 정리를 적용하지 않습니다.</span>`);
  }
  if (regime === "on-threshold") {
    add("판정", `<span class="st-unres">threshold 위: 수렴 정리 적용 안 함</span>`);
    return { regime, report, attractor, bracket, rows };
  }

  // While a bracket is computed the same rows stay, with empty values, so the panel keeps its height.
  const EMPTY = "—";
  const TWO = "two-line";
  if (regime === "bistable") {
    if (!status.bracket) {
      add("ε<sub>c</sub>(q)", status.state === "failed"
        ? `<span class="st-unres">계산 실패: 인증 구간 없음</span>`
        : `<span class="st-unres">계산 중…</span>`);
      add("Certified interval", EMPTY, TWO);
      add("Width", EMPTY);
      add("F<sub>c</sub> = 1 − 3ε<sub>c</sub>/4", EMPTY, TWO);
    } else {
      bracket = status.bracket;
      const { lo, hi } = bracket;
      const w = hi - lo;
      const running = status.state === "running";
      const midD = Math.min(14, Math.max(3, Math.floor(-Math.log10(w))));
      add("ε<sub>c</sub>(q)", `≈ ${num(((lo + hi) / 2).toFixed(midD))}`);
      add("Certified interval", num(intervalText(lo, hi)), TWO);
      add("Width", `${num(w.toExponential(1))}` + (running
        ? ` <span class="st-unres">· 계산 중</span>`
        : isDetermined(bracket) ? ` · N = ${num(`${bracket.roundsLo}, ${bracket.roundsHi}`)}`
        : ` <span class="st-unres">· 계산 미완료: ${bracket.stop === "budget" ? "라운드 한도 도달" : "float64 ball 해상도 한계"}</span>`));
      add("F<sub>c</sub> = 1 − 3ε<sub>c</sub>/4", num(intervalText(initialFidelity(hi), initialFidelity(lo))), TWO);
    }
  }

  add("Current ε", `${num(String(eps))}, F<sub>0</sub> = ${num(F0.toFixed(6))}`, "sep");

  if (!epsOk) {
    add("Predicted attractor", `<span class="st-unres">ε = ${eps}: 정리 범위(0 &lt; ε &lt; 1) 밖, 판정 없음</span>`);
    return { regime, report, attractor, bracket, rows };
  }

  const pbName = "Bell-connected fixed point P<sub>B</sub>";
  const zName = "separable boundary attractor Z(q)";
  if (regime === "unique-purified") {
    attractor = "P_B";
    add("Global attractor", `${pbName} (모든 0 &lt; ε &lt; 1)`);
  } else if (regime === "unique-unpurified") {
    attractor = "Z";
    add("Global attractor", `${zName} (모든 0 &lt; ε &lt; 1)`);
  } else if (!bracket) {
    attractor = "pending";
    add("Comparison", EMPTY);
    add("Predicted attractor", `<span class="st-unres">아직 판정 불가</span>`);
  } else {
    attractor = basinVerdict(bracket, eps);
    const cmp = { P_B: "ε ≤ ε<sub>L</sub> &lt; ε<sub>c</sub>", Z: "ε ≥ ε<sub>U</sub> &gt; ε<sub>c</sub>", undetermined: "ε<sub>L</sub> &lt; ε &lt; ε<sub>U</sub>" }[attractor];
    add("Comparison", cmp);
    add("Predicted attractor", attractor === "undetermined"
      ? `<span class="st-unres">basin 미확정: ε가 인증 구간 안에 있음${status.state === "running" ? " (계산 중)" : ""}</span>`
      : attractor === "P_B" ? pbName : zName);
  }

  // Values at the attractor: closed forms where the notes give them, else the labelled orbit.
  if (attractor === "Z") {
    add("Z(q)", zText(protocol, coords, q), TWO);
    add("F(Z) − F<sub>0</sub>", num((Z.F - F0).toExponential(3)));
  } else if (attractor === "P_B" && fp) {
    const [a, b] = protocol === "q4" ? [fp.PB.v, fp.PB.u] : [fp.PB.u, fp.PB.v];
    const mid = (x: readonly [number, number]) => ((x[0] + x[1]) / 2).toFixed(6);
    add("P<sub>B</sub>(q)", `(${coords[0]}, ${coords[1]}) = (${num(mid(a))}, ${num(mid(b))}), F = ${num(FB.toFixed(6))}`, TWO);
    add("F(P<sub>B</sub>) − F<sub>0</sub>", num((FB - F0).toExponential(3)));
  } else if (regime === "bistable") {
    add("Attractor", EMPTY, TWO);
    add("F(attractor) − F<sub>0</sub>", EMPTY);
  }

  // The uncertified orbit, only as a comparison, never as the verdict.
  const orbitF = num(orbit.F.toFixed(6));
  let check: string;
  if (!orbit.converged) {
    check = `<span class="st-unres">미수렴, 현재 F = ${orbitF}</span>`;
  } else if (regime === "bistable") {
    const near = (x: number) => Math.abs(orbit.F - x) < 1e-6;
    const found: Attractor | "S" | "?" = near(FB) ? "P_B" : near(Z.F) ? "Z" : near(FS) ? "S" : "?";
    const agree = found === attractor;
    check = `F = ${orbitF}, ${orbit.rounds.toLocaleString("en-US")} 라운드`;
    if (found === "S") check += ` <span class="st-unres">· saddle S 근처에 머문 상태</span>`;
    else if (attractor === "undetermined" || attractor === "pending") check += "";
    else if (!agree) check += ` <span class="st-unres">· 인증 판정과 다름: float64 orbit 한계</span>`;
    else check += " · 판정과 일치";
  } else if (regime === "unique-unpurified") {
    check = `F = ${orbitF}, ${orbit.rounds.toLocaleString("en-US")} 라운드` +
      (Math.abs(orbit.F - Z.F) < 1e-6 ? " · 판정과 일치" : ` <span class="st-unres">· Z(q)와 다름: float64 orbit 한계</span>`);
  } else {
    check = `F = ${orbitF}, ${orbit.rounds.toLocaleString("en-US")} 라운드`;
    add("F(orbit) − F<sub>0</sub>", num((orbit.F - F0).toExponential(3)) + " (비인증)");
    // below q_bi the transverse component can leave the float64 range on the way
    if (orbit.belowFloat64 || chartTrapped) {
      check += ` <span class="st-unres">· transverse 성분 float64 underflow: 극한은 log 추적 값${chartTrapped ? ", 차트 궤도는 불변선에 갇힌 값" : ""}</span>`;
    }
  }
  add("float64 orbit (비인증)", check, TWO);
  return { regime, report, attractor, bracket, rows };
}
