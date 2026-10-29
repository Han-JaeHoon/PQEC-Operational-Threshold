import uPlot from "uplot";
import "uplot/dist/uPlot.min.css";
import "./style.css";
import {
  allocTrajectory, computeTrajectory, summarize, Step, type Trajectory,
} from "./type3map";
import { qAbs, qRel } from "./thresholds";

const MAX_N = 500;
/** |dF| within this factor of its rounding bound is flagged "near precision limit". */
const NEAR_LIMIT = 1e3;
/** Linear half-width of the signed-log dF axis. */
const ASINH_LIN = 1e-13;

// ---------------------------------------------------------------- state
const state = {
  eps: 0.1,
  q: 0.0327,
  N: 30,
  sel: 0,
  qFull: false,
  fAxis: "auto" as "auto" | "fixed",
  dAxis: "lin" as "lin" | "asinh",
};
const traj: Trajectory = allocTrajectory(MAX_N);
/** Always MAX_N rounds, independent of the displayed N: feeds the monotonicity verdict. */
const trajLong: Trajectory = allocTrajectory(MAX_N);
/** First resolved decrease F_{n+1} < F_n over n = 0..MAX_N-1, or -1. */
let firstDecLong = -1;

// ---------------------------------------------------------------- DOM
const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;
const epsRange = $<HTMLInputElement>("eps-range");
const epsNum = $<HTMLInputElement>("eps-num");
const qRange = $<HTMLInputElement>("q-range");
const qNum = $<HTMLInputElement>("q-num");
const qFull = $<HTMLInputElement>("q-full");
const nRange = $<HTMLInputElement>("n-range");
const nNum = $<HTMLInputElement>("n-num");
const selRange = $<HTMLInputElement>("sel-range");
const selNum = $<HTMLInputElement>("sel-num");
const inspector = $<HTMLTableElement>("inspector");
const summaryTbl = $<HTMLTableElement>("summary");
const precWarn = $<HTMLParagraphElement>("precision-warn");
const timing = $<HTMLSpanElement>("timing");
const qMarkAbs = $<HTMLSpanElement>("qmark-abs");
const qMarkRel = $<HTMLSpanElement>("qmark-rel");
const qThLine = $<HTMLParagraphElement>("q-thresholds");
const monoBadge = $<HTMLSpanElement>("mono");

// ---------------------------------------------------------------- formatting
const fmtF = (x: number) => x.toFixed(12);
const fmtE = (x: number) => (x === 0 ? "0" : x.toExponential(6));
const fmtQ = (x: number) => String(+x.toPrecision(10));

function stepLabel(n: number): [string, string] {
  const s = traj.step[n];
  if (s === Step.Increase) return ["increase", "st-inc"];
  if (s === Step.Decrease) return ["decrease", "st-dec"];
  return ["unresolved (≤ float64 error bound)", "st-unres"];
}

function isNearLimit(n: number): boolean {
  return traj.step[n] !== Step.Unresolved && Math.abs(traj.dF[n]) < NEAR_LIMIT * traj.dFerr[n];
}

// ---------------------------------------------------------------- charts
function css(name: string): string {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

const SYNC = uPlot.sync("rounds");

type AxisFmt = "int" | "fixed" | "exp";

/** Tick labels with enough digits that adjacent ticks stay distinguishable. */
function tickValues(fmt: AxisFmt) {
  return (_u: uPlot, splits: number[], _ai: number, _space: number, incr: number) =>
    splits.map((v) => {
      if (v == null) return "";
      if (fmt === "int") return Number.isInteger(v) ? String(v) : "";
      if (fmt === "exp") return v === 0 ? "0" : v.toExponential(Math.abs(v) >= 1e-3 && incr >= 1e-4 ? 1 : 2);
      const dec = Math.min(15, Math.max(1, Math.ceil(-Math.log10(incr) + 1e-9)));
      return v.toFixed(dec);
    });
}

function axis(label: string, size: number, fmt: AxisFmt): uPlot.Axis {
  return {
    label,
    values: tickValues(fmt),
    labelSize: 18,
    size,
    stroke: css("--axis"),
    grid: { stroke: css("--grid"), width: 1 },
    ticks: { stroke: css("--grid"), width: 1 },
    font: "11px ui-monospace, SFMono-Regular, Menlo, monospace",
    labelFont: "12px -apple-system, BlinkMacSystemFont, Helvetica, Arial, sans-serif",
  };
}

/** Vertical marker at the selected round, drawn without rebuilding any series path. */
function selMarker(u: uPlot) {
  const x = u.valToPos(state.sel, "x", true);
  if (!Number.isFinite(x) || x < u.bbox.left - 1 || x > u.bbox.left + u.bbox.width + 1) return;
  const ctx = u.ctx;
  ctx.save();
  ctx.strokeStyle = css("--c-sel");
  ctx.lineWidth = 1.5 * devicePixelRatio;
  ctx.setLineDash([4 * devicePixelRatio, 3 * devicePixelRatio]);
  ctx.beginPath();
  ctx.moveTo(x, u.bbox.top);
  ctx.lineTo(x, u.bbox.top + u.bbox.height);
  ctx.stroke();
  ctx.restore();
}

function zeroLine(u: uPlot) {
  const y = u.valToPos(0, "y", true);
  if (!Number.isFinite(y) || y < u.bbox.top || y > u.bbox.top + u.bbox.height) return;
  const ctx = u.ctx;
  ctx.save();
  ctx.strokeStyle = css("--axis");
  ctx.lineWidth = devicePixelRatio;
  ctx.beginPath();
  ctx.moveTo(u.bbox.left, y);
  ctx.lineTo(u.bbox.left + u.bbox.width, y);
  ctx.stroke();
  ctx.restore();
}

/** Dashed F_0 reference, drawn as an overlay so it never widens the auto-zoomed y range. */
function f0Line(u: uPlot) {
  const F0 = traj.pPhiP[0];
  const y = u.valToPos(F0, "y", true);
  if (!Number.isFinite(y) || y < u.bbox.top - 1 || y > u.bbox.top + u.bbox.height + 1) return;
  const ctx = u.ctx;
  const dpr = devicePixelRatio;
  ctx.save();
  ctx.strokeStyle = css("--c-ref");
  ctx.fillStyle = css("--c-ref");
  ctx.lineWidth = dpr;
  ctx.setLineDash([6 * dpr, 4 * dpr]);
  ctx.beginPath();
  ctx.moveTo(u.bbox.left, y);
  ctx.lineTo(u.bbox.left + u.bbox.width, y);
  ctx.stroke();
  ctx.font = `${11 * dpr}px -apple-system, BlinkMacSystemFont, Helvetica, Arial, sans-serif`;
  ctx.textAlign = "right";
  ctx.fillText(`F₀ = ${F0.toFixed(6)}`, u.bbox.left + u.bbox.width - 4 * dpr, y - 4 * dpr);
  ctx.restore();
}

/** Hollow ring on the first round whose successor is lower (F_{n+1} < F_n). */
function decMarker(u: uPlot) {
  const n = firstDecLong;
  if (n < 0 || n > traj.N) return;
  const x = u.valToPos(n, "x", true);
  const y = u.valToPos(traj.pPhiP[n], "y", true);
  if (!Number.isFinite(x) || !Number.isFinite(y)) return;
  if (x < u.bbox.left - 1 || x > u.bbox.left + u.bbox.width + 1) return;
  const ctx = u.ctx;
  const dpr = devicePixelRatio;
  ctx.save();
  ctx.strokeStyle = css("--c-dec");
  ctx.lineWidth = 1.5 * dpr;
  ctx.beginPath();
  ctx.arc(x, y, 6 * dpr, 0, 2 * Math.PI);
  ctx.stroke();
  ctx.restore();
}

/** Sparse decade splits for the signed-log (asinh) dF axis. */
function asinhSplits(_u: uPlot, _ai: number, min: number, max: number): number[] {
  const top = Math.max(Math.abs(min), Math.abs(max), ASINH_LIN);
  const kMax = Math.floor(Math.log10(top));
  const kMin = Math.round(Math.log10(ASINH_LIN));
  const step = Math.max(1, Math.ceil((kMax - kMin + 1) / 3));
  const out = [0];
  for (let k = kMax; k >= kMin; k -= step) {
    const v = 10 ** k;
    if (v <= max) out.push(v);
    if (-v >= min) out.push(-v);
  }
  return out.sort((a, b) => a - b);
}

function padRange(min: number, max: number, includeZero = false): uPlot.Range.MinMax {
  if (includeZero) { min = Math.min(min, 0); max = Math.max(max, 0); }
  const span = max - min;
  const pad = Math.max(span * 0.08, Math.max(Math.abs(min), Math.abs(max)) * 1e-13, 1e-300);
  return [min - pad, max + pad];
}

function baseOpts(el: HTMLElement): Pick<uPlot.Options, "width" | "height" | "cursor" | "legend"> {
  return {
    width: el.clientWidth,
    height: el.clientHeight - 30,
    cursor: {
      sync: { key: SYNC.key },
      drag: { x: true, y: false, setScale: false },
      points: { size: 7 },
    },
    legend: { live: true },
  };
}

// ---------------------------------------------------------------- shared round-window zoom
/** Zoomed round window shared by all charts; null = all rounds. y re-ranges to the window. */
let xWin: [number, number] | null = null;

function applyZoom() {
  const N = traj.N;
  if (xWin) {
    const lo = Math.max(0, xWin[0]);
    const hi = Math.min(N, xWin[1]);
    if (hi - lo < 1) xWin = null;
  }
  for (const u of [chartF, chartP, chartD]) {
    const hi = u.data[0].length - 1;
    const min = xWin ? xWin[0] : 0;
    const max = xWin ? Math.min(xWin[1], N) : Math.max(hi, 1);
    u.setScale("x", { min, max });
  }
}

function onSelect(u: uPlot) {
  const { left, width } = u.select;
  if (width > 4) {
    const a = u.posToVal(left, "x");
    const b = u.posToVal(left + width, "x");
    xWin = [Math.max(0, Math.floor(a)), Math.ceil(b)];
    applyZoom();
  }
  u.setSelect({ left: 0, top: 0, width: 0, height: 0 }, false);
}

const elF = $<HTMLDivElement>("chart-f");
const elP = $<HTMLDivElement>("chart-p");
const elD = $<HTMLDivElement>("chart-d");
let chartF: uPlot, chartP: uPlot, chartD: uPlot;

function roundValue(_u: uPlot, v: number | null) {
  return v == null ? "–" : String(v);
}

function buildCharts() {
  chartF?.destroy();
  chartP?.destroy();
  chartD?.destroy();

  chartF = new uPlot({
    ...baseOpts(elF),
    scales: {
      x: { time: false },
      y: {
        range: (_u, min, max) => (state.fAxis === "fixed" ? [0, 1] : padRange(min, max)),
      },
    },
    axes: [axis("round n", 40, "int"), axis("F_n", 110, "fixed")],
    series: [
      { label: "n", value: roundValue },
      {
        label: "F_n", stroke: css("--c-f"), width: 2,
        points: { show: (_u, _si, i0, i1) => i1 - i0 <= 120, size: 5, fill: css("--c-f") },
        value: (_u, v) => (v == null ? "–" : fmtF(v)),
      },
    ],
    hooks: { draw: [f0Line, decMarker, selMarker], setSelect: [onSelect] },
  }, [[], []], elF);

  const pop = (label: string, c: string, dash?: number[]): uPlot.Series => ({
    label, stroke: css(c), width: 1.75, dash,
    points: { show: false },
    value: (_u, v) => (v == null ? "–" : fmtF(v)),
  });
  chartP = new uPlot({
    ...baseOpts(elP),
    scales: { x: { time: false }, y: { range: (_u, min, max) => padRange(min, max) } },
    axes: [axis("round n", 40, "int"), axis("p_β", 110, "fixed")],
    series: [
      { label: "n", value: roundValue },
      pop("Φ⁺", "--c-phip"),
      pop("Φ⁻", "--c-phim"),
      pop("Ψ⁺", "--c-psip"),
      pop("Ψ⁻", "--c-psim", [5, 4]),
    ],
    hooks: { draw: [selMarker], setSelect: [onSelect] },
  }, [[], [], [], [], []], elP);

  const dot = (label: string, c: string): uPlot.Series => ({
    label, stroke: css(c), width: 0,
    paths: () => null,
    points: { show: true, size: 6, fill: css(c), stroke: css(c) },
    value: (_u, v) => (v == null ? "–" : fmtE(v)),
  });
  chartD = new uPlot({
    ...baseOpts(elD),
    scales: {
      x: { time: false },
      y: state.dAxis === "asinh"
        ? { distr: 4, asinh: ASINH_LIN }
        : { range: (_u, min, max) => padRange(min, max, true) },
    },
    axes: [
      axis("round n", 40, "int"),
      state.dAxis === "asinh" ? { ...axis("ΔF_n", 110, "exp"), splits: asinhSplits, filter: (_u: uPlot, s: number[]) => s } : axis("ΔF_n", 110, "exp"),
    ],
    series: [
      { label: "n", value: roundValue },
      {
        label: "ΔF_n", stroke: css("--c-ref"), width: 1, points: { show: false },
        value: (_u, v) => (v == null ? "–" : fmtE(v)),
      },
      dot("increase", "--c-inc"),
      dot("decrease", "--c-dec"),
      dot("unresolved", "--c-unres"),
    ],
    hooks: { draw: [zeroLine, selMarker], setSelect: [onSelect] },
  }, [[], [], [], [], []], elD);

  for (const u of [chartF, chartP, chartD]) {
    let downX = 0;
    u.over.addEventListener("mousedown", (e) => { downX = e.clientX; });
    u.over.addEventListener("click", (e) => {
      if (Math.abs(e.clientX - downX) > 4) return; // that was a zoom drag
      const i = u.cursor.idx;
      if (i != null) setSelected(i);
    });
    u.over.addEventListener("dblclick", () => { xWin = null; applyZoom(); });
  }
  pushData();
}

// ---------------------------------------------------------------- data push
const xs = Array.from({ length: MAX_N + 2 }, (_, i) => i);
const dInc: (number | null)[] = new Array(MAX_N + 1);
const dDec: (number | null)[] = new Array(MAX_N + 1);
const dUnr: (number | null)[] = new Array(MAX_N + 1);

function pushData() {
  const N = traj.N;
  const n1 = N + 1;
  const x = xs.slice(0, n1);
  chartF.setData([x, traj.pPhiP.subarray(0, n1)] as unknown as uPlot.AlignedData);
  chartP.setData([
    x, traj.pPhiP.subarray(0, n1), traj.pPhiM.subarray(0, n1),
    traj.pPsiP.subarray(0, n1), traj.pPsiM.subarray(0, n1),
  ] as unknown as uPlot.AlignedData);

  // dF_n for n = 0..N-1 (both F_n and F_{n+1} displayed)
  const nd = Math.max(N, 1);
  for (let i = 0; i < nd; i++) {
    const s = traj.step[i];
    const d = traj.dF[i];
    dInc[i] = s === Step.Increase ? d : null;
    dDec[i] = s === Step.Decrease ? d : null;
    dUnr[i] = s === Step.Unresolved ? d : null;
  }
  chartD.setData([
    xs.slice(0, nd), traj.dF.subarray(0, nd),
    dInc.slice(0, nd), dDec.slice(0, nd), dUnr.slice(0, nd),
  ] as unknown as uPlot.AlignedData);
  if (xWin) applyZoom();
}

// ---------------------------------------------------------------- panels
function row(k: string, v: string, cls = "", txt = false): string {
  return `<tr class="${cls}"><th>${k}</th><td class="${txt ? "txt" : ""}">${v}</td></tr>`;
}

function renderInspector() {
  const n = state.sel;
  const [lab, cls] = stepLabel(n);
  const near = isNearLimit(n);
  inspector.innerHTML =
    row("u<sub>n</sub> = ⟨XX⟩", fmtF(traj.u[n])) +
    row("v<sub>n</sub> = ⟨ZZ⟩", fmtF(traj.v[n])) +
    row("p<sub>Φ⁺</sub> = F<sub>n</sub>", fmtF(traj.pPhiP[n]), "sep") +
    row("p<sub>Φ⁻</sub>", fmtF(traj.pPhiM[n])) +
    row("p<sub>Ψ⁺</sub>", fmtF(traj.pPsiP[n])) +
    row("p<sub>Ψ⁻</sub>", fmtF(traj.pPsiM[n])) +
    row("ΔF<sub>n</sub> (F<sub>n+1</sub>−F<sub>n</sub>)", fmtE(traj.dF[n]), "sep") +
    row("ΔF<sub>n</sub> (drift numerator)", fmtE(traj.dFmap[n])) +
    row("error bounds", `${traj.dFerr[n].toExponential(1)} / ${traj.dFmapErr[n].toExponential(1)}`) +
    row("step n→n+1",
      `<span class="${cls}">${lab}</span>` +
      (near ? ` <span class="st-unres">· near precision limit</span>` : ""), "", true);
}

function renderSummary() {
  const s = summarize(traj);
  const N = traj.N;
  const dec = s.firstDecrease >= 0
    ? `<span class="st-dec">n = ${s.firstDecrease} (F<sub>${s.firstDecrease + 1}</sub> &lt; F<sub>${s.firstDecrease}</sub>)</span>`
    : `No decrease detected up to round ${N}`;
  const unres = s.unresolved === 0
    ? "none"
    : `<span class="st-unres">${s.unresolved} step(s), first at n = ${s.firstUnresolved}</span>`;
  const [a, b, c, d] = s.finalPops;
  summaryTbl.innerHTML =
    row("F<sub>0</sub>", fmtF(s.F0)) +
    row(`F<sub>N</sub> (N = ${N})`, fmtF(s.FN)) +
    row(`max F (n = ${s.maxFRound})`, fmtF(s.maxF)) +
    row(`min F (n = ${s.minFRound})`, fmtF(s.minF)) +
    row("first resolved decrease", dec, "sep", true) +
    row("unresolved sign", unres, "", true) +
    row(`p<sub>Φ⁺</sub><sup>(N)</sup>`, fmtF(a), "sep") +
    row(`p<sub>Φ⁻</sub><sup>(N)</sup>`, fmtF(b)) +
    row(`p<sub>Ψ⁺</sub><sup>(N)</sup>`, fmtF(c)) +
    row(`p<sub>Ψ⁻</sub><sup>(N)</sup>`, fmtF(d));

  if (s.unresolved > 0) {
    precWarn.hidden = false;
    precWarn.textContent =
      `${s.unresolved} of ${N} steps have |ΔF| at or below the float64 error bound ` +
      `(or the two evaluation routes disagree): no sign is claimed for them (amber).`;
  } else {
    let near = 0;
    for (let i = 0; i < N; i++) if (isNearLimit(i)) near++;
    precWarn.hidden = near === 0;
    precWarn.textContent = near
      ? `${near} step(s) resolved but within ${NEAR_LIMIT}× of the float64 error bound.`
      : "";
  }
}

/** Analytic q_abs / q_rel (calculation notes) as ticks on the q slider + one text line. */
function renderThresholds() {
  const a = qAbs(state.eps);
  const r = qRel(state.eps);
  const qmax = +qRange.max;
  for (const [el, v] of [[qMarkAbs, a], [qMarkRel, r]] as const) {
    const ok = Number.isFinite(v) && v <= qmax;
    el.hidden = !ok;
    if (ok) el.style.setProperty("--f", String(v / qmax));
  }
  if (!Number.isFinite(a)) {
    qThLine.innerHTML =
      `<b class="k-abs">q<sub>abs</sub></b>, <b class="k-rel">q<sub>rel</sub></b>: ` +
      `defined only for entangled input, 0 &lt; ε &lt; 2/3`;
    return;
  }
  const where = state.q < a ? "below both" : state.q < r ? "between q<sub>abs</sub> and q<sub>rel</sub>" : "above both";
  qThLine.innerHTML =
    `<b class="k-abs">q<sub>abs</sub>(ε)</b> = <span class="num">${a.toFixed(6)}</span> · ` +
    `<b class="k-rel">q<sub>rel</sub>(ε)</b> = <span class="num">${r.toFixed(6)}</span> · ` +
    `q is ${where}`;
}

/** Monotonicity verdict for F_n over a fixed horizon of MAX_N rounds. */
function renderMonotone() {
  computeTrajectory(trajLong, state.eps, state.q, MAX_N);
  const s = summarize(trajLong);
  firstDecLong = s.firstDecrease;
  if (s.firstDecrease >= 0) {
    const n = s.firstDecrease;
    monoBadge.className = "mono dec";
    monoBadge.textContent = `↓ decreases at n = ${n}` + (n >= traj.N ? " (beyond N)" : "");
    monoBadge.title = `First resolved decrease: F_${n + 1} < F_${n} (both evaluation routes agree).`;
  } else {
    monoBadge.className = "mono ok";
    monoBadge.textContent = `↑ no decrease, n ≤ ${MAX_N}`;
    monoBadge.title =
      `No resolved decrease F_{n+1} < F_n for n = 0…${MAX_N - 1}.` +
      (s.unresolved > 0
        ? ` ${s.unresolved} step(s) near the fixed point are below the float64 error bound; no sign is claimed for them.`
        : "") +
      ` Finite-round numerical check, not a proof.`;
  }
}

// ---------------------------------------------------------------- update loop
let pending = false;
let dataDirty = false;
const times: number[] = [];

function schedule(recompute: boolean) {
  dataDirty ||= recompute;
  if (pending) return;
  pending = true;
  requestAnimationFrame(frame);
}

function frame() {
  pending = false;
  const t0 = performance.now();
  if (dataDirty) {
    dataDirty = false;
    computeTrajectory(traj, state.eps, state.q, state.N);
    renderMonotone();
    const t1 = performance.now();
    renderThresholds();
    pushData();
    renderSummary();
    renderInspector();
    const t2 = performance.now();
    times.push(t2 - t0);
    if (times.length > 120) times.shift();
    const avg = times.reduce((a, b) => a + b, 0) / times.length;
    const mx = Math.max(...times);
    timing.textContent =
      `update: compute ${(t1 - t0).toFixed(2)} ms · total ${(t2 - t0).toFixed(2)} ms ` +
      `(avg ${avg.toFixed(2)}, max ${mx.toFixed(2)} over last ${times.length})`;
  } else {
    renderInspector();
    for (const u of [chartF, chartP, chartD]) u.redraw(false, false);
  }
}

// ---------------------------------------------------------------- controls
function clamp(x: number, lo: number, hi: number) {
  return Math.min(hi, Math.max(lo, x));
}

function setQRange(full: boolean) {
  state.qFull = full;
  qFull.checked = full;
  qRange.max = full ? "1" : "0.2";
  qRange.value = String(state.q);
  renderThresholds();
}

function setSelected(n: number) {
  state.sel = clamp(Math.round(n), 0, state.N);
  selRange.value = selNum.value = String(state.sel);
  schedule(false);
}

function setN(n: number) {
  state.N = clamp(Math.round(n), 1, MAX_N);
  nRange.value = nNum.value = String(state.N);
  selRange.max = selNum.max = String(state.N);
  if (state.sel > state.N) setSelected(state.N);
  schedule(true);
}

epsRange.addEventListener("input", () => {
  state.eps = +epsRange.value;
  epsNum.value = fmtQ(state.eps);
  schedule(true);
});
epsNum.addEventListener("input", () => {
  const x = parseFloat(epsNum.value);
  if (!Number.isFinite(x) || x < 0 || x > 1) return;
  state.eps = x;
  epsRange.value = String(x);
  schedule(true);
});
epsNum.addEventListener("change", () => { epsNum.value = fmtQ(state.eps); });

qRange.addEventListener("input", () => {
  state.q = +qRange.value;
  qNum.value = fmtQ(state.q);
  schedule(true);
});
qNum.addEventListener("input", () => {
  const x = parseFloat(qNum.value);
  if (!Number.isFinite(x) || x < 0 || x > 1) return;
  state.q = x;
  if (x > 0.2 && !state.qFull) setQRange(true);
  qRange.value = String(x);
  schedule(true);
});
qNum.addEventListener("change", () => { qNum.value = fmtQ(state.q); });
qFull.addEventListener("change", () => {
  if (!qFull.checked && state.q > 0.2) { qFull.checked = true; return; }
  setQRange(qFull.checked);
});

nRange.addEventListener("input", () => setN(+nRange.value));
nNum.addEventListener("input", () => {
  const x = parseInt(nNum.value, 10);
  if (Number.isFinite(x) && x >= 1 && x <= MAX_N) setN(x);
});
nNum.addEventListener("change", () => { nNum.value = String(state.N); });

selRange.addEventListener("input", () => setSelected(+selRange.value));
selNum.addEventListener("input", () => {
  const x = parseInt(selNum.value, 10);
  if (Number.isFinite(x)) setSelected(x);
});

for (const r of document.querySelectorAll<HTMLInputElement>('input[name="faxis"]')) {
  r.addEventListener("change", () => { state.fAxis = r.value as "auto" | "fixed"; schedule(true); });
}
for (const r of document.querySelectorAll<HTMLInputElement>('input[name="daxis"]')) {
  r.addEventListener("change", () => {
    state.dAxis = r.value as "lin" | "asinh";
    buildCharts(); // scale distribution is fixed at construction
  });
}

// ---------------------------------------------------------------- layout / theme
const ro = new ResizeObserver(() => {
  for (const [u, el] of [[chartF, elF], [chartP, elP], [chartD, elD]] as const) {
    u.setSize({ width: el.clientWidth, height: el.clientHeight - 30 });
  }
});
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => buildCharts());

// ---------------------------------------------------------------- init
computeTrajectory(traj, state.eps, state.q, state.N);
renderMonotone();
renderThresholds();
buildCharts();
renderSummary();
renderInspector();
for (const el of [elF, elP, elD]) ro.observe(el);
