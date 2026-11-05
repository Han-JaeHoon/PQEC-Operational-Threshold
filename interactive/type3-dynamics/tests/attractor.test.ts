import { describe, expect, it } from "vitest";
import { attractorReport, fmtOutward, zText } from "../src/attractor";
import { computeEpsc } from "../src/epsc";
import { EpscClient, type EpscReply, type WorkerLike } from "../src/epscClient";
import { orbitLimit, THRESHOLDS, type RegimeProtocol } from "../src/regimes";

const COORDS: Record<RegimeProtocol, [string, string]> = { type3: ["u", "v"], type4: ["u", "v"], q4: ["a", "b"] };
const PROTOS: RegimeProtocol[] = ["type3", "type4", "q4"];

function report(p: RegimeProtocol, q: number, eps: number) {
  const b = computeEpsc(p, q);
  const status = b ? { state: "final" as const, bracket: b } : { state: "failed" as const, bracket: null };
  return attractorReport(p, COORDS[p], q, eps, status, orbitLimit(p, eps, q));
}
const text = (r: ReturnType<typeof report>) => r.rows.map((x) => `${x.k}: ${x.v}`).join("\n");

describe("attractor panel", () => {
  for (const p of PROTOS) {
    const t = THRESHOLDS[p];
    const mid = (t.qBi.q + t.qSN.q) / 2;
    it(`${p}: entering and leaving the bistable window`, () => {
      expect(report(p, t.qBi.q - 0.01, 0.5).attractor).toBe("P_B");
      expect(report(p, t.qSN.q + 0.01, 0.5).attractor).toBe("Z");
      expect(text(report(p, t.qBi.q - 0.01, 0.5))).toContain("Global attractor");
      expect(text(report(p, t.qSN.q + 0.01, 0.5))).toContain("Global attractor");
      expect(report(p, mid, 0.1).attractor).toBe("P_B");
      expect(report(p, mid, 0.999).attractor).toBe("Z");
      expect(text(report(p, mid, 0.1))).toContain("Bistable regime");
      expect(text(report(p, mid, 0.1))).toContain("Certified interval");
    });
    it(`${p}: around eps_c`, () => {
      const b = computeEpsc(p, mid)!;
      expect(report(p, mid, b.lo).attractor).toBe("P_B");
      expect(report(p, mid, b.hi).attractor).toBe("Z");
      const inside = report(p, mid, (b.lo + b.hi) / 2);
      expect(inside.attractor).toBe("undetermined");
      expect(text(inside)).toContain("basin 미확정");
    });
    it(`${p}: no theorem at eps = 0, eps = 1, q = 0, q = 1`, () => {
      for (const eps of [0, 1]) {
        for (const q of [t.qBi.q - 0.01, mid, t.qSN.q + 0.01]) {
          const r = report(p, q, eps);
          expect(r.attractor).toBe("out-of-range");
          expect(text(r)).toContain("정리 범위");
        }
      }
      for (const q of [0, 1]) expect(report(p, q, 0.5).attractor).toBe("out-of-range");
    });
    it(`${p}: while nothing is known it predicts nothing`, () => {
      const r = attractorReport(p, COORDS[p], mid, 0.5, { state: "pending", bracket: null }, orbitLimit(p, 0.5, mid));
      expect(r.attractor).toBe("pending");
      expect(r.bracket).toBeNull();
      const f = attractorReport(p, COORDS[p], mid, 0.5, { state: "failed", bracket: null }, orbitLimit(p, 0.5, mid));
      expect(f.attractor).toBe("pending");
      expect(f.rows.map((x) => x.v).join()).toContain("계산 실패");
    });
    it(`${p}: at the quoted threshold values it shows the boundary state`, () => {
      for (const q of [t.qBi.q, t.qSN.q]) expect(text(report(p, q, 0.5))).toContain("경계");
    });
    it(`${p}: never calls an attractor a purification or a gain by itself`, () => {
      for (const q of [t.qBi.q - 0.01, mid, t.qSN.q + 0.01]) {
        expect(text(report(p, q, 0.5))).not.toMatch(/purif/i);
      }
    });
  }
  it("names the line and transverse coordinate correctly in Z(q)", () => {
    const q = THRESHOLDS.type3.qBi.q - 0.01;
    expect(zText("type3", COORDS.type3, q)).toContain("v* =");
    expect(zText("type3", COORDS.type3, q)).toContain("u = 0");
    expect(zText("q4", COORDS.q4, 0.17)).toContain("a* =");
    expect(zText("q4", COORDS.q4, 0.17)).toContain("b = 0");
    expect(zText("type4", COORDS.type4, 0.3)).toContain("혼합점");
  });
  it("prints intervals outward", () => {
    expect(Number(fmtOutward(0.123456, 3, "down"))).toBeLessThan(0.123456);
    expect(Number(fmtOutward(0.123456, 3, "up"))).toBeGreaterThan(0.123456);
    expect(fmtOutward(0.5, 3, "down")).toBe("0.499");
    expect(fmtOutward(0.5, 3, "up")).toBe("0.501");
  });
});

describe("eps_c client", () => {
  function setup() {
    const sent: { id: number; protocol: RegimeProtocol; q: number }[] = [];
    const w: WorkerLike = { postMessage: (m) => sent.push(m as never), onmessage: null, onerror: null };
    let updates = 0;
    const c = new EpscClient(w, () => { updates++; });
    const reply = (r: Omit<EpscReply, "bracket"> & { lo?: number; hi?: number }) => w.onmessage!({
      data: { ...r, bracket: { protocol: r.protocol, q: r.q, lo: r.lo ?? 0.4, hi: r.hi ?? 0.6, roundsLo: 0, roundsHi: 0, totalRounds: 0, done: r.final, stop: r.final ? "resolution" : "" } },
    });
    return { c, w, sent, reply, updates: () => updates };
  }
  it("a late result for an earlier q is never shown for the current q", () => {
    const { c, sent, reply } = setup();
    c.request("q4", 0.178);
    c.request("q4", 0.179);
    expect(sent.length).toBe(2);
    reply({ id: sent[0].id, protocol: "q4", q: 0.178, final: true });
    expect(c.status("q4", 0.179).state).toBe("pending");
    expect(c.status("q4", 0.178).state).toBe("final");
    expect(c.status("type4", 0.178).state).toBe("pending");
  });
  it("caches finals and does not ask twice, re-asks if a job was abandoned", () => {
    const { c, sent, reply } = setup();
    c.request("q4", 0.178);
    c.request("q4", 0.178);
    expect(sent.length).toBe(1);
    reply({ id: 1, protocol: "q4", q: 0.178, final: true });
    c.request("q4", 0.178);
    expect(sent.length).toBe(1);
    c.request("q4", 0.179);
    reply({ id: 2, protocol: "q4", q: 0.179, final: false });
    expect(c.status("q4", 0.179).state).toBe("running");
    c.request("q4", 0.18);
    c.request("q4", 0.179); // superseded job: computed again
    expect(sent.length).toBe(4);
  });
  it("a worker failure reads as failed, not as a value", () => {
    const { c, w } = setup();
    w.onerror!({});
    expect(c.status("q4", 0.178).state).toBe("failed");
  });
});
