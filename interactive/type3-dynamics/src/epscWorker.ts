/**
 * Runs the eps_c certifier off the main thread.  Only the newest request is worked on: each one
 * supersedes the previous at the next tested input.  Every message names its own (protocol, q),
 * and each bracket it carries is already certified, so a partial bracket is safe to display.
 */
import { epscBrackets } from "./epsc";
import type { RegimeProtocol } from "./regimes";

export interface EpscRequest { id: number; protocol: RegimeProtocol; q: number }

const PROGRESS_MS = 60;
let latest = 0;

// Yield to the event loop quickly, so that a newer request can arrive between tested inputs.
const channel = new MessageChannel();
let resume: (() => void) | null = null;
channel.port1.onmessage = () => { resume?.(); };
const tick = () => new Promise<void>((r) => { resume = r; channel.port2.postMessage(0); });

async function run({ id, protocol, q }: EpscRequest) {
  const g = epscBrackets(protocol, q);
  let last = -Infinity;
  for (;;) {
    const r = g.next();
    if (r.done) {
      postMessage({ id, protocol, q, bracket: r.value, final: true });
      return;
    }
    const now = performance.now();
    if (now - last > PROGRESS_MS) {
      last = now;
      postMessage({ id, protocol, q, bracket: r.value, final: false });
    }
    await tick();
    if (id !== latest) return;
  }
}

onmessage = (e: MessageEvent<EpscRequest>) => {
  latest = e.data.id;
  void run(e.data);
};
