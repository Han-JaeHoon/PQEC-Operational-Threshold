/**
 * Asynchronous access to eps_c brackets, with a cache.  Results are stored under the (protocol, q)
 * named in the worker's own message and read back only under the current (protocol, q), so a late
 * result for another q can fill the cache but can never be shown for the current one.
 */
import type { EpscBracket } from "./epsc";
import type { RegimeProtocol } from "./regimes";

export interface WorkerLike {
  postMessage(msg: unknown): void;
  onmessage: ((e: { data: EpscReply }) => void) | null;
  onerror: ((e: unknown) => void) | null;
}

export interface EpscReply {
  id: number;
  protocol: RegimeProtocol;
  q: number;
  bracket: EpscBracket | null;
  final: boolean;
}

export type EpscStatus =
  /** No certified bracket yet. */
  | { state: "pending"; bracket: null }
  /** Certified so far, still improving. */
  | { state: "running"; bracket: EpscBracket }
  | { state: "final"; bracket: EpscBracket }
  /** The fixed points could not be enclosed, or the worker failed. */
  | { state: "failed"; bracket: null };

const CACHE_MAX = 4000;
const key = (p: RegimeProtocol, q: number) => `${p}|${q}`;

export class EpscClient {
  private finals = new Map<string, EpscBracket | null>();
  private partial = new Map<string, EpscBracket>();
  private nextId = 1;
  private runningKey: string | null = null;
  private broken = false;

  constructor(private worker: WorkerLike, private onUpdate: () => void) {
    worker.onmessage = (e) => this.receive(e.data);
    worker.onerror = () => { this.broken = true; this.runningKey = null; this.onUpdate(); };
  }

  /** Ask for (protocol, q) unless it is known or already being computed. */
  request(protocol: RegimeProtocol, q: number) {
    const k = key(protocol, q);
    if (this.broken || this.finals.has(k) || this.runningKey === k) return;
    this.runningKey = k;
    this.worker.postMessage({ id: this.nextId++, protocol, q });
  }

  /** What is known for exactly this (protocol, q). */
  status(protocol: RegimeProtocol, q: number): EpscStatus {
    const k = key(protocol, q);
    if (this.finals.has(k)) {
      const b = this.finals.get(k)!;
      return b ? { state: "final", bracket: b } : { state: "failed", bracket: null };
    }
    if (this.broken) return { state: "failed", bracket: null };
    const b = this.partial.get(k);
    return b ? { state: "running", bracket: b } : { state: "pending", bracket: null };
  }

  private receive({ protocol, q, bracket, final }: EpscReply) {
    const k = key(protocol, q);
    if (final) {
      if (this.finals.size >= CACHE_MAX) this.finals.delete(this.finals.keys().next().value as string);
      this.finals.set(k, bracket);
      this.partial.delete(k);
      if (this.runningKey === k) this.runningKey = null;
    } else if (bracket && !this.finals.has(k)) {
      if (this.partial.size >= CACHE_MAX) this.partial.delete(this.partial.keys().next().value as string);
      this.partial.set(k, bracket);
    }
    this.onUpdate();
  }
}
