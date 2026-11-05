/** Writes the brackets of src/epsc.ts for the q list on stdin, for generate_reference_epsc.py. */
import { computeEpsc } from "../src/epsc";
import type { RegimeProtocol } from "../src/regimes";
import { readFileSync } from "node:fs";

const list = JSON.parse(readFileSync(0, "utf8")) as [RegimeProtocol, number][];
console.log(JSON.stringify(list.map(([protocol, q]) => {
  const b = computeEpsc(protocol, q);
  return b && { protocol, q, lo: b.lo, hi: b.hi, stop: b.stop };
})));
