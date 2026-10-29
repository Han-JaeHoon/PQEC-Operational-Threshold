# Type-3 / Type-4 interactive dynamics

A local, browser-only simulator for the round-by-round dynamics of the **Type-3
(textbook 16-CNOT)** and **Type-4 (resynthesized 14-CNOT)** noisy PQEC gadgets. Pick the
protocol at the top, move the input noise `ε` and the per-CNOT noise `q`, and the Bell
fidelity, the four Bell populations and the fidelity change `ΔF_n` update immediately.
Everything is computed in the page; there is no server or Python process. (The folder keeps
its original name, `type3-dynamics`.)

This is a visualization tool for exact maps that are already in the research repositories.
It does not compute new results or thresholds. It evaluates a finite number of rounds in IEEE-754
double precision, so it is **not** an analytic proof or an infinite-round certificate.

## Run (macOS)

Requires Node.js 20 or newer, 22 LTS recommended (`node -v`; install with `brew install node` if needed).
The folder carries its own `.npmrc` (`legacy-peer-deps=true`) to work around an npm 10 crash
when it resolves vitest's optional peer dependencies; no extra flag is needed.

```bash
cd interactive/type3-dynamics
npm install
npm run dev
```

Vite prints a local URL (default `http://localhost:5173/`). Open it in a browser.

Other commands:

| command | what it does |
|--|--|
| `npm test` | runs the numerical tests (Tests A–G, see below) |
| `npm run build` | type-checks and builds a static bundle into `dist/` |
| `npm run preview` | serves the built bundle (`http://localhost:4173/`) |

## The model

The input is `(1−ε)|Φ⁺⟩⟨Φ⁺| + ε I/4 = ¼[II + (1−ε)(XX − YY + ZZ)]`. After every CNOT of
the 16-CNOT gadget, a two-qubit replacement depolarizing channel of strength `q` acts.
Single-qubit gates are ideal. On the invariant plane
`ρ = ¼[II + u XX − u YY + v ZZ]`, with `r = 1 − q`:

```
D_n     = 1 + r⁴ (2u_n² + v_n²)
u_{n+1} = 2 r⁶ u_n (1 + v_n) / D_n
v_{n+1} = r⁵ (1 + r)(u_n² + v_n) / D_n          u_0 = v_0 = 1 − ε

p_Φ⁺ = (1 + 2u + v)/4 = F_n,   p_Φ⁻ = (1 − 2u + v)/4,   p_Ψ⁺ = p_Ψ⁻ = (1 − v)/4
```

This is the exact Type-3 Bell-diagonal map of `RESEARCH_SUMMARY.md` §3.3 (note 01), which
`iterated_noisy_pqec.py` on the `iterated-noisy-pqec` branch verifies against the dense
circuit. Each round is computed from the full map. There is no fidelity-only recurrence and
no assumption that the state stays isotropic.

The map is implemented in [`src/type3map.ts`](src/type3map.ts) (`computeTrajectory`).

### Type 4

Same input, same noise channel after each of the 14 CNOTs of `pqec_resynth_noise.GATES`
(Qiskit-resynthesized `u`/`cx` sequence, unitarily identical to the 16-CNOT gadget), same
parity-weighted read-out. Its exact Bell-diagonal map is the `k = 2` member of the family
whose `k = 4` member is Type 3 (Bell-state-distillation `analytic_exact_map.py`; the Type-4
calculation notes, October 2026). It keeps `y = −x`, so the same `(u, v)` plane applies:

```
D_n     = 1 + r² (2u_n² + v_n²)
u_{n+1} = 2 r⁴ u_n (1 + v_n) / D_n
v_{n+1} = r³ (1 + r)(u_n² + v_n) / D_n          u_0 = v_0 = 1 − ε
```

with the same Bell populations. On general Bell-diagonal inputs this map matches the frozen
circuit-only reference of Bell-state-distillation (`circuit_reference.csv`, 957 Type-4 rows)
to 2.6e-15. It is implemented separately in [`src/type4map.ts`](src/type4map.ts)
(`computeTrajectory4`). Switching protocol keeps ε, q, N, the selected round and the zoom
window.

## Reading `ΔF_n` near the fixed point

`ΔF_n = F_{n+1} − F_n` is evaluated in two ways:

1. **Direct difference** of the two fidelities.
2. **Drift numerator** from the map itself: `ΔF_n = N_n / (4 D_n)` with
   `N_n = 2u(2r⁶(1+v) − D) + r⁵(1+r)(u² + v) − vD` (Type 4: `r⁶, r⁵ → r⁴, r³` and
   `D` with `r²`).

Each route carries a floating-point error bound. A step is labelled **increase** or
**decrease** only when both routes exceed their bounds and agree in sign. Otherwise the step
is marked **unresolved** (amber) and no sign is claimed for it. The inspector shows both
values and both bounds. Steps that are resolved but lie within 10³ of the bound are flagged
as near the precision limit.

As the trajectory converges to its fixed point, `ΔF_n` eventually falls below about `10⁻¹⁵`.
Beyond that round its sign cannot be decided in double precision, and the tool says so.

Both bounds are checked against exact rational arithmetic (`tests/type4map.test.ts`, for
both protocols). For very strong noise (q ≳ 0.5 over hundreds of rounds) the state reaches
I/4 so closely that u, v underflow below 1e-308; there the relative-error model alone does
not bound route 2, so its bound also carries an absolute term of `65 · 2⁻¹⁰⁷⁴`. This changes
no value, no sign and no chart, only that displayed bound in the underflow regime.

## Controls

- **Protocol**: Type 3 — 16 CNOT (default) or Type 4 — 14 CNOT, next to the title.
- **ε, q, N**: sliders with linked number inputs. The q slider spans 0–0.2 by default;
  tick **0–1** for the full range. Typing a value above 0.2 switches the range automatically.
- **q_abs / q_rel ticks**: two short ticks under the q slider mark the analytic thresholds
  of the selected protocol for the current ε (amber `q_abs`, cyan `q_rel`), with their values
  and where q sits relative to them on the line below. Both are the unique root in `(0, 1)`
  of the boundary polynomials of the October 2026 calculation notes (`r = 1−q`, `t = 1−ε`):
  - Type 3: `q_abs`: `3(1+t)r⁶ − (1+t)r⁵ − 3t²r⁴ − 1 = 0`;
    `q_rel`: `(1+t)(2+t)r⁶ − t(1+t)r⁵ − 3t²r⁴ − 1 = 0`
    ([`src/thresholds.ts`](src/thresholds.ts)).
  - Type 4: the same with every power of r lowered by two,
    `q_abs`: `3(1+t)r⁴ − (1+t)r³ − 3t²r² − 1 = 0`;
    `q_rel`: `(1+t)(2+t)r⁴ − t(1+t)r³ − 3t²r² − 1 = 0`
    ([`src/thresholds4.ts`](src/thresholds4.ts)).

  They are defined only for `0 < ε < 2/3`; outside that range the ticks are hidden.
- **Monotonicity badge** (next to the F_n title): `↑ no decrease, n ≤ 500` when no resolved
  step `F_{n+1} < F_n` occurs in rounds 0…500, otherwise `↓ decreases at n = k` and a ring
  on that point of the F_n chart. The check always runs over 500 rounds, independent of N,
  using the same two-route sign rule as the ΔF chart (unresolved steps never count as a
  decrease). It is a finite-round numerical check, not an analytic `q_th`.
- **F axis**: auto-zoom (fits the data) or fixed 0–1.
- **ΔF axis**: linear, or signed log (asinh), so that `10⁻⁴` and `10⁻¹⁰` steps are visible
  together.
- **Zoom rounds**: drag across any chart. All three charts share the round window, and y
  re-fits to it. Double-click to reset.
- **Inspect a round**: click a point, or use the inspector slider.
- **Show/hide populations**: click the legend entries.
- The footer shows the time of the last update (compute + chart update).

## Verification

`npm test` runs `tests/type3map.test.ts`, `tests/reference.test.ts`,
`tests/thresholds.test.ts` (Type 3) and `tests/type4map.test.ts`, `tests/reference4.test.ts`,
`tests/thresholds4.test.ts` (Type 4).

| test | what is checked |
|--|--|
| A | `u₀ = v₀ = 1−ε`, `F₀ = 1 − 3ε/4` for 104 values of ε |
| B | `Σ p_β = 1` (error below 1e-15) for every round, N = 500, on an (ε, q) grid |
| C | every Bell population ≥ −1e-15 on the same grid |
| D | `p_Ψ⁺ = p_Ψ⁻` exactly |
| E | at `q = 0`: populations follow `p → p²/Σp²` (`ρ → ρ²/Tr ρ²`), and `F₁ = (1+3t)²/(4(1+3t²))` |
| F | at `ε = 1`: `F_n = 1/4` exactly, for all q |
| G | Bell populations agree with independent full-circuit references (below), and the ε = 0.1, q = 0.0327 rise-then-decrease appears in the references at the same round |
| T1 | `q_abs`, `q_rel` reproduce Table 1 of the notes to 12 digits (ε = 0.02 … 0.65, and ε → 2/3) |
| T2 | NaN outside `0 < ε < 2/3`; `0 < q_abs < q_rel`, both increasing in ε |
| T3 | against the map itself: at `q_abs − 10⁻⁶` every non-target population (at `q_rel − 10⁻⁶` every ratio `p_i/p_Φ⁺`) decreases over 200 rounds; at `+10⁻⁶` the Φ⁻ quantity rises in round 1 |
| extra | `F₁` vs the repository closed form `pqec_cnot_threshold.F_dec`; the iterates converge to the note-01 fixed point `(u*, v*)`; the two `ΔF` routes agree within their bounds |
| **Type 4** | |
| A–F | as for Type 3, on the Type-4 map; round 0 is bitwise identical to Type 3; the general `(x, y, z)` map keeps `y = −x` and reduces to the `(u, v)` step; first-round `u₁, v₁` of the notes |
| sign | both `ΔF` bounds hold against exact BigInt rational arithmetic, including the underflow regime, for Type 4 and for Type 3; every resolved sign equals the exact sign |
| G4 | Bell populations agree with independent full-circuit references of the 14-CNOT circuit (below) to 1e-12, F₁ with `pqec_resynth_noise.F_resynth`; the references are unit-trace, Bell-diagonal (≤ 1e-9), positive; the circuit normalisation `⟨Z_a⟩` equals `D_n r¹⁰/4` every round; the ε = 0.1, q = 0.0405 rise-then-decrease appears in the references at the same round |
| T1-4 | Type-4 `q_abs`, `q_rel` reproduce the notes' tables to 12 digits (and ε → 2/3, ε = 0.15) |
| T2-4 | NaN outside `0 < ε < 2/3`; `0 < q_abs < q_rel`, both increasing; `r* > 5/6` and `> 4/5` as the proofs require; differ from Type 3 |
| T3-4 | as T3, on the Type-4 map; Φ⁻ is the limiting component |

The references in `verification/reference_type3.json` come from
[`verification/generate_reference.py`](verification/generate_reference.py). It does **not**
use the `(u, v)` closed form. Instead it iterates the full 4×4 two-qubit state through the
noisy 16-CNOT circuit by two routes:

- **R1**: dense 32×32 propagation (`iterated_noisy_pqec.py`, extracted read-only from the
  `iterated-noisy-pqec` branch).
- **R2**: PennyLane `default.mixed` with the repository's own Step-3 executor
  (`verify_analytic_decomposed._fred`). The next state is rebuilt from all 16 parity-weighted
  Pauli correlators.

To regenerate the references (needs numpy, pennylane and matplotlib; nothing in the
repository is modified):

```bash
python interactive/type3-dynamics/verification/generate_reference.py
```

The Type-4 references in `verification/reference_type4.json` come from
[`verification/generate_reference_type4.py`](verification/generate_reference_type4.py), the
same two routes on the 14-CNOT circuit: R1 is `iterated_noisy_pqec` with `"step4"`, R2 is
PennyLane with `pqec_resynth_noise._apply`, and R3 is `pqec_resynth_noise.F_resynth` for F₁.
The cases cover q = 0, small, moderate and strong noise (up to 0.6), ε = 0 and 1, a separable
input, and q just below Type-4 `q_abs(0.15)` and just above `q_rel(0.15)`.

```bash
python interactive/type3-dynamics/verification/generate_reference_type4.py
```

Note on `q → 1`: in the circuit, the normalisation `⟨Z_a⟩ = D_n r¹⁰/4` vanishes as `q → 1`.
It underflows double precision for `q ≳ 0.97` and is exactly 0 at `q = 1`, where the
parity-weighted estimator is 0/0. The simulator shows the closed-form limit there (the
state goes to I/4), for both protocols. The full-circuit references therefore stop at
`q = 0.6`.
