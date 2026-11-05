# Type-3 / Type-4 / 4Q / BBPSSW interactive dynamics

A local, browser-only simulator for the round-by-round dynamics of the **Type-3
(textbook 16-CNOT)** and **Type-4 (resynthesized 14-CNOT)** noisy PQEC gadgets and of the
**4Q (four-qubit, 5-CNOT, postselected)** noisy Bell-state distillation protocol, and of **Bennett's
BBPSSW** recurrence purification. Pick the
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

## BBPSSW

The fourth protocol is Bennett et al. (1996) Protocol A, with `Φ⁺` as the target: on two identical
pairs `[A1, B1, A2, B2]`, input twirl, `CNOT(A1→A2)`, `D_q` on `(A1, A2)`, `CNOT(B1→B2)`, `D_q` on
`(B1, B2)`, measure `(A2, B2)` in Z, keep `(A1, B1)` on `00` or `11`, output twirl. Everything except the
two `D_q` is ideal. Because of the twirls the state is `(F, h, h, h)` with `h = (1 − F)/3`, so one
scalar recurrence is exact ([`src/bbpsswmap.ts`](src/bbpsswmap.ts)), with `s = (1 − q)²`:

```
F_{n+1} = [s(F² + h²) + (1 − s)/8] / P_succ,   P_succ = s[(F + h)² + 4h²] + (1 − s)/2,   F_0 = 1 − 3ε/4
```

`P_succ^(n)` is shown per round in the inspector and is not accumulated. No threshold is shown for this
protocol. `F_n` belongs to the physical output after every earlier round succeeded, while Type 3/4 show
a reconstructed effective state, so equal `(ε, q)` does not mean equal resources or equal total noise.

`verification/bbpssw_stage1/` is the Stage 1 derivation and frozen circuit, copied unchanged. It is a
standalone reference, not one of this repository's earlier frozen baselines. Its `verify.py` rewrites
its own result files, so run it on a copy. `verification/generate_reference_bbpssw.py` turns it and
a second circuit route into `verification/reference_bbpssw.json`, which `tests/referencebbpssw.test.ts`
checks. `q_gain^(1)` appears there as a check on the sign of the first step, never in the page.

## Regimes: q_bi and q_SN

Each protocol has two proved constants, independent of `ε`, with `q_m` used only to identify the
low attractor. They live in [`src/regimes.ts`](src/regimes.ts) with their polynomials in
`s = 1 − q`; `tests/regimes.test.ts` checks each against bisection to `1e-9`.

| | `q_bi` | `q_SN` | `q_m` |
|--|--|--|--|
| Type 3 | 0.1182907629 | 0.1305792680 | 0.1187285384 |
| Type 4 | 0.1789326522 | 0.1894166559 | 0.1808274866 |
| 4Q | 0.1758332653 | 0.1806697260 | 0.1808274866 |

- `q < q_bi`: every `ε ∈ (0,1)` converges to the same Bell-connected fixed point `P_B`.
- `q_bi < q < q_SN`: bistable. `ε < ε_c(q)` goes to `P_B`, `ε > ε_c(q)` goes to `Z(q)`. Just below
  `q_SN`, `P_B` itself has `F < 1/2`, so reaching `P_B` is not the same as a fidelity gain.
- `q > q_SN`: every `ε ∈ (0,1)` converges to the same separable point `Z(q)`. `Z(q)` is the
  boundary point `B*` for `q < q_m`, which only 4Q reaches beyond `q_SN`, and the mixed point
  (`F = 1/4`) otherwise.

The q slider draws both constants as lines over three tinted bands, and the line below it names
the regime of the current `q` with the limit of the current orbit. That limit comes from a separate
orbit of up to 200,000 rounds, not from the charts, which stop at 500. It carries the transverse
component as a logarithm, so underflow cannot trap it on the invariant line, and a tiny component
that is still growing is never read as convergence. Running out of rounds shows "미수렴". This
orbit is uncertified and is labelled so; it never decides the attractor.

The regime shown comes from `regimeReport` in `src/regimes.ts`: the exact sign of the polynomial at
the float64 `q`, in BigInt arithmetic, not a comparison with the rounded constants. A float64 `q` is
a dyadic rational and none of these polynomials has a rational root in its bracket, so `q` is never
exactly on a threshold. Within `1e-9` of `q_bi` or `q_SN` the panel adds a boundary note.

## The threshold ε_c(q)

The "Asymptotic attractor" panel states what the theorems give for the current `(q, ε)`, separately
from the round-by-round charts. In the bistable window it shows a certified interval
`ε_L < ε_c(q) < ε_U`, compares `ε` with it, and marks it on the ε slider.

`ε ≤ ε_L` predicts `P_B`, `ε ≥ ε_U` predicts `Z(q)`, and in between the panel says the basin is not
determined. No theorem is applied at `ε = 0`, `ε = 1`, `q = 0` or `q = 1`.

**Algorithm** ([`src/epsc.ts`](src/epsc.ts)). Theorem 2.1 of the Nov 5 calculation notes: if
`T^n p_ε` is strictly above the saddle `S` in the order (`t_n > t_S`, `m_n < m_S`) then `ε < ε_c`, and
if strictly below then `ε > ε_c`. `S` is enclosed from its closed form, Eqs. 16-22. A bisection between
the closed-form bounds `L_0`, `U_0` of Eqs. 39-42 tests one input at a time. Exact root isolation
of `t_n(ε) = t_S` is not used, because the degree grows like `2^n`.

**What is certified.** Every float64 result is rounded outward. `T` preserves the order, so a box maps
inside the images of its two corners. Each endpoint of the interval carries its own certificate, and
so does every intermediate interval shown while the computation runs. The certificate relies on the
proved order property and global convergence (Facts 1-4 of the notes). The bisection only chooses
which inputs to test.

**What is not certified.** The `float64 orbit` rows, and the charts. The midpoint `ε_c ≈` is only the
centre of the interval. Point values of `P_B` and `Z(q)` are float64 values printed to 6 digits.
No table, interpolation or fit is used; `tests/epsc.test.ts` checks that `src/` reads no table.

**Accuracy and speed.** Measured widths are about `1e-13` through most of the window and `1e-12` to
`1e-11` near its ends, with a typical run of 1 to 40 ms in a Web Worker. Closer than about `1e-7` to
`q_bi`, convergence takes more than the 2,000,000 rounds allowed per input: the panel then keeps the
valid wider interval and says "계산 미완료". Going further would need more than float64.

**Checks.**

```bash
npm test
```

```bash
python verification/verify_epsc_polys.py
```

```bash
python verification/generate_reference_epsc.py
```

`npm test` includes all 57 rows of the research tables (`verification/epsc_research/`, regression
only). `verify_epsc_polys.py` checks the polynomials and the map identity exactly. The last
script needs python-flint and mpmath: it certifies the simulator's endpoints at 18 off-grid `q`, using
a different certificate (Proposition 4 of the `ε_c` note) in 256-bit ball arithmetic.

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

### 4Q

Four qubits, two copies of the input: the 5 CNOTs of Bell-state-distillation
`pqec_distill.circuit.CNOT_SEQUENCE`, each followed by the same replacement channel `D_q`,
then an ideal H on q3; (m₃, m₄) = (0, 0) is postselected and the first pair kept. Its exact
Bell-diagonal map (`analytic_exact_map.analytic_map_4q`; the 4Q calculation notes, October
2026) is not of the Type-3/4 form, and the Bell-isotropic input keeps `y = −z` instead of
`y = −x`. On that plane `ρ = ¼[II + a XX − b YY + b ZZ]`:

```
D_n     = 1 + r⁵ a_n² + (r⁵ + r³) b_n²          P_succ = D_n / 4
a_{n+1} = r³ [(1 + r) a_n + 2 r² b_n²] / D_n
b_{n+1} = r⁴ b_n [(1 + r) + 2 r a_n] / D_n       a_0 = b_0 = 1 − ε

p_Φ⁺ = (1 + a + 2b)/4 = F_n,   p_Ψ⁺ = (1 + a − 2b)/4,   p_Φ⁻ = p_Ψ⁻ = (1 − a)/4
```

So in 4Q mode the coinciding pair is **Φ⁻ = Ψ⁻** (the dashed Ψ⁻ line lies on Φ⁻), not
Ψ⁺ = Ψ⁻. The inspector names the coordinates `a_n = ⟨XX⟩`, `b_n = ⟨ZZ⟩`. The map matches
the frozen circuit-only reference (`circuit_reference.csv`, 957 4Q rows). It is implemented
in [`src/q4map.ts`](src/q4map.ts) (`computeTrajectory4Q`). Unlike Type 3/4, the normalisation
`P_succ ≥ 1/4` never vanishes, so `q = 1` is regular (one round gives I/4).

## Reading `ΔF_n` near the fixed point

`ΔF_n = F_{n+1} − F_n` is evaluated in two ways:

1. **Direct difference** of the two fidelities.
2. **Drift numerator** from the map itself: `ΔF_n = N_n / (4 D_n)` with
   `N_n = 2u(2r⁶(1+v) − D) + r⁵(1+r)(u² + v) − vD` (Type 4: `r⁶, r⁵ → r⁴, r³` and
   `D` with `r²`; 4Q: `N_n = A + 2B − (a + 2b)D` with `A = r³(1+r)a + 2r⁵b²`,
   `B = r⁴b(1 + r + 2ra)`).

Each route carries a floating-point error bound. A step is labelled **increase** or
**decrease** only when both routes exceed their bounds and agree in sign. Otherwise the step
is marked **unresolved** (amber) and no sign is claimed for it. The inspector shows both
values and both bounds. Steps that are resolved but lie within 10³ of the bound are flagged
as near the precision limit.

As the trajectory converges to its fixed point, `ΔF_n` eventually falls below about `10⁻¹⁵`.
Beyond that round its sign cannot be decided in double precision, and the tool says so.

Both bounds are checked against exact rational arithmetic (`tests/type4map.test.ts` for
Type 3 and Type 4, `tests/q4map.test.ts` for 4Q). The 4Q bounds come from an explicit
operation count (header of `src/q4map.ts`) and use `32·u` for both routes, with the
absolute route-2 term `128·2⁻¹⁰⁷⁴/(4D) + 2⁻¹⁰⁷⁴`. For very strong noise (q ≳ 0.5 over hundreds of rounds) the state reaches
I/4 so closely that u, v underflow below 1e-308; there the relative-error model alone does
not bound route 2, so its bound also carries an absolute term of `65 · 2⁻¹⁰⁷⁴`. This changes
no value, no sign and no chart, only that displayed bound in the underflow regime.

## Controls

- **Protocol**: Type 3 — 16 CNOT (default), Type 4 — 14 CNOT or 4Q — 5 CNOT, next to the
  title.
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
  - 4Q (limiting component Ψ⁺ for both):
    `q_abs`: `2(1+t−t²)r⁵ + r⁴ − (1+t²)r³ − 1 = 0`;
    `q_rel`: `1 − r⁴ + t(t+1)r³ + (2t²−3t−1)r⁵ = 0`
    ([`src/thresholds4q.ts`](src/thresholds4q.ts)).

  For Type 3 and Type 4 both are defined only for `0 < ε < 2/3`. For 4Q, `q_abs` is defined
  for `0 < ε < 2/3` and `q_rel` for the whole interior `0 < ε < 1` (the q_rel note proves it
  there), so for `2/3 ≤ ε < 1` only the `q_rel` tick and value are shown. Outside a
  threshold's range its tick is hidden.
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
`tests/thresholds4.test.ts` (Type 4) and `tests/q4map.test.ts`, `tests/reference4q.test.ts`,
`tests/thresholds4q.test.ts` (4Q).

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
| **4Q** | |
| A–F | same input state as Type 3; Σp = 1, p ≥ 0, `p_Φ⁻ = p_Ψ⁻` exactly and `a, b ≥ 0`, N = 500; the general `(x, y, z)` map keeps `y = −z` and reduces to the `(a, b)` step; `q = 0` gives `ρ²/Tr ρ²`; `ε = 1` stays at 1/4; `q = 1` gives I/4 in one round; first-round `a₁, b₁` of the q_abs note |
| sign | both `ΔF` bounds hold against exact BigInt rational arithmetic (including underflow, q up to 1); every resolved sign equals the exact sign |
| G4Q | Bell populations agree with independent full-circuit references of the 5-CNOT circuit (below) to 1e-12; the references are unit-trace, Bell-diagonal, positive, with `p_Φ⁻ = p_Ψ⁻`; the circuit's `P_succ` equals `D_n/4` every round, including q = 1 |
| T1-4Q | `q_abs` (Table 1 of the q_abs note) and `q_rel` (Sec. 17 of the q_rel note, including ε = 0.75, 0.9) to 12 digits; ε → 2/3 and ε → 1 limits; the q_rel maximum 0.158732684144 at ε = 0.767499473760; the small-ε series; the q-form quintic (rel Eq. 108) |
| T2-4Q | `q_abs` NaN outside `0 < ε < 2/3`, `q_rel` defined on `0 < ε < 1`; `0 < q_abs < q_rel`; `q_abs` increasing with `r* > 6/7`; `q_rel` not monotone on `(0, 1)`; above Type 3 everywhere, above Type 4 only for ε < 0.451 (abs) / 0.483 (rel) |
| T3-4Q | as T3, on the 4Q map, `q_rel` also for ε up to 0.99; Ψ⁺ is the limiting component |

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

The 4Q references in `verification/reference_4q.json` come from
[`verification/generate_reference_4q.py`](verification/generate_reference_4q.py), from a
checkout of Bell-state-distillation (default: next to this repository, `--bsd` to override).
R1 is that package's dense 16×16 simulator (`noisy_full_channel_dm` + `postselect_branch`);
R2 is PennyLane `default.mixed` written in the script from the gate list only, with the next
state rebuilt from all 16 expectations `⟨P ⊗ |00⟩⟨00|⟩`. Both record `P_succ`. The 18 cases
include q just below `q_abs(0.15)`, just above `q_rel(0.15)`, between the two at ε = 0.1,
separable inputs ε = 0.75 and 0.9 around `q_rel`, and q = 1.

```bash
python interactive/type3-dynamics/verification/generate_reference_4q.py
```

Note on `q → 1`: in the circuit, the normalisation `⟨Z_a⟩ = D_n r¹⁰/4` vanishes as `q → 1`.
It underflows double precision for `q ≳ 0.97` and is exactly 0 at `q = 1`, where the
parity-weighted estimator is 0/0. The simulator shows the closed-form limit there (the
state goes to I/4), for both protocols. The full-circuit references therefore stop at
`q = 0.6`. This does not apply to 4Q, whose `P_succ = D_n/4 ≥ 1/4`.
