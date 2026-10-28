# Type-3 interactive dynamics

A local, browser-only simulator for the round-by-round dynamics of the **Type-3
(textbook 16-CNOT) noisy PQEC gadget**. Move the input noise `ε` and the per-CNOT noise `q`,
and the Bell fidelity, the four Bell populations and the fidelity change `ΔF_n` update
immediately. Everything is computed in the page; there is no server or Python process.

This is a visualization tool for an exact map that is already in the repository. It does
not compute new results or thresholds. It evaluates a finite number of rounds in IEEE-754
double precision, so it is **not** an analytic proof or an infinite-round certificate.

## Run (macOS)

Requires Node.js 18 or newer (`node -v`; install with `brew install node` if needed).

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

## Reading `ΔF_n` near the fixed point

`ΔF_n = F_{n+1} − F_n` is evaluated in two ways:

1. **Direct difference** of the two fidelities.
2. **Drift numerator** from the map itself: `ΔF_n = N_n / (4 D_n)` with
   `N_n = 2u(2r⁶(1+v) − D) + r⁵(1+r)(u² + v) − vD`.

Each route carries a floating-point error bound. A step is labelled **increase** or
**decrease** only when both routes exceed their bounds and agree in sign. Otherwise the step
is marked **unresolved** (amber) and no sign is claimed for it. The inspector shows both
values and both bounds. Steps that are resolved but lie within 10³ of the bound are flagged
as near the precision limit.

As the trajectory converges to its fixed point, `ΔF_n` eventually falls below about `10⁻¹⁵`.
Beyond that round its sign cannot be decided in double precision, and the tool says so.

## Controls

- **ε, q, N**: sliders with linked number inputs. The q slider spans 0–0.2 by default;
  tick **0–1** for the full range. Typing a value above 0.2 switches the range automatically.
- **F axis**: auto-zoom (fits the data) or fixed 0–1.
- **ΔF axis**: linear, or signed log (asinh), so that `10⁻⁴` and `10⁻¹⁰` steps are visible
  together.
- **Zoom rounds**: drag across any chart. All three charts share the round window, and y
  re-fits to it. Double-click to reset.
- **Inspect a round**: click a point, or use the inspector slider.
- **Show/hide populations**: click the legend entries.
- The footer shows the time of the last update (compute + chart update).

## Verification

`npm test` runs `tests/type3map.test.ts` and `tests/reference.test.ts`.

| test | what is checked |
|--|--|
| A | `u₀ = v₀ = 1−ε`, `F₀ = 1 − 3ε/4` for 104 values of ε |
| B | `Σ p_β = 1` (error below 1e-15) for every round, N = 500, on an (ε, q) grid |
| C | every Bell population ≥ −1e-15 on the same grid |
| D | `p_Ψ⁺ = p_Ψ⁻` exactly |
| E | at `q = 0`: populations follow `p → p²/Σp²` (`ρ → ρ²/Tr ρ²`), and `F₁ = (1+3t)²/(4(1+3t²))` |
| F | at `ε = 1`: `F_n = 1/4` exactly, for all q |
| G | Bell populations agree with independent full-circuit references (below), and the ε = 0.1, q = 0.0327 rise-then-decrease appears in the references at the same round |
| extra | `F₁` vs the repository closed form `pqec_cnot_threshold.F_dec`; the iterates converge to the note-01 fixed point `(u*, v*)`; the two `ΔF` routes agree within their bounds |

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
