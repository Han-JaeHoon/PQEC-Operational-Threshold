# CNOT-only Noise Threshold — Analysis and Implementation

One round of PQEC on a Bell register, with each controlled-SWAP replaced by its
textbook 8-CNOT decomposition, **single-qubit gates ideal**, and a two-qubit
depolarizing channel after **each CNOT**.  This note separates the **theoretical
analysis** (Part I) from the **numerical implementation and verification**
(Part II).  Theory and circuit agree to `~1e-14`.

---

## 0. Notation and definitions

**Qubits / registers** (5 qubits, wire order `(a, A1, A2, B1, B2) = (0,1,2,3,4)`):

| symbol | definition |
|--------|------------|
| `a` | SWAP-test ancilla (wire 0) |
| `A = (A1, A2)` | **retained** register (wires 1,2) — the purified output lives here |
| `B = (B1, B2)` | **discarded** register (wires 3,4) — traced out |
| `Φ⁺` | `|Φ⁺⟩⟨Φ⁺|`, `|Φ⁺⟩ = (|00⟩+|11⟩)/√2` |

**Input state:**

| symbol | definition |
|--------|------------|
| `ε` | input Bell-state noise (global-depolarizing strength): `ρ_ε = (1−ε)Φ⁺ + ε I₄/4` |
| `p` | *alternative* input parametrization: local depolarizing of strength `p` on **each** Bell qubit |
| `t` | **input Bell-correlation strength**, `ρ_t = ¼[II + t(XX − YY + ZZ)]`; `t = 1−ε` (global input) `= (1−4p/3)²` (local input); `t=1` ideal Bell, `t=0` maximally mixed |
| `R` | two-copy input `ρ_t^A ⊗ ρ_t^B` (registers A and B hold identical copies) |

**Gate noise:**

| symbol | definition |
|--------|------------|
| `q ≡ ε₂` | per-CNOT **two-qubit replacement depolarizing** strength: `D_q(σ) = (1−q)σ + q·I₄/4` (q = replacement probability) |
| `s ≡ v` | **survival factor** `s = 1−q = 1−ε₂`; on a noisy CNOT pair every non-identity Pauli component is scaled by `s` (`v` is the same quantity used in the code) |

Single-qubit gates (`H, T, T†`, and the two ancilla Hadamards) are **ideal** throughout.

**Derived / output quantities:**

| symbol | definition |
|--------|------------|
| `σ_out` | final 5-qubit state after the two decomposed CSWAPs and the last ancilla `H` |
| `τ_A` | parity-weighted register-A operator `τ_A = Tr_{a,B}[(Z_a ⊗ I_A ⊗ I_B) σ_out]` |
| `Q` | parity **denominator** `Q = Tr τ_A = ⟨Z_a ⊗ I_A⟩` |
| `N_Φ` | Bell-projector **numerator** `N_Φ = Tr(Φ⁺ τ_A) = ⟨Z_a ⊗ Φ_A⟩` |
| `F_PQEC` | purified Bell fidelity `= N_Φ / Q` |
| `F_bare` | bare (no-QEC) fidelity `= (1+3t)/4` |
| `ρ_eff` | normalized effective state `= τ_A / Q` |
| `c_⊥, c_z` | effective-state Bell-diagonal correlators: `ρ_eff = ¼[II + c_⊥(XX−YY) + c_z ZZ]` |
| `ε₂* ≡ q_th` | **CNOT-noise threshold** (`F_PQEC = F_bare`) |
| `K₂ = α(t)` | small-noise slope `F_PQEC ≈ F_ideal − K₂ ε₂` |

---

## Part I — Theoretical analysis

### I.1 Circuit and input

Each Fredkin uses the textbook decomposition (Clifford+T Toffoli):

```
CSWAP(c; x, y) = CNOT(x→y) · Toffoli(c, y; x) · CNOT(x→y)     (8 CNOTs each)
```

The full one-round SWAP test is
`H_a → CSWAP(a; A1,B1) → CSWAP(a; A2,B2) → H_a`, i.e. **16 CNOTs** total.
After the first Hadamard the state is `σ₀ = |+⟩⟨+|_a ⊗ R`.

### I.2 Noise model

A two-qubit **replacement** depolarizing channel acts right after each CNOT on the
two qubits `(i,j)` it touched:

```
D_q^{(ij)}(σ) = (1−q) σ + q [ I_{ij}/4 ⊗ Tr_{ij}(σ) ] .
```

Single-qubit gates (`H, T, T†`) and the two ancilla Hadamards are ideal.

### I.3 Pauli-propagation rule

Writing states in the 5-qubit Pauli basis, the channel acts diagonally:

```
D_q^{(ij)}(P) = P            if P_i = P_j = I,
             = s·P           otherwise      (s = 1−q).
```

Only Pauli components that are identity on the noisy pair survive unscaled; every
component with a non-identity Pauli on that pair is damped by `s`. Propagating the
16 noisy CNOTs with the ideal Clifford+T conjugation rules gives the exact output.

### I.4 Parity-weighted operator (exact)

Keeping only the sector with `Z` on the ancilla and identity on `B`,

```
τ_A(t,q) = s¹⁰/16 · [ (1 + 3 s⁴ t²) II
                    + 2 s⁶ t(1+t) (XX − YY)
                    + s⁵ (1+s) t(1+t) ZZ ] .
```

The common leading factor is `s¹⁰` (not `s⁸`): the two data–data CNOTs also damp
the parity-carrying components.

### I.5 Denominator, numerator, fidelity (exact closed forms)

```
Q(t,q)   = Tr τ_A       = s¹⁰/4 · (1 + 3 s⁴ t²)
N_Φ(t,q) = Tr(Φ⁺ τ_A)   = s¹⁰/16 · [ 1 + 3 s⁴ t² + s⁵(1+5s) t(1+t) ]
```

The common `s¹⁰` cancels in the ratio:

```
F_PQEC(t,q) = N_Φ/Q = [ 1 + 3 s⁴ t² + s⁵(1+5s) t(1+t) ] / [ 4 (1 + 3 s⁴ t²) ]
            = ¼ [ 1 + s⁵(1+5s) t(1+t) / (1 + 3 s⁴ t²) ] ,     s = 1−q .
```

(Equivalently, with `v = s = 1−ε₂`, this is `¼[1 + t(1+t)(v⁵+5v⁶)/(1+3v⁴t²)]`, the
form used in `pqec_cnot_threshold.py`.)

### I.6 Effective state is anisotropic

```
ρ_eff = ¼ [ II + c_⊥ (XX − YY) + c_z ZZ ],
c_⊥ = 2 s⁶ t(1+t) / (1 + 3 s⁴ t²),      c_z = s⁵(1+s) t(1+t) / (1 + 3 s⁴ t²),
c_z − c_⊥ = s⁵(1−s) t(1+t) / (1 + 3 s⁴ t²) ≥ 0 .
```

So a **noisy CNOT turns the isotropic (Werner) Bell input into an anisotropic
Bell-diagonal state**: the `Z` correlation is preserved better than `X/Y`. At `q=0`,
`c_⊥ = c_z` and the state is isotropic again. This is a genuine new feature of the
decomposition-level model.

### I.7 Threshold and small-noise slope

`F_bare(t) = (1+3t)/4`. The CNOT-noise threshold `s_th = 1−q_th` solves

```
s_th⁵ (1 + 5 s_th)(1 + t) = 3 (1 + 3 s_th⁴ t²) .
```

Small-`q` expansion:

```
F_PQEC ≈ F_ideal(t) − α(t) q ,   α(t) = t(1+t)(35 + 33 t²) / (4(1 + 3 t²)²) ,   α(1) = 17/8 .
```

### I.8 Sampling overhead

The parity signal is `Q(t,q) = s¹⁰/4 (1 + 3 s⁴ t²)`, so the shot overhead is

```
N_shot(t,q)/N_shot(t,0) ~ [Q(t,0)/Q(t,q)]² = s⁻²⁰ [ (1 + 3 t²)/(1 + 3 s⁴ t²) ]² .
```

### I.9 Consistency limits

- `q=0` (`s=1`): `F_PQEC = (1+3t)²/(4(1+3t²))` — ideal one-round PQEC.
- `t=0`: `F_PQEC = ¼` — no recoverable Bell information.

### I.10 Contrast with the 3-qubit global-depolarizing model (Step 3a)

| | 3-qubit global depol after each whole CSWAP | two-qubit depol after each **CNOT** |
|--|---------------------------------------------|-------------------------------------|
| numerator/denominator | both `∝ (1−g)²` → **cancel** | `∝ s¹⁰ ×` different `t`-dependent factors → **do not cancel** |
| effective state | stays `∝ ρ_t²` (isotropic) | anisotropic (`c_z > c_⊥`) |
| operational threshold | **none** (`F` independent of `g`) | **finite** `ε₂*` |
| cost | sampling only, `N ~ (1−g)⁻⁴` | bias **and** sampling |

---

## Part II — Numerical implementation and verification

### II.1 Circuit (genuine PennyLane `default.mixed`)

- `verify_analytic_decomposed.py` — native gates: `_c2` (CNOT + 2-qubit depol `e2`),
  the Clifford+T `_tof` (single-qubit gates ideal), and `_fred(q,a,b,e2)` =
  `CNOT(b→a)·Toffoli(q,a;b)·CNOT(b→a)`. The full gadget:
  `QubitDensityMatrix → H_a → _fred → _fred → H_a`, reading `(⟨Z_a⊗Φ_A⟩, ⟨Z_a⊗I_A⟩)`.
- `pqec_cnot_threshold.py` — the analytic closed forms `F_dec`, `Q_denom`, `N_num`,
  `c_perp`, `c_z`, the threshold `eps2_star`, `eff_correlators_circuit`, and the
  figure.

### II.2 Noise-channel convention (Kraus)

The per-CNOT channel is `global_depol_kraus(ε₂)` (from `noisy_bell_state.py`):

```
K_I = √(1 − 15 ε₂/16) · I ,     K_P = √(ε₂/16) · P   (15 non-identity 2-qubit Paulis) ,
```

which is exactly `D_q(σ) = (1−ε₂)σ + ε₂ I₄/4`, i.e. **ε₂ is the replacement
probability** — the same convention as the 3-qubit global-depolarizing channel of
Step 3a.

### II.3 Verification (circuit vs Part I), all to `~1e-14`

| quantity | check |
|----------|-------|
| `F_PQEC` | circuit `= ¼[1 + s⁵(1+5s)t(1+t)/(1+3s⁴t²)]` |
| `Q`, `N_Φ` | match `s¹⁰/4(1+3s⁴t²)` and `s¹⁰/16[…]` |
| `c_⊥, c_z` | match; `c_z − c_⊥ > 0` for `q>0`, `= 0` for `q=0` |

Worked point `t = 0.6 (ε=0.4)`, `q = 0.12 (s=0.88)`:
`Q = 0.11471949`, `N_Φ = 0.07629940`, `F_PQEC = 0.66509535`,
`c_⊥ = 0.541161`, `c_z = 0.578059` (anisotropy `+0.036897`) — circuit = analytic.

### II.4 Threshold results

| input `ε` | 0.05 | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 |
|-----------|------|------|------|------|------|------|------|
| `ε₂*` | 0.017 | 0.033 | 0.061 | 0.085 | 0.103 | 0.117 | 0.126 |
| `16·ε₂*` (budget) | 0.27 | 0.53 | 0.98 | 1.35 | 1.65 | 1.87 | 2.01 |

`ε₂*` grows with input noise and stays **above realistic hardware CNOT error
(`~10⁻²`)** for `ε ≳ 0.03`, so one PQEC round with a genuinely noisy decomposed
CSWAP still yields a net fidelity gain. Figure: `pqec_cnot_threshold.png`.

### II.5 Circuit diagrams (`draw_cnot_noise.py`)

- `circuit_cswap_decomp.png` — one CSWAP decomposition
  (`CSWAP(0;1,2) = CNOT(2→1)·Toffoli(0,1;2)·CNOT(2→1)`).
- `circuit_swaptest_decomp.png` — full SWAP test (ideal), barrier-separated
  `state prep | CSWAP₁ | CSWAP₂ | final H`.
- `circuit_swaptest_cnot_noise.png` — same, with a 2-qubit depolarizing channel
  after each CNOT.

### II.6 Run

```bash
python pqec_cnot_threshold.py       # closed forms, threshold table, anisotropy check, figure
python verify_analytic_decomposed.py  # A/B/F_dec + K2 slope vs circuit
python draw_cnot_noise.py           # the three circuit diagrams
```
