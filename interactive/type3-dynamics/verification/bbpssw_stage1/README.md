# BBPSSW noisy one-round map: independent circuit-level verification

Date: 2026-11-03. All file contents are standalone. This environment did **not**
contain the user's Git repository: the immutable/frozen circuit specification is
stored in `frozen_circuit_spec.json` and realized in `frozen_circuit_reference.py`.
No preexisting repository freeze, commit, or tests were modified or compared.

## Scope and conventions

Reference: Bennett, Brassard, Popescu, Schumacher, Smolin and Wootters,
*Purification of Noisy Entanglement and Faithful Teleportation via Noisy
Channels*, Phys. Rev. Lett. **76**, 722 (1996), Protocol A and Eq. (7).
https://doi.org/10.1103/PhysRevLett.76.722

The paper uses a `Psi-` singlet and maps it to `Phi+` by a unilateral Pauli before
BXOR. Here `Phi+` is the target *from the start*, so the trivial relabeling
operations at the beginning/end are omitted. Everything is expressed in the
ordered Bell basis `[Phi+, Phi-, Psi+, Psi-]`.

Frozen noisy circuit specification:

1. Optionally pre-twirl each input pair (mandatory for literal BBPSSW on arbitrary
   input; the initial epsilon-isotropic state is already pre-twirled).
2. Prepare two identical input pairs in the 4-qubit order `[A1,B1,A2,B2]`.
3. Apply `CNOT(A1 -> A2)` followed by replacement depolarization `D_q` acting
   on the same local qubits `(A1,A2)`.
4. Apply `CNOT(B1 -> B2)` followed by the same `D_q` on `(B1,B2)`.
5. Measure the target pair `(A2,B2)` in the computational (Z) basis, retaining
   **both** equal-outcome events `00` and `11`. Discard `01` and `10`.
6. Retain the unmeasured source pair `(A1,B1)`, condition on success, and
   perform ideal `Phi+`-preserving twirling on the retained pair.

Twirling is the uniform single-qubit Clifford average of `U (x) U*` on a pair.
It is ideal, as are state preparation and readout. The noise is only on the two
local CNOTs, with no cross-party two-qubit gate. Each `D_q` acts as

`D_q(sigma) = (1-q)*sigma + q * [ (I_pair/4) (x) Tr_pair(sigma) ]`.

The independent simulator instead implements this channel using the equal
average of *sixteen* two-qubit Pauli conjugations. The original twirl is
implemented using *twenty-four* explicitly generated Clifford conjugations.

## Analytic derivation

Take an arbitrary Bell-diagonal input with population vector
`p=(a,b,c,d)`, entries nonnegative and summing to one. Define

```
w = (a*a+b*b, 2*a*b, c*c+d*d, 2*c*d)
P0 = (a+b)**2 + (c+d)**2
s = (1-q)**2
```

The ideal BXOR postselection gives the unnormalized source-pair state with Bell
population weights `w` and success probability `P0`. This follows from the
Bell-label bit/phase-error propagation through bilateral CNOT:
`bit_source' = bit_source`, `bit_target' = bit_target XOR bit_source`,
`phase_source' = phase_source XOR phase_target`.

For Bell-diagonal two-qubit input, the Alice and Bob local reduced states of
the *two-copy* input are both `I_4/4`. After either local CNOT, replacing its
two qubits gives `I_16/16`: the unreplaced party is already locally maximally
mixed and the replaced party becomes locally maximally mixed. Any subsequent
CNOT or replacement channel preserves `I_16/16`. Therefore the full noisy
four-qubit post-gate state is exactly

```
sigma_q = s * sigma_ideal + (1-s) * I_16/16.
```

A maximally mixed 4-qubit state passes the equal-Z-readout test with probability
`1/2` and yields unnormalized retained state `I_4/8`. Consequently,

```
unnormalized populations w_q = s*w + ((1-s)/8)*(1,1,1,1)
Psuccess = s*P0 + (1-s)/2
pre-output-twirl populations = w_q/Psuccess.
```

The *output* twirl preserves `Phi+` fidelity and redistributes all other Bell
populations equally. Define

```
F' = [s*(a*a+b*b)+(1-s)/8] / [s*((a+b)**2+(c+d)**2)+(1-s)/2]
post-output-twirl populations = (F', (1-F')/3, (1-F')/3, (1-F')/3).
```

**Important: original BBPSSW includes input twirling for arbitrary inputs.**
For an input not already isotropic, the initial twirl first transforms
`(a,b,c,d)` into `(a, (1-a)/3, (1-a)/3, (1-a)/3)`; the formula above must
be evaluated on these pre-twirled probabilities to recover the literal
Bennett Protocol A. The original complete noisy map is therefore

```
F = a
h = (1-F)/3
G_q(F) = [s*(F*F+h*h) + (1-s)/8] /
         [s*((F+h)**2 + (2*h)**2) + (1-s)/2]
p' = (G_q(F), (1-G_q(F))/3, (1-G_q(F))/3, (1-G_q(F))/3)
Psuccess = s*((F+h)**2+4*h*h) + (1-s)/2.
```

When the input is already isotropic, the initial twirl is the identity. For
repeated BBPSSW every post-round output is isotropic; the next round can use
the scalar `G_q(F)`. The successful raw state is physical (unlike a virtual
ratio estimator), but discarded failed attempts incur a resource overhead.

### Epsilon-isotropic initial input

Start from `rho_e = (1-epsilon) Phi+ + epsilon I_4/4` so that
`F0 = 1 - 3*epsilon/4`, `h=epsilon/4`. The exact noisy one-round map is

```
s = (1-q)**2
A = 1 - 3*epsilon/2 + 5*epsilon**2/8
P0 = 1 - epsilon + epsilon**2/2
Psuccess = s*P0 + (1-s)/2
F1 = (s*A+(1-s)/8)/Psuccess
epsilon1 = 4*(1-F1)/3
```

At `q=0` this is the original Bennett Eq. (7). At `q=1`, `Psuccess=1/2`
and `F1=1/4` for every input. The ideal gain region begins at `F0>1/2`,
which is equivalent to `0<epsilon<2/3`.

The one-round *strict* fidelity-gain threshold for `0<epsilon<2/3` is

```
q_gain_first(epsilon) = 1 - sqrt( 3/(3 + 2*epsilon - 3*epsilon**2) ).
```

This comes from factoring

```
Psuccess*(F1-F0) =
    (1-epsilon)/8 * (s*(3+2*epsilon-3*epsilon**2)-3).
```

This is **not** presented as the Rohde noise threshold nor as any all-round
monotonicity theorem. Degenerate endpoints are handled separately.

## Verification

Run from this folder:

```
python verify.py
python -m pytest -q test_bbpssw.py
```

The independent reference uses explicit 16x16 complex matrices, bit-permuted
CNOTs, independent Pauli-averaged depolarizing channels, the 00/11 projective
readout, and a physical Clifford twirl. The analytic module does not call
any of these circuit routines.

- 146 test Bell-diagonal input probability vectors across 13 q strengths:
  **1,898 independent one-round comparisons**.
- 672 additional repeated-round steps on isotropic inputs.
- `q=0`: agreement with Bennett's original scalar recurrence.
- `q=1`: maximally mixed output and postselection probability `1/2`.
- independent verification of twirl invariance and channel equivalence on
  generic full-rank four-qubit states (unit tests).
- 3 negative controls deliberately alter CNOT direction, CNOT noise, or
  postselection accept condition; all cause large detected disagreements.
- no physics assumptions or fitted coefficients enter the independent
  circuit simulator beyond the frozen circuit specification.

See `verification_results.json` for complete measured max errors and the
frozen specification SHA256. All comparison errors are near float64 roundoff
(about 1e-15, below 1e-12 tolerances).

## Notes and limits

- **No modification of the user's Git repository** was performed. The
  present folder is ready to be copied into it as a reproducible Stage-1
  reference but is not itself an existing Git-frozen baseline.
- The formula `sigma_q = s*sigma_ideal+(1-s)*I_16/16` assumes Bell-diagonal
  identical input pairs (maximally mixed local marginals). It is not claimed
  for arbitrary non-Bell-diagonal two-qubit inputs.
- The noisy map is a *conditional*, nonlinear, trace-normalized map; its
  success probability matters for resource comparisons.
- An ideal twirl uses correlated local single-qubit operations/classical
  randomness. Their imperfections are not included here.
- Running an analytic recurrence after each success describes the input to
  subsequent successful rounds, not a deterministic fixed-depth circuit with
  an automatically guaranteed successful branch.
