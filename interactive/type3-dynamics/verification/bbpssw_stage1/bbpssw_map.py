"""Analytical BBPSSW recurrence under per-local-CNOT replacement depolarization.

Convention: two identical Bell-diagonal Phi+-target pairs, BXOR A1->A2 and
B1->B2, local computational-basis measurements of A2/B2, accept equal outcomes.
Each local CNOT is followed by the replacement depolarizing channel on its two
wires. Twirling after success is noiseless and preserves Phi+ fidelity.
No pretwirl is applied to the arbitrary Bell-diagonal input supplied to a round.
"""

from __future__ import annotations
from math import isfinite
from typing import Sequence


def _check_probabilities(p: Sequence[float], q: float) -> tuple[float, float, float, float]:
    if len(p) != 4:
        raise ValueError('Bell population vector must contain 4 probabilities')
    v = tuple(float(z) for z in p)
    if not all(isfinite(z) for z in v) or any(z < -1e-13 for z in v):
        raise ValueError('invalid Bell probabilities')
    if abs(sum(v) - 1.0) > 1e-12:
        raise ValueError('Bell probabilities must sum to 1')
    if not isfinite(q) or not 0.0 <= q <= 1.0:
        raise ValueError('q must lie in [0,1]')
    return v


def round_untwirled(p: Sequence[float], q: float) -> tuple[tuple[float, ...], float]:
    """Return normalized pre-twirl populations and probability of success."""
    a, b, c, d = _check_probabilities(p, q)
    survival = (1.0 - q) ** 2
    w = (a*a + b*b, 2*a*b, c*c + d*d, 2*c*d)
    p0 = (a+b)**2 + (c+d)**2
    p_success = survival*p0 + (1-survival)*0.5
    background = (1-survival)/8.0
    return tuple((survival*z+background)/p_success for z in w), p_success


def round_twirl(p: Sequence[float], q: float) -> tuple[tuple[float, ...], float]:
    """Return normalized post-twirl Phi+-isotropic populations and success prob."""
    untwirled, ps = round_untwirled(p, q)
    f = untwirled[0]
    e = (1-f)/3.0
    return (f, e, e, e), ps


def round_original_bbpssw(p: Sequence[float], q: float) -> tuple[tuple[float, ...], float]:
    """Literal BBPSSW: input twirl, noisy BXOR+postselect, output twirl.

    For arbitrary Bell-diagonal input, the input twirl keeps a = fidelity and
    maps other populations to (1-a)/3 before running the two CNOTs.
    """
    a, b, c, d = _check_probabilities(p, q)
    isotropic = (a, (1-a)/3.0, (1-a)/3.0, (1-a)/3.0)
    return round_twirl(isotropic, q)


def isotropic_input(epsilon: float) -> tuple[float, float, float, float]:
    if not 0 <= epsilon <= 1:
        raise ValueError('epsilon must lie in [0, 1]')
    return (1-3*epsilon/4, epsilon/4, epsilon/4, epsilon/4)


def first_round_isotropic(epsilon: float, q: float) -> tuple[float, float]:
    """Closed scalar recurrence: return F_1 and P_success."""
    if not 0 <= epsilon <= 1 or not 0 <= q <= 1:
        raise ValueError('epsilon and q must lie in [0,1]')
    survival=(1-q)**2
    pure_w=1-1.5*epsilon+5*epsilon**2/8
    pure_ps=1-epsilon+epsilon**2/2
    success=survival*pure_ps+(1-survival)/2
    return (survival*pure_w+(1-survival)/8)/success, success


def trajectory(epsilon: float, q: float, rounds: int) -> list[dict[str, float | int | tuple]]:
    if not isinstance(rounds, int) or rounds < 0:
        raise ValueError('rounds must be a nonnegative integer')
    p = isotropic_input(epsilon)
    out = []
    for k in range(rounds+1):
        rec = {'round': k, 'populations': p, 'fidelity': p[0]}
        if k < rounds:
            p_new, ps = round_original_bbpssw(p, q)
            rec['success_probability'] = ps
            rec['delta_fidelity'] = p_new[0]-p[0]
            p=p_new
        out.append(rec)
    return out
