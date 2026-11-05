"""
epsc_core.py -- certified computation of the separating noise eps_c(q) in the bistable window.

Maps in (t,m) coordinates (see the q_SN note, Section 1):
  4Q   : Delta = 4 s^4 t (2m(1-s)+3s+1),  Sigma = 2[4 m^2 al + 4 m be + t^2 (3 s^5 + s^3) + alp],
         L3 = 4 m^2 alp + 4 m be + t^2 (s^3 - s^5) + al,   al = 1-s^3-s^4+s^5, alp = 1+s^3+s^4+s^5, be = 1-s^5
  T3/4 : e1=s^B, e2=s^(A-1), e3=s^A, ga = 1+e1-e2-e3, gap = 1+e1+e2+e3,
         Delta = 16 e3 t, Sigma = 2[4 m^2 ga + 4 m (1-e1) + t^2 (2e1+e2+e3) + gap],
         L3 = 4 m^2 gap + 4 m (1-e1) + t^2 (2e1-e2-e3) + ga
  T(t,m) = (Delta/Sigma, L3/Sigma).
The same function works with Fraction (exact), flint.arb (rigorous balls) and mpmath (floats).
"""
from fractions import Fraction as Fr
from flint import arb, ctx

AB = {'T4': (4, 2), 'T3': (6, 4)}


def T(kind, s, t, m):
    if kind == '4Q':
        al = 1 - s**3 - s**4 + s**5; alp = 1 + s**3 + s**4 + s**5; be = 1 - s**5
        D = 4*s**4*t*(2*m*(1 - s) + 3*s + 1)
        S = 2*(4*m*m*al + 4*m*be + t*t*(3*s**5 + s**3) + alp)
        L = 4*m*m*alp + 4*m*be + t*t*(s**3 - s**5) + al
    else:
        A, B = AB[kind]
        e1, e2, e3 = s**B, s**(A - 1), s**A
        ga = 1 + e1 - e2 - e3; gap = 1 + e1 + e2 + e3
        D = 16*e3*t
        S = 2*(4*m*m*ga + 4*m*(1 - e1) + t*t*(2*e1 + e2 + e3) + gap)
        L = 4*m*m*gap + 4*m*(1 - e1) + t*t*(2*e1 - e2 - e3) + ga
    return D/S, L/S


def quad(kind, s):
    """Coefficients A2, A1, A0 of Phi (q_SN note, Section 4) and tau(m)."""
    if kind == '4Q':
        A2 = 8*s**2*(5*s**5 - 2*s**4 + 2*s**3 + s + 2)
        A1 = -16*s**2*(3*s**5 + s**4 + s**3 - 1)
        A0 = 2*s**2*(1 - s**2)*(s**3 + 2*s**2 - s + 2)
        al = 1 - s**3 - s**4 + s**5
        tau = lambda m: (5*s**5 + s**4 - s**3 - 1 - 4*m*(1 - s**4) - 4*m*m*al)/(s**3*(3*s**2 + 1))
    else:
        A, B = AB[kind]
        e1, e2, e3 = s**B, s**(A - 1), s**A
        ga = 1 + e1 - e2 - e3; gap = 1 + e1 + e2 + e3
        if kind == 'T4':
            A2 = 8*s*(1 + s)*(3*s**2 + 1); A1 = -8*s*(2*s**5 + 2*s**4 + 5*s**3 + s**2 - s - 1)
            A0 = 2*s*(1 - s)*(4*s**4 + 8*s**3 + s**2 + 2*s + 1)
        else:
            A2 = 8*s*(1 + s)*(3*s**4 + 1); A1 = -8*s*(2*s**7 + 2*s**6 + 5*s**5 + s**4 - s - 1)
            A0 = 2*s*(1 - s)*(4*s**6 + 8*s**5 + s**4 + 2*s**3 + 2*s**2 + 2*s + 1)
        tau = lambda m: (8*e3 - gap - 4*m*(1 - e1) - 4*m*m*ga)/(2*e1 + e2 + e3)
    return A2, A1, A0, tau


def line_cb(kind, s):
    if kind == '4Q':
        return s**3*(1 + s), s**5
    A, B = AB[kind]
    return s**(A - 1)*(1 + s), s**B


def fixed_points_arb(kind, s):
    """Rigorous enclosures of S=(t_S,m_S) (saddle, root m_+), P_B (root m_-) and m_Z."""
    A2, A1, A0, tau = quad(kind, s)
    disc = A1*A1 - 4*A2*A0
    assert disc > 0
    mp_ = (-A1 + disc.sqrt())/(2*A2); mm_ = (-A1 - disc.sqrt())/(2*A2)
    tS = tau(mp_).sqrt(); tB = tau(mm_).sqrt()
    c, b = line_cb(kind, s)
    if c > 1:
        w = ((c - 1)/b).sqrt(); mZ = (1 - w)/(2*(1 + w))
    elif c < 1:
        mZ = arb(1)/2
    else:
        raise ValueError("q = q_m not handled (c undecided)")
    return (tS, mp_), (tB, mm_), mZ
