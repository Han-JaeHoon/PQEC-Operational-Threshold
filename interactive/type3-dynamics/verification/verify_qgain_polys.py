"""
Exact checks of the Nov 6 note "Exact Algebraic Determination of the Asymptotic Fidelity-Gain
Threshold" used by src/qgain.ts (sympy, rational arithmetic only).

  G1  E_{s,eps}(v) (Eqs. 14-17) is the elimination of u_B = (3 theta - v)/2 from Eq. (9).
  G2  the closed form (18) equals the Sylvester resultant Res_v(H_s, E_{s,eps}), symbolically in eps.
  G3  the note's own certificate (Appendix A): R_P(s, 2/3) factorization, one root in (0, 1),
      isolating interval, and the decisive signs of Q_P and V_P.
  G4  the isolating intervals of s_ent lie strictly between those of s_SN and s_bi (Eq. 25).

Run: python interactive/type3-dynamics/verification/verify_qgain_polys.py
"""
import sympy as S

s, v, eps = S.symbols("s v epsilon")
Q = S.Rational
th = 1 - eps

ITEMS = [
    ("T3", (2*s**6, 2*s**6, s**5*(1+s), s**5*(1+s), 2*s**4, s**4),
     112*s**8 - 16*s**7 - 19*s**6 + 10*s**5 - 3*s**4 - 49*s**2 + 14*s - 1,
     24*s**7 + 24*s**6 - 11*s**5 + 5*s**4 - 9*s - 9, 28*s**6 - 4*s**5 - 3*s**4 - 9,
     Q(8708, 10000), Q(8709, 10000), -s**20*(s+1)**2*(3*s**4+1),
     (Q(869, 1000), Q(870, 1000)), (Q(881, 1000), Q(882, 1000))),
    ("T4", (2*s**4, 2*s**4, s**3*(1+s), s**3*(1+s), 2*s**2, s**2),
     112*s**6 - 16*s**5 - 19*s**4 + 10*s**3 - 52*s**2 + 14*s - 1,
     24*s**5 + 24*s**4 - 11*s**3 + 5*s**2 - 9*s - 9, 28*s**4 - 4*s**3 - 3*s**2 - 9,
     Q(815, 1000), Q(816, 1000), -s**12*(s+1)**2*(3*s**2+1),
     (Q(810, 1000), Q(811, 1000)), (Q(821, 1000), Q(822, 1000))),
    ("4Q", (2*s**5, s**4*(1+s), 2*s**5, s**3*(1+s), s**3 + s**5, s**5),
     271*s**11 + 176*s**10 - 128*s**9 - 76*s**8 + 10*s**7 - 184*s**6 + 76*s**5 + 80*s**4 - s**3 + 8*s**2 - 36*s - 4,
     23*s**5 + 18*s**4 - 2*s**3 + 3*s - 18, 16*s**5 + 9*s**4 - 4*s**3 - 9,
     Q(820, 1000), Q(821, 1000), -4*s**19,
     (Q(819, 1000), Q(820, 1000)), (Q(824, 1000), Q(825, 1000))),
]


def bound(poly, a, b, lower):
    return sum(c*(a**k if (c > 0) == lower else b**k) for (k,), c in S.Poly(poly, s).terms())


for name, (a, be, g, d, et, b), R, QQ, VV, lo, hi, factor, sn, bi in ITEMS:
    H = (g*b + et*a)*v**2 + (et*(be - d) - g*a)*v + g*(1 - be)
    E = (g - 4*a)*v**2 + (-6*g*th - 4*(be - d))*v + 9*g*th**2
    # G1
    uB = (3*th - v)/2
    assert S.cancel(4*g*uB**2 - 4*v*(a*v + be - d) - E) == 0, name
    # G2
    h2, h1, h0 = H.coeff(v, 2), H.coeff(v, 1), H.coeff(v, 0)
    e2, e1, e0 = E.coeff(v, 2), E.coeff(v, 1), E.coeff(v, 0)
    closed = (h2*e0 - h0*e2)**2 - (h2*e1 - h1*e2)*(h1*e0 - h0*e1)
    assert S.expand(S.resultant(H, E, v) - closed) == 0, name
    # G3
    at23 = S.factor(closed.subs(eps, Q(2, 3)))
    assert S.cancel(at23 - factor*R) == 0, name
    assert S.Poly(R, s).count_roots(0, 1) == 1, name
    assert R.subs(s, lo) * R.subs(s, hi) < 0, name
    assert bound(QQ, lo, hi, True) > Q(1, 100), name
    assert bound(VV, lo, hi, False) < -Q(1, 10), name
    # G4
    assert sn[1] <= lo and hi <= bi[0], name
    print(name, "OK")
print("ALL OK")
