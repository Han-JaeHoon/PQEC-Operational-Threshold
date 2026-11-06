"""
Exact checks of the Nov 6 analytic note "Exact Algebraic Calculation of the Asymptotic
Target-Fidelity Noise Threshold" used by src/qtarget.ts (sympy, rational arithmetic only).

  T1  H_s = s^{k+1} (A v^2 + B v + C) and E_{s,lambda} = s^{k+1} (D v^2 + E v + G) (Eqs. 21-22),
      and the 4Q normalization H_s = s^5 (...), E = 2 s^3 (...) (Eqs. 32-33).
  T2  (AG - CD)^2 - (AE - BD)(BG - CE) is the Sylvester resultant of the normalized quadratics,
      and P_k = -N / (1+s)^2, P_4Q = N / s^3 are exact polynomial divisions (Eqs. 24, 34).
  T3  the explicit coefficient expansions c_0..c_4 (Eqs. 26-30) and d_0..d_4 (Eqs. 36-40).
  T4  degrees: 14, 10, 11 in s; 4 in lambda (Eqs. 31, 41).
  T5  lambda = 1: P(s, 1) = (3 s^4 + 1) R_T3, (3 s^2 + 1) R_T4, -R_4Q (Eq. 46), with R_P the q_ent
      polynomials of src/qgain.ts; s = 1: Eqs. (51)-(52).
  T6  dF_B/ds at s = 1: 17/8, 7/4, 1 (Section 6.3).
  T7  the raw resultant (Eq. 19) of the coefficient table equals the normalized one up to the
      positive factors used in src/qtarget.ts: Type 3/4: s^{4(k+1)} * N, 4Q: 4 s^16 * N.

Run: python interactive/type3-dynamics/verification/verify_qtarget_polys.py
"""
import sympy as S

s, v, L = S.symbols("s v L")

RAW = {
    "T3": (2*s**6, 2*s**6, s**5*(1+s), s**5*(1+s), 2*s**4, s**4),
    "T4": (2*s**4, 2*s**4, s**3*(1+s), s**3*(1+s), 2*s**2, s**2),
    "4Q": (2*s**5, s**4*(1+s), 2*s**5, s**3*(1+s), s**3 + s**5, s**5),
}
R = {
    "T3": 112*s**8 - 16*s**7 - 19*s**6 + 10*s**5 - 3*s**4 - 49*s**2 + 14*s - 1,
    "T4": 112*s**6 - 16*s**5 - 19*s**4 + 10*s**3 - 52*s**2 + 14*s - 1,
    "4Q": 271*s**11 + 176*s**10 - 128*s**9 - 76*s**8 + 10*s**7 - 184*s**6 + 76*s**5 + 80*s**4
          - s**3 + 8*s**2 - 36*s - 4,
}


def explicit(name):
    """Eqs. (25)-(30) and (35)-(40), typed in from the note."""
    if name == "4Q":
        d4 = 4*s**7*(2*s**2 + 1)**2
        d3 = -4*s**5*(2*s**2 + 1)*(3*s**4 + 1)
        d2 = -s**3*(131*s**8 + 80*s**7 + 28*s**6 + 40*s**5 + 14*s**4 - 80*s**3 - 4*s**2 - 40*s - 1)
        d1 = -4*s**2*(25*s**4 - 8*s**2 - 5)*(s**5 + s**4 - 1)
        d0 = -4*(s**5 + s**4 - 1)*(8*s**6 - 9*s**5 - 21*s**4 + 3*s**2 + 9*s + 1)
        return d4*L**4 + d3*L**3 + d2*L**2 + d1*L + d0
    h = s**(4 if name == "T3" else 2)
    c4 = -h**2*(5*s + 1)**2
    c3 = 4*h**2*(5*s + 1)*(s**3 + s**2 - s + 1)
    c2 = 2*h*(14*h*s**6 + 12*h*s**5 + 80*h*s**4 + 44*h*s**3 - 12*h*s**2 + 12*h*s - 6*h - 45*s**2 - 14*s - 1)
    c1 = 4*h*(2*h*s**2 - 1)*(7*s**4 + 6*s**3 + 12*s**2 - 8*s - 5)
    c0 = -(2*h*s**2 - 1)*(42*h*s**4 + 36*h*s**3 - 30*h*s**2 + 48*h*s - 24*h - 49*s**2 + 14*s - 1)
    return c4*L**4 + c3*L**3 + c2*L**2 + c1*L + c0


def main():
    for name, (al, be, ga, de, et, b) in RAW.items():
        H = (ga*b + et*al)*v**2 + (et*(be - de) - ga*al)*v + ga*(1 - be)
        E = ga*(L - v)**2 - 4*v*(al*v + be - de)
        if name == "4Q":
            A, B, C = 2*s**3*(1 + 2*s**2), -s*(1 + 3*s**4), 2*(1 - s**4 - s**5)
            D, Ee, G = -3*s**2, 2*(1 - s**2) - 2*s**2*L, s**2*L**2
            hs, es, fac, raw_fac = s**5, 2*s**3, s**3, 4*s**16
            deg = 11
        else:
            k = 4 if name == "T3" else 2
            A, B, C = s**k*(5*s + 1), 2*s**k*(s - 1 - s**2 - s**3), (1 + s)*(1 - 2*s**(k + 2))
            D, Ee, G = 1 - 7*s, 4*(1 - s) - 2*(1 + s)*L, (1 + s)*L**2
            hs, es, fac, raw_fac = s**(k + 1), s**(k + 1), -(1 + s)**2, s**(4*(k + 1))
            deg = 14 if name == "T3" else 10
        Hn, En = A*v**2 + B*v + C, D*v**2 + Ee*v + G
        assert S.expand(H - hs*Hn) == 0 and S.expand(E - es*En) == 0, (name, "T1")
        N = (A*G - C*D)**2 - (A*Ee - B*D)*(B*G - C*Ee)
        assert S.expand(N - S.resultant(Hn, En, v)) == 0, (name, "T2 resultant")
        P = S.cancel(N/fac)
        assert S.expand(N - fac*P) == 0 and S.Poly(P, s, L).domain == S.ZZ, (name, "T2 division")
        assert S.expand(P - explicit(name)) == 0, (name, "T3")
        assert S.degree(P, s) == deg and S.degree(P, L) == 4, (name, "T4", S.degree(P, s))
        lam1 = (3*s**4 + 1)*R[name] if name == "T3" else (3*s**2 + 1)*R[name] if name == "T4" else -R[name]
        assert S.expand(P.subs(L, 1) - lam1) == 0, (name, "T5 lambda=1")
        endsign = 12 if name == "4Q" else -12
        assert S.expand(P.subs(s, 1) - endsign*(L - 3)*(L + 1)**2*(3*L - 1)) == 0, (name, "T5 s=1")
        dv = -S.diff(H, s).subs({s: 1, v: 1})/S.diff(H, v).subs({s: 1, v: 1})
        W = v*(al*v + be - de)/ga
        du = (S.diff(W, s).subs({s: 1, v: 1}) + S.diff(W, v).subs({s: 1, v: 1})*dv)/2
        slope = {"T3": S.Rational(17, 8), "T4": S.Rational(7, 4), "4Q": S.Integer(1)}[name]
        assert S.simplify((dv + 2*du)/4 - slope) == 0, (name, "T6")
        # Eq. (19) on the raw coefficients
        h2, h1, h0 = ga*b + et*al, et*(be - de) - ga*al, ga*(1 - be)
        e2, e1, e0 = ga - 4*al, -2*ga*L - 4*(be - de), ga*L**2
        Rraw = (h2*e0 - h0*e2)**2 - (h2*e1 - h1*e2)*(h1*e0 - h0*e1)
        assert S.expand(Rraw - raw_fac*N) == 0, (name, "T7")
        print(name, "deg_s", deg, "deg_L 4: PASS")
    print("verify_qtarget_polys: all checks PASS")


if __name__ == "__main__":
    main()
