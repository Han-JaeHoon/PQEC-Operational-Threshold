"""
Independent reference values of q_target(F_target) for tests/qtarget.test.ts.

Different arithmetic from src/qtarget.ts: the decimal F_target is an exact rational, F_B(q) is
evaluated from the closed form (larger root v_B of H_s, positive u_B) in 256-bit arb balls, and the
bisection runs on exact rationals until the width is below 1e-30.  Each final interval is then
checked by an exact Sturm count to hold exactly one root of the note's normalized polynomial
P_P(s, 4 F_target - 1) in s = 1 - q (verify_qtarget_polys.py checks those polynomials).

Run from interactive/type3-dynamics (needs python-flint and sympy):
    python verification/generate_reference_qtarget.py
"""
import json
import os
from fractions import Fraction as Fr

import sympy as S
from flint import arb, ctx

ctx.prec = 256
HERE = os.path.dirname(os.path.abspath(__file__))
TARGETS = ["0.5000001", "0.51", "0.6", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95", "0.99",
           "0.999", "0.9999", "0.99999"]
# lower end of the isolating interval of s_ent (Nov 6 note, Appendix A.1)
S_ENT_LO = {"type3": Fr(8708, 10000), "type4": Fr(815, 1000), "q4": Fr(820, 1000)}

s_, L_ = S.symbols("s L")


def coeffs(p, s):
    if p == "q4":
        return 2*s**5, s**4*(1 + s), 2*s**5, s**3*(1 + s), s**3 + s**5, s**5
    a, b = (6, 4) if p == "type3" else (4, 2)
    c = s**(a - 1)*(1 + s)
    return 2*s**a, 2*s**a, c, c, 2*s**b, s**b


def fb(p, q):
    s = 1 - arb(q.numerator)/q.denominator
    al, be, ga, de, et, b = coeffs(p, s)
    h2, h1, h0 = ga*b + et*al, et*(be - de) - ga*al, ga*(1 - be)
    vb = (-h1 + (h1*h1 - 4*h2*h0).sqrt())/(2*h2)
    ub = (vb*(al*vb + be - de)/ga).sqrt()
    return (1 + vb + 2*ub)/4


def sign(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    raise ArithmeticError("ball contains 0")


def normalized(p):
    if p == "q4":
        A, B, C = 2*s_**3*(1 + 2*s_**2), -s_*(1 + 3*s_**4), 2*(1 - s_**4 - s_**5)
        D, E, G = -3*s_**2, 2*(1 - s_**2) - 2*s_**2*L_, s_**2*L_**2
        fac = s_**3
    else:
        k = 4 if p == "type3" else 2
        A, B, C = s_**k*(5*s_ + 1), 2*s_**k*(s_ - 1 - s_**2 - s_**3), (1 + s_)*(1 - 2*s_**(k + 2))
        D, E, G = 1 - 7*s_, 4*(1 - s_) - 2*(1 + s_)*L_, (1 + s_)*L_**2
        fac = -(1 + s_)**2
    return S.cancel(((A*G - C*D)**2 - (A*E - B*D)*(B*G - C*E))/fac)


def main():
    out = []
    for p in ("type3", "type4", "q4"):
        P = normalized(p)
        for t in TARGETS:
            F = Fr(t)
            Fa = arb(F.numerator)/F.denominator
            lo, hi = Fr(0), 1 - S_ENT_LO[p]
            assert sign(fb(p, hi) - Fa) == -1
            while hi - lo > Fr(1, 10**30):
                mid = (lo + hi)/2
                if sign(fb(p, mid) - Fa) > 0:
                    lo = mid
                else:
                    hi = mid
            lam = S.Rational(4*F.numerator - F.denominator, F.denominator)
            Ps = S.Poly(P.subs(L_, lam), s_)
            n = Ps.count_roots(S.Rational((1 - hi).numerator, (1 - hi).denominator),
                               S.Rational((1 - lo).numerator, (1 - lo).denominator))
            assert n == 1, (p, t, n)
            out.append({"protocol": p, "F": t, "lo": f"{float(lo):.17g}", "hi": f"{float(hi):.17g}",
                        "lo30": str(S.Float(S.Rational(lo.numerator, lo.denominator), 35)),
                        "sturm_roots": int(n)})
            print(p, t, out[-1]["lo30"])
    with open(os.path.join(HERE, "reference_qtarget.json"), "w") as f:
        json.dump({"source": "generate_reference_qtarget.py (arb 256-bit, exact rational bisection, "
                             "Sturm count on P_P(s, 4F-1))", "rows": out}, f, indent=1)


if __name__ == "__main__":
    main()
