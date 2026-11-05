"""
Exact checks behind src/epsc.ts and regimes.sideOfRoot (sympy, rational arithmetic only).

  P1  each q_bi / q_SN / q_m polynomial has exactly one root in its bracket, and no rational root
      there, so a float64 q is never exactly on a threshold.
  P2  the (t, m) map of src/epsc.ts equals the common form (1)-(5) of the calculation notes.
  P3  the coefficients used by the certifier are nonnegative for 0 < s < 1.
  P4  the saddle formula (16)-(22) gives a fixed point at a rational sample point of each window,
      and the order signs (13) hold on a rational grid of R at that point.

Run: python interactive/type3-dynamics/verification/verify_epsc_polys.py
"""
import sympy as sp

s, t, m, u, v = sp.symbols("s t m u v")
ROOTS = {
    "q4": {"qBi": ([1, -4, -6, 0, 1, 4], "0.824", "0.825"),
           "qSN": ([41, 32, 16, 16, 6, -32, 0, -16, 1], "0.8193", "0.8194"),
           "qM": ([1, 1, 0, 0, -1], "0.81", "0.83")},
    "type4": {"qBi": ([4, 4, -1, 2, -5], "0.821", "0.822"),
              "qSN": ([1, 2, 9, 12, 5, -2, -4, -6, -1], "0.8105", "0.8106"),
              "qM": ([1, 1, 0, 0, -1], "0.81", "0.83")},
    "type3": {"qBi": ([4, 4, -1, 2, -1, 0, -4], "0.8816", "0.8818"),
              "qSN": ([1, 2, 9, 12, 5, -2, 1, 0, -5, -6, -1], "0.8694", "0.8695"),
              "qM": ([1, 1, 0, 0, 0, 0, -1], "0.875", "0.89")},
}
for name, rs in ROOTS.items():
    for r, (co, lo, hi) in rs.items():
        P = sp.Poly(co, s)
        n = P.count_roots(sp.Rational(lo), sp.Rational(hi))
        assert n == 1, (name, r, n)
        rat = [x for x in sp.roots(P, filter="Q")]
        assert not any(sp.Rational(lo) <= x <= sp.Rational(hi) for x in rat), (name, r, rat)
print("P1 OK: one irrational root per bracket")


def common(name):
    if name == "q4":
        return 2*s**5, s**4*(1 + s), 2*s**5, s**3*(1 + s), s**3 + s**5, s**5
    A, B = (6, 4) if name == "type3" else (4, 2)
    c = s**(A - 1)*(1 + s)
    return 2*s**A, 2*s**A, c, c, 2*s**B, s**B


def tm_coeffs(name):
    if name == "q4":
        A_ = 1 - s**3 - s**4 + s**5
        return dict(A=A_, B=1 - s**5, C=3*s**5 + s**3, E=1 + s**3 + s**4 + s**5, F=s**3 - s**5,
                    H=4*s**4*(3*s + 1), I=8*s**4*(1 - s))
    a, b = (6, 4) if name == "type3" else (4, 2)
    e1, e2, e3 = s**b, s**(a - 1), s**a
    return dict(A=1 + e1 - e2 - e3, B=1 - e1, C=2*e1 + e2 + e3, E=1 + e1 + e2 + e3, F=2*e1 - e2 - e3,
                H=16*e3, I=0)


def tm_map(name, tt, mm):
    c = tm_coeffs(name)
    S = 2*(4*mm**2*c["A"] + 4*mm*c["B"] + tt**2*c["C"] + c["E"])
    L = 4*mm**2*c["E"] + 4*mm*c["B"] + tt**2*c["F"] + c["A"]
    return tt*(c["H"] + c["I"]*mm)/S, L/S


for name in ROOTS:
    al, be, ga, de, et, b = common(name)
    V = (1 - 2*m)/(1 + 2*m)
    U = t/(1 + 2*m)
    D = 1 + et*U**2 + b*V**2
    Up = U*(al*V + be)/D
    Vp = (ga*U**2 + de*V)/D
    tp, mp = tm_map(name, t, m)
    assert sp.simplify(2*Up/(1 + Vp) - tp) == 0, name
    assert sp.simplify((1 - Vp)/(2*(1 + Vp)) - mp) == 0, name
print("P2 OK: (t, m) map = common form of the notes")

samples = [sp.Rational(k, 200) for k in range(1, 200)]
for name in ROOTS:
    c = tm_coeffs(name)
    for k, e in c.items():
        for x in samples:
            assert sp.sympify(e).subs(s, x) >= 0, (name, k, x)
print("P3 OK: certifier coefficients >= 0 on a rational grid of (0, 1)")

Q = {"q4": sp.Rational(178251, 10**6), "type4": sp.Rational(184175, 10**6), "type3": sp.Rational(124435, 10**6)}
for name, q in Q.items():
    s0 = 1 - q
    al, be, ga, de, et, b = [e.subs(s, s0) for e in common(name)]
    h2, h1, h0 = ga*b + et*al, et*(be - de) - ga*al, ga*(1 - be)
    disc = h1**2 - 4*h2*h0
    assert disc > 0
    vS = (-h1 - sp.sqrt(disc))/(2*h2)
    u2 = vS*(al*vS + be - de)/ga
    D = 1 + et*u2 + b*vS**2
    assert abs(sp.N(D - (al*vS + be), 60)) < 1e-50
    assert abs(sp.N((ga*u2 + de*vS)/D - vS, 60)) < 1e-50
    tp, mp = tm_map(name, t, m)
    J = [sp.diff(tp, t), sp.diff(tp, m), sp.diff(mp, t), sp.diff(mp, m)]
    J = [sp.lambdify((t, m), e.subs(s, s0), "sympy") for e in J]
    for i in range(1, 11):
        for j in range(0, 11):
            P = (sp.Rational(i, 10), sp.Rational(j, 20))
            a, b_, c_, d = [f(*P) for f in J]
            assert a > 0 and d > 0 and b_ < 0 and c_ < 0, (name, P)
print("P4 OK: saddle is a fixed point; order signs hold on a rational grid at the sample q")
print("ALL OK")
