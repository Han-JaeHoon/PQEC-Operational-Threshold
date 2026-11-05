"""
Independent check of the eps_c brackets of src/epsc.ts (Theorem 2.1 of the calculation notes, float64
outward rounding, order-preserving corners) by a different certificate and a different arithmetic:
Proposition 4 of "The separating noise eps_c(q)" with 256-bit arb balls, mean-value enclosures and
exactly verified sub/super-solutions u, l (epsc_research/epsc_core.py and jac_gen.py, the research code).

For each q (off the research grid, near q_m, near q_bi, near q_SN) the exact float64 values of q, lo and
hi from the simulator are used: lo must be certified into the P_B basin and hi into the Z basin.

Run from interactive/type3-dynamics (needs node_modules, and python-flint, mpmath):
    python verification/generate_reference_epsc.py
"""
import json, os, subprocess, sys
from fractions import Fraction as Fr

import mpmath as mp
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "epsc_research"))
from epsc_core import T, fixed_points_arb  # noqa: E402
from jac_gen import JAC  # noqa: E402

ctx.prec = 256
mp.mp.dps = 50
KIND = {"q4": "4Q", "type4": "T4", "type3": "T3"}
LIST = [
    ["q4", 0.1762], ["q4", 0.17712345], ["q4", 0.1785], ["q4", 0.1801], ["q4", 0.1806], ["q4", 0.17584],
    ["type4", 0.17955], ["type4", 0.1808], ["type4", 0.180828], ["type4", 0.183], ["type4", 0.1871], ["type4", 0.1892],
    ["type3", 0.1185], ["type3", 0.1187], ["type3", 0.11873], ["type3", 0.1205], ["type3", 0.1277], ["type3", 0.1303],
]


def to_arb(x):
    return arb(x.numerator)/arb(x.denominator)


def setup(kind, q):
    s = 1 - q
    sa = to_arb(s)
    (tS, mS), _, mZ = fixed_points_arb(kind, sa)
    sm = mp.mpf(s.numerator)/s.denominator
    tf, mf = mp.mpf(tS.mid().str(40, radius=False)), mp.mpf(mS.mid().str(40, radius=False))
    h = mp.mpf("1e-20")
    f = lambda a, b: T(kind, sm, a, b)
    A = mp.matrix([[(f(tf + h, mf)[0] - f(tf - h, mf)[0])/(2*h), (f(tf, mf + h)[0] - f(tf, mf - h)[0])/(2*h)],
                   [(f(tf + h, mf)[1] - f(tf - h, mf)[1])/(2*h), (f(tf, mf + h)[1] - f(tf, mf - h)[1])/(2*h)]])
    ev, V = mp.eig(A)
    k = 0 if mp.re(ev[0]) > mp.re(ev[1]) else 1
    w = [mp.re(V[0, k]), mp.re(V[1, k])]
    if w[0] < 0:
        w = [-w[0], -w[1]]
    assert w[0] > 0 and w[1] < 0
    for r in [Fr(1, 10**3), Fr(1, 10**4), Fr(1, 10**5), Fr(1, 10**6)]:
        den = 10**14
        u = (Fr(int(mp.nint((tf + mp.mpf(r.numerator)/r.denominator*w[0])*den)), den),
             Fr(int(mp.nint((mf + mp.mpf(r.numerator)/r.denominator*w[1])*den)), den))
        l = (Fr(int(mp.nint((tf - mp.mpf(r.numerator)/r.denominator*w[0])*den)), den),
             Fr(int(mp.nint((mf - mp.mpf(r.numerator)/r.denominator*w[1])*den)), den))
        Tu, Tl = T(kind, s, *u), T(kind, s, *l)
        if (Tu[0] >= u[0] and Tu[1] <= u[1] and to_arb(u[0]) > tS and u[0] <= 1 and
                Tl[0] <= l[0] and Tl[1] >= l[1] and to_arb(l[1]) > mS and to_arb(l[1]) < mZ and l[0] >= 0):
            return sa, u, l
    raise RuntimeError(f"no u, l for {kind} {q}")


def side(kind, sa, u, l, eps, nmax=400000):
    th = 1 - eps
    t, m = to_arb(2*th/(1 + th)), to_arb((1 - th)/(2*(1 + th)))
    U = [to_arb(x) for x in u]
    L = [to_arb(x) for x in l]
    for n in range(nmax):
        if t >= U[0] and m <= U[1]:
            return +1, n
        if t <= L[0] and m >= L[1]:
            return -1, n
        tm, mm = arb(t.mid()), arb(m.mid())
        rt, rm = arb(0, t.rad()), arb(0, m.rad())
        T0 = T(kind, sa, tm, mm)
        a, b, c, d = JAC[kind](sa, t, m)
        t, m = T0[0] + a*rt + b*rm, T0[1] + c*rt + d*rm
        if t.rad() > 1e-6 or m.rad() > 1e-6:
            return 0, n
    return 0, nmax


here = os.path.dirname(HERE)
bundle = os.path.join(os.environ.get("TMPDIR", "/tmp"), "dump_epsc_brackets.mjs")
subprocess.check_call(["node_modules/.bin/esbuild", "verification/dump_epsc_brackets.ts", "--bundle",
                       "--platform=node", "--format=esm", f"--outfile={bundle}", "--log-level=warning"], cwd=here)
out = subprocess.run(["node", bundle], input=json.dumps(LIST), capture_output=True, text=True, check=True, cwd=here)
rows, ok = [], True
for r in json.loads(out.stdout):
    kind = KIND[r["protocol"]]
    sa, u, l = setup(kind, Fr(r["q"]))
    a = side(kind, sa, u, l, Fr(r["lo"]))
    b = side(kind, sa, u, l, Fr(r["hi"]))
    good = a[0] == +1 and b[0] == -1
    ok &= good
    rows.append(dict(protocol=r["protocol"], q=r["q"], lo=r["lo"], hi=r["hi"], width=r["hi"] - r["lo"],
                     arb_side_lo=a[0], arb_rounds_lo=a[1], arb_side_hi=b[0], arb_rounds_hi=b[1]))
    print(f"{r['protocol']:6s} q={r['q']:<11} [{r['lo']:.14f}, {r['hi']:.14f}] w={r['hi']-r['lo']:.1e}  "
          f"arb: lo->{a[0]:+d} (n={a[1]}) hi->{b[0]:+d} (n={b[1]})  {'OK' if good else 'FAIL'}", flush=True)
json.dump(dict(note="lo, hi from src/epsc.ts; sides certified independently with arb (+1: P_B, -1: Z)", rows=rows),
          open(os.path.join(HERE, "reference_epsc.json"), "w"), indent=1)
print("ALL OK" if ok else "FAILED")
sys.exit(0 if ok else 1)
