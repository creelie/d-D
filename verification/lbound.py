# Rigorous (ball arithmetic) explicit Littlewood bounds from a zero-free half-plane Re s > B.
from flint import arb, arb_series, ctx
from sympy import primerange
ctx.prec = 128
PI = arb.pi(); EULER = arb.const_euler()
B1 = arb('0.26149721284764278375542683860869585905') + arb('0 +/- 1e-35')  # Meissel-Mertens
P_SUM = 10**6
_G = None
def Gsum():
    global _G
    if _G is None:
        s = arb(0)
        for p in primerange(2, P_SUM + 1):
            ip = arb(1)/p
            s += ip - (1 + ip).log()
        _G = s
    return _G
def Z(u):
    z = arb_series([1 + u, 1], prec=2).zeta()
    return -z[1]/z[0]
def S_bound(q, u, B):
    return (1 + u - B)/(1 - B)**2 * ((q/PI).log()/2 + ((2 + u)/2).digamma()/2 + Z(u))
def F(v):
    return v*v.log() - v + B1*v - 1/(2*v)
def bounds(q, B, logy, ell, u):
    """returns (lower, upper) for log|L(1,chi)|, chi primitive mod q, zeros in Re s <= B."""
    q = arb(q); B = arb(B); u = arb(u); logy = arb(logy); ell = arb(ell)
    logx = logy + ell
    y = logy.exp(); x = logx.exp()
    S = S_bound(q, u, B)
    E = ((((B - 1)*logx).exp()/logx + ((B - 1)*logy).exp()/logy)*S
         + (1/(x*logx) + 1/(y*logy))*PI**2/8)/ell
    W = (F(logx) - F(logy))/ell
    tau = 3/y.sqrt()
    lo = -W + Gsum() - tau - E
    hi = W + (EULER - B1) + tau + E
    return lo, hi, dict(S=S, E=E, W=W)
if __name__ == '__main__':
    import itertools, math
    B = arb(7)/8
    q0 = 2**40
    best = None
    for logy10 in [x/4 for x in range(28, 61)]:
        for ell in [2,3,4,5,6,7,8,9,10,12]:
            for u in ['0.06','0.08','0.1','0.12','0.15']:
                lo, hi, d = bounds(q0, B, arb(logy10)*arb(10).log(), ell, arb(u))
                v = float(lo.lower())
                if best is None or v > best[0]: best = (v, logy10, ell, u, d)
    v = best[0]
    print('q0=2^40 best logL >=', v, 'L>=', math.exp(v), 'params log10 y, ell, u =', best[1:4])
    print('  S,E,W =', [float(best[4][k].mid()) for k in ('S','E','W')])
    print('  h >=', math.sqrt(q0)*math.exp(v)/math.pi)
