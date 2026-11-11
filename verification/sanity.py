# numerical check of Corollary 2.3 and Theorem 5.1 for real characters chi_d, d=-q
import numpy as np, math, sys
from sympy import jacobi_symbol, primerange
from lbound import bounds, arb
def kron(d,n):  # Kronecker symbol (d/n), d<0 fundamental, n prime power base p
    if n==2: return 0 if d%2==0 else (1 if d%8 in (1,7) else -1)
    return jacobi_symbol(d% n, n)
def check(D, h, logy, ell):
    d=-D; y=math.exp(logy); x=y*math.exp(ell)
    M=0.0
    for p in primerange(2,int(x)+1):
        c=kron(d,p); 
        if c==0: continue
        pk=p; k=1
        while pk<=x:
            w=1 if pk<=y else math.log(x/pk)/ell
            M+= (c**k)*w/(k*pk); pk*=p; k+=1
    logL=math.log(math.pi*h/math.sqrt(D))
    lo,hi,dd=bounds(arb(D),arb('0.875'),arb(logy),arb(ell),arb('0.1'))
    print(D,'logL=%.6f M=%.6f diff=%.2e  E(7/8)=%.4f  bounds [%.4f, %.4f]'%(logL,M,logL-M,float(dd['E'].mid()),float(lo.mid()),float(hi.mid())))
for D,h in [(991027,63),(59249107,512),(227932027,996),(163,1)]:
    check(D,h,math.log(3e5),3.0)
