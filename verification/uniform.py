# Certified constants of the uniform theorem (ball arithmetic).
from lbound import *
import math
def uniform(B, q0, y0, ell, u):
    B=arb(B); q0=arb(q0); y0=arb(y0); ell=arb(ell); u=arb(u)
    k=1/(1-B)
    logy0=y0.log(); logx0=logy0+ell; x0=logx0.exp()
    S0=S_bound(q0,u,B)
    c0=(((2+u)/2).digamma()-PI.log())/2+Z(u)
    assert c0>0
    E0=((((B-1)*logx0).exp()/logx0+((B-1)*logy0).exp()/logy0)*S0+(1/(x0*logx0)+1/(y0*logy0))*PI**2/8)/ell
    tau0=3/y0.sqrt()
    kappa0=logy0+ell/2-k*q0.log().log()
    Cm=(Gsum()-B1-1/(2*logy0**2)-tau0-E0).exp()
    Cp=(EULER+1/(2*logy0**2)+tau0+E0).exp()
    return dict(k=k,kappa0=kappa0,Cm=Cm,Cp=Cp,E0=E0,S0=S0,c0=c0,G=Gsum(),tau0=tau0)
if __name__=='__main__':
    q0=2**40
    for name,B,y0,ell,u in [('7/8','0.875','5e6',14,'0.1'),('11/12',arb(11)/12,'1.8e9',20,'0.08')]:
        r=uniform(B,q0,arb(y0),ell,u)
        print(name,{k:(float(v.lower()),float(v.upper())) for k,v in r.items()})
        lq=math.log(math.log(q0))
        denom=float(r['k'].mid())*lq+float(r['kappa0'].upper())
        Lb=float(r['Cm'].lower())/denom
        print('   L(1)>=',Lb,' at q0; h>=',math.sqrt(q0)*Lb/math.pi,' Lmax<=',float(r['Cp'].upper())*denom)
