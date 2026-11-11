import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from uniform import uniform; from lbound import arb
q0=2**40
P={}
for name,B,y0,ell,u in [('7/8','0.875','5e6',14,'0.1'),('11/12',arb(11)/12,'1.8e9',20,'0.08')]:
    r=uniform(B,q0,arb(y0),ell,u); P[name]=(float(r['k'].mid()),float(r['kappa0'].upper()),float(r['Cm'].lower()))
t=np.linspace(np.log10(q0),20,400); D=10**t
fig,ax=plt.subplots(figsize=(6.4,4.0))
for name,c,ls in [('7/8','#1f4e9c','-'),('11/12','#b8461b','--')]:
    k,ka,C=P[name]; h=np.sqrt(D)*C/(np.pi*(k*np.log(np.log(D))+ka))
    ax.semilogy(t,h,ls,color=c,lw=2,label=r'zero-free for Re $s>%s$'%name)
ax.axhline(1e4,color='k',lw=0.8,ls=':'); ax.text(16.6,1.12e4,r'$h=10^4$',fontsize=9)
ax.axvspan(9,np.log10(q0),color='#9bbf85',alpha=0.35)
ax.text(9.25,3e5,'tabulated\nunconditionally\n$|d|<2^{40}$',fontsize=8.5)
ax.axvline(np.log10(q0),color='#3a6b2a',lw=1)
ax.set_xlim(9,20); ax.set_ylim(3e3,1e9)
ax.set_xlabel(r'$\log_{10}|d|$'); ax.set_ylabel(r'certified lower bound for $h(d)$')
ax.legend(loc='upper left',bbox_to_anchor=(0.36,0.98),fontsize=8.5,frameon=False)
fig.tight_layout()
for e in ['png','pdf']: fig.savefig('../paper/figures/fig_hbound.'+e,dpi=200)
print(P)
