import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpmath import zeta, mpc
OUT='/tmp/claude-0/-home-user/71e36f3f-5336-59e8-93a9-66a59a73ec3f/scratchpad/paper/figures/'
plt.rcParams.update({'font.family':'serif','font.size':10})
# Figure 1: log|zeta| over the strip with the zero-free half-plane
sig=np.linspace(0.05,1.3,90); t=np.linspace(2,42,240)
S,T=np.meshgrid(sig,t)
Zv=np.vectorize(lambda s,tt: float(abs(zeta(mpc(s,tt)))))(S,T)
LZ=np.log(np.maximum(Zv,1e-3))
fig=plt.figure(figsize=(7.2,5.0)); ax=fig.add_subplot(111,projection='3d')
ax.plot_surface(S,T,LZ,cmap='viridis',rstride=2,cstride=2,linewidth=0,alpha=0.85)
yy,zz=np.meshgrid(np.linspace(2,42,2),np.linspace(-6,3,2))
ax.plot_surface(0.875*np.ones_like(yy),yy,zz,color='crimson',alpha=0.25)
ax.plot_surface(1.0*np.ones_like(yy),yy,zz,color='gray',alpha=0.15)
g=[14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719]
ax.scatter([0.5]*len(g),g,[-6]*len(g),color='k',s=12)
ax.set_xlabel(r'$\sigma$'); ax.set_ylabel(r'$t$'); ax.set_zlabel(r'$\log|\zeta(\sigma+it)|$')
ax.set_zlim(-6,3); ax.view_init(28,-60)
ax.text(0.9,40,3.2,r'$\mathrm{Re}\,s=7/8$',color='crimson')
plt.tight_layout(); plt.savefig(OUT+'fig_strip.png',dpi=200); plt.savefig(OUT+'fig_strip.pdf'); plt.close()
# Figure 2: certified lower bound surface at q=2^40
from lbound import bounds, arb
ly=np.linspace(4.5,9.5,36); el=np.linspace(4,26,34)
LY,EL=np.meshgrid(ly,el)
f=lambda a,b: float(bounds(2**40,arb(7)/8,arb(a)*arb(10).log(),arb(b),arb('0.1'))[0].mid())
V=np.vectorize(f)(LY,EL)
fig=plt.figure(figsize=(7.2,5.0)); ax=fig.add_subplot(111,projection='3d')
ax.plot_surface(LY,EL,np.exp(V),cmap='plasma',linewidth=0,alpha=0.9)
i=np.unravel_index(np.argmax(V),V.shape)
ax.scatter([LY[i]],[EL[i]],[np.exp(V[i])],color='k',s=30)
ax.set_xlabel(r'$\log_{10} y$'); ax.set_ylabel(r'$\ell=\log(x/y)$'); ax.set_zlabel(r'certified lower bound for $L(1,\chi)$')
ax.view_init(25,-130)
plt.tight_layout(); plt.savefig(OUT+'fig_surface.png',dpi=200); plt.savefig(OUT+'fig_surface.pdf'); plt.close()
print('max at',LY[i],EL[i],np.exp(V[i]))
