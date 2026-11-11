import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fn, timestr = sys.argv[1], sys.argv[2]
L=open(fn).read().split('\n'); hdr=L[0]
kv=dict(t.split('=') for t in hdr[2:].replace(' at d=',' d=').split())
rows=[tuple(map(int,l.split())) for l in L[1:] if l.strip()]
H=len(rows); F={h:f for h,f,m in rows}; M={h:m for h,f,m in rows}
ntot=sum(F.values()); dmax=max(M.values()); hd=[h for h in M if M[h]==dmax][0]
g=lambda n: f'{n:,}'.replace(',','\\,')
P='../paper/'
open(P+'data_macros.tex','w').write(
 f"\\newcommand{{\\NBOUND}}{{{g(int(kv['N']))}}}\n\\newcommand{{\\NTOT}}{{{g(ntot)}}}\n\\newcommand{{\\DMAX}}{{{g(dmax)}}}\n"
 f"\\newcommand{{\\HDMAX}}{{{hd}}}\n\\newcommand{{\\NFUND}}{{{g(int(kv['fundamental']))}}}\n\\newcommand{{\\MINL}}{{{kv['minL']}}}\n"
 f"\\newcommand{{\\MINLD}}{{-{g(-int(kv['d']))}}}\n\\newcommand{{\\TIMECOUNT}}{{{timestr}}}\n")
sel=list(range(1,41))+[h for h in range(50,H+1,25)]
sel=[h for h in sel if h<=H]
# three-column-pair table
per=(len(sel)+2)//3; cols=[sel[i*per:(i+1)*per] for i in range(3)]
out=['\\begin{center}\\small','\\begin{tabular}{@{}rrr@{\\qquad}rrr@{\\qquad}rrr@{}}','\\toprule',
     ' & '.join(['$h$ & $F(h)$ & max $|d|$']*3)+r'\\','\\midrule']
for i in range(per):
    cells=[]
    for c in cols:
        if i<len(c): h=c[i]; cells+= [str(h),g(F[h]),g(M[h])]
        else: cells+=['','','']
    out.append(' & '.join(cells)+r'\\')
out+=['\\bottomrule','\\end{tabular}','\\end{center}']
open(P+'data_table.tex','w').write('\n'.join(out)+'\n')
hs=np.array(sorted(F)); f=np.array([F[h] for h in hs]); m=np.array([M[h] for h in hs])
fig,(a,b)=plt.subplots(1,2,figsize=(8,3.3))
v2=np.array([(int(h) & -int(h)).bit_length()-1 for h in hs])
sc=a.scatter(hs,f,c=np.minimum(v2,6),s=5,cmap='viridis'); a.set_xlabel('$h$'); a.set_ylabel('$F(h)$')
cb=fig.colorbar(sc,ax=a,pad=0.01); cb.set_label('2-adic valuation of $h$',fontsize=8)
a.set_title('(a) number of fields',fontsize=9)
b.loglog(hs,m,'.',ms=3,color='#1f4e9c',label='largest $|d|$')
hh=np.linspace(1,H,400); b.loglog(hh,(np.pi*hh/0.1706)**2,'-',color='#b8461b',lw=1.2,label=r'$(\pi h/0.1706)^2$')
b.set_xlabel('$h$'); b.set_ylabel('$|d|$'); b.legend(fontsize=8,frameon=False); b.set_title('(b) largest discriminant',fontsize=9)
fig.tight_layout()
for e in ['png','pdf']: fig.savefig(P+'figures/fig_counts.'+e,dpi=200)
print(ntot,dmax,hd,kv)
