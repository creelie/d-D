from lbound import *
import math
def opt(q, B, grid_y, grid_l, grid_u):
    best=None
    for ly in grid_y:
        for ell in grid_l:
            for u in grid_u:
                lo,hi,d = bounds(q,B,arb(ly)*arb(10).log(),ell,arb(u))
                v=float(lo.lower())
                if best is None or v>best[0]: best=(v,ly,ell,u,d,float(hi.upper()))
    return best
if __name__=="__main__":
  B=arb(7)/8
  b=opt(2**40,B,[x/4 for x in range(12,48)],list(range(4,41,2)),['0.06','0.08','0.1','0.12','0.15'])
  print(b[:4], [float(b[4][k].mid()) for k in ('S','E','W')], 'L>=',math.exp(b[0]),'h>=',math.sqrt(2**40)*math.exp(b[0])/math.pi)
  ly,ell=b[1],b[2]
  b=opt(2**40,B,[ly+x/20 for x in range(-5,6)],[ell+x/2 for x in range(-4,5)],[str(0.06+0.01*i) for i in range(10)])
  print(b[:4], 'L>=',math.exp(b[0]),'h>=',math.sqrt(2**40)*math.exp(b[0])/math.pi)
