from opt import opt
from lbound import arb
import math
for name,B in [('11/12',arb(11)/12),('7/8',arb(7)/8),('Liu',arb('0.874957019420098946'))]:
  b=opt(2**40,B,[x/4 for x in range(12,60)],list(range(4,41,2)),['0.06','0.08','0.1','0.12','0.15'])
  print(name,b[:4],'L>=',math.exp(b[0]),'h>=',math.sqrt(2**40)*math.exp(b[0])/math.pi, 'Lmax<=',math.exp(b[5]))
