import numpy as np, sys
def h(D):
    n=0; amax=int((D/3)**0.5)+1
    for a in range(1,amax+1):
        b=np.arange(-a+1,a+1,dtype=np.int64)
        num=b*b+D; ok=(num%(4*a)==0); c=num//(4*a)
        ok&=(c>=a); ok&=~((c==a)&(b<0))
        n+=int(ok.sum())
    return n
for D in map(int,sys.argv[1:]): print(D,h(D))
