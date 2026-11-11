/* Exact class numbers h(d) of all fundamental discriminants -N <= d < 0 by
   counting reduced forms (a,b,c): b^2-4ac=d, |b|<=a<=c, b>=0 if |b|=a or a=c.
   Output: for each h<=H: F(h)=#fundamental d with h(d)=h, and max|d|; plus
   the minimum of h*pi/sqrt|d| (= L(1,chi_d) for d<-4) over the range.
   Usage: classcount N H outfile                                   */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#define BLK (1<<22)
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: N H out\n");return 1;}
  int64_t N=atoll(argv[1]); int H=atoi(argv[2]);
  uint16_t *cnt=calloc(N+1,sizeof(uint16_t));
  uint8_t *sqf=malloc(N/4+2); /* squarefree flag for m<=N/4+1 ... we need up to N */
  free(sqf);
  /* squarefree sieve up to N (bit array) */
  uint8_t *nsq=calloc(N/8+2,1);
  for(int64_t p=2;p*p<=N;p++){int64_t s=p*p; for(int64_t m=s;m<=N;m+=s) nsq[m>>3]|=1<<(m&7);}
  int64_t nblk=(N+BLK)/BLK;
  #pragma omp parallel for schedule(dynamic,1)
  for(int64_t bi=0;bi<nblk;bi++){
    int64_t D0=bi*(int64_t)BLK, D1=D0+BLK; if(D1>N+1) D1=N+1; if(D0<3) D0=3;
    if(D0>=D1) continue;
    int64_t A=(int64_t)sqrt((double)(D1-1)/3.0)+1;
    for(int64_t a=1;a<=A;a++){
      int64_t fa=4*a;
      for(int64_t b=-a+1;b<=a;b++){
        int64_t b2=b*b;
        int64_t cmin=(D0+b2+fa-1)/fa; int64_t c0=a; if(b<0) c0=a+1;
        if(cmin<c0) cmin=c0;
        int64_t cmax=(D1-1+b2)/fa;
        if(cmin>cmax) continue;
        int64_t D=fa*cmin-b2;
        for(int64_t c=cmin;c<=cmax;c++,D+=fa) cnt[D]++;
      }
    }
  }
  /* collect */
  int64_t *F=calloc(H+1,sizeof(int64_t)), *mx=calloc(H+1,sizeof(int64_t));
  double minL=1e9; int64_t argmin=0; int64_t nfund=0; uint64_t sumh=0; int maxh=0;
  for(int64_t D=3;D<=N;D++){
    int fund=0;
    if((D&3)==3){ if(!(nsq[D>>3]&(1<<(D&7)))) fund=1; }
    else if((D&3)==0){ int64_t m=D>>2; if(((m&3)==1||(m&3)==2) && !(nsq[m>>3]&(1<<(m&7)))) fund=1; }
    if(!fund) continue;
    nfund++; int h=cnt[D]; sumh+=h; if(h>maxh) maxh=h;
    if(h==0){fprintf(stderr,"zero count at %lld\n",(long long)D);return 2;}
    if(D>4){ double L=M_PI*h/sqrt((double)D); if(L<minL){minL=L;argmin=D;} }
    if(h<=H){F[h]++; if(D>mx[h]) mx[h]=D;}
  }
  FILE*o=fopen(argv[3],"w");
  fprintf(o,"# N=%lld H=%d fundamental=%lld sum_h=%llu max_h=%d minL=%.6f at d=-%lld\n",(long long)N,H,(long long)nfund,(unsigned long long)sumh,maxh,minL,(long long)argmin);
  for(int h=1;h<=H;h++) fprintf(o,"%d %lld %lld\n",h,(long long)F[h],(long long)mx[h]);
  fclose(o); return 0;
}
