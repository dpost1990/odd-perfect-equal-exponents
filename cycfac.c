// cycfac: find all q <= X with q | Phi_l(Q), for prime l and integer Q (decimal, any size).
// Every prime factor q != l of Phi_l(Q) satisfies q == 1 (mod 2l) (q odd) and Q^l == 1, Q != 1 (mod q).
// We test every q == 1 (mod 2l) coprime to 3*5*7 (composite q passing the test are reported too;
// the caller only keeps primes).  Usage: cycfac l Q lo hi
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>
typedef unsigned __int128 u128;
static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m){ return (uint64_t)(((u128)a*b)%m); }
static inline uint64_t powmod(uint64_t a, uint64_t e, uint64_t m){ uint64_t r=1; a%=m; while(e){ if(e&1) r=mulmod(r,a,m); a=mulmod(a,a,m); e>>=1;} return r; }
int main(int argc,char**argv){
  if(argc<5){fprintf(stderr,"usage: cycfac l Q lo hi\n");return 1;}
  uint64_t l=strtoull(argv[1],0,10); mpz_t Q; mpz_init_set_str(Q,argv[2],10);
  uint64_t lo=strtoull(argv[3],0,10), hi=strtoull(argv[4],0,10);
  uint64_t step=2*l; uint64_t q = (lo/step)*step+1; if(q<lo) q+=step; if(q==1) q+=step;
  int small = mpz_sizeinbase(Q,2)<=64; uint64_t Qs = small? mpz_get_ui(Q):0;
  // for big Q: precompute Q = sum limbs; use mpz_fdiv_ui
  for(; q<=hi; q+=step){
    if(q%3==0||q%5==0||q%7==0||q%11==0||q%13==0) { if(q>13) continue; }
    uint64_t r = small? Qs%q : mpz_fdiv_ui(Q,q);
    if(r<=1) continue;
    if(powmod(r,l,q)==1) printf("%llu\n",(unsigned long long)q);
  }
  return 0;
}
