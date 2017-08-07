#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "point.h"

void bsgs(mpz_t p,PROJECTIVE_POINT P,const unsigned long int B1,const unsigned long int B2,const int window_size,const mpz_t N,const FILE *fp,const mpz_t D,const mpz_t ma, const mpz_t mb){

  unsigned long int f=210;
  unsigned long int f2 = f/2;
  unsigned long int f4 = f/4;
  unsigned long int i,s,v,v2,u;
  mpz_t tB1,ts,inv,d,product,a1,gcd,ti,tf;

  mpz_inits(tB1,ts,inv,d,product,a1,gcd,ti,tf,NULL);
  mpz_set_ui(tB1,B1);
  mpz_nextprime(ts,tB1);
  s=mpz_get_ui(ts);
  mpz_set_ui(d,1);
  mpz_set_ui(tf,f);
 
  EXTENDED_POINT eP;
  MONTGOMERY_POINT mP;
  MONTGOMERY_POINT baby_step[f4];
  MONTGOMERY_POINT double_mP;
  MONTGOMERY_POINT add_mP;
  MONTGOMERY_POINT giant_step;
  MONTGOMERY_POINT Giant;
  extended_point_init(eP);
  montgomery_point_init(mP);
  fprintf(stdout,"X=%ld Y=%ld Z=%ld\n",mpz_get_ui(mP->X),mpz_get_ui(mP->Y),mpz_get_ui(mP->Z));
  protoext(eP,P,N);
  fprintf(stdout,"X=%ld Y=%ld Z=%ld\n",mpz_get_ui(eP->X),mpz_get_ui(eP->Y),mpz_get_ui(eP->Z));
  exttomon(mP,eP,N);
  fprintf(stdout,"X=%ld Y=%ld Z=%ld\n",mpz_get_ui(mP->X),mpz_get_ui(mP->Y),mpz_get_ui(mP->Z));
  montgomery_point_init(giant_step);
  montgomery_point_init(Giant);
  fprintf(stdout,"1\n");
  montgomery_point_init(double_mP);
  fprintf(stdout,"2\n");
  montgomery_double(double_mP,mP,ma,N);///
  for(i=0;i<f4;i++){
    montgomery_point_init(baby_step[i]);
  }
  fprintf(stdout,"3\n");
  montgomery_point_init(add_mP);
  montgomery_point_set(add_mP,mP);
  for(i=0;i<f4;i++){
    mpz_set_ui(ti,2*i+1);
    mpz_gcd(gcd,ti,tf);
    if(mpz_cmp_ui(gcd,1)==0){
      montgomery_point_set(baby_step[i],add_mP);
    }
    montgomery_add(add_mP,add_mP,double_mP,N);
  }
  fprintf(stdout,"4\n");
  v2=(s-f/2)/f;
  /*動かない
  mscalar(giant_step,mP,f,window_size,ma,mb,N);
  mscalar(Giant,giant_step,v2,window_size,ma,mb,N);
  */
  fprintf(stdout,"5\n");
  while(s<=B2){
    v=(s-f/2)/f;
    while(v!=v2){
      v2++;
      montgomery_add(Giant,Giant,giant_step,N);
    }
    fprintf(stdout,"6\n");
    u=abs((s-f/2)%f-f/2);
    mpz_mul_mod(product,Giant->X,baby_step[(u-1)/2]->Z,N);
    mpz_mul_mod(a1,Giant->Z,baby_step[(u-1)/2]->X,N);
    mpz_sub(product,product,a1);
    mpz_mul_mod(d,d,product,N);
    mpz_nextprime(ts,ts);
    s=mpz_get_ui(ts);
    fprintf(stdout,"7\n");
  }
  mpz_gcd(p,d,N);
  mpz_clears(tB1,ts,inv,d,product,a1,gcd,ti,tf,NULL);
  extended_point_clear(eP);
  montgomery_point_clear(mP);
  for(i=0;i<f4;i++){
    montgomery_point_clear(baby_step[i]);
  }
  montgomery_point_clear(giant_step);
  montgomery_point_clear(Giant);
  montgomery_point_clear(double_mP);
  montgomery_point_clear(add_mP);
  fprintf(stdout,"8\n");
}
