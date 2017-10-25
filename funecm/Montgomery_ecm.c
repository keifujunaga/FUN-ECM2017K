#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "gmp.h"
#include "point.h"
//aaaaaaaaaaaaaaa

//モンゴメリー曲線の公式（射影座標）
//b*y^2=x^3+ax^2+x

void Montgomery_ecm(mpz_t f, const mpz_t N, const mpz_t X, const mpz_t Y, mpz_t d, const mpz_t montgomery_a , const unsigned long int B1, const unsigned long int B2, FILE *fp, const int window_size)
{
  PROJECTIVE_POINT P;
  int e;
  int i;
  mpz_t tmp;
  mpz_t tmp2;
  mpz_t inv;
  mpz_t montgomery_b;
	
  mpz_init(tmp);
  mpz_init(tmp2);
  mpz_init(inv);
  mpz_init(montgomery_b);

  projective_point_init(P);
		
  /* set P */
  if (X == NULL)
    mpz_set_ui(P->X, 2);
  else
    mpz_set(P->X, X);
  mpz_set(P->Y, Y);
  mpz_set_ui(P->Z, 1);

  /* calcurate d if atkin_flag isn't set*/
 //モンゴメリー曲線の計算
 //  calc montgomery_b = (x^3+ax^2+x)/y^2
  if (X == NULL) {	
    mpz_pow_ui(tmp,P->X,2); //tmp = x^2
    mpz_mod(tmp,tmp,N);
    mpz_mul_mod(tmp2,P->X,montgomery_a,N); // tmp2 => a*x  
    mpz_add(tmp,tmp,tmp2);	//tmp= x^2+a*x
    mpz_add_ui(tmp,tmp,1);
    mpz_mul_mod(tmp,P->X,tmp,N); // tmp = x^3+ax^2+x	
    mpz_pow_ui(tmp2,P->Y,2); //tmp2 = y^2
    mpz_mod(tmp2,tmp2,N);
    mpz_invert(tmp2,tmp2,N);
    mpz_mul_mod(montgomery_b,tmp,tmp2,N); //montgomery_b = (x^3+ax^2+x)/y^2
  }

  /* set prime number */
  unsigned long int p = 2;
  mpz_t prime;
  mpz_init(prime);
  mpz_set_ui(prime, p);

  double stage1_time, stage2_time = -1;
  double start, end;
  start = omp_get_wtime();
  /* stage1 */
  while (p <= B1) {
    /* e = log p k */
    e = (int)(log(B1) / log(p));
    for (i = 1; i <= e; i++) {
      /* P_Z <- 1 */
      montgomery_scalar(P, P, p, montgomery_a, N);
      mpz_gcd(f, P->Z, N);
     /*	
      fprintf(fp,"prime number = %d time= %d \n",p,i);
      gmp_fprintf(fp,"montgomery_scalar_X = %Zd\n",P->X); 
      gmp_fprintf(fp,"montgomery_scalar_Y = %Zd\n",P->Y); 
      gmp_fprintf(fp,"montgomery_scalar_Z = %Zd\n",P->Z); 
      */
      if(mpz_cmp_ui(f,1) != 0&&mpz_cmp(f,N)!=0 ){
	gmp_printf("value of f = %Zd \n",f);        
	end = omp_get_wtime();
	stage1_time = end - start;
        //[WORNING]
        mpz_set_ui(f,1);
        //[WORNIG]
	goto FACTOR_FOUND;
      }
    }

    mpz_nextprime(prime, prime);
    p = mpz_get_ui(prime);
  }
  end = omp_get_wtime();
  stage1_time = end - start;
 
  start = omp_get_wtime();
  // stage2 
  printf("transition to the stage2");
  montgomery_bsgs(f,P,B1,B2,window_size,N,fp,d,montgomery_a,montgomery_b);
  end = omp_get_wtime();
  stage2_time = end - start;
 

 FACTOR_FOUND:
  fprintf(fp, "Stage1 time: %f seconds\n", stage1_time);
  if (stage2_time != -1){
    fprintf(fp, "Stage2 time: %f seconds\n", stage2_time);
    fprintf(fp, "Stage2 gcd: %ld\n", mpz_get_ui(f));
  }
  else
    fprintf(fp, "Stage2 time: ----\n");

  projective_point_clear(P);
  mpz_clear(tmp);
  mpz_clear(tmp2);
  mpz_clear(inv);
  mpz_clear(prime);
}

