#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "point.h"

int main(){
  PROJECTIVE_POINT P;
  PROJECTIVE_POINT R;
  mpz_t tmp;
  mpz_t tmp2;
  mpz_t inv;
  mpz_t a;//モンゴメリー曲線の座標 
  mpz_t b;
  mpz_t N;
	
  mpz_init(tmp);
  mpz_init(tmp2);
  mpz_init(inv); 
  mpz_init(a);
  mpz_init(b);
  mpz_init (N);

  projective_point_init(P);
  projective_point_init(R);
  mpz_set_ui(a,11);
  mpz_set_ui(P->Y,3);
  mpz_set_ui(P->X, 2);
  mpz_set_ui(P->Z, 1);
  mpz_set_ui(N,3925865844970303603);   


  mpz_set_ui(R->Y,3);
  mpz_set_ui(R->X, 2);
  mpz_set_ui(R->Z, 1);

  printf("Montgomery_a = ");
  mpz_out_str(stdout,10,a);
  printf("\n");
  gmp_printf("P->X = %Zd\n",P->X);
  gmp_printf("P->Y = %Zd\n",P->Y);
  gmp_printf("P->Z = %Zd\n",P->Z);

  montgomery_double(P,P,a,N);
  printf("after_double(double)\n");
  gmp_printf("P->X = %Zd\n",P->X);
  gmp_printf("P->Y = %Zd\n",P->Y);
  gmp_printf("P->Z = %Zd\n",P->Z);
  
  /* 
  montgomery_add(P,P,P,P);
  printf("after_add\n");
  gmp_printf("P->X = %Zd\n",P->X);
  gmp_printf("P->Y = %Zd\n",P->Y);
  gmp_printf("P->Z = %Zd\n",P->Z);
  */

  montgomery_double(P,P,a,N);
  printf("after_double(fourfold)\n");
  gmp_printf("P->X = %Zd\n",P->X);
  gmp_printf("P->Y = %Zd\n",P->Y);
  gmp_printf("P->Z = %Zd\n",P->Z);

  montgomery_scalar(R,P,4,a,N);
  printf("after_scalar\n");
  gmp_printf("P->X = %Zd\n",P->X);
  gmp_printf("P->Y = %Zd\n",P->Y);
  gmp_printf("P->Z = %Zd\n",P->Z);

return 0;


}
