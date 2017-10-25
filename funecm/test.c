#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "point.h"



int
main ()
{
  PROJECTIVE_POINT P;
  PROJECTIVE_POINT R;
  PROJECTIVE_POINT T;
  PROJECTIVE_POINT I;
  AFFINE_POINT A;
 
  mpz_t tmp;
  mpz_t tmp2;
  mpz_t inv;
  mpz_t a;			//モンゴメリー曲線の座標 
  mpz_t b;
  mpz_t N;
  mpz_t s;

  mpz_init (tmp);
  mpz_init (tmp2);
  mpz_init (inv);
  mpz_init (a);
  mpz_init (b);
  mpz_init (N);
  mpz_init (s);

  projective_point_init (P);
  projective_point_init (R);
  projective_point_init (T);
  projective_point_init (I);
  affine_point_init(A);  

  mpz_set_ui (N, 311119);
  mpz_set_ui (tmp, 13);
  mpz_set_ui (inv, 16);
  mpz_invert (inv, inv, N);
  mpz_mul (a, tmp, inv);
  mpz_set_ui (P->X, 4);
  mpz_set_ui (P->Z, 1);
  mpz_set_ui (R->X, 4);
  mpz_set_ui (R->Z, 1);

  //printf("Montgomery_a = ");
  //mpz_out_str(stdout,10,a);//10は10進数表記をするという意味
  //printf("\n");

  gmp_printf ("P->X = %Zd\n", P->X);
  gmp_printf ("P->Z = %Zd\n", P->Z);
/*
  printf ("\n");
  printf ("after montgomery_double\n");
  montgomery_double (R, R, a, N);
  print_p(R);
  //montgomery_double (I,R,a,N);//1043435
  //montgomery_double (I,I,a,N);
  montgomery_add(T,R,P,P,N);
  printf("T\n");
  print_p(T);
  montgomery_add(I,T,R,P,N);
  printf("I\n");
  print_p(I);
*/
  printf ("\n");
  printf ("after_scalar\n");
  montgomery_scalar (P, P,5, a, N);
  //printf ("after_scalar\n");
  print_p(P);
  protoaff(A,P,N);
  gmp_printf("X=%Zd\n",A->x);   
  return 0;


}
