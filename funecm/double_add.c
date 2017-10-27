#include <stdio.h>
#include "gmp.h"
#include "point.h"
 void
double_add (PROJECTIVE_POINT R, PROJECTIVE_POINT P, const mpz_t N) 
{
  mpz_t B, C, D, E, F, H, J;
  mpz_inits (B, C, D, E, F, H, J, NULL);
  mpz_add (B, P->X, P->Y);	//B = X1+Y1
  mpz_pow_ui (B, B, 2);		//B = (X1+Y1)^2
  mpz_mod (B, B, N);
  mpz_pow_ui (C, P->X, 2);	//C = X1^2
  mpz_mod (C, C, N);
  mpz_pow_ui (D, P->Y, 2);	//D = Y1^2 
  mpz_mod (D, D, N);
  mpz_mul_si (E, C, -1);	//E=-C
  mpz_add (F, E, D);		//F = E+D
  mpz_pow_ui (H, P->Z, 2);	//H = Z1^2
  mpz_mod (H, H, N);
  mpz_sub (J, F, H);		//J = F-H
  mpz_sub (J, J, H);		//J = F-2H
  
    /* calcurate X_3 */ 
    mpz_sub (R->X, B, C);	//X3=B-C
  mpz_sub (R->X, R->X, D);	//X3=B-C-D
  mpz_mul_mod (R->X, R->X, J, N);	//X3
  
    /* calcurate Y_3 */ 
    mpz_sub (R->Y, E, D);	//Y3 = E-D
  mpz_mul_mod (R->Y, R->Y, F, N);	//Y3
  
    /* calcurate Z_3 */ 
    mpz_mul_mod (R->Z, F, J, N);	//Z3
  mpz_clears (B, C, D, E, F, H, J, NULL);
} void

dedicated_doubling (EXTENDED_POINT R, const EXTENDED_POINT P, const mpz_t N) 
{
  mpz_t A, B, C, D, E, G, F, H;
  mpz_inits (A, B, C, D, E, G, F, H, NULL);
  mpz_pow_ui (A, P->X, 2);
  mpz_mod (A, A, N);
  mpz_pow_ui (B, P->Y, 2);
  mpz_mod (B, B, N);
  mpz_pow_ui (C, P->Z, 2);
  mpz_mod (C, C, N);
  mpz_mul_ui (C, C, 2);
  mpz_mod (C, C, N);
  mpz_mul_si (D, A, -1);
  mpz_add (E, P->X, P->Y);
  mpz_pow_ui (E, E, 2);
  mpz_mod (E, E, N);
  mpz_sub (E, E, A);
  mpz_sub (E, E, B);
  mpz_add (G, D, B);
  mpz_sub (F, G, C);
  mpz_sub (H, D, B);
  mpz_mul_mod (R->X, E, F, N);
  mpz_mul_mod (R->Y, G, H, N);
  mpz_mul_mod (R->T, E, H, N);
  mpz_mul_mod (R->Z, F, G, N);
  mpz_clears (A, B, C, D, E, G, F, H, NULL);
} 

//calc montgomery double （2倍算）
  void
montgomery_double (PROJECTIVE_POINT R, PROJECTIVE_POINT P,
		   const mpz_t montgomery_a, const mpz_t N)
{
  
    //PがX1とx2の座標を示す
    mpz_t A, AA, B, BB, C, D, E, F, G, four, inv;
  mpz_inits (A, AA, B, BB, C, D, E, F, G, four, inv, NULL);
  mpz_set_ui (four, 4);
  
    //printf("\ndouble\n");
    // A = P->X + P->Z = xp+zp
    mpz_add (A, P->X, P->Z);
  
    //gmp_printf("Xp+Zp=%Zd\n",A); 
    // AA = (A * A) mod N = (xp+zp)^2
    mpz_mul_mod (AA, A, A, N);
  
    //gmp_printf("(Xp+Zp)^2=%Zd\n",AA); 
    // B = P->X - P->Z = x1-z1
    mpz_sub (B, P->X, P->Z);
  
    //gmp_printf("Xp-Zp=%Zd\n",B); 
    // BB = (B * B) mod N = (xp-zp)^2
    mpz_mul_mod (BB, B, B, N);
  
    //gmp_printf("(Xp-Zp)^2=%Zd\n",BB); 
    // C = AA - BB = (Xp+Zp)^2-(Xp-Zp)^2
    //4xz
    mpz_sub (C, AA, BB);
  
    //gmp_printf("4XpZp=%Zd\n",C); 
    // R->X3 = (AA * BB) mod N = (Xp+Zp)^2*(Xp-Zp)^2
    //calc X_point
    mpz_mul_mod (R->X, AA, BB, N);
  
    //gmp_printf("X=%Zd\n",R->X); 
    // D = a + 2 
    mpz_add_ui (D, montgomery_a, 2);
  
    //gmp_printf("a+2=%F\n",D); 
    // E = D / 4 = (a+2)/4
    mpz_invert (inv, four, N);
  mpz_mul (E, D, inv);
  
    //gmp_printf("(a+2)/4=%.*Ff\n",E); 
    // F = E * C = E * (AA - BB)
    //           = {(a+2)/4}*{(x1+z1)^2-(x1-z1)^2}
    // a24*4xz
    mpz_mul_mod (F, E, C, N);
  
    //gmp_printf("a24*4XpZp=%Ff\n",F); 
    // G = F + BB = E * C + BB
    //            = E * (AA - BB) + BB 
    //            = {(a+2)/4}*{(x1+z1)^2-(x1-z1)^2}+(x1-z1)^2
    // (x-z)^2+a24*4xz
    mpz_add (G, F, BB);
  
    //gmp_printf("(Xp-Zp)^2+a24*4XpZp=%Zd\n",G); 
    // R->Z = C * G = C * (F + BB)
    //              = C * (E * C + BB)
    //              = {(x1+z1)^2-(x1-z1)^2} *
    //                        [{(a+2)/4}*{(x1+z1)^2-(x1-z1)^2}+(x1-z1)^2]
    mpz_mul_mod (R->Z, C, G, N);
  
    //gmp_printf("Z=%Zd\n",R->Z); 
    mpz_clears (A, AA, B, BB, C, D, E, F, G, four, inv, NULL);
} 
