#include <stdio.h>
#include "gmp.h"
#include <time.h>
#include "point.h"

/* convert affine point to projective point */ 
  void
afftopro (PROJECTIVE_POINT R, const AFFINE_POINT P) 
{
  mpz_set (R->X, P->x);
  mpz_set (R->Y, P->y);
  mpz_set_ui (R->Z, 1);
} 

/* convert projective point to affine point */ 
  void
protoaff (AFFINE_POINT R, const PROJECTIVE_POINT P, const mpz_t N) 
{
  mpz_t inv;
  mpz_init (inv);
  mpz_invert (inv, P->Z, N);
  mpz_mul (R->x, P->X, inv);
  mpz_mul (R->y, P->Y, inv);
  mpz_mod (R->x, R->x, N);
  mpz_mod (R->y, R->y, N);
  mpz_clear (inv);
} 

/* convert projective point to extended point */ 
  void
protoext (EXTENDED_POINT R, const PROJECTIVE_POINT P, const mpz_t N) 
{
  
    /* (X:Y:Z) -> (XZ:YZ:XY:Z^2) */ 
    mpz_mul (R->X, P->X, P->Z);
  mpz_mod (R->X, R->X, N);
  mpz_mul (R->Y, P->Y, P->Z);
  mpz_mod (R->Y, R->Y, N);
  mpz_mul (R->T, P->X, P->Y);
  mpz_mod (R->T, R->T, N);
  mpz_pow_ui (R->Z, P->Z, 2);
  mpz_mod (R->Z, R->Z, N);
} 

/* convert extended point to projective point */ 
  void
exttopro (PROJECTIVE_POINT R, const EXTENDED_POINT P) 
{
  mpz_set (R->X, P->X);
  mpz_set (R->Y, P->Y);
  mpz_set (R->Z, P->Z);
} void

montgomery_coefficient (mpz_t A, mpz_t B, const mpz_t d, const mpz_t N) 
{				//montgomery曲線の係数A,Bをedward曲線のdから算出
  mpz_t C, D, E, F, inv;
  mpz_inits (C, D, E, F, inv, NULL);
  
    //C = d - 1
    mpz_sub_ui (C, d, 1);
  
    //D = d + 1
    mpz_add_ui (D, d, 1);
  
    //E = D * -1 = -d-1
    mpz_mul_si (E, D, -1);
  
    //inv = 1 / E = 1/(-d-1)
    mpz_invert (inv, E, N);
  
    //F = C * 2 = 2*(d-1)
    mpz_mul_ui (F, C, 2);
  
    //A = F * inv = 2*(d-1)/(-d-1)
    mpz_mul_mod (A, F, inv, N);
  
    //B = E * 4 = 4*(-d-1)
    mpz_mul_ui (B, inv, 4);
  mpz_clears (C, D, E, F, inv, NULL);
} void

exttomon (MONTGOMERY_POINT R, const EXTENDED_POINT P, const mpz_t N) 
{				//Edwards曲線でのX,Y座標をmontgomery曲線でのX,Z座標に変換
  mpz_t A, B, C, inv;
  mpz_inits (A, B, C, inv, NULL);
  
    //A = y + 1 
    mpz_add_ui (A, P->Y, 1);
  
    //B = y *(-1)
    mpz_mul_si (B, P->Y, -1);
  
    //C = - y + 1
    mpz_add_ui (C, B, 1);
  
    //inv = 1 / C = 1 / (1-y)
    mpz_invert (inv, C, N);
  
    //R->X = A * inv = (y+1)*{1/(1-y)} = (y+1)/(1-y)
    mpz_mul_mod (R->X, A, inv, N);
  
    //R->Z = 1
    mpz_set_ui (R->Z, 1);
  mpz_clears (A, B, C, inv, NULL);
} 
