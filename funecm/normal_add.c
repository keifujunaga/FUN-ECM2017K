#include <stdio.h>
#include "gmp.h"
#include "point.h"
 void
extended_dedicated_add (EXTENDED_POINT R, EXTENDED_POINT P, EXTENDED_POINT Q,
			const mpz_t N) 
{
  mpz_t A, B, C, D, E, F, G, H, tmp;
  mpz_inits (A, B, C, D, E, F, G, H, tmp, NULL);
  
    /* A<-(Y1-X1)*(Y2+X2) */ 
    mpz_sub (A, P->Y, P->X);
  mpz_add (tmp, Q->Y, Q->X);
  mpz_mul_mod (A, A, tmp, N);
  
    /* B<-(Y1+X1)*(Y2-X2) */ 
    mpz_add (B, P->Y, P->X);
  mpz_sub (tmp, Q->Y, Q->X);
  mpz_mul_mod (B, B, tmp, N);
  
    /* C<-2*Z1*T2 */ 
    mpz_mul_ui (C, P->Z, 2);
  mpz_mul_mod (C, C, Q->T, N);
  
    /* D<-2*T1*Z2 */ 
    mpz_mul_ui (D, P->T, 2);
  mpz_mod (D, D, N);
  mpz_mul_mod (D, D, Q->Z, N);
  
    /* E<-D+C */ 
    mpz_add (E, D, C);
  mpz_mod (E, E, N);
  
    /* F<-B-A */ 
    mpz_sub (F, B, A);
  mpz_mod (F, F, N);
  
    /* G<-B+A */ 
    mpz_add (G, B, A);
  mpz_mod (G, G, N);
  
    /* H<-D-C */ 
    mpz_sub (H, D, C);
  mpz_mod (H, H, N);
  
    /* X3<-E*F */ 
    mpz_mul_mod (R->X, E, F, N);
  
    /* Y3<-G*H */ 
    mpz_mul_mod (R->Y, G, H, N);
  
    /* T3<-E*H */ 
    mpz_mul_mod (R->T, E, H, N);
  
    /* Z3<-F*G */ 
    mpz_mul_mod (R->Z, F, G, N);
  mpz_clears (A, B, C, D, E, F, G, H, tmp, NULL);
} void

montgomery_add (PROJECTIVE_POINT R, PROJECTIVE_POINT P, PROJECTIVE_POINT Q,
		PROJECTIVE_POINT P0, const mpz_t N) 
{				//montgomery曲線の加算関数
  //P0は初期座標を表している
  //Qは足す座標を示している
  //Rは格納する座標
  //printf("\nnormal\n"); 
  mpz_t A, B, C, D, E, F, G, H, GG, HH;
  mpz_inits (A, B, C, D, E, F, G, H, GG, HH, NULL);
  
    /* A = P->X + P->Z 
       = xp+zp */ 
    mpz_add (A, P->X, P->Z);
  
    //gmp_printf("Xp+Zp=%Zd\n",A);
    /* B = P->X - P->Z 
       = xp-zp */ 
    mpz_sub (B, P->X, P->Z);
  
    //gmp_printf("Xp-Zp=%Zd\n",B);
    /* C = Q->X + Q->Z 
       = xQ+zQ */ 
    mpz_add (C, Q->X, Q->Z);
  
    //gmp_printf("XQ+ZQ=%Zd\n",C);
    /* D = Q->X - Q->Z 
       = xQ-zQ */ 
    mpz_sub (D, Q->X, Q->Z);
  
    //gmp_printf("XQ-ZQ=%Zd\n",D);
    /* E = D * A 
       = (xQ-zQ)*(xP+zP) */ 
    mpz_mul_mod (E, D, A, N);
  
    //gmp_printf("(XQ-ZQ)*(Xp+Zp)=%Zd\n",E);
    /* F = C * B 
       = (xQ+zQ)*(xP-zP) */ 
    mpz_mul_mod (F, C, B, N);
  
    //gmp_printf("(XQ+ZQ)*(Xp-Zp)=%Zd\n",F);
    /* G = E + F = DA + CB 
       = (xQ-zQ)*(xP+zP)+(xQ+zQ)*(xP-zP) */ 
    mpz_add (G, E, F);
  
    //gmp_printf("(XQ-ZQ)*(Xp+Zp)+(XQ+ZQ)*(Xp-Zp)=%Zd\n",G);
    /* H = E - F = DA - CB 
       = (xQ-zQ)*(xP+zP)-(xQ+zQ)*(xP-zP) */ 
    mpz_sub (H, E, F);
  
    //gmp_printf("(xQ-zQ)*(xP+zP)-(xQ+zQ)*(xP-zP)=%Zd\n",H);
    /* GG = G ^ 2 = (DA + CB) ^ 2 
       = {(xQ-zQ)*(xP+zP)+(xQ+zQ)*(xP-zP)}^2 */ 
    mpz_pow_ui (GG, G, 2);
  
    //gmp_printf("{(xQ-zQ)*(xP+zP)+(xQ+zQ)*(xP-zP)}^2=%Zd\n",GG);
    mpz_mod (GG, GG, N);
  
    /* HH = H ^ 2 = (DA - CB) ^ 2 
       = {(xQ-zQ)*(xP+zP)-(xQ+zQ)*(xP-zP)}^2 */ 
    mpz_pow_ui (HH, H, 2);
  
    //gmp_printf("{(xQ-zQ)*(xP+zP)-(xQ+zQ)*(xP-zP)}^2=%Zd\n",HH);
    mpz_mod (HH, HH, N);
  
    /* 
       R->X = O->Z * GG = O->Z * (DA + CB) ^ 2 
       = z1*{(xQ-zQ)*(xP+zP)+(xQ+zQ)*(xP-zP)}^2
     */ 
    mpz_mul_mod (R->X, P0->Z, GG, N);
  
    //gmp_printf("Xp-q=%Zd\n",P0->X);
    //gmp_printf("X=%Zd\n",R->X);
    /*　
       R->Z = O->X * HH = O->X * (DA - CB) ^ 2
       =x1*{(xQ-zQ)*(xP+zP)-(xQ+zQ)*(xP-zP)}^2
     */ 
    mpz_mul_mod (R->Z, P0->X, HH, N);
  
    //gmp_printf("Zp-q=%Zd\n",P0->Z);
    //gmp_printf("Z=%Zd\n",R->Z);
    mpz_clears (A, B, C, D, E, F, G, H, GG, HH, NULL);
} 
