#include <stdio.h>
#include "gmp.h"
#include "point.h"

void extended_dedicated_add(EXTENDED_POINT R, EXTENDED_POINT P, EXTENDED_POINT Q, const mpz_t N)
{
	mpz_t A,B,C,D,E,F,G,H,tmp;
	
	mpz_inits(A,B,C,D,E,F,G,H,tmp,NULL);
	
	/* A<-(Y1-X1)*(Y2+X2) */
	mpz_sub(A, P->Y, P->X);
	mpz_add(tmp, Q->Y, Q->X);
	mpz_mul_mod(A, A, tmp, N);

	/* B<-(Y1+X1)*(Y2-X2) */
	mpz_add(B, P->Y, P->X);
	mpz_sub(tmp, Q->Y, Q->X);
	mpz_mul_mod(B, B, tmp, N);

	/* C<-2*Z1*T2 */
	mpz_mul_ui(C, P->Z, 2);
	mpz_mul_mod(C, C, Q->T, N);

	/* D<-2*T1*Z2 */
	mpz_mul_ui(D, P->T, 2);
	mpz_mod(D, D, N);
	mpz_mul_mod(D, D, Q->Z, N);

	/* E<-D+C */
	mpz_add(E, D, C);
	mpz_mod(E, E, N);
	
	/* F<-B-A */
	mpz_sub(F, B, A);
	mpz_mod(F, F, N);

	/* G<-B+A */
	mpz_add(G, B, A);
	mpz_mod(G, G, N);
	
	/* H<-D-C */
	mpz_sub(H, D, C);
	mpz_mod(H, H, N);

	/* X3<-E*F */
	mpz_mul_mod(R->X, E, F, N);
	/* Y3<-G*H */
	mpz_mul_mod(R->Y, G, H, N);
	/* T3<-E*H */
	mpz_mul_mod(R->T, E, H, N);
	/* Z3<-F*G */
	mpz_mul_mod(R->Z, F, G, N);

	mpz_clears(A,B,C,D,E,F,G,H,tmp,NULL);
}

void montgomery_add (PROJECTIVE_POINT R, PROJECTIVE_POINT P, 
		     PROJECTIVE_POINT Q, PROJECTIVE_POINT O, const mpz_t N)
{//montgomery曲線の加算関数
  mpz_t A,B,C,D,E,F,G,H,GG,HH;
  mpz_inits(A,B,C,D,E,F,G,H,GG,HH,NULL);
  /* A = P->X + P->Z 
       = x2+z2 */
  mpz_add(A,P->X,P->Z);
  /* B = P->X - P->Z 
       = x2-z2 */
  mpz_sub(B,P->X,P->Z);
  /* C = Q->X + Q->Z 
       = x3+z3 */
  mpz_add(C,Q->X,Q->Z);
  /* D = Q->X - Q->Z 
       = x3-z3 */
  mpz_sub(D,Q->X,Q->Z);
  /* E = D * A 
       = (x3-z3)*(x2+z2) */
  mpz_mul_mod(E,D,A,N);
  /* F = C * B 
       = (x3+z3)*(x2-z2) */
  mpz_mul_mod(F,C,B,N);
  /* G = E + F = DA + CB 
               = (x3-z3)*(x2+z2)+(x3+z3)*(x2-z2) */
  mpz_add(G,E,F);
  /* H = E - F = DA - CB 
               = (x3-z3)*(x2+z2)-(x3+z3)*(x2-z2) */
  mpz_sub(H,E,F);
  /* GG = G ^ 2 = (DA + CB) ^ 2 
                = {(x3-z3)*(x2+z2)+(x3+z3)*(x2-z2)}^2 */
  mpz_pow_ui(GG,G,2);
  mpz_mod(GG,GG,N);
  /* HH = H ^ 2 = (DA - CB) ^ 2 
                = {(x3-z3)*(x2+z2)-(x3+z3)*(x2-z2)}^2 */
  mpz_pow_ui(HH,H,2);
  mpz_mod(HH,HH,N);
  /* 
     R->X = O->Z * GG = O->Z * (DA + CB) ^ 2 
                    = z1*{(x3-z3)*(x2+z2)+(x3+z3)*(x2-z2)}^2
   */
  mpz_mul_mod(R->X,O->Z,GG,N);
  /*　
    R->Z = O->X * HH = O->X * (DA - CB) ^ 2
                   =x1*{(x3-z3)*(x2+z2)-(x3+z3)*(x2-z2)}^2
  */
  mpz_mul_mod(R->Z,O->X,HH,N);
  mpz_clears(A,B,C,D,E,F,G,H,GG,HH,NULL);
} 
  
