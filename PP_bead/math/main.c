#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#define M_E  2.7182818284590452354 /* e */
#define M_PI 3.1415926535897932384 /* pi */


complex double gmm;
complex double bessel;
complex double sqrtval;
complex double eadi;
complex double eadpi;
complex double ipw;
float flt = 0.1111111111111111111111111111;
double dbl = 0.1111111111111111111111111111;
long double  ld = 0.1111111111111111111111111111L;




int main(){
  gmm = gamma(5);
  bessel = j1( 0.5);
  sqrtval = csqrt(-2);
  eadi = cexp(I);
  eadpi = cexp(I*M_PI);
  ipw = cpow(I, M_E);
  printf("Gamma(5) = %g + %g i\n" ,creal(gmm),cimag(gmm));
  printf("bessel (0.5) = %g + %g i\n",creal(bessel), cimag(bessel));
  printf("sqrt(2) = %g + %g i\n",creal(sqrtval), cimag(sqrtval));
  printf("e^i = %g + %g i\n",creal(eadi), cimag(eadi));
  printf("e^i*pi = %g + %g i\n",creal(eadpi), cimag(eadpi));
  printf("i^e = %g +%g i\n", creal(ipw), cimag(ipw));

  printf("this is a float %.25g\n",flt);
  printf("this is a double %.25lg\n",dbl);
  printf("this is a long double %.25Lg\n",ld);

  return 0;
}
