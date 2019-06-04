#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bases.h"
#include <assert.h>

double adapt24(double f(double), double a, double b, double acc, double eps, double *error, double f2, double f3, int nrec){ //fro m the lecture notes
  assert(nrec < 1e6);
  double f1=f(a+(b-a)/6.);
  double f4=f(a+5.*(b-a)/6.);
  double Q=(2*f1+ f2+f3+2*f4)/6.0*(b-a);
  double q=(f1+f2+f3+f4)/4.0*(b-a);
  double tolerance = acc + eps*fabs(Q);
  *error = fabs(Q-q);
  if (*error < tolerance) {
    return Q;
  }
  else{
    double Q1 = adapt24(f,a,(a+b)/2.,acc/sqrt(2.),eps,error,f1,f2,nrec+1);
    double Q2 = adapt24(f,(a+b)/2.,b,acc/sqrt(2.),eps,error,f3,f4,nrec+1);
    return Q1 + Q2;
  }
}

double adapt(double f(double), double a, double b, double acc, double eps, double* err){ //from the lecture notes
  double f2 = f(a+2.*(b-a)/6.);
  double f3 = f(a+4.*(b-a)/6.);
  int nrec = 0;
  return adapt24(f,a,b,acc,eps, err,f2,f3,nrec);
}



double adapt_CC(double f(double), double a, double b, double acc, double eps, double* err)
{
  //let's try this nested function and it might just work
  double cc(double fi){
  return (f((a+b)/2+(a-b)/2*cos(fi))) * sin(fi)*(b-a)/2;
  }
  double a2 = 0;
  double b2 = 3.141592653589793238462;
  return adapt(cc, a2, b2, acc, eps, err);
}
