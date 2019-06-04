#include "integration.h"
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_blas.h>

double my_function(double x){return x*x;}//according to the excersise, we should use  nested functions

int main(int argc, char** argv){
  double f1(double x){return sqrt(x);};
  double err=0;//we would need this to be added to be able to get it out
  double a=0;
  double b=1;
  double acc=0.0001;
  double eps=0.0001;
  double sol1=adapt(f1, a, b, acc, eps, &err);
  printf("First integral from 0-1 sqrt(x) supposed to be = 2/3\n" );
  printf("The integral's value is: %g +/- %g \n",sol1, err );
  printf("now try it with the CC transformation:\n" );
  err=0;//we would need this to be added to be able to get it out
  a=0;
  b=1;
  acc=0.0001;
  eps=0.0001;
  sol1=adapt_CC(f1, a, b, acc, eps, &err);
  printf("The integral's value now is: %g +/- %g \n",sol1, err );
  printf("with GSL routines: \n" );

  int lt=100;
  gsl_integration_workspace * w;
  w = gsl_integration_workspace_alloc (lt);
  double f1_gsl(double x,void * params){ return sqrt(x);};
  gsl_function  F;
  F.function = f1_gsl;
  gsl_integration_qags(&F, 0, 1, acc, eps, lt, w, &sol1, &err);
  gsl_integration_workspace_free(w);
  printf("The integral's value now is: %g +/- %g \n",sol1, err );

  double f2(double x){return 1.0/sqrt(x);};
   err=0;//we would need this to be added to be able to get it out
   a=0;
   b=1;
   acc=0.0001;
   eps=0.0001;
  double sol2=adapt(f2, a, b, acc, eps, &err);
  printf("Second integral from 0-1 1/sqrt(x) is supposed to be = 2\n" );
  printf("The integral's value is: %g +/- %g \n",sol2, err );
  printf("now try it with the CC transformation:\n" );
  err=0;//we would need this to be added to be able to get it out
  a=0;
  b=1;
  acc=0.0001;
  eps=0.0001;
  sol2=adapt_CC(f2, a, b, acc, eps, &err);
  printf("The integral's value now is: %g +/- %g \n",sol2, err );
  printf("with GSL routines: \n" );

  lt=100;
  gsl_integration_workspace * w2;
  w2 = gsl_integration_workspace_alloc (lt);
  double f2_gsl(double x,void * params){ return 1.0/sqrt(x);};
  F.function = f2_gsl;
  gsl_integration_qags(&F, 0, 1, acc, eps, lt, w2, &sol2, &err);
  gsl_integration_workspace_free(w2);
  printf("The integral's value now is: %g +/- %g \n",sol2, err );




  double f3(double x){return log(x)/sqrt(x);};
  err=0;//we would need this to be added to be able to get it out
  a=0;
  b=1;
  acc=0.0001;
  eps=0.0001;
  double sol3=adapt(f3, a, b, acc, eps, &err);
  printf("third integral from 0-1 ln(x)/sqrt(x) supposed to be = -4\n" );
  printf("The integral's value is: %g +/- %g \n",sol3, err );
  printf("now try it with the CC transformation:\n" );
  err=0;//we would need this to be added to be able to get it out
  a=0;
  b=1;
  acc=0.0001;
  eps=0.0001;
  sol3=adapt_CC(f3, a, b, acc, eps, &err);
  printf("The integral's value now is: %g +/- %g \n",sol3, err );
  printf("with GSL routines: \n" );

  lt=100;
  gsl_integration_workspace * w3;
  w3 = gsl_integration_workspace_alloc (lt);
  double f3_gsl(double x,void * params){ return log(x)/sqrt(x);};
  F.function = f3_gsl;
  gsl_integration_qags(&F, 0, 1, acc, eps, lt, w, &sol3, &err);
  gsl_integration_workspace_free(w3);
  printf("The integral's value now is: %g +/- %g \n",sol3, err );




  double f4(double x){return 4.*sqrt(1.0-(1.0-x)*(1.0-x));};
  err=0;//we would need this to be added to be able to get it out
  a=0;
  b=1;
  acc=0.0001;
  eps=0.0001;
  double sol4=adapt(f4, a, b, acc, eps, &err);
  printf("our last integral from 0-1  4sqrt(1-(1-x)2) supposed to be = pi (3.141592653589793238462...)\n" );
  printf("The integral's value is: %g +/- %g \n",sol4, err );
  printf("now try it with the CC transformation:\n" );
  err=0;//we would need this to be added to be able to get it out
  a=0;
  b=1;
  acc=0.0001;
  eps=0.0001;
  sol4=adapt_CC(f4, a, b, acc, eps, &err);
  printf("The integral's value now is: %g +/- %g \n",sol4, err );
  printf("with GSL routines: \n" );

  lt=100;
  gsl_integration_workspace * w4;
  w4 = gsl_integration_workspace_alloc (lt);
  double f4_gsl(double x,void * params){return 4.*sqrt(1.0-(1.0-x)*(1.0-x));};
  F.function = f4_gsl;
  gsl_integration_qags(&F, 0, 1, acc, eps, lt, w, &sol4, &err);
  gsl_integration_workspace_free(w4);
  printf("The integral's value now is: %g +/- %g \n",sol4, err );

printf("this is succesfull so far... Maybe i'll return to implement excersise C, since this doesn't seem to be as time consuming\n ----------------------------------------------------------------------\n" );

//I don't need to free anything since I didn't allocate anything nor did I open files...


  return 0;
}
