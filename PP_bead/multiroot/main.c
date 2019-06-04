#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
#include"funcs.h"




int main(){
  double x = 1.5;
  double y = 1.5;

  double *px = &x;
  double *py = &y;

  printf("starting point: %g \t %g\n",*px,*py);
  int success=root(px,py);
  printf("result should be: 1 \t 1 \n it is: %g \t %g\n",*px,*py);



double xmax = 10;
double estart = -0.5;


gsl_multiroot_function F;
F.f=equation;
F.n=1;
F.params=(void*)&xmax;

gsl_multiroot_fsolver * S;
S = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids,F.n);

gsl_vector* start = gsl_vector_alloc(F.n);
gsl_vector_set(start,0,estart);
gsl_multiroot_fsolver_set(S,&F,start);

int flag;
do{
	gsl_multiroot_fsolver_iterate(S);
	flag=gsl_multiroot_test_residual(S->f,EPS);
}while(flag==GSL_CONTINUE);

double result=gsl_vector_get(S->x,0);
gsl_multiroot_fsolver_free(S);
gsl_vector_free(start);
printf("Now the hydrogen: \n expected:  -0.5 \n");
printf("result: %g\n",result);

for (size_t i = 1; i < 100; i++) {
  fprintf(stderr,"%g %g %g\n",xmax*i/100,fe(result,xmax*i/100),xmax*i/100*exp(-xmax*i/100));
}

return 0;
}
