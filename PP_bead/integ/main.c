#include<stdio.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>
#include<math.h>

double integrand (double t, void * params) {
  double f = log(t)/sqrt(t);
	return f;
}

double integrandnorm (double t, void * params) {
  double alpha = *(double*)params;
  double f = exp(-alpha * t*t);
  return f;
}
double integrandham (double t, void * params) {
  double alpha = *(double*)params;
  double f = (-alpha * alpha * t*t/2 + alpha/2 + t * t/2) * exp(-alpha * t*t);
  return f;
}


int main(){
    int lim=100;
  	gsl_integration_workspace * w;
  	w = gsl_integration_workspace_alloc (lim);
    gsl_function F;
  	 F.function = integrand;
    double result,error;
    int flag = gsl_integration_qags
    		(&F, 0, 1, 1e-8, 1e-8, lim, w, &result, &error);
    gsl_integration_workspace_free(w);
//    if(flag!=GSL_SUCCESS)  printf("something is wrong?");
    printf("first result = %g\n",result);
//now b part

gsl_integration_workspace * w2;
w2 = gsl_integration_workspace_alloc (lim);
gsl_function F2;
F2.function = integrandnorm;

gsl_integration_workspace * w3;
w3 = gsl_integration_workspace_alloc (lim);
gsl_function F3;
F3.function = integrandham;


for(double alpha=0.1;alpha<20;alpha+=0.1){

  F2.params = (void*)&alpha;
	double result,error,acc2=1e-8,eps2=1e-8;
  flag = gsl_integration_qagi
		(&F2, acc2, eps2, lim, w2, &result, &error);
//flag not needed, it works

  F3.params = (void*)&alpha;
  double result2,error2,acc3=1e-8,eps3=1e-8;
  flag = gsl_integration_qagi
    (&F3, acc3, eps3, lim, w3, &result2, &error2);
//flag not needed, it works

  fprintf(stderr,"%g %g %g %g\n",alpha,result,result2,result2/result);
}
gsl_integration_workspace_free(w2);
gsl_integration_workspace_free(w3);
printf("The minimum we found is the ground state\n");
return 0;
}
