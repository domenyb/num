#include<math.h>
#include<gsl/gsl_vector.h>


double rosen_func(const gsl_vector *v, void *params)
{
  double x, y;

  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);

  return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}


struct experimental_data {int n; double *t, *y, *e;};

double decayfunction(double tim, double A, double T, double B){
	return A*exp(-tim/T)+B;
}


double function_to_minimize(const gsl_vector *x, void *params){
	double A = gsl_vector_get(x,0);
	double T = gsl_vector_get(x,1);
	double B = gsl_vector_get(x,2);
	struct experimental_data p = *(struct experimental_data *)params;
	double chisq=0;
	for(int i=0;i<p.n;i++){
		double tim=p.t[i];
		double sig=p.y[i];
		double err=p.e[i];
		chisq+=pow(decayfunction(tim,A,T,B)-sig,2)/err/err;
		}
	return chisq;
}
