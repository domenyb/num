#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_odeiv2.h>
#include<math.h>
#include<assert.h>
#define EPS 1e-12

int ode_H(double x, const double y[], double dydx[], void *params)
{
	double e=*(double*)params;
	dydx[0] = y[1];
	dydx[1] = 2*(-1/x-e)*y[0];
	return GSL_SUCCESS;
}

double fe(double e,double maxx)
{ assert(maxx>=0);
	gsl_odeiv2_system system;
	system.function = ode_H;
	system.jacobian = NULL;
	system.dimension = 2;
	system.params = (void*)&e;

	gsl_odeiv2_driver *driver;
	double hstart = 0.001;
	driver = gsl_odeiv2_driver_alloc_y_new(&system,
					       gsl_odeiv2_step_rkf45,
					       hstart, 1e-8, 1e-8);
	double x0 = 0.0004;
	double y[] = { x0-x0*x0, 1-2*x0};
	gsl_odeiv2_driver_apply(driver, &x0, maxx, y);
	gsl_odeiv2_driver_free(driver);
	return y[0];
}

int equation (const gsl_vector * x, void * params, gsl_vector * f)
{
	double maxx = *(double*)params;
	double e = gsl_vector_get(x,0);
	double mm = fe(e,maxx);
	gsl_vector_set(f,0,mm);
	return GSL_SUCCESS;
}

int root_equation(const gsl_vector * x, void * params, gsl_vector * f){
  const double r = gsl_vector_get(x,0);
  const double r2 = gsl_vector_get(x,1);
	double mm =  2* ( 200 * r*r*r - 200 * r * r2 + r -1 );
	const double mm2 =  200*(r2 - r*r);


	gsl_vector_set(f,0,mm);
	gsl_vector_set(f,1,mm2);
	return GSL_SUCCESS;
}

int root(double * z,double * z2){
	gsl_multiroot_function F;
	F.f=root_equation;
	F.n=2;
	F.params=NULL;

	gsl_multiroot_fsolver * S;
	S = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids,F.n);
	gsl_vector* start = gsl_vector_alloc(F.n);
  gsl_vector_set (start, 0, *z);
  gsl_vector_set (start, 1, *z2);
	gsl_multiroot_fsolver_set(S,&F,start);


	int flag,iter=0;
	flag=GSL_CONTINUE;
	while(flag==GSL_CONTINUE){
		iter++;
		gsl_multiroot_fsolver_iterate(S);
		flag=gsl_multiroot_test_residual(S->f,EPS);
  fprintf(stdout,"x= %g, \t  y= %g \t  iter: %i \t   the gradient: %g \t %g\n",gsl_vector_get(S->x,0),gsl_vector_get(S->x,1),iter,gsl_vector_get(S->f,0),gsl_vector_get(S->f,1));
	};
  *z = gsl_vector_get(S->x,0);
  *z2 = gsl_vector_get(S->x,1);



	gsl_multiroot_fsolver_free(S);
	gsl_vector_free(start);
	return GSL_SUCCESS;
}
