#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>


int ode_err(double t, const double y[], double dydt[], void* params){
  dydt[0]= (2/sqrt(M_PI))*exp(-t*t);
return GSL_SUCCESS;
}

double err_func(double t){
  //after the documentation
	if (t<0) {
    return -err_func(-t);
  }
  gsl_odeiv2_system sys;
  sys.function=ode_err;
  sys.jacobian=NULL;
  sys.dimension=1;
  sys.params=NULL;
  gsl_odeiv2_driver* driver;
  double hstart=0.0001;
  driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,hstart,1e-8,1e-8);

  double t0=0;
  double y[1]={0.0};
  gsl_odeiv2_driver_apply(driver,&t0,t,y);

  gsl_odeiv2_driver_free(driver);
  return y[0];
}



int main(int argc, char* argv[]){
  if (argc!=4){
    fprintf(stderr, "I need 3 command line arguments!\n" );
    return 1;
  }
	double a = atof(argv[1]);
	double b = atof(argv[2]);
	double dx = atof(argv[3]);

	for (double x = a; x <= b+1e-6; x += dx)
		fprintf(stdout, "%g %g %g\n", x, err_func(x), gsl_sf_erf(x));

return 0;
}
